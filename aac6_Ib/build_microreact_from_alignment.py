#!/usr/bin/env python3

import argparse
import base64
import csv
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


API_URL = "https://microreact.org/api/projects/create"
DEFAULT_TOKEN_ENV = "MICROREACT_ACCESS_KEY"
TOKEN_ENV_FALLBACKS = ("MICROREACT_ACCESS_KEY", "MICROREACT_ACCESS_TOKEN")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract variable positions from an aligned FASTA file, write a CSV, "
            "build a Microreact payload, and optionally create the project via API.\n\n"
            "To upload directly to Microreact, first put your access key in an "
            "environment variable, for example:\n"
            "  export MICROREACT_ACCESS_KEY='your-key-here'\n"
            "Then run this script normally. If you only want the CSV and "
            ".microreact file without uploading, add --skip-upload."
        )
        ,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("alignment", type=Path, help="Aligned FASTA file")
    parser.add_argument("tree", type=Path, help="Newick tree file")
    parser.add_argument(
        "--project-name",
        required=True,
        help="Project name to use in Microreact",
    )
    parser.add_argument(
        "--csv-output",
        type=Path,
        default=Path("variant_positions.csv"),
        help="Output CSV path",
    )
    parser.add_argument(
        "--payload-output",
        type=Path,
        default=Path("project.microreact"),
        help="Output .microreact path",
    )
    parser.add_argument(
        "--id-column",
        default="id",
        help="Identifier column name for the CSV",
    )
    parser.add_argument(
        "--reference-id",
        help=(
            "Optional reference sequence ID from the alignment. Columns where the "
            "reference has a gap are removed, and remaining columns are numbered "
            "by reference position."
        ),
    )
    parser.add_argument(
        "--token-env",
        default=DEFAULT_TOKEN_ENV,
        help=(
            "Environment variable containing the Microreact API token "
            f"(default: {DEFAULT_TOKEN_ENV}).\n"
            "Example:\n"
            f"  export {DEFAULT_TOKEN_ENV}='your-key-here'"
        ),
    )
    parser.add_argument(
        "--access",
        choices=("private", "public"),
        default="private",
        help="Project access mode when uploading (default: private)",
    )
    parser.add_argument(
        "--skip-upload",
        action="store_true",
        help="Only write the CSV and .microreact JSON locally",
    )
    parser.add_argument(
        "--ignore-case",
        action="store_true",
        help="Treat uppercase/lowercase alignment characters as equivalent when finding variant positions",
    )
    parser.add_argument(
        "--secondary-table",
        type=Path,
        help="Optional second CSV/TSV file to add as a linked metadata table",
    )
    parser.add_argument(
        "--secondary-master-field",
        default="id",
        help="Field in the primary dataset used to link to the secondary table (default: id)",
    )
    parser.add_argument(
        "--secondary-link-field",
        default="#Allele",
        help="Field in the secondary table used to link back to the primary dataset (default: #Allele)",
    )
    parser.add_argument(
        "--secondary-title",
        default="Metadata",
        help="Title for the optional secondary table",
    )
    parser.add_argument(
        "--secondary-delimiter",
        choices=("auto", "comma", "tab"),
        default="auto",
        help="Delimiter for the optional secondary table (default: auto)",
    )
    return parser.parse_args()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    current_name: str | None = None
    current_chunks: list[str] = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records.append((current_name, "".join(current_chunks)))
                current_name = line[1:].strip().split()[0]
                current_chunks = []
                continue
            if current_name is None:
                raise SystemExit(f"Invalid FASTA: sequence data found before a header in {path}")
            current_chunks.append(line)

    if current_name is not None:
        records.append((current_name, "".join(current_chunks)))

    if not records:
        raise SystemExit(f"No FASTA records found in {path}")

    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise SystemExit(
            f"Alignment contains sequences with differing lengths: {sorted(lengths)}"
        )

    seen_ids: set[str] = set()
    duplicate_ids = [name for name, _ in records if name in seen_ids or seen_ids.add(name)]
    if duplicate_ids:
        raise SystemExit(f"Duplicate FASTA IDs found: {', '.join(duplicate_ids[:5])}")

    return records


def normalise_base(value: str, ignore_case: bool) -> str:
    return value.upper() if ignore_case else value


def build_column_map(
    records: list[tuple[str, str]], reference_id: str | None
) -> list[tuple[int, int]]:
    alignment_length = len(records[0][1])
    if reference_id is None:
        return [(index, index + 1) for index in range(alignment_length)]

    record_map = dict(records)
    if reference_id not in record_map:
        raise SystemExit(f"Reference sequence {reference_id!r} was not found in the alignment")

    reference_sequence = record_map[reference_id]
    column_map: list[tuple[int, int]] = []
    reference_position = 0
    for alignment_index, value in enumerate(reference_sequence):
        if value in {"-", "."}:
            continue
        reference_position += 1
        column_map.append((alignment_index, reference_position))

    if not column_map:
        raise SystemExit(f"Reference sequence {reference_id!r} contains no non-gap positions")

    return column_map


def find_variant_positions(
    records: list[tuple[str, str]], column_map: list[tuple[int, int]], ignore_case: bool
) -> list[tuple[int, int]]:
    sequences = [sequence for _, sequence in records]
    variant_positions: list[tuple[int, int]] = []

    for alignment_index, label_position in column_map:
        column_values = [sequence[alignment_index] for sequence in sequences]
        observed = {
            normalise_base(value, ignore_case)
            for value in column_values
            if value not in {"?", "."}
        }
        if len(observed) > 1:
            variant_positions.append((alignment_index, label_position))

    if not variant_positions:
        raise SystemExit("No variant positions were found in the supplied alignment")

    return variant_positions


def build_rows(
    records: list[tuple[str, str]], variant_positions: list[tuple[int, int]], id_column: str
) -> tuple[list[str], list[dict[str, str]]]:
    header = [id_column] + [f"pos_{position}" for _, position in variant_positions]
    rows: list[dict[str, str]] = []

    for sample_id, sequence in records:
        row = {id_column: sample_id}
        for alignment_index, label_position in variant_positions:
            row[f"pos_{label_position}"] = sequence[alignment_index]
        rows.append(row)

    return header, rows


def write_csv(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, lineterminator="\r\n")
        writer.writeheader()
        writer.writerows(rows)


def detect_delimiter(path: Path, delimiter_mode: str) -> str:
    if delimiter_mode == "comma":
        return ","
    if delimiter_mode == "tab":
        return "\t"
    if path.suffix.lower() in {".tsv", ".tab"}:
        return "\t"
    return ","


def read_delimited_header(path: Path, delimiter_mode: str) -> list[str]:
    delimiter = detect_delimiter(path, delimiter_mode)
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise SystemExit(f"Secondary table is empty: {path}") from exc
    if not header:
        raise SystemExit(f"Secondary table has an empty header row: {path}")
    return header


def build_microreact_payload(
    project_name: str,
    csv_path: Path,
    csv_bytes: bytes,
    tree_path: Path,
    tree_bytes: bytes,
    id_column: str,
    table_fields: list[str],
    reference_id: str | None = None,
    secondary_table_path: Path | None = None,
    secondary_table_bytes: bytes | None = None,
    secondary_table_fields: list[str] | None = None,
    secondary_master_field: str | None = None,
    secondary_link_field: str | None = None,
    secondary_title: str = "Metadata",
) -> dict:
    csv_blob = "data:text/csv;base64," + base64.b64encode(csv_bytes).decode("ascii")
    tree_blob = (
        "data:application/octet-stream;base64,"
        + base64.b64encode(tree_bytes).decode("ascii")
    )
    position_title = "Variant positions"
    if reference_id:
        position_title = f"Variant positions ({reference_id} reference)"

    payload = {
        "schema": "https://microreact.org/schema/v1.json",
        "meta": {
            "name": project_name,
            "description": (
                f"Variant positions numbered relative to reference {reference_id}."
                if reference_id
                else "Variant positions numbered by alignment column."
            ),
        },
        "files": {
            "data-file-1": {
                "id": "data-file-1",
                "format": "text/csv",
                "name": csv_path.name,
                "type": "data",
                "size": len(csv_bytes),
                "blob": csv_blob,
            },
            "tree-file-1": {
                "id": "tree-file-1",
                "format": "text/x-nh",
                "name": tree_path.name,
                "type": "tree",
                "size": len(tree_bytes),
                "blob": tree_blob,
            },
        },
        "datasets": {
            "dataset-1": {
                "id": "dataset-1",
                "file": "data-file-1",
                "idFieldName": id_column,
            }
        },
        "tables": {
            "table-1": {
                "dataset": "dataset-1",
                "title": position_title,
                "columns": [{"field": field} for field in table_fields],
            }
        },
        "trees": {
            "tree-1": {
                "title": "Tree",
                "file": "tree-file-1",
            }
        },
    }

    if secondary_table_path is not None:
        assert secondary_table_bytes is not None
        assert secondary_table_fields is not None
        assert secondary_master_field is not None
        assert secondary_link_field is not None

        secondary_blob = (
            "data:text/csv;base64,"
            + base64.b64encode(secondary_table_bytes).decode("ascii")
        )
        payload["files"]["data-file-2"] = {
            "id": "data-file-2",
            "format": "text/csv",
            "name": secondary_table_path.name,
            "type": "data",
            "size": len(secondary_table_bytes),
            "blob": secondary_blob,
        }
        payload["datasets"]["dataset-2"] = {
            "id": "dataset-2",
            "file": "data-file-2",
            "masterFieldName": secondary_master_field,
            "linkFieldName": secondary_link_field,
        }
        payload["tables"]["table-2"] = {
            "displayMode": "cosy",
            "hideUnselected": False,
            "title": secondary_title,
            "paneId": "table-2",
            "file": "data-file-2",
            "columns": [{"field": field, "fixed": False} for field in secondary_table_fields],
        }

    return payload


def build_request(url: str, payload_bytes: bytes, token: str) -> urllib.request.Request:
    return urllib.request.Request(
        url,
        data=payload_bytes,
        headers={
            "Content-Type": "application/json; charset=utf-8",
            "Access-Token": token,
            "User-Agent": "curl/8.7.1",
            "Accept": "*/*",
        },
        method="POST",
    )


def upload_project(payload: dict, token: str, access: str) -> dict:
    query = ""
    if access == "private":
        query = "?" + urllib.parse.urlencode({"access": "private"})

    payload_bytes = json.dumps(payload).encode("utf-8")
    request = build_request(API_URL + query, payload_bytes=payload_bytes, token=token)

    try:
        with urllib.request.urlopen(request) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        if exc.code in {301, 302, 307, 308}:
            location = exc.headers.get("Location")
            if location:
                redirect_url = urllib.parse.urljoin(API_URL + query, location)
                redirect_request = build_request(
                    redirect_url, payload_bytes=payload_bytes, token=token
                )
                with urllib.request.urlopen(redirect_request) as response:
                    return json.loads(response.read().decode("utf-8"))
        error_body = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code}: {error_body}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(str(exc)) from exc


def get_token(token_env: str) -> tuple[str | None, str]:
    token = os.environ.get(token_env)
    if token:
        return token, token_env

    for fallback in TOKEN_ENV_FALLBACKS:
        if fallback == token_env:
            continue
        token = os.environ.get(fallback)
        if token:
            return token, fallback

    return None, token_env


def main() -> int:
    args = parse_args()
    if args.payload_output.suffix != ".microreact":
        args.payload_output = args.payload_output.with_suffix(".microreact")

    records = read_fasta(args.alignment)
    column_map = build_column_map(records, args.reference_id)
    variant_positions = find_variant_positions(
        records, column_map=column_map, ignore_case=args.ignore_case
    )
    header, rows = build_rows(records, variant_positions, args.id_column)

    write_csv(args.csv_output, header, rows)

    csv_bytes = args.csv_output.read_bytes()
    tree_bytes = args.tree.read_text().strip().encode("utf-8") + b"\n"
    secondary_table_fields = None
    secondary_table_bytes = None
    if args.secondary_table:
        if args.secondary_master_field not in header:
            raise SystemExit(
                f"Primary field {args.secondary_master_field!r} was not found in the variant CSV"
            )
        secondary_table_fields = read_delimited_header(
            args.secondary_table, args.secondary_delimiter
        )
        if args.secondary_link_field not in secondary_table_fields:
            raise SystemExit(
                f"Secondary field {args.secondary_link_field!r} was not found in {args.secondary_table}"
            )
        secondary_table_bytes = args.secondary_table.read_bytes()
    payload = build_microreact_payload(
        project_name=args.project_name,
        csv_path=args.csv_output,
        csv_bytes=csv_bytes,
        tree_path=args.tree,
        tree_bytes=tree_bytes,
        id_column=args.id_column,
        table_fields=header,
        reference_id=args.reference_id,
        secondary_table_path=args.secondary_table,
        secondary_table_bytes=secondary_table_bytes,
        secondary_table_fields=secondary_table_fields,
        secondary_master_field=args.secondary_master_field,
        secondary_link_field=args.secondary_link_field,
        secondary_title=args.secondary_title,
    )
    args.payload_output.write_text(json.dumps(payload, indent=2) + "\n")

    summary = {
        "project_name": args.project_name,
        "samples": len(records),
        "alignment_length": len(records[0][1]),
        "reference_id": args.reference_id,
        "variant_positions": len(variant_positions),
        "csv_output": str(args.csv_output.resolve()),
        "payload_output": str(args.payload_output.resolve()),
    }

    if args.skip_upload:
        print(json.dumps({"created": False, **summary}, indent=2))
        return 0

    token, token_source = get_token(args.token_env)
    if not token:
        print(
            (
                "Missing Microreact API token in environment variables "
                f"{', '.join(TOKEN_ENV_FALLBACKS)}. "
                "Use --skip-upload to only build local files."
            ),
            file=sys.stderr,
        )
        return 2

    try:
        response = upload_project(payload, token=token, access=args.access)
    except RuntimeError as exc:
        print(f"Microreact API request failed: {exc}", file=sys.stderr)
        return 1

    print(
        json.dumps(
            {"created": True, **summary, "token_env": token_source, "project": response},
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
