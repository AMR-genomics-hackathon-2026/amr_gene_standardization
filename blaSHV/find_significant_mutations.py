from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
from statsmodels.tools import add_constant

# Load mutation data
mutations = pd.read_csv('variant-positions.csv')


# Load phenotype data (read header row, strip whitespace from columns)
genes = pd.read_csv('SHV_ncbirefgenes_noBla.tsv', sep='\t', header=0)
genes.columns = genes.columns.str.strip()




# Merge on 'id' (mutations) and '#Allele' (genes), keeping binary columns and BETA-LACTAM
binary_cols = ['CEFIDEROCOL/CEPHALOSPORIN', 'CEPHALOSPORIN', 'BETA-LACTAM']
cols_to_keep = ['#Allele'] + [col for col in binary_cols if col in genes.columns]
merged = pd.merge(mutations, genes[cols_to_keep], left_on='id', right_on='#Allele')

print('Merged dataframe preview:')
print(merged.head(20))
print(f"Merged dataframe shape: {merged.shape}")

# Generate wildtype pseudo-sequence from BETA-LACTAM=1 rows
mutation_cols = [col for col in mutations.columns if col != 'id']
wildtype_df = merged[merged['BETA-LACTAM'] == 1][mutation_cols]
wildtype_pseudo = wildtype_df.mode().iloc[0]  # most common amino acid per position
print('\nWildtype pseudo-sequence:')
print(wildtype_pseudo)


# Encode mutations as 1 (different from wildtype), 0 (same as wildtype) for all rows
mutation_matrix = merged[mutation_cols].apply(lambda col: col != wildtype_pseudo[col.name]).astype(int)

# Prepare X and y for logistic regression
X = mutation_matrix
y = merged['CEPHALOSPORIN']

print('\nMutation matrix preview:')
print(X.head())

# Remove columns with zero variance
X = X.loc[:, X.var() > 0]

print(f'Number of mutation columns after removing zero-variance: {X.shape[1]}')



# Scoary-like association analysis with sensitivity/specificity
print('\nRunning Scoary-like association analysis (Fisher exact test for each mutation) with sensitivity/specificity...')
mutation_names = X.columns
odds_ratios = []
p_values = []


sensitivity = []
specificity = []
table_sizes = []
n_00 = []  # mutation=0, phenotype=0
n_01 = []  # mutation=0, phenotype=1
n_10 = []  # mutation=1, phenotype=0
n_11 = []  # mutation=1, phenotype=1
for col in mutation_names:
    mut = X[col]
    pheno = y
    table = pd.crosstab(mut, pheno)
    # Sensitivity: fraction of R (y==1) samples with mutation
    sens = mut[pheno == 1].mean() if (pheno == 1).sum() > 0 else float('nan')
    # Specificity: fraction of S (y==0) samples without mutation
    spec = 1 - mut[pheno == 0].mean() if (pheno == 0).sum() > 0 else float('nan')
    sensitivity.append(sens)
    specificity.append(spec)
    # Number of samples in the contingency table (non-NA)
    table_sizes.append(table.values.sum())
    # Individual cell counts
    n_00.append(table.at[0,0] if (0 in table.index and 0 in table.columns) else 0)
    n_01.append(table.at[0,1] if (0 in table.index and 1 in table.columns) else 0)
    n_10.append(table.at[1,0] if (1 in table.index and 0 in table.columns) else 0)
    n_11.append(table.at[1,1] if (1 in table.index and 1 in table.columns) else 0)
    # Fisher's exact test
    if table.shape == (2,2):
        oddsratio, p = fisher_exact(table)
    else:
        oddsratio, p = float('nan'), float('nan')
    odds_ratios.append(oddsratio)
    p_values.append(p)

# Multiple testing correction (Benjamini-Hochberg FDR)
reject, pvals_corrected, _, _ = multipletests(p_values, method='fdr_bh')

# Collect results
results = pd.DataFrame({
    'mutation': mutation_names,
    'odds_ratio': odds_ratios,
    'p_value': p_values,
    'fdr_bh': pvals_corrected,
    'significant': reject,
    'sensitivity': sensitivity,
    'specificity': specificity,
    'n_samples_in_table': table_sizes,
    'n_00': n_00,
    'n_01': n_01,
    'n_10': n_10,
    'n_11': n_11
})
results = results.sort_values('p_value')
results.to_csv('mutation_scoary_results.csv', index=False)
print('Top significant mutations:')
print(results.head(10))
print('Full results saved to mutation_scoary_results.csv')

# Random Forest for combinatorial prediction
print('\nTraining Random Forest to find combinations of mutations predictive of resistance...')
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X, y)
y_pred = rf.predict(X)
print(classification_report(y, y_pred, target_names=['Sensitive', 'Resistant']))

# Feature importances
importances = pd.Series(rf.feature_importances_, index=mutation_names)
importances = importances.sort_values(ascending=False)
importances.to_csv('mutation_rf_importances.csv')
print('Top mutations by Random Forest importance:')
print(importances.head(10))
print('Full importances saved to mutation_rf_importances.csv')


# Save merged dataframe to TSV file
merged.to_csv('merged_data.tsv', sep='\t', index=False)
print('Merged dataframe saved to merged_data.tsv')

