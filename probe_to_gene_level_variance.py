import pandas as pd
import numpy as np
import re

expr = pd.read_csv("/Users/yingjiawu/Desktop/AI_in_med_final_project/AI_in_med/annotated_expression_wk0_ust_only.csv")

annotation_cols = ["ID_REF", "Gene_symbol"]

sample_cols = [c for c in expr.columns if c.startswith("GSM")]
if len(sample_cols) == 0:
    sample_cols = [c for c in expr.columns if c not in annotation_cols]

print("Original shape:", expr.shape)
print("Number of samples:", len(sample_cols))

# Remove missing gene_symbol
expr = expr.dropna(subset=["Gene_symbol"]).copy()

expr["Gene_symbol"] = expr["Gene_symbol"].astype(str).str.strip()

# Remove empty strings
expr = expr[expr["Gene_symbol"] != ""]

# Remove ambiguous multi-gene mappings e.g. "GENE1 /// GENE2"
ambiguous_pattern = r"///|//|;|,"
expr = expr[~expr["Gene_symbol"].str.contains(ambiguous_pattern, regex=True)].copy()

# Fix gene symbols corrupted by Excel auto date conversion
def fix_excel_date_gene(symbol):
    m = re.match(r"^(\d{1,2})-Sep$", symbol)
    if m:
        return f"SEPT{m.group(1)}"
    m = re.match(r"^(\d{1,2})-Mar$", symbol)
    if m:
        return f"MARCH{m.group(1)}"
    m = re.match(r"^(\d{1,2})-Dec$", symbol)
    if m:
        return f"DEC{m.group(1)}"
    return symbol

date_pattern = r"^\d{1,2}-[A-Z][a-z]{2}$"
corrupted_mask = expr["Gene_symbol"].str.match(date_pattern, na=False)

expr["gene_symbol"] = expr["Gene_symbol"].apply(
    lambda x: fix_excel_date_gene(x) if re.match(date_pattern, str(x)) else x
)

print("After cleaning gene symbols:", expr.shape)

sample_cols = [c for c in sample_cols if c in expr.columns]
expr[sample_cols] = expr[sample_cols].apply(pd.to_numeric, errors="coerce")

# Remove probes where more than 20% of samples are missing
missing_rate = expr[sample_cols].isnull().mean(axis=1)
expr = expr[missing_rate < 0.2].copy()
print(f"After removing probes with >20% missing values: {expr.shape}")

# Collapse probes to gene-level. Keep the probe with the highest variance per gene
expr["variance"] = expr[sample_cols].var(axis=1)

expr = expr.sort_values(
    by=["gene_symbol", "variance"],
    ascending=[True, False]
)

gene_level = expr.drop_duplicates(subset="gene_symbol", keep="first").copy()
gene_level = gene_level[["gene_symbol"] + sample_cols].copy()

print("Gene-level matrix shape:", gene_level.shape)
print("Number of unique genes:", gene_level["gene_symbol"].nunique())

labels = pd.read_excel("/Users/yingjiawu/Desktop/AI_in_med_final_project/data/DEG_trial/response_label.xlsx")

print("\nLabels shape:", labels.shape)
print(labels.head())

labels = labels.rename(columns={"i_wk8_response": "response"})
labels = labels[["sample_id", "response"]].copy()
labels["sample_id"] = labels["sample_id"].astype(str).str.strip()

# Transpose: rows=genes, cols=samples  ->  rows=samples, cols=genes
X = gene_level.set_index("gene_symbol").T
X.index.name = "sample_id"
X.reset_index(inplace=True)
X["sample_id"] = X["sample_id"].astype(str).str.strip()

print("\nModel matrix shape before merge:", X.shape)

model_df = X.merge(labels, on="sample_id", how="inner")

lost = len(X) - len(model_df)

print(model_df["response"].value_counts())
response_rate = (model_df["response"] == "Y").mean()
print(f"Response rate: {response_rate:.2%}")

print("\nFinal modeling dataframe shape:", model_df.shape)
print(model_df[["sample_id", "response"]].head())

model_df.to_csv("/Users/yingjiawu/Desktop/model_input_gene_level_variance.csv", index=False)
gene_level.to_csv("/Users/yingjiawu/Desktop/gene_level_expression_bestprobe_variance.csv", index=False)