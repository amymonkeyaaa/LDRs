from pathlib import Path
import pandas as pd

# Input raw hpep enrichment file
input_file = Path("Raw Data/tertiles.hpep.GO-term.enrichments.txt")

# Output cleaned file
output_file = Path("hpep.same_sequence_dedup.cleanedv2.txt")

rows = []
current_GO = None

with input_file.open("r") as f:
    for line in f:
        line = line.strip()

        if not line:
            continue

        # Header line
        if line.startswith("#TERTILE_GOTERM_SET"):
            parts = line.split("\t")

            current_GO = {
                "tertile_type": parts[1],          # hpep
                "tertile": parts[2],               # H / I / L
                "bias": parts[3],                  # amino acid bias
                "parameter_set_count": parts[4],   # e.g. H3I0L0
                "GO_term": parts[5],               # GO term
            }

        # Protein / ID-CBR line
        else:
            parts = line.split("\t")

            if len(parts) < 8:
                continue

            rows.append({
                **current_GO,
                "protein": parts[0],
                "bias_line": parts[1],
                "tertile_line": parts[2],
                "parameter_set": parts[3],
                "IDCBR_id": parts[4],
                "hpep_value": float(parts[5]),
                "sequence": parts[6],
                "GO_YES_NO": parts[7],
            })

df = pd.DataFrame(rows)

print("Parsed hpep data:")
print(df.head())
print("Raw parsed shape:", df.shape)

# Extract cov percentage
df["cov_percent"] = (
    df["parameter_set"]
    .str.extract(r"cov=(\d+)%")
    .astype(float)
)

# Count how many unique parameter sets support each exact same sequence signal
df["n_parameter_sets"] = (
    df.groupby(["GO_term", "bias", "tertile", "protein", "sequence"])["parameter_set"]
    .transform("nunique")
)

# Mark cov=5 rows
df["is_cov5"] = df["cov_percent"] == 5

# Sort so cov=5 is preferentially kept when exact sequence duplicates exist
df = df.sort_values(
    by=["GO_term", "bias", "tertile", "protein", "sequence", "is_cov5"],
    ascending=[True, True, True, True, True, False]
)

# Count unique exact-sequence signals before deduplication
n_unique_sequences_before = (
    df[["GO_term", "bias", "tertile", "protein", "sequence"]]
    .drop_duplicates()
    .shape[0]
)

# Collapse ONLY exact same-sequence duplicates
df_cleaned = df.drop_duplicates(
    subset=["GO_term", "bias", "tertile", "protein", "sequence"],
    keep="first"
).copy()

# Count unique exact-sequence signals after deduplication
n_unique_sequences_after = (
    df_cleaned[["GO_term", "bias", "tertile", "protein", "sequence"]]
    .drop_duplicates()
    .shape[0]
)

# Save in original-like hierarchical format
with output_file.open("w") as out:
    for header_values, group in df_cleaned.groupby(
        ["tertile_type", "tertile", "bias", "parameter_set_count", "GO_term"],
        sort=False
    ):
        tertile_type, tertile, bias, parameter_set_count, GO_term = header_values

        out.write(
            f"#TERTILE_GOTERM_SET\t{tertile_type}\t{tertile}\t{bias}\t{parameter_set_count}\t{GO_term}\n"
        )

        for _, row in group.iterrows():
            out.write(
                f"{row['protein']}\t"
                f"{row['bias_line']}\t"
                f"{row['tertile_line']}\t"
                f"{row['parameter_set']}\t"
                f"{row['IDCBR_id']}\t"
                f"{row['hpep_value']}\t"
                f"{row['sequence']}\t"
                f"{row['GO_YES_NO']}\t"
                f"{row['n_parameter_sets']}\n"
            )

        out.write("\n")

print("Saved original-format cleaned hpep file:")
print(output_file)
print("Unique exact-sequence signals before cleaning:", n_unique_sequences_before)
print("Unique exact-sequence signals after cleaning:", n_unique_sequences_after)
print("These should be equal if only exact duplicate sequence rows were removed.")