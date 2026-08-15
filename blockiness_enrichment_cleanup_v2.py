from pathlib import Path
import pandas as pd

# Input raw blockiness enrichment file
input_file = Path("Raw Data/tertiles.blockiness.GO-term.enrichments.txt")

# Output cleaned file
output_file = Path("blockiness.same_sequence_dedup.cleanedv2.txt")

# Store parsed protein/ID-CBR rows here
rows = []

# Store the current GO-term header information
# Protein rows inherit this information until the next header appears
current_GO = None

# Read the raw text file line by line
with input_file.open("r") as f:
    for line in f:
        # Remove leading/trailing whitespace and newline characters
        line = line.strip()

        # Skip empty lines
        if not line:
            continue

        # Header lines contain the GO term, bias, tertile type, and tertile category
        if line.startswith("#TERTILE_GOTERM_SET"):
            parts = line.split("\t")

            current_GO = {
                "tertile_type": parts[1],          # blockiness
                "tertile": parts[2],               # H / I / L
                "bias": parts[3],                  # amino acid bias, e.g. A, QN, EK
                "parameter_set_count": parts[4],   # header-level parameter support, e.g. H0I3L0
                "GO_term": parts[5],               # GO term for this section
            }

        # Protein/ID-CBR rows inherit GO information from the most recent header
        else:
            parts = line.split("\t")

            # Skip malformed rows with too few columns
            if len(parts) < 8:
                continue

            rows.append({
                **current_GO,                       # copy GO_term, bias, tertile context into row
                "protein": parts[0],               # UniProt accession
                "bias_line": parts[1],             # row-level bias annotation
                "tertile_line": parts[2],          # row-level tertile annotation
                "parameter_set": parts[3],         # parameter set, e.g. len=20_cov=5%
                "IDCBR_id": parts[4],              # ID-CBR identifier
                "blockiness_std": float(parts[5]), # numeric blockiness value
                "sequence": parts[6],              # ID-CBR amino-acid sequence
                "GO_YES_NO": parts[7],             # YES/NO: whether protein is assigned this GO term
            })

# Convert parsed rows into a pandas dataframe
df = pd.DataFrame(rows)

print("Parsed data:")
print(df.head())
print("Raw parsed shape:", df.shape)

# Extract cov percentage from parameter_set.
# Example: len=20_cov=5% -> cov_percent = 5
# If cov is missing, value becomes NaN.
df["cov_percent"] = (
    df["parameter_set"]
    .str.extract(r"cov=(\d+)%")
    .astype(float)
)

# Count how many unique parameter sets support each exact same-sequence signal.
# This does NOT collapse different sequences from the same protein.
# Grouping key = same GO term + same bias + same tertile + same protein + exact same sequence.
df["n_parameter_sets"] = (
    df.groupby(["GO_term", "bias", "tertile", "protein", "sequence"])["parameter_set"]
    .transform("nunique")
)

# Mark rows where the row itself is cov=5.
df["is_cov5"] = df["cov_percent"] == 5

# Sort so that if exact duplicate sequences exist across parameter sets,
# the cov=5 row is kept first when possible.
df = df.sort_values(
    by=["GO_term", "bias", "tertile", "protein", "sequence", "is_cov5"],
    ascending=[True, True, True, True, True, False]
)

# Count unique exact-sequence signals before deduplication.
# This should equal the number after deduplication if we only remove exact duplicates.
n_unique_sequences_before = (
    df[["GO_term", "bias", "tertile", "protein", "sequence"]]
    .drop_duplicates()
    .shape[0]
)

# Collapse ONLY exact same-sequence duplicates caused by multiple parameter sets.
# This keeps:
# - all proteins
# - all different ID-CBR sequences
# - all different GO terms
# - all different bias/tertile categories
# It only removes repeated rows where the exact same sequence appears multiple times.
df_cleaned = df.drop_duplicates(
    subset=["GO_term", "bias", "tertile", "protein", "sequence"],
    keep="first"
).copy()

# Count unique exact-sequence signals after deduplication.
n_unique_sequences_after = (
    df_cleaned[["GO_term", "bias", "tertile", "protein", "sequence"]]
    .drop_duplicates()
    .shape[0]
)

# Save cleaned data in an original-like hierarchical text format.
# This keeps one header per GO-term/tertile/bias section, just like the raw file.
# Extra columns are appended to each protein row:
# n_parameter_sets
with output_file.open("w") as out:
    for header_values, group in df_cleaned.groupby(
        ["tertile_type", "tertile", "bias", "parameter_set_count", "GO_term"],
        sort=False
    ):
        tertile_type, tertile, bias, parameter_set_count, GO_term = header_values

        # Write the section header in the same style as the original file.
        out.write(
            f"#TERTILE_GOTERM_SET\t{tertile_type}\t{tertile}\t{bias}\t{parameter_set_count}\t{GO_term}\n"
        )

        # Write each cleaned protein/ID-CBR row underneath the header.
        for _, row in group.iterrows():
            out.write(
                f"{row['protein']}\t"
                f"{row['bias_line']}\t"
                f"{row['tertile_line']}\t"
                f"{row['parameter_set']}\t"
                f"{row['IDCBR_id']}\t"
                f"{row['blockiness_std']}\t"
                f"{row['sequence']}\t"
                f"{row['GO_YES_NO']}\t"
                f"{row['n_parameter_sets']}\n"
            )

        # Add a blank line between sections for readability.
        out.write("\n")

print("Saved original-format cleaned file:")
print(output_file)
print("Unique exact-sequence signals before cleaning:", n_unique_sequences_before)
print("Unique exact-sequence signals after cleaning:", n_unique_sequences_after)
print("These should be equal. If equal, only exact duplicate sequence rows were removed.")

