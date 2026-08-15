#v3: same 9-way tertile pattern grouping as v1/v2, but sourced from the
#3-way merged file (cluster + blockiness + hpep) so each matched row carries
#a compositional ClusterID. Also computes, per cluster occurrence within a
#GO+Bias block, the YES:NO (go_annotated) ratio among its member rows.
from pathlib import Path

input_file = Path(
    "/Users/monkeyaaa/Documents/GitHub/LDRs/merged files/"
    "cluster_hpep_blockiness_merged_v1_same_sequence+parameter.txt"
)
output_file = Path(
    "/Users/monkeyaaa/Documents/GitHub/LDRs/key pattern tertile/"
    "pattern_combination_tertile_v3.txt"
)
ratio_summary_file = Path(
    "/Users/monkeyaaa/Documents/GitHub/LDRs/key pattern tertile/"
    "cluster_yes_no_ratio_v1.txt"
)

DATA_COLUMNS = [
    "accession", "bias", "IDCBR", "sequence", "go_annotated",
    "cluster_id", "cluster_GO_stats",
    "b_tertile", "b_deviation",
    "h_tertile", "h_value",
    "matched_params", "n_cluster_param_sets",
    "n_blockiness_param_sets", "n_hpep_param_sets",
]

#parse into: blocks = [ (header_fields, [ (cluster_id, [rows]) , ... ]) ]
blocks = []
current_header = None
current_clusters = []   #list of (cluster_id, rows) for current block
current_cluster_id = None
current_rows = []

def flush_cluster():
    global current_cluster_id, current_rows
    if current_cluster_id is not None:
        current_clusters.append((current_cluster_id, current_rows))
    current_cluster_id = None
    current_rows = []

def flush_block():
    global current_header, current_clusters
    flush_cluster()
    if current_header is not None:
        blocks.append((current_header, current_clusters))
    current_header = None
    current_clusters = []

with input_file.open("r") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith("#MERGED_GOTERM_SET"):
            flush_block()
            current_header = line.split("\t")
        elif line.startswith("##CLUSTER_ID"):
            flush_cluster()
            current_cluster_id = line.split("\t")[1]
        else:
            current_rows.append(line.split("\t"))

flush_block()

#header indices: 0 tag,1 merged,2 bias,3 c_n_clusters,4 c_n_param_sets,
#5 c_size_range,6 b_tertile,7 b_counts,8 h_tertile,9 h_counts,
#10 n_matched,11 GO_term

pattern_dict = {
    "HH": [], "HI": [], "HL": [],
    "IH": [], "II": [], "IL": [],
    "LH": [], "LI": [], "LL": [],
}

for header, clusters in blocks:
    b_tertile = header[6]
    h_tertile = header[8]
    combo = b_tertile + h_tertile
    pattern_dict[combo].append((header, clusters))

GO_ANNOTATED_COL = DATA_COLUMNS.index("go_annotated")

verify_total_rows = 0
verify_total_yes = 0
verify_unique_proteins = set()

with output_file.open("w") as f, ratio_summary_file.open("w") as rf:
    rf.write("pattern\tGO_term\tbias\tcluster_id\tn_yes\tn_no\tn_total\tyes_no_ratio\n")

    for pattern, block_list in pattern_dict.items():
        f.write(f"#Pattern: {pattern}\n")
        for header, clusters in block_list:
            bias = header[2]
            go_term = header[11]
            f.write("\t".join(header) + "\n")

            for cluster_id, rows in clusters:
                n_yes = sum(1 for r in rows if r[GO_ANNOTATED_COL] == "YES")
                n_no = sum(1 for r in rows if r[GO_ANNOTATED_COL] == "NO")
                n_total = len(rows)
                ratio = (n_yes / n_no) if n_no > 0 else (float("inf") if n_yes > 0 else 0.0)
                ratio_str = f"{ratio:.3f}" if ratio != float("inf") else "inf"

                verify_total_rows += n_total
                verify_total_yes += n_yes
                for r in rows:
                    verify_unique_proteins.add(r[0])

                f.write(f"##CLUSTER_ID\t{cluster_id}\tn_yes={n_yes}\tn_no={n_no}\tyes_no_ratio={ratio_str}\n")
                f.write("\t".join(DATA_COLUMNS) + "\n")
                for row in rows:
                    f.write("\t".join(row) + "\n")

                rf.write(f"{pattern}\t{go_term}\t{bias}\t{cluster_id}\t{n_yes}\t{n_no}\t{n_total}\t{ratio_str}\n")

            f.write("\n")
        f.write("\n")

print(f"Blocks: {sum(len(v) for v in pattern_dict.values())}")
print(f"Total matched rows: {verify_total_rows}")
print(f"Total YES rows: {verify_total_yes}")
print(f"Unique proteins (row-level accession set): {len(verify_unique_proteins)}")
