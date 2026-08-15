"""
Deduplicate compositional_clusters.GO-term.enrichments.txt

Produces compositional_clusters.GO-term.enrichments.cleaned.txt in the same
format as the original file, with three extra tab-separated columns appended
to each data line:
  col 8: n_parameter_sets_found  — how many parameter sets detected this protein
  col 9: all_parameter_sets      — semicolon-separated list (provenance)
  col 10: all_idcbrs_collapsed   — IDCBRs merged into this row (semicolon-separated)

Two-stage deduplication within each #CLUSTER_GOTERM_SET block:
  Stage 1 (Case 5): Same protein, same cluster, multiple IDCBRs
    → keep longest; tiebreak: IDCBR bias matching consensus bias
  Stage 2 (Cases 1+2): Same protein, multiple parameter sets in same block
    → keep lowest cov%; tiebreak: shortest len; then smallest cluster_id
    (NOT sorted by p-value)

Cross-block cases (3, 4, 6, 7) are preserved naturally since each block is
a unique (go_term, consensus_bias) pair.
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

INPUT_FILE = Path(__file__).parent.parent / "Raw Data" / "compositional_clusters.GO-term.enrichments.txt"
OUTPUT_FILE = Path(__file__).parent.parent / "compositional_clusters.GO-term.enrichments.cleanedV1.txt"

HEADER_RE = re.compile(r"^#CLUSTER_GOTERM_SET\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)")
GO_STATS_RE = re.compile(r"\(([^,]+),(\d+),(\d+),(\d+),([\d.e+\-]+),([ed])\)")
PARAM_RE = re.compile(r"len=(\d+)_cov=(\d+)%")
IDCBR_RE = re.compile(r"^[A-Z0-9]+_(\d+)_(\d+)_(.+)$")


def parse_idcbr(idcbr: str):
    m = IDCBR_RE.match(idcbr)
    if not m:
        return None, None, ""
    start, end, bias = int(m.group(1)), int(m.group(2)), m.group(3)
    return start, end, bias


def parse_param_set(param_set: str):
    m = PARAM_RE.match(param_set)
    if not m:
        return 9999, 9999
    return int(m.group(1)), int(m.group(2))


def parse_data_line(line: str, header_meta: dict):
    """Parse one data line. Returns a dict of fields, or None on failure."""
    parts = line.split()
    if len(parts) < 7:
        return None

    # Reconstruct the GO stats field (it may be split by spaces within parens)
    # Fields: acc bias cluster_id param_set idcbr go_stats has_go_term
    # go_stats is wrapped in () and may contain commas but not spaces — safe
    acc = parts[0]
    bias = parts[1]
    cluster_id = parts[2]
    param_set = parts[3]
    idcbr = parts[4]

    # GO stats: find the field starting with '('
    go_stats_str = ""
    has_go_term = ""
    for i, p in enumerate(parts[5:], 5):
        if p.startswith("("):
            go_stats_str = p
            if i + 1 < len(parts):
                has_go_term = parts[i + 1]
            break

    m = GO_STATS_RE.match(go_stats_str)
    if not m:
        return None

    go_term_count_orig = m.group(2)
    cluster_member_count = m.group(3)
    proteome_count = m.group(4)
    pvalue_str = m.group(5)
    enriched_depleted = m.group(6)

    try:
        pvalue = float(pvalue_str)
    except ValueError:
        pvalue = float("inf")

    len_val, cov_pct = parse_param_set(param_set)
    idcbr_start, idcbr_end, idcbr_bias = parse_idcbr(idcbr)
    idcbr_length = (idcbr_end - idcbr_start) if idcbr_start is not None else 0

    return {
        **header_meta,
        "uniprot_accession": acc,
        "protein_bias": bias,
        "cluster_id": cluster_id,
        "parameter_set": param_set,
        "idcbr_identifier": idcbr,
        "go_term_count_orig": go_term_count_orig,
        "cluster_member_count": cluster_member_count,
        "proteome_count": proteome_count,
        "pvalue": pvalue,
        "pvalue_str": pvalue_str,
        "enriched_depleted": enriched_depleted,
        "has_go_term": has_go_term,
        "_len_val": len_val,
        "_cov_pct": cov_pct,
        "_idcbr_length": idcbr_length,
        "_idcbr_bias": idcbr_bias,
    }


def stage1_within_cluster(entries: list, consensus_bias: str) -> dict:
    """
    Stage 1: For each (cluster_id, acc) group, collapse multiple IDCBRs to one.
    Returns dict keyed by (cluster_id, acc) → best entry dict with all_idcbrs added.
    """
    groups = defaultdict(list)
    for e in entries:
        groups[(e["cluster_id"], e["uniprot_accession"])].append(e)

    result = {}
    for key, group in groups.items():
        if len(group) == 1:
            best = group[0]
            best["all_idcbrs_this_cluster"] = best["idcbr_identifier"]
        else:
            # Pick longest; tiebreak: exact bias match with consensus_bias
            def sort_key(e):
                exact_match = 0 if e["_idcbr_bias"] == consensus_bias else 1
                return (-e["_idcbr_length"], exact_match)

            group_sorted = sorted(group, key=sort_key)
            best = group_sorted[0]
            best["all_idcbrs_this_cluster"] = ";".join(e["idcbr_identifier"] for e in group)
        result[key] = best
    return result


def stage2_cross_cluster(stage1_result: dict) -> list:
    """
    Stage 2: For each acc appearing in multiple clusters, keep lowest cov% entry.
    Returns list of final rows.
    """
    by_acc = defaultdict(list)
    for (cluster_id, acc), entry in stage1_result.items():
        by_acc[acc].append(entry)

    rows = []
    for acc, group in by_acc.items():
        if len(group) == 1:
            best = group[0]
            best["all_parameter_sets"] = best["parameter_set"]
            best["all_cluster_ids"] = best["cluster_id"]
            best["n_parameter_sets_found"] = 1
            best["single_param_set"] = True
        else:
            # Collect all parameter sets and cluster ids
            all_param_sets = ";".join(sorted(set(e["parameter_set"] for e in group)))
            all_cluster_ids = ";".join(sorted(set(e["cluster_id"] for e in group)))
            n_param = len(set(e["parameter_set"] for e in group))

            # Pick best: lowest cov_pct, then shortest len_val, then smallest cluster_id
            def sort_key(e):
                return (e["_cov_pct"], e["_len_val"], int(e["cluster_id"]) if e["cluster_id"].isdigit() else 999999)

            best = sorted(group, key=sort_key)[0]
            best["all_parameter_sets"] = all_param_sets
            best["all_cluster_ids"] = all_cluster_ids
            best["n_parameter_sets_found"] = n_param
            best["single_param_set"] = False

        rows.append(best)
    return rows


def process_block(header_meta: dict, raw_entries: list) -> list:
    entries = []
    for line in raw_entries:
        parsed = parse_data_line(line, header_meta)
        if parsed:
            entries.append(parsed)

    if not entries:
        return []

    stage1 = stage1_within_cluster(entries, header_meta["consensus_bias"])
    return stage2_cross_cluster(stage1)


FILE_HEADER = """\
# Cleaned version of compositional_clusters.GO-term.enrichments.txt
# Deduplication: one entry per protein per #CLUSTER_GOTERM_SET block.
#
# Data lines have 3 extra columns appended (cols 8-10):
#   col 8  n_parameter_sets_found  Number of parameter sets in which this protein was detected
#   col 9  all_parameter_sets      All parameter sets (semicolon-separated)
#   col 10 all_idcbrs_collapsed    All IDCBRs merged into this row (semicolon-separated)
#
# Selection rules:
#   Cross-parameter (same protein, multiple clusters in same block):
#     keep lowest cov% (highest confidence); tiebreak: shortest len; then smallest cluster_id
#   Within-cluster (same protein, multiple IDCBRs in same cluster):
#     keep longest IDCBR; tiebreak: bias matching cluster consensus bias
"""


def format_data_line(row: dict) -> str:
    """Reconstruct the original 7-field line and append 3 extra columns."""
    go_stats = (
        f"({row['go_term']},{row['go_term_count_orig']},{row['cluster_member_count']},"
        f"{row['proteome_count']},{row['pvalue_str']},{row['enriched_depleted']})"
    )
    fields = [
        row["uniprot_accession"],
        row["protein_bias"],
        row["cluster_id"],
        row["parameter_set"],
        row["idcbr_identifier"],
        go_stats,
        row["has_go_term"],
        str(row["n_parameter_sets_found"]),
        row["all_parameter_sets"],
        row["all_idcbrs_this_cluster"],
    ]
    return "\t".join(fields)


def main():
    print(f"Reading {INPUT_FILE}", file=sys.stderr)

    current_header_line = None
    current_meta = None
    current_raw = []
    total_input_rows = 0
    total_output_rows = 0
    blocks_processed = 0

    with (
        open(INPUT_FILE, encoding="utf-8") as fin,
        open(OUTPUT_FILE, "w", encoding="utf-8") as fout,
    ):
        fout.write(FILE_HEADER)

        def flush_block():
            nonlocal total_input_rows, total_output_rows, blocks_processed
            if current_meta is None or not current_raw:
                return
            total_input_rows += len(current_raw)
            rows = process_block(current_meta, current_raw)
            # Sort by cluster_id then accession to keep cluster groups together
            rows.sort(key=lambda r: (r["cluster_id"], r["uniprot_accession"]))
            fout.write("\n\n")
            fout.write(current_header_line + "\n")
            for row in rows:
                fout.write(format_data_line(row) + "\n")
            total_output_rows += len(rows)
            blocks_processed += 1
            if blocks_processed % 100 == 0:
                print(
                    f"  Processed {blocks_processed} blocks, "
                    f"{total_input_rows} input rows → {total_output_rows} output rows",
                    file=sys.stderr,
                )

        for line in fin:
            line = line.rstrip("\n\r")
            m = HEADER_RE.match(line)
            if m:
                flush_block()
                current_header_line = line
                current_meta = {
                    "go_term": m.group(1),
                    "consensus_bias": m.group(2),
                    "n_clusters_in_set": m.group(3),
                    "n_parameter_sets_in_set": m.group(4),
                    "cluster_size_range": m.group(5),
                }
                current_raw = []
            elif current_meta and line.strip() and not line.startswith("#"):
                current_raw.append(line)

        flush_block()

    print(
        f"\nDone. {total_input_rows} input rows → {total_output_rows} output rows "
        f"across {blocks_processed} blocks.",
        file=sys.stderr,
    )
    print(f"Output: {OUTPUT_FILE}", file=sys.stderr)


if __name__ == "__main__":
    main()
