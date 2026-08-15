#!/usr/bin/env python3
"""
3-way merge: compositional clusters + blockiness + hpep tertiles.  (v2)

v2 extends v1 (merge_cluster_hpep_blockiness_v1.py) to recover matches for
cluster GO+Bias blocks whose bias is a multi-amino-acid combination
(e.g. "EK", "GN", "PQSMY"). Blockiness/hpep only ever report single-letter
biases (plus the one two-letter "EK"). Under v1's exact-bias-string match,
a cluster block like "GN" could never line up with any tertile block, even
though individual windows inside that cluster block often carry a genuine
single-letter signature (e.g. IDCBR ...(_G) that IS present in the
blockiness/hpep raw files under bias "G").

Extension (v2 only, on top of everything v1 already does):
  For a cluster row whose bias is multi-letter, if no exact-bias match is
  found in blockiness/hpep, retry the lookup using just the FIRST LETTER
  of the cluster bias as the tertile bias. GO term, accession, IDCBR, and
  parameter-set overlap must still all match — only the bias string used
  to find the tertile row is relaxed. A new column, tertile_bias_matched,
  records which bias value (exact cluster bias, or its first letter)
  actually produced the match for that row.

Sources:
  Raw files  — full parameter set lists per (go_term, bias, accession, IDCBR)
  Cleaned files — representative metric values after deduplication

Matching logic:
  1. GO_term + Bias    (block key; cluster blocks are kept if the exact bias,
                        OR — when the cluster bias is multi-letter — its
                        first letter, is present as a block in BOTH
                        blockiness and hpep)
  2. UniprotAccession  (protein must appear in all three for the same GO+Bias)
  3. IDCBR             (same protein region across all three; same IDCBR ->
                        same sequence. Row-level bias used for the
                        blockiness/hpep lookup is the exact cluster bias if
                        that finds a hit, else its first letter)
  4. Parameter overlap (cluster (intersect) blockiness (intersect) hpep must be non-empty)

Block header:
  #MERGED_GOTERM_SET  merged  bias  c_n_clusters  c_n_param_sets  c_size_range
                      b_tertile  b_counts  h_tertile  h_counts  n_matched  GO_term
  (b_tertile/b_counts/h_tertile/h_counts are reported per tertile-bias
  actually used by matched rows in this block, e.g. "EK=I;E=L" when rows in
  a single cluster block matched under two different tertile biases)

Data row columns (16):
  accession  bias  IDCBR  sequence  go_annotated
  cluster_id  cluster_GO_stats
  tertile_bias_matched
  b_tertile  b_deviation
  h_tertile  h_value
  matched_params  n_cluster_param_sets  n_blockiness_param_sets  n_hpep_param_sets
"""

RAW_C   = "/Users/monkeyaaa/Documents/GitHub/LDRs/Raw Data/compositional_clusters.GO-term.enrichments.txt"
RAW_B   = "/Users/monkeyaaa/Documents/GitHub/LDRs/Raw Data/tertiles.blockiness.GO-term.enrichments.txt"
RAW_H   = "/Users/monkeyaaa/Documents/GitHub/LDRs/Raw Data/tertiles.hpep.GO-term.enrichments.txt"
CLEAN_C = "/Users/monkeyaaa/Documents/GitHub/LDRs/corrected cluster file/compositional_clusters.GO-term.enrichments.cleanedV1.txt"
CLEAN_B = "/Users/monkeyaaa/Documents/GitHub/LDRs/corrected tertile file/blockiness.same_sequence_dedup.cleanedv2.txt"
CLEAN_H = "/Users/monkeyaaa/Documents/GitHub/LDRs/corrected tertile file/hpep.same_sequence_dedup.cleanedv2.txt"
OUT     = "/Users/monkeyaaa/Documents/GitHub/LDRs/merged files/cluster_hpep_blockiness_merged_v2.txt"


# ---------------------------------------------------------------------------
# Parsers (identical to v1)
# ---------------------------------------------------------------------------

def parse_raw_cluster(path):
    """
    Returns:
      data:       {(go_term, bias, acc, idcbr): set of param_sets}
      block_info: {(go_term, bias): {n_clusters, n_param_sets, size_range}}
    """
    data       = {}
    block_info = {}
    current    = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#CLUSTER"):
                p       = line.split("\t")
                go_term = p[1]
                bias    = p[2]
                current = (go_term, bias)
                block_info[current] = {
                    "n_clusters":   p[3],
                    "n_param_sets": p[4],
                    "size_range":   p[5],
                }
            elif not line.startswith("#"):
                p     = line.split("\t")
                acc   = p[0]
                ps    = p[3]   # param_set
                idcbr = p[4]   # IDCBR
                key   = (*current, acc, idcbr)
                data.setdefault(key, set()).add(ps)
    return data, block_info


def parse_raw_tertile_by_idcbr(path):
    """
    Returns:
      data: {(go_term, bias, acc, idcbr): set of param_sets}
      seqs: {(go_term, bias, acc, idcbr): sequence}   (first seen)
    """
    data    = {}
    seqs    = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                p       = line.split("\t")
                current = (p[5], p[3])   # (go_term, bias)
            else:
                p     = line.split("\t")
                acc   = p[0]
                ps    = p[3]
                idcbr = p[4]
                seq   = p[6]
                key   = (*current, acc, idcbr)
                data.setdefault(key, set()).add(ps)
                seqs.setdefault(key, seq)
    return data, seqs


def parse_clean_cluster(path):
    """
    Returns: {(go_term, bias, acc): {cluster_id, GO_stats, yes_no}}
    Keeps first row per protein per block (dedup already applied).
    """
    data    = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#CLUSTER"):
                p       = line.split("\t")
                current = (p[1], p[2])   # (go_term, bias)
            elif not line.startswith("#"):
                p   = line.split("\t")
                acc = p[0]
                key = (*current, acc)
                if key not in data:
                    data[key] = {
                        "cluster_id": p[2],
                        "GO_stats":   p[5],
                        "yes_no":     p[6],
                    }
    return data


def parse_clean_tertile_by_idcbr(path):
    """
    Returns:
      data:       {(go_term, bias, acc, idcbr): {tertile, metric, yes_no}}
      block_meta: {(go_term, bias): {tertile, block_counts}}
    """
    data       = {}
    block_meta = {}
    current    = None
    cur_meta   = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                p       = line.split("\t")
                current = (p[5], p[3])   # (go_term, bias)
                cur_meta = {"tertile": p[2], "block_counts": p[4]}
                block_meta.setdefault(current, cur_meta)
            else:
                p     = line.split("\t")
                acc   = p[0]
                idcbr = p[4]
                key   = (*current, acc, idcbr)
                if key not in data:
                    data[key] = {
                        "tertile": cur_meta["tertile"],
                        "metric":  p[5],
                        "yes_no":  p[7],
                    }
    return data, block_meta


# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------

def tertile_bias_candidates(bias):
    """Exact bias first; for a multi-amino-acid bias, also try its first
    letter as a fallback tertile bias."""
    candidates = [bias]
    if len(bias) > 1 and bias[0] not in candidates:
        candidates.append(bias[0])
    return candidates


def build_merged_blocks(raw_c, block_info_c,
                        raw_b, seqs_b,
                        raw_h,
                        clean_c, clean_b, block_meta_b,
                        clean_h, block_meta_h):

    c_blocks = set(k[:2] for k in raw_c)
    b_blocks = set(k[:2] for k in raw_b)
    h_blocks = set(k[:2] for k in raw_h)

    # A cluster GO+Bias block is kept if the exact bias, or (for multi-letter
    # biases) its first letter, forms a block present in BOTH blockiness and
    # hpep for that GO term.
    shared = []
    for (go_term, bias) in sorted(c_blocks):
        valid = [c for c in tertile_bias_candidates(bias)
                 if (go_term, c) in b_blocks and (go_term, c) in h_blocks]
        if valid:
            shared.append((go_term, bias, valid))

    output_blocks = []
    for (go_term, bias, valid_tbias) in shared:
        c_keys = {k for k in raw_c if k[:2] == (go_term, bias)}

        out_rows = []
        for key in sorted(c_keys):
            (gt, bi, acc, idcbr) = key
            c_params = raw_c[key]

            # Find which tertile bias (exact, or first-letter fallback)
            # actually has this same accession+IDCBR in both blockiness
            # and hpep.
            matched_tbias = None
            for tbias in valid_tbias:
                tkey = (go_term, tbias, acc, idcbr)
                if tkey in raw_b and tkey in raw_h:
                    matched_tbias = tbias
                    break
            if matched_tbias is None:
                continue

            tkey     = (go_term, matched_tbias, acc, idcbr)
            b_params = raw_b[tkey]
            h_params = raw_h[tkey]

            matched = sorted(c_params & b_params & h_params)
            if not matched:
                continue   # require at least one param shared across all three

            # Sequence from raw blockiness (guaranteed same as hpep since same IDCBR)
            seq = seqs_b.get(tkey, "")

            # Representative values from cleaned files
            cc = clean_c.get((go_term, bias, acc), {})
            bc = clean_b.get(tkey, {})
            hc = clean_h.get(tkey, {})

            go_annotated = "YES" if (
                cc.get("yes_no") == "YES" or
                bc.get("yes_no") == "YES" or
                hc.get("yes_no") == "YES"
            ) else "NO"

            out_rows.append([
                acc,
                bias,
                idcbr,
                seq,
                go_annotated,
                cc.get("cluster_id", ""),
                cc.get("GO_stats",   ""),
                matched_tbias,
                bc.get("tertile",    ""),
                bc.get("metric",     ""),
                hc.get("tertile",    ""),
                hc.get("metric",     ""),
                ";".join(matched),
                str(len(c_params)),
                str(len(b_params)),
                str(len(h_params)),
            ])

        if not out_rows:
            continue

        # Sort by cluster_id (numeric), then accession within the same cluster
        out_rows.sort(key=lambda r: (int(r[5]), r[0]))

        # tertile bias(es) actually used by matched rows in this block
        tbiases_used = sorted(set(r[7] for r in out_rows))

        bi_c = block_info_c.get((go_term, bias), {})

        def fmt_meta(block_meta, field):
            parts = []
            for tb in tbiases_used:
                meta = block_meta.get((go_term, tb), {})
                val  = meta.get(field, "")
                parts.append(val if len(tbiases_used) == 1 else f"{tb}={val}")
            return ";".join(parts)

        output_blocks.append({
            "go_term":      go_term,
            "bias":         bias,
            "c_n_clusters":   bi_c.get("n_clusters",   ""),
            "c_n_param_sets": bi_c.get("n_param_sets", ""),
            "c_size_range":   bi_c.get("size_range",   ""),
            "b_tertile":      fmt_meta(block_meta_b, "tertile"),
            "b_counts":       fmt_meta(block_meta_b, "block_counts"),
            "h_tertile":      fmt_meta(block_meta_h, "tertile"),
            "h_counts":       fmt_meta(block_meta_h, "block_counts"),
            "n_rows":  len(out_rows),
            "rows":    out_rows,
        })

    return output_blocks


def write_output(output_blocks, out_path):
    with open(out_path, "w") as f:
        for blk in output_blocks:
            header = "\t".join([
                "#MERGED_GOTERM_SET",
                "merged",
                blk["bias"],
                blk["c_n_clusters"],
                blk["c_n_param_sets"],
                blk["c_size_range"],
                blk["b_tertile"],
                blk["b_counts"],
                blk["h_tertile"],
                blk["h_counts"],
                f"n_matched={blk['n_rows']}",
                blk["go_term"],
            ])
            f.write(header + "\n")
            current_cluster = None
            for row in blk["rows"]:
                if row[5] != current_cluster:
                    current_cluster = row[5]
                    f.write(f"##CLUSTER_ID\t{current_cluster}\n")
                f.write("\t".join(row) + "\n")


def print_stats(output_blocks):
    total    = sum(b["n_rows"] for b in output_blocks)
    prots    = len(set(row[0] for b in output_blocks for row in b["rows"]))
    go_yes   = sum(1 for b in output_blocks for row in b["rows"] if row[4] == "YES")
    multi    = sum(1 for b in output_blocks if len(b["bias"]) > 1)
    fallback = sum(1 for b in output_blocks for row in b["rows"] if row[7] != row[1])
    print(f"GO+Bias blocks       : {len(output_blocks)}")
    print(f"  multi-letter bias  : {multi}")
    print(f"Matched rows         : {total}")
    print(f"  via first-letter   : {fallback}")
    print(f"Unique proteins      : {prots}")
    print(f"go_annotated=YES     : {go_yes}")


if __name__ == "__main__":
    print("Parsing raw cluster...")
    raw_c, block_info_c = parse_raw_cluster(RAW_C)
    print(f"  {len(set(k[:2] for k in raw_c))} GO+Bias blocks, {len(raw_c)} (go+bias+acc+idcbr) keys")

    print("Parsing raw blockiness...")
    raw_b, seqs_b = parse_raw_tertile_by_idcbr(RAW_B)
    print(f"  {len(set(k[:2] for k in raw_b))} GO+Bias blocks, {len(raw_b)} keys")

    print("Parsing raw hpep...")
    raw_h, _ = parse_raw_tertile_by_idcbr(RAW_H)
    print(f"  {len(set(k[:2] for k in raw_h))} GO+Bias blocks, {len(raw_h)} keys")

    print("Parsing cleaned cluster...")
    clean_c = parse_clean_cluster(CLEAN_C)

    print("Parsing cleaned blockiness...")
    clean_b, block_meta_b = parse_clean_tertile_by_idcbr(CLEAN_B)

    print("Parsing cleaned hpep...")
    clean_h, block_meta_h = parse_clean_tertile_by_idcbr(CLEAN_H)

    print("Building merged blocks...")
    output_blocks = build_merged_blocks(
        raw_c, block_info_c,
        raw_b, seqs_b,
        raw_h,
        clean_c, clean_b, block_meta_b,
        clean_h, block_meta_h,
    )

    print(f"Writing {OUT} ...")
    write_output(output_blocks, out_path=OUT)

    print_stats(output_blocks)
    print("Done.")
