#!/usr/bin/env python3
"""Summarize multi-amino-acid bias module patterns from the cleaned module file."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path


INPUT = Path("core protein:bias:module analysis/multi_amino_acid_bias_modules.cleanedV1.txt")
CLEAN_CLUSTER_INPUT = Path("corrected cluster file/compositional_clusters.GO-term.enrichments.cleanedV1.txt")
OUT_DIR = Path("core protein:bias:module analysis")
MULTIPLE_MODULES_OUT = OUT_DIR / "case1_multiple_modules_same_GO_terms.txt"
SAME_MODULE_OUT = OUT_DIR / "case2_same_module+clusterID_across_multiple_GO_terms.txt"
SAME_BIAS_ANY_CLUSTER_OUT = OUT_DIR / "case2.2_same_bias_across_multiple_GO_terms.txt"
COMMON_PROTEIN_SAME_BIAS_OUT = OUT_DIR / "case3_common_protein_same_bias_multiple_GO_terms.txt"
COMMON_PROTEIN_SAME_MODULE_OUT = OUT_DIR / "case4_common_protein_same_module_multiple_GO_terms.txt"
CORE_PROTEIN_GROUP_SAME_BIAS_OUT = OUT_DIR / "case5_core_protein_groups_same_bias_multiple_GO_terms.txt"
CORE_PROTEIN_GROUP_SAME_MODULE_OUT = OUT_DIR / "case6_core_protein_groups_same_module_multiple_GO_terms.txt"
GO_TERMS_REPEATED_BIAS_OUT = OUT_DIR / "case7_GO_terms_repeated_bias_class.txt"

GO_STATS_RE = re.compile(r"\(([^,]+),(\d+),(\d+),(\d+),([\d.eE+\-]+),([ed])\)")


def parse_module_file(path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    modules = []
    protein_rows = []
    current_go = None
    current_consensus_bias = None
    current_module = None

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if fields[0] == "#CLUSTER_GOTERM_SET":
                current_go = fields[1]
                current_consensus_bias = fields[2]
            elif fields[0] == "##MODULE":
                current_module = {
                    "GO_term": current_go,
                    "module_bias": current_consensus_bias,
                    "ClusterID": fields[1],
                    "row_biases": fields[2],
                    "n_proteins": fields[3],
                    "n_go_annotated_proteins": fields[4],
                    "n_parameter_sets": fields[5],
                    "parameter_sets": fields[6],
                    "p_value": fields[7],
                    "enrichment_type": fields[8],
                    "GO_stats": fields[9],
                }
                modules.append(current_module)
            elif not line.startswith("#") and current_module is not None:
                protein_rows.append(
                    {
                        **current_module,
                        "accession": fields[0],
                        "protein_bias": fields[1],
                        "parameter_set": fields[3],
                        "IDCBR": fields[4],
                        "protein_GO_stats": fields[5],
                        "GO_annotated": fields[6],
                        "n_parameter_sets_found": fields[7] if len(fields) > 7 else "",
                        "all_parameter_sets": fields[8] if len(fields) > 8 else "",
                        "all_idcbrs_collapsed": fields[9] if len(fields) > 9 else "",
                    }
                )

    return modules, protein_rows


def parse_modules(path: Path) -> list[dict[str, str]]:
    modules, _protein_rows = parse_module_file(path)
    return modules


def parse_p_value(go_stats: str) -> str:
    match = GO_STATS_RE.search(go_stats)
    return match.group(5) if match else ""


def parse_cleaned_cluster_rows(path: Path) -> list[dict[str, str]]:
    rows = []
    current_go = None
    current_consensus_bias = None

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if fields[0] == "#CLUSTER_GOTERM_SET":
                current_go = fields[1]
                current_consensus_bias = fields[2]
                continue

            if line.startswith("#") or current_go is None:
                continue

            rows.append(
                {
                    "GO_term": current_go,
                    "consensus_bias": current_consensus_bias,
                    "accession": fields[0],
                    "protein_bias": fields[1],
                    "ClusterID": fields[2],
                    "parameter_set": fields[3],
                    "IDCBR": fields[4],
                    "GO_stats": fields[5],
                    "GO_annotated": fields[6],
                    "n_parameter_sets_found": fields[7] if len(fields) > 7 else "",
                    "all_parameter_sets": fields[8] if len(fields) > 8 else "",
                    "all_idcbrs_collapsed": fields[9] if len(fields) > 9 else "",
                    "p_value": parse_p_value(fields[5]),
                }
            )

    return rows


def write_multiple_modules_same_go(modules: list[dict[str, str]], output: Path) -> int:
    by_go = defaultdict(list)
    for module in modules:
        by_go[module["GO_term"]].append(module)

    rows_written = 0
    with output.open("w") as handle:
        handle.write("# Case 1: GO terms with multiple multi-amino-acid bias modules\n")
        handle.write("# Source: multi_amino_acid_bias_modules.cleanedV1.txt\n")
        handle.write("# A module is counted as GO_term + module_bias + ClusterID.\n")
        handle.write("#\n")

        for go_term, go_modules in sorted(
            by_go.items(), key=lambda item: (-len(item[1]), item[0])
        ):
            if len(go_modules) <= 1:
                continue

            handle.write(f"#GO_TERM\t{go_term}\tn_modules={len(go_modules)}\n")
            handle.write(
                "module_bias\tClusterID\tn_parameter_sets\tparameter_sets\t"
                "n_proteins\tn_go_annotated_proteins\tp_value\n"
            )

            for module in sorted(
                go_modules,
                key=lambda item: (
                    item["module_bias"],
                    int(item["ClusterID"]) if item["ClusterID"].isdigit() else item["ClusterID"],
                ),
            ):
                handle.write(
                    "\t".join(
                        [
                            module["module_bias"],
                            module["ClusterID"],
                            module["n_parameter_sets"],
                            module["parameter_sets"],
                            module["n_proteins"],
                            module["n_go_annotated_proteins"],
                            module["p_value"],
                        ]
                    )
                    + "\n"
                )
                rows_written += 1
            handle.write("\n")

    return rows_written


def write_same_module_across_go(modules: list[dict[str, str]], output: Path) -> int:
    by_module = defaultdict(list)
    for module in modules:
        key = (module["module_bias"], module["ClusterID"])
        by_module[key].append(module)

    rows_written = 0
    with output.open("w") as handle:
        handle.write("# Case 2: same multi-amino-acid bias module across multiple GO terms\n")
        handle.write("# Source: multi_amino_acid_bias_modules.cleanedV1.txt\n")
        handle.write("# Same module is defined as module_bias + ClusterID.\n")
        handle.write("#\n")

        for (module_bias, cluster_id), module_rows in sorted(
            by_module.items(), key=lambda item: (-len(item[1]), item[0][0], item[0][1])
        ):
            go_terms = sorted({row["GO_term"] for row in module_rows})
            if len(go_terms) <= 1:
                continue

            parameter_sets = sorted(
                {param for row in module_rows for param in row["parameter_sets"].split(";") if param}
            )
            handle.write(
                "\t".join(
                    [
                        "##MODULE",
                        module_bias,
                        f"ClusterID={cluster_id}",
                        f"n_GO_terms={len(go_terms)}",
                        f"n_parameter_sets={len(parameter_sets)}",
                        "parameter_sets=" + ";".join(parameter_sets),
                    ]
                )
                + "\n"
            )
            handle.write(
                "GO_term\tmodule_bias\tClusterID\tn_parameter_sets\t"
                "parameter_sets\tn_proteins\tn_go_annotated_proteins\tp_value\n"
            )

            for row in sorted(module_rows, key=lambda item: item["GO_term"]):
                handle.write(
                    "\t".join(
                        [
                            row["GO_term"],
                            row["module_bias"],
                            row["ClusterID"],
                            row["n_parameter_sets"],
                            row["parameter_sets"],
                            row["n_proteins"],
                            row["n_go_annotated_proteins"],
                            row["p_value"],
                        ]
                    )
                    + "\n"
                )
                rows_written += 1
            handle.write("\n")

    return rows_written


def write_same_bias_across_go_any_cluster(modules: list[dict[str, str]], output: Path) -> int:
    by_bias = defaultdict(list)
    for module in modules:
        by_bias[module["module_bias"]].append(module)

    rows_written = 0
    with output.open("w") as handle:
        handle.write("# Case 2.2: same multi-amino-acid bias across multiple GO terms, regardless of ClusterID\n")
        handle.write("# Source: multi_amino_acid_bias_modules.cleanedV1.txt\n")
        handle.write("# Unlike Case 2, a module here is defined by module_bias alone, ignoring ClusterID.\n")
        handle.write("# ClusterID is still reported per row for reference.\n")
        handle.write("#\n")

        for module_bias, bias_rows in sorted(
            by_bias.items(),
            key=lambda item: (-len({row["GO_term"] for row in item[1]}), item[0]),
        ):
            go_terms = sorted({row["GO_term"] for row in bias_rows})
            if len(go_terms) <= 1:
                continue

            cluster_ids = sorted(
                {row["ClusterID"] for row in bias_rows},
                key=lambda x: int(x) if x.isdigit() else x,
            )
            parameter_sets = sorted(
                {param for row in bias_rows for param in row["parameter_sets"].split(";") if param}
            )
            handle.write(
                "\t".join(
                    [
                        "##MODULE_BIAS",
                        module_bias,
                        f"n_GO_terms={len(go_terms)}",
                        f"n_clusters={len(cluster_ids)}",
                        "ClusterIDs=" + ";".join(cluster_ids),
                        f"n_parameter_sets={len(parameter_sets)}",
                        "parameter_sets=" + ";".join(parameter_sets),
                    ]
                )
                + "\n"
            )
            handle.write(
                "GO_term\tmodule_bias\tClusterID\tn_parameter_sets\t"
                "parameter_sets\tn_proteins\tn_go_annotated_proteins\tp_value\n"
            )

            for row in sorted(
                bias_rows,
                key=lambda item: (
                    item["GO_term"],
                    int(item["ClusterID"]) if item["ClusterID"].isdigit() else item["ClusterID"],
                ),
            ):
                handle.write(
                    "\t".join(
                        [
                            row["GO_term"],
                            row["module_bias"],
                            row["ClusterID"],
                            row["n_parameter_sets"],
                            row["parameter_sets"],
                            row["n_proteins"],
                            row["n_go_annotated_proteins"],
                            row["p_value"],
                        ]
                    )
                    + "\n"
                )
                rows_written += 1
            handle.write("\n")

    return rows_written


def write_common_protein_same_bias(cluster_rows: list[dict[str, str]], output: Path) -> int:
    by_protein_bias = defaultdict(list)
    for row in cluster_rows:
        by_protein_bias[(row["accession"], row["consensus_bias"])].append(row)

    rows_written = 0
    with output.open("w") as handle:
        handle.write("# Case 3: common protein with the same compositional bias across multiple GO terms\n")
        handle.write("# Source: compositional_clusters.GO-term.enrichments.cleanedV1.txt\n")
        handle.write("# This includes all consensus biases, not only multi-amino-acid biases.\n")
        handle.write("# Group key: accession + ConsensusBias.\n")
        handle.write("#\n")

        grouped_by_protein = defaultdict(list)
        for (accession, consensus_bias), rows in by_protein_bias.items():
            go_terms = sorted({row["GO_term"] for row in rows})
            if len(go_terms) > 1:
                grouped_by_protein[accession].append((consensus_bias, rows))

        for accession, bias_groups in sorted(
            grouped_by_protein.items(),
            key=lambda item: (
                -len({row["GO_term"] for _bias, rows in item[1] for row in rows}),
                item[0],
            ),
        ):
            handle.write(f"##PROTEIN\t{accession}\tn_biases={len(bias_groups)}\n")

            for consensus_bias, rows in sorted(
                bias_groups,
                key=lambda item: (-len({row["GO_term"] for row in item[1]}), item[0]),
            ):
                go_terms = sorted({row["GO_term"] for row in rows})
                cluster_ids = sorted(
                    {row["ClusterID"] for row in rows},
                    key=lambda x: int(x) if x.isdigit() else x,
                )
                parameter_sets = sorted(
                    {
                        param
                        for row in rows
                        for param in row["all_parameter_sets"].split(";")
                        if param
                    }
                )
                handle.write(
                    "\t".join(
                        [
                            "###BIAS",
                            consensus_bias,
                            f"n_GO_terms={len(go_terms)}",
                            f"n_clusters={len(cluster_ids)}",
                            "ClusterIDs=" + ";".join(cluster_ids),
                            "parameter_sets=" + ";".join(parameter_sets),
                        ]
                    )
                    + "\n"
                )
                handle.write(
                    "GO_term\tConsensusBias\tClusterID\tprotein_bias\tIDCBR\tGO_annotated\t"
                    "n_parameter_sets_found\tall_parameter_sets\tp_value\n"
                )
                for row in sorted(rows, key=lambda item: (item["GO_term"], item["ClusterID"])):
                    handle.write(
                        "\t".join(
                            [
                                row["GO_term"],
                                row["consensus_bias"],
                                row["ClusterID"],
                                row["protein_bias"],
                                row["IDCBR"],
                                row["GO_annotated"],
                                row["n_parameter_sets_found"],
                                row["all_parameter_sets"],
                                row["p_value"],
                            ]
                        )
                        + "\n"
                    )
                    rows_written += 1
                handle.write("\n")

            handle.write("\n")

    return rows_written


def write_common_protein_same_multi_aa_bias_legacy(
    protein_rows: list[dict[str, str]], output: Path
) -> int:
    by_protein_bias = defaultdict(list)
    for row in protein_rows:
        by_protein_bias[(row["accession"], row["module_bias"])].append(row)

    rows_written = 0
    with output.open("w") as handle:
        handle.write("# Legacy: common protein with the same multi-amino-acid bias across multiple GO terms\n")
        handle.write("# Source: multi_amino_acid_bias_modules.cleanedV1.txt\n")
        handle.write("# Group key: accession + module_bias.\n")
        handle.write("# Superseded by case3_common_protein_same_bias_multiple_GO_terms.txt for the full-bias case.\n")
        handle.write("#\n")

        for (accession, module_bias), rows in sorted(
            by_protein_bias.items(),
            key=lambda item: (-len({row["GO_term"] for row in item[1]}), item[0][1], item[0][0]),
        ):
            go_terms = sorted({row["GO_term"] for row in rows})
            if len(go_terms) <= 1:
                continue

            cluster_ids = sorted({row["ClusterID"] for row in rows}, key=lambda x: int(x) if x.isdigit() else x)
            handle.write(
                "\t".join(
                    [
                        "##PROTEIN_BIAS",
                        accession,
                        module_bias,
                        f"n_GO_terms={len(go_terms)}",
                        f"n_clusters={len(cluster_ids)}",
                        "ClusterIDs=" + ";".join(cluster_ids),
                    ]
                )
                + "\n"
            )
            handle.write(
                "GO_term\tmodule_bias\tClusterID\tprotein_bias\tIDCBR\tGO_annotated\t"
                "n_parameter_sets_found\tall_parameter_sets\tp_value\n"
            )
            for row in sorted(rows, key=lambda item: (item["GO_term"], item["ClusterID"])):
                handle.write(
                    "\t".join(
                        [
                            row["GO_term"],
                            row["module_bias"],
                            row["ClusterID"],
                            row["protein_bias"],
                            row["IDCBR"],
                            row["GO_annotated"],
                            row["n_parameter_sets_found"],
                            row["all_parameter_sets"],
                            row["p_value"],
                        ]
                    )
                    + "\n"
                )
                rows_written += 1
            handle.write("\n")

    return rows_written


def write_common_protein_same_module(
    protein_rows: list[dict[str, str]], output: Path
) -> int:
    by_protein_module_bias = defaultdict(list)
    for row in protein_rows:
        key = (row["accession"], row["module_bias"])
        by_protein_module_bias[key].append(row)

    rows_written = 0
    with output.open("w") as handle:
        handle.write("# Case 4: common protein with the same multi-amino-acid bias module across multiple GO terms\n")
        handle.write("# Source: multi_amino_acid_bias_modules.cleanedV1.txt\n")
        handle.write("# Group key: accession + module_bias. ClusterIDs are listed as supporting detail.\n")
        handle.write("#\n")

        for (accession, module_bias), rows in sorted(
            by_protein_module_bias.items(),
            key=lambda item: (-len({row["GO_term"] for row in item[1]}), item[0][1], item[0][0]),
        ):
            go_terms = sorted({row["GO_term"] for row in rows})
            if len(go_terms) <= 1:
                continue

            cluster_ids = sorted(
                {row["ClusterID"] for row in rows},
                key=lambda x: int(x) if x.isdigit() else x,
            )
            parameter_sets = sorted(
                {
                    param
                    for row in rows
                    for param in row["all_parameter_sets"].split(";")
                    if param
                }
            )
            handle.write(
                "\t".join(
                    [
                        "##PROTEIN_MULTI_AA_BIAS",
                        accession,
                        module_bias,
                        f"n_GO_terms={len(go_terms)}",
                        f"n_clusters={len(cluster_ids)}",
                        "ClusterIDs=" + ";".join(cluster_ids),
                        "parameter_sets=" + ";".join(parameter_sets),
                    ]
                )
                + "\n"
            )
            handle.write(
                "GO_term\tmodule_bias\tClusterID\tprotein_bias\tIDCBR\tGO_annotated\t"
                "n_parameter_sets_found\tall_parameter_sets\tp_value\n"
            )
            for row in sorted(rows, key=lambda item: (item["GO_term"], item["ClusterID"])):
                handle.write(
                    "\t".join(
                        [
                            row["GO_term"],
                            row["module_bias"],
                            row["ClusterID"],
                            row["protein_bias"],
                            row["IDCBR"],
                            row["GO_annotated"],
                            row["n_parameter_sets_found"],
                            row["all_parameter_sets"],
                            row["p_value"],
                        ]
                    )
                    + "\n"
                )
                rows_written += 1
            handle.write("\n")

    return rows_written


def write_core_protein_groups_same_bias(
    cluster_rows: list[dict[str, str]], output: Path, min_group_size: int = 2
) -> int:
    protein_bias_to_go = defaultdict(set)
    protein_bias_to_clusters = defaultdict(set)
    protein_bias_to_parameter_sets = defaultdict(set)
    for row in cluster_rows:
        key = (row["consensus_bias"], row["accession"])
        protein_bias_to_go[key].add(row["GO_term"])
        protein_bias_to_clusters[key].add(row["ClusterID"])
        for param in row["all_parameter_sets"].split(";"):
            if param:
                protein_bias_to_parameter_sets[key].add(param)

    by_signature = defaultdict(list)
    for (consensus_bias, accession), go_terms in protein_bias_to_go.items():
        if len(go_terms) <= 1:
            continue
        signature = (consensus_bias, tuple(sorted(go_terms)))
        by_signature[signature].append(accession)

    groups = [
        (signature, sorted(accessions))
        for signature, accessions in by_signature.items()
        if len(accessions) >= min_group_size
    ]
    groups.sort(key=lambda item: (-len(item[0][1]), -len(item[1]), item[0][0], item[1][0]))

    rows_written = 0
    with output.open("w") as handle:
        handle.write("# Case 5: core group of proteins with the same bias across multiple GO terms\n")
        handle.write("# Source: compositional_clusters.GO-term.enrichments.cleanedV1.txt\n")
        handle.write("# Core group definition: same ConsensusBias and exact same multi-GO signature, with at least 2 proteins.\n")
        handle.write("# This includes all consensus biases, not only multi-amino-acid biases.\n")
        handle.write("#\n")

        for (consensus_bias, go_terms), accessions in groups:
            cluster_ids = sorted(
                {
                    cluster_id
                    for accession in accessions
                    for cluster_id in protein_bias_to_clusters[(consensus_bias, accession)]
                },
                key=lambda x: int(x) if x.isdigit() else x,
            )
            parameter_sets = sorted(
                {
                    parameter_set
                    for accession in accessions
                    for parameter_set in protein_bias_to_parameter_sets[(consensus_bias, accession)]
                }
            )
            handle.write(
                "\t".join(
                    [
                        "#BIAS_GO_TERMS",
                        f"bias={consensus_bias}",
                        "GO_terms=" + ";".join(go_terms),
                        f"n_proteins={len(accessions)}",
                        f"n_GO_terms={len(go_terms)}",
                        f"n_clusters={len(cluster_ids)}",
                        "ClusterIDs=" + ";".join(cluster_ids),
                        "parameter_sets=" + ";".join(parameter_sets),
                    ]
                )
                + "\n"
            )
            handle.write("proteins\t" + ";".join(accessions) + "\n\n")
            rows_written += 1

    return rows_written


def write_core_protein_groups_same_module(
    protein_rows: list[dict[str, str]], output: Path, min_group_size: int = 2
) -> int:
    module_to_go = defaultdict(set)
    module_to_clusters = defaultdict(set)
    module_to_parameter_sets = defaultdict(set)
    for row in protein_rows:
        key = (row["module_bias"], row["accession"])
        module_to_go[key].add(row["GO_term"])
        module_to_clusters[key].add(row["ClusterID"])
        for param in row["all_parameter_sets"].split(";"):
            if param:
                module_to_parameter_sets[key].add(param)

    by_signature = defaultdict(list)
    for (module_bias, accession), go_terms in module_to_go.items():
        if len(go_terms) <= 1:
            continue
        signature = (module_bias, tuple(sorted(go_terms)))
        by_signature[signature].append(accession)

    groups = [
        (signature, sorted(accessions))
        for signature, accessions in by_signature.items()
        if len(accessions) >= min_group_size
    ]
    groups.sort(key=lambda item: (-len(item[0][1]), -len(item[1]), item[0][0], item[1][0]))

    rows_written = 0
    with output.open("w") as handle:
        handle.write("# Case 6: core group of proteins with the same multi-amino-acid bias module across multiple GO terms\n")
        handle.write("# Source: multi_amino_acid_bias_modules.cleanedV1.txt\n")
        handle.write("# Core group definition: same module_bias (the multi-amino-acid compositional bias itself, e.g. ED, HS) and exact same multi-GO signature, with at least 2 proteins.\n")
        handle.write("# A module here is identified by module_bias alone, independent of ClusterID; ClusterIDs are listed as supporting detail (same convention as Case 4).\n")
        handle.write("#\n")

        for (module_bias, go_terms), accessions in groups:
            cluster_ids = sorted(
                {
                    cluster_id
                    for accession in accessions
                    for cluster_id in module_to_clusters[(module_bias, accession)]
                },
                key=lambda x: int(x) if x.isdigit() else x,
            )
            parameter_sets = sorted(
                {
                    parameter_set
                    for accession in accessions
                    for parameter_set in module_to_parameter_sets[(module_bias, accession)]
                }
            )
            handle.write(
                "\t".join(
                    [
                        "#MODULE_GO_TERMS",
                        f"module_bias={module_bias}",
                        "GO_terms=" + ";".join(go_terms),
                        f"n_proteins={len(accessions)}",
                        f"n_GO_terms={len(go_terms)}",
                        f"n_clusters={len(cluster_ids)}",
                        "ClusterIDs=" + ";".join(cluster_ids),
                        "parameter_sets=" + ";".join(parameter_sets),
                    ]
                )
                + "\n"
            )
            handle.write("proteins\t" + ";".join(accessions) + "\n\n")
            rows_written += 1

    return rows_written


def write_go_terms_repeated_bias_class(
    cluster_rows: list[dict[str, str]], output: Path, min_repeats: int = 2
) -> int:
    cluster_agg: dict[tuple[str, str, str], dict[str, object]] = defaultdict(
        lambda: {"accessions": set(), "parameter_sets": set(), "p_value": ""}
    )
    for row in cluster_rows:
        key = (row["GO_term"], row["consensus_bias"], row["ClusterID"])
        agg = cluster_agg[key]
        agg["accessions"].add(row["accession"])
        for param in row["all_parameter_sets"].split(";"):
            if param:
                agg["parameter_sets"].add(param)
        if not agg["p_value"] and row["p_value"]:
            agg["p_value"] = row["p_value"]

    by_go_bias = defaultdict(set)
    for go_term, consensus_bias, cluster_id in cluster_agg:
        by_go_bias[(go_term, consensus_bias)].add(cluster_id)

    by_go = defaultdict(list)
    for (go_term, consensus_bias), cluster_ids in by_go_bias.items():
        if len(cluster_ids) >= min_repeats:
            by_go[go_term].append(consensus_bias)

    rows_written = 0
    with output.open("w") as handle:
        handle.write("# Case 7: GO terms that repeatedly use the same bias class across multiple clusters\n")
        handle.write("# Source: compositional_clusters.GO-term.enrichments.cleanedV1.txt\n")
        handle.write("# Includes all consensus biases, both single- and multi-amino-acid.\n")
        handle.write(
            "# A bias class is 'repeated' for a GO term when it appears in more than one distinct "
            "ClusterID within that GO term.\n"
        )
        handle.write("# Ranked by number of distinct repeated bias classes per GO term.\n")
        handle.write("#\n")

        for go_term, biases in sorted(
            by_go.items(), key=lambda item: (-len(item[1]), item[0])
        ):
            handle.write(f"#GO_TERM\t{go_term}\tn_repeated_bias_classes={len(biases)}\n")

            for consensus_bias in sorted(
                biases, key=lambda b: (-len(by_go_bias[(go_term, b)]), b)
            ):
                cluster_ids = sorted(
                    by_go_bias[(go_term, consensus_bias)],
                    key=lambda x: int(x) if x.isdigit() else x,
                )
                total_proteins = {
                    accession
                    for cluster_id in cluster_ids
                    for accession in cluster_agg[(go_term, consensus_bias, cluster_id)]["accessions"]
                }
                total_parameter_sets = sorted(
                    {
                        param
                        for cluster_id in cluster_ids
                        for param in cluster_agg[(go_term, consensus_bias, cluster_id)]["parameter_sets"]
                    }
                )
                handle.write(
                    "\t".join(
                        [
                            "##BIAS",
                            consensus_bias,
                            f"n_clusters={len(cluster_ids)}",
                            "ClusterIDs=" + ";".join(cluster_ids),
                            f"n_proteins={len(total_proteins)}",
                            "parameter_sets=" + ";".join(total_parameter_sets),
                        ]
                    )
                    + "\n"
                )
                handle.write("ClusterID\tn_proteins\tparameter_sets\tp_value\n")
                for cluster_id in cluster_ids:
                    agg = cluster_agg[(go_term, consensus_bias, cluster_id)]
                    handle.write(
                        "\t".join(
                            [
                                cluster_id,
                                str(len(agg["accessions"])),
                                ";".join(sorted(agg["parameter_sets"])),
                                agg["p_value"],
                            ]
                        )
                        + "\n"
                    )
                    rows_written += 1
                handle.write("\n")

            handle.write("\n")

    return rows_written


def main() -> None:
    modules, protein_rows = parse_module_file(INPUT)
    cluster_rows = parse_cleaned_cluster_rows(CLEAN_CLUSTER_INPUT)
    case1_rows = write_multiple_modules_same_go(modules, MULTIPLE_MODULES_OUT)
    case2_rows = write_same_module_across_go(modules, SAME_MODULE_OUT)
    case2_2_rows = write_same_bias_across_go_any_cluster(modules, SAME_BIAS_ANY_CLUSTER_OUT)
    case3_rows = write_common_protein_same_bias(cluster_rows, COMMON_PROTEIN_SAME_BIAS_OUT)
    case4_rows = write_common_protein_same_module(protein_rows, COMMON_PROTEIN_SAME_MODULE_OUT)
    case5_rows = write_core_protein_groups_same_bias(cluster_rows, CORE_PROTEIN_GROUP_SAME_BIAS_OUT)
    case6_rows = write_core_protein_groups_same_module(protein_rows, CORE_PROTEIN_GROUP_SAME_MODULE_OUT)
    case7_rows = write_go_terms_repeated_bias_class(cluster_rows, GO_TERMS_REPEATED_BIAS_OUT)
    print(f"Parsed modules: {len(modules)}")
    print(f"Parsed protein rows: {len(protein_rows)}")
    print(f"Parsed cleaned cluster rows: {len(cluster_rows)}")
    print(f"Case 1 rows written: {case1_rows}")
    print(f"Case 1 output: {MULTIPLE_MODULES_OUT}")
    print(f"Case 2 rows written: {case2_rows}")
    print(f"Case 2 output: {SAME_MODULE_OUT}")
    print(f"Case 2.2 rows written: {case2_2_rows}")
    print(f"Case 2.2 output: {SAME_BIAS_ANY_CLUSTER_OUT}")
    print(f"Case 3 rows written: {case3_rows}")
    print(f"Case 3 output: {COMMON_PROTEIN_SAME_BIAS_OUT}")
    print(f"Case 4 rows written: {case4_rows}")
    print(f"Case 4 output: {COMMON_PROTEIN_SAME_MODULE_OUT}")
    print(f"Case 5 groups written: {case5_rows}")
    print(f"Case 5 output: {CORE_PROTEIN_GROUP_SAME_BIAS_OUT}")
    print(f"Case 6 groups written: {case6_rows}")
    print(f"Case 6 output: {CORE_PROTEIN_GROUP_SAME_MODULE_OUT}")
    print(f"Case 7 rows written: {case7_rows}")
    print(f"Case 7 output: {GO_TERMS_REPEATED_BIAS_OUT}")


if __name__ == "__main__":
    main()
