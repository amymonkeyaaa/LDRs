#!/usr/bin/env python3
"""Extract multi-amino-acid compositional bias modules from a cleaned cluster file."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


GO_STATS_RE = re.compile(
    r"\(([^,]+),(\d+),(\d+),(\d+),([\d.eE+\-]+),([ed])\)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract multi-amino-acid compositional bias modules from a cleaned "
            "#CLUSTER_GOTERM_SET file. Default output is a clean block-style "
            "text file that resembles the cleaned cluster input."
        )
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="corrected cluster file/compositional_clusters.GO-term.enrichments.cleanedV1.txt",
        help="Cleaned compositional cluster GO enrichment file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="core protein:bias:module analysis/multi_amino_acid_bias_modules.cleanedV1.txt",
        help="Output TXT path.",
    )
    parser.add_argument(
        "--include-single-aa",
        action="store_true",
        help="Also include single-amino-acid consensus biases. Default: only biases like EK, SN, EDK.",
    )
    return parser.parse_args()


def parse_go_stats(value: str) -> dict[str, str]:
    match = GO_STATS_RE.search(value)
    if not match:
        return {
            "go_term_count": "",
            "cluster_member_count": "",
            "proteome_count": "",
            "p_value": "",
            "enrichment_type": "",
        }

    go_term, go_count, member_count, proteome_count, p_value, enrich_type = match.groups()
    return {
        "go_term_count": go_count,
        "cluster_member_count": member_count,
        "proteome_count": proteome_count,
        "p_value": p_value,
        "enrichment_type": enrich_type,
    }


def is_multi_amino_acid_bias(consensus_bias: str) -> bool:
    return len(consensus_bias) > 1


def parse_cluster_modules(
    path: Path, include_single_aa: bool = False
) -> dict[tuple[str, str, str], dict]:
    modules: dict[tuple[str, str, str], dict] = {}
    current_header = None

    with path.open(newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")

            if fields[0] == "#CLUSTER_GOTERM_SET":
                current_header = {
                    "GO_term": fields[1],
                    "consensus_bias": fields[2],
                    "n_clusters_in_go_bias_block": fields[3],
                    "n_parameter_sets_in_go_bias_block": fields[4],
                    "cluster_size_range_in_go_bias_block": fields[5],
                    "include_block": include_single_aa or is_multi_amino_acid_bias(fields[2]),
                }
                continue

            if (
                line.startswith("#")
                or current_header is None
                or not current_header["include_block"]
            ):
                continue

            if len(fields) < 7:
                continue

            accession = fields[0]
            row_bias = fields[1]
            cluster_id = fields[2]
            parameter_set = fields[3]
            idcbr = fields[4]
            go_stats = fields[5]
            go_annotated = fields[6]
            n_parameter_sets_found = fields[7] if len(fields) > 7 else ""
            all_parameter_sets = fields[8].split(";") if len(fields) > 8 and fields[8] else []
            all_idcbrs = fields[9].split(";") if len(fields) > 9 and fields[9] else [idcbr]

            key = (
                current_header["GO_term"],
                current_header["consensus_bias"],
                cluster_id,
            )

            if key not in modules:
                modules[key] = {
                    **current_header,
                    "ClusterID": cluster_id,
                    "row_biases": set(),
                    "parameter_sets": set(),
                    "accessions": set(),
                    "go_annotated_accessions": set(),
                    "idcbrs": set(),
                    "representative_idcbrs": set(),
                    "n_parameter_sets_found_values": [],
                    "go_stats": go_stats,
                    "rows": [],
                    **parse_go_stats(go_stats),
                }

            module = modules[key]
            module["row_biases"].add(row_bias)
            module["parameter_sets"].add(parameter_set)
            module["parameter_sets"].update(all_parameter_sets)
            module["accessions"].add(accession)
            module["representative_idcbrs"].add(idcbr)
            module["idcbrs"].update(all_idcbrs)
            if go_annotated == "YES":
                module["go_annotated_accessions"].add(accession)
            if n_parameter_sets_found:
                module["n_parameter_sets_found_values"].append(int(n_parameter_sets_found))
            module["rows"].append(fields)

    return modules


def module_sort_key(item: tuple[tuple[str, str, str], dict]) -> tuple:
    (go_term, consensus_bias, cluster_id), _module = item
    cluster_sort = int(cluster_id) if cluster_id.isdigit() else cluster_id
    return go_term, consensus_bias, cluster_sort


def write_modules(modules: dict[tuple[str, str, str], dict], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    sorted_modules = sorted(modules.items(), key=module_sort_key)

    with output.open("w", newline="") as handle:
        handle.write("# Extracted multi-amino-acid compositional bias modules from cleanedV1\n")
        handle.write("# Module definition here: ConsensusBias contains more than one amino acid, e.g. EK, SN, EDK.\n")
        handle.write("# Format preserves the cleaned cluster file style, with one ##MODULE section per ClusterID.\n")
        handle.write("# Data rows are the original cleaned cluster rows for that module.\n")
        handle.write("#\n")
        handle.write(
            "# GO block columns:\t#CLUSTER_GOTERM_SET\tGO_term\tConsensusBias\t"
            "ModulesInThisOutput\tOriginalNumberOfClusters\t"
            "NumberOfParameterSets\tClusterSizeRange\n"
        )
        handle.write(
            "# Module columns:\t##MODULE\tClusterID\trow_biases\tn_proteins\t"
            "n_go_annotated_proteins\tn_parameter_sets\tparameter_sets\t"
            "p_value\tenrichment_type\tGO_stats\n"
        )
        handle.write(
            "# Data columns:\tUniprotAccession\tBias\tClusterID\tParameter_Set\t"
            "IDCBR_Identifier\tGO-term statistics\tYES/NO\t"
            "n_parameter_sets_found\tall_parameter_sets\tall_idcbrs_collapsed\n\n"
        )

        current_block = None
        modules_in_block = {}
        for key, _module in sorted_modules:
            block_key = key[:2]
            modules_in_block[block_key] = modules_in_block.get(block_key, 0) + 1

        for (go_term, consensus_bias, _cluster_id), module in sorted_modules:
            block_key = (go_term, consensus_bias)
            if block_key != current_block:
                if current_block is not None:
                    handle.write("\n")
                handle.write(
                    "\t".join(
                        [
                            "#CLUSTER_GOTERM_SET",
                            go_term,
                            consensus_bias,
                            str(modules_in_block[block_key]),
                            module["n_clusters_in_go_bias_block"],
                            module["n_parameter_sets_in_go_bias_block"],
                            module["cluster_size_range_in_go_bias_block"],
                        ]
                    )
                    + "\n"
                )
                current_block = block_key

            handle.write(
                "\t".join(
                    [
                        "##MODULE",
                        module["ClusterID"],
                        ";".join(sorted(module["row_biases"])),
                        str(len(module["accessions"])),
                        str(len(module["go_annotated_accessions"])),
                        str(len(module["parameter_sets"])),
                        ";".join(sorted(module["parameter_sets"])),
                        module["p_value"],
                        module["enrichment_type"],
                        module["go_stats"],
                    ]
                )
                + "\n"
            )

            for row in module["rows"]:
                handle.write("\t".join(row) + "\n")
            handle.write("\n")


def main() -> None:
    args = parse_args()
    modules = parse_cluster_modules(Path(args.input), include_single_aa=args.include_single_aa)
    write_modules(modules, Path(args.output))
    print(f"Modules written: {len(modules)}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
