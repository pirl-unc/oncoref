"""Current antigen catalog commands; no evidence or curation rules live here."""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import pandas as pd

from . import antigen_sets, target_filters, target_selection


def emit(frame, output_format):
    if output_format == "json":
        print(frame.to_json(orient="records", indent=2))
    elif output_format == "csv":
        print(frame.to_csv(index=False), end="")
    else:
        print(frame.to_string(index=False))


def list_antigens(args):
    rows = antigen_sets.proteins(args.set)
    overrides = {
        name: getattr(args, name)
        for name in ("min_tumor_fold", "min_reproductive_fold", "min_tumor_patients")
        if getattr(args, name) is not None
    }
    if args.set in antigen_sets.EXCLUSION_CATEGORIES:
        if args.mode != "all" or overrides:
            raise ValueError("Candidate filters do not apply to Unsafe/Other backgrounds")
    else:
        evaluated = target_filters.apply(
            replace(target_filters.FILTER_MODES[args.mode], **overrides), rows
        )
        ca, ra = evaluated.filtered_cancer_antigen, evaluated.filtered_reproductive_antigen
        member = {"CTA": ca & ra, "CA": ca, "RA": ra, "CA|RA": ca | ra}[args.set]
        rows = evaluated[evaluated.keep & member].copy()
        rows["category"] = rows.filtered_category
    if args.ids:
        ids = sorted({gene for members in rows.member_gene_ids for gene in json.loads(members)})
        if args.format == "table":
            print("\n".join(ids))
        else:
            emit(pd.DataFrame({"gene_id": ids}), args.format)
    else:
        emit(
            rows[
                [
                    "protein_group_label",
                    "entity_id",
                    "category",
                    "tumor_fold",
                    "reproductive_fold",
                    "independent_patients_at_tenfold",
                ]
            ],
            args.format,
        )
    return 0


def explain(args):
    rows = target_selection.results(args.gene)
    if rows.empty:
        # Resolve exact member symbols as well as grouped display labels.
        from .gene_matching import gene_identity_catalog

        identities = gene_identity_catalog()
        ids = identities.loc[identities.canonical_symbol.eq(args.gene), "canonical_gene_id"]
        pieces = [target_selection.results(gene) for gene in ids]
        rows = pd.concat(pieces).drop_duplicates("entity_id") if pieces else rows
    if rows.empty:
        raise ValueError(f"No current antigen catalog record for {args.gene!r}")
    fields = [
        "protein_group_label",
        "entity_id",
        "classification",
        "decision_reason",
        "normal_reference_tpm",
        "normal_reference_contexts",
        "disposable_max_tpm",
        "tumor_fold",
        "samples_at_tenfold",
        "independent_patients_at_tenfold",
        "cancer_types_at_tenfold",
        "reproductive_fold",
        "reproductive_supported_tpm",
        "reproductive_supported_context",
        "reproductive_top_context",
        "rare_cell_status",
        "rare_cell_basis",
        "rare_cell_mapping_notes",
        "literature_recorded_tier",
    ]
    output = rows[fields].copy()
    if args.loss:
        output["loss_consequences"] = [
            antigen_sets.loss_consequences(label) for label in rows.protein_group_label
        ]
    for mode in target_filters.FILTER_MODES:
        evaluated = target_filters.apply(target_filters.FILTER_MODES[mode], rows)
        output[mode] = [
            "pass" if row.keep else "fail: " + row.filter_failed for row in evaluated.itertuples()
        ]
    if args.format == "table":
        for record in output.to_dict("records"):
            for name, value in record.items():
                print(f"{name}: {value}")
    else:
        emit(output, args.format)
    return 0


def funnel(args):
    table = target_selection.results()
    if args.published_only:
        table = table[table.published_nomination]
    categories = table.category.value_counts().rename_axis("stage").reset_index(name="remaining")
    categories["kind"] = "catalog category"
    alive = pd.Series(True, index=table.index)
    shared = [{"stage": "All assessed proteins", "remaining": len(table), "removed": 0}]
    gates = [("No accepted normal exclusion", ~table.category.eq("Unsafe"))]
    gates.extend((name, table[name]) for name in target_selection.SHARED_CHECKS)
    gates.append(("Qualifying CA or RA; home-tissue rule", table.is_target_candidate))
    for name, passed in gates:
        before = int(alive.sum())
        alive &= passed
        shared.append(
            {"stage": name, "remaining": int(alive.sum()), "removed": before - int(alive.sum())}
        )
    shared = pd.DataFrame(shared).assign(kind="sequential catalog")
    filtered = target_filters.funnel(args.mode, table).rename(columns={"layer": "stage"})
    filtered["kind"] = "sequential filter"
    emit(pd.concat([shared, categories, filtered], ignore_index=True), args.format)
    return 0


def query_sequences(value):
    try:
        is_file = Path(value).is_file()
    except OSError:
        is_file = False
    if not is_file:
        yield "query", value.upper()
        return
    label, parts, seen = None, [], set()
    with Path(value).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if label is not None:
                    yield label, "".join(parts).upper()
                label, parts = line[1:].strip(), []
                if not label or label in seen:
                    raise ValueError("FASTA records need unique, nonempty identifiers")
                seen.add(label)
            elif label is None:
                raise ValueError("Sequence files must be FASTA")
            else:
                parts.append(line)
    if label is None:
        raise ValueError("Empty FASTA")
    yield label, "".join(parts).upper()


def subtract(args):
    from .peptide_subtraction import MatchPolicy, PeptideIndex

    costs = {}
    for value in args.substitution_cost:
        try:
            query, background, cost = value.split(":")
            costs[(query, background)] = float(cost)
        except ValueError as exc:
            raise ValueError("Substitution costs use QUERY:BACKGROUND:COST, e.g. I:L:0.5") from exc
    policy = MatchPolicy(
        max_distance=args.max_distance,
        equivalence_groups=tuple(args.equiv),
        substitution_costs=costs,
        insertion_cost=args.insertion_cost,
        deletion_cost=args.deletion_cost,
    )
    if not 1 <= args.min_length <= args.max_length:
        raise ValueError("Require 1 <= min-length <= max-length")
    lengths = range(args.min_length, args.max_length + 1)
    queries = list(query_sequences(args.query))
    for label, sequence in queries:
        if not sequence or not set(sequence) <= set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*"):
            raise ValueError(f"Invalid protein sequence for {label}")
    records = antigen_sets.sequences(args.background)
    if records.attrs.get("missing_gene_ids"):
        raise ValueError(
            "Background has missing gene translations; subtraction cannot claim complete coverage"
        )
    index = PeptideIndex(records, policy=policy, lengths=lengths)
    parts = [
        index.subtract(sequence, lengths=lengths).assign(query_id=label)
        for label, sequence in queries
    ]
    output = pd.concat(parts, ignore_index=True)
    emit(output, args.format)
    print(
        json.dumps(
            {
                "background": args.background,
                "records": len(records),
                "match_policy": policy.as_dict(),
                "lengths": list(lengths),
                "interpretation": "Retained sequence windows, not validated antigen presentation or safety",
            }
        ),
        file=sys.stderr,
    )
    return 0


def register(subparsers):
    group = subparsers.add_parser(
        "antigens", help="Current antigen catalog, filters and peptide subtraction"
    )
    commands = group.add_subparsers(dest="antigen_command", required=True)
    listing = commands.add_parser(
        "list", help="List current catalog members after optional filtering"
    )
    listing.add_argument("--set", choices=list(antigen_sets.SETS), default="CTA")
    listing.add_argument("--mode", choices=list(target_filters.FILTER_MODES), default="all")
    listing.add_argument("--min-tumor-fold", type=float)
    listing.add_argument("--min-reproductive-fold", type=float)
    listing.add_argument("--min-patients", dest="min_tumor_patients", type=int)
    listing.add_argument("--ids", action="store_true")
    listing.set_defaults(func=list_antigens)
    explanation = commands.add_parser(
        "explain", help="Explain every gate and named filter for a gene"
    )
    explanation.add_argument("gene")
    explanation.add_argument(
        "--loss", action="store_true", help="Include reviewed tissue-loss consequences and evidence"
    )
    explanation.set_defaults(func=explain)
    counts = commands.add_parser(
        "funnel", help="Catalog categories and sequential filter attrition"
    )
    counts.add_argument("--mode", choices=list(target_filters.FILTER_MODES), default="all")
    counts.add_argument("--published-only", action="store_true")
    counts.set_defaults(func=funnel)
    peptides = commands.add_parser(
        "subtract", help="Subtract background matches from a protein sequence or FASTA"
    )
    peptides.add_argument("query")
    peptides.add_argument(
        "--background", choices=list(antigen_sets.EXCLUSION_CATEGORIES), default="Unsafe|Other"
    )
    peptides.add_argument("--max-distance", type=float, default=0)
    peptides.add_argument(
        "--equiv",
        action="append",
        default=[],
        help="Explicit equivalence group, e.g. ILV; repeat for disjoint groups",
    )
    peptides.add_argument(
        "--substitution-cost",
        action="append",
        default=[],
        help="Directed QUERY:BACKGROUND:COST; repeat as needed",
    )
    peptides.add_argument("--insertion-cost", type=float, default=1)
    peptides.add_argument("--deletion-cost", type=float, default=1)
    peptides.add_argument("--min-length", type=int, default=8)
    peptides.add_argument("--max-length", type=int, default=15)
    peptides.set_defaults(func=subtract)
    for command in (listing, explanation, counts, peptides):
        command.add_argument("--format", choices=["table", "csv", "json"], default="table")
