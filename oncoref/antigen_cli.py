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


def _host_overrides(args) -> dict:
    sex = getattr(args, "sex", None)
    if sex is not None and getattr(args, "profile", "all") != "all":
        raise ValueError("Use --sex or legacy --profile, not both")
    if sex is not None:
        from . import antigens

        host = antigens._build_host_context(
            species=getattr(args, "species", "human"),
            sex=sex,
            pregnant=bool(getattr(args, "pregnant", False)),
            lactating=bool(getattr(args, "lactating", False)),
            protect=tuple(getattr(args, "protect_tissues", None) or ()),
            accept_loss=tuple(getattr(args, "allow_tissues", None) or ()),
            absent=tuple(getattr(args, "absent_tissues", None) or ()),
        )
        return {"host": host} if host is not None else {}

    overrides = {}
    if getattr(args, "profile", "all") != "all":
        overrides["host_profile"] = args.profile
    if getattr(args, "lactating", False):
        overrides["lactating"] = True
    if getattr(args, "pregnant", False):
        overrides["pregnant"] = True
    if getattr(args, "allow_tissues", None):
        overrides["allow_tissues"] = tuple(args.allow_tissues)
    if getattr(args, "absent_tissues", None):
        overrides["absent_tissues"] = tuple(args.absent_tissues)
    if getattr(args, "protect_tissues", None):
        overrides["protect_tissues"] = tuple(args.protect_tissues)
    return overrides


def summary_antigens(args):
    from . import antigens

    emit(antigens.summary(species=getattr(args, "species", "all")), args.format)
    return 0


def list_antigens(args):
    from .species import resolve_species

    resolved = resolve_species(getattr(args, "species", "human"))
    overrides = {
        name: getattr(args, name)
        for name in ("min_tumor_fold", "min_reproductive_fold", "min_tumor_patients")
        if getattr(args, name) is not None
    }
    overrides.update(_host_overrides(args))
    ectopic_only = bool(getattr(args, "ectopic_only", False))
    replicated = bool(getattr(args, "replicated", False))
    if args.set in antigen_sets.EXCLUSION_CATEGORIES:
        if args.mode != "all" or overrides or ectopic_only or replicated:
            raise ValueError("Candidate filters do not apply to Unsafe/Other backgrounds")
        rows = antigen_sets.proteins(args.set, species=resolved)
    else:
        selection = antigen_sets.select(args.set, mode=args.mode, species=resolved, **overrides)
        rows = selection.proteins()
        if ectopic_only and "tumor_activation_mode" in rows.columns:
            rows = rows[
                rows["tumor_activation_mode"].isin(["ectopic", "ectopic_and_lineage_retained"])
            ].reset_index(drop=True)
        if replicated:
            if "tumor_support" in rows.columns:
                rows = rows[rows["tumor_support"].eq("replicated")].reset_index(drop=True)
            else:
                p_cnt = rows.get(
                    "independent_patients_at_tenfold", pd.Series(0, index=rows.index)
                ).fillna(0)
                s_cnt = rows.get("samples_at_tenfold", pd.Series(0, index=rows.index)).fillna(0)
                rows = rows[(p_cnt >= 2) | (s_cnt >= 2)].reset_index(drop=True)
        if args.format == "json":
            print(json.dumps(selection.metadata), file=sys.stderr)
    if args.ids:
        ids = sorted({gene for members in rows.member_gene_ids for gene in json.loads(members)})
        if args.format == "table":
            print("\n".join(ids))
        else:
            emit(pd.DataFrame({"gene_id": ids}), args.format)
    else:
        cols = [
            c
            for c in (
                "protein_group_label",
                "entity_id",
                "category",
                "sex_restriction_code",
                "tumor_activation_mode",
                "tumor_support",
                "reproductive_qualifying_compartments",
                "tumor_fold",
                "reproductive_fold",
                "independent_patients_at_tenfold",
                *(
                    ("host_status", "effective_off_target_tpm", "effective_tumor_fold")
                    if "host_status" in rows
                    else ()
                ),
            )
            if c in rows.columns
        ]
        emit(rows[cols], args.format)
    return 0


def explain(args):
    from .species import resolve_species
    from .species_target_selection import (
        packaged_species_host_exposure,
        packaged_species_results,
    )

    resolved = resolve_species(getattr(args, "species", "human"))
    rows = packaged_species_results(resolved, args.gene)
    if rows.empty:
        raise ValueError(f"No current antigen catalog record for {args.gene!r}")
    fields = [
        f
        for f in (
            "protein_group_label",
            "entity_id",
            "classification",
            "sex_restriction_code",
            "tumor_activation_mode",
            "tumor_support",
            "decision_reason",
            "normal_reference_tpm",
            "normal_reference_contexts",
            "disposable_max_tpm",
            "tumor_fold",
            "samples_at_tenfold",
            "independent_patients_at_tenfold",
            "cancer_types_at_tenfold",
            "reproductive_fold",
            "reproductive_qualifying_compartments",
            "reproductive_supported_tpm",
            "reproductive_supported_context",
            "reproductive_top_context",
            "cell_line_max_tpm",
            "qualifying_cell_lines",
            "rare_cell_status",
            "rare_cell_basis",
            "rare_cell_mapping_notes",
            "literature_recorded_tier",
            "paper_source_count",
            "paper_source_ids",
        )
        if f in rows
    ]
    output = rows[fields].copy()
    if "paper_source_count" in rows:
        counts = target_filters.paper_counts()
        for column in counts:
            output[column] = rows.entity_id.map(counts[column]).fillna(0).astype(int)
        if args.papers:
            from .cta_literature import paper_gene_claims

            output["paper_claims"] = [
                json.loads(paper_gene_claims(label).to_json(orient="records"))
                for label in rows.protein_group_label
            ]
        if args.loss:
            output["loss_consequences"] = [
                antigen_sets.loss_consequences(label) for label in rows.protein_group_label
            ]
    elif args.papers or args.loss:
        raise ValueError("--papers and --loss are currently human-only curation views")
    species_host_evidence = packaged_species_host_exposure(resolved)
    for mode in target_filters.FILTER_MODES:
        evaluated = target_filters.apply(
            target_filters.FILTER_MODES[mode], rows, host_evidence=species_host_evidence
        )
        output[mode] = [
            "pass" if row.keep else "fail: " + row.filter_failed for row in evaluated.itertuples()
        ]
    for label, filt in (
        (
            "profile:male-patient",
            target_filters.TargetFilter(membership="cta", host_profile="male-patient"),
        ),
        (
            "profile:female-patient",
            target_filters.TargetFilter(membership="cta", host_profile="female-patient"),
        ),
        (
            "profile:female-lactating",
            target_filters.TargetFilter(
                membership="cta", host_profile="female-patient", lactating=True
            ),
        ),
    ):
        evaluated = target_filters.apply(filt, rows, host_evidence=species_host_evidence)
        output[label] = [
            "pass" if row.keep else "fail: " + row.filter_failed for row in evaluated.itertuples()
        ]
    if overrides := _host_overrides(args):
        selected = target_filters.apply(
            replace(target_filters.TargetFilter(), **overrides),
            rows,
            host_evidence=species_host_evidence,
        )
        for column in (
            "host_status",
            "host_reasons",
            "host_unknown_compartments",
            "host_limiting_compartment",
            "effective_off_target_tpm",
            "effective_tumor_fold",
            "keep",
        ):
            output[column] = selected[column].to_numpy()
        if args.format == "json":
            print(json.dumps(selected.attrs), file=sys.stderr)
    if args.format == "table":
        for record in output.to_dict("records"):
            for name, value in record.items():
                print(f"{name}: {value}")
    else:
        emit(output, args.format)
    return 0


def funnel(args):
    from .species import resolve_species
    from .species_target_selection import packaged_species_results

    resolved = resolve_species(getattr(args, "species", "human"))
    table = packaged_species_results(resolved)
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
    filt = replace(target_filters.FILTER_MODES[args.mode], **_host_overrides(args))
    filtered = target_filters.funnel(filt, table).rename(columns={"layer": "stage"})
    filtered["kind"] = "sequential filter"
    extra_steps = []
    if getattr(args, "ectopic_only", False) or getattr(args, "replicated", False):
        cur = target_filters.apply(filt, table)
        cur = cur[cur["keep"]]
        if getattr(args, "ectopic_only", False) and "tumor_activation_mode" in cur.columns:
            nxt = cur[
                cur["tumor_activation_mode"].isin(["ectopic", "ectopic_and_lineage_retained"])
            ]
            extra_steps.append(
                {
                    "stage": "ectopic_only",
                    "remaining": len(nxt),
                    "removed": len(cur) - len(nxt),
                    "kind": "sequential filter",
                }
            )
            cur = nxt
        if getattr(args, "replicated", False) and "tumor_support" in cur.columns:
            nxt = cur[cur["tumor_support"].eq("replicated")]
            extra_steps.append(
                {
                    "stage": "replicated",
                    "remaining": len(nxt),
                    "removed": len(cur) - len(nxt),
                    "kind": "sequential filter",
                }
            )
    parts = [shared, categories, filtered]
    if extra_steps:
        parts.append(pd.DataFrame(extra_steps))
    emit(pd.concat(parts, ignore_index=True), args.format)
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
    host_overrides = _host_overrides(args)
    species = getattr(args, "species", "human")
    sequence_source = getattr(args, "sequence_source", "auto")
    if host_overrides:
        records = antigen_sets.select(
            "CA|RA",
            species=species,
            sequence_source=sequence_source,
            **host_overrides,
        ).background_sequences(args.background)
    elif species == "human" and sequence_source == "auto":
        records = antigen_sets.sequences(args.background)
    else:
        records = antigen_sets.sequences(
            args.background,
            species=species,
            sequence_source=sequence_source,
        )
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
                "selection": records.attrs,
                "match_policy": policy.as_dict(),
                "lengths": list(lengths),
                "interpretation": "Retained sequence windows, not validated antigen presentation or safety",
            }
        ),
        file=sys.stderr,
    )
    return 0


def _add_host_options(parser):
    parser.add_argument(
        "--sex",
        choices=["female", "male"],
        default=None,
        help="Patient biological sex: protects all adult reproductive & resting mammary organs of that sex",
    )
    parser.add_argument(
        "--profile",
        choices=list(target_filters.HOST_PROFILES),
        default="all",
        help="Legacy host profile: all (no host screen), male-patient, or female-patient",
    )
    parser.add_argument(
        "--lactating",
        action="store_true",
        help="Protect lactating mammary tissue (excludes milk-program antigens such as LALBA)",
    )
    parser.add_argument(
        "--pregnant",
        action="store_true",
        help="Protect embryo, placenta and milk-secretory state; developmental coverage remains incomplete",
    )
    parser.add_argument(
        "--accept-loss",
        "--allow-tissue",
        dest="allow_tissues",
        action="append",
        default=[],
        help="Explicitly accept whole-compartment loss; acute injury remains a separate concern",
    )
    parser.add_argument(
        "--absent-tissue",
        dest="absent_tissues",
        action="append",
        default=[],
        help="Declare a verified absent compartment (distinct from accepted injury)",
    )
    parser.add_argument(
        "--protect-tissue",
        dest="protect_tissues",
        action="append",
        default=[],
        help="Explicitly protect a compartment from toxicity (e.g. testis, ovary)",
    )


def register(subparsers):
    group = subparsers.add_parser(
        "antigens", help="Current antigen catalog, filters and peptide subtraction"
    )
    commands = group.add_subparsers(dest="antigen_command", required=True)
    overview = commands.add_parser(
        "summary", help="Cross-species antigen catalog summary (human, dog, mouse)"
    )
    overview.add_argument(
        "--species",
        default="all",
        help="Species name, taxonomy ID, or 'all' (default: all)",
    )
    overview.add_argument("--format", choices=["table", "csv", "json"], default="table")
    overview.set_defaults(func=summary_antigens)
    listing = commands.add_parser(
        "list", help="List current catalog members after optional filtering"
    )
    listing.add_argument("--set", choices=list(antigen_sets.ALL_SETS), default="CTA")
    listing.add_argument("--mode", choices=list(target_filters.FILTER_MODES), default="all")
    _add_host_options(listing)
    listing.add_argument(
        "--ectopic-only",
        action="store_true",
        help="Keep only antigens with confirmed ectopic tumor activation outside homologous source organs",
    )
    listing.add_argument(
        "--replicated",
        action="store_true",
        help="Keep only antigens supported by >=2 independent patients or tumor samples",
    )
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
        "--papers",
        action="store_true",
        help="Include each paper's exact gene claims and allowed uses",
    )
    _add_host_options(explanation)
    explanation.add_argument(
        "--loss", action="store_true", help="Include reviewed tissue-loss consequences and evidence"
    )
    explanation.set_defaults(func=explain)
    counts = commands.add_parser(
        "funnel", help="Catalog categories and sequential filter attrition"
    )
    counts.add_argument("--mode", choices=list(target_filters.FILTER_MODES), default="all")
    _add_host_options(counts)
    counts.add_argument("--ectopic-only", action="store_true")
    counts.add_argument("--replicated", action="store_true")
    counts.add_argument("--published-only", action="store_true")
    counts.set_defaults(func=funnel)
    peptides = commands.add_parser(
        "subtract", help="Subtract background matches from a protein sequence or FASTA"
    )
    peptides.add_argument("query")
    _add_host_options(peptides)
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
    peptides.add_argument(
        "--sequence-source",
        choices=list(antigen_sets.SEQUENCE_SOURCES),
        default="auto",
        help="Translation source: auto (prefer local sequence_bundle, else pyensembl), sequence_bundle, or pyensembl",
    )
    peptides.set_defaults(func=subtract)
    for command in (listing, explanation, counts, peptides):
        command.add_argument(
            "--species",
            default="human",
            help="Species name or taxonomy ID (human/9606, dog/canine/9615, mouse/10090)",
        )
        command.add_argument("--format", choices=["table", "csv", "json"], default="table")
