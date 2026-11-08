#!/usr/bin/env python
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Derive the proteoform-group registry: sets of CGA genes encoding identical
proteins.

Some cancer-germline antigens are duplicated to two (or more) distinct genomic
loci that encode a *byte-identical* protein — e.g. CTAG1A/CTAG1B (NY-ESO-1),
XAGE1A/XAGE1B, SSX4/SSX4B, MAGEA2/MAGEA2B. RNA-seq reads multi-map between such
loci, so each gene's individual TPM under-counts the protein; the biologically
meaningful unit is the *proteoform* (the sum across member genes).

By default this scans the cancer-testis-antigen universe (plus an explicit seed of
known identical-protein partners that may not be in the CTA table yet), matches
the approved MANE Select or Ensembl canonical proteins from the pinned selection
registry, and writes one row per member locus for every shared sequence.

``--scope genome`` instead scans *all* protein-coding genes in the release,
producing the full identical-protein grouping (histone clusters, tubulins,
ubiquitins, … — any family with byte-identical loci whose RNA-seq reads multi-map
and individually under-count). The CTA-scoped registry is a refinement of it: each
CTA group's members fall within one genome group, which may merge in extra non-CTA
paralogs (so the genome label can be larger than the CTA label).
The genome registry is the opt-in variant (issue #12): genome-wide summation
shifts many more genes' expression than the focused CTA subset, so it is offered,
not shipped as the default for the percentile/within-sample artifacts.

The default generator uses packaged offline evidence. Only historical reproduction
with ``--legacy-longest`` requires the optional pyensembl dependency. Run:

    python scripts/generate_proteoform_groups.py               # CTA-scoped (default)
    python scripts/generate_proteoform_groups.py --scope genome  # genome-wide

Output: ``oncoref/data/proteoform-groups.csv`` (cta) or
``…/proteoform-groups-genome.csv`` (genome). The label is the slash-joined, sorted
member symbols (``SSX4/SSX4B``), matching the grouped-CTA labels used by the
downstream target-selection layer.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_CTA_CSV = _REPO_ROOT / "oncoref" / "data" / "cancer-testis-antigens.csv"
_DEFAULT_OUTPUT = _REPO_ROOT / "oncoref" / "data" / "proteoform-groups.csv"
_DEFAULT_GENOME_OUTPUT = _REPO_ROOT / "oncoref" / "data" / "proteoform-groups-genome.csv"

#: Known identical-protein partners that may not yet be curated into the CTA
#: table (e.g. the ``*B`` paralog whose ``*A``/base sibling is present). Seeding
#: their gene IDs lets the group form regardless of CTA-table freshness; the
#: scan still confirms protein identity, so a wrong seed simply drops out.
_SEED_GENE_IDS: tuple[str, ...] = (
    "ENSG00000269791",  # SSX4B  (partner of SSX4)
    "ENSG00000183305",  # MAGEA2B (partner of MAGEA2)
)

_OUTPUT_COLUMNS = (
    "proteoform_id",
    "member_symbol",
    "member_gene_id",
    "protein_length",
    "n_members",
)

#: Groups that MUST be derivable from any sane Ensembl release. If a release
#: renames or alters one of these gene IDs / sequences they would silently drop
#: out, leaving the registry on a basis that no longer matches the expression
#: data — so fail the build loudly instead.
_ANCHOR_GROUPS: frozenset[str] = frozenset(
    {"CTAG1A/CTAG1B", "XAGE1A/XAGE1B", "SSX4/SSX4B", "MAGEA2/MAGEA2B", "RBMY1F/RBMY1J"}
)


def _canonical_protein(gene) -> str | None:
    """Historical length-only choice; the name does not imply a canonical tag.

    This legacy registry builder does not check transcript biotype or CDS
    qualification. See audit_cta_transcripts.py before migrating its identities.
    """
    seqs = [t.protein_sequence for t in gene.transcripts if t.protein_sequence]
    if not seqs:
        return None
    return max(seqs, key=len)


def _candidate_gene_ids(cta_csv: Path) -> list[str]:
    ids: list[str] = []
    if cta_csv.exists():
        df = pd.read_csv(cta_csv, low_memory=False)
        if "Ensembl_Gene_ID" in df.columns:
            ids = df["Ensembl_Gene_ID"].dropna().astype(str).str.split(".").str[0].tolist()
    ids.extend(_SEED_GENE_IDS)
    # Stable de-dup, preserving first occurrence.
    return list(dict.fromkeys(ids))


def _genes_for_scope(data, gene_ids: list[str] | None):
    """Resolve the genes to scan: a specific id list (CTA scope) or every gene in
    the release (``gene_ids is None`` -> genome scope)."""
    if gene_ids is None:
        yield from data.genes()
        return
    for gene_id in gene_ids:
        try:
            yield data.gene_by_id(gene_id)
        except Exception:
            continue


def build_legacy_proteoform_groups(
    gene_ids: list[str] | None, ensembl_release: int, min_members: int
):
    import pyensembl

    data = pyensembl.EnsemblRelease(ensembl_release)
    by_protein: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for gene in _genes_for_scope(data, gene_ids):
        if gene.biotype != "protein_coding":
            continue
        protein = _canonical_protein(gene)
        if not protein:
            continue
        by_protein[protein].append((gene.gene_name, gene.gene_id))

    rows: list[dict[str, object]] = []
    for protein, members in by_protein.items():
        # Distinct gene IDs only (a gene can't group with itself).
        members = sorted(set(members))
        if len(members) < min_members:
            continue
        # Genome scope includes genes with no HGNC symbol; fall back to the gene id
        # for both the label and the member symbol so neither is empty/NaN (CTA
        # members all have symbols, so this is a no-op for the shipped default).
        members = [(symbol or gene_id, gene_id) for symbol, gene_id in members]
        label = "/".join(sorted(symbol for symbol, _ in members))
        for symbol, gene_id in members:
            rows.append(
                {
                    "proteoform_id": label,
                    "member_symbol": symbol,
                    "member_gene_id": gene_id,
                    "protein_length": len(protein),
                    "n_members": len(members),
                }
            )
    out = pd.DataFrame(rows, columns=list(_OUTPUT_COLUMNS))
    return out.sort_values(["proteoform_id", "member_symbol"]).reset_index(drop=True)


def build_proteoform_groups(gene_ids: list[str] | None, ensembl_release: int, min_members: int):
    """Group only qualified, approved default sequences from the shared registry."""
    from oncoref import transcript_selection, translation_qualification

    if ensembl_release != 115:
        raise ValueError("Current transcript selection is pinned to Ensembl 115")
    decisions = transcript_selection.records()
    table = translation_qualification.records()
    metadata = (
        table.loc[table.contig.isin(translation_qualification.PRIMARY_CONTIGS)]
        .drop_duplicates("gene_id")
        .set_index("gene_id")
    )
    rows = decisions.loc[decisions.selection_status.eq("selected")].copy()
    if gene_ids is not None:
        rows = rows.loc[rows.gene_id.isin(gene_ids)]
    rows = rows.loc[rows.gene_id.map(metadata.gene_biotype).eq("protein_coding")]
    result = []
    for sequence_hash, members in rows.groupby("selected_sequence_sha256", sort=True):
        if len(members) < min_members:
            continue
        symbols = {g: metadata.loc[g].symbol or g for g in members.gene_id}
        label = "/".join(sorted(symbols.values()))
        lengths = {n for values in members.candidate_lengths for n in values}
        if len(lengths) != 1:
            raise ValueError("Identical sequence hashes have inconsistent lengths")
        for gene, symbol in sorted(symbols.items()):
            result.append(
                {
                    "proteoform_id": label,
                    "member_symbol": symbol,
                    "member_gene_id": gene,
                    "protein_length": next(iter(lengths)),
                    "n_members": len(members),
                    "sequence_sha256": sequence_hash,
                    "selection_policy_id": transcript_selection.POLICY_ID,
                }
            )
    return (
        pd.DataFrame(result)
        .sort_values(["proteoform_id", "member_symbol", "member_gene_id"])
        .reset_index(drop=True)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ensembl-release", type=int, default=115)
    parser.add_argument(
        "--legacy-longest", action="store_true", help="Explicit historical reproduction only"
    )
    parser.add_argument("--cta-csv", type=Path, default=_DEFAULT_CTA_CSV)
    parser.add_argument(
        "--scope",
        choices=("cta", "genome"),
        default="cta",
        help="cta: scan the CTA universe (shipped default); genome: scan all "
        "protein-coding genes (opt-in, issue #12)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="defaults to proteoform-groups.csv (cta) or proteoform-groups-genome.csv (genome)",
    )
    parser.add_argument("--min-members", type=int, default=2)
    args = parser.parse_args()

    if args.output is None:
        args.output = _DEFAULT_OUTPUT if args.scope == "cta" else _DEFAULT_GENOME_OUTPUT
    gene_ids = None if args.scope == "genome" else _candidate_gene_ids(args.cta_csv)
    builder = build_legacy_proteoform_groups if args.legacy_longest else build_proteoform_groups
    groups = builder(gene_ids, args.ensembl_release, args.min_members)

    missing_anchors = _ANCHOR_GROUPS - set(groups["proteoform_id"])
    if missing_anchors:
        raise SystemExit(
            f"expected anchor proteoform groups missing: {sorted(missing_anchors)} — "
            f"Ensembl release {args.ensembl_release} may have changed these gene IDs or "
            f"sequences. Refusing to write a registry on a drifted basis."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    groups.to_csv(args.output, index=False)

    n_groups = groups["proteoform_id"].nunique()
    print(f"Wrote {len(groups)} member rows across {n_groups} proteoform groups to {args.output}")
    # CTA scope is small enough to enumerate; genome scope would print hundreds.
    if args.scope == "cta":
        for label, sub in groups.groupby("proteoform_id"):
            members = ", ".join(sub["member_symbol"].astype(str))
            print(f"  {label}: {members} ({sub['protein_length'].iloc[0]} aa)")


if __name__ == "__main__":
    main()
