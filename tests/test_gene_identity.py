import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import oncoref
from oncoref import gene_identity

PRAME = "ENSG00000185686"
ALTERNATE = "ENSG00000275013"


class Annotation:
    species = SimpleNamespace(latin_name="homo_sapiens")
    reference_name = "GRCh38"
    release = 93

    def __init__(self, genes):
        self.genes = {g.gene_id: g for g in genes}

    def gene_by_id(self, gid):
        try:
            return self.genes[gid]
        except KeyError as exc:
            raise ValueError(gid) from exc

    def gene_ids(self):
        return list(self.genes)


def gene(gid, name="PRAME", contig="22"):
    return SimpleNamespace(gene_id=gid, gene_name=name, contig=contig, start=1, end=100, strand="+")


def test_real_shipped_prame_alias_has_auditable_annotation_identity():
    genome = Annotation([gene(PRAME), gene(ALTERNATE, contig="CHR_HSCHR22_1_CTG3")])
    record = oncoref.resolve_gene_identity(f" {ALTERNATE}.2 ", genome=genome)
    assert record.verified and record.status == "alias"
    assert record.input_gene_id == f" {ALTERNATE}.2 "
    assert record.source_gene_id == ALTERNATE
    assert record.canonical_gene_id == PRAME
    assert record.mapping_method == "ensembl_alt_contig"
    assert record.annotation_release == 93
    assert record.canonical_reference_releases == ("115",)
    assert record.source_contig == "CHR_HSCHR22_1_CTG3"
    assert record.alias_rows == ((PRAME, "PRAME", "ensembl_alt_contig"),)
    assert len(record.aliases_sha256) == len(record.canonical_gene_space_sha256) == 64
    assert json.loads(json.dumps(record.as_dict()))["input_gene_id"] == record.input_gene_id
    assert oncoref.cta_annotation_gene_ids(genome) == {PRAME, ALTERNATE}
    assert oncoref.cta_annotation_gene_ids(genome, unfiltered=True) == {PRAME, ALTERNATE}


@pytest.mark.parametrize("name,contig", [("PRAMEF1", "ALT"), ("PRAME", "22")])
def test_conflicting_alias_not_admitted_or_excluded(name, contig):
    genome = Annotation([gene(PRAME), gene(ALTERNATE, name, contig)])
    record = gene_identity.resolve_gene_identity(ALTERNATE, genome=genome)
    assert not record.verified and record.status == "annotation_conflict"
    assert record.canonical_gene_id is None and record.candidate_gene_ids == (PRAME,)
    assert gene_identity.cta_annotation_gene_ids(genome, unfiltered=True) == {PRAME}
    assert len(gene_identity.cta_annotation_gene_identities(genome)) == 2


def test_same_symbol_and_unknown_locus_never_creates_alias():
    unknown = "ENSG99999999999"
    genome = Annotation([gene(PRAME), gene(unknown, contig="ALT")])
    assert gene_identity.resolve_gene_identity(unknown, genome=genome).status == "unmapped"
    assert gene_identity.cta_annotation_gene_ids(genome, unfiltered=True) == {PRAME}
    assert gene_identity.resolve_gene_identity(ALTERNATE, genome=genome).status == "absent"


def test_missing_primary_still_verifies_explicit_alternate():
    genome = Annotation([gene(ALTERNATE, contig="ALT")])
    assert gene_identity.cta_annotation_gene_ids(genome) == {ALTERNATE}


def test_canonical_locus_conflict_is_retained():
    genome = Annotation([gene(PRAME, contig="1")])
    assert gene_identity.resolve_gene_identity(PRAME, genome=genome).status == "annotation_conflict"
    assert gene_identity.cta_annotation_gene_ids(genome, unfiltered=True) == set()


@pytest.fixture
def mapping_tables(monkeypatch):
    actual = gene_identity.get_data
    aliases = actual("ensembl-id-aliases").copy()
    canonical = actual("canonical-gene-space").copy()

    def get_data(name, **kwargs):
        return {"ensembl-id-aliases": aliases, "canonical-gene-space": canonical}[name]

    monkeypatch.setattr(gene_identity, "get_data", get_data)
    gene_identity._reference.cache_clear()
    yield aliases, canonical
    gene_identity._reference.cache_clear()


def test_ambiguous_alias_targets_are_audited_without_exclusion(mapping_tables):
    aliases, _ = mapping_tables
    aliases.loc[len(aliases)] = [ALTERNATE, "ENSG00000141510", "TP53", "ensembl_alt_contig"]
    genome = Annotation([gene(ALTERNATE, contig="ALT")])
    record = gene_identity.resolve_gene_identity(ALTERNATE, genome=genome)
    assert record.status == "ambiguous" and not record.verified
    assert set(record.candidate_gene_ids) == {PRAME, "ENSG00000141510"}
    assert gene_identity.cta_annotation_gene_ids(genome, unfiltered=True) == set()
    assert len(gene_identity.cta_annotation_gene_identities(genome)) == 1


def test_nonunique_primary_symbol_does_not_validate_alias(mapping_tables):
    _, canonical = mapping_tables
    extra = canonical.loc[canonical.ensembl_gene_id.eq(PRAME)].iloc[0].copy()
    extra["ensembl_gene_id"] = "ENSG99999999999"
    canonical.loc[len(canonical)] = extra
    record = gene_identity.resolve_gene_identity(
        ALTERNATE, genome=Annotation([gene(ALTERNATE, contig="ALT")])
    )
    assert record.status == "annotation_conflict" and not record.verified


def test_historical_alias_requires_absent_successor_and_matching_chromosome():
    old, current = "ENSG00000005955", "ENSG00000278311"
    annotation = Annotation([gene(old, "GGNBP2", "17")])
    assert gene_identity.resolve_gene_identity(old, genome=annotation).status == "alias"
    annotation.genes[current] = gene(current, "GGNBP2", "17")
    assert (
        gene_identity.resolve_gene_identity(old, genome=annotation).status == "annotation_conflict"
    )


@pytest.mark.parametrize("attribute,value", [("reference_name", "GRCm39"), ("release", None)])
def test_annotation_identity_required(attribute, value):
    genome = Annotation([])
    setattr(genome, attribute, value)
    with pytest.raises(ValueError, match="requires an explicit"):
        gene_identity.resolve_gene_identity(PRAME, genome=genome)


def test_query_failures_propagate():
    genome = Annotation([])
    genome.gene_by_id = lambda gid: (_ for _ in ()).throw(RuntimeError("broken annotation"))
    with pytest.raises(RuntimeError, match="broken annotation"):
        gene_identity.resolve_gene_identity(PRAME, genome=genome)


def test_default_and_full_candidate_universe_stay_distinct(monkeypatch):
    from oncoref import cta

    monkeypatch.setattr(cta, "cta_gene_ids", lambda: set())
    genome = Annotation([gene(PRAME), gene(ALTERNATE, contig="ALT")])
    assert not gene_identity.cta_annotation_gene_ids(genome)
    assert gene_identity.cta_annotation_gene_ids(genome, unfiltered=True) == {PRAME, ALTERNATE}


@pytest.mark.parametrize("release", [93, 112])
def test_prame_real_annotation_sources(release):
    from pyensembl import EnsemblRelease

    genome = EnsemblRelease(release, species="human")
    if not genome.required_local_files_exist():
        pytest.skip(f"Ensembl {release} not installed")
    primary = gene_identity.resolve_gene_identity(PRAME, genome=genome)
    assert primary.verified and primary.canonical_gene_id == PRAME
    record = gene_identity.resolve_gene_identity(ALTERNATE, genome=genome)
    if ALTERNATE not in genome.gene_ids():
        # The legacy CI mirror contains primary-assembly genes only. An absent
        # locus must not be invented from the shipped alias table.
        assert record.status == "absent" and not record.verified
        assert ALTERNATE not in gene_identity.cta_annotation_gene_ids(genome, unfiltered=True)
        return
    assert record.verified and record.canonical_gene_id == PRAME
    assert ALTERNATE in gene_identity.cta_annotation_gene_ids(genome, unfiltered=True)
    transcripts = {t.transcript_id: t.protein_id for t in genome.gene_by_id(ALTERNATE).transcripts}
    assert transcripts["ENST00000539862"] == "ENSP00000445097"
    assert transcripts["ENST00000617728"] == "ENSP00000484066"


def test_unsupported_mapping_method_is_not_verified(mapping_tables):
    aliases, _ = mapping_tables
    aliases.loc[aliases.alt_haplotype_id.eq(ALTERNATE), "source"] = "ensembl_archive_replacement"
    record = gene_identity.resolve_gene_identity(
        ALTERNATE, genome=Annotation([gene(ALTERNATE, contig="ALT")])
    )
    assert record.status == "unverified" and not record.verified


def test_annotation_species_is_required():
    genome = Annotation([])
    genome.species = SimpleNamespace(latin_name="mus_musculus")
    with pytest.raises(ValueError, match="human"):
        gene_identity.resolve_gene_identity(PRAME, genome=genome)


def test_cta_curation_and_bundle_versions_unchanged():
    from oncoref.version import DATA_VERSION, SOURCE_MATRIX_VERSION

    assert len(oncoref.cta_gene_ids()) == 624
    assert len(oncoref.cta_unfiltered_gene_ids()) == 2532
    assert DATA_VERSION == "5.23.25"
    assert SOURCE_MATRIX_VERSION == "5.22.14"


@pytest.mark.parametrize("release", [93, 112])
def test_pinned_real_alternate_annotation_and_genuine_non_cta_source(release):
    path = Path(__file__).parent / "fixtures" / "gene-identity" / f"prame-ensembl-{release}.json"
    snapshot = json.loads(path.read_text())
    genome = Annotation([SimpleNamespace(**row) for row in snapshot["genes"]])
    genome.release = snapshot["annotation_release"]
    genome.reference_name = snapshot["annotation_assembly"]
    genome.species = SimpleNamespace(latin_name=snapshot["annotation_species"])
    record = gene_identity.resolve_gene_identity(ALTERNATE, genome=genome)
    assert record.verified and record.canonical_gene_id == PRAME
    assert record.source_contig == ("CHR_HSCHR22_1_CTG3" if release == 93 else "HSCHR22_1_CTG3")
    assert len(snapshot["original_file_identities"]) == 2
    assert all(len(row["sha256"]) == 64 for row in snapshot["original_file_identities"])
    assert gene_identity.cta_annotation_gene_ids(genome, unfiltered=True) == {PRAME, ALTERNATE}
    transcripts = {
        t["transcript_id"]: t["protein_id"] for t in genome.gene_by_id(ALTERNATE).transcripts
    }
    assert transcripts["ENST00000539862"] == "ENSP00000445097"
    assert transcripts["ENST00000617728"] == "ENSP00000484066"
    sequences = {
        row.gene_id: max((t["protein_sequence"] for t in row.transcripts), key=len)
        for row in genome.genes.values()
    }

    def kmers(sequence):
        return {sequence[i : i + 9] for i in range(len(sequence) - 8)}

    primary_kmers = kmers(sequences[PRAME])
    assert len(primary_kmers) == 501
    assert kmers(sequences[ALTERNATE]) == primary_kmers
    # PRAMEF10 is a distinct non-CTA gene whose seven shared peptides stay in
    # the background even though the alternate PRAME assembly copy is excluded.
    assert (
        gene_identity.resolve_gene_identity("ENSG00000187545", genome=genome).canonical_gene_id
        != PRAME
    )
    assert len(primary_kmers & kmers(sequences["ENSG00000187545"])) == 7
