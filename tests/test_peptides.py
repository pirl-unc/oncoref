# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

from types import SimpleNamespace

import pandas as pd
import pytest

from oncoref import peptides
from oncoref.gene_identity import resolve_gene_identity as peptides_gene_identity


class _FakeTr:
    def __init__(self, gene_id, biotype, protein_sequence, contig="22", name="PRAME"):
        self.gene_id = gene_id
        self.contig = contig
        self.gene_name = name
        self.biotype = biotype
        self.protein_sequence = protein_sequence


class _FakeGenome:
    release = 999
    reference_name = "GRCh38"
    species = SimpleNamespace(latin_name="homo_sapiens")

    def __init__(self, transcripts):
        self._trs = transcripts

    def transcripts(self):
        return self._trs

    def gene_ids(self):
        return sorted({t.gene_id for t in self._trs})

    def gene_by_id(self, gid):
        for t in self._trs:
            if t.gene_id == gid:
                return t
        raise ValueError(gid)


def test_kmers():
    assert peptides._kmers("AAACCC", 3) == {"AAA", "AAC", "ACC", "CCC"}
    assert peptides._kmers("AB", 3) == set()  # shorter than k


def test_longest_protein_per_gene_keeps_longest():
    genome = _FakeGenome(
        [
            _FakeTr("G1", "protein_coding", "AAAA"),
            _FakeTr("G1", "protein_coding", "AAAAAA"),  # longer -> wins
            _FakeTr("G1", "lncRNA", "ZZZZZZZZ"),  # non-coding -> ignored
            _FakeTr("G2", "protein_coding", "CC"),  # shorter than k -> dropped
            _FakeTr("G3", "protein_coding", "MKLP*"),  # stop codon stripped
        ]
    )
    longest = peptides._longest_protein_per_gene(genome, 3)
    assert longest["G1"] == "AAAAAA"
    assert "G2" not in longest
    assert longest["G3"] == "MKLP"


@pytest.fixture
def fake_proteome(monkeypatch, tmp_path):
    # CTA protein "AAACCC" (3-mers: AAA,AAC,ACC,CCC); background "CCCDDD" shares CCC.
    genome = _FakeGenome(
        [
            _FakeTr("ENSG00000185686", "protein_coding", "AAACCC"),
            _FakeTr("ENSG00000141510", "protein_coding", "CCCDDD", "17", "TP53"),
            _FakeTr("ENSG_NC", "lncRNA", "EEEEEE"),
        ]
    )
    monkeypatch.setattr(peptides, "_usable_genome", lambda: genome)
    monkeypatch.setattr(peptides, "_derived_cache_dir", lambda: tmp_path)
    monkeypatch.setattr(peptides, "cta_gene_ids", lambda: ["ENSG00000185686"])
    monkeypatch.setattr(peptides, "cta_unfiltered_gene_ids", lambda: ["ENSG00000185686"])
    monkeypatch.setattr(peptides, "cta_gene_id_to_name", lambda: {"ENSG00000185686": "CTAX"})
    monkeypatch.setattr(peptides, "_MIN_PROTEOME_GENES", 1)  # tiny fake proteome
    peptides._COUNTS_CACHE.clear()
    monkeypatch.setattr(
        peptides,
        "cta_annotation_gene_identities",
        lambda g: (peptides_gene_identity("ENSG00000185686", genome=g),),
    )
    yield
    peptides._COUNTS_CACHE.clear()


def test_specific_9mer_counts_subtracts_background(fake_proteome):
    df = peptides.cta_specific_9mer_counts(k=3)
    assert list(df["Symbol"]) == ["CTAX"]
    row = df.iloc[0]
    assert row["n_9mers"] == 4  # AAA, AAC, ACC, CCC
    assert row["n_specific_9mers"] == 3  # CCC also in background -> excluded


def test_specific_9mer_counts_caches_with_fingerprint(fake_proteome, tmp_path):
    peptides.cta_specific_9mer_counts(k=3)
    # filename is keyed by release AND a CTA-set fingerprint
    cached = list(tmp_path.glob("cta_specific_3mers_r999_*.csv"))
    assert len(cached) == 1


def test_specific_9mer_counts_returns_fresh_copy(fake_proteome):
    a = peptides.cta_specific_9mer_counts(k=3)
    a.loc[0, "n_specific_9mers"] = -999  # mutating the result must not corrupt the cache
    b = peptides.cta_specific_9mer_counts(k=3)
    assert b.loc[0, "n_specific_9mers"] == 3


def test_specific_9mer_counts_refresh_rebuilds(fake_proteome, tmp_path):
    peptides.cta_specific_9mer_counts(k=3)
    # corrupt the on-disk cache; refresh=True must drop it and rebuild correctly
    (next(tmp_path.glob("cta_specific_3mers_r999_*.csv"))).write_text("garbage\n")
    peptides._COUNTS_CACHE.clear()
    df = peptides.cta_specific_9mer_counts(k=3, refresh=True)
    assert int(df.loc[0, "n_specific_9mers"]) == 3


def test_specific_9mer_weights_keyed_by_proteoform_key_by_default(fake_proteome):
    # Default key is the proteoform_key; for the singleton ENSG00000185686 that is its ENSG.
    assert peptides.cta_specific_9mer_weights(k=3) == {"ENSG00000185686": 3}
    assert peptides.cta_specific_9mer_weights(k=3, by="ensembl_gene_id") == {"ENSG00000185686": 3}
    assert peptides.cta_specific_9mer_weights(k=3, by="symbol") == {"CTAX": 3}
    with pytest.raises(ValueError, match="by must be"):
        peptides.cta_specific_9mer_weights(k=3, by="nonsense")


def test_weight_by_proteoform_key_covers_group_via_any_member(monkeypatch):
    # A real CGB3/5/8 group whose canonical min-ENSG (CGB3) is unexpressed: the counts
    # table has only the expressed member (CGB8), yet keying by proteoform_key gives the
    # group its weight — no canonical-member ambiguity.
    import pandas as pd

    fake_counts = pd.DataFrame(
        {
            "Ensembl_Gene_ID": ["ENSG00000213030"],  # CGB8 (expressed member of CGB3/5/8)
            "Symbol": ["CGB8"],
            "n_9mers": [40],
            "n_specific_9mers": [31],
        }
    )
    monkeypatch.setattr(
        peptides, "cta_specific_9mer_counts", lambda *, k=peptides.DEFAULT_K: fake_counts
    )
    weights = peptides.cta_specific_9mer_weights(by="proteoform_key")
    assert weights == {"CGB3/5/8": 31}


def test_specific_9mer_load_joins_on_proteoform_key(fake_proteome, monkeypatch):
    # The load joins on proteoform_key (the uniform key), NOT Symbol/ENSG: a stub whose
    # Symbol/ENSG don't match the weight key still resolves via proteoform_key.
    from oncoref import coverage

    def fake_fractions(code, *, threshold_tpm):
        return pd.DataFrame(
            {
                "proteoform_key": ["ENSG00000185686"],  # the weight key (singleton -> ENSG)
                "Ensembl_Gene_ID": ["ENSG00000185686"],
                "Symbol": ["SOMETHING_ELSE"],  # deliberately not the join key
                "fraction_expressing": [0.5],
            }
        )

    monkeypatch.setattr(coverage, "cta_patient_fractions", fake_fractions)
    # load = fraction (0.5) * n_specific_9mers (3) = 1.5
    assert peptides.cta_specific_9mer_load("X", threshold_tpm=5, k=3) == pytest.approx(1.5)


def test_specific_9mer_load_empty_cohort(fake_proteome, monkeypatch):
    from oncoref import coverage

    monkeypatch.setattr(coverage, "cta_patient_fractions", lambda code, **k: pd.DataFrame())
    assert peptides.cta_specific_9mer_load("X", k=3) == 0.0


def test_verified_alt_copy_does_not_mask_specificity_but_real_background_does(monkeypatch):
    primary, alternate, background = "ENSG00000185686", "ENSG00000275013", "ENSG00000141510"
    genome = _FakeGenome(
        [
            _FakeTr(primary, "protein_coding", "AAACCC"),
            _FakeTr(alternate, "protein_coding", "AAACCC", "ALT"),
            _FakeTr(background, "protein_coding", "CCCDDD", "17", "TP53"),
        ]
    )
    monkeypatch.setattr(peptides, "cta_gene_ids", lambda: {primary})
    monkeypatch.setattr(peptides, "cta_unfiltered_gene_ids", lambda: {primary})
    monkeypatch.setattr(peptides, "_MIN_PROTEOME_GENES", 1)
    counts = peptides._build_counts(genome, 3).set_index("Ensembl_Gene_ID")
    assert counts.loc[primary, "n_specific_9mers"] == 3
    genome._trs[-1].protein_sequence = "AAACCC"
    assert peptides._build_counts(genome, 3).iloc[0].n_specific_9mers == 0
    genome._trs = genome._trs[1:2]  # verified alternate is the only annotated copy
    counts = peptides._build_counts(genome, 3)
    assert counts.iloc[0].n_9mers == counts.iloc[0].n_specific_9mers == 4


def test_conflicting_alias_stays_in_background(monkeypatch):
    primary, alternate = "ENSG00000185686", "ENSG00000275013"
    genome = _FakeGenome(
        [
            _FakeTr(primary, "protein_coding", "AAACCC"),
            _FakeTr(alternate, "protein_coding", "AAACCC", "ALT", "OTHER"),
        ]
    )
    monkeypatch.setattr(peptides, "cta_gene_ids", lambda: {primary})
    monkeypatch.setattr(peptides, "cta_unfiltered_gene_ids", lambda: {primary})
    monkeypatch.setattr(peptides, "_MIN_PROTEOME_GENES", 1)
    assert peptides._build_counts(genome, 3).iloc[0].n_specific_9mers == 0


def test_identity_change_invalidates_existing_count_cache(fake_proteome, monkeypatch, tmp_path):
    from dataclasses import replace

    peptides.cta_specific_9mer_counts(k=3)
    original = peptides.cta_annotation_gene_identities
    monkeypatch.setattr(
        peptides,
        "cta_annotation_gene_identities",
        lambda g: tuple(replace(record, aliases_sha256="f" * 64) for record in original(g)),
    )
    peptides.cta_specific_9mer_counts(k=3)
    assert len(list(tmp_path.glob("cta_specific_3mers_r999_*.csv"))) == 2
