import json
from types import SimpleNamespace

import pandas as pd
import pytest

from oncoref import antigen_sets


def test_current_memberships_and_forbidden_background():
    sets = {g: set(antigen_sets.proteins(g).entity_id) for g in antigen_sets.SETS}
    assert sets["CTA"] == sets["CA"] & sets["RA"]
    assert sets["CA|RA"] == sets["CA"] | sets["RA"]
    assert not sets["CA|RA"] & sets["Unsafe"]
    assert sets["Unsafe|Other"] == sets["Unsafe"] | sets["Other"]
    assert not sets["Unsafe"] & sets["Other"]
    assert not sets["CA|RA"] & sets["Unsafe|Other"]
    assert sets["CA|RA"] | sets["Unsafe|Other"] == set(
        antigen_sets.target_selection.results().entity_id
    )
    assert antigen_sets.gene_ids("CTA") == {
        g for ids in antigen_sets.proteins("CTA").member_gene_ids for g in json.loads(ids)
    }
    with pytest.raises(ValueError):
        antigen_sets.proteins("Neither")


def test_sequences_include_every_annotated_translation_and_report_gaps(monkeypatch):
    monkeypatch.setattr(
        antigen_sets,
        "proteins",
        lambda group: pd.DataFrame(
            {
                "entity_id": ["p1", "p2"],
                "member_gene_ids": ['["g1"]', '["g2"]'],
            }
        ),
    )
    transcripts = [
        SimpleNamespace(id="t1", biotype="protein_coding", protein_sequence="ACDEFGHIK"),
        SimpleNamespace(id="t2", biotype="protein_coding_LoF", protein_sequence="ACDEFGHIK"),
        SimpleNamespace(id="t3", biotype="nonsense_mediated_decay", protein_sequence="LMNPQRSTV"),
        SimpleNamespace(id="t4", biotype="retained_intron", protein_sequence=None),
    ]
    genome = SimpleNamespace(
        release=115,
        species=SimpleNamespace(latin_name="homo_sapiens"),
        gene_by_id=lambda gene: SimpleNamespace(transcripts=transcripts if gene == "g1" else []),
    )
    result = antigen_sets.sequences("Unsafe", genome=genome)
    assert list(result.transcript_id) == ["t1", "t2", "t3"]
    assert set(result.transcript_biotype) == {
        "protein_coding",
        "protein_coding_LoF",
        "nonsense_mediated_decay",
    }
    assert result.attrs["missing_gene_ids"] == ["g2"]
    assert result.attrs["sequence_scope"] == "all annotated protein translations"
    assert result.sequence_sha256.iloc[0] == result.sequence_sha256.iloc[1]
