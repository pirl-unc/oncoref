import json

import pytest

from oncoref import antigen_sets


def test_current_memberships_and_forbidden_background():
    sets = {g: set(antigen_sets.proteins(g).entity_id) for g in antigen_sets.SETS}
    assert sets["CTA"] == sets["CA"] & sets["RA"]
    assert sets["CA|RA"] == sets["CA"] | sets["RA"]
    assert not sets["CA|RA"] & sets["Unsafe"]
    assert antigen_sets.gene_ids("CTA") == {
        g for ids in antigen_sets.proteins("CTA").member_gene_ids for g in json.loads(ids)
    }
    with pytest.raises(ValueError):
        antigen_sets.proteins("Neither")
