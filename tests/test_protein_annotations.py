from oncoref.protein_annotations import _gene_rows


def test_shared_entry_keeps_exact_loci_and_evidence():
    entry = {
        "primaryAccession": "P00001",
        "genes": [{"geneName": {"value": "DO_NOT_MAP_BY_SYMBOL"}}],
        "uniProtKBCrossReferences": [
            {
                "database": "Ensembl",
                "properties": [
                    {"key": "GeneId", "value": "ENSG00000000001.3"},
                ],
            },
            {
                "database": "Ensembl",
                "properties": [
                    {"key": "GeneId", "value": "ENSG00000000002.4"},
                ],
            },
        ],
        "comments": [
            {
                "commentType": "FUNCTION",
                "texts": [
                    {
                        "value": "May participate in repair (By similarity).",
                        "evidences": [{"source": "PubMed", "id": "123"}],
                    }
                ],
            }
        ],
    }
    rows = _gene_rows([entry])
    assert {row["gene_id"] for row in rows} == {"ENSG00000000001", "ENSG00000000002"}
    assert all(row["shared_entry"] for row in rows)
    assert all(row["function"].endswith("(By similarity).") for row in rows)
    assert all(row["function_pubmed_ids"] == "123" for row in rows)
    entry["uniProtKBCrossReferences"] = []
    assert _gene_rows([entry]) == []
