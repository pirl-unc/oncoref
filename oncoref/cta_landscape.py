"""Reproducible intake of complete published CTA landscape nominations.

Paper membership is nomination evidence. Shared probes retain every annotated
member and their assay scope; no membership asserts unique peptide presentation.
"""

from __future__ import annotations

import hashlib
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

from .gene_ids import canonical_gene_id, canonical_gene_space, resolve_symbol

INPUTS = {
    "PMC4737856-supp/ncomms10499-s4.xlsx": "4d5a8a68dba0c5499074b79796429f1c61f7b043130c0add6f629bb9449e5e87",
    "PMC6193945-supp/41388_2018_357_MOESM4_ESM.xlsx": "0d60e63830e622c3bcdfdeb4152767904ae09eb06ce144ffae1af9b32579e586",
    "PMC5696236-supp/oncotarget-08-92966-s002.xls": "4b83d4b1ec8102a3719e79e24e3cbb037a439264f130b62a6e45ddbd32b9fdb5",
    "PMC5696236-supp/oncotarget-08-92966-s003.xls": "5198424def9b66c8b24f662490b557126ed334da99667d028426298f2e7e2839",
    "PMC8564638-supp/MOL2-15-3003-s007.xlsx": "f1620f3a66b288fcb8e0b06d13579437a6c287c8e5e2041f3b458867c110d9bb",
    "PMC8564638-supp/MOL2-15-3003-s003.xlsx": "27c88594f2b3ad9a2a4fcae10ab1f28169500041802b8cc5f133faf1604d8f48",
    "PMC6601584-supp/CAM4-8-3511-s002.xlsx": "19c1cb9825eda2c42ca727b08ea4b1d928d11dde01893087cfadf8e438548941",
    "PMC10749072-supp/jitc-2023-007935supp001.pdf": "87b3c5d3bb8ef245b436aae186235081d42b220309fa170fce5373f6b0b8ca33",
    "PMC10851610.xml": "ea59a269180b6bd9c40cd74951c042850f95560331ae9ee638cffa4de2e511eb",
    "PMC5696236-supp/oncotarget-08-92966-s004.xlsx": "8f70d3ab237c7df0f1b1885705998e6a6c97d4d0c9261814750e5fef3bedc31f",
    "gong.xlsx": "b91d126e9e0a270309cb3f6ad5253f1153db3a21f0a8c19d8ef55080410a2007",
    "loriot-s1.xlsx": "5d5eb8777d16d48eee3c50b556c805f3efd75d2a6de7c8884e40bf374e51ef65",
    "loriot-s2.xlsx": "03188e479ca52419b30162004344105f1241c44cebf37f2f11d6bdd74dbd5a77",
    "PMC5120866-bioc.xml": "e3217c6848f79582427787ca4c64fd1852d97a4c67ee6310ee9a6e8973d544a9",
}

INPUT_URLS = {
    "gong.xlsx": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-021-22695-y/MediaObjects/41467_2021_22695_MOESM2_ESM.xlsx",
    "loriot-s1.xlsx": "https://doi.org/10.1371/journal.pgen.1011734.s004",
    "loriot-s2.xlsx": "https://doi.org/10.1371/journal.pgen.1011734.s005",
    "PMC5120866-bioc.xml": "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC5120866/unicode",
}

# A paper can contribute both its broad nomination and a stricter nested set.
PAPERS = {
    "Wang2016_CT": (
        "Wang 2016",
        "10.1038/ncomms10499",
        "Systematic identification of genes with a cancer-testis expression pattern in 19 cancer types",
        1019,
    ),
    "Bruggeman2018_GC": (
        "Bruggeman 2018",
        "10.1038/s41388-018-0357-2",
        "Massive expression of germ cell-specific genes is a hallmark of cancer and a potential target for novel treatment development",
        756,
    ),
    "daSilva2017_testis_biased": (
        "da Silva 2017 · testis-biased",
        "10.18632/oncotarget.21715",
        "Genome-wide identification of cancer/testis genes and their association with prognosis in a pan-cancer analysis",
        1103,
    ),
    "daSilva2017_CT": (
        "da Silva 2017 · CT",
        "10.18632/oncotarget.21715",
        "Genome-wide identification of cancer/testis genes and their association with prognosis in a pan-cancer analysis",
        745,
    ),
    "Jamin2021_CT": (
        "Jamin 2021",
        "10.1002/1878-0261.12900",
        "Combined RNA/tissue profiling identifies novel Cancer/testis genes",
        478,
    ),
    "Jamin2021_core": (
        "Jamin 2021 · core",
        "10.1002/1878-0261.12900",
        "Combined RNA/tissue profiling identifies novel Cancer/testis genes",
        124,
    ),
    "Chang2019_TGCT": (
        "Chang 2019 · TGCT",
        "10.1002/cam4.2223",
        "Comprehensive characterization of cancer-testis genes in testicular germ cell tumor",
        1036,
    ),
    "Carter2023_CT": (
        "Carter 2023",
        "10.1136/jitc-2023-007935",
        "Identification of pan-cancer/testis genes and validation of therapeutic targeting in triple-negative breast cancer",
        103,
    ),
    "Seager2024_CTA": (
        "Seager 2024 · panel",
        "10.1186/s12967-024-04918-0",
        "Cancer testis antigen burden (CTAB): a novel biomarker of tumor-associated antigens in lung cancer",
        17,
    ),
}

PAPERS.update(
    {
        "Loriot2025_S1": (
            "Loriot 2025 · strict",
            "10.1371/journal.pgen.1011734",
            "A survey of human cancer-germline genes: Linking X chromosome localization, DNA methylation and sex-biased expression in early embryos",
            146,
        ),
        "Loriot2025_S2": (
            "Loriot 2025 · preferential",
            "10.1371/journal.pgen.1011734",
            "A survey of human cancer-germline genes: Linking X chromosome localization, DNA methylation and sex-biased expression in early embryos",
            134,
        ),
        "Gong2021_reproductive_PC": (
            "Gong 2021 · reproductive coding",
            "10.1038/s41467-021-22695-y",
            "The RNA landscape of the human placenta in health and disease",
            744,
        ),
        "Gong2021_reproductive_ncRNA": (
            "Gong 2021 · reproductive noncoding",
            "10.1038/s41467-021-22695-y",
            "The RNA landscape of the human placenta in health and disease",
            1067,
        ),
        "daSilva2017_tumor_proteomics": (
            "da Silva 2017 · tumor proteomics",
            "10.18632/oncotarget.21715",
            "Genome-wide identification of cancer/testis genes and their association with prognosis in a pan-cancer analysis",
            136,
        ),
        "Bai2016_EGFL6": (
            "Bai 2016",
            "10.1158/0008-5472.CAN-16-0225",
            "EGFL6 regulates the asymmetric division, maintenance and metastasis of ALDH+ ovarian cancer cells",
            1,
        ),
    }
)


def _text(value):
    return "" if value is None else str(value).strip()


def _xlsx(path, sheet):
    import openpyxl

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        return list(enumerate(workbook[sheet].values, 1))
    finally:
        workbook.close()


def _xls(path, sheet):
    import xlrd

    workbook = xlrd.open_workbook(path)
    tab = (
        workbook.sheet_by_index(sheet) if isinstance(sheet, int) else workbook.sheet_by_name(sheet)
    )
    return [(i + 1, tab.row_values(i)) for i in range(tab.nrows)]


def resolve_landscape_membership(rows):
    """Use pinned ID/synonym maps only; ambiguous/unknown identities stay visible.

    Historical Ensembl and Entrez identifiers never fall back to a matching symbol.
    Symbol-only papers use unambiguous canonical symbols or shipped NCBI synonyms;
    no live annotation download or silently chosen paralog is permitted.
    """
    result = rows.copy()
    space = canonical_gene_space().set_index("ensembl_gene_id")
    symbols = defaultdict(set)
    for gid, symbol in space.symbol.items():
        if pd.notna(symbol) and symbol:
            symbols[str(symbol).upper()].add(gid)
    mapped = []
    for row in result.itertuples(index=False):
        identifier = _text(row.source_gene_id)
        if row.source_id_type in {"ensembl", "entrez"}:
            gid = canonical_gene_id(identifier)
            method = f"{row.source_id_type}_reference_map"
        elif row.source_id_type == "symbol":
            direct = symbols.get(row.source_symbol.upper(), set())
            official = resolve_symbol(row.source_symbol)
            ids = direct or symbols.get(official.upper(), set())
            gid = next(iter(ids)) if len(ids) == 1 else None
            method = "exact_symbol" if direct else "NCBI_symbol_synonym"
        else:
            gid, method = None, "unannotated_probe"
        label = space.loc[gid, "symbol"] if gid else ""
        if gid and (pd.isna(label) or not label):
            label = gid  # A canonical locus can have no approved symbol.
        mapped.append(
            (
                gid or "",
                label,
                space.loc[gid, "biotype"] if gid else "",
                method if gid else "unmapped",
            )
        )
    result[["Ensembl_Gene_ID", "Symbol", "biotype", "mapping_method"]] = pd.DataFrame(
        mapped, index=result.index
    )
    result["mapping_status"] = result.Ensembl_Gene_ID.ne("").map(
        {True: "mapped", False: "unmapped"}
    )
    result["protein_candidate_eligible"] = result.biotype.eq("protein_coding")
    return result


def extract_landscapes(input_dir):
    """Return complete source memberships and source metadata from pinned inputs.

    The CTA intake uses Wang's 1,019 coding CT genes, Bruggeman's 756 GC genes,
    da Silva's broad 1,103 and tumor-positive 745 sets, all 602 Jamin CT probes
    (every annotated member plus unresolved probes), Chang's 1,036 expressed C1
    coding genes, Carter's final 103, and Seager's 17-target panel. Noncoding-only
    CT-RNA screens are not protein antigen nominations and are outside this intake.
    """
    root = Path(input_dir)
    for name, expected in INPUTS.items():
        actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"Published input checksum mismatch: {name}")
    rows = []

    def add(
        tag,
        file,
        sheet,
        row,
        symbol,
        gid="",
        id_type="symbol",
        *,
        scope="Cancer-associated RNA nomination; no antigen validation asserted",
        subset="selected_CT",
        biotype="not specified",
        assay="gene-labelled RNA; paralog resolution not established",
    ):
        rows.append(
            {
                "source_tag": tag,
                "source_gene_id": _text(gid),
                "source_symbol": _text(symbol),
                "source_id_type": id_type,
                "source_biotype": biotype,
                "evidence_type": "published_CTA_nomination",
                "validation_scope": scope,
                "source_file": Path(file).name,
                "source_sheet": sheet,
                "source_rows": str(row),
                "source_subset": subset,
                "assay_resolution": assay,
            }
        )

    file = "PMC4737856-supp/ncomms10499-s4.xlsx"
    for n, r in _xlsx(root / file, "Supplementary Data 3A"):
        if re.fullmatch(r"ENSG\d{11}", _text(r[0])):
            add(
                "Wang2016_CT",
                file,
                "Supplementary Data 3A",
                n,
                r[1],
                r[0],
                "ensembl",
                biotype="protein_coding",
            )

    file = "PMC6193945-supp/41388_2018_357_MOESM4_ESM.xlsx"
    tgct = {_text(r[0]) for n, r in _xlsx(root / file, "2") if n > 2}
    for n, r in _xlsx(root / file, "1D"):
        if n > 2 and r[0]:
            scope = (
                "GC RNA nomination; inclusion depends on TGCT under this study's thresholds"
                if r[0] in tgct
                else "GC RNA nomination across cancer types; no antigen validation asserted"
            )
            add(
                "Bruggeman2018_GC",
                file,
                "1D",
                n,
                r[0],
                scope=scope,
                subset="TGCT_dependent" if r[0] in tgct else "non_TGCT",
            )

    file = "PMC5696236-supp/oncotarget-08-92966-s002.xls"
    broad = {
        r[1]: (n, str(int(r[0])))
        for n, r in _xls(root / file, "1103 CTs predicted")
        if n > 1 and r[1]
    }
    broad_upper = {symbol.upper(): value for symbol, value in broad.items()}
    if len(broad_upper) != len(broad):
        raise ValueError("Ambiguous case-insensitive symbols in da Silva Table 1")
    for symbol, (n, gid) in broad.items():
        add(
            "daSilva2017_testis_biased",
            file,
            "1103 CTs predicted",
            n,
            symbol,
            gid,
            "entrez",
            scope="At least 90% normal-tissue RNA attributable to testis; tumor positivity not required",
            subset="testis_biased_nomination",
        )
    file = "PMC5696236-supp/oncotarget-08-92966-s003.xls"
    for n, r in _xls(root / file, 0):
        if n <= 3:
            continue
        percentages = [v for v in r[2::2] if isinstance(v, (int, float))]
        if percentages and max(percentages) >= 10:
            # Table 2 contains TRIM49D2P whereas Table 1 labels TRIM49D2.
            # Preserve the discrepancy; do not invent an Entrez join.
            matched = broad_upper.get(r[0].upper())
            add(
                "daSilva2017_CT",
                file,
                "Sheet1",
                n,
                r[0],
                matched[1] if matched else "",
                "entrez" if matched else "symbol",
                scope="RSEM >1 in at least 10% of samples in a tumor type; RNA nomination",
            )

    file = "PMC8564638-supp/MOL2-15-3003-s007.xlsx"
    probes = []
    for n, r in _xlsx(root / file, "Expression data and classes"):
        if n == 1 or r[6] not in {"SET", "SEHET", "PET", "PEHET"} or _text(r[11]) in {"", "-"}:
            continue
        probes.append(r[0])
        members = _text(r[1]).split(";")
        for symbol in members:
            unnamed = symbol == "-"
            add(
                "Jamin2021_CT",
                file,
                "Expression data and classes",
                n,
                f"unannotated:{r[0]}" if unnamed else symbol,
                r[0] if unnamed else "",
                "probe" if unnamed else "symbol",
                subset=r[6],
                scope="CT-pattern microarray nomination; tumor upregulation with no corresponding normal-tissue detection",
                assay=f"Affymetrix {r[0]}: {_text(r[1])}; shared probe"
                if len(members) > 1
                else f"Affymetrix {r[0]}: {_text(r[1])}",
            )
    if len(probes) != 602:
        raise ValueError(f"Expected 602 Jamin CT probes, found {len(probes)}")
    file = "PMC8564638-supp/MOL2-15-3003-s003.xlsx"
    for n, r in _xlsx(root / file, "Annotation"):
        if re.fullmatch(r"ENSG\d{11}", _text(r[0])):
            add(
                "Jamin2021_core",
                file,
                "Annotation",
                n,
                r[1],
                r[0],
                "ensembl",
                subset="core_CT",
                scope="Core CT annotation; underlying microarray probes may recognize several loci",
            )

    file = "PMC6601584-supp/CAM4-8-3511-s002.xlsx"
    for n, r in _xlsx(root / file, "Suppl Table 2"):
        if (
            r[2] == "C1"
            and r[5] == "protein_coding"
            and r[3] not in (None, "NA")
            and float(str(r[3]).split("%")[0]) >= 1
        ):
            add(
                "Chang2019_TGCT",
                file,
                "Suppl Table 2",
                n,
                r[1],
                r[0],
                "ensembl",
                biotype="protein_coding",
                scope="Testis-specific coding gene expressed in at least 1% of TGCT; no broad somatic-cancer claim",
                subset="TGCT_coding",
            )

    from pypdf import PdfReader

    file = "PMC10749072-supp/jitc-2023-007935supp001.pdf"
    pattern = re.compile(
        r"\s*(ENSG\d+)\s+(.+?)\s+(\d[\d.]*)\s+(True|False)\s+(True|False)\s+(True|False)\s*$"
    )
    for page, content in enumerate(PdfReader(root / file).pages, 1):
        if not 17 <= page <= 27:
            continue
        for n, line in enumerate(content.extract_text().splitlines(), 1):
            match = pattern.fullmatch(line)
            if match:
                gid, symbol, _, thymus, _, selected = match.groups()
                if selected == "True" and thymus == "False":
                    add(
                        "Carter2023_CT",
                        file,
                        "Supplementary Table 2",
                        f"PDF page {page}, text line {n}",
                        symbol,
                        gid,
                        "ensembl",
                        scope="Final CT-expression set without detectable thymic expression; functional vaccination evidence is not assigned to every gene",
                        biotype="protein_coding",
                    )

    # Literal panel from Methods; preserve the historical GAGE2 label and flag
    # group-level measurement rather than pretending it is a unique-locus assay.
    file = "PMC10851610.xml"
    panel = (
        "BAGE",
        "CTAG1B",
        "CTAG2",
        "GAGE1",
        "GAGE10",
        "GAGE12J",
        "GAGE13",
        "GAGE2",
        "MAGEA1",
        "MAGEA10",
        "MAGEA12",
        "MAGEA3",
        "MAGEA4",
        "MAGEC2",
        "MLANA",
        "SSX2",
        "XAGE1B",
    )
    for n, symbol in enumerate(panel, 1):
        add(
            "Seager2024_CTA",
            file,
            "Methods: Data processing and statistical analysis; Figure 1",
            n,
            symbol,
            scope="Member of the published 17-target CTA profiling panel; includes differentiation antigen MLANA",
            assay="Reported RNA-panel target label; related-gene assays may be shared",
            subset="profiled_panel",
        )

    for number, sheet in [(1, "List of 146 CG genes"), (2, "List of 134 CG-Preferential gen")]:
        file = f"loriot-s{number}.xlsx"
        for n, r in _xlsx(root / file, sheet):
            if n > 2:
                add(
                    f"Loriot2025_S{number}",
                    file,
                    sheet,
                    n,
                    r[2],
                    r[1],
                    "ensembl",
                    subset="strict_CG" if number == 1 else "preferential_CG",
                    scope="Cancer-germline RNA nomination; strict and preferential classes remain distinct",
                )

    file = "gong.xlsx"
    for sheet, kind in [
        ("Data 5 - tissue-enriched PC", "PC"),
        ("Data 6 tissue-enriched lncR", "ncRNA"),
    ]:
        for n, r in _xlsx(root / file, sheet):
            if r[0] in {"Testis", "Placenta", "Ovary"}:
                add(
                    f"Gong2021_reproductive_{kind}",
                    file,
                    sheet,
                    n,
                    r[2],
                    r[1],
                    "ensembl",
                    subset=r[0],
                    biotype="protein_coding" if kind == "PC" else "lncRNA",
                    scope="Normal reproductive-tissue enrichment; no tumor validation",
                )

    file = "PMC5696236-supp/oncotarget-08-92966-s004.xlsx"
    for n, r in _xlsx(root / file, "Supplementary_Table3"):
        if n > 2 and r[0]:
            add(
                "daSilva2017_tumor_proteomics",
                file,
                "Supplementary_Table3",
                n,
                r[0],
                subset="tumor_proteomics",
                scope="Gene-labelled tumor proteomics; unique-locus peptide assignment is not established",
                assay="Gene-labelled mass spectrometry; not HLA immunopeptidomics",
            )

    add(
        "Bai2016_EGFL6",
        "PMC5120866-bioc.xml",
        "Abstract",
        1,
        "EGFL6",
        subset="targeted_cancer_expression",
        scope="Ovarian-tumor and vascular expression/function; no normal-tissue exclusivity claim",
        assay="Targeted EGFL6 RNA/protein/functional experiments; no HLA or T-cell validation",
    )

    # A microarray gene can have many selected probes. Preserve every original
    # location/annotation but count its source identity only once in intersections.
    raw = pd.DataFrame(rows)
    grouped = []
    for _, group in raw.groupby(
        ["source_tag", "source_id_type", "source_gene_id", "source_symbol"],
        sort=False,
        dropna=False,
    ):
        row = group.iloc[0].to_dict()
        for col in ("source_rows", "source_subset", "assay_resolution"):
            row[col] = " | ".join(dict.fromkeys(group[col]))
        grouped.append(row)
    membership = resolve_landscape_membership(pd.DataFrame(grouped))
    sizes = membership.groupby("source_tag").size()
    if sizes["Jamin2021_CT"] != 607 or sizes["Jamin2021_core"] != 125:
        raise ValueError("Unexpected Jamin supplement annotation counts")
    for tag, (_, _, _, expected) in PAPERS.items():
        if tag not in {"Jamin2021_CT", "Jamin2021_core"} and sizes[tag] != expected:
            raise ValueError(f"Unexpected {tag} count: {sizes[tag]} != {expected}")
    sources = []
    canonical = canonical_gene_space()
    reference_hash = hashlib.sha256(canonical.to_csv(index=False).encode()).hexdigest()
    for tag, (label, doi, title, reported) in PAPERS.items():
        sub = membership[membership.source_tag.eq(tag)]
        names = list(sub.source_file.unique())
        input_name = next(k for k in INPUTS if Path(k).name == names[0])
        pmc = input_name.split("-")[0].replace(".xml", "")
        sources.append(
            {
                "source_tag": tag,
                "citation": label,
                "title": title,
                "doi": doi,
                "source_url": INPUT_URLS.get(
                    input_name,
                    f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/"
                    + ("fullTextXML" if input_name.endswith(".xml") else "supplementaryFiles"),
                ),
                "source_table": "; ".join(sub.source_sheet.unique()),
                "published_rows": len(sub),
                "paper_reported_genes": reported,
                "input_sha256": INPUTS[input_name],
                "hash_scope": "Original downloaded file bytes",
                "evidence_scope": "; ".join(sub.validation_scope.unique()),
                "canonical_reference_sha256": reference_hash,
                "canonical_ensembl_release": ";".join(
                    sorted(canonical.ensembl_release.astype(str).unique())
                ),
                "retrieved_date": "2026-09-24",
            }
        )
    return membership, pd.DataFrame(sources)
