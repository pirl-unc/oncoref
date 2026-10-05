#!/usr/bin/env python3
"""Rebuild issue #524's reference inventory from the frozen bibliographic snapshot.

This is an offline export of audit flags, citation metadata and manually traced
source chains. It does not revalidate scientific claims or change reference data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DIRECTORY = ROOT / "docs/audits/issue-524-references"
FLAGS = {
    "tmb_source_not_revalidated",
    "numeric_source_location_not_verified",
    "response_denominator_missing",
    "small_denominator",
}
SOURCE_FILES = (
    "docs/audits/tmb-ici-apd1.csv",
    "docs/audits/source-citations.csv",
    "oncoref/data/cancer-tmb.csv",
    "oncoref/data/cancer-tmb-source-audit.csv",
    "oncoref/data/cancer-ici-response-estimates.csv",
    "oncoref/data/cancer-ici-source-locator-audit.csv",
)


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def csv_text(rows, fields=None):
    if fields is None:
        fields = list(dict.fromkeys(key for row in rows for key in row))
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def json_text(value):
    return json.dumps(value, indent=2) + "\n"


def bibliography_html(references, records, traces, summary):
    escape = html.escape
    grouped = defaultdict(list)
    for row in records:
        grouped[row["reference"]].append(row)
    cards = []
    for reference in sorted(
        references,
        key=lambda row: (-row["tmb_record_count"], -row["response_record_count"], row["reference"]),
    ):
        rows = grouped[reference["reference"]]
        groups = " ".join(
            group for group in ("tmb", "response") if reference[f"{group}_record_count"]
        )
        links = [(reference["source_url"], "Citation")]
        if reference["doi"]:
            links.append(("https://doi.org/" + reference["doi"], "Publisher / DOI"))
        if reference["pmcid"]:
            links.append(
                (
                    "https://pmc.ncbi.nlm.nih.gov/articles/" + reference["pmcid"] + "/",
                    "PMC full text",
                )
            )
        linked_urls = {url for url, _ in links}
        for url in reference["followup_source_urls"].split(" | "):
            if url and url not in linked_urls:
                links.append((url, "Traced source"))
                linked_urls.add(url)
        link_html = " · ".join(
            f'<a href="{escape(url, quote=True)}">{escape(label)}</a>' for url, label in links
        )
        fields = (
            "dataset",
            "row_key",
            "cancer_code",
            "metric",
            "value",
            "findings",
            "source_locator",
        )
        table_rows = "".join(
            "<tr>" + "".join(f"<td>{escape(str(row[key]))}</td>" for key in fields) + "</tr>"
            for row in rows
        )
        trace_notes = reference["source_trace_notes"]
        trace_html = f"<p><b>Source trace:</b> {escape(trace_notes)}</p>" if trace_notes else ""
        cards.append(
            f'<article data-group="{groups}" id="{reference["reference_id"]}">'
            f"<h2>{escape(reference['reference_id'])} · {escape(reference['title'])}</h2>"
            f"<p>{escape(reference['reference'])} · {escape(reference['journal'])} · "
            f"{escape(reference['year'])} · {len(rows)} flagged records</p><p>{link_html}</p>"
            f"<p><b>Cancer codes:</b> {escape(reference['cancer_codes'])}</p>"
            f"<p><b>Flags:</b> {escape(reference['review_flags'])}</p>{trace_html}"
            f"<details><summary>Show {len(rows)} mapped records and existing locators</summary>"
            '<div class="table"><table><thead><tr><th>Dataset</th><th>Row key</th>'
            "<th>Cancer</th><th>Metric</th><th>Value</th><th>Flags</th><th>Existing locator</th>"
            f"</tr></thead><tbody>{table_rows}</tbody></table></div></details></article>"
        )
    trace_rows = "".join(
        "<li>" + escape(trace["record_key"] + ": " + trace["trace_notes"]) + "</li>"
        for trace in traces
    )
    missing = "".join(
        "<li>"
        + escape(row["row_key"] + " · " + row["cancer_code"] + " · " + row["notes"])
        + "</li>"
        for row in records
        if not row["reference"]
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Issue 524 reference inventory</title>
<style>
body{{max-width:1200px;margin:40px auto;padding:0 24px;font:16px/1.55 system-ui;color:#182b3a;background:#f7f9fa}}
h1{{font-size:32px}} h2{{font-size:20px}}
article{{background:white;border:1px solid #d7e1e7;border-radius:8px;padding:20px;margin:18px 0}}
a{{color:#085f91}} input,select{{font:inherit;padding:10px;margin:8px 8px 8px 0}}
input{{width:min(600px,70%)}} .table{{overflow-x:auto}}
table{{border-collapse:collapse;font-size:13px}}
th,td{{border:1px solid #d7e1e7;padding:8px;text-align:left;vertical-align:top}}
summary{{cursor:pointer}} #shown{{font-weight:600}}
</style>
</head>
<body>
<h1>Issue #524 · Complete flagged-reference inventory</h1>
<p>{summary["flagged_records"]} flagged records · {len(references)} cited reference identifiers ·
{summary["tmb_references"]} TMB references · {summary["response_references"]} response references ·
{summary["uncited_records"]} records without a direct citation.</p>
<p>Bibliographic snapshot acquired {escape(summary["retrieved_on"])}.
PubMed and the DOI registry establish citation identity and available document links.
Existing numeric locators retain their prior status; exact values, populations and denominators
still require scientific review.</p>
<p><a href="references.csv">Reference CSV</a> ·
<a href="flagged-records.csv">All mapped records</a> ·
<a href="uncited-records.csv">Uncited records</a> ·
<a href="lookup-audit.json">Lookup and coverage audit</a></p>
<details><summary>Recovered source chains and remaining gaps</summary>
<ul>{trace_rows}</ul>
<p><a href="source-traces.csv">Source-chain mapping</a> ·
<a href="additional-primary-reference.csv">Recovered original chordoma reference</a></p>
</details>
<input id="search" aria-label="Search references" placeholder="Search paper, cancer code, PMID, DOI or flag">
<select id="group" aria-label="Reference group"><option value="all">All references</option>
<option value="tmb">TMB</option><option value="response">Response</option></select>
<p id="shown"></p>
{"".join(cards)}
<h2>Records without a direct citation</h2><ul>{missing}</ul>
<script>
const query = document.querySelector('#search');
const group = document.querySelector('#group');
const cards = [...document.querySelectorAll('article')];
function filter() {{
  let count = 0;
  for (const card of cards) {{
    const matchesText = card.textContent.toLowerCase().includes(query.value.toLowerCase());
    const matchesGroup = group.value === 'all' || card.dataset.group.split(' ').includes(group.value);
    card.hidden = !(matchesText && matchesGroup);
    if (!card.hidden) count++;
  }}
  document.querySelector('#shown').textContent = count + ' references shown';
}}
query.addEventListener('input', filter);
group.addEventListener('change', filter);
filter();
</script>
</body>
</html>
"""


def build_inventory(root=ROOT, directory=DEFAULT_DIRECTORY):
    ledger = read_csv(root / SOURCE_FILES[0])
    selected = [row for row in ledger if set(row["findings"].split(";")) & FLAGS]
    if len(selected) != len({(row["dataset"], row["row_key"]) for row in selected}):
        raise ValueError("Duplicate dataset/row keys in the selected audit records")
    grouped = defaultdict(list)
    for row in selected:
        grouped[row["reference"]].append(row)
    identifiers = sorted(set(grouped) - {""})
    snapshot = json.loads((directory / "bibliographic-metadata.json").read_text())
    metadata = {row["reference"]: row for row in snapshot}
    if len(metadata) != len(snapshot):
        raise ValueError("Duplicate citation identifiers in the bibliographic snapshot")
    missing_metadata = set(identifiers) - metadata.keys()
    if missing_metadata:
        raise ValueError(f"Missing bibliographic metadata: {sorted(missing_metadata)}")
    traces = read_csv(directory / "source-traces.csv")
    additional = read_csv(directory / "additional-primary-reference.csv")
    citations = {row["reference"]: row for row in read_csv(root / SOURCE_FILES[1])}
    tmb = {row["cancer_code"]: row for row in read_csv(root / SOURCE_FILES[2])}
    reviews = {row["cancer_code"]: row for row in read_csv(root / SOURCE_FILES[3])}
    estimates = read_csv(root / SOURCE_FILES[4])
    by_id = {row["estimate_id"]: row for row in estimates}
    primary = {
        (row["cancer_code"], row["regimen"]): row
        for row in estimates
        if row["role"] == "primary" and row["metric"] == "ORR"
    }
    locators = {row["estimate_id"]: row for row in read_csv(root / SOURCE_FILES[5])}

    def source_for(row):
        if row["dataset"] == "cancer-tmb":
            return tmb[row["cancer_code"]]
        if row["dataset"] == "cancer-ici-response-estimates":
            return by_id[row["row_key"]]
        return primary[(row["cancer_code"], row["regimen"])]

    references = []
    for index, identifier in enumerate(identifiers, 1):
        reference = metadata[identifier].copy()
        if not reference["title"] or not reference["source_url"].startswith("https://"):
            raise ValueError(f"Missing citation title or HTTPS source URL: {identifier}")
        if identifier.startswith("PMID:") and identifier[5:] != reference["pmid"]:
            raise ValueError(f"PubMed identity mismatch: {identifier}")
        if identifier.startswith("DOI:") and identifier[4:].lower() != reference["doi"].lower():
            raise ValueError(f"DOI identity mismatch: {identifier}")
        rows = grouped[identifier]
        links = set()
        if reference["pmcid"]:
            links.add("https://pmc.ncbi.nlm.nih.gov/articles/" + reference["pmcid"] + "/")
        for row in rows:
            source = source_for(row)
            link = locators.get(source.get("estimate_id"), {}).get("source_document_url", "")
            if link:
                links.add(link)
        matches = [trace for trace in traces if trace["cited_reference"] == identifier]
        reference.update(
            reference_id=f"R{index:03d}",
            full_text_or_existing_document_urls=" | ".join(sorted(links)),
            numeric_claim_review_status="unresolved_audit_flags",
            flagged_record_count=len(rows),
            tmb_record_count=sum(row["dataset"] == "cancer-tmb" for row in rows),
            response_record_count=sum(row["dataset"] != "cancer-tmb" for row in rows),
            cancer_codes="; ".join(sorted({row["cancer_code"] for row in rows})),
            review_flags="; ".join(
                sorted(
                    {flag for row in rows for flag in row["findings"].split(";") if flag in FLAGS}
                )
            ),
            existing_locators=" | ".join(
                sorted({row["source_locator"] for row in rows if row["source_locator"]})
            ),
            inventory_title=citations.get(identifier, {}).get("title", ""),
            followup_source_references="; ".join(
                trace["related_references"] for trace in matches if trace["related_references"]
            ),
            followup_source_urls=" | ".join(trace["related_urls"] for trace in matches),
            source_trace_notes=" | ".join(trace["trace_notes"] for trace in matches),
        )
        references.append(reference)
    reference_map = {row["reference"]: row for row in references}
    records = []
    for row in selected:
        source = source_for(row)
        reference = reference_map.get(row["reference"], {})
        is_tmb = row["dataset"] == "cancer-tmb"
        review = reviews[row["cancer_code"]] if is_tmb else {}
        locator = locators.get(source.get("estimate_id"), {})
        records.append(
            {
                **row,
                "reference_id": reference.get("reference_id", ""),
                "title": reference.get("title", ""),
                "source_url": reference.get("source_url", ""),
                "full_text_or_existing_document_urls": reference.get(
                    "full_text_or_existing_document_urls", ""
                ),
                "citation_status": "cited_reference_found"
                if row["reference"]
                else "no_direct_reference",
                "underlying_estimate_id": source.get("estimate_id", ""),
                "source_document_url": locator.get("source_document_url", ""),
                "source_document_kind": locator.get("source_document_kind", ""),
                "metric_n": source.get("metric_n", ""),
                "responders": source.get("responders", ""),
                "source_n": source.get("source_n", source.get("n_samples", "")),
                "unit": source.get("unit", "mutations_per_megabase" if is_tmb else ""),
                "value_basis": source.get("value_basis", ""),
                "trial_name": source.get("trial_name", ""),
                "trial_nct": source.get("trial_nct", ""),
                "source_population_label": source.get("source_population_label", ""),
                "setting": source.get("setting", ""),
                "tmb_assay": review.get("tmb_assay", ""),
                "source_review_notes": review.get("source_review_notes", ""),
                "notes": source.get("note", source.get("notes", "")),
            }
        )
    summary = json.loads((directory / "lookup-audit.json").read_text())
    summary.update(
        flagged_records=len(records),
        cited_reference_identifiers=len(references),
        distinct_canonical_references=len(
            {row["pmid"] or row["doi"].lower() or row["reference"] for row in references}
        ),
        tmb_references=sum(bool(row["tmb_record_count"]) for row in references),
        response_references=sum(bool(row["response_record_count"]) for row in references),
        uncited_records=sum(not row["reference"] for row in records),
        flags=dict(
            Counter(flag for row in records for flag in row["findings"].split(";") if flag in FLAGS)
        ),
        metadata_statuses=dict(Counter(row["metadata_status"] for row in references)),
        references_with_pmc=sum(bool(row["pmcid"]) for row in references),
        source_files={
            name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in SOURCE_FILES
        },
        snapshot_inputs={
            name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
            for name in (
                "bibliographic-metadata.json",
                "source-traces.csv",
                "additional-primary-reference.csv",
            )
        },
        additional_primary_references=len(additional),
        source_trace_records=len(traces),
        modeled_blends_without_direct_citation=sum(
            not row["reference"] and row["value_basis"] == "derived_blend" for row in records
        ),
    )
    return {
        "references.csv": csv_text(references),
        "flagged-records.csv": csv_text(records),
        "uncited-records.csv": csv_text(
            [row for row in records if not row["reference"]], list(records[0])
        ),
        "lookup-audit.json": json_text(summary),
        "bibliography.html": bibliography_html(references, records, traces, summary),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--directory", type=Path, default=DEFAULT_DIRECTORY)
    parser.add_argument(
        "--output", type=Path, help="Write elsewhere, preserving the checked-in export"
    )
    parser.add_argument("--check", action="store_true", help="Fail if generated files are stale")
    args = parser.parse_args()
    if args.check and args.output:
        parser.error("--check and --output cannot be used together")
    artifacts = build_inventory(args.root, args.directory)
    if args.check:
        stale = [
            name
            for name, content in artifacts.items()
            if not (args.directory / name).is_file()
            or (args.directory / name).read_bytes() != content.encode("utf-8")
        ]
        if stale:
            parser.exit(1, f"Stale reference-inventory artifacts: {', '.join(stale)}\n")
        print(
            "Reference inventory reproduces from the checked-in snapshot and current audit tables."
        )
    else:
        output = args.output or args.directory
        output.mkdir(parents=True, exist_ok=True)
        for name, content in artifacts.items():
            (output / name).write_text(content, encoding="utf-8")
        print(f"Wrote {len(artifacts)} reference-inventory artifacts to {output}")


if __name__ == "__main__":
    main()
