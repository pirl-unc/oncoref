#!/usr/bin/env python3
"""Write the Bradley nomination audit from oncoref's HPA v23 gates and inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from oncoref import hpa, reference_data
from oncoref.cta_gate_audit import bradley_candidate_audit
from oncoref.cta_tissues import CORE_REPRODUCTIVE_TISSUES


def save(fig, path):
    with matplotlib.rc_context({"pdf.fonttype": 42}):
        fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    audit = bradley_candidate_audit()
    audit.to_csv(out / "bradley-gate-audit.csv", index=False)
    ids = set(audit.Ensembl_Gene_ID)
    rna = hpa.hpa_rna_consensus()
    rna = rna[rna.Gene.isin(ids)].copy()
    ihc = hpa.hpa_normal_tissue()
    ihc = ihc[ihc.Gene.isin(ids)].copy()
    rna.to_csv(out / "bradley-normal-rna.csv", index=False)
    ihc.to_csv(out / "bradley-normal-ihc.csv", index=False)
    inputs = {}
    for name in ["hpa_rna_consensus", "hpa_normal_tissue"]:
        path = reference_data.local_path(name)
        inputs[name] = {"release": "v23", "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    (out / "source-audit.json").write_text(
        json.dumps(
            {
                "inputs": inputs,
                "paper_doi": "10.1038/s41467-020-19141-w",
                "n_candidates": len(audit),
                "default_retained": int(audit.default_panel.sum()),
                "cohort_context": "Normal-tissue audit; not cancer IHC or patient eligibility",
            },
            indent=2,
        )
        + "\n"
    )
    y = np.arange(len(audit))
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(
        y,
        audit.rna_deflated_reproductive_frac * 100,
        color=np.where(audit.rna_gate_pass, "#267b53", "#748ea5"),
        height=0.6,
        label="observed RNA fraction",
    )
    ax.scatter(
        audit.required_rna_fraction * 100,
        y,
        marker="|",
        s=230,
        color="#a63b32",
        linewidths=2,
        label="required for protein tier",
    )
    for i, row in audit.iterrows():
        ax.text(
            102,
            i,
            f"{row.rna_deflated_reproductive_frac:.2%} / {row.required_rna_fraction:.0%}",
            va="center",
            fontsize=9,
        )
    ax.set_yticks(y, audit.Symbol)
    ax.invert_yaxis()
    ax.set_xlim(0, 127)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xlabel("Deflated RNA fraction in testis, ovary and placenta (%)")
    ax.set_title("Bradley cancer-placenta candidates: the RNA restriction gate")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, out / "bradley-rna-gates")
    columns = ["coding_gate_pass", "protein_gate_pass", "rna_gate_pass", "default_panel"]
    values = audit[columns].astype(int).to_numpy()
    fig, ax = plt.subplots(figsize=(8.5, 6))
    from matplotlib.colors import ListedColormap

    ax.imshow(values, cmap=ListedColormap(["#eee5df", "#83b79a"]), vmin=0, vmax=1, aspect="auto")
    for (i, j), val in np.ndenumerate(values):
        ax.text(j, i, "pass" if val else "fail", ha="center", va="center", fontsize=10)
    ax.set_xticks(
        range(4), ["Protein-coding", "IHC restriction*", "RNA threshold", "Default panel"]
    )
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.set_yticks(y, audit.Symbol)
    ax.spines[:].set_visible(False)
    ax.set_title("Independent gates explain all ten outcomes", pad=35)
    fig.text(
        0.5,
        0.015,
        "* No IHC evidence is not an IHC veto; it requires the stricter 97% RNA gate.\nAn IHC pass is a normal-tissue rule, not proof of tumor antigenicity or clinical safety.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.075, 1, 1))
    save(fig, out / "bradley-gate-matrix")
    vgll = rna[rna.Gene.isin(audit.loc[audit.Symbol.eq("VGLL1"), "Ensembl_Gene_ID"])].copy()
    vgll = vgll[vgll.Tissue.ne("thymus")].sort_values("nTPM", ascending=False).head(12)
    fig, ax = plt.subplots(figsize=(9.2, 6))
    colors = ["#267b53" if t in CORE_REPRODUCTIVE_TISSUES else "#748ea5" for t in vgll.Tissue]
    ax.barh(np.arange(len(vgll)), vgll.nTPM, color=colors, height=0.65)
    for i, n in enumerate(vgll.nTPM):
        ax.text(n + 2, i, f"{n:g}", va="center", fontsize=9)
    ax.set_yticks(np.arange(len(vgll)), vgll.Tissue)
    ax.invert_yaxis()
    ax.set_xlim(0, vgll.nTPM.max() * 1.13)
    ax.set_xlabel("HPA v23 normal-tissue RNA (nTPM)")
    ax.set_title(
        "VGLL1: strongest normal-tissue RNA signals\nPlacental enrichment coexists with expression outside the core tissues"
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, out / "vgll1-normal-rna")
    return audit


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("docs/audits/cta-bradley"))
    plot(parser.parse_args().out)
