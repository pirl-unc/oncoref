"""Parse reported staining into categorical observations, separate from confidence.

The parser records what a source says. It does not establish reagent specificity,
normal-tissue coverage or eligibility to override another assay. Original text
and numeric ranges remain available for review; absent numbers stay absent.
"""

from __future__ import annotations

import re
from typing import Any


def parse_staining_statement(statement: str) -> dict[str, Any]:
    """Extract result, intensity, localization and explicitly reported percentages.

    Mixed descriptions retain both positive and negative findings. Cross-reactive
    staining is a source claim, never silently translated into protein absence.
    A recommendation about a control is not an observed staining measurement.
    """
    text = " ".join(statement.split())
    lower = text.lower()
    control = bool(re.search(r"\b(should|must|expected to)\b", lower))
    cross = bool(re.search(r"cross[- ]?react|non[- ]?specific|background staining", lower))
    negative = bool(
        re.search(
            r"\bnegative\b|not (?:detect|observ|seen|stain)|no (?:\w+ ){0,3}(?:staining|expression)|"
            r"(?:staining|expression|positivity) (?:is |was |are )?(?:largely )?absent|"
            r"absence of|(?:staining|expression) is (?:usually )?lacking",
            lower,
        )
    )
    # Remove explicitly negated staining clauses before detecting positive words.
    positive_text = re.sub(
        r"(?:not|no) (?:\w+ ){0,3}(?:detected|observed|seen|staining|expression)", "", lower
    )
    positive = bool(
        re.search(
            r"\bpositive\b|\bpositivity\b|\b(?:weak|faint|moderate|strong|intense|intensive)\b|"
            r"(?:staining|expression) (?:is |was )?(?:seen|observed|present|found|occurs)|"
            r"(?:cells?|fibres?|fibers?) (?:\w+ ){0,2}stain|"
            r"(?:staining|immunostaining) (?:of|occurs|in)|"
            r"(?:cytoplasmic|nuclear|membranous|distinct) (?:\w+ ){0,2}staining|"
            r"expression is high|(?:is |are )?(?:found|expressed|occurs) in",
            positive_text,
        )
    )
    if control:
        result = "not_assessed"
    elif cross and not positive:
        result = "reported_cross_reactivity"
    elif cross:
        result = (
            "reported_cross_reactivity"
            if not re.search(r"but|while|whereas|however|except", lower)
            else "mixed"
        )
    elif negative and positive:
        result = "mixed"
    elif negative:
        result = "not_detected"
    elif positive:
        result = "detected"
    else:
        result = "not_assessed"
    intensity = []
    for name, pattern in [
        ("low", r"weak|faint"),
        ("moderate", r"moderate"),
        ("high", r"strong|intense|intensive"),
    ]:
        if re.search(pattern, lower):
            intensity.append(name)
    localization = [
        name
        for name, pattern in [
            ("nuclear", r"\bnuclear\b"),
            ("cytoplasmic", r"cytoplasm"),
            ("membranous", r"membran"),
            ("extracellular", r"extracellular|fibres|fibers"),
        ]
        if re.search(pattern, lower)
    ]
    ranges = []
    for match in re.finditer(
        r"(?P<operator><=|>=|<|>|≤|≥|~)?\s*(?P<low>\d+(?:\.\d+)?)(?:\s*[-–]\s*(?P<high>\d+(?:\.\d+)?))?\s*%",
        text,
    ):
        ranges.append(
            {
                "operator": match["operator"] or ("range" if match["high"] else "reported"),
                "low": float(match["low"]),
                "high": float(match["high"] or match["low"]),
            }
        )
    return {
        "result": result,
        "intensity": ";".join(intensity) or ("none" if result == "not_detected" else "unspecified"),
        "localization": ";".join(localization) or "unspecified",
        "reported_cross_reactivity": cross,
        "is_control_recommendation": control,
        "reported_percentages": ranges,
        "source_statement": text,
    }
