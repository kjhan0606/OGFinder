#!/usr/bin/env python3
"""Summarize verified independent-support ICL stacks by mask and filter."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_WORKS = (
    DEFAULT_ROOT / "icl_2d_independent_support_k2p5_2026-10-19_02",
    DEFAULT_ROOT / "icl_2d_independent_support_k3p5_2026-10-19_02",
    DEFAULT_ROOT / "icl_2d_independent_support_k5p0_2026-10-19_02",
)
DEFAULT_OUTPUT = HERE / "icl_2d_independent_support_assessment_2026-10-19.md"
RADIAL_FILE = "independent_support_radial_profiles.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def verified_profile(work: Path) -> tuple[dict, list[dict[str, str]]]:
    summary_path = work / "summary.json"
    verification_path = work / "independent_verification.json"
    radial_path = work / RADIAL_FILE
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    verification = json.loads(verification_path.read_text(encoding="utf-8"))
    if (summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False
            or verification.get("status") != "verified"
            or verification.get("science_measurement") is not False
            or verification.get("summary_sha256") != sha256(summary_path)):
        raise ValueError(f"run is not independently verified: {work}")
    if verification.get("n_pairs") != summary.get("n_pairs"):
        raise ValueError(f"verification pair count differs: {work}")
    if sha256(radial_path) != summary.get("outputs", {}).get("radial_profiles_sha256"):
        raise ValueError(f"radial-profile checksum differs: {work}")
    rows = read_tsv(radial_path)
    if len(rows) != 40:
        raise ValueError(f"expected 40 radial bins: {work}")
    return summary, rows


def value(row: dict[str, str], key: str) -> float | None:
    raw = row.get(key, "")
    return None if raw in ("", "None") else float(raw)


def fmt(x: float | None, digits: int = 4) -> str:
    return "—" if x is None or not math.isfinite(x) else f"{x:.{digits}g}"


def fmt_percent(x: float | None) -> str:
    return "—" if x is None or not math.isfinite(x) else f"{x:.3g}%"


def pct_change(new: float | None, old: float | None) -> float | None:
    if new is None or old is None or old == 0:
        return None
    return 100.0 * (new - old) / abs(old)


def render(works: list[Path], output: Path) -> str:
    loaded = [(work, *verified_profile(work)) for work in works]
    mask_labels = []
    for work, summary, _ in loaded:
        label = f"K={float(summary['mask_factor']):g}"
        if label in mask_labels:
            raise ValueError(f"duplicate mask factor: {label}")
        mask_labels.append(label)
    if set(mask_labels) != {"K=2.5", "K=3.5", "K=5"}:
        raise ValueError(f"expected K=2.5, 3.5, and 5.0; got {mask_labels}")

    lines = [
        "# Independent-support ICL estimator sensitivity — 2026-10-19",
        "",
        "## Input and method",
        "",
        "This assessment compares the independent-support estimator—target and matched-control means on each frame's own valid pixels—with the strict common-support estimator on identical target/control pairs and identical 25 kpc annuli. The three source runs use mask factors K=2.5, 3.5, and 5.0.",
        "",
        "Each input run passed `verify_icl_independent_support.py`. All outputs remain exploratory and keep `science_measurement=false`. The tabulated pair SEM is descriptive; it is not a spatial-block covariance interval.",
        "",
        "## 50–200 kpc radial bins",
        "",
        "The support-gain columns count pairs with a valid independent target-minus-control mean versus pairs that also have a valid common-support mean. Coverage columns are median geometric pixel fractions across available pairs.",
        "",
        "| Mask | Filter | R (kpc) | N independent | N common | Added | Independent Δ ± pair SEM | Independent Δ on common pairs | Common Δ on same pairs | Shift (%) | Median target / control / common coverage |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    shifts: dict[tuple[str, str], list[float]] = {}
    for work, summary, rows in loaded:
        mask = f"K={float(summary['mask_factor']):g}"
        for row in rows:
            lo, hi = int(row["r_min_kpc"]), int(row["r_max_kpc"])
            if lo < 50 or hi > 200:
                continue
            band = row["filter"]
            ind_common = value(row, "mean_independent_delta_on_common_pairs_njy_arcsec2")
            common = value(row, "mean_common_delta_on_same_pairs_njy_arcsec2")
            shift = pct_change(ind_common, common)
            if shift is not None:
                shifts.setdefault((mask, band), []).append(abs(shift))
            lines.append(
                f"| {mask} | {band} | {lo}–{hi} | {row['n_pairs_independent_delta']} | "
                f"{row['n_pairs_common_support']} | {row['n_pairs_added_vs_common']} | "
                f"{fmt(value(row, 'mean_delta_independent_njy_arcsec2'))} ± "
                f"{fmt(value(row, 'pair_sem_delta_independent_njy_arcsec2'))} | "
                f"{fmt(ind_common)} | {fmt(common)} | {fmt_percent(shift)} | "
                f"{fmt(value(row, 'median_target_coverage'), 3)} / "
                f"{fmt(value(row, 'median_control_coverage'), 3)} / "
                f"{fmt(value(row, 'median_common_coverage'), 3)} |"
            )

    lines.extend(["", "## Reading the comparison", ""])
    for mask, summary, _ in loaded:
        lines.append(
            f"- {mask}: {summary['n_pairs']} matched pairs "
            f"(i={summary['n_by_filter']['i']}, i2={summary['n_by_filter']['i2']}); "
            "all source-input paths exist and the saved annular aggregates reproduce from per-pair profiles."
        )
    for key in sorted(shifts):
        vals = shifts[key]
        lines.append(
            f"- {key[0]}, {key[1]}: maximum absolute same-pair estimator shift over 50–200 kpc is {max(vals):.2f}% across the six annuli."
        )
    lines.extend([
        "",
        "The independent estimator recovers more pairs or valid pixels only where the target and control have separate support. The shift can therefore reflect both extra area and differing azimuthal/source-mask environments. It is a sensitivity diagnostic, not a correction to the common-support stack.",
        "",
        "These runs do not resolve residual sky transfer, position-dependent extended PSF wings, extinction gradients, coadd noise covariance, or satellite-selection completeness. They do not measure ICL–satellite association significance. Keep the common-support estimator as the reference and retain `science_measurement=false` pending those gates and covariance-aware analysis.",
        "",
        "## Reproducibility",
        "",
        "- Astrafek runner: [`stack_hsc_icl_independent_support.py`](stack_hsc_icl_independent_support.py)",
        "- Independent verifier: [`verify_icl_independent_support.py`](verify_icl_independent_support.py)",
        "- Assessment generator: [`summarize_icl_independent_support.py`](summarize_icl_independent_support.py)",
    ])
    for work, summary, _ in loaded:
        lines.append(f"- {work.name}: `summary.json` SHA-256 `{sha256(work / 'summary.json')}`; verifier SHA-256 `{sha256(work / 'independent_verification.json')}`.")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, action="append",
                        help="verified run directory; may be repeated, defaults to all three K factors")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    works = args.work or list(DEFAULT_WORKS)
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite assessment: {output}")
    report = render(works, output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(report, encoding="utf-8")
    print(f"wrote {output}")
    print(f"sha256={sha256(output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
