#!/usr/bin/env python3
"""Compare pooled and within-field ICL–satellite slopes on the same pairs."""

import argparse
import json
from pathlib import Path

import numpy as np

from analyze_s23b_icl_satellite_relation import BLOCK_WIDTHS, GRID_PHASES
from audit_s23b_icl_satellite_field_omission import (
    DATA, FLOORS, N_BOOTSTRAP, ORIGIN, RELATION, RELATION_SHA,
    connected_pair_blocks, model_matrix, read_tsv, sha256, write_tsv,
)

FIELD_AUDIT = DATA / "s23b_icl_satellite_field_omission_2026-11-18_01"
FIELD_AUDIT_SHA = "c06deb42fcd060b262df31f968cca3ebf2deef7db1226f5dc09f47ee18cba0fd"


def load_pairs():
    inputs = {
        RELATION / "summary.json": RELATION_SHA,
        FIELD_AUDIT / "summary.json": FIELD_AUDIT_SHA,
    }
    for path, digest in list(inputs.items()):
        if sha256(path) != digest:
            raise ValueError(f"Changed input summary: {path}")
        name = ("paired_icl_satellite_profiles.tsv" if path.parent == RELATION
                else "pair_field_membership.tsv")
        inputs[path.parent / name] = json.loads(path.read_text())["output_sha256"][name]
    for path, digest in inputs.items():
        if sha256(path) != digest:
            raise ValueError(f"Changed input: {path}")
    ledger = read_tsv(FIELD_AUDIT / "pair_field_membership.tsv")
    keys = ("target_id", "control_id", "filter")
    fields = {tuple(row[k] for k in keys): row["target_field"] for row in ledger}
    if len(fields) != len(ledger) or any(r["target_field"] != r["control_field"] for r in ledger):
        raise ValueError("Field ledger has duplicate or cross-field pairs")
    pairs = [r for r in read_tsv(RELATION / "paired_icl_satellite_profiles.tsv")
             if r["icl_relative_contrast_njy_arcsec2"]]
    pair_keys = [tuple(row[k] for k in keys) for row in pairs]
    if len(pairs) != 988 or len(set(pair_keys)) != len(pairs) or set(pair_keys) != fields.keys():
        raise ValueError("Field membership differs from fitted pairs")
    for row, key in zip(pairs, pair_keys, strict=True):
        row["field"] = fields[key]
    return pairs, {str(path): digest for path, digest in inputs.items()}


def fitted_slopes(design, outcome, field_columns, weights):
    use = weights > 0
    roots = np.sqrt(weights[use])
    field_columns = field_columns[use]
    # Empty fields in a resample contribute no intercept or information.
    present = field_columns.any(axis=0)
    pooled = design[use] * roots[:, None]
    within = np.column_stack((pooled[:, 1:], field_columns[:, present] * roots[:, None]))
    response = outcome[use] * roots
    slopes = []
    for matrix, index in ((pooled, 1), (within, 0)):
        coef, _, rank, singular = np.linalg.lstsq(matrix, response, rcond=None)
        if rank != matrix.shape[1] or singular[-1] <= singular[0] * 1e-5:
            return None
        slopes.append(float(coef[index]))
    return slopes


def paired_intervals(design, outcome, field_columns, labels, seed):
    groups, row_group = np.unique(labels, return_inverse=True)
    draws = np.random.default_rng(seed).multinomial(
        len(groups), np.full(len(groups), 1 / len(groups)), size=N_BOOTSTRAP,
    )
    valid = [slopes for draw in draws
             if (slopes := fitted_slopes(design, outcome, field_columns, draw[row_group])) is not None]
    if len(valid) < 0.9 * N_BOOTSTRAP:
        raise ValueError(f"Only {len(valid)} valid matched bootstrap draws")
    slopes = np.asarray(valid)
    contrasts = np.column_stack((slopes, slopes[:, 1] - slopes[:, 0]))
    return np.percentile(contrasts, [2.5, 97.5], axis=0).T, len(valid)


def run(output):
    if output.exists():
        raise FileExistsError(output)
    pairs, input_digests = load_pairs()
    results = []
    for band_index, band in enumerate(("i", "i2")):
        for floor_index, floor in enumerate(FLOORS):
            rows = [r for r in pairs if r["filter"] == band and min(
                float(r["inner_common_support_fraction"]),
                float(r["outer_common_support_fraction"]),
            ) >= floor]
            fields, membership = np.unique([r["field"] for r in rows], return_inverse=True)
            indicators = membership[:, None] == np.arange(len(fields))
            controls = {r["target_id"]: dict(ra_deg=r["control_ra_deg"],
                                            dec_deg=r["control_dec_deg"]) for r in rows}
            design, outcome, _ = model_matrix(rows)
            pooled, within = fitted_slopes(design, outcome, indicators, np.ones(len(rows)))
            # Verify the categorical-intercept slope independently by partialling out field means.
            joined = np.column_stack((design[:, 1:], outcome))
            centered = joined.copy()
            for mask in indicators.T:
                centered[mask] -= joined[mask].mean(axis=0)
            check = np.linalg.lstsq(centered[:, :-1], centered[:, -1], rcond=None)[0][0]
            if not np.isclose(within, check, atol=1e-10, rtol=0):
                raise ValueError("Within-field point estimate differs from field demeaning")
            for grid_index, (width, phase) in enumerate(
                (width, phase) for width in BLOCK_WIDTHS for phase in GRID_PHASES
            ):
                labels, blocks = connected_pair_blocks(rows, controls, width, phase, ORIGIN)
                seed = 2026111810 + 100 * band_index + 10 * floor_index + grid_index
                intervals, n_valid = paired_intervals(design, outcome, indicators, labels, seed)
                result = dict(filter=band, support_floor=floor, width_deg=width, phase=phase,
                              n_pairs=len(rows), n_fields=len(fields), n_blocks=len(blocks),
                              valid_draws=n_valid, seed=seed)
                for key, point, (lower, upper) in zip(
                    ("pooled", "within_field", "difference"), (pooled, within, within-pooled),
                    intervals, strict=True,
                ):
                    result.update({key: point, f"{key}_low": lower, f"{key}_high": upper})
                results.append(result)
    output.mkdir(parents=True, exist_ok=False)
    table = output / "within_field_sensitivity.tsv"
    write_tsv(table, results)
    summary = dict(
        science_measurement=False, inputs_sha256=input_digests, n_fits=len(results),
        draws_per_layout=N_BOOTSTRAP, runner_sha256=sha256(Path(__file__)),
        output_sha256={table.name: sha256(table)}, slope_unit="nJy arcsec^-2 dex^-1",
        method="Same covariates and pairs. Compare a shared intercept with one intercept per "
               "observed field. Matched connected-block resamples refit both models. Differences "
               "are within-field minus pooled slopes.",
        limitation="Exploratory conditional slopes and percentile intervals on saved profiles. "
                   "Field adjustment does not correct cluster selection, within-field systematics, "
                   "or provide a survey covariance.",
    )
    with (output / "summary.json").open("x") as stream:
        json.dump(summary, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output.resolve())
