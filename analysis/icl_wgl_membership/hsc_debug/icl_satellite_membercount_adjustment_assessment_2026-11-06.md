# Contributing-member count sensitivity for the ICL–satellite relation

## Sample and method

The fixed common-support sample contains 94 matched cluster/filter pairs,
including 26 `i` pairs and 68 `i2` pairs. Each pair relates the target-minus-
control ICL surface brightness at 200–300 proper kpc to the equal-member
satellite excess at 50–75 proper kpc. The analysis compares a partial Pearson
correlation controlling redshift and linear CAMIRA richness with a correlation adding
the log number of detected members in the satellite stack. The member count
was recomputed from the contributing per-member profile rows for every pair.

Uncertainty comes from 10,000 requested joint resamples of connected 5-degree
target and control sky blocks. A total of 9,989 resamples produced finite
estimates in both filters. Both filters share each resampled block draw.
Reported 95% intervals are pointwise and exploratory. The analysis did not
test a new confirmatory hypothesis.

## Results

For `i2`, the baseline partial correlation is 0.293. Adding log member count
changes the estimate to 0.295. The paired adjustment is 0.002 with a 95%
block-bootstrap interval of −0.023 to 0.044. The partial correlations between
member count and ICL excess, and between member count and satellite excess,
are 0.057 (−0.185, 0.313) and −0.032 (−0.226, 0.136), respectively.

For `i`, the baseline correlation is 0.079 and the member-count-adjusted
correlation is 0.088. The paired change is 0.009 with a 95% interval of
−0.081 to 0.192. Member-count correlations with ICL excess and satellite
excess are −0.074 (−0.506, 0.267) and 0.113 (−0.389, 0.427). All intervals
for member-count associations and model changes include zero.

The measured member count does not explain the positive `i2` point estimate
within the selected 94-pair sample. The broad `i` intervals permit substantial
changes in either direction. The result does not correct the detected-satellite
selection, member incompleteness, source-Poisson variance, PSF wings, sky
transfer, or coadd covariance. Member count can be a precision proxy, mediator,
or collider, so the adjusted correlation is a sensitivity estimate rather than
a causal effect. Keep `science_measurement=false`.

## Reproducibility

- Astrafek runner: [`analyze_icl_satellite_membercount_adjustment.py`](analyze_icl_satellite_membercount_adjustment.py).
- Runner SHA-256: `13a6c6b6168a65e5dc0f71e40a963f0298c62c04a14266683aae3b466a39f16b`.
- Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_membercount_adjustment_2026-11-06_02/`.
- Summary SHA-256: `12420feadba8d83badf3c1cb8446652c4bf3c3041bcbc70eae84e2f94011b485`.
- Pair table SHA-256: `df9484f4da89de4e57cdcaeb2cc7a98fdc82e1d0a9dbd84cbb666397dae45ccf`.
- Filter summary SHA-256: `7c7cadb7cdbc54edaa355a486f3329dd96dfc44b62e52cf15eb9a0d392e8692f`.
- The run used existing local profile tables and coordinate manifests. No HSC query, download, or image-pixel read occurred.
