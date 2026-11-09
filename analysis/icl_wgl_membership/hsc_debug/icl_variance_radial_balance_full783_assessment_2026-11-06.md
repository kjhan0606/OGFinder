# Full-sample radial coadd-variance balance (2026-11-06)

## Scope and method

The 94-pair supported analysis was extended to every existing same-filter
target/control coadd pair: 302 `i` and 481 `i2` pairs (783 total). Astrafek
verified saved image hashes, IMAGE/MASK/VARIANCE HDU mappings, filter labels,
and center WCS before reading the delivered variance and mask planes. It
measured median variance in ten proper-radius bins from 0 to 900 kpc. Pixels
with non-finite or non-positive variance or reviewed bad/bright-star mask bits
were excluded. Bins with fewer than 50 usable pixels in either member of a
pair were omitted for that pair. Median variance and valid-pixel fractions are
recorded per pair; these are not full annulus coverage fractions.

The initial paired contrasts used 5,000 bootstrap resamples of connected
5-degree sky blocks, separately by filter. To account for the 20 correlated
filter/radius contrasts, the follow-up uses a single joint resampling of the
connected blocks across both filters, 10,000 replicates, and centered
max-absolute-t simultaneous 95% intervals. The joint graph contains 38
approximate blocks; its largest block has 107 of the 783 pairs.

## Findings

The median target/control variance ratio is near one in `i`. In `i2`, it is
about 0.96 from the center through 300 kpc and approaches one at larger radii.
After familywise adjustment across all 20 radial contrasts, three `i2`
intervals still lie below equal variance:

| Filter | Radius (kpc) | Median variance ratio | Simultaneous 95% interval |
|---|---:|---:|---:|
| `i2` | 75–100 | 0.961 | 0.924–0.999 |
| `i2` | 100–150 | 0.955 | 0.915–0.996 |
| `i2` | 200–300 | 0.961 | 0.934–0.988 |

The other 17 simultaneous intervals include one. This identifies a small
filter-dependent difference in the *delivered variance planes* for this
selected target/control sample. It does not establish that HSC targets have
lower total noise: the approximate block partition has only 38 components and
one large component, the target/control mask coverage differs, and the
variance plane omits the unknown effective-gain term for source Poisson noise.
All ICL-versus-variance partial-correlation intervals within the 94 supported
subset include zero. No variance correction was applied to the ICL profile.

An earlier familywise output directory (`icl_variance_radial_balance_familywise_2026-11-06_01`)
was invalidated: it pooled both filters into each labelled metric. The
corrected, filter-separated result is in `_02`; its summary records the
invalidated predecessor.

## Limits

The max-t intervals adjust only these 20 variance contrasts. Connected 5°
equal-area cells are approximate spatial blocks, not validated HSC fields.
This diagnostic does not model the covariance of the masked ICL estimator,
source-Poisson noise, PSF wings, sky transfer, or astrophysical covariance.
Keep `science_measurement=false`.

## Reproduction and artifacts

- Main runner: [`analyze_icl_variance_radial_balance.py`](analyze_icl_variance_radial_balance.py), SHA-256 `446f34b327f47ed04dbc658911b1a1306913515265f6bea59a9ab316abdcf9d5`.
- Familywise runner: [`analyze_icl_variance_radial_familywise.py`](analyze_icl_variance_radial_familywise.py), SHA-256 `e39e2ce8cc6bed3af0f2aa21832461202b3885bd2bcabfbc6bd3bc277777d7dd`.
- Full pair results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_radial_balance_full783_2026-11-06_01/` (`summary.json` SHA-256 `f155075dfc5bfdc978fdcff9d3252d369ee81d5ef1af23ae0371aef1d4d04632`).
- Corrected simultaneous intervals: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_radial_balance_familywise_2026-11-06_02/` (`summary.json` SHA-256 `bee18c612211a672059108ee6d59f537dd4ce27638f4d9eccaef5c26f677ccfb`; TSV SHA-256 `436286d95f341629ad7d6be17382abaf8fea780d7709c08c300f91ac9e332fe8`).

All 1,566 coadd file hashes matched their saved acquisition records. The
analysis read only already downloaded FITS data; no HSC query or download was
made.
