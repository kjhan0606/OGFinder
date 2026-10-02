# Multiwavelength cross-match (plugin `xmatch`, tab Measure, chip "X-match")

Match the catalog with a reference table by sky position (or pixel position) and copy information over. Engine `ogfkit/xmatch.py` (pure functions), CLI `plugins/xmatch/xmatch.py`.
Reference sources: a local file (FITS table, CSV, TSV/ASCII, VOTable; RA/Dec in degrees or sexagesimal, hours for RA) or **one cone query** over the whole field through an IVOA TAP
service (the existing `ai_bridge` `tap_query` adapter, `<base>/sync`, ADQL template with `{ra} {dec} {radius_deg}`), matched locally. No network access happens unless
"Allow network access" is on (CLI `--allow-network`); the raw answer is kept in `<work>/xmatch/xmatch_tap.csv`.

## Output
Step **Cross-match** adds, per catalog object: `XM_SEP` (arcsec, or pixels in pixel mode), `XM_N` (candidates within the radius), `XM_FLAG`, `XM_DRA`/`XM_DDEC` (reference − object, arcsec, dRA includes cos Dec),
`XM_ID` (text of the ID column, exact, e.g. 19-digit Gaia ids), `XM_V1..XM_V4` (up to four copied columns of the nearest counterpart).
`XM_FLAG` bits: 1 = matched, 2 = ambiguous (>1 candidate inside the radius), 4 = mutual nearest (the object is also the nearest catalog object of its counterpart), 8 = the counterpart is claimed by several objects.
Files in `<work>/xmatch/`: `xmatch_summary.json` (counts, median offset and scatter of the clean pairs, chance-match estimate, TAP info), `xmatch_pairs.tsv`, `xmatch_tap.csv`.
Windows: *Match Summary…*, *Pair Table…*. Meta key `xmatch` in the catalog metadata.

## Systematic offset and chance matches
* The median RA/Dec offset (3σ-clipped, 5 iterations) of the matched pairs is always reported; "Remove the median offset" repeats the match after subtracting it from the reference.
* Chance matches are estimated by shifting the reference by 30–120″ in 8 directions and counting objects with a counterpart (edge-corrected for the footprint); `chance.mean` / `n_matched` is the false-match rate.

## Validation (tests in `plugins/xmatch/tests`, numbers measured)
* Synthetic, σ = 0.15″ per axis per catalog (pair separation Rayleigh 0.21″): completeness 0.62 / 0.98 / 1.00 at radius 0.3 / 0.6 / 1.0″ (analytic 0.98 at 0.6″), purity 1.00 with 300 unrelated reference sources.
* Chance-match estimate vs truth (12 random fields, 500 vs 800 sources, 3″): estimated 79.2, true 75.8 (+4.5 %; without the edge correction it was 18 % low).
* Systematic shift (+0.45″, −0.30″) recovered as +0.436″, −0.284″ (381 pairs); matching at 0.6″ goes from 61.8 % to 99.8 % correct when the shift is removed.
* CLI file mode 200/200 ids correct with an offset of 0.3″/0.2″ recovered within 0.05″; pixel mode, FITS/CSV/ECSV/VOTable/TSV readers and sexagesimal input tested; fake local TAP server: 100/100 ids correct, request is a POST to `/tap/sync` with LANG=ADQL.
* **Live**, real Gaia DR3 (ESA TAP, M51 field): 294 real sources perturbed by 0.1″ noise and +0.35″/−0.25″ offset, queried and matched through the TAP path: 284/294 correct ids (the rest are close pairs in crowded
  spots at the 2″ radius), offset recovered +0.352″/−0.241″, chance matches 14.5 of 294 (909 reference rows). Skipped when offline.
* Live on the GUI M51 catalog (206 objects) vs Gaia DR3 at 2″: 105 matched, 100 unique, 5 ambiguous, chance 1.2 (1.1 %), median offset (+0.38″, −0.58″) with scatter 0.3″ — the offset is a property of that image's WCS, not a defect.

## Limits
* Positions only; no proper-motion/epoch propagation, no magnitude/colour priors or likelihood-ratio matching, nearest-neighbour resolution (not a global one-to-one assignment: flag 8 marks conflicts).
* The TAP path sends one cone query for the field (radius from the catalog extent); services with row limits truncate large fields (`TOP 50000` in the default ADQL) — use your own ADQL for other tables/surveys. The answer is not cached between runs, so a replay of an exported session queries again.
* Sky mode needs ALPHA_J2000/DELTA_J2000 in the catalog (or set the column names); only simple tangent-plane flat separations are used for the shift (fine below ~1°).
* TAP columns are numeric when all values parse, otherwise text; copied values that are integers are written as text.

## CLI
```
python3 plugins/xmatch/xmatch.py --catalog cat.tsv --work out --source file --ref-file ref.fits --radius-arcsec 1.0 --ref-id-col objid --copy-columns gmag,plx [--apply-shift] [--mode pixel --ref-x-col X --ref-y-col Y]
python3 plugins/xmatch/xmatch.py --catalog cat.tsv --work out --source tap --allow-network --tap-url https://gea.esac.esa.int/tap-server/tap --ref-id-col source_id --copy-columns phot_g_mean_mag,parallax
```
