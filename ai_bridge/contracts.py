"""Task contracts: the standard input record and the standard output columns per task type.

This file is the single source of truth for the docs table (docs/ai_services.md), for
`validate-profile`, for the TSV the CLI writes, and for the deterministic `mock` service.

Input record (one per catalog object; every adapter sees the same dict, see `INPUT_RECORD`):
    id, x, y, ra, dec, mags{band}, mag_errs{band}, wcs{...}, pixel_scale_arcsec,
    psf_fwhm_arcsec, cutouts{band: {...}} (lazy, produced on demand), catalog_row{...}

Output: named columns (below) + provenance columns  AI_<P>_SERVICE / _MODEL / _REQID / _TIME / _ERROR
where <P> is the task's provenance prefix.  Column types: float | int | str | bool | json.
"""

INPUT_RECORD = [
    ("id", "str", "catalog object identifier (NUMBER column by default)"),
    ("x, y", "float", "pixel position in the first image (1-based, SExtractor convention: X_IMAGE, Y_IMAGE)"),
    ("ra, dec", "float", "degrees, ICRS/J2000 (ALPHA_J2000, DELTA_J2000) or null"),
    ("mags", "dict band->float", "from MAG_<band> columns (MAG_AUTO -> band 'AUTO'); null when 99/blank"),
    ("mag_errs", "dict band->float", "from MAGERR_<band> columns"),
    ("cutouts", "dict band->file", "per-band cutout, FITS / PNG / npy, produced lazily (see cutouts.py)"),
    ("wcs", "dict", "FITS WCS keywords of the cutout (CTYPE, CRVAL, CRPIX, CD/CDELT, ...)"),
    ("pixel_scale_arcsec", "float", "arcsec/pixel (WCS, or --pixel-scale)"),
    ("psf_fwhm_arcsec", "float", "from --psf-fwhm, else null"),
    ("catalog_row", "dict col->str", "the full catalog row, as text"),
]

# provenance columns (values come from the service response, never invented)
PROVENANCE = ["SERVICE", "MODEL", "REQID", "TIME", "ERROR"]

# mock generators: ("float", lo, hi) | ("choice", [..]) | ("str", prefix)
TASKS = {
    "photoz": {
        "prefix": "PHOTOZ", "needs": ["mags"],
        "doc": "photometric redshift from catalogue magnitudes (and/or cutouts)",
        "outputs": [
            ("PHOTOZ", "float", "", "best / point-estimate redshift", ("float", 0.0, 3.0)),
            ("PHOTOZ_ERR", "float", "", "1-sigma uncertainty", ("float", 0.01, 0.3)),
            ("PHOTOZ_P16", "float", "", "16th percentile of the PDF", ("float", 0.0, 3.0)),
            ("PHOTOZ_P50", "float", "", "median of the PDF", ("float", 0.0, 3.0)),
            ("PHOTOZ_P84", "float", "", "84th percentile of the PDF", ("float", 0.0, 3.0)),
        ]},
    "sed_fit": {
        "prefix": "SED", "needs": ["mags"],
        "doc": "stellar-population / SED fit; names match the local SED-fit columns",
        "outputs": [
            ("LOG_MASS", "float", "log10(Msun)", "stellar mass", ("float", 7.0, 12.0)),
            ("LOG_MASS_ERR", "float", "dex", "mass uncertainty", ("float", 0.05, 0.4)),
            ("LOG_AGE", "float", "log10(yr)", "mass-weighted age", ("float", 7.5, 10.1)),
            ("LOG_AGE_ERR", "float", "dex", "age uncertainty", ("float", 0.05, 0.4)),
            ("LOG_Z", "float", "log10(Z/Zsun)", "metallicity", ("float", -2.0, 0.4)),
            ("AV", "float", "mag", "V-band dust attenuation", ("float", 0.0, 2.0)),
            ("SFR", "float", "Msun/yr", "star-formation rate", ("float", 0.0, 50.0)),
            ("SED_CHI2", "float", "", "reduced chi-square of the fit", ("float", 0.5, 3.0)),
        ]},
    "morphology": {
        "prefix": "MORPH", "needs": ["cutouts"],
        "doc": "galaxy morphology from cutouts; names match the local morphology columns",
        "outputs": [
            ("MORPH_TYPE", "str", "", "morphological class label", ("choice", ["E", "S0", "Sa", "Sb", "Sc", "Irr"])),
            ("MORPH_CONF", "float", "", "confidence / probability of MORPH_TYPE, 0..1", ("float", 0.3, 1.0)),
            ("MORPH_DESC", "str", "", "free-text description (optional)", ("str", "mock morphology ")),
        ]},
    "star_galaxy": {
        "prefix": "SG", "needs": ["cutouts|mags"],
        "doc": "star / galaxy separation",
        "outputs": [
            ("STAR_PROB", "float", "", "probability that the object is a star, 0..1", ("float", 0.0, 1.0)),
            ("CLASS_LABEL", "str", "", "STAR | GALAXY | (QSO ...) as returned by the service", ("choice", ["STAR", "GALAXY"])),
        ]},
    "real_bogus": {
        "prefix": "RB", "needs": ["cutouts"],
        "doc": "transient real/bogus score on (science, template, difference) cutouts",
        "outputs": [
            ("REALBOGUS_SCORE", "float", "", "1 = real, 0 = bogus (service convention noted in the profile)", ("float", 0.0, 1.0)),
            ("CLASS_LABEL", "str", "", "REAL | BOGUS (optional)", ("choice", ["REAL", "BOGUS"])),
        ]},
    "transient_classification": {
        "prefix": "TRANS", "needs": ["mags|cutouts"],
        "doc": "transient type classification",
        "outputs": [
            ("CLASS_LABEL", "str", "", "most probable class", ("choice", ["SN Ia", "SN II", "AGN", "other"])),
            ("CLASS_PROB", "float", "", "probability of CLASS_LABEL", ("float", 0.2, 1.0)),
            ("CLASS_PROBS", "json", "", "all class probabilities as a JSON object (optional)", ("str", "{}")),
        ]},
    "moving_object_classification": {
        "prefix": "MOV", "needs": ["cutouts|catalog_row"],
        "doc": "moving-object / asteroid / artifact classification",
        "outputs": [
            ("CLASS_LABEL", "str", "", "class label", ("choice", ["MOVING", "STATIONARY", "ARTIFACT"])),
            ("CLASS_PROB", "float", "", "probability of CLASS_LABEL", ("float", 0.2, 1.0)),
            ("MOVING_SCORE", "float", "", "moving-object score 0..1 (optional)", ("float", 0.0, 1.0)),
        ]},
    "anomaly_detection": {
        "prefix": "ANOM", "needs": ["cutouts|mags"],
        "doc": "anomaly / outlier score",
        "outputs": [
            ("ANOMALY_SCORE", "float", "", "larger = more anomalous (convention noted in the profile)", ("float", 0.0, 1.0)),
            ("ANOMALY_FLAG", "bool", "", "service-defined anomaly decision (optional)", ("choice", ["False", "True"])),
        ]},
    "embedding_similarity": {
        "prefix": "EMB", "needs": ["cutouts|mags"],
        "doc": "embedding / similarity search; the vector itself stays with the service, a reference is stored",
        "outputs": [
            ("EMBEDDING_ID", "str", "", "identifier of the stored embedding", ("str", "mock-emb-")),
            ("SIMILAR_IDS", "str", "", "comma separated ids of the nearest neighbours (optional)", ("str", "mock-nn-")),
            ("SIMILARITY_TOP", "float", "", "similarity of the nearest neighbour (optional)", ("float", 0.0, 1.0)),
        ]},
    "captioning": {
        "prefix": "CAP", "needs": ["cutouts"],
        "doc": "image captioning / object description (free text)",
        "outputs": [
            ("CAPTION", "str", "", "text description of the cutout", ("str", "mock caption for object ")),
        ]},
    "image_retrieval": {
        "prefix": "RET", "needs": ["cutouts"],
        "doc": "image retrieval: ids of similar images in a remote archive",
        "outputs": [
            ("RETRIEVED_IDS", "str", "", "comma separated ids of retrieved images", ("str", "mock-img-")),
            ("RETRIEVAL_SCORE", "float", "", "score of the best hit (optional)", ("float", 0.0, 1.0)),
        ]},
    "generic": {
        "prefix": "GEN", "needs": [],
        "doc": "anything else: the response mapping defines the column names (also used by tap_query)",
        "outputs": [],
    },
}

TASK_NAMES = list(TASKS)


def prov_prefix(task, profile=None):
    """Provenance prefix: profile key "provenance_prefix" (A-Z0-9) overrides the task default."""
    p = (profile or {}).get("provenance_prefix")
    return p if p else TASKS[task]["prefix"]


def prov_columns(task, profile=None):
    p = prov_prefix(task, profile)
    return ["AI_%s_%s" % (p, k) for k in PROVENANCE]


def output_names(task):
    return [o[0] for o in TASKS[task]["outputs"]]


def output_spec(task, name):
    for o in TASKS[task]["outputs"]:
        if o[0] == name:
            return {"type": o[1], "units": o[2], "doc": o[3], "mock": o[4]}
    return None


def contract_table_markdown():
    lines = ["| task | provenance prefix | output columns (type) |", "|---|---|---|"]
    for t, d in TASKS.items():
        cols = ", ".join("`%s` (%s)" % (o[0], o[1]) for o in d["outputs"]) or "defined by the profile's response mapping"
        lines.append("| `%s` | `AI_%s_*` | %s |" % (t, d["prefix"], cols))
    return "\n".join(lines)

# columns whose values must lie in [0, 1]; anything else is rejected as a mapping error
UNIT_INTERVAL = {"STAR_PROB", "REALBOGUS_SCORE", "MORPH_CONF", "CLASS_PROB", "MOVING_SCORE", "ANOMALY_SCORE"}
