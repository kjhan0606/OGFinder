"""Arc-length-aware dynamical classification (pure numpy, no ephemeris)."""
import numpy as np
from moving import kepler as K


def test_short_arc_is_unclassified_even_with_a_confident_looking_fit():
    r = K.classify_arc(0.03, 2.7, 0.1, 5.0, sig=dict(a=0.01, e=0.01, i=0.1), determined=True)
    assert r["status"] == "unclassified" and r["label"].startswith("unclassified")
    assert "0.03" in r["label"]
    assert r["probs"]                                  # uncertainty is still reported


def test_short_arc_uses_ranging_probabilities():
    rp = {"MBA": 0.6, "NEO": 0.3, "TNO": 0.1}
    r = K.classify_arc(0.05, 2.7, 0.1, 5.0, determined=False, ranging_probs=rp)
    assert r["status"] == "unclassified" and r["most_probable"] == "MBA"
    assert abs(sum(r["probs"].values()) - 1.0) < 1e-9 and abs(r["p_most_probable"] - 0.6) < 1e-9


def test_degenerate_hyperbolic_fit_is_unclassified_on_long_arc():
    r = K.classify_arc(30.0, -5.0, 1.8, 5.0, determined=False)
    assert r["status"] == "unclassified" and "degenerate" in r["label"]
    r2 = K.classify_arc(30.0, 2.7, 1.3, 5.0, determined=True)
    assert r2["status"] == "unclassified"                    # e >= 1 is never a bound class


def test_long_precise_arc_is_classified():
    r = K.classify_arc(400.0, 2.7, 0.1, 5.0, sig=dict(a=0.001, e=0.001, i=0.01), determined=True)
    assert r["status"] == "classified" and r["label"] == "MBA" and r["p_most_probable"] > 0.99


def test_long_arc_straddling_a_boundary_is_ambiguous_with_probabilities():
    # a = 3.3 AU is the MBA / Cybele boundary; sigma_a = 0.2 -> ~50/50
    r = K.classify_arc(60.0, 3.3, 0.1, 5.0, sig=dict(a=0.2, e=0.01, i=0.5), determined=True)
    assert r["status"] == "ambiguous" and r["label"].startswith("ambiguous")
    assert 0.3 < r["probs"].get("MBA", 0) < 0.7 and r["probs"].get("Hilda/Cybele/Jupiter-region", 0) > 0.3


def test_non_psd_covariance_is_handled():
    cov = np.array([[1e-4, 0, 0], [0, -1e-6, 0], [0, 0, 1e-2]])      # negative eigenvalue (degenerate fit)
    p = K.class_probs_from_covariance(2.7, 0.1, 5.0, cov)
    assert p and abs(sum(p.values()) - 1.0) < 1e-6
    assert K.class_probs_from_covariance(np.nan, 0.1, 5.0, cov) is None
