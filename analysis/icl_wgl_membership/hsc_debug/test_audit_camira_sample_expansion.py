import pytest

from audit_camira_sample_expansion import evaluate, greedy_additions, sky_cell


def row(ident, ra, dec, z=0.4, richness=20.0, z_bcg=0.4):
    return {"id": str(ident), "name": f"cluster-{ident}", "ra_deg": ra,
            "dec_deg": dec, "z": z, "richness": richness, "z_bcg": z_bcg}


def test_greedy_pool_respects_spherical_separation_and_current_ids():
    current = [row(i, (130 + i * 0.5) % 360, 0) for i in range(800)]
    candidates = [row("near", 130, 0, richness=100),
                  row("new", 350, 20, richness=50),
                  row("duplicate", 0, 0, richness=25)]
    additions = greedy_additions(current, candidates, 0.1)
    assert [r["id"] for r in additions] == ["new"]


def test_expansion_scenarios_separate_strict_and_relaxed_redshift_rules():
    current = [row(i, (i * 0.5) % 360, (i % 120) - 60) for i in range(800)]
    catalog = current + [
        row("spec_ok", 140, 20, z=.6, z_bcg=.6),
        row("photoz_only", 200, 20, z=.6, z_bcg=-1),
        row("bad_spec", 260, 20, z=.4, z_bcg=.45),
    ]
    result = evaluate(catalog, current)
    assert result["selection_scenarios"]["broader_z_with_bcg_spec"]["greedy_additions_after_6arcmin_separation"] == 1
    assert result["selection_scenarios"]["broader_z_cluster_redshift"]["greedy_additions_after_6arcmin_separation"] == 3


def test_sky_cell_rejects_invalid_coordinates():
    with pytest.raises(ValueError):
        sky_cell(row("bad", 360, 0))
