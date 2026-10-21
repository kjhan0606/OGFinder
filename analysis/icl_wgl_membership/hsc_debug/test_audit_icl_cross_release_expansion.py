from audit_icl_cross_release_expansion import evaluate


def candidate(ident, ra, richness, z=0.4, z_bcg=0.4):
    return {
        "id": str(ident), "name": str(ident), "ra_deg": ra, "dec_deg": 0.0,
        "z": z, "richness": richness, "z_bcg": z_bcg,
    }


def power_plan():
    scenarios = []
    for effect in (0.20, 0.293):
        for design_effect in (1.0, 1.5):
            scenarios.append({
                "planning_partial_r": effect, "target_power": 0.80,
                "multiplicity_scenario": "bonferroni_168",
                "design_effect": design_effect, "inflated_usable_rows": 500,
            })
    return {
        "input": {
            "current_parent_lenses": 800,
            "current_joint_support_rows": {"i": 26, "i2": 68},
            "observed_support_yield": {"i": 0.0325, "i2": 0.085},
        },
        "scenarios": scenarios,
    }


def test_cross_release_pool_applies_redshift_spec_current_and_separation_cuts():
    current = [candidate(f"current-{index}", 0.0, 1.0) for index in range(800)]
    s20a = [
        candidate("a", 180.0, 10.0),
        candidate("b", 180.05, 9.0),  # within 6 arcmin of a
        candidate("bad-z", 181.0, 50.0, z=0.6),
        candidate("bad-spec", 182.0, 50.0, z_bcg=0.5),
        candidate("current-match", 0.0001, 50.0),  # within one arcsec of current
    ]
    s23b = [
        candidate("a-new", 180.0, 10.0),  # exact cross-release duplicate; tie sorts to S20A
        candidate("c", 180.2, 8.0),
    ]
    result, additions = evaluate(current, s20a, s23b, power_plan())
    assert len(additions) == 2
    assert [row["release"] for row in additions] == ["S20A", "S23B"]
    assert result["catalog_counts"]["s20a_s23b_eligible_cross_pairs_within_1arcsec"] == 1
    assert result["catalog_counts"]["eligible_rows_within_1arcsec_of_current_targets"] == 1
    assert result["catalog_counts"]["greedy_additions"] == 2
    assert result["cas_query_submitted"] is False
    assert result["science_measurement"] is False


def test_cross_release_output_projects_support_as_assumption_only():
    current = [candidate(f"current-{index}", 0.0, 1.0) for index in range(800)]
    result, _ = evaluate(current, [candidate("one", 180.0, 10.0)], [], power_plan())
    assert result["support_projection"]["projected_total_support_rows_if_yield_transfers"] == {
        "i": 26.0325, "i2": 68.085,
    }
    assert result["support_projection"]["sample_size_scenarios"][0]["meets_scenario"] is False
