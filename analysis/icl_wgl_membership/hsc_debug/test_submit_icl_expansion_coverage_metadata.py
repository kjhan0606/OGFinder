import pytest

from submit_icl_expansion_coverage_metadata import check_submission_authorization


PLAN_SHA = "e615ee7c2414109af6d3c156104d71766355d3b4e2a5dba3636a0142995557d7"
S23B_SCHEMA = "astrafek-s23b-pdr3-coverage-query-plan-v1"
S20A_SCHEMA = "astrafek-icl-pdr3-coverage-query-plan-v1"


def test_s23b_submission_requires_exact_plan_approval():
    with pytest.raises(ValueError, match="exact plan SHA-256"):
        check_submission_authorization(S23B_SCHEMA, PLAN_SHA, None)
    with pytest.raises(ValueError, match="exact plan SHA-256"):
        check_submission_authorization(S23B_SCHEMA, PLAN_SHA, "0" * 64)
    check_submission_authorization(S23B_SCHEMA, PLAN_SHA, PLAN_SHA)


def test_existing_s20a_approval_flow_remains_unchanged():
    check_submission_authorization(S20A_SCHEMA, "a" * 64, None)
    with pytest.raises(ValueError, match="does not match"):
        check_submission_authorization(S20A_SCHEMA, "a" * 64, "b" * 64)
