import json

import pytest

from audit_icl_s23b_expansion_feasibility import run


def test_rejects_unpinned_or_modified_catalogue_inputs(tmp_path):
    catalog = tmp_path / "catalog"
    catalog.mkdir()
    (catalog / "summary.json").write_text(json.dumps({"schema": "wrong"}))
    (catalog / "greedy_candidate_additions.tsv").write_text("id\n1\n")
    power = tmp_path / "power.json"
    power.write_text("{}")
    with pytest.raises(ValueError, match="summary is not the reviewed pinned result"):
        run(catalog, power)
