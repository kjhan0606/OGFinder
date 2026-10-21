from __future__ import annotations

import pytest

import check_wl_cas_job_status as checker


def test_queries_only_existing_job_status_and_clears_password(monkeypatch):
    calls = []

    class FakeCas:
        def __init__(self, user, password):
            self.credential = {"user": user, "password": password}

        def json(self, action, **kwargs):
            calls.append((action, kwargs))
            return {"id": kwargs["id"], "status": "queued"}

    monkeypatch.setattr(checker, "Cas", FakeCas)
    monkeypatch.setattr(checker.getpass, "getpass", lambda prompt: "hidden")
    result = checker.query_status("12345")
    assert calls == [("status", {"id": "12345"})]
    assert result["status"] == "queued"
    assert result["error"] is None
    assert result["query_submitted"] is False
    assert result["rows_downloaded"] is False
    assert result["credentials_saved"] is False


def test_rejects_non_numeric_job_id_before_prompt(monkeypatch):
    monkeypatch.setattr(checker.getpass, "getpass", lambda prompt: pytest.fail("should not prompt"))
    with pytest.raises(ValueError, match="positive integer"):
        checker.query_status("not-a-job")
