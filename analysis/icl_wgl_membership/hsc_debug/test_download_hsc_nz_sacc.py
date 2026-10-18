import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from download_hsc_nz_sacc import (
    FILENAME,
    OFFICIAL_URL,
    download,
    validate_payload,
    validate_url,
)


class HscNzSaccDownloadTests(unittest.TestCase):
    def test_only_pinned_official_url_is_accepted(self):
        validate_url(OFFICIAL_URL)
        for url in (
            "http://hsc-release.mtk.nao.ac.jp/archive/filetree/shape_catalog_y3/dalal23/hsc_y3_fourier_space_data_vector.sacc",
            OFFICIAL_URL + "?other=1",
            OFFICIAL_URL.replace("dalal23", "other"),
            OFFICIAL_URL.replace("hsc-release", "example"),
        ):
            with self.subTest(url=url), self.assertRaises(ValueError):
                validate_url(url)

    def test_rejects_html_and_empty_payloads(self):
        for payload in (b"", b"<!doctype html><title>login</title>"):
            with self.subTest(payload=payload), self.assertRaises(ValueError):
                validate_payload(payload)

    def test_download_records_hash_without_credentials_and_refuses_overwrite(self):
        payload = b"SIMPLE  = T\nTEST SACC CONTENT\n"
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / "calibration"
            with mock.patch("download_hsc_nz_sacc.fetch", return_value=(payload, OFFICIAL_URL)):
                result = download(out, "researcher", "hidden")
            artifact = out / FILENAME
            manifest = json.loads((out / "download_manifest.json").read_text())
            self.assertEqual(artifact.read_bytes(), payload)
            self.assertEqual(result["sha256"], hashlib.sha256(payload).hexdigest())
            self.assertEqual(manifest["sha256"], result["sha256"])
            self.assertFalse(manifest["credential_values_saved"])
            with self.assertRaises(FileExistsError):
                download(out, "researcher", "hidden")


if __name__ == "__main__":
    unittest.main()
