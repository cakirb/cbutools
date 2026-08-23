import builtins

import pytest

import cbutools as cbu


def test_report_has_clear_error_without_papermill(monkeypatch):
    original_import = builtins.__import__

    def import_without_papermill(name, *args, **kwargs):
        if name == "papermill":
            raise ImportError("papermill is unavailable")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", import_without_papermill)

    with pytest.raises(ImportError, match=r"cbutools\[report\]"):
        cbu.get_barcodes_report()
