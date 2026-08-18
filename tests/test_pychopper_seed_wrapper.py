from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

SCRIPT = Path(__file__).parents[1] / "scripts" / "run_pychopper_seeded.py"


def _load_wrapper():
    spec = importlib.util.spec_from_file_location("run_pychopper_seeded", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_seed_wrapper_forwards_arguments_and_seeds_numpy(monkeypatch) -> None:
    observed: dict[str, object] = {}
    fake_numpy = types.ModuleType("numpy")

    class FakeRandom:
        @staticmethod
        def seed(value: int) -> None:
            observed["seed"] = value

    fake_numpy.random = FakeRandom()  # type: ignore[attr-defined]
    fake_package = types.ModuleType("pychopper")
    fake_package.__path__ = []  # type: ignore[attr-defined]
    fake_scripts = types.ModuleType("pychopper.scripts")
    fake_scripts.__path__ = []  # type: ignore[attr-defined]
    fake_command = types.ModuleType("pychopper.scripts.pychopper")

    def fake_main() -> int:
        observed["argv"] = list(sys.argv)
        return 0

    fake_command.main = fake_main  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "numpy", fake_numpy)
    monkeypatch.setitem(sys.modules, "pychopper", fake_package)
    monkeypatch.setitem(sys.modules, "pychopper.scripts", fake_scripts)
    monkeypatch.setitem(sys.modules, "pychopper.scripts.pychopper", fake_command)
    monkeypatch.setattr(sys, "argv", ["wrapper.py"])

    wrapper = _load_wrapper()
    assert wrapper.main(["20261005", "-k", "PCS109", "input.fq", "-"]) == 0
    assert observed == {
        "seed": 20261005,
        "argv": ["wrapper.py", "-k", "PCS109", "input.fq", "-"],
    }


def test_seed_wrapper_rejects_invalid_seed() -> None:
    wrapper = _load_wrapper()
    try:
        wrapper.main(["not-an-integer", "input.fq", "-"])
    except SystemExit as error:
        assert str(error) == "RANDOM_SEED must be an integer"
    else:
        raise AssertionError("invalid seed was accepted")
