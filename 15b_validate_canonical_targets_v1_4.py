#!/usr/bin/env python3
"""Validate the Pass-10 semantic canonical-target vocabulary without writes."""
from __future__ import annotations

import importlib.util
from pathlib import Path


SOURCE = Path(__file__).with_name("15b_validate_canonical_targets_v1_3.py")
SPEC = importlib.util.spec_from_file_location("canonical_targets_v1_3_base", SOURCE)
base = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(base)

base.VALIDATOR_VERSION = "canonical-target-browser-validation-pass10-v1.4"
base.AUTHORITY_VERSION = "canonical-target-authority-pass10-v1"
base.EXPECTED_TARGET_BROWSER_GROUPS = 34


if __name__ == "__main__":
    base.main()
