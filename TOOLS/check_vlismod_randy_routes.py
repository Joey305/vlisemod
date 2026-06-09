#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from typing import Any

import requests


DEFAULT_VIRUS = os.environ.get("VLISMOD_SAMPLE_VIRUS", "Human immunodeficiency virus 1")
DEFAULT_PDB = os.environ.get("VLISMOD_SAMPLE_PDB", "3EKY")
DEFAULT_LIGAND = os.environ.get("VLISMOD_SAMPLE_LIGAND", "DR7")
DEFAULT_CHAIN = os.environ.get("VLISMOD_SAMPLE_CHAIN", "A")
DEFAULT_LIGAND_ID = os.environ.get("VLISMOD_SAMPLE_LIGAND_ID", "100")
DEFAULT_PROTEIN_TYPE = os.environ.get("VLISMOD_SAMPLE_PROTEIN_TYPE", "Capsid Protein")
TIMEOUT_SECONDS = int(os.environ.get("VLISMOD_ROUTE_CHECK_TIMEOUT", "60"))


@dataclass
class CheckResult:
    name: str
    ok: bool
    status_code: int | None
    detail: str


def require_env(name: str) -> str:
    value = os.environ.get(name, "").strip()
    if not value:
        raise SystemExit(f"Missing required environment variable: {name}")
    return value


def call_json(
    session: requests.Session,
    method: str,
    url: str,
    *,
    token: str | None = None,
    timeout: int = TIMEOUT_SECONDS,
    **kwargs: Any,
) -> requests.Response:
    headers = dict(kwargs.pop("headers", {}) or {})
    if token:
        headers["Authorization"] = f"Bearer {token}"
    return session.request(method, url, headers=headers, timeout=timeout, **kwargs)


def check_response(name: str, response: requests.Response, expected_status: int) -> CheckResult:
    ok = response.status_code == expected_status
    detail = f"expected {expected_status}, got {response.status_code}"
    content_type = response.headers.get("content-type", "")
    if content_type.startswith("application/json"):
        try:
            payload = response.json()
            if isinstance(payload, dict):
                detail = ", ".join(sorted(payload.keys())[:8]) or "json"
            elif isinstance(payload, list):
                detail = f"list[{len(payload)}]"
        except Exception:
            detail = "json-parse-failed"
    return CheckResult(name=name, ok=ok, status_code=response.status_code, detail=detail)


def main() -> int:
    backup_url = require_env("VLISMOD_BACKUP_URL").rstrip("/")
    token = require_env("RANDY_API_TOKEN")
    app_base = os.environ.get("VLISMOD_APP_BASE_URL", "").strip().rstrip("/")

    session = requests.Session()
    checks: list[CheckResult] = []

    direct_gets = [
        ("db-health", f"{backup_url}/db-health", 200),
        ("ligands-with-synonyms", f"{backup_url}/ligands/with-synonyms", 200),
        ("ligand-info", f"{backup_url}/ligand-info?ligand_code={DEFAULT_LIGAND}", 200),
        ("ligand-options", f"{backup_url}/ligand-options", 200),
        ("ligand-smiles", f"{backup_url}/ligand-smiles?ligand_id={DEFAULT_LIGAND}", 200),
        ("pdb-mapping", f"{backup_url}/pdb-mapping?ligand_code={DEFAULT_LIGAND}", 200),
        ("pdb-residues-by-ligand", f"{backup_url}/pdb-residues/by-ligand?ligand_code={DEFAULT_LIGAND}", 200),
        ("virus-names", f"{backup_url}/virus-proteins/virus-names", 200),
        ("protein-types", f"{backup_url}/virus-proteins/protein-types", 200),
        ("protacability-source", f"{backup_url}/protacability/source", 200),
    ]
    for name, url, expected in direct_gets:
        checks.append(check_response(name, call_json(session, "GET", url, token=token), expected))

    direct_posts = [
        (
            "ligand-images-data",
            f"{backup_url}/ligand-images-data",
            {"virus_name": DEFAULT_VIRUS, "pdb_code": DEFAULT_PDB, "ligand_name": DEFAULT_LIGAND, "chain": DEFAULT_CHAIN},
            200,
        ),
        (
            "pymol-session-data",
            f"{backup_url}/pymol-session-data",
            {
                "pdb_code": DEFAULT_PDB,
                "ligand_name": DEFAULT_LIGAND,
                "chain": DEFAULT_CHAIN,
                "options": {
                    "functional_groups": True,
                    "binding_pocket": True,
                    "distal_atoms": True,
                    "solvent_exposed_atoms": True,
                    "hydrated_atoms": True,
                    "rupley_sasa": True,
                },
            },
            200,
        ),
        (
            "compare-interactions",
            f"{backup_url}/ligand-interactions/compare",
            {"ligand": DEFAULT_LIGAND, "pdb_ids": [f"{DEFAULT_PDB}-{DEFAULT_LIGAND_ID}-{DEFAULT_CHAIN}"]},
            200,
        ),
        (
            "virus-proteins-pdbs",
            f"{backup_url}/virus-proteins/pdbs",
            {"virus_name": DEFAULT_VIRUS, "protein_types": [DEFAULT_PROTEIN_TYPE]},
            200,
        ),
        (
            "export-data",
            f"{backup_url}/export-data",
            {"pdb_codes": [DEFAULT_PDB], "data_sets": ["Interatomic Interactions"]},
            200,
        ),
    ]
    for name, url, payload, expected in direct_posts:
        checks.append(check_response(name, call_json(session, "POST", url, token=token, json=payload), expected))

    checks.append(
        check_response(
            "unauthorized-pymol-session-data",
            call_json(
                session,
                "POST",
                f"{backup_url}/pymol-session-data",
                json={"pdb_code": DEFAULT_PDB, "ligand_name": DEFAULT_LIGAND, "chain": DEFAULT_CHAIN, "options": {"functional_groups": True}},
            ),
            401,
        )
    )

    if app_base:
        app_gets = [
            ("app-get-viruses", f"{app_base}/get_viruses", 200),
            ("app-get-ligands-with-synonyms", f"{app_base}/get_ligands_with_synonyms", 200),
            ("app-get-pdb-codes", f"{app_base}/get_pdb_codes/{DEFAULT_VIRUS}", 200),
            ("app-get-ligands", f"{app_base}/get_ligands/{DEFAULT_PDB}", 200),
            ("app-check-functional-groups", f"{app_base}/check_functional_groups/{DEFAULT_PDB}", 200),
            ("app-protacability-filters", f"{app_base}/api/protacability/filter_options", 200),
        ]
        for name, url, expected in app_gets:
            checks.append(check_response(name, call_json(session, "GET", url), expected))

        app_posts = [
            (
                "app-generate-ligand-images",
                f"{app_base}/generate_ligand_images",
                {"virus": DEFAULT_VIRUS, "pdb_code": DEFAULT_PDB, "ligand": DEFAULT_LIGAND, "chain": DEFAULT_CHAIN},
                200,
                False,
            ),
            (
                "app-generate-pymol-session",
                f"{app_base}/generate_pymol_session",
                {"pdb_code": DEFAULT_PDB, "ligand": DEFAULT_LIGAND, "chain": DEFAULT_CHAIN, "functional_groups": "on"},
                200,
                False,
            ),
            (
                "app-compare-ligand-interactions",
                f"{app_base}/compare_ligand_interactions",
                {"ligand": DEFAULT_LIGAND, "pdb_ids": [f"{DEFAULT_PDB}-{DEFAULT_LIGAND_ID}-{DEFAULT_CHAIN}"]},
                200,
                True,
            ),
        ]
        for name, url, payload, expected, send_json in app_posts:
            if send_json:
                response = call_json(session, "POST", url, json=payload)
            else:
                response = call_json(session, "POST", url, data=payload)
            checks.append(check_response(name, response, expected))

    failures = [check for check in checks if not check.ok]
    for check in checks:
        prefix = "PASS" if check.ok else "FAIL"
        status = check.status_code if check.status_code is not None else "ERR"
        print(f"[{prefix}] {check.name}: {status} ({check.detail})")

    if failures:
        print(f"\n{len(failures)} check(s) failed.", file=sys.stderr)
        return 1

    print(f"\nAll {len(checks)} checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
