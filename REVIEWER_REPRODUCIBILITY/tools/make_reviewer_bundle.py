#!/usr/bin/env python3
"""Build and audit the distributable reviewer ZIP without derived run artifacts."""
from __future__ import annotations

import argparse
import re
import shutil
import tempfile
import zipfile
from pathlib import Path


PACKAGE_NAME = "VLiSEMOD_Reviewer_Reproducibility_v1.0"
EXCLUDED_PARTS = {"__pycache__", ".git", ".DS_Store"}
EXCLUDED_SUFFIXES = {".pyc", ".db", ".sqlite", ".zip", ".log"}
PRIVATE_PATTERNS = (
    re.compile(b"/" + b"Users/"),
    re.compile(b"jxs" + b"794"),
    re.compile(b"Her" + b"oku", re.I),
    re.compile(b"Bear" + b"er\\s+", re.I),
)
AUTHOR_ONLY_TOOLS = {"extract_frozen_inputs.py", "sanitize_reference_artifacts.py"}


def included(root: Path, path: Path) -> bool:
    relative = path.relative_to(root)
    if relative.parts[:1] == ("tools",) and path.name in AUTHOR_ONLY_TOOLS:
        return False
    if any(part in EXCLUDED_PARTS for part in relative.parts):
        return False
    if path.suffix.lower() in EXCLUDED_SUFFIXES:
        return False
    if ".download_cache" in relative.parts or path.suffix.lower() == ".part":
        return False
    # Only the deliberately small fixture coordinates are distributable here.
    if path.suffix.lower() == ".cif" and relative.parts[:1] != ("fixture",):
        return False
    if relative.parts and relative.parts[0] == "outputs":
        return False
    return path.is_file()


def audit(paths: list[Path]) -> list[str]:
    findings = []
    for path in paths:
        try:
            data = path.read_bytes()
        except OSError as exc:
            findings.append(f"unreadable: {path}: {exc}")
            continue
        for pattern in PRIVATE_PATTERNS:
            if pattern.search(data):
                findings.append(f"private/deployment marker {pattern.pattern!r}: {path}")
                break
    return findings


def build(root: Path, destination: Path) -> None:
    files = [path for path in root.rglob("*") if included(root, path)]
    findings = audit(files)
    if findings:
        raise RuntimeError("Package privacy audit failed:\n" + "\n".join(findings))
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(files):
            archive.write(path, (Path(PACKAGE_NAME) / path.relative_to(root)).as_posix())
    with tempfile.TemporaryDirectory(prefix="vlisemod-reviewer-audit-") as tmp:
        with zipfile.ZipFile(destination) as archive:
            archive.extractall(tmp)
            bad = zipfile.ZipFile(destination).testzip()
            members = archive.namelist()
        extracted = Path(tmp) / PACKAGE_NAME
        if bad or not (extracted / "reproduce.py").is_file() or (extracted / "outputs").exists():
            raise RuntimeError(f"ZIP extraction audit failed: bad_member={bad}")
    prefix = PACKAGE_NAME + "/"
    cif_members = [member for member in members if member.lower().endswith(".cif")]
    general_cifs = [member for member in cif_members if not member.startswith(prefix + "fixture/")]
    if general_cifs:
        raise RuntimeError("ZIP includes non-fixture CIF files: " + ", ".join(general_cifs))
    print(f"bundle: {destination}")
    print(f"files: {len(files)}")
    print(f"general release CIF files: {len(general_cifs)}")
    print("fixture CIF files: " + (", ".join(Path(member).relative_to(PACKAGE_NAME).as_posix() for member in cif_members) or "none"))
    print("privacy audit: PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", default=Path(__file__).resolve().parents[2] / f"{PACKAGE_NAME}.zip")
    args = parser.parse_args()
    root, output = Path(args.root).resolve(), Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    build(root, output)


if __name__ == "__main__":
    main()
