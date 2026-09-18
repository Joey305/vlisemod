"""Move the two introductory DOI handoff notes to their exact sentence ends.

This is a minimal OOXML edit: it only replaces the affected paragraph's two
runs, preserving every other document part (including Mendeley fields).
"""

from copy import deepcopy
from pathlib import Path
from tempfile import NamedTemporaryFile
from zipfile import ZIP_DEFLATED, ZipFile
import os

from lxml import etree


PACKAGE = Path(
    "/Users/jxs794/Library/CloudStorage/OneDrive-UniversityofMiami/"
    "VLiSEMOD/Supplementary/Re-Submission-Package"
)
W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
NS = {"w": W_NS}

POCKET_NOTE = (
    " [DOI: 10.1093/nar/28.1.235; 10.1038/s41573-021-00371-6; "
    "10.1016/j.cbpa.2019.02.022; 10.1038/nchembio.2329.]"
)
ACTIVITY_NOTE = (
    " [DOI: 10.1021/acs.jcim.0c00589; 10.1021/acs.jmedchem.8b01413; "
    "10.1016/j.bmcl.2024.129676; 10.1038/s41467-018-08027-7.]"
)


def tag(name):
    return f"{{{W_NS}}}{name}"


def run_with_text(template, text, red=False):
    run = etree.Element(tag("r"))
    r_pr = template.find(tag("rPr"))
    if r_pr is not None:
        run.append(deepcopy(r_pr))
    if red:
        r_pr = run.find(tag("rPr"))
        if r_pr is None:
            r_pr = etree.Element(tag("rPr"))
            run.insert(0, r_pr)
        color = etree.SubElement(r_pr, tag("color"))
        color.set(tag("val"), "C00000")
        bold = etree.SubElement(r_pr, tag("b"))
        size = etree.SubElement(r_pr, tag("sz"))
        size.set(tag("val"), "18")  # 9 pt, matching the existing red notes
    text_element = etree.SubElement(run, tag("t"))
    if text.startswith(" ") or text.endswith(" "):
        text_element.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    text_element.text = text
    return run


def reposition_notes(document_xml):
    root = etree.fromstring(document_xml)
    target = None
    for paragraph in root.xpath(".//w:body/w:p", namespaces=NS):
        text = "".join(paragraph.xpath(".//w:t/text()", namespaces=NS))
        if text.startswith("For a medicinal chemist"):
            target = paragraph
            break
    if target is None:
        raise ValueError("The medicinal-chemist paragraph was not found.")

    runs = target.findall(tag("r"))
    red_note = next(
        (
            run
            for run in runs
            if "[ADD AFTER \u2018...pocket.\u2019" in "".join(run.xpath(".//w:t/text()", namespaces=NS))
        ),
        None,
    )
    if red_note is None:
        raise ValueError("The existing paragraph-end DOI note was not found.")
    source_runs = [run for run in runs if run is not red_note]
    if len(source_runs) != 1:
        raise ValueError("Expected one original text run in the target paragraph.")

    source_run = source_runs[0]
    source_text = "".join(source_run.xpath(".//w:t/text()", namespaces=NS))
    before_pocket, remainder = source_text.split("pocket.", 1)
    before_activity, after_activity = remainder.rsplit("biological activity.", 1)

    replacements = [
        run_with_text(source_run, before_pocket + "pocket."),
        run_with_text(source_run, POCKET_NOTE, red=True),
        run_with_text(source_run, before_activity + "biological activity."),
        run_with_text(source_run, ACTIVITY_NOTE, red=True),
    ]
    if after_activity:
        replacements.append(run_with_text(source_run, after_activity))

    insert_at = target.index(source_run)
    target.remove(source_run)
    target.remove(red_note)
    for offset, replacement in enumerate(replacements):
        target.insert(insert_at + offset, replacement)
    return etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)


def update_docx(path):
    with ZipFile(path) as source:
        updated_xml = reposition_notes(source.read("word/document.xml"))
        with NamedTemporaryFile(dir=path.parent, suffix=".docx", delete=False) as temp_file:
            temp_path = Path(temp_file.name)
        try:
            with ZipFile(temp_path, "w", ZIP_DEFLATED) as destination:
                for info in source.infolist():
                    data = updated_xml if info.filename == "word/document.xml" else source.read(info.filename)
                    destination.writestr(info, data)
            os.replace(temp_path, path)
        finally:
            if temp_path.exists():
                temp_path.unlink()


def main():
    for filename in (
        "ACS_Infectious_Disease_2026_Schulz_etAL_Revised_MENDELEY_REPAIRED.docx",
        "ACS_Infectious_Disease_2026_Schulz_etAL_MARKUP_MENDELEY_REPAIRED.docx",
    ):
        update_docx(PACKAGE / filename)


if __name__ == "__main__":
    main()
