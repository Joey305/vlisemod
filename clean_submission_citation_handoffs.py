from copy import deepcopy
from pathlib import Path
import shutil
import zipfile

from lxml import etree


NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

PACKAGE = Path(
    "/Users/jxs794/Library/CloudStorage/OneDrive-UniversityofMiami/VLiSEMOD/"
    "Supplementary/Re-Submission-Package"
)
RESPONSE_SOURCE = PACKAGE / "Response_to_Reviewer_Comments_MENDELEY_REPAIRED.docx"
RESPONSE_OUTPUT = PACKAGE / "Response_to_Reviewer_Comments_SUBMISSION_CLEAN.docx"
S5_SOURCE = PACKAGE / "Supplementary_File_5_Supplemental_MENDELEY_REPAIRED.docx"
S5_OUTPUT = PACKAGE / "Supplementary_File_5_Supplemental_SUBMISSION_CLEAN.docx"


def paragraph_text(paragraph):
    return "".join(paragraph.xpath(".//w:t/text()", namespaces=NS))


def replace_paragraph_text(paragraph, text):
    """Keep paragraph properties and the first run's base formatting."""
    runs = paragraph.xpath("./w:r", namespaces=NS)
    first_rpr = None
    if runs:
        rpr = runs[0].find("w:rPr", namespaces=NS)
        if rpr is not None:
            first_rpr = deepcopy(rpr)
    for child in list(paragraph):
        if child.tag != W + "pPr":
            paragraph.remove(child)
    run = etree.SubElement(paragraph, W + "r")
    if first_rpr is not None:
        run.append(first_rpr)
    t = etree.SubElement(run, W + "t")
    if text.startswith(" ") or text.endswith(" "):
        t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    t.text = text


def rewrite_docx(source, output, transform):
    with zipfile.ZipFile(source, "r") as zin, zipfile.ZipFile(
        output, "w", zipfile.ZIP_DEFLATED
    ) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "word/document.xml":
                root = etree.fromstring(data)
                transform(root)
                data = etree.tostring(
                    root, xml_declaration=True, encoding="UTF-8", standalone=True
                )
            zout.writestr(item, data)


def clean_response(root):
    for p in root.xpath(".//w:p", namespaces=NS):
        text = paragraph_text(p)
        if "[REMOVE CURRENT CITATION HERE]" in text:
            replace_paragraph_text(
                p,
                "Response. We removed the residual Southan citation from the affected "
                "opening-design sentence because it did not directly support that statement. "
                "The finalized manuscript retains directly relevant structural and PROTAC "
                "literature only at the appropriate sentences.",
            )
        elif "[FINAL RDKit DOI:" in text:
            replace_paragraph_text(
                p,
                "Response. We replaced the generic RDKit concept DOI with the "
                "version-specific DOI for RDKit 2026_03_2 (10.5281/zenodo.19922430) "
                "in the Methods sentence. The reference list has been updated accordingly.",
            )


def clean_s5(root):
    reference_label_updates = {
        "RCSB PDB / wwPDB Chemical Component Dictionary (refs. 1, 10, 40)":
            "RCSB PDB / wwPDB Chemical Component Dictionary (refs. 1, 9, 19)",
        "PLIP / Arpeggio (refs. 13, 14)": "PLIP / Arpeggio (refs. 12, 13)",
        "Solvent-accessibility tools / FreeSASA (ref. 15)":
            "Solvent-accessibility tools / FreeSASA (ref. 14)",
        "ChEMBL (ref. 11)": "ChEMBL (ref. 10)",
        "PROTAC-DB 3.0 (ref. 34)": "PROTAC-DB 3.0 (ref. 18)",
    }
    for p in root.xpath(".//w:p", namespaces=NS):
        text = paragraph_text(p)
        for current, corrected in reference_label_updates.items():
            if text.startswith(current):
                replace_paragraph_text(p, corrected)
                break
        else:
            current = None
        if current is not None:
            continue
        if "[DOI TO ADD MANUALLY: 10.1038/s41589-020-00689-z]" in text:
            clean = text.replace(
                " [DOI TO ADD MANUALLY: 10.1038/s41589-020-00689-z]", " (DOI: 10.1038/s41589-020-00689-z)"
            )
            clean = clean.replace(
                "  [FINAL CITE DOI: 10.1038/s41589-020-00689-z; REMOVE TEMPORARY DOI MARKER.]", ""
            )
            replace_paragraph_text(p, clean)
            continue
        for run in list(p.xpath(".//w:r", namespaces=NS)):
            run_text = "".join(run.xpath(".//w:t/text()", namespaces=NS))
            if "[FINAL CITE DOI:" in run_text:
                run.getparent().remove(run)


if __name__ == "__main__":
    rewrite_docx(RESPONSE_SOURCE, RESPONSE_OUTPUT, clean_response)
    rewrite_docx(S5_SOURCE, S5_OUTPUT, clean_s5)
    print(RESPONSE_OUTPUT)
    print(S5_OUTPUT)
