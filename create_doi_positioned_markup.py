"""Create fresh, non-stale clean and redline copies with inline DOI notes.

The markup copy follows the supplied Pass 10 convention: existing deletions
remain black and struck through; each new DOI addition is ordinary red text.
"""

from pathlib import Path
from shutil import copy2
from tempfile import NamedTemporaryFile
from zipfile import ZIP_DEFLATED, ZipFile

from lxml import etree

from move_sentence_doi_notes import update_docx


CLEAN_SOURCE = Path(
    "/Users/jxs794/Library/CloudStorage/OneDrive-UniversityofMiami/VLiSEMOD/"
    "Supplementary/Re-Submission-Package/"
    "ACS_Infectious_Disease_2026_Schulz_etAL_Revised_MENDELEY_REPAIRED.docx"
)
CLEAN_OUTPUT = CLEAN_SOURCE.with_name(
    "ACS_Infectious_Disease_2026_Schulz_etAL_Revised_MENDELEY_REPAIRED_DOI_POSITIONED.docx"
)
MARKUP_SOURCE = Path(
    "/Users/jxs794/Library/CloudStorage/OneDrive-UniversityofMiami/VLiSEMOD/"
    "Supplementary/ACS_Infectious_Disease_2026_Schulz_etAL_Reviewer_Pass10_METHOD_MARKUP.docx"
)
MARKUP_OUTPUT = MARKUP_SOURCE.with_name(
    "ACS_Infectious_Disease_2026_Schulz_etAL_Reviewer_Pass10_METHOD_MARKUP_DOI_POSITIONED.docx"
)
W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
NS = {"w": W_NS}

POCKET_DOIS = (
    " [DOI: 10.1093/nar/28.1.235; 10.1038/s41573-021-00371-6; "
    "10.1016/j.cbpa.2019.02.022; 10.1038/nchembio.2329.]"
)
ACTIVITY_DOIS = (
    " [DOI: 10.1021/acs.jcim.0c00589; 10.1021/acs.jmedchem.8b01413; "
    "10.1016/j.bmcl.2024.129676; 10.1038/s41467-018-08027-7.]"
)


def tag(name):
    return f"{{{W_NS}}}{name}"


def red_run(text):
    run = etree.Element(tag("r"))
    r_pr = etree.SubElement(run, tag("rPr"))
    # Mirrors the supplied Pass 10 markup's addition styling.
    strike = etree.SubElement(r_pr, tag("strike"))
    strike.set(tag("val"), "0")
    color = etree.SubElement(r_pr, tag("color"))
    color.set(tag("val"), "FF0000")
    text_element = etree.SubElement(run, tag("t"))
    text_element.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    text_element.text = text
    return run


def add_inline_dois(document_xml):
    root = etree.fromstring(document_xml)
    target = next(
        (
            paragraph
            for paragraph in root.xpath(".//w:body/w:p", namespaces=NS)
            if "".join(paragraph.xpath(".//w:t/text()", namespaces=NS)).startswith(
                "For a medicinal chemist"
            )
        ),
        None,
    )
    if target is None:
        raise ValueError("The medicinal-chemist paragraph was not found.")

    citation_sdts = [
        child
        for child in target.findall(tag("sdt"))
        if "".join(child.xpath(".//w:t/text()", namespaces=NS)) in {"1\u20135", "6\u20139"}
    ]
    if len(citation_sdts) != 2:
        raise ValueError("Expected the two existing citation fields after the target sentences.")
    for citation, doi_note in reversed(list(zip(citation_sdts, (POCKET_DOIS, ACTIVITY_DOIS)))):
        target.insert(target.index(citation) + 1, red_run(doi_note))
    return etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)


def write_markup_output():
    with ZipFile(MARKUP_SOURCE) as source:
        document_xml = add_inline_dois(source.read("word/document.xml"))
        with NamedTemporaryFile(dir=MARKUP_OUTPUT.parent, suffix=".docx", delete=False) as tmp:
            temp_path = Path(tmp.name)
        try:
            with ZipFile(temp_path, "w", ZIP_DEFLATED) as destination:
                for info in source.infolist():
                    content = document_xml if info.filename == "word/document.xml" else source.read(info.filename)
                    destination.writestr(info, content)
            temp_path.replace(MARKUP_OUTPUT)
        finally:
            if temp_path.exists():
                temp_path.unlink()


def main():
    # A fresh filename avoids Word displaying an already-open, stale document buffer.
    copy2(CLEAN_SOURCE, CLEAN_OUTPUT)
    # The open Word document has retained the earlier paragraph-end handoff text.
    # Move those DOI groups into their exact sentences in the fresh copy.
    update_docx(CLEAN_OUTPUT)
    write_markup_output()


if __name__ == "__main__":
    main()
