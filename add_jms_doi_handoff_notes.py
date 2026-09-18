"""Create a red DOI handoff copy without changing Mendeley citation fields."""

from pathlib import Path
from tempfile import NamedTemporaryFile
from zipfile import ZIP_DEFLATED, ZipFile

from lxml import etree


SOURCE = Path(
    "/Users/jxs794/Library/CloudStorage/OneDrive-UniversityofMiami/VLiSEMOD/"
    "Supplementary/Re-Submission-Package/"
    "ACS_Infectious_Disease_2026_Schulz_etAL_Revised_JMS.docx"
)
OUTPUT = SOURCE.with_name("ACS_Infectious_Disease_2026_Schulz_etAL_Revised_JMS_DOI_HANDOFF.docx")
W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
NS = {"w": W_NS}


def tag(name):
    return f"{{{W_NS}}}{name}"


def text_of(element):
    return "".join(element.xpath(".//w:t/text()", namespaces=NS))


def red_run(note):
    run = etree.Element(tag("r"))
    r_pr = etree.SubElement(run, tag("rPr"))
    color = etree.SubElement(r_pr, tag("color"))
    color.set(tag("val"), "FF0000")
    bold = etree.SubElement(r_pr, tag("b"))
    size = etree.SubElement(r_pr, tag("sz"))
    size.set(tag("val"), "18")
    text = etree.SubElement(run, tag("t"))
    text.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    text.text = " " + note
    return run


def paragraph_starting(root, phrase):
    for paragraph in root.xpath(".//w:body/w:p", namespaces=NS):
        if text_of(paragraph).startswith(phrase):
            return paragraph
    raise ValueError(f"Paragraph not found: {phrase}")


def sdt_with_display(paragraph, display):
    for sdt in paragraph.findall(tag("sdt")):
        if text_of(sdt) == display:
            return sdt
    raise ValueError(f"Citation field {display!r} not found")


def add_note_after(paragraph, child, note):
    paragraph.insert(paragraph.index(child) + 1, red_run(note))


def update_xml(data):
    root = etree.fromstring(data)

    # Original manuscript evidence: no citation follows the Figure 5B contact-count
    # paragraph. The present four-source field supports the Agreement paragraph.
    agreement = paragraph_starting(root, "4,28,33,34Agreement across several structures")
    agreement.append(
        red_run(
            "[MOVE 4,28,33,34 HERE. DOI: 10.1038/nchembio.2329; "
            "10.1038/nrmicro2747; 10.1038/s41467-023-39904-5; "
            "10.1016/j.coph.2016.08.014.]"
        )
    )

    # The generic RDKit concept DOI remains in two fields. Its version-specific
    # replacement is already correctly present in the Methods RDKit statement.
    resources = paragraph_starting(root, "Answering that question manually")
    add_note_after(
        resources,
        sdt_with_display(resources, "9\u201317"),
        "[REPLACE ONLY THE RDKit ITEM WITH DOI: 10.5281/zenodo.19922430; "
        "REMOVE OBSOLETE DOI: 10.5281/zenodo.591637.]",
    )
    cheminformatics = paragraph_starting(root, "Cheminformatics-derived ligand depictions")
    add_note_after(
        cheminformatics,
        sdt_with_display(cheminformatics, "17"),
        "[REPLACE WITH RDKit DOI: 10.5281/zenodo.19922430; "
        "REMOVE OBSOLETE DOI: 10.5281/zenodo.591637.]",
    )

    # These two sources already occur after their respective source sentences;
    # only the final duplicated group needs to be removed.
    retrospective = paragraph_starting(root, "Two published HIV-1 protease series")
    retrospective.append(
        red_run(
            "[REMOVE DUPLICATE FINAL 45,46. Hidaka DOI 10.1021/acs.bioconjchem.9b00195 "
            "and Kawasaki DOI 10.1111/j.1747-0285.2009.00921.x are already cited "
            "after their respective sentences.]"
        )
    )
    return etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)


def main():
    with ZipFile(SOURCE) as source:
        updated = update_xml(source.read("word/document.xml"))
        with NamedTemporaryFile(dir=OUTPUT.parent, suffix=".docx", delete=False) as temp:
            temp_path = Path(temp.name)
        try:
            with ZipFile(temp_path, "w", ZIP_DEFLATED) as destination:
                for info in source.infolist():
                    payload = updated if info.filename == "word/document.xml" else source.read(info.filename)
                    destination.writestr(info, payload)
            temp_path.replace(OUTPUT)
        finally:
            if temp_path.exists():
                temp_path.unlink()


if __name__ == "__main__":
    main()
