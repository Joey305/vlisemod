from copy import deepcopy
from pathlib import Path
import shutil
import zipfile

from lxml import etree


NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
TARGET = Path(
    "/Users/jxs794/Library/CloudStorage/OneDrive-UniversityofMiami/VLiSEMOD/"
    "Supplementary/Re-Submission-Package/Resubmission/"
    "ACS_Infectious_Disease_2026_Schulz_etAL_Revised.docx"
)
NEEDLE = "V-LiSEMOD complements the resources that supply these individual data types."
REVISED_TEXT = (
    "V-LiSEMOD complements the resources that supply these individual data types. "
    "The PDB and Chemical Component Dictionary provide structures and component records "
    "(DOIs: 10.1093/nar/28.1.235; 10.1093/bioinformatics/btu789; 10.1093/nar/gkae1091); "
    "PLIP and Arpeggio describe protein-ligand contacts "
    "(DOIs: 10.1093/nar/gkv315; 10.1016/j.jmb.2016.12.004); "
    "solvent-accessibility tools calculate surface exposure "
    "(DOI: 10.12688/f1000research.7931.1); "
    "ChEMBL organizes bioactivity and medicinal-chemistry data "
    "(DOI: 10.1093/nar/gkad1004); "
    "and PROTAC-DB catalogs reported degraders and their components "
    "(DOI: 10.1093/nar/gkae768). "
    "V-LiSEMOD keeps these pieces connected to the same viral structure and ligand atom so that "
    "a chemist can move from a bound ligand to a specific modification-site hypothesis in one workflow. "
    "PLIP and LigPlot+ are complementary interaction-analysis and visualization tools "
    "(DOIs: 10.1093/nar/gkv315; 10.1021/ci200227u); "
    "the protein–ligand contact records reported and analyzed in V-LiSEMOD were generated with "
    "pdbe-arpeggio 1.4.4 (DOI: 10.1016/j.jmb.2016.12.004)."
)


def text_of(element):
    return "".join(element.xpath(".//w:t/text()", namespaces=NS))


def replace_paragraph_keep_citation(paragraph):
    runs = paragraph.xpath("./w:r", namespaces=NS)
    first_rpr = None
    if runs:
        rpr = runs[0].find("w:rPr", namespaces=NS)
        if rpr is not None:
            first_rpr = deepcopy(rpr)
    sdt_children = paragraph.xpath("./w:sdt", namespaces=NS)
    for child in list(paragraph):
        if child.tag not in {W + "pPr", W + "sdt"}:
            paragraph.remove(child)
    run = etree.Element(W + "r")
    if first_rpr is not None:
        run.append(first_rpr)
    text = etree.SubElement(run, W + "t")
    text.text = REVISED_TEXT
    if sdt_children:
        paragraph.insert(paragraph.index(sdt_children[0]), run)
    else:
        paragraph.append(run)


def update_document():
    temporary = TARGET.with_suffix(".tmp.docx")
    changed = 0
    with zipfile.ZipFile(TARGET, "r") as zin, zipfile.ZipFile(
        temporary, "w", zipfile.ZIP_DEFLATED
    ) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "word/document.xml":
                root = etree.fromstring(data)
                for paragraph in root.xpath(".//w:p", namespaces=NS):
                    if NEEDLE in text_of(paragraph):
                        replace_paragraph_keep_citation(paragraph)
                        changed += 1
                data = etree.tostring(
                    root, xml_declaration=True, encoding="UTF-8", standalone=True
                )
            zout.writestr(item, data)
    if changed != 1:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"Expected one target paragraph, changed {changed}.")
    shutil.move(temporary, TARGET)
    print(f"Updated {TARGET}")


if __name__ == "__main__":
    update_document()
