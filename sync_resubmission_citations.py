from copy import deepcopy
from pathlib import Path
import base64
import json
import shutil
import zipfile

from lxml import etree


NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
PREFIX = "MENDELEY_CITATION_v3_"
RESUB = Path(
    "/Users/jxs794/Library/CloudStorage/OneDrive-UniversityofMiami/VLiSEMOD/"
    "Supplementary/Re-Submission-Package/Resubmission"
)
FINAL = RESUB / "ACS_Infectious_Disease_2026_Schulz_etAL_Revised.docx"
MARKUP = RESUB / "ACS_Infectious_Disease_2026_Schulz_etAL_Revised_MARKUP.docx"
S5 = RESUB / "Supplementary_File_5_Supplemental_SUBMISSION_CLEAN.docx"


def text_of(element):
    return "".join(element.xpath(".//w:t/text()", namespaces=NS))


def paragraph_of(element):
    current = element
    while current is not None and current.tag != W + "p":
        current = current.getparent()
    return current


def citation_tag(sdt):
    tag = sdt.find("w:sdtPr/w:tag", namespaces=NS)
    return tag.get(W + "val") if tag is not None else ""


def decode_tag(value):
    encoded = value[len(PREFIX) :]
    encoded += "=" * (-len(encoded) % 4)
    return json.loads(base64.urlsafe_b64decode(encoded).decode("utf-8"))


def citation_key(sdt):
    data = decode_tag(citation_tag(sdt))
    items = []
    for item in data.get("citationItems", []):
        source = item.get("itemData", {})
        items.append((source.get("DOI") or source.get("title") or "").lower())
    return tuple(items)


def field_text(sdt):
    return text_of(sdt)


def set_field_from_template(destination, source):
    old_tag = destination.find("w:sdtPr/w:tag", namespaces=NS)
    new_tag = source.find("w:sdtPr/w:tag", namespaces=NS)
    if old_tag is None or new_tag is None:
        raise RuntimeError("Mendeley tag not found")
    old_tag.set(W + "val", new_tag.get(W + "val"))
    text_nodes = destination.xpath(".//w:sdtContent//w:t", namespaces=NS)
    if not text_nodes:
        raise RuntimeError("Citation display text not found")
    text_nodes[0].text = field_text(source)
    for node in text_nodes[1:]:
        node.text = ""


def rewrite_docx(path, transform):
    temporary = path.with_suffix(".sync.tmp.docx")
    with zipfile.ZipFile(path, "r") as zin, zipfile.ZipFile(
        temporary, "w", zipfile.ZIP_DEFLATED
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
    shutil.move(temporary, path)


def final_field_map():
    return {citation_key(sdt): sdt for sdt in final_fields_all()}


def final_fields_all():
    with zipfile.ZipFile(FINAL) as z:
        root = etree.fromstring(z.read("word/document.xml"))
    result = []
    for sdt in root.xpath(".//w:sdt", namespaces=NS):
        if citation_tag(sdt).startswith(PREFIX):
            result.append(sdt)
    return result


def final_bibliography_content():
    with zipfile.ZipFile(FINAL) as z:
        root = etree.fromstring(z.read("word/document.xml"))
    for sdt in root.xpath(".//w:sdt", namespaces=NS):
        if citation_tag(sdt) == "MENDELEY_BIBLIOGRAPHY":
            content = sdt.find("w:sdtContent", namespaces=NS)
            if content is not None:
                return [deepcopy(child) for child in content]
    raise RuntimeError("Final manuscript Mendeley bibliography not found")


def find_final_field(fields, sentence_start, display=None):
    for field in fields:
        paragraph = paragraph_of(field)
        if paragraph is not None and sentence_start in text_of(paragraph):
            if display is None or field_text(field) == display:
                return field
    raise RuntimeError(f"Missing final citation field for {sentence_start!r}")


def sync_markup(root):
    final_fields = final_field_map()
    final_field_list = final_fields_all()
    exact_updates = 0
    special_updates = 0
    for sdt in root.xpath(".//w:sdt", namespaces=NS):
        tag = citation_tag(sdt)
        if not tag.startswith(PREFIX):
            continue
        key = citation_key(sdt)
        if key in final_fields:
            if field_text(sdt) != field_text(final_fields[key]):
                set_field_from_template(sdt, final_fields[key])
                exact_updates += 1
            continue
        paragraph = paragraph_of(sdt)
        context = text_of(paragraph) if paragraph is not None else ""
        if "Answering that question manually usually requires" in context and "10.5281/zenodo.591637" in key:
            source = find_final_field(final_field_list, "Answering that question manually usually requires", "9–17")
            set_field_from_template(sdt, source)
            special_updates += 1
        elif "Cheminformatics-derived ligand depictions" in context and "10.5281/zenodo.591637" in key:
            source = find_final_field(final_field_list, "Cheminformatics-derived ligand depictions")
            set_field_from_template(sdt, source)
            special_updates += 1

    resource_text = "V-LiSEMOD complements the resources that supply these individual data types."
    for paragraph in root.xpath(".//w:p", namespaces=NS):
        if resource_text in text_of(paragraph):
            red_runs = [
                run
                for run in paragraph.xpath("./w:r", namespaces=NS)
                if run.find("w:rPr/w:color", namespaces=NS) is not None
                and run.find("w:rPr/w:color", namespaces=NS).get(W + "val") == "FF0000"
                and text_of(run) == "1,10,11,13–15,34,40"
            ]
            if len(red_runs) == 1:
                node = red_runs[0].find("w:t", namespaces=NS)
                node.text = "1,9,18; 12,13; 14; 10; 19–21"
                special_updates += 1
            elif len(red_runs) != 0:
                raise RuntimeError("Expected one active resource-citation run in markup")

    destination_bibliography = None
    for sdt in root.xpath(".//w:sdt", namespaces=NS):
        if citation_tag(sdt) == "MENDELEY_BIBLIOGRAPHY":
            destination_bibliography = sdt.find("w:sdtContent", namespaces=NS)
            break
    if destination_bibliography is None:
        raise RuntimeError("Markup Mendeley bibliography not found")
    for child in list(destination_bibliography):
        destination_bibliography.remove(child)
    for child in final_bibliography_content():
        destination_bibliography.append(child)


def set_plain_paragraph_text(paragraph, text):
    runs = paragraph.xpath("./w:r", namespaces=NS)
    rpr = deepcopy(runs[0].find("w:rPr", namespaces=NS)) if runs and runs[0].find("w:rPr", namespaces=NS) is not None else None
    for child in list(paragraph):
        if child.tag != W + "pPr":
            paragraph.remove(child)
    run = etree.SubElement(paragraph, W + "r")
    if rpr is not None:
        run.append(rpr)
    node = etree.SubElement(run, W + "t")
    node.text = text


def sync_s5(root):
    updates = {
        "RCSB PDB / wwPDB Chemical Component Dictionary (refs. 1, 9, 19)":
            "RCSB PDB / wwPDB Chemical Component Dictionary (refs. 1, 9, 18)",
        "PROTAC-DB 3.0 (ref. 18)": "PROTAC-DB 3.0 (refs. 19–21)",
    }
    count = 0
    for paragraph in root.xpath(".//w:p", namespaces=NS):
        current = text_of(paragraph)
        if current in updates:
            set_plain_paragraph_text(paragraph, updates[current])
            count += 1
    if count not in (0, len(updates)):
        raise RuntimeError(f"Expected zero or {len(updates)} S5 label updates, found {count}")


if __name__ == "__main__":
    rewrite_docx(MARKUP, sync_markup)
    rewrite_docx(S5, sync_s5)
    print(MARKUP)
    print(S5)
