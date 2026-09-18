import base64
import json
import re
import sys
import unicodedata
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

from lxml import etree


NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
PREFIX = "MENDELEY_CITATION_v3_"


def norm(value):
    value = unicodedata.normalize("NFKD", value or "")
    value = "".join(c for c in value if not unicodedata.combining(c))
    return re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()


def element_text(element):
    return "".join(element.xpath(".//w:t/text()", namespaces=NS))


def decode_tag(value):
    encoded = value[len(PREFIX) :]
    encoded += "=" * (-len(encoded) % 4)
    return json.loads(base64.urlsafe_b64decode(encoded).decode("utf-8"))


def source_item(data):
    issued = data.get("issued", {}).get("date-parts", [[]])
    year = issued[0][0] if issued and issued[0] else None
    authors = data.get("author", [])
    first_author = authors[0].get("family", "") if authors else ""
    return {
        "doi": (data.get("DOI") or "").lower(),
        "title": data.get("title", "").replace("\n", " ").strip(),
        "title_norm": norm(data.get("title", "")),
        "author": first_author,
        "year": year,
        "identity": (data.get("DOI") or norm(data.get("title", "")) or f"{first_author} {year}").lower(),
    }


def containing_paragraph(node):
    current = node
    while current is not None and current.tag != W + "p":
        current = current.getparent()
    return current


def inspect(path):
    with zipfile.ZipFile(path) as archive:
        root = etree.fromstring(archive.read("word/document.xml"))
        names = archive.namelist()
    citations = []
    bibliography = []
    for sdt in root.xpath(".//w:sdt", namespaces=NS):
        props = sdt.find("w:sdtPr", namespaces=NS)
        tag = props.find("w:tag", namespaces=NS) if props is not None else None
        value = tag.get(W + "val") if tag is not None else ""
        p = containing_paragraph(sdt)
        location = (
            {
                "index": int(p.xpath("count(preceding::w:p)", namespaces=NS)) + 1,
                "text": element_text(p),
            }
            if p is not None
            else {"index": None, "text": ""}
        )
        if value.startswith(PREFIX):
            data = decode_tag(value)
            items = [source_item(x.get("itemData", {})) for x in data.get("citationItems", [])]
            citations.append(
                {
                    "paragraph": location["index"],
                    "context": location["text"],
                    "display": element_text(sdt),
                    "items": items,
                    "tag": value,
                }
            )
        elif value == "MENDELEY_BIBLIOGRAPHY":
            content = sdt.find("w:sdtContent", namespaces=NS)
            for entry in content.xpath(".//w:p", namespaces=NS) if content is not None else []:
                text = element_text(entry)
                if text:
                    bibliography.append(text)
    citations_by_source = defaultdict(list)
    for citation in citations:
        for item in citation["items"]:
            citations_by_source[item["identity"]].append(
                {"paragraph": citation["paragraph"], "display": citation["display"], "context": citation["context"]}
            )
    return {
        "file": str(path),
        "sdt_count": len(root.xpath(".//w:sdt", namespaces=NS)),
        "citation_count": len(citations),
        "bibliography_sdt": sum(1 for s in root.xpath(".//w:sdt", namespaces=NS) if (s.find("w:sdtPr/w:tag", namespaces=NS) is not None and s.find("w:sdtPr/w:tag", namespaces=NS).get(W+"val") == "MENDELEY_BIBLIOGRAPHY")),
        "citations": citations,
        "bibliography_paragraphs": bibliography,
        "cited_sources": citations_by_source,
        "parts": names,
    }


def main():
    reports = {Path(p).name: inspect(Path(p)) for p in sys.argv[1:]}
    print(json.dumps(reports, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
