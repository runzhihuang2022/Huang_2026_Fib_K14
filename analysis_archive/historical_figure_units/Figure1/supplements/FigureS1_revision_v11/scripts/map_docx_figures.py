from __future__ import annotations

from pathlib import Path
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
WORD = ROOT / "inputs" / "old_supplement_docx" / "word"
REL_NS = "{http://schemas.openxmlformats.org/package/2006/relationships}"
WP_NS = "{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}"
A_NS = "{http://schemas.openxmlformats.org/drawingml/2006/main}"
R_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"

rels = ET.parse(WORD / "_rels" / "document.xml.rels")
targets = {
    node.attrib["Id"]: node.attrib["Target"]
    for node in rels.getroot().findall(f"{REL_NS}Relationship")
}

for _event, node in ET.iterparse(WORD / "document.xml", events=("end",)):
    if node.tag != f"{WP_NS}inline":
        continue
    prop = node.find(f"{WP_NS}docPr")
    blip = node.find(f".//{A_NS}blip")
    if prop is not None and blip is not None:
        relation = blip.attrib.get(f"{R_NS}embed", "")
        print(f"{prop.attrib.get('descr', '')}\t{relation}\t{targets.get(relation, '')}")
    node.clear()

print("\n--- Figure S31 legend context ---")
W_NS = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
paragraphs = []
for _event, node in ET.iterparse(WORD / "document.xml", events=("end",)):
    if node.tag == f"{W_NS}p":
        value = "".join(text.text or "" for text in node.iter(f"{W_NS}t")).strip()
        if value:
            paragraphs.append(value)
        node.clear()
indices = [index for index, value in enumerate(paragraphs) if "Figure S31" in value]
if indices:
    for value in paragraphs[indices[0]:indices[0] + 8]:
        print(value)
