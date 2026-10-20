"""Read-only audit of supplied dual-reporter PPT/AI; no RGB cell quantification."""
from pathlib import Path
import csv, hashlib, io, json, re, zipfile
import xml.etree.ElementTree as ET
from PIL import Image
import pymupdf

ROOT = Path(r'data/external/local_87055123a971')
OUT = ROOT / '双同源比例分析_20261002'
for folder in ['00_manifest', '01_extracted', '04_per_slide', '08_audit']:
    (OUT / folder).mkdir(exist_ok=True, parents=True)

def sha(b):
    return hashlib.sha256(b).hexdigest()

def table(path, rows, cols=None):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols or list(rows[0]), delimiter='\t' if path.suffix == '.tsv' else ',')
        w.writeheader()
        w.writerows(rows)

manifest = []
for p in sorted(ROOT.iterdir()):
    if not p.is_file() or p.name.startswith('~$'):
        continue
    b = p.read_bytes()
    manifest.append(dict(source_path=str(p), bytes=len(b), sha256=sha(b), format=p.suffix,
                         status='DISPLAY_SOURCE_ONLY_NO_CELL_COUNTS'))
table(OUT / '00_manifest' / 'input_manifest.tsv', manifest)

ns = {'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
      'p':'http://schemas.openxmlformats.org/presentationml/2006/main',
      'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
audit = []
with zipfile.ZipFile(ROOT / '20260930.pptx') as z:
    slides = sorted([n for n in z.namelist() if re.fullmatch(r'ppt/slides/slide\d+\.xml', n)],
                    key=lambda n:int(re.search(r'slide(\d+)',n).group(1)))
    for sl in slides:
        num = int(re.search(r'slide(\d+)', sl).group(1))
        xml = ET.fromstring(z.read(sl))
        text = ' | '.join(t.text or '' for t in xml.findall('.//a:t', ns))
        rels = {e.attrib['Id']:e.attrib['Target'] for e in ET.fromstring(z.read(f'ppt/slides/_rels/slide{num}.xml.rels'))}
        for blip in xml.findall('.//p:pic//a:blip', ns):
            target = rels[blip.attrib['{'+ns['r']+'}embed']]
            media = 'ppt/' + target.removeprefix('../')
            b = z.read(media)
            fname = Path(media).name
            (OUT / '01_extracted' / fname).write_bytes(b)
            im = Image.open(io.BytesIO(b))
            match = re.search(r'^D(\d+)', text, re.I)
            audit.append(dict(slide_number=num, source_label=text, embedded_file=media,
                              extracted_file=str(OUT/'01_extracted'/fname), width_px=im.width,
                              height_px=im.height, image_mode=im.mode, frame_count=getattr(im,'n_frames',1),
                              sha256=sha(b), candidate_day=int(match.group(1)) if match else '',
                              time_semantics='UNCONFIRMED', subject_id='', group='',
                              analysis_scope='displayed_section_image' if num < 9 else 'montage_with_repeated_zoom_insets',
                              status='BLOCKED_SOURCE_CHANNELS_AND_METADATA_MISSING'))
table(OUT / '00_manifest' / 'ppt_image_audit.tsv', audit)
table(OUT / '00_manifest' / 'sample_metadata_pending.csv', audit)

pending = []
for row in audit:
    if not row['source_label']:
        continue
    pending.append(dict(source_label=row['source_label'], candidate_day=row['candidate_day'],
                        subject_id='', group='', analysis_region='', numerator_double_positive_n='',
                        denominator_valid_DAPI_n='', ratio_pct='',
                        source_counts_or_channels='', status='PENDING_INPUT_NOT_ZERO'))
table(OUT / '04_per_slide' / '待补充细胞计数.csv', pending)

doc = pymupdf.open(ROOT / '未标题-1.ai')
ai_text = '\n'.join(page.get_text() for page in doc)
(OUT / '08_audit' / 'ai_visible_labels.txt').write_text(ai_text, encoding='utf-8')
ai_images = []
for page_idx, page in enumerate(doc):
    for item in page.get_images(full=True):
        info = doc.extract_image(item[0])
        ai_images.append(dict(page=page_idx+1, xref=item[0], width_px=item[2], height_px=item[3],
                              bits_per_component=item[4], colorspace=info.get('cs-name'),
                              purpose='display_raster_no_source_channel_metadata'))
table(OUT / '00_manifest' / 'ai_image_audit.tsv', ai_images)

notes = {
    'status':'BLOCKED_NOT_QUANTIFIED', 'ppt_images':len(audit), 'labeled_section_images':len(pending),
    'reason':'Only 8-bit RGB display images and a montage; no raw marker channels or cell counts.',
    'requested_endpoint':'same-cell Pdgfra-tdT and Krt14-GFP double-positive fraction over time',
    'default_denominator_pending_confirmation':'valid DAPI-positive cells',
    'formula':'100 * double_positive_cell_count / valid_denominator_cell_count',
    'replicate_unit':'animal; fields and sections nested; identity not confirmed',
    'candidate_times':[3,5,7,12], 'candidate_times_confirmed':False,
    'label_conflicts':['PPT D5-G2-y2 vs AI d5-g2-y21', 'PPT D12-G-y2 vs existing review crop d12-g2-y2'],
    'normal_label_interpretation':'Nor is unconfirmed and must not be silently pooled with wound samples',
    'no_numeric_results':True, 'line_plot_created':False,
    'required_next_input':['raw registered channels or traceable cell counts', 'time/group/animal mapping',
                           'denominator and valid analysis region', 'positivity/threshold validation']}
(OUT / '08_audit' / 'analysis_status.json').write_text(json.dumps(notes, ensure_ascii=False, indent=2), encoding='utf-8')
animal_ppt = ROOT.parent / '动物实验' / '20260908双同源动物图组图D0-D12.pptx'
with zipfile.ZipFile(animal_ppt) as animal_zip:
    schedule_text = [t.text for t in ET.fromstring(animal_zip.read('ppt/slides/slide1.xml')).findall('.//a:t', ns)]
(OUT / '08_audit' / 'animal_sampling_schedule.json').write_text(json.dumps({
    'source':str(animal_ppt), 'sha256':sha(animal_ppt.read_bytes()), 'slide1_text':schedule_text,
    'interpretation':'Experiment record lists sampling D3, D5, D7, D12 after modeling; stained-image linkage, groups and animal identities still require confirmation.'
}, ensure_ascii=False, indent=2), encoding='utf-8')
readme = '''# 双同源标记细胞比例：输入核查（2026-10-02）

状态：待补充原始输入；尚未计算细胞比例，尚未生成数值折线图。

## 当前文件

原目录有 20260930.pptx、未标题-1.ai 和 JPG 拼图。PPT 提取出 8 张有编号的组织图、1 张无编号的拼图；均为单帧、8 位 RGB 展示图，无独立荧光通道。AI 的 PDF 兼容内容也只包含展示用栅格图，没有细胞计数表。所有原始文件保持不变。

PPT 标签：D3-G6-Nor、D3-G6-y2、D5-G2-y2、D5-G2-z1、D5-G4-y1、D7-G5-nor、D12-G5-m2、D12-G-y2。邻近动物实验目录的原始取材记录列出造模后 D3、D5、D7、D12 取材，与染色标签的时间编码相符；逐张染色图与该取材记录的关联仍待确认。

AI/JPG 拼图标签与 PPT 有差异，并包含 D7-G3-Z2；拼图的放大插图重复展示同一区域，不能增加细胞分母或独立样本数。D3/D5/D7/D12 暂记为候选时间编码，未确认伤后天数。Nor、G 编号和 y/z/m 后缀的含义待确认；同一 G 编号不可自行认定为不同动物。

## 为什么暂缓比例和折线图

双阳性必须归属同一个细胞。RGB 合成图的黄色像素受显示映射、亮度及相邻信号重叠影响，不能直接转换为双标细胞数。当前没有原始通道强度、阳性阈值依据、有效细胞分母或核/细胞分割核查。缺失数值在 CSV 中留空，空白不代表 0%。

## 补充后计算

请提供同一切片注册的 DAPI、tdT、GFP 原始通道文件，或已经逐视野核查的细胞计数表，并确认动物编号、时间/组别、计数范围及分母。

默认待确认的比例 = 同一细胞 tdT+ 且 GFP+ 的数量 / 分析范围内有效 DAPI+ 细胞数 × 100%。若要求以 tdT+ 或 GFP+ 细胞作分母，请明确。

每只动物多视野的主要比例先汇总分子和分母再相除，不能把视野当作独立动物。折线图以独立动物为点、时间为横轴；只有重复设计和独立 n 明确后才添加误差线或统计检验。正常对照是否独立展示须按元数据确认。

## 文件入口

- 00_manifest/input_manifest.tsv：输入文件 SHA-256。
- 00_manifest/ppt_image_audit.tsv：页码、标签、尺寸、格式、提取路径。
- 00_manifest/sample_metadata_pending.csv：待确认元数据。
- 04_per_slide/待补充细胞计数.csv：待补充计数；未填任何估算值。
- 01_extracted/：PPT 内嵌图像的逐字节副本。
- source_contact_sheet.jpg：仅用于图像核查的缩略图。
- 08_audit/analysis_status.json：当前完成范围及缺失项。
'''
(OUT / 'README_输入核查与待补充项.md').write_text(readme, encoding='utf-8')
assert len(audit) == 9 and len(pending) == 8
assert all(r['image_mode'] == 'RGB' and r['frame_count'] == 1 for r in audit)
assert all(r['ratio_pct'] == '' for r in pending)
for row in audit:
    assert sha(Path(row['extracted_file']).read_bytes()) == row['sha256']
print(json.dumps({'output':str(OUT), 'status':notes['status'], 'ppt_images':len(audit), 'pending_count_rows':len(pending)}, ensure_ascii=False))
