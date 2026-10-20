from pathlib import Path
import hashlib,json,re,sys,xml.etree.ElementTree as ET
import numpy as np,pandas as pd,tifffile
from PIL import Image,ImageDraw,ImageFont
sys.stdout.reconfigure(encoding='utf-8')

ROOT=Path(r'data/external/local_4a28549fddd8')
SRC=ROOT/'20260930爬片'
OUT=ROOT/'统计分析_20261013'
for s in ['00_manifest','02_qc','03_per_cell','04_per_field','05_statistics','06_figures','07_excel','08_audit']:(OUT/s).mkdir(parents=True,exist_ok=True)
def dump(obj,p):p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def csv(rows,p):pd.DataFrame(rows).to_csv(p,index=False,encoding='utf-8-sig')
def local(tag):return tag.split('}')[-1]
rows=[]; detail=[]; protocols=[]; pdetail=[]; thumbs=[]
font=ImageFont.truetype('data/external/local_f014c5000f62',20)
for folder in sorted(SRC.iterdir()):
 if not folder.is_dir():continue
 print('Auditing',folder.name,flush=True)
 field=folder.name.split('.2026')[0]
 group,mag,rep=field.rsplit('-',2)
 for p in sorted(folder.iterdir()):
  b=p.read_bytes(); sha=hashlib.sha256(b).hexdigest()
  base=dict(field_id=field,filename=p.name,relative_path=str(p.relative_to(ROOT)),source_path=str(p),bytes=len(b),mtime=p.stat().st_mtime,sha256=sha)
  if p.suffix.lower()=='.scanprotocol':
   rows.append(dict(base,kind='protocol'))
   rt=ET.fromstring(b)
   for e in rt.iter():e.tag=local(e.tag)
   channels=[]
   for par in rt.iter('Channels'):
    channels.extend([dict(index=i,**{local(c.tag):c.text for c in e if len(c)==0},filter=e.findtext('./InstrumentChannelSettings/FilterConfiguration')) for i,e in enumerate(par) if local(e.tag)=='Channel'])
   acq=[]
   for par in rt.iter('PassChannelSettings'):
    for kv in par:
     ch=kv.findtext('Key');v=kv.find('Value')
     if v is None:continue
     vals={k:v.findtext('.//'+k) for k in ['ExposureTime','Gain','LightIntensity','AcquireAsGrayscale']}
     acq.append(dict(channel_index=ch,**vals))
   run=[]
   for i,e in enumerate(rt.iter('RunSequenceV2')):
    run.append(dict(pass_index=i,objective=e.findtext('.//Objective/PartNumber'),capture=[(c.findtext('ChannelIndex'),c.findtext('Capture')) for c in e.findall('./RunSequenceChannels/ChannelSequence')]))
   pdetail.append(dict(field_id=field,channels=channels,acquisition=acq,runs=run))
   protocols.append(dict(field_id=field,channels=json.dumps(channels,ensure_ascii=False),acquisition=json.dumps(acq,ensure_ascii=False)))
   continue
  if p.suffix.lower() not in ['.tif','.tiff']:continue
  print(' reading',p.name,flush=True)
  with tifffile.TiffFile(p) as t:
   pg=t.pages[0]
   with Image.open(p) as pim:a=np.array(pim)
   tags={x.name:str(x.value) for x in pg.tags.values()}
   kind=re.search('_Plate_([DRM])_',p.name)[1];ch=re.search(r'd(\d+)\.TIF',p.name)[1]
   row=dict(base,kind=kind,channel_file=ch,group_label=group,magnification_label=mag,field_number=int(rep),shape=str(a.shape),axes=pg.axes,dtype=str(a.dtype),pages=len(t.pages),min=float(a.min()),max=float(a.max()),mean=float(a.mean()),p01=float(np.percentile(a,1)),p50=float(np.percentile(a,50)),p99=float(np.percentile(a,99)),saturation_dtype_pct=float(100*np.mean(a==np.iinfo(a.dtype).max)))
   rows.append(row);detail.append(dict(field_id=field,filename=p.name,tags=tags))
   if kind=='M':
    im=Image.fromarray(a.astype('uint8'));im.thumbnail((400,310))
    tile=Image.new('RGB',(410,350),'#fafafa');tile.paste(im,((410-im.width)//2,35));ImageDraw.Draw(tile).text((8,5),field,fill='black',font=font);thumbs.append(tile)
csv(rows,OUT/'00_manifest/input_manifest.csv');csv(protocols,OUT/'00_manifest/protocol_summary.csv')
dump(detail,OUT/'00_manifest/tiff_tags.json');dump(pdetail,OUT/'00_manifest/protocol_details.json')
canvas=Image.new('RGB',(410*4,350*((len(thumbs)+3)//4)), 'white')
for i,im in enumerate(thumbs):canvas.paste(im,((i%4)*410,(i//4)*350))
canvas.save(OUT/'02_qc/all_fields_contact_sheet.jpg',quality=92)
print('OUT',OUT)
df=pd.DataFrame(rows);print(df.groupby(['group_label','magnification_label','kind'],dropna=True).size().to_string())
print('unique image specs',df[['kind','channel_file','shape','axes','dtype']].drop_duplicates().to_dict('records'))
print('first tags',json.dumps(detail[0],ensure_ascii=False)[:12000])
print('first protocol',json.dumps(pdetail[0],ensure_ascii=False))
