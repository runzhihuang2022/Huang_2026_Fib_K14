from pathlib import Path
import openpyxl,json,csv,hashlib
R=Path(__file__).resolve().parents[1];D=R/'data/processed/table_s1';O=R/'validation/table_s1';O.mkdir(parents=True,exist_ok=True)
p=D/'Table_S1_source.xlsx';wb=openpyxl.load_workbook(p,read_only=True,data_only=False);cached=openpyxl.load_workbook(p,read_only=True,data_only=True)
rows=[];errors=[]
for i,ws in enumerate(wb.worksheets,1):
 vals=list(ws.values);dest=R/'outputs/tables/table_s1'/f'sheet_{i:02}.csv';dest.parent.mkdir(parents=True,exist_ok=True)
 with dest.open('w',newline='',encoding='utf-8-sig') as f:csv.writer(f).writerows(vals)
 formula=0
 for row in ws:
  for c in row:
   if c.data_type=='e':errors.append({'sheet':ws.title,'cell':c.coordinate,'value':c.value})
   if c.data_type=='f':formula+=1
 for row in cached[ws.title]:
  for c in row:
   if c.data_type=='e':errors.append({'sheet':ws.title,'cell':c.coordinate,'cached_error':c.value})
 existing=D/f'sheet_{i:02}.csv'
 # Both exports are values/string representations from the preserved source.
 same=list(csv.reader(existing.open(encoding='utf-8-sig')))==list(csv.reader(dest.open(encoding='utf-8-sig')))
 rows.append(dict(sheet=ws.title,rows=ws.max_row,columns=ws.max_column,formula_cells=formula,export_matches_input_csv=same))
assert len(rows)==11 and not errors and all(x['export_matches_input_csv'] for x in rows),(rows,errors)
report=dict(workbook_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),sheets=rows,spreadsheet_error_cells=errors,scope='Read-only content audit and CSV re-export; source formatting preserved; no formulas added')
(O/'report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');print(json.dumps(report,ensure_ascii=False))
