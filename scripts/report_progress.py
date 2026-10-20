"""Write an explicitly provisional HTML status page for one immutable run ID."""
import argparse
from datetime import datetime, timezone
from html import escape
import json
from pathlib import Path
import pandas as pd


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir',type=Path,required=True)
    p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();state=json.loads((a.run_dir/'progress.json').read_text())
    coverage=pd.read_csv(a.inputs/'coverage.csv')
    parts=[a.run_dir/v['file'] for v in state['parts']]
    columns=['band','status','C','A_raw','A_corrected','Gini_aperture_v1','M20_aperture_v1']
    if parts:
        frame=pd.concat([pd.read_parquet(q).reindex(columns=columns) for q in parts],ignore_index=True)
        counts=frame.groupby('band')[columns[2:]].count().to_html()
        outcomes=frame.groupby(['band','status']).size().unstack(fill_value=0).to_html()
    else:counts='<p>No measurement partition has been committed yet.</p>';outcomes=''
    fraction=state['checkpointed']/state['count']
    html=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Full catalog morphology — processing status</title><style>body{{font:16px system-ui;max-width:1000px;margin:2rem auto;padding:0 1rem;line-height:1.55;color:#222}}table{{border-collapse:collapse;display:block;overflow:auto}}td,th{{padding:.5rem;border:1px solid #aaa}}progress{{width:100%}}code{{overflow-wrap:anywhere}}</style></head><body><h1>Full upstream catalog: processing status</h1><p><strong>Status: {escape(state['status'])}. This is a provisional progress page, not the completed analysis.</strong></p><p>84,200 catalog IDs × three bands = 252,600 requested rows. All upstream IDs are selected, without an availability cut.</p><progress value="{state['checkpointed']}" max="{state['count']}"></progress><p>{state['checkpointed']:,} / {state['count']:,} rows checkpointed ({fraction:.1%}). Error rows: {state['failed']:,}; invalid CAS rows: {state['invalid']:,}.</p><h2>Input availability for the complete catalog</h2>{coverage.to_html(index=False)}<h2>Measurements already committed</h2><p>These counts describe the completed prefix only. They must not be used as final catalog statistics.</p>{counts}{outcomes}<h2>Metrics being evaluated</h2><p>C; raw, sky and corrected A; both named raw, sky and corrected S profiles; aperture Gini and M20; positive-flux aperture axis ratio, ellipticity and image-coordinate position angle.</p><p>The pipeline automatically writes the full HTML analysis after the final table is committed. No user decision or scientific approval is needed to keep processing. Calibration limits remain explicit in the final report.</p><p>Run: <code>{escape(state['run_id'])}</code><br>Updated: {datetime.now(timezone.utc).isoformat()}</p></body></html>'''
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(html);print(a.output.resolve())
