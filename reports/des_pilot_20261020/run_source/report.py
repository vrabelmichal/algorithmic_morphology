"""Self-contained, run-bound exploratory HTML report; never infer input from latest."""
from __future__ import annotations

import base64
import hashlib
from html import escape
from io import BytesIO
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .metrics import PROFILES


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(2 ** 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _picture(frame, series, label, *, against=None):
    fig, ax = plt.subplots(figsize=(6.4, 3.7), constrained_layout=True)
    if against:
        both = frame[[series, against]].apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        ax.scatter(both[against], both[series], s=7, alpha=0.35)
        ax.set_xlabel(against); ax.set_ylabel(series)
        ax.set_title(f"{label}: n={len(both)}")
        if len(both) > 0:
            lo = float(both.min().min()); hi = float(both.max().max())
            ax.plot([lo,hi],[lo,hi],ls='--',lw=1,color='grey')
    else:
        values = pd.to_numeric(frame[series], errors="coerce").replace([np.inf,-np.inf],np.nan).dropna()
        ax.hist(values, bins=min(40,max(5,int(np.sqrt(len(values))) if len(values) else 5)), color='#617e8e')
        ax.set_xlabel(series);ax.set_ylabel('Galaxies')
        ax.set_title(f"{label}: n={len(values)} / {len(frame)}")
    ax.grid(alpha=.2)
    with BytesIO() as buffer:
        fig.savefig(buffer, format="png", dpi=110)
        plt.close(fig)
        return base64.b64encode(buffer.getvalue()).decode("ascii")


def summarize(frame):
    """Finite-only per-band statistics, retaining denominators and flag counts."""
    metrics = ['C','rpetro','A_raw','A_sky','A_corrected','mask_fraction']
    metrics += [f'S_{p.name}_{suffix}' for p in PROFILES for suffix in ['raw','sky','corrected']]
    stats = []
    for band, part in frame.groupby('band',sort=True):
        for metric in metrics:
            vals = pd.to_numeric(part.get(metric,pd.Series(dtype=float)),errors='coerce')
            vals = vals.replace([np.inf,-np.inf],np.nan).dropna()
            stats.append(dict(band=str(band),metric=metric,total=len(part),valid=len(vals),
                              negative=int((vals<0).sum()),median=float(vals.median()) if len(vals) else None,
                              p16=float(vals.quantile(.16)) if len(vals) else None,
                              p84=float(vals.quantile(.84)) if len(vals) else None))
    flags = frame.get('quality_flags',pd.Series(dtype=str)).dropna().str.split('|').explode()
    flags = flags[flags.ne('')].value_counts().to_dict()
    errors = frame.get('error_message',pd.Series(dtype=str)).dropna().value_counts().to_dict()
    return {'statistics':stats,'quality_flags':{str(k):int(v) for k,v in flags.items()},
            'errors':{str(k):int(v) for k,v in errors.items()}}


def render_report(table_path, run_id, config_hash, *, inputs=None):
    """Return HTML string; caller publishes atomically with this table digest."""
    df = pd.read_parquet(table_path)
    digest = sha256_file(table_path)
    summary = summarize(df)
    status = df.get("status",pd.Series(dtype="str")).value_counts().to_dict()
    plots = []
    keys = ["C", "A_raw", "A_corrected"]
    for p in PROFILES:
        keys.extend([f"S_{p.name}_raw",f"S_{p.name}_corrected"])
    for band, part in df.groupby("band",sort=True):
        for key in keys:
            if key in part:
                plots.append((f"{band}: {key}",_picture(part,key,str(band))))
        for y,x in [("A_corrected","A_raw"),
                    (f"S_{PROFILES[1].name}_raw",f"S_{PROFILES[0].name}_raw"),
                    (f"S_{PROFILES[1].name}_corrected",f"S_{PROFILES[0].name}_corrected")]:
            if y in part and x in part:
                plots.append((f"{band}: {y} vs {x}",_picture(part,y,str(band),against=x)))
    rows = []
    for band, part in df.groupby("band",sort=True):
        def count_numeric(key):
            if key not in part:
                return 0
            return int(pd.to_numeric(part[key],errors="coerce").replace([np.inf,-np.inf],np.nan).notna().sum())
        valid_raw = count_numeric('A_raw')
        valid_corr = count_numeric('A_corrected')
        ss = [count_numeric(f'S_{p.name}_corrected') for p in PROFILES]
        rows.append(f"<tr><td>{escape(str(band))}</td><td>{len(part)}</td><td>{valid_raw}</td><td>{valid_corr}</td><td>{ss[0]}</td><td>{ss[1]}</td></tr>")
    stats_table = pd.DataFrame(summary['statistics']).to_html(index=False,float_format=lambda v:f'{v:.5g}',na_rep='missing')
    flag_table = escape(json.dumps({'quality_flags':summary['quality_flags'],'errors':summary['errors']},indent=2))
    title = 'Algorithmic morphology — profile comparison'
    status_display = escape(json.dumps({str(k):int(v) for k,v in status.items()},sort_keys=True))
    figures = ''.join(f'<figure><img src="data:image/png;base64,{src}" alt="{escape(label)}"><figcaption>{escape(label)}</figcaption></figure>' for label,src in plots)
    metadata = json.dumps({"run_id":run_id,"parquet_file":table_path.name,"parquet_sha256":digest,
                           "config_sha256":config_hash,"inputs":inputs or {}},sort_keys=True,default=str,indent=2)
    return f'''<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>{escape(title)}</title>
<style>body{{font:16px system-ui,sans-serif;max-width:1150px;margin:2rem auto;padding:0 1rem;color:#252525;line-height:1.5}}figure{{display:inline-block;vertical-align:top;width:min(100%,540px);margin:1rem .4rem}}img{{width:100%}}figcaption{{font-size:.85rem}}table{{border-collapse:collapse}}td,th{{padding:.5rem;border:1px solid #aaa}}pre{{overflow:auto;background:#f1f1f1;padding:1rem;white-space:pre-wrap}}.note{{border-left:3px solid #777;padding-left:1rem}}</style></head>
<body><h1>{escape(title)}</h1><p><strong>Run:</strong> {escape(run_id)}<br><strong>Parquet:</strong> {escape(table_path.name)}<br><strong>SHA-256:</strong> <code>{digest}</code><br><strong>Rows:</strong> {len(df)}<br><strong>Statuses:</strong> {status_display}</p>
<p class="note">Exploratory comparison, not a benchmark or proof of physical accuracy. Both clumpiness profiles differ in smoothing, central exclusion and ×10 scaling: their numerical values are not directly interchangeable. Corrected A/S require a clean sky tile; a missing correction is shown as missing, not as raw or zero. The surviving valid subset can be biased.</p>
<h2>Valid measurement counts by band</h2><table><thead><tr><th>Band</th><th>Rows</th><th>A raw</th><th>A corrected</th><th>S Lotz corrected</th><th>S Conselice corrected</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
<h2>Quality flags and failures (all rows retained)</h2><pre>{flag_table}</pre>
<h2>Finite metric summaries</h2><p>p16–p84 describes the measured sample, not uncertainty on individual galaxies.</p>{stats_table}
<h2>Distributions and direct comparisons</h2>{figures}
<h2>Exact input / output provenance</h2><pre>{escape(metadata)}</pre>
</body></html>'''
