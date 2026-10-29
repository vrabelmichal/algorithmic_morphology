"""Publish report v4: embed the reproducible association explorer in section 6.4."""
import argparse
from datetime import datetime,timezone
from html import escape
import json
from pathlib import Path
import subprocess,sys
from algorithmic_morphology.report import sha256_file

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report-dir',type=Path,default=Path('reports/full_catalog_20261020'));a=p.parse_args();root=a.report_dir
if (root/'report_v4.html').exists():raise FileExistsError('Preserve report v4; choose a new revision')
subprocess.run([sys.executable,str(Path(__file__).with_name('build_association_explorer.py')),'--report-dir',str(root)],check=True)
explorer=(root/'association_explorer.html').read_text()
# Isolate styles and global JavaScript names; srcdoc keeps the explorer available offline.
explorer=explorer.replace('</style>','body{padding:16px;box-sizing:border-box;background:#fff;color:#222}</style>')
explorer=explorer.replace('</body>',"<script>new ResizeObserver(()=>parent.postMessage({type:'morphology-explorer-height',height:Math.ceil(document.body.getBoundingClientRect().height)},'*')).observe(document.body);</script></body>")
section='<div id="correlation-explorer"><h3>6.4 Interactive correlation explorer</h3><p>Choose a human answer or DECaLS10 class, photometric band and comparison cohort. Each bar uses every eligible finite pair; counts and corrected/uncorrected quantities remain explicit.</p><p><a href="association_explorer.html">Open the explorer separately</a></p><iframe id="association-explorer-frame" title="Interactive morphology and classification rank correlations" style="display:block;width:100%;height:1200px;border:0;border-radius:8px;background:white" srcdoc="'+escape(explorer,quote=True)+'"></iframe></div>'
s=(root/'report_v3.html').read_text().replace('scientific analysis · revision 3','scientific analysis · revision 4')
s=s.replace('<h3>6.4 Correlation and class-separation atlas</h3>',section+'<h3>6.4a Correlation and class-separation atlas</h3>',1)
s=s.replace('<a href="#classification">6 · Human / DECaLS10</a>','<a href="#classification">6 · Human / DECaLS10</a><a href="#correlation-explorer">Interactive correlations</a>',1)
s=s.replace('</body>',"<script>const explorerFrame=document.getElementById('association-explorer-frame');function resizeExplorer(){try{explorerFrame.style.height=Math.ceil(explorerFrame.contentDocument.body.getBoundingClientRect().height)+16+'px';}catch(e){}}explorerFrame.addEventListener('load',resizeExplorer);window.addEventListener('resize',resizeExplorer);resizeExplorer();window.addEventListener('message',event=>{const frame=document.getElementById('association-explorer-frame');if(event.source===frame.contentWindow&&event.data?.type==='morphology-explorer-height'&&Number.isFinite(event.data.height)){frame.style.height=Math.max(500,Math.min(4000,event.data.height+16))+'px';}});</script></body>")
(root/'report_v4.html').write_text(s)
manifest=dict(report_version='scientific_report_v4',built_utc=datetime.now(timezone.utc).isoformat(),html_sha256=sha256_file(root/'report_v4.html'),
 previous_report_sha256=sha256_file(root/'report_v3.html'),explorer_sha256=sha256_file(root/'association_explorer.html'),
 association_table_sha256=sha256_file(root/'paper_v3_data/associations.csv'),morphology_sha256=sha256_file(root/'measurements_metrics.parquet'),
 source_code={str(Path(f).resolve()):sha256_file(f) for f in [__file__,Path(__file__).with_name('build_association_explorer.py')]},
 note='Presentation-only addition; report v3 data, figures and validation are preserved.')
(root/'report_v4_manifest.json').write_text(json.dumps(manifest,indent=2))
(root/'index.html.tmp').write_text(s);(root/'index.html.tmp').replace(root/'index.html')
history=json.loads((root/'report_history.json').read_text());history['v4']={'html':'report_v4.html','sha256':manifest['html_sha256'],'validation':'report_v4_validation.json'};history['default']='index.html';(root/'report_history.json').write_text(json.dumps(history,indent=2))
print('Explorer embedded in section 6.4 and linked in report navigation')
