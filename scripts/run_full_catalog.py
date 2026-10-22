"""Run the whole catalog and publish its HTML analysis automatically on completion."""
import argparse
import logging
from pathlib import Path
import subprocess
import sys
import time
from algorithmic_morphology.runner import run


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--source-project',type=Path,default=Path('../nucleated_lsbg'))
    a=p.parse_args()
    logging.basicConfig(level=logging.INFO,stream=sys.stdout,
                        format='%(asctime)sZ %(levelname)s %(message)s',datefmt='%Y-%m-%dT%H:%M:%S')
    logging.Formatter.converter=time.gmtime
    run_dir=run(a.config,chunk_size=1000)
    subprocess.run([sys.executable,str(Path(__file__).with_name('report_full_catalog.py')),
                    '--run-dir',str(run_dir),'--inputs',str(a.inputs),'--output',str(a.output)],check=True)
    subprocess.run([sys.executable,str(Path(__file__).with_name('build_paper_report.py')),
                    '--report-dir',str(a.output),'--source-project',str(a.source_project)],check=True)
    subprocess.run([sys.executable,str(Path(__file__).with_name('validate_paper_report.py')),
                    '--report-dir',str(a.output)],check=True)
    logging.info('Full-catalog HTML published: %s',a.output/'index.html')
