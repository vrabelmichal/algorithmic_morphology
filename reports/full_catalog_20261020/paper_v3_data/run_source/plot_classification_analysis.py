"""Render class/probability figures from frozen matched inputs."""
import argparse,json
from pathlib import Path
import pandas as pd
from algorithmic_morphology.classification_plots import create_plots

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);a=p.parse_args()
contract=json.loads((a.input/'feature_contract.json').read_text())
figures=create_plots(pd.read_parquet(a.input/'joined.parquet'),pd.read_csv(a.input/'associations.csv'),contract,a.input/'figures')
(a.input/'figures.json').write_text(json.dumps(figures,indent=2));print(f'Created {len(figures)} classification figures',flush=True)
