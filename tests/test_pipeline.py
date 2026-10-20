import hashlib
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip('pyarrow')
fits = pytest.importorskip('astropy.io.fits')

from algorithmic_morphology.runner import run


def test_pipeline_checkpoint_report_and_nonoverwrite(tmp_path):
    data=tmp_path/'images'; masks=tmp_path/'masks'
    data.mkdir(); masks.mkdir()
    y,x=np.indices((96,96))
    img=40*np.exp(-((x-47.5)**2+(y-47.5)**2)/(2*9**2))
    fits.writeto(data/'LSB_ID1.fits',img,overwrite=False)
    fits.writeto(masks/'mask1.fits',np.zeros((96,96),dtype=np.uint8),overwrite=False)
    fits.BinTableHDU.from_columns([fits.Column(name='COADD_OBJECT_ID',format='K',array=[1,2])]).writeto(tmp_path/'catalog.fits')
    c=tmp_path/'example.py'
    c.write_text(f'''from pathlib import Path
lsb_catalog_fits_path=Path({str(tmp_path/'catalog.fits')!r})
band_order=['g']
band_sky_subtracted_dir_paths={{'g':Path({str(data)!r})}}
band_mask_dir_paths={{'g':Path({str(masks)!r})}}
allow_geometric_center=True
morphology_output_dir_path=Path({str(tmp_path/'output')!r})
''')
    first=run(c,chunk_size=1)
    assert len(list((first/'parts').glob('*.parquet')))==2
    table=next(first.glob('measurements__*.parquet'))
    report=next(first.glob('report__*.html'))
    import pandas as pd
    df=pd.read_parquet(table)
    assert list(df['status'])==['ok','error']
    text=report.read_text()
    assert hashlib.sha256(table.read_bytes()).hexdigest() in text
    assert 'lotz2004_statmorph_like_v1' in text
    assert 'conselice2003_like_v1' in text
    assert (tmp_path/'output'/'latest').resolve()==first.resolve()
    second=run(c,chunk_size=1)
    assert second!=first
    assert table.exists() and report.exists()
    assert (tmp_path/'output'/'latest').resolve()==second.resolve()
