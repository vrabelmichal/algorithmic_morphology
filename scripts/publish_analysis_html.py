"""Combine the written pilot analysis and immutable run plots in one HTML file."""
import argparse
import base64
from html import escape
import mimetypes
from pathlib import Path
import re

import markdown


def publish(report_dir):
    report_dir = Path(report_dir)
    narrative = markdown.markdown((report_dir/'README.md').read_text(), extensions=['fenced_code','tables'])
    def inline_image(match):
        source = (report_dir/match.group(1)).resolve()
        if not source.is_relative_to(report_dir.resolve()):
            raise ValueError('Image must be inside report directory')
        mime = mimetypes.guess_type(source.name)[0]
        data = base64.b64encode(source.read_bytes()).decode('ascii')
        return f'src="data:{mime};base64,{data}"'
    narrative = re.sub(r'src="([^"]+)"',inline_image,narrative)
    reports = list(report_dir.glob('report__*.html'))
    if len(reports) != 1:
        raise ValueError('Expected exactly one immutable run report')
    plots = reports[0].read_text().split('<body>',1)[1].rsplit('</body>',1)[0]
    style = '''body{font:16px system-ui,sans-serif;max-width:1150px;margin:2rem auto;padding:0 1rem;line-height:1.55;color:#252525}img{max-width:100%;height:auto}figure{display:inline-block;vertical-align:top;width:min(100%,540px);margin:1rem .4rem}figcaption{font-size:.85rem}table{border-collapse:collapse;display:block;overflow:auto}td,th{padding:.5rem;border:1px solid #aaa}pre{overflow:auto;background:#f1f1f1;padding:1rem;white-space:pre-wrap}code{overflow-wrap:anywhere}nav{padding:1rem;background:#edf3f6}hr{margin:3rem 0}'''
    html = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>DES CAS pilot — analysis and measurements</title><style>{style}</style></head><body><nav><a href="#analysis">Written analysis</a> · <a href="#measurements">Full measurement plots and provenance</a></nav><main id="analysis">{narrative}</main><hr><section id="measurements">{plots}</section></body></html>'''
    output = report_dir/'index.html'
    output.write_text(html)
    print(output.resolve())
    return output


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report-dir',type=Path,required=True)
    publish(parser.parse_args().report_dir)
