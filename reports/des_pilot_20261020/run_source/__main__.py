"""Command-line entry point."""
import argparse
import logging
import sys

def main():
    parser=argparse.ArgumentParser(description='Algorithmic morphology CAS pilot')
    sub=parser.add_subparsers(dest='command',required=True)
    run=sub.add_parser('run',help='Process catalog and generate bound HTML report')
    run.add_argument('--config-file-name',required=True)
    run.add_argument('--limit',type=int)
    run.add_argument('--chunk-size',type=int,default=50)
    args=parser.parse_args()
    logging.basicConfig(level=logging.INFO,stream=sys.stdout,
        format='%(asctime)sZ %(levelname)s %(name)s %(message)s',datefmt='%Y-%m-%dT%H:%M:%S')
    logging.Formatter.converter=__import__('time').gmtime
    if args.command=='run':
        from .runner import run as execute
        print(execute(args.config_file_name,limit=args.limit,chunk_size=args.chunk_size))

if __name__=='__main__': main()
