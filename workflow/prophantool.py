#!/usr/bin/env python

import pathlib
from snakemake.cli import main
import argparse
import sys
import json
import yaml

import warnings
warnings.filterwarnings("ignore")

__version__ = '1.0.0'

parser = argparse.ArgumentParser(description="ProphAnTool's main script")
parser.add_argument("-s", "--snakefile", default=f'{sys.path[0]}/Snakefile', help="Path to Snakefile")
parser.add_argument(
    "-c", "--configfile", required=True,
    help="Configuration file for ProphAnTool (JSON)")
parser.add_argument(
    '--unlock', action='store_true', default=False,
    help='If user forced termination of workflow, this might be required')
parser.add_argument('-v', '--version', action='version', version=f'ProphAnTool {__version__}')
args = parser.parse_args()


def read_config(filename):
    if filename.split('.')[-1] == 'yaml':
        with open(filename) as stream:
            try:
                return yaml.safe_load(stream), 'yaml'
            except yaml.YAMLError as exc:
                print(exc)
    elif filename.split('.')[-1] == 'json':
        with open(filename) as f:
            return json.load(f), 'json'
    else:
        sys.exit('ERROR: Config file must end in either ".json" or ".yaml"')


def save_config(config_data, filename, output_format):
    with open(filename, 'w', encoding='utf-8') as f:
        if output_format == 'json':
            json.dump(config_data, f, ensure_ascii=False, indent=2)
        elif output_format == 'yaml':
            yaml.dump(config_data, f, indent=2)
        else:
            return NotImplementedError


user_config, config_format = read_config(args.configfile)
config = read_config(f'{sys.path[0]}/default_config.json')[0]       # default configurations
for key in config.keys():
    if key in user_config.keys():
        config[key] = user_config[key]                              # set default values
#validate_config(config)
pathlib.Path(config["output"]).mkdir(parents=True, exist_ok=True)
save_config(config, f'{config["output"]}/config.json', output_format=config_format)

command = (
    f"-s {args.snakefile} --printshellcmds --cores {config['threads']} --configfile {config['output']}/config.json "
    f"--use-conda{' --unlock' if args.unlock else ''}")

print(f"ProphAnTool command: snakemake {command}")
main(command)