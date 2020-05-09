import os
import sys
from os.path import *

import click

sys.path.insert(0, dirname(dirname(__file__)))
from toolkit.parse_config import modify_json_from_config
from toolkit import IO_json, process_path

def main(input_json, input_config, output_file):
    output_file = process_path(output_file)
    if not exists(dirname(output_file)):
        os.makedirs(dirname(output_file))
    ori_json_obj = dict(IO_json(file=input_json, way='r'))
    new_json_obj = modify_json_from_config(ori_json_obj, input_config)
    IO_json(json_obj=new_json_obj,
            file=output_file,
            way='w')

@click.command()
@click.option("-i", "input_json", help="input json file")
@click.option("-c", "input_config", help="input configuration")
@click.option("-o", "output_file", help="output json file")
def cli(input_json, input_config, output_file):
    main(input_json, input_config, output_file)


if __name__ == '__main__':
    cli()
