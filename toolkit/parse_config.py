"""
This script is designed to parse the given configuration file.

"""

import configparser


def modify_json_from_config(json_obj_default, config_file):
    json_obj = json_obj_default
    if config_file is not None:
        new_json_obj = read_config(config_file, default_json=json_obj)
        return new_json_obj
    return json_obj

def read_config(config_file,default_json=None):
    config = configparser.ConfigParser()
    config.read(config_file)
    # if success
    json_obj = parse_config2json(config,default_json=default_json)
    return json_obj


def process_value(v):
    if v.isnumeric():
        return int(v)
    elif v.lower() == 'true':
        return True
    elif v.lower() == 'false':
        return False
    return v


def parse_config2json(config,default_json=None):
    if default_json is None:
        collect_dict = {_: {}
                    for _ in config.sections()}
    else:
        collect_dict = default_json

    for key in config.sections():
        sto_dict = collect_dict[key]
        section_obj = config[key]
        for _key, value in section_obj.items():
            value = process_value(value)
            for idx, step_key in enumerate(_key.split('.')):
                if step_key not in sto_dict:
                    sto_dict[step_key] = {}
                if not idx == _key.count('.'):
                    sto_dict = sto_dict[step_key]
            sto_dict[step_key] = value
            sto_dict = collect_dict[key]
    return collect_dict


if __name__ == '__main__':
    pass
