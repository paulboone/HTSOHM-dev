import os
import yaml

def load_config_file(file_name):
    """Reads input file.

    Args:
        file_name (str): auto-created config filename (ex. htsohm.sample.yaml).

    Returns:
        config (dict): parameters specified in config.

    """
    with open(file_name) as config_file:
         config = yaml.load(config_file)

    if not 'override_restart_errors' in config:
        config['override_restart_errors'] = False
    if not 'keep_configs' in config:
        config['keep_configs'] = False
    if not "output_dir" in config or not config["output_dir"]:
        config["output_dir"] = os.getcwd()
    else:
        os.makedirs(config['output_dir'], exist_ok=True)
    return config
