"""Test validate map"""
import yaml

from validate_map import *


def test__function_map():
    """Test that all functions referenced in the config file are listed in the function map."""
    config = read_config("config.yaml")
    fxn_map = create_function_map()

    fxn_config = []
    for check in config["check"]:
        fxn_config.append(config["check"][check]["function"])

    config_not_map = set(fxn_config) - set(fxn_map.keys())

    assert len(config_not_map) == 0
