"""
    Contains settings related to the simulator, the data, and the repository
"""
import os
import json

act_covid_key = None
binary_path = None
cache_folder = None
default_input_file = None
default_output_file = None


def _load_credentials():
    global act_covid_key
    global default_input_file
    global default_output_file
    global binary_path
    global cache_folder
    settings_file = "settings.gitignore.json" if os.path.exists("settings.gitignore.json") else "settings.json"
    with open(settings_file, "r") as handle:
        settings_dict = json.load(handle)
        act_covid_key = settings_dict["act_covid_key"]
        default_input_file = settings_dict["default_input_file"]
        default_output_file = settings_dict["default_output_file"]
        binary_path = settings_dict["binary_path"]
        cache_folder = settings_dict["cache_folder"]


_load_credentials()
