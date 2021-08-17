"""
    This script will fetch data from online services and save it locally in the cache directory.

    Current sources:
        * Covid Act Now (covidactnow.org), for which you will need a free API key
        * Covidestim.org

    You will need a "settings.json" file in this directory with your Covid Act Now API key

"""

import json
import os
import pickle
from typing import List, Dict
from datetime import date as Date

import requests
import settings
from states import state_abbrevs


def main():
    if not os.path.exists(settings.cache_folder):
        os.makedirs(settings.cache_folder)

    fetch_covid_act_now()
    fetch_covid_estim()
    extract_variant_data()


def fetch_covid_estim():
    print("Fetching covidestim.org data")
    url = "https://covidestim.s3.us-east-2.amazonaws.com/latest/state/estimates.csv"
    output_path = os.path.join(settings.cache_folder, "covid-estim.csv")
    with requests.get(url, stream=True) as result:
        with open(output_path, "wb") as handle:
            for chunk in result.iter_content(chunk_size=1024 ** 2):
                handle.write(chunk)


def fetch_covid_act_now():
    for state in state_abbrevs():
        print(f"Fetching covidactnow.org data for {state}")
        url = f"https://api.covidactnow.org/v2/state/{state}.timeseries.json?apiKey={settings.act_covid_key}"
        response = requests.get(url)
        data = response.json()

        output_path = os.path.join(settings.cache_folder, f"covid_act_now-{state}.json")
        with open(output_path, "w") as handle:
            json.dump(data, handle)


def extract_variant_data():
    """
    The best source of variant data I've found so far is https://covid.cdc.gov/covid-data-tracker/#variant-proportions,
    which has moderately detailed variant proportions for different regions. Unfortunately I don't see any easy way
    to fetch the data automatically.

    If you wish to have the improved, updated data, you will need to go to that website, go to the first dashboard, and
    select the button at the bottom of the dashboard which says "Download" and then select "Tableau Workbook". Use the
    most recent version and save the file in the "cache" folder.

    It should be named "Variant_Proportions_Plus_Nowcasting.twbx"
    """
    expected_path = os.path.join(settings.cache_folder, "Variant_Proportions_Weekly.twbx")
    if not os.path.exists(expected_path):
        print("No tableau file for the variant proportions, if you're interested in using updated variant data"
              "please follow the instructions in the fetch_data.py script source code.")
        return
    print(f"Extracting variant data from {expected_path}...")

    from zipfile import ZipFile

    extract_path = os.path.join(settings.cache_folder, "variant_proportions.hyper")
    with ZipFile(expected_path, "r") as handle:
        hyper_files = [f for f in handle.infolist() if f.filename.endswith(".hyper")]
        if len(hyper_files) != 1:
            print(f" > error, wrong number of hyper files found ({len(hyper_files)} but was expecting 1)")
            return
        with handle.open(hyper_files[0], "r") as zip_file:
            with open(extract_path, "wb") as target_file:
                target_file.write(zip_file.read())

    if not os.path.exists(extract_path):
        print(" > error extracting hyper database, file not found")
        return

    try:
        variant_dictionary = _extract_variants_from_tableau(extract_path)
    except RuntimeError as e:
        return

    with open(os.path.join(settings.cache_folder, "variant_data.pickle"), "wb") as handle:
        pickle.dump(variant_dictionary, handle)


def _extract_variants_from_tableau(extract_path: str):
    from tableauhyperapi import HyperProcess, Telemetry, Connection, TableDefinition
    with HyperProcess(telemetry=Telemetry.DO_NOT_SEND_USAGE_DATA_TO_TABLEAU) as hyper:
        with Connection(endpoint=hyper.endpoint, database=extract_path) as connection:
            table_names = connection.catalog.get_table_names(schema="Extract")

            def find_table(name: str) -> TableDefinition:
                target_table = [t for t in table_names if str(t.name.unescaped).startswith(name)]
                if len(target_table) != 1:
                    print(f" > error, wrong number of {name} tables found ({len(target_table)} but was expecting 1)")
                    raise RuntimeError()
                return connection.catalog.get_table_definition(name=target_table[0])

            def get_with_keys(table: TableDefinition) -> List[dict]:
                keys = [column.name.unescaped for column in table.columns]
                rows = connection.execute_list_query(query=f"SELECT * FROM {table.table_name}")
                return [{k: v for k, v in zip(keys, row)} for row in rows]

            variant_info = get_with_keys(find_table("variant_definitions"))
            region_info = get_with_keys(find_table("hhs_regions"))
            proportion_info = get_with_keys(find_table("proportion_modeling"))

            delta_lineages = {v["lineage"] for v in variant_info if v["who_label"] == "Delta"}

            regional = {}

            region_info.append({"hhs_region": "USA"})
            for region in region_info:
                region_key = region["hhs_region"]
                proportions = [p for p in proportion_info if p["usa_or_hhsregion"] == region_key]
                weeks = {p["week_ending"] for p in proportions}
                weekly_data = []
                for week in sorted(weeks):
                    proportions_this_week = [p for p in proportions if p["week_ending"] == week]
                    variants = {p["variant"]: p["share"] for p in proportions_this_week}
                    delta = sum(v for k, v in variants.items() if k in delta_lineages)
                    weekly_data.append({"date": Date(week.year, week.month, week.day),
                                        "variants": {"alpha": 1-delta, "delta": delta}})

                regional["USA" if region_key == "USA" else region['state_abv']] = weekly_data

            return regional




if __name__ == '__main__':
    main()
