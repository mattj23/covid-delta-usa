"""
    This script will fetch data from online services and save it locally in the cache directory.

    Current sources:
        * Covid Act Now (covidactnow.org), for which you will need a free API key
        * Covidestim.org

    You will need a "settings.json" file in this directory with your Covid Act Now API key

"""

import json
import os
import requests
import settings
from states import state_abbrevs


def main():
    if not os.path.exists(settings.cache_folder):
        os.makedirs(settings.cache_folder)

    fetch_covid_act_now()
    fetch_covid_estim()


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


if __name__ == '__main__':
    main()
