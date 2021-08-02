"""
    Current workaround for not having a good, ingest-able source for covid variant data.  I collected these values on
    the delta variant by going to the CDC's website (https://covid.cdc.gov/covid-data-tracker/#variant-proportions) and
    typing in the weighted averages manually.  I'm also making the incorrect but hopefully functional assumption that
    all of the non-delta variants have, in aggregate, similar enough mechanics that it's not necessary to model them
    separate.  Thus the "alpha" variant here is a catch-all for everything that isn't delta.

"""
from typing import List, Tuple, Dict
from datetime import date as Date


def load_variant_history() -> List[Dict]:
    variants = [
        {"date": Date(2021, 5, 8), "variants": {"alpha": None, "delta": 0.014}},
        {"date": Date(2021, 5, 22), "variants": {"alpha": None, "delta": 0.031}},
        {"date": Date(2021, 6, 5), "variants": {"alpha": None, "delta": 0.103}},
        {"date": Date(2021, 6, 19), "variants": {"alpha": None, "delta": 0.319}},
        {"date": Date(2021, 7, 3), "variants": {"alpha": None, "delta": 0.627}},
        {"date": Date(2021, 7, 17), "variants": {"alpha": None, "delta": 0.822}},
    ]

    for row in variants:
        vs = row['variants']
        vs["alpha"] = 1.0 - vs["delta"]

    return variants
