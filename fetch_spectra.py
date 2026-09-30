# Download excitation and emission spectra from FPbase (https://www.fpbase.org)
# and save them to data/spectra.csv. Run once; the app only reads the CSV.
import json
import time
import urllib.error
import urllib.request
import pandas as pd

URL = "https://www.fpbase.org/graphql/"

# name shown in the app : name on FPbase
DYES = {
    "BUV395": "BD Horizon BUV395", "BUV496": "BD Horizon BUV496", "BUV563": "BD Horizon BUV563",
    "BUV661": "BD Horizon BUV661", "BUV737": "BD Horizon BUV737", "BUV805": "BD Horizon BUV805",
    "DAPI": "DAPI", "Hoechst 33342": "Hoechst 33342", "Zombie UV": "Zombie UV",
    "LIVE/DEAD Blue": "LIVE/DEAD Fixable Blue",
    "BV421": "Brilliant Violet 421", "BV480": "Brilliant Violet 480", "BV510": "Brilliant Violet 510",
    "BV570": "Brilliant Violet 570", "BV605": "Brilliant Violet 605", "BV650": "Brilliant Violet 650",
    "BV711": "Brilliant Violet 711", "BV750": "Brilliant Violet 750", "BV785": "Brilliant Violet 785",
    "Pacific Blue": "Pacific Blue", "eFluor 450": "eFluor 450", "V450": "BD Horizon V450",
    "V500": "BD Horizon V500", "Super Bright 436": "Super Bright 436", "Super Bright 600": "Super Bright 600",
    "Super Bright 645": "Super Bright 645", "Super Bright 702": "Super Bright 702",
    "Super Bright 780": "Super Bright 780", "Spark Violet 538": "Spark Violet 538",
    "Zombie Violet": "Zombie Violet", "Zombie Aqua": "Zombie Aqua", "LIVE/DEAD Aqua": "LIVE/DEAD Fixable Aqua",
    "FITC": "Fluorescein (FITC)", "Alexa Fluor 488": "Alexa Fluor 488", "BB515": "BD Horizon BB515",
    "BB700": "BD Horizon - BB700", "PerCP": "PerCP", "PerCP-Cy5.5": "PerCP-Cy5.5",
    "PerCP-eFluor 710": "PerCP-eFluor 710", "GFP (EGFP)": "EGFP", "Spark Blue 550": "Spark Blue 550",
    "Zombie Green": "Zombie Green",
    "PE": "PE (R-PE / R-phycoerythrin)", "PE-CF594": "PE-CF594", "PE-Dazzle 594": "PE/Dazzle 594",
    "PE-Texas Red": "PE-Texas Red", "PE-eFluor 610": "PE-eFluor 610", "PE-Cy5": "PE-Cy5",
    "PE-Cy5.5": "PE-Cy5.5", "PE-Cy7": "PE-Cy7", "PE-Fire 780": "PE/Fire 780",
    "PE-Alexa Fluor 700": "PE-Alexa Fluor 700", "7-AAD": "7-AAD", "PI": "Propidium Iodide",
    "mCherry": "mCherry", "tdTomato": "tdTomato", "Zombie Red": "Zombie Red", "Zombie Yellow": "Zombie Yellow",
    "APC": "APC (allophycocyanin)", "Alexa Fluor 647": "Alexa Fluor 647", "Alexa Fluor 700": "Alexa Fluor 700",
    "APC-Cy7": "APC/Cy7", "APC-Fire 750": "APC/Fire-750", "APC-eFluor 780": "APC-eFluor 780",
    "APC-H7": "APC/H7", "eFluor 660": "eFluor 660", "RealBlue 670": "RealBlue 670",
    "Zombie NIR": "Zombie NIR", "LIVE/DEAD Near-IR": "LIVE/DEAD Fixable Near-IR",
    "FVD eFluor 780": "Fixable Viability Dye eFluor 780",
}


def gql(query):
    req = urllib.request.Request(URL, data=json.dumps({"query": query}).encode(),
                                 headers={"Content-Type": "application/json", "User-Agent": "flow-panel-picker"})
    for wait in [0, 10, 30, 60, 120]:
        time.sleep(wait)
        try:
            return json.load(urllib.request.urlopen(req))["data"]
        except urllib.error.HTTPError as e:
            if e.code != 429:
                raise
            print("rate limited, waiting...")
    raise RuntimeError("FPbase kept rate limiting, try again later")


spectra = gql("{ spectra { id subtype owner { name } } }")["spectra"]

# find the excitation and emission spectrum id of each dye
wanted = []
for short, fp_name in DYES.items():
    ids = {s["subtype"]: s["id"] for s in spectra if s["owner"] and s["owner"]["name"] == fp_name}
    ex_id = ids.get("EX") or ids.get("AB")
    em_id = ids.get("EM")
    if ex_id is None or em_id is None:
        print("missing:", short, ids)
        continue
    wanted += [(short, "ex", ex_id), (short, "em", em_id)]

# download them 20 at a time
rows = []
for start in range(0, len(wanted), 20):
    batch = wanted[start:start + 20]
    query = " ".join(f"s{i}: spectrum(id: {sid}) {{ data }}" for i, (_, _, sid) in enumerate(batch))
    result = gql("{ " + query + " }")
    for i, (short, kind, _) in enumerate(batch):
        for wl, val in result[f"s{i}"]["data"]:
            rows.append({"fluor": short, "kind": kind, "wavelength": wl, "value": val})
    time.sleep(3)

pd.DataFrame(rows).to_csv("data/spectra.csv", index=False)
print("saved", len(set(r["fluor"] for r in rows)), "dyes")
