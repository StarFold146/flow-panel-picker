"""Spectral signatures, similarity and panel search for spectral flow cytometry."""
import heapq
import itertools

import numpy as np
import pandas as pd

# BD FACSDiscover S8 lasers (nm)
LASERS = {"UV": 349, "V": 405, "B": 488, "YG": 561, "R": 637}

# Detector bins are approximate: 20 nm wide, from just above each laser to 850 nm.
BIN_WIDTH = 20
LAST_NM = 850


# Rough relative brightness of each dye on antibodies, 1 (dim) to 5 (brightest).
# Taken loosely from vendor brightness charts; real brightness depends on the clone,
# the antigen and the instrument, so use it for ranking only.
BRIGHTNESS = {
    "PE": 5, "BV421": 5, "PE-Cy5": 5, "APC": 4, "BB515": 4, "BB700": 4,
    "BUV395": 4, "BUV563": 4, "BUV661": 4, "BUV737": 4, "RealBlue 670": 4,
    "BV605": 4, "BV650": 4, "BV711": 4, "PE-CF594": 4, "PE-Dazzle 594": 4,
    "PE-eFluor 610": 4, "PE-Cy7": 4, "PE-Fire 780": 4,
    "Super Bright 436": 4, "Super Bright 600": 4, "Super Bright 645": 4,
    "Alexa Fluor 647": 3, "Alexa Fluor 488": 3, "APC-Cy7": 3, "APC-Fire 750": 3,
    "BUV496": 3, "BV480": 3, "BV750": 3, "BV785": 3, "PE-Cy5.5": 3,
    "PE-Texas Red": 3, "PE-Alexa Fluor 700": 3, "PerCP-eFluor 710": 3,
    "Spark Blue 550": 3, "Spark Violet 538": 3, "Super Bright 702": 3, "eFluor 660": 3,
    "APC-H7": 2, "APC-eFluor 780": 2, "BUV805": 2, "BV510": 2, "BV570": 2,
    "FITC": 2, "PerCP-Cy5.5": 2, "V450": 2, "eFluor 450": 2, "Super Bright 780": 2,
    "PerCP": 1, "Pacific Blue": 1, "V500": 1, "Alexa Fluor 700": 1,
    "Kiravia Blue 520": 4, "PE-Alexa Fluor 610": 4,
    "Alexa Fluor 532": 3, "Alexa Fluor 546": 3, "Alexa Fluor 555": 3, "Alexa Fluor 568": 3,
    "Alexa Fluor 594": 3, "APC-Alexa Fluor 750": 3, "eFluor 570": 3, "eFluor 615": 3,
    "DyLight 488": 3, "DyLight 550": 3, "DyLight 650": 3, "Cy3": 3, "Cy5": 3,
    "iFluor 488": 3, "iFluor 647": 3,
    "Alexa Fluor 680": 2, "Alexa Fluor 750": 2, "Pacific Green": 2, "eFluor 710": 2,
    "DyLight 680": 2, "DyLight 755": 2, "Cy5.5": 2,
    "Alexa Fluor 405": 1, "Alexa Fluor 430": 1, "Pacific Orange": 1, "eFluor 506": 1,
    "DyLight 405": 1, "Cy7": 1,
}

# Dyes that stain dead cells or DNA, and fluorescent proteins: never suggested for an antibody.
VIABILITY = {"7-AAD", "DAPI", "Hoechst 33342", "PI", "LIVE/DEAD Fixable Aqua", "LIVE/DEAD Fixable Blue",
             "LIVE/DEAD Fixable Near-IR", "LIVE/DEAD Fixable Violet", "LIVE/DEAD Fixable Green",
             "LIVE/DEAD Fixable Yellow", "LIVE/DEAD Fixable Red", "LIVE/DEAD Fixable Far Red",
             "FVD eFluor 450", "FVD eFluor 455UV", "FVD eFluor 506", "FVD eFluor 520",
             "FVD eFluor 660", "FVD eFluor 780", "SYTOX Blue", "SYTOX Green", "SYTOX Red",
             "SYTOX AADvanced", "Zombie Aqua", "Zombie Green",
             "Zombie NIR", "Zombie Red", "Zombie UV", "Zombie Violet", "Zombie Yellow"}
# Fluorescent proteins and cell tracking dyes: also not antibody colours.
PROTEINS = {"GFP (EGFP)", "mCherry", "tdTomato", "CFSE", "CellTrace Violet", "CellTrace Far Red",
            "CellTrace Yellow", "CPD eFluor 450", "CPD eFluor 670"}

# Importance level -> lowest brightness worth suggesting
LEVELS = {"high": 4, "medium": 3, "low": 1}


def brightness(dye):
    return BRIGHTNESS.get(dye, 3)


def is_viability_marker(name):
    n = name.lower()
    return any(w in n for w in ("live", "dead", "viab", "zombie", "dapi", "7-aad"))


def suggest_pool(marker, level, fixed, sim, n=12):
    """Dyes worth trying for a marker the user has no antibody for yet.

    fixed: dyes that are certain to be in the panel (markers with a single option).
    High-level markers get bright dyes; low-level markers prefer dimmer ones,
    leaving the bright dyes for the markers that need them.
    """
    dyes = list(sim.index)
    if is_viability_marker(marker):
        pool = [d for d in dyes if d in VIABILITY]
    else:
        pool = [d for d in dyes if d not in VIABILITY and d not in PROTEINS]
    bright = [d for d in pool if brightness(d) >= LEVELS[level]]
    if len(bright) >= 3:
        pool = bright

    def score(d):
        clash = max((sim.loc[d, f] for f in fixed if f != d), default=0.0)
        tilt = {"high": -0.02, "medium": -0.01, "low": 0.01}[level]
        return clash + tilt * brightness(d)

    pool = [d for d in pool if d not in fixed]
    return sorted(pool, key=score)[:n]


def make_detectors():
    """List of (detector name, laser, start nm, end nm)."""
    detectors = []
    for laser, nm in LASERS.items():
        start = nm + 15
        i = 1
        while start + BIN_WIDTH <= LAST_NM:
            detectors.append((f"{laser}{i}", laser, start, start + BIN_WIDTH))
            start += BIN_WIDTH
            i += 1
    return detectors


def signatures_from_spectra(spectra):
    """Simulate each dye's signature across the detectors.

    spectra: long table with columns fluor, kind ('ex' or 'em'), wavelength, value.
    Returns a table: one row per dye, one column per detector, max value = 1.
    """
    detectors = make_detectors()
    grid = np.arange(300, 901)
    rows = {}
    for fluor, d in spectra.groupby("fluor"):
        ex = d[d["kind"] == "ex"].sort_values("wavelength")
        em = d[d["kind"] == "em"].sort_values("wavelength")
        ex_curve = np.interp(grid, ex["wavelength"], ex["value"], left=0, right=0)
        em_curve = np.interp(grid, em["wavelength"], em["value"], left=0, right=0)
        ex_curve = ex_curve / ex_curve.max()
        em_curve = em_curve / em_curve.max()

        sig = []
        for name, laser, a, b in detectors:
            excitation = ex_curve[grid == LASERS[laser]][0]
            emission = em_curve[(grid >= a) & (grid < b)].sum()
            sig.append(excitation * emission)
        sig = np.array(sig)
        rows[fluor] = sig / sig.max()
    return pd.DataFrame(rows, index=[d[0] for d in detectors]).T


def similarity_matrix(signatures):
    """Cosine similarity between every pair of dyes (1 = identical signature)."""
    x = signatures.to_numpy()
    x = x / np.linalg.norm(x, axis=1, keepdims=True)
    return pd.DataFrame(x @ x.T, index=signatures.index, columns=signatures.index)


def complexity_index(signatures, fluors):
    """Condition number of the signature matrix. Lower is easier to unmix."""
    return float(np.linalg.cond(signatures.loc[fluors].to_numpy()))


def worst_pairs(sim, fluors, n=5):
    """The n most similar pairs among the chosen dyes."""
    pairs = [(a, b, sim.loc[a, b]) for a, b in itertools.combinations(fluors, 2)]
    pairs.sort(key=lambda p: -p[2])
    return pd.DataFrame(pairs[:n], columns=["dye 1", "dye 2", "similarity"])


def find_panels(options, key_markers, sim, top=10, max_steps=2_000_000):
    """Search for the best dye for each marker.

    options: {marker: [dyes the user has for it]}
    key_markers: markers that matter most
    Panels are ranked by
      1. the most similar pair that involves a key marker (lower is better)
      2. the most similar pair in the whole panel (lower is better)
    Returns (list of (key_max, all_max, {marker: dye}), stopped_early).
    """
    # key markers first, then markers with fewer choices; this makes the search faster
    order = sorted(options, key=lambda m: (m not in key_markers, len(options[m])))
    # plain arrays: much faster than sim.loc, which matters when running in the browser
    names = list(sim.index)
    pos = {d: i for i, d in enumerate(names)}
    s = sim.to_numpy()
    best = []            # heap of (-key_max, -all_max, counter, panel)
    counter = itertools.count()
    steps = 0
    stopped = False

    def worst_kept():
        if len(best) < top:
            return (np.inf, np.inf)
        k, a, _, _ = best[0]
        return (-k, -a)

    def search(i, chosen, key_max, all_max):
        nonlocal steps, stopped
        if stopped:
            return
        steps += 1
        if steps > max_steps:
            stopped = True
            return
        if (key_max, all_max) >= worst_kept():
            return
        if i == len(order):
            panel = dict(chosen)
            item = (-key_max, -all_max, next(counter), panel)
            if len(best) < top:
                heapq.heappush(best, item)
            else:
                heapq.heapreplace(best, item)
            return
        marker = order[i]
        used = {d for _, d in chosen}
        tries = []
        for dye in options[marker]:
            if dye in used:
                continue
            new_key, new_all = key_max, all_max
            for other_marker, other_dye in chosen:
                v = s[pos[dye], pos[other_dye]]
                new_all = max(new_all, v)
                if marker in key_markers or other_marker in key_markers:
                    new_key = max(new_key, v)
            tries.append((new_key, new_all, dye))
        # most promising dye first, so good panels are found early and the rest is cut off
        for new_key, new_all, dye in sorted(tries, key=lambda t: (t[0], t[1])):
            chosen.append((marker, dye))
            search(i + 1, chosen, new_key, new_all)
            chosen.pop()

    search(0, [], 0.0, 0.0)
    results = sorted([(-k, -a, p) for k, a, _, p in best], key=lambda r: (r[0], r[1]))
    return results, stopped
