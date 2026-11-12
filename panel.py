"""Spectral signatures, similarity and panel search for spectral flow cytometry."""
import heapq
import re
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

# Viability dyes that survive fixation (needed for intracellular staining such as FOXP3)
FIXABLE = {d for d in VIABILITY if d.startswith(("Zombie", "LIVE/DEAD Fixable", "FVD"))}

# Dyes most vendors (BioLegend, Thermo, BD) sell on most clones at ordinary prices.
# BD-only dyes (BUV, BB, RealBlue, V450/V500), Thermo's Super Bright (same spectra as BV, rarer)
# and other rarer dyes are left out of suggestions by default.
COMMON = {"FITC", "Alexa Fluor 488", "PE", "PE-Dazzle 594", "PE-CF594", "PE-Texas Red",
          "PE-eFluor 610", "PE-Cy5", "PE-Cy5.5", "PE-Cy7", "PerCP", "PerCP-Cy5.5",
          "PerCP-eFluor 710", "APC", "Alexa Fluor 647", "Alexa Fluor 700", "APC-Cy7",
          "APC-Fire 750", "APC-eFluor 780", "APC-H7", "eFluor 660", "Pacific Blue",
          "eFluor 450", "BV421", "BV510", "BV570", "BV605", "BV650", "BV711",
          "BV785", "DAPI", "7-AAD", "PI", "SYTOX Blue", "SYTOX AADvanced"} | FIXABLE

# The classic colours nearly every lab knows and every vendor sells on most clones:
# no vendor-specific variants (eFluor, Fire, PE-Texas Red, PerCP-eFluor 710 ...)
CLASSIC = {"FITC", "Alexa Fluor 488", "PE", "PE-Dazzle 594", "PE-CF594", "PE-Cy5", "PE-Cy7",
           "PerCP-Cy5.5", "APC", "Alexa Fluor 647", "Alexa Fluor 700", "APC-Cy7", "Pacific Blue",
           "BV421", "BV510", "BV605", "BV650", "BV711", "BV785",
           "DAPI", "7-AAD", "PI", "Zombie NIR", "Zombie Aqua", "Zombie UV", "Zombie Violet",
           "LIVE/DEAD Fixable Blue", "LIVE/DEAD Fixable Aqua", "LIVE/DEAD Fixable Near-IR",
           "LIVE/DEAD Fixable Violet", "FVD eFluor 780", "FVD eFluor 506"}

# Which dyes suggestions may use
TIERS = {"classic": CLASSIC, "common": COMMON, "all": None}

# Dyes that are the same colour under different vendors' names (signatures within ~0.90).
# The built-in list shows each group as one dye, so you can buy whichever vendor has your clone.
GROUPS = {
    "PE-CF594 / PE-Dazzle 594 / PE-eFluor 610": ["PE-CF594", "PE-Dazzle 594", "PE-eFluor 610", "PE-Texas Red"],
    "APC-Cy7 / APC-H7 / APC-Fire 750 / APC-eFluor 780": ["APC-Cy7", "APC-H7", "APC-Fire 750", "APC-eFluor 780"],
    "PerCP-Cy5.5 / PerCP-eFluor 710": ["PerCP-Cy5.5", "PerCP-eFluor 710"],
    "FITC / Alexa Fluor 488": ["FITC", "Alexa Fluor 488"],
    "APC / Alexa Fluor 647 / eFluor 660": ["APC", "Alexa Fluor 647", "eFluor 660"],
    "Pacific Blue / eFluor 450 / V450": ["Pacific Blue", "eFluor 450", "V450"],
    "BV421 / Super Bright 436": ["BV421", "Super Bright 436"],
    "BV510 / V500": ["BV510", "V500"],
    "BV605 / Super Bright 600": ["BV605", "Super Bright 600"],
    "BV650 / Super Bright 645": ["BV650", "Super Bright 645"],
    "BV785 / Super Bright 780": ["BV785", "Super Bright 780"],
    "Zombie Aqua / LIVE/DEAD Fixable Aqua": ["Zombie Aqua", "LIVE/DEAD Fixable Aqua"],
    "Zombie UV / LIVE/DEAD Fixable Blue / FVD eFluor 455UV": ["Zombie UV", "LIVE/DEAD Fixable Blue", "FVD eFluor 455UV"],
    "Zombie Green / LIVE/DEAD Fixable Green / FVD eFluor 520": ["Zombie Green", "LIVE/DEAD Fixable Green", "FVD eFluor 520"],
    "Zombie Yellow / LIVE/DEAD Fixable Yellow": ["Zombie Yellow", "LIVE/DEAD Fixable Yellow"],
    "Zombie Red / LIVE/DEAD Fixable Red": ["Zombie Red", "LIVE/DEAD Fixable Red"],
}
GROUP_OF = {m: g for g, ms in GROUPS.items() for m in ms}

# a group belongs to a set (viability, a tier ...) when one of its dyes does
for _set in (VIABILITY, PROTEINS, FIXABLE, CLASSIC, COMMON):
    _set.update(g for g, ms in GROUPS.items() if _set & set(ms))
for g, ms in GROUPS.items():
    BRIGHTNESS[g] = max(BRIGHTNESS.get(m, 3) for m in ms)


def merge_groups(signatures):
    """One row per group of same-colour dyes: the mean of its members' signatures."""
    sig = signatures.copy()
    for g, ms in GROUPS.items():
        present = [m for m in ms if m in sig.index]
        if not present:
            continue
        row = sig.loc[present].mean()
        sig = sig.drop(present)
        sig.loc[g] = row / row.max()
    return sig


def resolve(dye, names):
    """The name a dye goes by in this list (its group's name if it was merged)."""
    if dye in names:
        return dye
    g = GROUP_OF.get(dye)
    return g if g in names else None


# Importance level -> lowest brightness worth suggesting
LEVELS = {"high": 4, "medium": 3, "low": 1}


def brightness(dye):
    return BRIGHTNESS.get(dye, 3)


def is_viability_marker(name):
    n = name.lower()
    return any(w in n for w in ("live", "dead", "viab", "zombie", "dapi", "7-aad"))


# Markers usually stained after fixing and permeabilising (cytoplasm or nucleus)
INSIDE = ("foxp3", "helios", "ikzf", "ki67", "ki-67", "t-bet", "tbet", "gata3", "gata-3", "ror",
          "eomes", "bcl6", "bcl-6", "blimp", "tcf1", "tcf-1", "tox", "irf4", "batf", "nfat",
          "ifn", "il-", "il1", "il2", "il4", "il5", "il6", "il9", "il10", "il13", "il17", "il21",
          "il22", "tnf", "tgf", "gzm", "granzyme", "perforin", "prf", "ctla", "pstat", "p-stat",
          "phospho", "caspase", "h2ax")


def needs_perm(name):
    """Guess from the name whether a marker is stained inside the cell."""
    n = name.lower().replace(" ", "").replace("_", "-")
    n = re.sub(r"il-?\d+r", "", n)   # IL-7R, IL2RA and other receptors sit on the surface
    return any(w in n for w in INSIDE)


def suggest_pool(marker, level, fixed, sim, n=20, tier="classic", fixable=True):
    """Dyes worth trying for a marker the user has no antibody for yet.

    fixed: dyes that are certain to be in the panel (markers with a single option).
    fixable: the cells will be fixed and permeabilised, so live/dead dyes must survive that;
    otherwise live/dead markers get the cheaper DNA dyes that only work on unfixed cells.
    High-level markers get bright dyes; low-level markers prefer dimmer ones,
    leaving the bright dyes for the markers that need them.
    """
    dyes = list(sim.index)
    allowed = TIERS[tier]
    if allowed is not None and len([d for d in dyes if d in allowed]) >= 10:
        dyes = [d for d in dyes if d in allowed]
    if is_viability_marker(marker):
        pool = [d for d in dyes if d in VIABILITY]
        # fixed cells need a fixable dye; unfixed cells can use the cheap DNA dyes (DAPI, 7-AAD, ...)
        pool = [d for d in pool if (d in FIXABLE) == fixable] or pool
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


# Similarity bands: below CAREFUL is fine, CAREFUL to TOO_CLOSE needs care, above is very hard to unmix
CAREFUL, TOO_CLOSE = 0.90, 0.98


def band(v):
    return 0 if v < CAREFUL else 1 if v < TOO_CLOSE else 2


def rank(key_max, all_max):
    """Sort key for a panel, lower is better.

    1. the band of the most similar pair in the whole panel: no pair should be too close
    2. the band of the most similar pair that involves a key marker
    3. a balance of the two, in steps of 0.05: key pairs count fully, other pairs are
       given 0.15 of slack, so cleaner key markers win unless the rest of the panel
       gets clearly worse for it
    4. the most similar pair in the whole panel, then the key pairs
    Every part only grows as dyes are added, so the search can cut branches with it.
    """
    balance = max(key_max, all_max - 0.15)
    return (band(all_max), band(key_max), int(np.ceil(balance / 0.05 - 1e-9)), all_max, key_max)


def find_panels(options, key_markers, sim, top=10, max_steps=2_000_000):
    """Search for the best dye for each marker.

    options: {marker: [dyes the user has for it]}
    key_markers: markers that matter most
    Panels are ranked by rank(): no pair too close first, then clean key markers,
    then the most similar pair in the whole panel.
    Returns (list of (key_max, all_max, {marker: dye}), stopped_early).
    """
    # key markers first, then markers with fewer choices; this makes the search faster
    order = sorted(options, key=lambda m: (m not in key_markers, len(options[m])))
    # plain arrays: much faster than sim.loc, which matters when running in the browser
    names = list(sim.index)
    pos = {d: i for i, d in enumerate(names)}
    s = sim.to_numpy()
    best = []            # heap of (negated rank, counter, key_max, all_max, panel)
    counter = itertools.count()
    steps = 0
    stopped = False

    def worst_kept():
        if len(best) < top:
            return (np.inf,)
        return tuple(-x for x in best[0][0])

    def search(i, chosen, key_max, all_max):
        nonlocal steps, stopped
        if stopped:
            return
        steps += 1
        if steps > max_steps:
            stopped = True
            return
        r = rank(key_max, all_max)
        if r >= worst_kept():
            return
        if i == len(order):
            panel = dict(chosen)
            item = (tuple(-x for x in r), next(counter), key_max, all_max, panel)
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
        for new_key, new_all, dye in sorted(tries, key=lambda t: rank(t[0], t[1])):
            chosen.append((marker, dye))
            search(i + 1, chosen, new_key, new_all)
            chosen.pop()

    search(0, [], 0.0, 0.0)
    results = sorted([(k, a, p) for _, _, k, a, p in best], key=lambda r: rank(r[0], r[1]))
    return results, stopped


def panel_score(panel, key_markers, s, pos):
    """(key_max, all_max) of a {marker: dye} panel, using the plain similarity array."""
    items = list(panel.items())
    key_max = all_max = 0.0
    for (m1, d1), (m2, d2) in itertools.combinations(items, 2):
        v = s[pos[d1], pos[d2]]
        all_max = max(all_max, v)
        if m1 in key_markers or m2 in key_markers:
            key_max = max(key_max, v)
    return key_max, all_max


def local_panels(options, key_markers, sim, restarts=40, seed=0):
    """For many markers with many choices, where trying every combination is too slow.

    Each round builds a panel greedily (key markers first, in a shuffled order, picking
    among the best few dyes at random), then keeps changing one marker's dye, or swapping
    the dyes of two markers, while that makes the panel better. Different rounds end in
    different panels. Returns the same shape as find_panels.
    """
    rng = np.random.default_rng(seed)
    names = list(sim.index)
    pos = {d: i for i, d in enumerate(names)}
    s = sim.to_numpy()
    found = {}
    for _ in range(restarts):
        markers = list(options)
        rng.shuffle(markers)
        markers.sort(key=lambda m: (m not in key_markers, len(options[m])))
        panel = {}
        for m in markers:
            used = set(panel.values())
            tries = []
            for d in options[m]:
                if d in used:
                    continue
                k, a = panel_score({**panel, m: d}, key_markers, s, pos)
                tries.append((rank(k, a), d))
            if not tries:
                break
            tries.sort()
            panel[m] = tries[min(int(rng.integers(0, 3)), len(tries) - 1)][1]
        if len(panel) < len(options):
            continue
        best = rank(*panel_score(panel, key_markers, s, pos))
        better = True
        while better:
            better = False
            for m in options:
                used = set(panel.values())
                for d in options[m]:
                    if d in used:
                        continue
                    r = rank(*panel_score({**panel, m: d}, key_markers, s, pos))
                    if r < best:
                        panel, best, better = {**panel, m: d}, r, True
                        used = set(panel.values())
            for m1, m2 in itertools.combinations(options, 2):
                d1, d2 = panel[m1], panel[m2]
                if d2 in options[m1] and d1 in options[m2]:
                    r = rank(*panel_score({**panel, m1: d2, m2: d1}, key_markers, s, pos))
                    if r < best:
                        panel, best, better = {**panel, m1: d2, m2: d1}, r, True
        found[tuple(sorted(panel.items()))] = panel
    results = [(*panel_score(p, key_markers, s, pos), p) for p in found.values()]
    results.sort(key=lambda r: rank(r[0], r[1]))
    return results


def diverse(results, n=10, min_diff=2):
    """Keep up to n panels that each differ from the ones already kept in at least min_diff markers."""
    kept = []
    for r in results:
        if all(sum(r[2][m] != k[2][m] for m in r[2]) >= min_diff for k in kept):
            kept.append(r)
        if len(kept) == n:
            break
    return kept


def alternatives(panel, options, sim, n=3):
    """For each marker, other dyes from its options that would fit the rest of the panel.

    Returns {marker: [(dye, most similar to the rest), ...]}, best first, only dyes that
    stay below the 'be careful' band.
    """
    out = {}
    for m, dye in panel.items():
        rest = [d for k, d in panel.items() if k != m]
        alts = []
        for d in options.get(m, []):
            if d == dye or d in rest:
                continue
            worst = max(sim.loc[d, rest]) if rest else 0.0
            if worst < CAREFUL:
                alts.append((d, float(worst)))
        out[m] = sorted(alts, key=lambda a: a[1])[:n]
    return out
