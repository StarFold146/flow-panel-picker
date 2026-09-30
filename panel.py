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
        for dye in options[marker]:
            if dye in used:
                continue
            new_key, new_all = key_max, all_max
            for other_marker, other_dye in chosen:
                s = sim.loc[dye, other_dye]
                new_all = max(new_all, s)
                if marker in key_markers or other_marker in key_markers:
                    new_key = max(new_key, s)
            chosen.append((marker, dye))
            search(i + 1, chosen, new_key, new_all)
            chosen.pop()

    search(0, [], 0.0, 0.0)
    results = sorted([(-k, -a, p) for k, a, _, p in best], key=lambda r: (r[0], r[1]))
    return results, stopped
