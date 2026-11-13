import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from export import workbook
from panel import (apart_pairs, panel_score, parse_cells, CO_EXPRESSED, EXPRESSION, FRAGILE, staining_notes, merge_groups, resolve, CAREFUL, LEVELS, TOO_CLOSE, alternatives,
                   brightness, diverse, local_panels, misfit, needs_perm, rank, complexity_index, find_panels,
                   similarity_matrix, suggest_pool, typical_expression, worst_pairs)

st.set_page_config(page_title="Spectral Panel Picker", layout="wide")


@st.cache_data
def load_builtin():
    # precomputed from data/spectra.csv by fetch_spectra.py, so the page loads faster
    return merge_groups(pd.read_csv("data/signatures.csv", index_col=0))


@st.cache_data
def load_uploaded(file):
    sig = pd.read_csv(file, index_col=0)
    sig = sig.div(sig.max(axis=1), axis=0)
    return sig


def color_similarity(v):
    if v >= TOO_CLOSE:
        return "background-color: #f4a6a6"
    if v >= CAREFUL:
        return "background-color: #fde2a7"
    return ""


def warn_pairs(fluors, labels, sim, markers=None, apart=()):
    """Say out loud which pairs are too close, instead of leaving it to the heatmap.

    markers / apart: the marker of each dye and the pairs of markers never on the same cell;
    close colours on such pairs are fine below 0.98.
    """
    name = dict(zip(fluors, labels))
    marker = dict(zip(fluors, markers or fluors))
    pairs = worst_pairs(sim, fluors, n=len(fluors) ** 2)
    separate = pairs.apply(lambda r: frozenset((marker[r["dye 1"]], marker[r["dye 2"]])) in apart, axis=1) \
        if len(pairs) else pairs["similarity"] > 2
    bad = pairs[pairs["similarity"] >= TOO_CLOSE]
    ok_apart = pairs[separate & (pairs["similarity"] >= CO_EXPRESSED) & (pairs["similarity"] < TOO_CLOSE)]
    pairs = pairs[~separate | (pairs["similarity"] >= TOO_CLOSE)]
    care = pairs[(pairs["similarity"] >= CAREFUL) & (pairs["similarity"] < TOO_CLOSE)]
    if len(bad):
        st.error("Too similar to unmix well, change one of each pair: " + "; ".join(
            f"{name[a]} and {name[b]} ({v:.2f})" for a, b, v in bad.itertuples(index=False)))
    if len(care):
        st.warning("Close, OK only if the two markers are not co-expressed and neither is dim: " + "; ".join(
            f"{name[a]} and {name[b]} ({v:.2f})" for a, b, v in care.itertuples(index=False)))
    if not len(bad) and not len(care):
        st.success(f"No pair above {CAREFUL:.2f} similarity" +
                   (" among markers on the same cells." if len(ok_apart) else "."))
    if len(ok_apart):
        st.caption("Close colours on markers that never share a cell, fine: " + "; ".join(
            f"{name[a]} and {name[b]} ({v:.2f})" for a, b, v in ok_apart.itertuples(index=False)))
    near = pairs[(pairs["similarity"] >= CO_EXPRESSED) & (pairs["similarity"] < CAREFUL)]
    if len(near):
        st.info(f"Fine in general, but keep apart if they are on the same cells (OMIP-105 kept "
                f"co-expressed markers below {CO_EXPRESSED:.2f}): " + "; ".join(
                    f"{name[a]} and {name[b]} ({v:.2f})" for a, b, v in near.itertuples(index=False)))


def show_panel(fluors, labels, sig, sim, markers=None, apart=()):
    """Heatmap, worst pairs and signature plot for one set of dyes."""
    warn_pairs(fluors, labels, sim, markers, apart)
    sub = sim.loc[fluors, fluors].copy()
    sub.index = labels
    sub.columns = labels

    c1, c2 = st.columns([3, 2])
    with c1:
        fig = px.imshow(sub.round(2), text_auto=True, zmin=0, zmax=1,
                        color_continuous_scale="Reds", aspect="auto")
        fig.update_layout(height=150 + 32 * len(fluors), margin=dict(t=20))
        st.plotly_chart(fig, width="stretch")
    with c2:
        st.metric("Complexity index", f"{complexity_index(sig, fluors):.1f}",
                  help="Condition number of the signature matrix. Lower is easier to unmix.")
        pairs = worst_pairs(sim, fluors)
        name = dict(zip(fluors, labels))
        pairs["dye 1"] = pairs["dye 1"].map(name)
        pairs["dye 2"] = pairs["dye 2"].map(name)
        st.write("Most similar pairs")
        st.dataframe(pairs.style.format({"similarity": "{:.2f}"})
                     .map(color_similarity, subset=["similarity"]), hide_index=True)

    long = sig.loc[fluors].copy()
    long.index = labels
    long = long.reset_index(names="dye").melt(id_vars="dye", var_name="detector", value_name="signal")
    fig = px.line(long, x="detector", y="signal", color="dye", markers=True)
    fig.update_layout(height=350, margin=dict(t=20), xaxis_tickangle=-90)
    st.plotly_chart(fig, width="stretch")


# ---------------- sidebar ----------------
st.sidebar.header("Spectra")
upload = st.sidebar.file_uploader(
    "Optional: your own reference signatures (.csv)",
    help="One row per dye, one column per detector, first column = dye name. "
         "For example exported from single-stain reference controls.")
if upload is not None:
    sig = load_uploaded(upload)
    st.sidebar.success(f"Using your signatures: {len(sig)} dyes")
else:
    sig = load_builtin()
    st.sidebar.info("Using built-in signatures simulated from public spectra (FPbase) "
                    "for BD FACSDiscover S8 lasers. Good for ranking, the exact numbers "
                    "will differ from the instrument.")
sim = similarity_matrix(sig)
all_dyes = sorted(sig.index)
tier = st.sidebar.radio(
    "Suggest dyes from", ["classic", "common", "all"],
    format_func={"classic": "Classic colours (BV, PE, APC, FITC ...)",
                 "common": "+ vendor variants (eFluor, Fire, PerCP-eFluor 710 ...)",
                 "all": "Every dye (incl. BD-only BUV/BB)"}.get,
    help="Only changes what is suggested for markers with no dyes; dyes you pick yourself are always used.")

st.sidebar.markdown("""
**Colour guide**
- similarity ≥ 0.98: very hard to unmix
- 0.90 – 0.98: be careful, avoid for dim or co-expressed markers
- 0.70 – 0.90: fine, unless the two markers are on the same cells

**Can I buy it?**
This tool doesn't know which clones are sold in which colour. Check there:
- [FluoroFinder Panel Builder](https://app.fluorofinder.com/panels/new) (all vendors, free academic account)
- [Biocompare Panel Builder](https://app.fluorofinder.com/biocompare/x) (same, no account)
- [Thermo Fisher Panel Builder](https://www.thermofisher.com/order/panel-builder)
- [BioLegend Multicolor Panel Selector](https://www.biolegend.com/multicolor_staining)
- [BD Panel Designer](https://www.bdbiosciences.com/en-us/resources/panel-design)

Dyes a vendor doesn't have for your clone: tick them under **To buy** and it picks again.
""")

# ---------------- main ----------------
st.title("Spectral Panel Picker")
st.write("Tell it which dyes you have for each marker. It picks one dye per marker so that "
         "the most similar pair in the panel is as different as possible. "
         "Set markers that matter most to **high** importance: they get the cleanest dyes first. "
         "Set how strongly each marker is **expressed**: dim markers get bright dyes, bright markers "
         "leave the bright dyes to others. Intracellular markers avoid PE/APC tandems, which fixation can break. "
         "Leave a marker's dyes empty and it will suggest a colour to buy.")

tab_pick, tab_check = st.tabs(["Pick a panel", "Check a panel"])

with tab_pick:
    ss = st.session_state
    listing = st.file_uploader(
        "Optional: load a marker list (.csv)", type="csv", key="listing",
        help="Columns: marker, dyes, importance, expression. Dyes separated by ';' (leave empty to get a "
             "suggestion), importance high / medium / low, expression dim / medium / bright "
             "(empty = guessed from the name). 'Download this list' below saves one.")
    if listing is not None and ss.get("listing_id") != listing.file_id:
        table0 = pd.read_csv(listing, dtype=str).fillna("")
        table0.columns = [c.strip().lower() for c in table0.columns]
        unknown = []
        for i, r in enumerate(table0.to_dict("records")):
            dyes = []
            for d in r.get("dyes", "").split(";"):
                if d.strip():
                    found = resolve(d.strip(), all_dyes)
                    if found and found not in dyes:
                        dyes.append(found)
                    elif not found:
                        unknown.append(d.strip())
            ss[f"name{i}"] = r.get("marker", "").strip()
            ss[f"dyes{i}"] = dyes
            level = r.get("importance", "").strip().lower()
            ss[f"level{i}"] = level if level in LEVELS else "medium"
            expr = r.get("expression", "").strip().lower()
            ss[f"expr{i}"] = expr if expr in EXPRESSION else typical_expression(ss[f"name{i}"])
            ss[f"cells{i}"] = r.get("cells", "").strip()
        if "cells" in table0.columns and table0["cells"].str.strip().any():
            ss["use_cells"] = True
        ss["n"] = max(2, len(table0))
        ss["listing_id"] = listing.file_id
        ss["listing_note"] = (f"Loaded {len(table0)} markers." +
                              (f" Not in the dye list, skipped: {', '.join(unknown)}" if unknown else ""))
    if ss.get("listing_note"):
        st.caption(ss["listing_note"])

    # starting values go into session state before the widgets are drawn, and the widgets
    # read only from there; giving a widget a value as well would make it a new widget
    # (and reset it) whenever that value stops being passed
    def first(key, value):
        ss.setdefault(key, value)
        return {}

    n = st.number_input("Number of markers", min_value=2, max_value=40, key="n", **first("n", 6))

    example = {
        0: ("CD4", ["BUV737", "BV605", "PE-Cy7"], False),
        1: ("CD8", ["BV711", "PerCP-Cy5.5", "BUV805"], False),
        2: ("CD25", ["BV421", "PE", "APC", "BB515"], True),
        3: ("FOXP3", ["PE", "APC", "Alexa Fluor 488", "eFluor 450"], True),
        4: ("Helios", ["PE-Cy7", "Alexa Fluor 647", "FITC"], False),
        5: ("Live/Dead", ["Zombie NIR", "Zombie Aqua", "LIVE/DEAD Fixable Blue"], False),
    }

    use_cells = st.toggle(
        "My markers are on different cell types", key="use_cells", **first("use_cells", False),
        help="Adds a 'Cells' column. Two markers with cell types in common (or one left empty) count as "
             "on the same cells. Markers that never share a cell can have closer colours: spill from one "
             "lands on cells without the other. Leave it off for a panel of one cell type.")
    widths = [1.2, 3.4, 0.9, 0.9, 1.2] if use_cells else [1.2, 4, 0.9, 0.9]
    heads = st.columns(widths)
    for h, text in zip(heads, ["Marker", "Dyes you have", "Importance", "Expression", "Cells"]):
        h.caption(text)
    rows = []
    for i in range(int(n)):
        name0, dyes0, key0 = example.get(i, (f"Marker {i + 1}", [], False))
        c1, c2, c3, c4, *c5 = st.columns(widths)
        name = c1.text_input("Marker", key=f"name{i}", label_visibility="collapsed",
                             **first(f"name{i}", name0))
        dyes = c2.multiselect("Dyes", all_dyes, key=f"dyes{i}", label_visibility="collapsed",
                              **first(f"dyes{i}", [r for r in (resolve(d, all_dyes) for d in dyes0) if r]),
                              placeholder="dyes you have (empty = suggest one to buy)")
        level = c3.selectbox("Importance", list(LEVELS), key=f"level{i}", label_visibility="collapsed",
                             **first(f"level{i}", "high" if key0 else "medium"),
                             help="high: key marker, its pairs must be cleanest. low: can live with some spread.")
        expr = c4.selectbox("Expression", EXPRESSION, key=f"expr{i}", label_visibility="collapsed",
                            **first(f"expr{i}", typical_expression(name0)),
                            help="dim: needs a bright dye. bright: fine on a dim dye. "
                                 "First guess comes from the name; change it if you know better.")
        cells = ""
        if use_cells:
            cells = c5[0].text_input("Cells", key=f"cells{i}", label_visibility="collapsed",
                                     **first(f"cells{i}", ""), placeholder="e.g. T, Treg (empty = all)")
        if name:
            rows.append((name, dyes, level, expr, cells))

    seen, dupes = set(), set()
    for r in rows:
        (dupes if r[0] in seen else seen).add(r[0])
    if dupes:
        st.warning("Same marker name twice, only the first is used: " + ", ".join(sorted(dupes)))
        rows = [r for k, r in enumerate(rows) if r[0] not in {x[0] for x in rows[:k]}]

    st.download_button("Download this list", pd.DataFrame(
        [{"marker": name, "dyes": "; ".join(dyes), "importance": level, "expression": expr,
          **({"cells": cells} if use_cells else {})}
         for name, dyes, level, expr, cells in rows]
    ).to_csv(index=False).encode(), "markers.csv", key="save_list")

    key_markers = {name for name, _, level, *_ in rows if level == "high"}
    options = {name: dyes for name, dyes, *_ in rows if dyes}
    apart = apart_pairs({name: parse_cells(cells) for name, *_, cells in rows}) if use_cells else set()
    fixed = {dyes[0] for dyes in options.values() if len(dyes) == 1}
    inside = [name for name, *_ in rows if needs_perm(name)]
    fixable = st.checkbox(
        "Fix and permeabilise (intracellular or nuclear staining)", value=bool(inside),
        key=f"perm{bool(inside)}",
        help="On: live/dead suggestions are fixable dyes only (Zombie, LIVE/DEAD Fixable, FVD). "
             "Off: the cheaper DAPI, 7-AAD, PI or SYTOX instead.")
    if inside:
        st.caption("Ticked because these look intracellular: " + ", ".join(inside))
    banned = ss.setdefault("banned", set())     # (marker, dye) the vendor doesn't sell
    info = {name: (expr, fixable and needs_perm(name)) for name, _, _, expr, _ in rows}
    suggested = set()
    for name, dyes, level, expr, _ in rows:
        if not dyes:
            options[name] = suggest_pool(name, level, fixed, sim, tier=tier, fixable=fixable,
                                         expression=expr, banned={d for m, d in banned if m == name})
            suggested.add(name)
    if banned:
        b1, b2 = st.columns([5, 1])
        b1.caption("Not sold, left out of suggestions: " + "; ".join(f"{m} {d}" for m, d in sorted(banned)))
        if b2.button("Clear", key="clear_banned"):
            ss["banned"] = set()
            st.rerun()

    # what the results were computed from, to notice when the inputs change afterwards
    inputs = repr((rows, tier, fixable, sorted(banned), use_cells))

    if st.button("Find panels", type="primary") or ss.pop("search_again", False):
        if len(options) < 2:
            st.warning("Add at least two markers.")
        else:
            combos = float(np.prod([len(v) for v in options.values()], dtype=float))
            with st.spinner("Searching..."):
                if combos <= 5e6:
                    results, stopped = find_panels(options, key_markers, sim, top=200, info=info, apart=apart)
                    if stopped:
                        results += local_panels(options, key_markers, sim, info=info, apart=apart)
                else:
                    # far too many to try one by one: build and improve panels instead
                    results = local_panels(options, key_markers, sim, info=info, apart=apart)
                results.sort(key=lambda r: (rank(*r[:3]), complexity_index(sig, list(r[-1].values()))))
                # panels that differ only in one or two low-priority markers are not real options
                results = diverse(results, n=10, min_diff=max(1, round(len(options) / 4)))
            if not results:
                st.error("No panel found. Two markers may only have the same dye.")
            else:
                table = []
                # results are already in rank order (with the cell-type relief); the table shows
                # the real similarities
                s_arr, pos = sim.to_numpy(), {d: i for i, d in enumerate(sim.index)}
                for order_no, (_, _, mis, panel) in enumerate(results):
                    key_max, all_max, _ = panel_score(panel, key_markers, s_arr, pos)
                    fl = list(panel.values())
                    row = {"#": order_no, "key pairs max": key_max, "all pairs max": all_max, "dye fit": mis,
                           "complexity": complexity_index(sig, fl)}
                    row.update(panel)
                    table.append(row)
                table = pd.DataFrame(table)
                table = table.sort_values("#").drop(columns="#").head(10).reset_index(drop=True)
                table.index = table.index + 1
                st.session_state["results"] = table
                st.session_state["suggested"] = suggested
                st.session_state["levels"] = {name: level for name, _, level, *_ in rows}
                st.session_state["apart"] = apart
                st.session_state["info"] = info
                st.session_state["key_markers"] = key_markers
                st.session_state["options"] = options
                st.session_state["inputs"] = inputs

    if "results" in st.session_state:
        table = st.session_state["results"]
        markers = [c for c in table.columns if c not in ["key pairs max", "all pairs max", "dye fit", "complexity"]]
        st.subheader("Best panels")
        if ss.get("inputs") != inputs:
            st.info("The markers or settings changed since this search. Press **Find panels** to update.")
        st.dataframe(table.style.format({"key pairs max": "{:.2f}", "all pairs max": "{:.2f}",
                                         "complexity": "{:.1f}"}), width="stretch")
        st.caption("dye fit: 0 means every dim marker has a bright dye and no intracellular marker sits on "
                   "a fragile tandem; higher means some compromise (see the notes under the panel).")
        info_r = ss.get("info", {})

        def note(m, dye):
            expression, inside_m = info_r.get(m, ("medium", False))
            out = []
            if inside_m and dye in FRAGILE:
                out.append("tandem on an intracellular marker, may degrade after fix/perm")
            if expression == "dim" and brightness(dye) < 4:
                out.append("dim marker on a dimmer dye")
            if expression == "bright" and brightness(dye) >= 5:
                out.append("bright marker on one of the brightest dyes")
            return "; ".join(out)

        d1, d2 = st.columns(2)
        d1.download_button("Download as Excel (one sheet per panel)", workbook(
            table, markers, sig, sim, key_markers=ss.get("key_markers", set()),
            to_buy=ss.get("suggested", set()),
            note=note, bench=lambda panel: staining_notes(panel, fixable)),
            "panels.xlsx", type="primary")
        d2.download_button("Download table (.csv)", table.to_csv().encode(), "panels.csv")

        pick = st.selectbox("Look at panel", table.index)
        row = table.loc[pick]
        fluors = [row[m] for m in markers]
        labels = [f"{m} ({row[m]})" for m in markers]
        to_buy = [m for m in markers if m in st.session_state.get("suggested", set())]
        alts = alternatives({m: row[m] for m in markers}, st.session_state.get("options", {}), sim)
        issues = {m: note(m, row[m]) for m in markers if note(m, row[m])}
        if issues:
            st.warning("Worth a second look: " + "; ".join(f"{m} ({row[m]}): {v}" for m, v in issues.items()))
        if to_buy:
            st.write("**To buy** for this panel (suggested colours)")
            levels = st.session_state.get("levels", {})
            buy = pd.DataFrame([{"marker": m, "dye": row[m], "importance": levels.get(m, ""),
                                 "brightness (1-5)": brightness(row[m]),
                                 "also fits": ", ".join(d for d, _ in alts.get(m, []))} for m in to_buy])
            st.dataframe(buy, hide_index=True)
            st.download_button("Download shopping list", buy.to_csv(index=False).encode(),
                               "to_buy.csv", key="buy")
            st.caption("Check the vendor actually sells your clone in that colour (links in the sidebar). "
                       "'Also fits' keeps every pair below 0.90 with the rest of this panel. "
                       "Brightness is a rough guide from vendor charts.")
            not_sold = st.multiselect("Not sold for my clone", [f"{m}: {row[m]}" for m in to_buy],
                                      key=f"not_sold{pick}",
                                      placeholder="tick colours you can't buy, then search again")
            if not_sold and st.button("Search again without these"):
                for item in not_sold:
                    m, d = item.split(": ", 1)
                    ss["banned"].add((m, d))
                ss["search_again"] = True
                st.rerun()
        show_panel(fluors, labels, sig, sim, markers, ss.get("apart", set()))
        if ss.get("apart"):
            st.caption("Ranked with markers on different cells allowed closer colours "
                       "(their similarity counts 0.10 less; pairs at 0.98 or above still count in full).")
        with st.expander("Bench notes for this panel (from OMIP-104, OMIP-105 and related papers)"):
            for line in staining_notes({m: row[m] for m in markers}, fixable):
                st.markdown("- " + line)

with tab_check:
    st.write("Already have a panel? Pick its dyes to see which pairs are too close.")
    chosen = st.multiselect("Dyes in your panel", all_dyes,
                            default=[resolve(d, all_dyes) for d in ["DAPI", "GFP (EGFP)", "PE", "APC"]])
    if len(chosen) >= 2:
        show_panel(chosen, chosen, sig, sim)
