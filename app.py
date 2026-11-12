import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from panel import (CAREFUL, LEVELS, TOO_CLOSE, alternatives, brightness, diverse, local_panels, needs_perm, rank, complexity_index, find_panels,
                   signatures_from_spectra, similarity_matrix, suggest_pool, worst_pairs)

st.set_page_config(page_title="Spectral Panel Picker", layout="wide")


@st.cache_data
def load_builtin():
    spectra = pd.read_csv("data/spectra.csv")
    return signatures_from_spectra(spectra)


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


def warn_pairs(fluors, labels, sim):
    """Say out loud which pairs are too close, instead of leaving it to the heatmap."""
    name = dict(zip(fluors, labels))
    pairs = worst_pairs(sim, fluors, n=len(fluors) ** 2)
    bad = pairs[pairs["similarity"] >= TOO_CLOSE]
    care = pairs[(pairs["similarity"] >= CAREFUL) & (pairs["similarity"] < TOO_CLOSE)]
    if len(bad):
        st.error("Too similar to unmix well, change one of each pair: " + "; ".join(
            f"{name[a]} and {name[b]} ({v:.2f})" for a, b, v in bad.itertuples(index=False)))
    if len(care):
        st.warning("Close, OK only if the two markers are not co-expressed and neither is dim: " + "; ".join(
            f"{name[a]} and {name[b]} ({v:.2f})" for a, b, v in care.itertuples(index=False)))
    if not len(bad) and not len(care):
        st.success(f"No pair above {CAREFUL:.2f} similarity.")


def show_panel(fluors, labels, sig, sim):
    """Heatmap, worst pairs and signature plot for one set of dyes."""
    warn_pairs(fluors, labels, sim)
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
""")

# ---------------- main ----------------
st.title("Spectral Panel Picker")
st.write("Tell it which dyes you have for each marker. It picks one dye per marker so that "
         "the most similar pair in the panel is as different as possible. "
         "Set markers that matter most (or are dim) to **high**: they get the cleanest and brightest dyes first. "
         "Leave a marker's dyes empty and it will suggest a colour to buy.")

tab_pick, tab_check = st.tabs(["Pick a panel", "Check a panel"])

with tab_pick:
    n = st.number_input("Number of markers", min_value=2, max_value=40, value=6)

    example = {
        0: ("CD4", ["BUV737", "BV605", "PE-Cy7"], False),
        1: ("CD8", ["BV711", "PerCP-Cy5.5", "BUV805"], False),
        2: ("CD25", ["BV421", "PE", "APC", "BB515"], True),
        3: ("FOXP3", ["PE", "APC", "Alexa Fluor 488", "eFluor 450"], True),
        4: ("Helios", ["PE-Cy7", "Alexa Fluor 647", "FITC"], False),
        5: ("Live/Dead", ["Zombie NIR", "Zombie Aqua", "LIVE/DEAD Fixable Blue"], False),
    }

    rows = []
    for i in range(int(n)):
        name0, dyes0, key0 = example.get(i, (f"Marker {i + 1}", [], False))
        c1, c2, c3 = st.columns([1.2, 4, 0.9])
        name = c1.text_input("Marker", value=name0, key=f"name{i}", label_visibility="collapsed")
        dyes = c2.multiselect("Dyes", all_dyes, default=[d for d in dyes0 if d in all_dyes],
                              key=f"dyes{i}", label_visibility="collapsed",
                              placeholder="dyes you have (empty = suggest one to buy)")
        level = c3.selectbox("Importance", list(LEVELS), index=0 if key0 else 1, key=f"level{i}",
                             label_visibility="collapsed",
                             help="high: key or dim marker, gets bright clean dyes. "
                                  "low: bright marker, happy with a dim dye.")
        if name:
            rows.append((name, dyes, level))

    key_markers = {name for name, _, level in rows if level == "high"}
    options = {name: dyes for name, dyes, _ in rows if dyes}
    fixed = {dyes[0] for dyes in options.values() if len(dyes) == 1}
    inside = [name for name, _, _ in rows if needs_perm(name)]
    fixable = st.checkbox(
        "Fix and permeabilise (intracellular or nuclear staining)", value=bool(inside),
        key=f"perm{bool(inside)}",
        help="On: live/dead suggestions are fixable dyes only (Zombie, LIVE/DEAD Fixable, FVD). "
             "Off: the cheaper DAPI, 7-AAD, PI or SYTOX instead.")
    if inside:
        st.caption("Ticked because these look intracellular: " + ", ".join(inside))
    suggested = set()
    for name, dyes, level in rows:
        if not dyes:
            options[name] = suggest_pool(name, level, fixed, sim, tier=tier,
                                         fixable=fixable)
            suggested.add(name)

    if st.button("Find panels", type="primary"):
        if len(options) < 2:
            st.warning("Add at least two markers.")
        else:
            combos = float(np.prod([len(v) for v in options.values()], dtype=float))
            with st.spinner("Searching..."):
                if combos <= 5e6:
                    results, stopped = find_panels(options, key_markers, sim, top=200)
                    if stopped:
                        results += local_panels(options, key_markers, sim)
                else:
                    # far too many to try one by one: build and improve panels instead
                    results = local_panels(options, key_markers, sim)
                results.sort(key=lambda r: rank(r[0], r[1]))
                # panels that differ only in one or two low-priority markers are not real options
                results = diverse(results, n=10, min_diff=max(1, round(len(options) / 4)))
            if not results:
                st.error("No panel found. Two markers may only have the same dye.")
            else:
                table = []
                for key_max, all_max, panel in results:
                    fl = list(panel.values())
                    row = {"key pairs max": key_max, "all pairs max": all_max,
                           "complexity": complexity_index(sig, fl)}
                    row.update(panel)
                    table.append(row)
                table = pd.DataFrame(table)
                order = sorted(table.index, key=lambda i: (rank(table.at[i, "key pairs max"],
                                                                table.at[i, "all pairs max"]),
                                                           table.at[i, "complexity"]))
                table = table.loc[order].head(10).reset_index(drop=True)
                table.index = table.index + 1
                st.session_state["results"] = table
                st.session_state["suggested"] = suggested
                st.session_state["levels"] = {name: level for name, _, level in rows}
                st.session_state["options"] = options

    if "results" in st.session_state:
        table = st.session_state["results"]
        markers = [c for c in table.columns if c not in ["key pairs max", "all pairs max", "complexity"]]
        st.subheader("Best panels")
        st.dataframe(table.style.format({"key pairs max": "{:.2f}", "all pairs max": "{:.2f}",
                                         "complexity": "{:.1f}"}), width="stretch")
        st.download_button("Download table", table.to_csv().encode(), "panels.csv")

        pick = st.selectbox("Look at panel", table.index)
        row = table.loc[pick]
        fluors = [row[m] for m in markers]
        labels = [f"{m} ({row[m]})" for m in markers]
        to_buy = [m for m in markers if m in st.session_state.get("suggested", set())]
        alts = alternatives({m: row[m] for m in markers}, st.session_state.get("options", {}), sim)
        if to_buy:
            st.write("**To buy** for this panel (suggested colours)")
            levels = st.session_state.get("levels", {})
            buy = pd.DataFrame([{"marker": m, "dye": row[m], "importance": levels.get(m, ""),
                                 "brightness (1-5)": brightness(row[m]),
                                 "also fits": ", ".join(d for d, _ in alts.get(m, []))} for m in to_buy])
            st.dataframe(buy, hide_index=True)
            st.download_button("Download shopping list", buy.to_csv(index=False).encode(),
                               "to_buy.csv", key="buy")
            st.caption("Check the vendor actually sells your clone in that colour; if not, try one "
                       "from 'also fits' (each keeps every pair below 0.90 with the rest of this panel). "
                       "Brightness is a rough guide from vendor charts.")
        show_panel(fluors, labels, sig, sim)

with tab_check:
    st.write("Already have a panel? Pick its dyes to see which pairs are too close.")
    chosen = st.multiselect("Dyes in your panel", all_dyes,
                            default=["DAPI", "GFP (EGFP)", "PE", "APC"])
    if len(chosen) >= 2:
        show_panel(chosen, chosen, sig, sim)
