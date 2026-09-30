import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from panel import (complexity_index, find_panels, signatures_from_spectra,
                   similarity_matrix, worst_pairs)

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
    if v >= 0.98:
        return "background-color: #f4a6a6"
    if v >= 0.90:
        return "background-color: #fde2a7"
    return ""


def show_panel(fluors, labels, sig, sim):
    """Heatmap, worst pairs and signature plot for one set of dyes."""
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

st.sidebar.markdown("""
**Colour guide**
- similarity ≥ 0.98: very hard to unmix
- 0.90 – 0.98: be careful, avoid for dim or co-expressed markers
""")

# ---------------- main ----------------
st.title("Spectral Panel Picker")
st.write("Tell it which dyes you have for each marker. It picks one dye per marker so that "
         "the most similar pair in the panel is as different as possible. "
         "Tick **key** for the markers that matter most: they get the cleanest dyes first.")

tab_pick, tab_check = st.tabs(["Pick a panel", "Check a panel"])

with tab_pick:
    n = st.number_input("Number of markers", min_value=2, max_value=40, value=6)

    example = {
        0: ("CD4", ["BUV737", "BV605", "PE-Cy7"], False),
        1: ("CD8", ["BV711", "PerCP-Cy5.5", "BUV805"], False),
        2: ("CD25", ["BV421", "PE", "APC", "BB515"], True),
        3: ("FOXP3", ["PE", "APC", "Alexa Fluor 488", "eFluor 450"], True),
        4: ("Helios", ["PE-Cy7", "Alexa Fluor 647", "FITC"], False),
        5: ("Live/Dead", ["Zombie NIR", "Zombie Aqua", "LIVE/DEAD Blue"], False),
    }

    options, key_markers = {}, set()
    for i in range(int(n)):
        name0, dyes0, key0 = example.get(i, (f"Marker {i + 1}", [], False))
        c1, c2, c3 = st.columns([1.2, 4, 0.6])
        name = c1.text_input("Marker", value=name0, key=f"name{i}", label_visibility="collapsed")
        dyes = c2.multiselect("Dyes", all_dyes, default=[d for d in dyes0 if d in all_dyes],
                              key=f"dyes{i}", label_visibility="collapsed",
                              placeholder="dyes you have for this marker")
        key = c3.checkbox("key", value=key0, key=f"key{i}")
        if dyes:
            options[name] = dyes
            if key:
                key_markers.add(name)

    if st.button("Find panels", type="primary"):
        if len(options) < 2:
            st.warning("Give at least two markers some dyes.")
        elif len(options) < int(n):
            st.warning("Some markers have no dyes yet.")
        else:
            with st.spinner("Searching..."):
                results, stopped = find_panels(options, key_markers, sim, top=20)
            if stopped:
                st.info("Too many combinations, stopped early. The results are good but may not be the best.")
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
                table = pd.DataFrame(table).sort_values(
                    ["key pairs max", "all pairs max", "complexity"]).head(10).reset_index(drop=True)
                table.index = table.index + 1
                st.session_state["results"] = table

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
        show_panel(fluors, labels, sig, sim)

with tab_check:
    st.write("Already have a panel? Pick its dyes to see which pairs are too close.")
    chosen = st.multiselect("Dyes in your panel", all_dyes,
                            default=["DAPI", "GFP (EGFP)", "PE", "APC"])
    if len(chosen) >= 2:
        show_panel(chosen, chosen, sig, sim)
