# Spectral Panel Picker

A small Streamlit app for choosing fluorochromes on a spectral flow cytometer
(set up for the BD FACSDiscover S8).

You list the dyes you have for each marker. The app picks one dye per marker so that
the most similar pair of dyes in the panel is as different as possible. Markers ticked
as **high** importance are placed first: the app keeps their dyes far away from everything else.

## Use it online

**https://starfold146.github.io/flow-panel-picker/**

Nothing to install. The app runs entirely in your browser (via [stlite](https://github.com/whitphx/stlite)),
so any file you upload stays on your computer. The first load takes about half a minute.

## Run it locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

## What it does

- **Pick a panel**: for every marker, choose the dyes you have and set its importance (high = key).
  The app searches all combinations (one dye per marker, no dye used twice) and ranks them by
  1. how close the most similar pair in the whole panel is (below 0.90, 0.90–0.98, above 0.98),
     so a panel with a pair that cannot be unmixed never comes first,
  2. the same band for pairs that involve a high-importance marker,
  3. a balance of the two in steps of 0.05: high-importance pairs count fully, other pairs get
     0.15 of slack, so cleaner key markers win unless the rest of the panel gets clearly worse,
  4. the most similar pair in the whole panel, then the complexity index.
  When you look at a panel it says which pairs are too close (above 0.98) or need care (0.90–0.98).
- **Suggest colours to buy**: leave a marker's dyes empty and the app suggests one.
  Each marker has an importance level: **high** (key or dim marker) gets bright, clean dyes;
  **medium** gets moderately bright ones; **low** (bright marker) is happy with a dim dye, which
  leaves the bright ones for the markers that need them. Live/dead markers get viability dyes.
  By default suggestions only use the classic colours every vendor sells (BV421-785, FITC/AF488, PE and
  its Cy/Dazzle tandems, PerCP-Cy5.5, APC, AF647, AF700, APC-Cy7, Pacific Blue); the sidebar can add
  vendor variants (eFluor, Fire, PerCP-eFluor 710 ...) or every dye including BD-only BUV/BB, and live/dead markers get fixable dyes when the panel is
  fixed and permeabilised (ticked automatically when a marker such as FOXP3, Ki67 or a cytokine is
  in the panel), otherwise the cheaper DAPI, 7-AAD, PI or SYTOX are suggested. With many markers the app builds panels and improves them one marker at a
  time instead of trying every combination, and shows panels that differ in at least a quarter of
  their markers. The chosen panel comes with a shopping list, with other dyes that would also fit
  each marker in case your clone is not sold in the suggested colour. Brightness is a rough 1–5 guide from vendor charts;
  check that your clone is sold in that colour.
- **Check a panel**: pick the dyes of an existing panel and see the similarity matrix,
  the most similar pairs and the signatures.

Similarity is the cosine similarity between two signatures (1 = identical).
Rough guide: ≥ 0.98 is very hard to unmix; 0.90–0.98 needs care, especially for dim or co-expressed markers.
The complexity index is the condition number of the signature matrix (lower is easier to unmix).

## Where the signatures come from

By default the signatures are **simulated**, not measured. Excitation and emission spectra of
120 dyes (common antibody conjugates including the cheaper Alexa Fluor, eFluor, DyLight and Cy dyes, viability dyes and cell tracking dyes) were downloaded from [FPbase](https://www.fpbase.org) (`fetch_spectra.py`).
For each S8 laser (349, 405, 488, 561, 637 nm) the excitation at the laser line is multiplied by
the emission collected in 20 nm bins up to 850 nm.

This is good enough to rank dyes against each other, but the numbers will not match the
instrument exactly: the real detector filters, laser powers and detector gains are not modelled,
and some pairs (for example APC vs Alexa Fluor 647) come out less similar than on a real machine.

If you have single-stain reference controls, export the signatures (one row per dye, one column
per detector) and upload the CSV in the sidebar. The app will use your data instead.

## Not included yet

- Real antigen density (brightness matching uses the importance level only)
- Spillover spreading from real data
- Tandem dye degradation

## License

- Code: MIT (see `LICENSE`).
- `data/spectra.csv`: spectra from [FPbase](https://www.fpbase.org), licensed under
  [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). Signatures simulated from it are shared under the same license.
