# Spectral Panel Picker

A small Streamlit app for choosing fluorochromes on a spectral flow cytometer
(set up for the BD FACSDiscover S8).

You list the dyes you have for each marker. The app picks one dye per marker so that
the most similar pair of dyes in the panel is as different as possible. Markers ticked
as **high** importance are placed first: the app keeps their dyes far away from everything else.
Each marker's **expression** (dim / medium / bright) decides how bright its dye should be.

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

- **Pick a panel**: for every marker, choose the dyes you have, set its importance (high = key)
  and its expression level. Expression is first guessed from the name, following OMIP-105's
  antigen classes: clear on/off markers (CD3, CD4, CD8, CD45 ...) are bright, continuous ones
  (CD25, PD-1, CD44, KLRG1 ...) medium, and low or uncertain ones (FOXP3, Helios, Ki-67, CTLA-4,
  cytokines, tags such as FLAG) dim.
  The app searches all combinations (one dye per marker, no dye used twice) and ranks them by
  1. how close the most similar pair in the whole panel is (below 0.90, 0.90–0.98, above 0.98),
     so a panel with a pair that cannot be unmixed never comes first,
  2. the same band for pairs that involve a high-importance marker,
  3. a balance of the two in steps of 0.05: high-importance pairs count fully, other pairs get
     0.15 of slack, so cleaner key markers win unless the rest of the panel gets clearly worse,
     Below 0.50 every pair is comfortably apart, so differences there don't count,
  4. the **dye fit**: dim markers on bright dyes, bright markers on dimmer ones, and no PE/APC/PerCP
     tandem on an intracellular marker (FOXP3 and methanol buffers can break tandems).
     Panels that use the same colours handed out differently count as one, the best fit kept,
  5. the most similar pair in the whole panel, then the complexity index.
  Pairs between 0.70 and 0.90 are pointed out too: fine in general, but OMIP-105 kept
  co-expressed markers below 0.70. Each panel comes with bench notes (Brilliant Stain Buffer
  when there are two or more polymer dyes, a monocyte blocker with cyanine tandems,
  markers to stain before fixing, fixation time).
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
  check that your clone is sold in that colour. The sidebar links to panel builders that know
  what is sold (FluoroFinder, Biocompare, Thermo Fisher, BioLegend, BD). Tick suggested colours
  you can't buy and the app searches again without them.
- **Download**: an Excel file with an overview of all panels and one sheet per panel, sorted
  and coloured by laser, with key markers starred, empty Clone and Vendor columns to fill in,
  notes and a disclaimer. The plain CSV is still there.
- **Load / save a marker list**: a CSV with columns `marker`, `dyes` (separated by `;`, empty = suggest)
  `importance` (high / medium / low) and `expression` (dim / medium / bright, empty = guess) fills in every row at once; *Download this list* saves the
  current rows in the same format. Handy for keeping the colours you know your vendors sell.
- **Check a panel**: pick the dyes of an existing panel and see the similarity matrix,
  the most similar pairs and the signatures.

Similarity is the cosine similarity between two signatures (1 = identical).
Rough guide: ≥ 0.98 is very hard to unmix; 0.90–0.98 needs care, especially for dim or co-expressed markers.
The complexity index is the condition number of the signature matrix (lower is easier to unmix).

## Where the signatures come from

By default the signatures are **simulated**, not measured. Excitation and emission spectra of
120 dyes (common antibody conjugates including the cheaper Alexa Fluor, eFluor, DyLight and Cy dyes, viability dyes and cell tracking dyes) were downloaded from [FPbase](https://www.fpbase.org) (`fetch_spectra.py`).
For each S8 laser (349, 405, 488, 561, 637 nm) the excitation at the laser line is multiplied by
the emission collected in 20 nm bins up to 850 nm. The result is saved as `data/signatures.csv`,
which is what the app loads.

This is good enough to rank dyes against each other, but the numbers will not match the
instrument exactly: the real detector filters, laser powers and detector gains are not modelled,
and some pairs (for example APC vs Alexa Fluor 647) come out less similar than on a real machine.

If you have single-stain reference controls, export the signatures (one row per dye, one column
per detector) and upload the CSV in the sidebar. The app will use your data instead.

Dyes that are the same colour under different vendors' names are shown as one entry, for
example `PE-CF594 / PE-Dazzle 594 / PE-eFluor 610 / PE-Texas Red` or `APC-Cy7 / APC-H7 / APC-Fire 750 / APC-eFluor 780`
(only where the signatures are about 0.90 similar or more), so you can buy whichever vendor sells your clone.

## Not included yet

- Spillover spreading from real data
- Which markers are on the same cells (the 0.70 note is shown for every pair)

## License

- Code: MIT (see `LICENSE`).
- `data/spectra.csv`: spectra from [FPbase](https://www.fpbase.org), licensed under
  [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). Signatures simulated from it are shared under the same license.
