# Spectral Panel Picker

A small Streamlit app for choosing fluorochromes on a spectral flow cytometer
(set up for the BD FACSDiscover S8).

You list the dyes you have for each marker. The app picks one dye per marker so that
the most similar pair of dyes in the panel is as different as possible. Markers ticked
as **key** are placed first: the app keeps their dyes far away from everything else.

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

- **Pick a panel**: for every marker, choose the dyes you have and tick the key markers.
  The app searches all combinations (one dye per marker, no dye used twice) and ranks them by
  1. the highest similarity among pairs that involve a key marker,
  2. the highest similarity in the whole panel,
  3. the complexity index.
- **Check a panel**: pick the dyes of an existing panel and see the similarity matrix,
  the most similar pairs and the signatures.

Similarity is the cosine similarity between two signatures (1 = identical).
Rough guide: ≥ 0.98 is very hard to unmix; 0.90–0.98 needs care, especially for dim or co-expressed markers.
The complexity index is the condition number of the signature matrix (lower is easier to unmix).

## Where the signatures come from

By default the signatures are **simulated**, not measured. Excitation and emission spectra of
70 common dyes were downloaded from [FPbase](https://www.fpbase.org) (`fetch_spectra.py`).
For each S8 laser (349, 405, 488, 561, 637 nm) the excitation at the laser line is multiplied by
the emission collected in 20 nm bins up to 850 nm.

This is good enough to rank dyes against each other, but the numbers will not match the
instrument exactly: the real detector filters, laser powers and detector gains are not modelled,
and some pairs (for example APC vs Alexa Fluor 647) come out less similar than on a real machine.

If you have single-stain reference controls, export the signatures (one row per dye, one column
per detector) and upload the CSV in the sidebar. The app will use your data instead.

## Not included yet

- Matching bright dyes to dim antigens (antigen density)
- Spillover spreading from real data
- Tandem dye degradation

## License

- Code: MIT (see `LICENSE`).
- `data/spectra.csv`: spectra from [FPbase](https://www.fpbase.org), licensed under
  [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). Signatures simulated from it are shared under the same license.
