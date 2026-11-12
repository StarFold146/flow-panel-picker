"""A tidy Excel file of the panels: one sheet per panel, sorted by laser, coloured by laser."""
import io
import re

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side

from panel import LASERS, brightness

LASER_NAME = {"UV": "UV 349", "V": "V 405", "B": "B 488", "YG": "YG 561", "R": "R 637"}
LASER_FILL = {"UV": "EDE7F6", "V": "E3E9F7", "B": "E2F2E6", "YG": "FBF1DC", "R": "F8E1E1"}
HEAD_FILL = "404040"
THIN = Side(style="thin", color="D0D0D0")


def peak(sig, dye):
    """(laser, detector number) where a dye's signature is highest."""
    det = sig.loc[dye].idxmax()
    laser, num = re.match(r"([A-Z]+)(\d+)", det).groups()
    return laser, int(num)


def _header(ws, row, names, widths):
    for col, (name, w) in enumerate(zip(names, widths), start=1):
        c = ws.cell(row=row, column=col, value=name)
        c.font = Font(bold=True, color="FFFFFF")
        c.fill = PatternFill("solid", fgColor=HEAD_FILL)
        c.alignment = Alignment(vertical="center")
        ws.column_dimensions[c.column_letter].width = w


FOOTER = [
    "Spectral similarities are simulated from public spectra (FPbase), not measured on the instrument.",
    "Please confirm the clone is sold in that colour before ordering.",
    "Titration, single-stain controls and FMOs are still needed.",
]


def workbook(table, markers, sig, sim, key_markers=(), to_buy=(), note=None, title="Panel",
             n_sheets=5, bench=None):
    """Excel bytes. table: the results table (one row per panel, a column per marker).

    note(marker, dye) -> text for the Note column; bench({marker: dye}) -> list of bench notes.
    """
    wb = Workbook()
    summary = wb.active
    summary.title = "Overview"
    summary["A1"] = f"{title}: best panels ({len(markers)} colours)"
    summary["A1"].font = Font(bold=True, size=14)
    cols = ["Panel", "Most similar pair", "Key marker pairs", "Complexity"] + list(markers)
    _header(summary, 3, cols, [8, 16, 16, 12] + [16] * len(markers))
    for r, (i, row) in enumerate(table.iterrows(), start=4):
        vals = [f"Panel {i}", round(row["all pairs max"], 2), round(row["key pairs max"], 2),
                round(row["complexity"], 1)] + [row[m] for m in markers]
        for col, v in enumerate(vals, start=1):
            c = summary.cell(row=r, column=col, value=v)
            c.border = Border(bottom=THIN)
            if col > 4:
                laser = peak(sig, v)[0]
                c.fill = PatternFill("solid", fgColor=LASER_FILL[laser])
    summary.freeze_panes = "B4"
    end = 4 + len(table) + 1
    summary.cell(row=end, column=1, value="Cells are coloured by the laser that excites the dye best. "
                 "Lower similarity and complexity are better.").font = Font(italic=True, color="808080")

    order = list(LASERS)
    for i, row in list(table.iterrows())[:n_sheets]:
        ws = wb.create_sheet(f"Panel {i}")
        dyes = {m: row[m] for m in markers}
        ws["A1"] = f"{title} {i}, BD FACSDiscover S8 ({len(markers)} colours)"
        ws["A1"].font = Font(bold=True, size=14)
        _header(ws, 3, ["Laser", "Marker", "Fluorochrome", "Brightness", "Clone", "Vendor", "Note"],
                [10, 20, 30, 11, 14, 12, 40])
        lines = sorted(markers, key=lambda m: (order.index(peak(sig, dyes[m])[0]), peak(sig, dyes[m])[1]))
        for r, m in enumerate(lines, start=4):
            laser = peak(sig, dyes[m])[0]
            key = m in key_markers
            notes = []
            if m in to_buy:
                notes.append("to buy (suggested)")
            if note and note(m, dyes[m]):
                notes.append(note(m, dyes[m]))
            vals = [LASER_NAME[laser], m + (" ★" if key else ""), dyes[m],
                    "●" * brightness(dyes[m]), "", "", "; ".join(notes)]
            for col, v in enumerate(vals, start=1):
                c = ws.cell(row=r, column=col, value=v)
                c.fill = PatternFill("solid", fgColor=LASER_FILL[laser])
                c.font = Font(bold=key)
                c.border = Border(bottom=THIN)
        r = 4 + len(lines) + 1
        ws.cell(row=r, column=1,
                value=f"★ key markers. Most similar pair {row['all pairs max']:.2f} (simulated); "
                      f"complexity index {row['complexity']:.1f}. Clone and vendor left blank to fill in.")
        for k, line in enumerate(FOOTER, start=r + 1):
            ws.cell(row=k, column=1, value=line).font = Font(italic=True, color="808080")
        if bench:
            k += 2
            ws.cell(row=k, column=1, value="Bench notes").font = Font(bold=True)
            for line in bench(dyes):
                k += 1
                ws.cell(row=k, column=1, value="• " + line)
        ws.freeze_panes = "A4"

    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()
