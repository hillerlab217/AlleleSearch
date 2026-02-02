#!/usr/bin/env python3
"""
colour_amino_acids.py
Colour the amino-acid alignment sheet, center cells, format columns,
and export as Excel, PDF, or both.

Usage:
    python colour_amino_acids.py ALIGNMENT.xlsx [--out {excel,pdf,both}] [--sheet SHEETNAME]

Defaults:
    --out both
    --sheet alignment
"""

import sys
import argparse
from pathlib import Path
from typing import List
from openpyxl import load_workbook
from openpyxl.styles import PatternFill, Alignment, Font, Border, Side

# PDF
from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape, A4
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle

# --------------------------- colour palette ----------------------------
AA_COLOUR = {
    "A": "f4b183", "G": "f4b183", "I": "f4b183", "L": "f4b183",
    "M": "f4b183", "V": "f4b183", "P": "f4b183",
    "R": "4f81bd", "K": "4f81bd", "H": "4f81bd",
    "N": "9bbb59", "C": "9bbb59", "Q": "9bbb59", "S": "9bbb59", "T": "9bbb59",
    "D": "e26b6b", "E": "e26b6b",
    "F": "fff475", "W": "fff475", "Y": "fff475",
}
DEFAULT_FILL = PatternFill("solid", fgColor="ffffff")
CENTER       = Alignment(horizontal="center", vertical="center")
HEADER_FONT  = Font(size=8)

LEFT = Side(style="thin")
RIGHT = Side(style="thin")
BOTTOM = Side(style="thin")

ALLELE_BORDER = Border(right=RIGHT, bottom=BOTTOM)
HEADER_BORDER = Border(left=LEFT, right=RIGHT, bottom=BOTTOM)
BLANK_BORDER = Border(bottom=BOTTOM)

# ----------------------------- excel styling ---------------------------------
def colour_sheet_inplace(xlsx_path: str, sheet_name: str = "alignment") -> None:
    wb = load_workbook(xlsx_path)
    if sheet_name not in wb.sheetnames:
        raise ValueError(f"Sheet '{sheet_name}' not found in {xlsx_path}")
    ws = wb[sheet_name]

    # hide gridlines
    ws.sheet_view.showGridLines = False

    # header row (row 1)
    for cell in ws[1]:
        if isinstance(cell.value, str):
            cell.value = cell.value.replace("_name", "")
        header_txt = str(cell.value).strip() if cell.value else ""
        # font size for B+ headers
        if cell.column_letter != "A":
            cell.font = HEADER_FONT
        cell.border = HEADER_BORDER if header_txt else BLANK_BORDER
        cell.alignment = CENTER

    cell_A1 = ws["A1"]
    cell_A1.border    = ALLELE_BORDER
    cell_A1.alignment = CENTER

    # -- column A (allele names)
    ws.column_dimensions["A"].width = 8.0
    for cell in ws["A"][1:]:
        if isinstance(cell.value, str):
            cell.value = cell.value.replace("_D39V", "")
        cell.alignment = CENTER

    # -- columns B onward
    for col_cells in ws.iter_cols(min_col=2, min_row=2):
        col_letter = col_cells[0].column_letter
        header_cell = ws[f"{col_letter}1"]
        is_gap = not (header_cell.value and str(header_cell.value).strip())

        # set width per column type
        ws.column_dimensions[col_letter].width = 0.5 if is_gap else 2.4

        # colour body cells; no borders in body
        for cell in col_cells:
            aa = str(cell.value).strip().upper() if cell.value else ""
            cell.fill = PatternFill("solid", fgColor=AA_COLOUR.get(aa, "ffffff"))
            cell.alignment = CENTER

    wb.save(xlsx_path)
    print(f"✓ Excel updated: {Path(xlsx_path).name}")

# ----------------------------- pdf export ------------------------------------
def xlsx_to_pdf(xlsx_path: str, sheet_name: str = "alignment") -> Path:
    wb = load_workbook(xlsx_path, data_only=True)
    if sheet_name not in wb.sheetnames:
        raise ValueError(f"Sheet '{sheet_name}' not found in {xlsx_path}")
    ws = wb[sheet_name]

    # convert to 2D array (strings)
    data: List[List[str]] = []
    for row in ws.iter_rows(values_only=True):
        data.append([str(c) if c is not None else "" for c in row])

    if not data:
        raise ValueError("Sheet is empty; nothing to render.")

    # figure out which columns are 'gaps' from header row (row 0)
    header = data[0]
    ncols = len(header)

    # Estimate column widths to roughly mirror Excel (A wider, gaps very thin)
    # A4 landscape width ~ 11.69in; ReportLab handles points. We'll allocate relative units.
    # Base units: allele col = 8, aa col = 2.4, gap col = 0.5 (same ratios as Excel widths)
    rel = []
    for j in range(ncols):
        if j == 0:
            rel.append(8.0)
        else:
            is_gap = (header[j].strip() == "")
            rel.append(0.5 if is_gap else 2.4)

    # Normalize to page width with margins
    page_w, page_h = landscape(A4)
    left_margin = right_margin = 24  # pts
    usable_w = page_w - left_margin - right_margin
    total_rel = sum(rel)
    colWidths = [usable_w * (r / total_rel) for r in rel]

    # Build table
    pdf_path = Path(xlsx_path).with_suffix(".pdf")
    doc = SimpleDocTemplate(str(pdf_path), pagesize=landscape(A4),
                            leftMargin=left_margin, rightMargin=right_margin,
                            topMargin=24, bottomMargin=24)
    table = Table(data, colWidths=colWidths, repeatRows=1)

    # Base style
    style = TableStyle([
        ("ALIGN", (0,0), (-1,-1), "CENTER"),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("FONTSIZE", (1,0), (-1,0), 6),   # header row small for B+ headers
        ("FONTSIZE", (0,0), (0,0), 8),    # A1 a bit larger
        ("GRID", (0,0), (-1,-1), 0.25, colors.black),
        # Excel-like emphasis lines:
        ("LINEAFTER", (0,0), (0,-1), 0.5, colors.black),   # right border after allele column A
        ("LINEBELOW", (0,0), (-1,0), 0.5, colors.black),   # bottom border under header
    ])

    # Header cleanup visuals: light grey background for header cells with text
    for c, val in enumerate(header):
        if val.strip():
            style.add("BACKGROUND", (c, 0), (c, 0), colors.whitesmoke)

    # Colour amino-acid cells
    for r in range(1, len(data)):
        row = data[r]
        for c in range(1, ncols):  # AA cells begin at column B
            aa = row[c].strip().upper()
            hexcol = AA_COLOUR.get(aa)
            if hexcol:
                style.add("BACKGROUND", (c, r), (c, r), colors.HexColor("#" + hexcol))

    table.setStyle(style)
    doc.build([table])
    print(f"✓ PDF written : {pdf_path.name}")
    return pdf_path

# ----------------------------- CLI ------------------------------------------
def main():
    p = argparse.ArgumentParser(description="Colour amino-acid alignment and export as Excel/PDF.")
    p.add_argument("xlsx", help="Path to the alignment .xlsx file (will be updated in place for Excel output).")
    p.add_argument("--out", choices=["excel", "pdf", "both"], default="both",
                   help="Choose output format (default: both).")
    p.add_argument("--sheet", default="alignment", help="Worksheet name (default: alignment).")
    args = p.parse_args()

    xlsx_path = args.xlsx
    out = args.out
    sheet = args.sheet

    if out in ("excel", "both"):
        colour_sheet_inplace(xlsx_path, sheet)

    if out in ("pdf", "both"):
        xlsx_to_pdf(xlsx_path, sheet)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print(__doc__.strip())
        sys.exit(1)
    main()