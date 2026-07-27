from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional
from datetime import date

logger = logging.getLogger(__name__)


def save_generalisation_report_pdf(
    report_entries: List[Dict[str, Any]],
    pdf_path: Path,
    dataset_label: Optional[str] = None,
) -> None:
    """
    Save a PDF report for Step5 generalisation.

    Each entry should contain:
      - row_index
      - base_chebi
      - path_labels (or path_ids)
      - mapped_chebi
      - mapped_chebi_name (optional)
      - pathways (list[str])
    """
    mapped_entries = [
        e for e in report_entries
        if e.get("mapped_chebi") and e.get("pathways")
    ]
    if not mapped_entries:
        logger.warning("No mapped Step5 entries to write to PDF.")
        return

    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.pdfgen import canvas
    except ImportError as e:
        raise RuntimeError(
            "reportlab is required to generate the Step5 PDF report. "
            "Install it with: pip install reportlab"
        ) from e

    pdf_path = Path(pdf_path)
    c = canvas.Canvas(str(pdf_path), pagesize=A4)
    width, height = A4

    def wrap_text(text: str, max_width: float, font_name: str = "Helvetica", font_size: int = 10):
        words = text.split()
        if not words:
            return [""]
        lines = []
        current = words[0]
        for w in words[1:]:
            cand = current + " " + w
            if c.stringWidth(cand, font_name, font_size) <= max_width:
                current = cand
            else:
                lines.append(current)
                current = w
        lines.append(current)
        return lines

    left = 50
    right = 50
    bottom = 80
    max_w = width - left - right

    # Header
    y = height - 50
    c.setFont("Helvetica-Bold", 16)
    c.drawString(left, y, "Reactome is_a Generalisation Report (Step5)")
    y -= 24

    c.setFont("Helvetica", 11)
    c.drawString(left, y, f"Dataset: {dataset_label or 'Unknown'}")
    y -= 14
    c.drawString(left, y, f"Entries: {len(mapped_entries)}")
    y -= 14
    c.drawString(left, y, f"Date: {date.today().isoformat()}")
    y -= 24

    current_row = None
    c.setFont("Helvetica", 10)

    for entry in mapped_entries:
        row_idx = entry.get("row_index")
        base_chebi = entry.get("base_chebi")
        mapped_chebi = entry.get("mapped_chebi")
        mapped_name = entry.get("mapped_chebi_name")
        pathways = entry.get("pathways", [])

        path_labels = entry.get("path_labels") or entry.get("path_ids") or []
        route = " -> ".join(path_labels) if path_labels else str(base_chebi)

        if current_row is None or row_idx != current_row:
            current_row = row_idx
            y -= 10
            if y < bottom:
                c.showPage()
                y = height - 50
                c.setFont("Helvetica", 10)
            c.setFont("Helvetica-Bold", 11)
            c.drawString(left, y, f"Row {row_idx}")
            y -= 18
            c.setFont("Helvetica", 10)

        if y < bottom:
            c.showPage()
            y = height - 50
            c.setFont("Helvetica", 10)

        c.drawString(left + 10, y, f"Start CHEBI: {base_chebi}")
        y -= 12

        for line in wrap_text("Path: " + route, max_w - 20):
            if y < bottom:
                c.showPage()
                y = height - 50
                c.setFont("Helvetica", 10)
            c.drawString(left + 20, y, line)
            y -= 12

        label = f"Mapped ancestor: {mapped_chebi}"
        if mapped_name:
            label += f" ({mapped_name})"
        c.drawString(left + 20, y, label)
        y -= 12

        pw_text = "Reactome pathways: " + ", ".join(pathways)
        for line in wrap_text(pw_text, max_w - 20):
            if y < bottom:
                c.showPage()
                y = height - 50
                c.setFont("Helvetica", 10)
            c.drawString(left + 20, y, line)
            y -= 12

        y -= 6

    c.save()
    logger.info(f"Step5 PDF report saved to: {pdf_path}")
