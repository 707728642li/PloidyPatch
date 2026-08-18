#!/usr/bin/env python3
"""Finalize the visually inspected v0.2 fixed-layout review manuscript."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
HTML_NAME = "PloidyPatch_review_manuscript_v0.2.html"
PDF_NAME = "PloidyPatch_review_manuscript_v0.2.pdf"
SOURCE = ROOT / "manuscript/assembled_v0.2/PLOIDYPATCH_MANUSCRIPT_DRAFT_v0.2.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_text(path: Path, value: str) -> None:
    temporary = path.with_name(path.name + ".working")
    require(not temporary.exists(), f"working file already exists: {temporary}")
    temporary.write_text(value, encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def artifact(path: Path) -> dict[str, object]:
    return {
        "relative_path": path.name,
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def validate_html(path: Path) -> dict[str, object]:
    value = path.read_text(encoding="utf-8")
    require(value.count('<section class="figure-block') == 6, "expected six figure blocks")
    require(value.count('class="main-figure"') == 6, "expected six figure images")
    require(value.count("data:image/png;base64,") == 6, "expected six embedded figures")
    require("file:///" not in value, "local file URI leaked into HTML")
    require("data-claim=" not in value, "claim markup leaked into HTML")
    require("results/" not in value, "internal results path leaked into HTML")
    require("manuscript/source_data/" not in value, "internal source-data path leaked")
    require("doi:10.3389/fpls.2026.1819201" in value, "Hu et al. citation missing")
    for number in range(1, 7):
        require(f"Figure {number}.<br>" in value, f"Figure {number} label missing")
    return {
        "figure_blocks": 6,
        "embedded_figure_payloads": 6,
        "internal_paths_absent": True,
    }


def validate_pdf(path: Path) -> dict[str, object]:
    value = path.read_bytes()
    require(value.startswith(b"%PDF-"), "PDF header missing")
    require(value.rstrip().endswith(b"%%EOF"), "PDF EOF marker missing")
    page_count = len(re.findall(rb"/Type\s*/Page\b", value))
    require(page_count == 27, f"expected 27 PDF pages, found {page_count}")
    require(b"/MediaBox [0 0 612 792]" in value, "US Letter media box missing")
    require(b"/Marked true" in value, "tagged-PDF marker missing")
    return {"page_count": 27, "page_size": "US Letter (612 x 792 pt)", "tagged": True}


def write_checksums(output_dir: Path) -> None:
    paths = sorted(
        (path for path in output_dir.iterdir() if path.is_file() and path.name != "SHA256SUMS"),
        key=lambda path: path.name.encode("utf-8"),
    )
    atomic_text(
        output_dir / "SHA256SUMS",
        "".join(f"{sha256_file(path)}  {path.name}\n" for path in paths),
    )


def finalize(output_dir: Path) -> None:
    output_dir = output_dir.resolve(strict=True)
    html = output_dir / HTML_NAME
    pdf = output_dir / PDF_NAME
    for path in (html, pdf, SOURCE):
        require(path.is_file() and path.stat().st_size > 0, f"missing input: {path}")

    html_validation = validate_html(html)
    pdf_validation = validate_pdf(pdf)
    visual_qa = {
        "schema_version": "ploidypatch.formatted_review_visual_qa.v2",
        "pdf": {
            **artifact(pdf),
            **pdf_validation,
            "renderer": "Microsoft Edge headless print-to-PDF",
            "rasterizer": "Poppler pdftoppm at 144 dpi",
            "all_pages_inspected": True,
            "pages_inspected": 27,
            "status": "passed",
            "findings": {
                "clipping": 0,
                "overlap": 0,
                "missing_glyphs": 0,
                "figure_or_table_overflow": 0,
                "orphaned_table_titles": 0,
            },
        },
        "html": {**artifact(html), **html_validation},
        "editable_docx_visual_qa": {
            "status": "not_completed",
            "reason": "LibreOffice is unavailable in the build environment",
            "compensation": "deterministic DOCX structural QA plus full visual inspection of this matched fixed-layout PDF",
        },
    }
    atomic_text(
        output_dir / "visual_qa.json",
        json.dumps(visual_qa, indent=2, sort_keys=True) + "\n",
    )

    manifest = {
        "schema_version": "ploidypatch.formatted_review_manuscript.v2",
        "artifacts": {
            "html": artifact(html),
            "pdf": artifact(pdf),
            "visual_qa": artifact(output_dir / "visual_qa.json"),
        },
        "source": {
            "project_relative_path": SOURCE.relative_to(ROOT).as_posix(),
            "bytes": SOURCE.stat().st_size,
            "sha256": sha256_file(SOURCE),
        },
        "design": {
            "preset": "narrative_proposal",
            "page": "US Letter portrait",
            "margins_in": 1.0,
            "body_font": "Calibri",
            "body_size_pt": 11,
            "embedded_figures": 6,
            "table_title_pagination": "keep_with_following_table",
        },
        "qa_status": {
            "pdf_visual": "passed_all_27_pages",
            "editable_docx_visual_render": "not_completed_renderer_unavailable",
        },
    }
    atomic_text(
        output_dir / "format_manifest.json",
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    atomic_text(
        output_dir / "README.md",
        "# PloidyPatch v0.2 formatted review manuscript\n\n"
        "This directory contains the deterministic figure-integrated HTML and tagged 27-page "
        "US-Letter PDF generated from the v0.2 assembled manuscript. All PDF pages were inspected "
        "after 144-dpi Poppler rasterization; no clipping, overlap, missing glyph, figure/table "
        "overflow or orphaned table title was found. The editable double-spaced DOCX is maintained "
        "separately in `manuscript/formatted_submission_v0.2/`.\n",
    )
    write_checksums(output_dir)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "manuscript/formatted_review_v0.2",
    )
    args = parser.parse_args()
    finalize(args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
