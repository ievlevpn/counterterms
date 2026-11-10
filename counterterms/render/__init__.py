"""Rendering: trees and the renormalized-equation report.

``render(eq, fmt)`` (fmt = text/markdown/json/latex/html) and the single-tree
drawers ``shorthand`` / ``ascii_art`` / ``forest`` / ``svg``.  See notes/output.md.
"""
from .tree import shorthand, ascii_art, forest, coord_names
from .svg import svg
from .report import render, text_report, markdown_report, json_report
from .html import html_document
from .export import export_structure, structure_json, tree_to_dict, tree_from_dict

__all__ = ["render", "text_report", "markdown_report", "json_report", "html_document",
           "shorthand", "ascii_art", "forest", "coord_names", "svg",
           "export_structure", "structure_json", "tree_to_dict", "tree_from_dict"]
