"""Rendering: trees and the renormalized-equation report.

``render(eq, fmt)`` (fmt = text/markdown/json/latex) and the single-tree
drawers ``shorthand`` / ``ascii_art`` / ``forest``.  See notes/output.md.
"""
from .tree import shorthand, ascii_art, forest, coord_names
from .report import render, text_report, markdown_report, json_report
from .export import export_structure, structure_json, tree_to_dict, tree_from_dict

__all__ = ["render", "text_report", "markdown_report", "json_report",
           "shorthand", "ascii_art", "forest", "coord_names",
           "export_structure", "structure_json", "tree_to_dict", "tree_from_dict"]
