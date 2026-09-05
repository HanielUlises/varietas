# Configuration file for the Sphinx documentation builder.
#
# The documentation is hand-written reStructuredText using the C++ domain
# rather than Doxygen output: the library is header-only and templated, and
# the parts of it worth documenting are the constructions, not the signatures.

project = "varietas"
copyright = "2026, Haniel Ulises"
author = "Haniel Ulises"

version = "0.1"
release = "0.1.0"

extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.githubpages",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "requirements.txt"]

# -- C++ domain ---------------------------------------------------------------

primary_domain = "cpp"
highlight_language = "cpp"
cpp_index_common_prefix = ["varietas::", "varietas::urdf_import::"]

# -- HTML output --------------------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_title = "varietas"
html_short_title = "varietas"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 3,
    "sticky_navigation": True,
    "titles_only": False,
    "prev_next_buttons_location": "both",
    "style_external_links": True,
}

html_context = {
    "display_github": True,
    "github_user": "HanielUlises",
    "github_repo": "varietas",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

# -- MathJax ------------------------------------------------------------------

mathjax3_config = {
    "tex": {
        "macros": {
            "Q": r"\mathbb{Q}",
            "R": r"\mathbb{R}",
            "C": r"\mathbb{C}",
            "F": r"\mathbb{F}",
            "A": r"\mathbb{A}",
            "V": r"\mathbf{V}",
            "I": r"\mathbf{I}",
            "LT": r"\mathrm{LT}",
            "p": r"\boldsymbol{p}",
        }
    }
}

# -- Other outputs ------------------------------------------------------------

latex_elements = {
    "papersize": "a4paper",
    "pointsize": "11pt",
    "preamble": r"\usepackage{amsmath}\usepackage{amssymb}",
}

intersphinx_mapping = {}
todo_include_todos = False
