#!/usr/bin/env python3
"""Regenerate the hardcoded API/crate reference tables from code sources.

Sources:
  - ``crates/*/Cargo.toml`` + ``bindings/*/Cargo.toml`` ([package] description
    and workspace crate dep edges) -> ``docs/reference/crate-overview.mdx``.
  - ``python/nucleide/_internal.pyi`` (typed signatures + class docstrings),
    ``python/nucleide/*.py`` facades (module docstrings + ``__all__`` +
    alias assignments), and ``python/nucleide/__init__.py`` (submodule list)
    -> ``docs/reference/python-api.mdx``.

The two MDX files keep hand-written prose outside marked regions. This
script only rewrites the bodies delimited by::

    {/* GEN:<name>:START */}
    ...generated...
    {/* GEN:<name>:END */}

The crate overview is per-crate card sections (``### <name>`` heading, path
line with a ``[source]`` tree link, one-line responsibility, then a tag row
with ``Depends on:`` pills and a ``Python:`` pill). Pills are Markdown links
whose text is raw inline HTML ``<span style="..."><code>...</code></span>``
(plain ``[<code>](url)`` links do not pill-style); they render static with no
hydration and use ``currentColor`` so they stay readable in light and dark
themes. Each pill carries ``white-space: nowrap`` so tokens never split
mid-pill, while spaces between pills allow wrapping. Markdown (not raw
``<a>``) is required so the docs sync rewrites ``./``/``../`` ``.mdx`` hrefs
to routes; raw ``<a href>`` HTML passes through verbatim and would 404. The
module overview table stays a docs-kit ``<DataTable>`` component (3 columns
fit, so linked ``<code>`` cells never squeeze); per-module symbol bullets
stay Markdown (they render fine). Linked DataTable cells use
``React.createElement("a", ...)`` (not MDX ``<a>`` JSX, which would compile
to Astro elements and break the React island build); ``python-api.mdx``
carries ``import * as React from "react"`` outside the GEN regions.

Generated bullets link each symbol name to its GitHub source line
(``_internal.pyi`` / ``data.py`` definitions, or the facade ``alias = real``
assignment for aliases); the module cards/tables and per-module backing lines
link the docs inward to each other. Cross-page links always start with ``./``
or ``../`` so Astro resolves them to routes (``./python-api.mdx#<slug>`` for
Python pills, ``../architecture/crate-responsibilities.md#<slug>`` for crate
pills, same-page ``#<slug>`` stays as-is); source/tree links stay absolute
GitHub URLs.

Usage:
    python3 scripts/gen-reference.py --write   # regenerate the docs
    python3 scripts/gen-reference.py --check   # exit 1 with a diff if stale

Stdlib only, runs on Python 3.10+ (no ``tomllib``; manifests are parsed with
regex/line parsing and signatures use ``ast.unparse``).
"""

import argparse
import ast
import difflib
import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PYTHON_DIR = REPO_ROOT / "python" / "nucleide"
PYI_PATH = PYTHON_DIR / "_internal.pyi"
CRATE_OVERVIEW_MD = REPO_ROOT / "docs" / "reference" / "crate-overview.mdx"
PYTHON_API_MD = REPO_ROOT / "docs" / "reference" / "python-api.mdx"

# Same file base the website docs sync uses
# (website/package.json `sync-docs --github-file-base`).
GITHUB_BASE = "https://github.com/nukehub-dev/nucleide/blob/main"
GITHUB_TREE_BASE = "https://github.com/nukehub-dev/nucleide/tree/main"

CRATE_ORDER = [
    "linalg",
    "nuclei",
    "material",
    "mcnp-io",
    "serpent-io",
    "fluka-io",
    "alara-io",
    "cccc-io",
    "fispact-io",
    "origen-io",
    "enrichment",
    "depletion",
    "vr-tools",
    "r2s",
    "nucleide-bindings",
    "nucleide-wasm",
]

MODULE_ORDER = [
    "nuclei",
    "material",
    "mcnp",
    "serpent",
    "fluka",
    "vr",
    "enrichment",
    "depletion",
    "alara",
    "cccc",
    "fispact",
    "origen",
    "r2s",
    "data",
]

LONG_THRESHOLD = 110

# Markdownlint `line-length` budget (.markdownlint-cli2.jsonc); code blocks
# and tables are exempt in `.md`, and `.mdx` files are not globbed at all, so
# only bullet/paragraph lines must fit (or at least carry no wrappable spaces
# past this column, which bare GitHub URLs satisfy). Keep DataTable JSX rows
# (one per line) sane anyway for readability.
MARKDOWN_LINE_LENGTH = 120

DESC_RE = re.compile(r'^\s*description\s*=\s*"((?:[^"\\]|\\.)*)"', re.M)
NAME_RE = re.compile(r'^\s*name\s*=\s*"([^"]+)"', re.M)
# `foo.workspace = true` edge inside [dependencies].
DOTTED_DEP_RE = re.compile(r"^\s*([A-Za-z0-9_-]+)\.workspace\s*=\s*true", re.M)
# `foo = { workspace = true, ... }` edge inside [dependencies].
INLINE_DEP_RE = re.compile(
    r"^\s*([A-Za-z0-9_-]+)\s*=\s*\{[^}\n]*workspace\s*=\s*true[^}\n]*\}",
    re.M,
)
SECTION_RE = re.compile(r"^\[(?P<name>[^\]]+)\]", re.M)
BACKED_BY_RE = re.compile(r"backed by the `([^`]+)` crate")
BACKED_PAREN_RE = re.compile(r"\s*\(backed by the `[^`]+` crate\)\.?\s*$")
# GEN markers: legacy HTML comments in `.md` and MDX `{/* */}` comments in
# `.mdx` (bare `<!-- -->` does not parse as MDX). Both styles are accepted on
# read; writes always emit the MDX style since both targets are `.mdx`.
GEN_RE = re.compile(
    r"(?:<!--\s*GEN:(?P<open_a>[A-Za-z0-9_-]+):START\s*-->"
    r"|\{\s*/\*\s*GEN:(?P<open_b>[A-Za-z0-9_-]+):START\s*\*/\s*\})"
    r"(?P<body>.*?)"
    r"(?:<!--\s*GEN:(?P<close_a>[A-Za-z0-9_-]+):END\s*-->"
    r"|\{\s*/\*\s*GEN:(?P<close_b>[A-Za-z0-9_-]+):END\s*\*/\s*\})",
    re.DOTALL,
)


def github_slug(text: str) -> str:
    """Slugify heading text the way GitHub anchors do.

    Lowercase, drop backticks and punctuation except ``-``/``_``, turn spaces
    into ``-`` (e.g. ```nucleide.alara``` -> ``nucleidealara``). Used for every
    same-page and cross-page anchor this script emits.
    """
    slug = text.lower().replace("`", "")
    slug = re.sub(r"[^a-z0-9 _-]", "", slug)
    return slug.strip().replace(" ", "-")


def _pyi_url(lineno: int) -> str:
    """GitHub source link for a ``_internal.pyi`` definition line."""
    return f"{GITHUB_BASE}/python/nucleide/_internal.pyi#L{lineno}"


def _data_url(lineno: int) -> str:
    """GitHub source link for a ``data.py`` definition line."""
    return f"{GITHUB_BASE}/python/nucleide/data.py#L{lineno}"


def _facade_url(mod: str, lineno: int) -> str:
    """GitHub source link for a facade ``alias = real`` assignment line."""
    return f"{GITHUB_BASE}/python/nucleide/{mod}.py#L{lineno}"


def _tree_url(rel_path: str) -> str:
    """GitHub tree link for a repo directory (crate sources)."""
    return f"{GITHUB_TREE_BASE}/{rel_path}"


def _responsibilities_url(crate: str) -> str:
    """Relative link to a crate's section in crate-responsibilities.md."""
    return f"../architecture/crate-responsibilities.md#{github_slug(crate)}"


def _section_text(text: str, section: str) -> str:
    """Return the raw text of a `[section]` block (to the next header)."""
    matches = list(SECTION_RE.finditer(text))
    for i, match in enumerate(matches):
        if match.group("name").strip() == section:
            start = match.end()
            end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
            return text[start:end]
    return ""


def parse_manifest(path: Path) -> tuple[str, str, list[str]]:
    """Parse one Cargo.toml into (package name, description, raw dep names)."""
    text = path.read_text(encoding="utf-8")
    package = _section_text(text, "package")
    name_match = NAME_RE.search(package)
    desc_match = DESC_RE.search(package)
    if name_match is None or desc_match is None:
        raise RuntimeError(f"{path}: missing [package] name/description")
    deps_block = _section_text(text, "dependencies")
    # Scan line-by-line so dep order follows the manifest file order.
    raw: list[str] = []
    for line in deps_block.splitlines():
        dotted = DOTTED_DEP_RE.match(line)
        if dotted is not None:
            raw.append(dotted.group(1))
            continue
        inline = INLINE_DEP_RE.match(line)
        if inline is not None:
            raw.append(inline.group(1))
    # Preserve order, drop duplicates.
    seen: set[str] = set()
    deps: list[str] = []
    for dep in raw:
        if dep not in seen:
            seen.add(dep)
            deps.append(dep)
    return name_match.group(1), desc_match.group(1), deps


def load_crates() -> list[dict[str, object]]:
    """Load every workspace crate manifest with workspace-only dep edges."""
    manifests = sorted(REPO_ROOT.glob("crates/*/Cargo.toml")) + sorted(
        REPO_ROOT.glob("bindings/*/Cargo.toml")
    )
    if not manifests:
        raise RuntimeError("no crate manifests found")
    parsed = [(m, *parse_manifest(m)) for m in manifests]
    known = {name for _, name, _, _ in parsed}
    crates: list[dict[str, object]] = []
    for manifest, name, desc, raw_deps in parsed:
        deps = [d for d in raw_deps if d in known]
        rel = manifest.parent.relative_to(REPO_ROOT).as_posix()
        crates.append({"name": name, "path": rel, "description": desc, "deps": deps})
    order = {name: i for i, name in enumerate(CRATE_ORDER)}
    missing = known - set(order)
    if missing:
        raise RuntimeError(f"CRATE_ORDER missing crates: {sorted(missing)}")
    crates.sort(key=lambda c: order[str(c["name"])])
    return crates


def escape_cell(text: str) -> str:
    """Escape table-cell text so markdownlint emphasis checks stay quiet.

    Retained for history: the overview tables are now docs-kit
    ``<DataTable>`` components (plain JS strings, no Markdown escaping).
    """
    return text.replace("|", "\\|").replace("_", "\\_")


def _js_string(text: str) -> str:
    """Render text as a double-quoted JS string literal for DataTable data."""
    return json.dumps(text, ensure_ascii=False)


def _link_cell(href: str, text: str) -> str:
    """Render ``<a><code>`` as React elements for DataTable data.

    MDX ``<a>`` JSX compiles to Astro elements (``{astro:jsx, ...}``), which
    React cannot render inside the docs-kit DataTable island (Astro-React
    interop). ``React.createElement`` builds real React elements instead, so
    the site build passes and the DOM is still ``<a><code>``. Both elements
    carry ``whiteSpace: nowrap`` to override docs-kit ``.prose a``
    ``overflow-wrap: anywhere`` (which leaks into ``not-prose`` islands and
    would otherwise split identifiers like ``linalg`` mid-token); plain-string
    cells (Path/Responsibility/Depends/Contents) wrap only at spaces. Requires
    ``import * as React from "react"`` at the top of the MDX file (outside GEN
    regions);     the MDX files also import the React ``DataTable`` directly so
    ``<DataTable>`` renders static (no island hydration, which cannot serialize
    element props).
    """
    href_js = _js_string(href)
    text_js = _js_string(text)
    nowrap = '{ whiteSpace: "nowrap" }'
    inner = f'React.createElement("code", {{ style: {nowrap} }}, {text_js})'
    return f'React.createElement("a", {{ href: {href_js}, style: {nowrap} }}, {inner})'


PILL_STYLE = (
    "display:inline-block;border:1px solid currentColor;"
    "border-radius:9999px;padding:0 .6em;margin:.1em .15em;"
    "font-size:.85em;white-space:nowrap;text-decoration:none"
)


def _pill(href: str, text: str) -> str:
    """Render a pill as a Markdown link with a styled inner span.

    Plain Markdown ``[<code>](url)`` links do not pill-style, so the link
    text carries raw inline HTML ``<span style><code></code></span>`` with
    ``currentColor`` (theme-safe in light and dark). ``inline-block`` plus
    ``white-space: nowrap`` keeps each pill whole while spaces between pills
    still allow wrapping. The link itself stays Markdown (not raw ``<a>``)
    because only Markdown ``[text](href)`` links are rewritten to routes by
    the docs sync (``sync-docs.mjs``) and Astro; raw ``<a href>`` HTML passes
    through verbatim and would 404 as ``*.mdx``. ``href`` must start with
    ``./`` or ``../`` for cross-page links; absolute GitHub URLs stay as-is
    for source links. Static output, no hydration.
    """
    return f'[<span style="{PILL_STYLE}"><code>{text}</code></span>]({href})'


def render_crate_cards(crates: list[dict[str, object]], crate_to_module: dict[str, str]) -> str:
    """Render the crate overview as per-crate card sections.

    Each card is a ``### <name>`` heading, a ```<path>`` line with a
    ``[source]`` tree link, a one-line responsibility, then a tag row with
    ``Depends on:`` pills and a ``Python:`` pill (``-`` when empty).
    """
    parts: list[str] = []
    for crate in crates:
        name = str(crate["name"])
        path = str(crate["path"])
        desc = str(crate["description"]).strip()
        deps = list(crate["deps"])  # type: ignore[arg-type]
        if deps:
            pills = " ".join(_pill(_responsibilities_url(str(d)), str(d)) for d in deps)
        else:
            pills = "-"
        mod = crate_to_module.get(name)
        if mod is None:
            python = "-"
        else:
            slug = github_slug(f"nucleide.{mod}")
            python = _pill(f"./python-api.mdx#{slug}", f"nucleide.{mod}")
        parts.append(
            f"### {name}\n\n`{path}` · [source]({_tree_url(path)})\n\n"
            f"{desc}\n\nDepends on: {pills} · Python: {python}"
        )
    return "\n\n".join(parts)


def parse_facade(mod: str) -> dict[str, object]:
    """Parse a facade module docstring, __all__, and alias assignments."""
    path = PYTHON_DIR / f"{mod}.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    doc = ast.get_docstring(tree) or ""
    first_line = doc.strip().splitlines()[0].strip() if doc.strip() else ""
    backed: str | None = None
    backed_match = BACKED_BY_RE.search(doc)
    if backed_match is not None:
        backed = backed_match.group(1)
    contents = BACKED_PAREN_RE.sub("", first_line).strip()
    all_names: list[str] = []
    aliases: dict[str, str] = {}
    alias_lines: dict[str, int] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if (
                isinstance(target, ast.Name)
                and isinstance(node.value, ast.Name)
                and target.id != node.value.id
            ):
                aliases[target.id] = node.value.id
                alias_lines[target.id] = node.lineno
            if (
                isinstance(target, ast.Name)
                and target.id == "__all__"
                and isinstance(node.value, ast.List)
            ):
                for elt in node.value.elts:
                    if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                        all_names.append(elt.value)
    if not all_names:
        raise RuntimeError(f"{path}: could not parse __all__")
    return {
        "doc": doc,
        "first_line": first_line,
        "backed": backed,
        "contents": contents,
        "all": all_names,
        "aliases": aliases,
        "alias_lines": alias_lines,
    }


def parse_pyi() -> tuple[dict[str, tuple[str, int]], dict[str, dict[str, object]]]:
    """Parse _internal.pyi into function (sig, lineno) and class info.

    Class entries carry ``doc``, ``members``, and the ``class`` statement
    ``line`` so bullets can link to the GitHub source line.
    """
    tree = ast.parse(PYI_PATH.read_text(encoding="utf-8"))
    funcs: dict[str, tuple[str, int]] = {}
    classes: dict[str, dict[str, object]] = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            funcs[node.name] = (_func_sig(node.name, node, strip_self=False), node.lineno)
        elif isinstance(node, ast.ClassDef):
            members: list[tuple[str, str, str]] = []
            for item in node.body:
                if not isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    continue
                decos = {_unparse(d) for d in item.decorator_list}
                if "property" in decos:
                    ret = _unparse(item.returns) if item.returns else "Any"
                    members.append(("property", item.name, f"{item.name}: {ret}"))
                elif "staticmethod" in decos or "classmethod" in decos:
                    members.append(
                        (
                            "method",
                            item.name,
                            _func_sig(item.name, item, strip_self=False),
                        )
                    )
                else:
                    members.append(
                        (
                            "method",
                            item.name,
                            _func_sig(item.name, item, strip_self=True),
                        )
                    )
            classes[node.name] = {
                "doc": ast.get_docstring(node) or "",
                "members": members,
                "line": node.lineno,
            }
    return funcs, classes


def parse_data_py() -> tuple[dict[str, tuple[str, int]], dict[str, tuple[str, int]]]:
    """Parse pure-Python data.py into (sig, lineno) funcs and const bullets."""
    path = PYTHON_DIR / "data.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    funcs: dict[str, tuple[str, int]] = {}
    consts: dict[str, tuple[str, int]] = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            funcs[node.name] = (_func_sig(node.name, node, strip_self=False), node.lineno)
        elif isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name) and target.id.isupper():
                consts[target.id] = (f"{target.id} = {_unparse(node.value)}", node.lineno)
    return funcs, consts


def _unparse(node: ast.AST) -> str:
    return ast.unparse(node)


def _format_args(args: ast.arguments, *, strip_self: bool) -> str:
    text = _unparse(args)
    # ast.unparse renders defaults without padding (`top: str | None=None`).
    # Pad bare `=` (none of our defaults contain `=` inside string literals).
    text = re.sub(r"(?<![\s=!<>:-])=(?![=\s])", " = ", text)
    if strip_self:
        if text == "self":
            return ""
        if text.startswith("self, "):
            return text[len("self, ") :]
    return text


def _func_sig(name: str, node: ast.FunctionDef | ast.AsyncFunctionDef, *, strip_self: bool) -> str:
    args = _format_args(node.args, strip_self=strip_self)
    ret = _unparse(node.returns) if node.returns else "None"
    return f"{name}({args}) -> {ret}"


def _long_bullet(
    name: str,
    code: str,
    suffix: str,
    indent: str,
    url: str | None = None,
    trailing: list[str] | None = None,
) -> list[str]:
    """Render a signature too long for one line as a fenced python block.

    The head links just the symbol ``name`` (a full linked signature plus a
    wrappable suffix would break the markdownlint line-length budget); an
    alias note that does not fit the head goes in ``trailing`` instead.
    """
    head = f"{indent}- `{name}`" if url is None else f"{indent}- [`{name}`]({url})"
    if suffix:
        head += suffix
    head += ":"
    fence = f"{indent}  ```python"
    body = f"{indent}  {code}"
    close = f"{indent}  ```"
    lines = ["", head, "", fence, body, close]
    lines.extend(trailing or [])
    lines.append("")
    return lines


def _sig_bullet(
    name: str,
    code: str,
    alias_of: str | None,
    url: str | None = None,
    indent: str = "",
) -> list[str]:
    """Render one signature bullet, linking the symbol name to its source.

    ``url`` is None for class-member sub-bullets, which stay unlinked (one
    link per symbol). Single-line linked bullets always fit the lint budget
    unless they carry a wrappable ``(alias of ...)`` suffix past column 120;
    those wrap the suffix onto a continuation line instead.
    """
    suffix = f" (alias of `{alias_of}`)" if alias_of else ""
    if len(code) + len(suffix) > LONG_THRESHOLD:
        if url is not None and alias_of is not None:
            trailing = [f"{indent}  (alias of `{alias_of}`)"]
            return _long_bullet(name, code, "", indent, url, trailing)
        return _long_bullet(name, code, suffix, indent, url)
    if url is None:
        return [f"{indent}- `{code}`{suffix}"]
    if alias_of is None:
        return [f"{indent}- [`{code}`]({url})"]
    single = f"{indent}- [`{code}`]({url}){suffix}"
    if len(single) <= MARKDOWN_LINE_LENGTH:
        return [single]
    return [f"{indent}- [`{code}`]({url})", f"{indent}  (alias of `{alias_of}`)"]


def render_module_table(facades: dict[str, dict[str, object]]) -> str:
    """Render the Python module map as a docs-kit DataTable.

    Submodule/Backing-crate cells render ``<a><code>`` links (same-page anchors
    and relative responsibility links) via ``React.createElement``; Contents
    stay plain strings. Tables render static; ``sortable`` stays off to
    preserve MODULE_ORDER.
    """
    lines = [
        "<DataTable",
        '  caption="Python submodules, their backing crates, and contents."',
        "  columns={[",
        '    { key: "submodule", header: "Submodule" },',
        '    { key: "backing", header: "Backing crate" },',
        '    { key: "contents", header: "Contents" },',
        "  ]}",
        "  data={[",
    ]
    for mod in MODULE_ORDER:
        facade = facades[mod]
        backed = facade["backed"]
        if backed:
            crate = str(backed)
            url = _responsibilities_url(crate)
            backing = _link_cell(url, crate)
        else:
            backing = _js_string("— (pure Python)")
        slug = github_slug(f"nucleide.{mod}")
        submodule = _link_cell(f"#{slug}", f"nucleide.{mod}")
        contents = _js_string(str(facade["contents"]).strip())
        row = f"    {{ submodule: {submodule}, backing: {backing}, contents: {contents} }},"
        lines.append(row)
    lines.extend(
        [
            "  ]}",
            "/>",
        ]
    )
    return "\n".join(lines)


def render_backing_line(mod: str, facade: dict[str, object], crate_dirs: dict[str, str]) -> str:
    """One line under each module `##` heading: backing crate + source dir."""
    backed = facade["backed"]
    if backed is None:
        return f"Backing crate: pure Python · [source file]({GITHUB_BASE}/python/nucleide/data.py)"
    crate = str(backed)
    if crate not in crate_dirs:
        raise RuntimeError(f"nucleide.{mod}: backing crate `{crate}` has no manifest dir")
    return (
        f"Backing crate: [`{crate}`]({_responsibilities_url(crate)}) · "
        f"[source dir]({_tree_url(crate_dirs[crate])})"
    )


def render_module(
    mod: str,
    facade: dict[str, object],
    pyi_funcs: dict[str, tuple[str, int]],
    pyi_classes: dict[str, dict[str, object]],
    data_funcs: dict[str, tuple[str, int]],
    data_consts: dict[str, tuple[str, int]],
) -> str:
    all_names = list(facade["all"])  # type: ignore[arg-type]
    aliases = dict(facade["aliases"])  # type: ignore[arg-type]
    alias_lines = dict(facade["alias_lines"])  # type: ignore[arg-type]
    lines: list[str] = []
    missing: list[str] = []
    for name in all_names:
        if name in aliases:
            real = aliases[name]
            if real in pyi_funcs:
                sig = real_to_alias_sig(pyi_funcs[real][0], real, name)
            elif real in data_funcs:
                sig = real_to_alias_sig(data_funcs[real][0], real, name)
            else:
                missing.append(f"{name} (alias of unknown {real})")
                continue
            url = _facade_url(mod, alias_lines[name])
            lines.extend(_sig_bullet(name, sig, real, url))
        elif name in pyi_classes:
            info = pyi_classes[name]
            url = _pyi_url(int(info["line"]))
            doc = str(info["doc"]).strip().splitlines()
            summary = doc[0].strip() if doc else ""
            if summary:
                single = f"- [`{name}`]({url}) — {summary}"
                if len(single) <= MARKDOWN_LINE_LENGTH:
                    lines.append(single)
                else:
                    lines.append(f"- [`{name}`]({url})")
                    lines.append(f"  — {summary}")
            else:
                lines.append(f"- [`{name}`]({url})")
            members = list(info["members"])  # type: ignore[arg-type]
            for kind, member_name, member_sig in members:
                if kind == "property":
                    lines.append(f"  - `{member_sig}` (property)")
                else:
                    lines.extend(_sig_bullet(member_name, member_sig, None, None, "  "))
            # Blank line after a class block keeps nested lists readable.
            lines.append("")
        elif name in pyi_funcs:
            sig, lineno = pyi_funcs[name]
            lines.extend(_sig_bullet(name, sig, None, _pyi_url(lineno)))
        elif name in data_funcs:
            sig, lineno = data_funcs[name]
            lines.extend(_sig_bullet(name, sig, None, _data_url(lineno)))
        elif name in data_consts:
            text, lineno = data_consts[name]
            lines.append(f"- [`{text}`]({_data_url(lineno)})")
        else:
            missing.append(name)
    if missing:
        raise RuntimeError(f"nucleide.{mod}: exports missing from sources: " + ", ".join(missing))
    # Drop the trailing blank line the class branch may leave behind.
    while lines and lines[-1] == "":
        lines.pop()
    return "\n".join(lines)


def real_to_alias_sig(sig: str, real: str, alias: str) -> str:
    if sig.startswith(real + "("):
        return alias + sig[len(real) :]
    return sig


def _gen_name(match: re.Match[str]) -> str:
    """Return the GEN region name from a dual-style marker match."""
    name = match.group("open_a") or match.group("open_b")
    close = match.group("close_a") or match.group("close_b")
    if name is None or close is None or name != close:
        raise RuntimeError(f"mismatched GEN markers: {match.group(0)[:60]!r}")
    return name


def replace_regions(text: str, regions: dict[str, str]) -> str:
    def _repl(match: re.Match[str]) -> str:
        name = _gen_name(match)
        if name not in regions:
            raise RuntimeError(f"unexpected GEN region: {name}")
        start = f"{{/* GEN:{name}:START */}}"
        end = f"{{/* GEN:{name}:END */}}"
        return f"{start}\n\n{regions[name]}\n\n{end}"

    result = GEN_RE.sub(_repl, text)
    found = {_gen_name(m) for m in GEN_RE.finditer(text)}
    wanted = set(regions)
    if found != wanted:
        raise RuntimeError(
            f"GEN region mismatch: file has {sorted(found)}, expected {sorted(wanted)}"
        )
    return result


def build_crate_regions() -> dict[str, str]:
    crates = load_crates()
    facades = {mod: parse_facade(mod) for mod in MODULE_ORDER}
    crate_to_module = {
        str(facade["backed"]): mod for mod, facade in facades.items() if facade["backed"]
    }
    return {"crate-table": render_crate_cards(crates, crate_to_module)}


def build_python_regions() -> dict[str, str]:
    facades = {mod: parse_facade(mod) for mod in MODULE_ORDER}
    check_init_submodules()
    crate_dirs = {str(c["name"]): str(c["path"]) for c in load_crates()}
    pyi_funcs, pyi_classes = parse_pyi()
    data_funcs, data_consts = parse_data_py()
    regions = {"module-table": render_module_table(facades)}
    for mod in MODULE_ORDER:
        body = render_module(mod, facades[mod], pyi_funcs, pyi_classes, data_funcs, data_consts)
        backing = render_backing_line(mod, facades[mod], crate_dirs)
        regions[f"module-{mod}"] = f"{backing}\n\n{body}"
    return regions


def check_init_submodules() -> None:
    """Assert MODULE_ORDER matches the submodule list in __init__.py."""
    path = PYTHON_DIR / "__init__.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    imported: list[str] = []
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module == "nucleide":
            for alias in node.names:
                if alias.name != "*":
                    imported.append(alias.name)
    if set(imported) != set(MODULE_ORDER):
        raise RuntimeError(
            f"{path}: submodule list {sorted(imported)} "
            f"does not match MODULE_ORDER {sorted(MODULE_ORDER)}"
        )


def file_update(path: Path, regions: dict[str, str]) -> tuple[str, str]:
    old = path.read_text(encoding="utf-8")
    new = replace_regions(old, regions)
    if not new.endswith("\n"):
        new += "\n"
    return old, new


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Regenerate GEN regions in the reference docs.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true", help="regenerate the markdown files")
    group.add_argument(
        "--check",
        action="store_true",
        help="exit 1 with a diff when the docs are stale",
    )
    args = parser.parse_args(argv)

    targets = [
        (CRATE_OVERVIEW_MD, build_crate_regions()),
        (PYTHON_API_MD, build_python_regions()),
    ]
    stale: list[Path] = []
    for path, regions in targets:
        old, new = file_update(path, regions)
        if old == new:
            continue
        stale.append(path)
        if args.write:
            path.write_text(new, encoding="utf-8")
        else:
            diff = difflib.unified_diff(
                old.splitlines(keepends=True),
                new.splitlines(keepends=True),
                fromfile=str(path),
                tofile=str(path),
            )
            sys.stdout.writelines(diff)
    if args.write:
        if stale:
            print("updated: " + ", ".join(str(p) for p in stale))
        else:
            print("reference docs already fresh")
        return 0
    if stale:
        print(
            "stale reference docs: " + ", ".join(str(p) for p in stale),
            file=sys.stderr,
        )
        print("run: python3 scripts/gen-reference.py --write", file=sys.stderr)
        return 1
    print("reference docs fresh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
