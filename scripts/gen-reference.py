#!/usr/bin/env python3
"""Regenerate the hardcoded API/crate reference tables from code sources.

Sources:
  - ``crates/*/Cargo.toml`` + ``bindings/*/Cargo.toml`` ([package] description
    and workspace crate dep edges) -> ``docs/reference/crate-overview.md``.
  - ``python/nucleide/_internal.pyi`` (typed signatures + class docstrings),
    ``python/nucleide/*.py`` facades (module docstrings + ``__all__`` +
    alias assignments), and ``python/nucleide/__init__.py`` (submodule list)
    -> ``docs/reference/python-api.md``.

The two markdown files keep hand-written prose outside marked regions. This
script only rewrites the bodies delimited by::

    <!-- GEN:<name>:START -->
    ...generated...
    <!-- GEN:<name>:END -->

Usage:
    python3 scripts/gen-reference.py --write   # regenerate the docs
    python3 scripts/gen-reference.py --check   # exit 1 with a diff if stale

Stdlib only, runs on Python 3.10+ (no ``tomllib``; manifests are parsed with
regex/line parsing and signatures use ``ast.unparse``).
"""

import argparse
import ast
import difflib
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PYTHON_DIR = REPO_ROOT / "python" / "nucleide"
PYI_PATH = PYTHON_DIR / "_internal.pyi"
CRATE_OVERVIEW_MD = REPO_ROOT / "docs" / "reference" / "crate-overview.md"
PYTHON_API_MD = REPO_ROOT / "docs" / "reference" / "python-api.md"

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
GEN_RE = re.compile(
    r"<!-- GEN:(?P<name>[A-Za-z0-9_-]+):START -->"
    r"(?P<body>.*?)"
    r"<!-- GEN:(?P=name):END -->",
    re.DOTALL,
)


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
    """Escape table-cell text so markdownlint emphasis checks stay quiet."""
    return text.replace("|", "\\|").replace("_", "\\_")


def render_crate_table(crates: list[dict[str, object]]) -> str:
    lines = [
        "| Crate | Path | Responsibility | Depends on |",
        "| --- | --- | --- | --- |",
    ]
    for crate in crates:
        name = str(crate["name"])
        path = str(crate["path"])
        desc = escape_cell(str(crate["description"]).strip())
        deps = list(crate["deps"])  # type: ignore[arg-type]
        depends = ", ".join(f"`{d}`" for d in deps) if deps else "-"
        lines.append(f"| `{name}` | `{path}` | {desc} | {depends} |")
    return "\n".join(lines)


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
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if (
                isinstance(target, ast.Name)
                and isinstance(node.value, ast.Name)
                and target.id != node.value.id
            ):
                aliases[target.id] = node.value.id
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
    }


def parse_pyi() -> tuple[dict[str, str], dict[str, dict[str, object]]]:
    """Parse _internal.pyi into function sigs and class info."""
    tree = ast.parse(PYI_PATH.read_text(encoding="utf-8"))
    funcs: dict[str, str] = {}
    classes: dict[str, dict[str, object]] = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            funcs[node.name] = _func_sig(node.name, node, strip_self=False)
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
            }
    return funcs, classes


def parse_data_py() -> tuple[dict[str, str], dict[str, str]]:
    """Parse pure-Python data.py into function sigs and constant bullets."""
    path = PYTHON_DIR / "data.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    funcs: dict[str, str] = {}
    consts: dict[str, str] = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            funcs[node.name] = _func_sig(node.name, node, strip_self=False)
        elif isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name) and target.id.isupper():
                consts[target.id] = f"{target.id} = {_unparse(node.value)}"
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


def _long_bullet(name: str, code: str, suffix: str, indent: str) -> list[str]:
    """Render a signature too long for one line as a fenced python block."""
    head = f"{indent}- `{name}`"
    if suffix:
        head += suffix
    head += ":"
    fence = f"{indent}  ```python"
    body = f"{indent}  {code}"
    close = f"{indent}  ```"
    return ["", head, "", fence, body, close, ""]


def _sig_bullet(name: str, code: str, alias_of: str | None, indent: str = "") -> list[str]:
    suffix = f" (alias of `{alias_of}`)" if alias_of else ""
    if len(code) + len(suffix) > LONG_THRESHOLD:
        short_suffix = f" (alias of `{alias_of}`)" if alias_of else ""
        return _long_bullet(name, code, short_suffix, indent)
    return [f"{indent}- `{code}`{suffix}"]


def render_module_table(facades: dict[str, dict[str, object]]) -> str:
    lines = [
        "| Submodule | Backing crate | Contents |",
        "| --- | --- | --- |",
    ]
    for mod in MODULE_ORDER:
        facade = facades[mod]
        backed = facade["backed"]
        backing = f"`{backed}`" if backed else "— (pure Python)"
        contents = escape_cell(str(facade["contents"]).strip())
        lines.append(f"| `nucleide.{mod}` | {backing} | {contents} |")
    return "\n".join(lines)


def render_module(
    mod: str,
    facade: dict[str, object],
    pyi_funcs: dict[str, str],
    pyi_classes: dict[str, dict[str, object]],
    data_funcs: dict[str, str],
    data_consts: dict[str, str],
) -> str:
    all_names = list(facade["all"])  # type: ignore[arg-type]
    aliases = dict(facade["aliases"])  # type: ignore[arg-type]
    lines: list[str] = []
    missing: list[str] = []
    for name in all_names:
        if name in aliases:
            real = aliases[name]
            sig = pyi_funcs.get(real, data_funcs.get(real, ""))
            if not sig:
                missing.append(f"{name} (alias of unknown {real})")
                continue
            # Re-target the signature at the alias name.
            sig = real_to_alias_sig(sig, real, name)
            lines.extend(_sig_bullet(name, sig, real))
        elif name in pyi_classes:
            info = pyi_classes[name]
            doc = str(info["doc"]).strip().splitlines()
            summary = doc[0].strip() if doc else ""
            head = f"- `{name}`"
            if summary:
                head += f" — {summary}"
            lines.append(head)
            members = list(info["members"])  # type: ignore[arg-type]
            for kind, member_name, member_sig in members:
                if kind == "property":
                    lines.append(f"  - `{member_sig}` (property)")
                else:
                    lines.extend(_sig_bullet(member_name, member_sig, None, indent="  "))
            # Blank line after a class block keeps nested lists readable.
            lines.append("")
        elif name in pyi_funcs:
            lines.extend(_sig_bullet(name, pyi_funcs[name], None))
        elif name in data_funcs:
            lines.extend(_sig_bullet(name, data_funcs[name], None))
        elif name in data_consts:
            lines.append(f"- `{data_consts[name]}`")
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


def replace_regions(text: str, regions: dict[str, str]) -> str:
    def _repl(match: re.Match[str]) -> str:
        name = match.group("name")
        if name not in regions:
            raise RuntimeError(f"unexpected GEN region: {name}")
        return f"<!-- GEN:{name}:START -->\n\n{regions[name]}\n\n<!-- GEN:{name}:END -->"

    result = GEN_RE.sub(_repl, text)
    found = {m.group("name") for m in GEN_RE.finditer(text)}
    wanted = set(regions)
    if found != wanted:
        raise RuntimeError(
            f"GEN region mismatch: file has {sorted(found)}, expected {sorted(wanted)}"
        )
    return result


def build_crate_regions() -> dict[str, str]:
    crates = load_crates()
    return {"crate-table": render_crate_table(crates)}


def build_python_regions() -> dict[str, str]:
    facades = {mod: parse_facade(mod) for mod in MODULE_ORDER}
    check_init_submodules()
    pyi_funcs, pyi_classes = parse_pyi()
    data_funcs, data_consts = parse_data_py()
    regions = {"module-table": render_module_table(facades)}
    for mod in MODULE_ORDER:
        regions[f"module-{mod}"] = render_module(
            mod, facades[mod], pyi_funcs, pyi_classes, data_funcs, data_consts
        )
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
