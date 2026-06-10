from __future__ import annotations

from typing import Iterable, Sequence

import pandas as pd

from .datamodels import MarkerProgram, MalignantProgram


def _as_sequence(value):
    if value is None:
        return []
    if isinstance(value, (str, MarkerProgram, MalignantProgram, dict)):
        return [value]
    return list(value)


def _iter_hierarchy_nodes(items, *, source: str):
    for item in _as_sequence(items):
        if isinstance(item, str):
            from .builtin import get_builtin_hierarchy
            yield from _iter_hierarchy_nodes(get_builtin_hierarchy(item), source=f"hierarchy:{item}")
        elif isinstance(item, MarkerProgram):
            yield item, source
            for child, child_source in _iter_hierarchy_nodes(item.children, source=source):
                yield child, child_source
        elif isinstance(item, dict):
            yield from _iter_hierarchy_nodes(MarkerProgram.from_dict(item), source=source)
        else:
            raise TypeError("hierarchy must contain MarkerProgram, dict, or built-in hierarchy name strings")


def _iter_flat_programs(items, *, source: str):
    for item in _as_sequence(items):
        if isinstance(item, str):
            from .builtin import get_builtin_marker_program_set
            yield from _iter_flat_programs(get_builtin_marker_program_set(item), source=f"program_set:{item}")
        elif isinstance(item, MalignantProgram):
            yield item, source
        elif isinstance(item, dict):
            yield MalignantProgram.from_dict(item), source
        else:
            raise TypeError("programs must contain MalignantProgram, dict, or built-in marker-program set name strings")


def _dedupe_preserve(values: Iterable[str]) -> list[str]:
    out = []
    seen = set()
    for value in values:
        if value not in seen:
            out.append(value)
            seen.add(value)
    return out


def collect_positive_marker_genes(
    hierarchy=None,
    program_sets=None,
    programs=None,
    *,
    uppercase: bool = True,
    deduplicate: bool = True,
    return_format: str = "list",
):
    """Collect positive marker genes from hierarchies and flat program sets.

    This small convenience helper is intended for sanity-check visualizations,
    such as constructing a marker heatmap from the same canonical markers used
    by HierAnnot. ``hierarchy`` may be a list of :class:`MarkerProgram` objects,
    dictionaries, or a built-in hierarchy name such as ``"breast_tme"``.
    ``program_sets`` may be one or more built-in flat marker-program set names
    such as ``"tumor_breast"``, ``"cell_cycle"`` or ``"stress"``. ``programs``
    may contain custom :class:`MalignantProgram` objects or dictionaries.

    Parameters
    ----------
    hierarchy
        Hierarchical marker program(s), dictionaries, or built-in hierarchy
        names from which to collect ``positive_markers``.
    program_sets
        Built-in flat marker-program set name or names. These are resolved with
        :func:`hierannot.builtin.get_builtin_marker_program_set`.
    programs
        Additional custom flat marker programs.
    uppercase
        If True, return marker names in upper case.
    deduplicate
        If True, remove duplicate genes while preserving first occurrence.
    return_format
        ``"list"`` returns a list of marker genes. ``"dataframe"`` returns one
        row per marker with source/program metadata.
    """
    rows = []

    for node, source in _iter_hierarchy_nodes(hierarchy, source="hierarchy"):
        for marker in list(getattr(node, "positive_markers", []) or []):
            gene = str(marker).strip()
            if not gene:
                continue
            rows.append({"source": source, "program": str(node.name), "program_type": "hierarchy", "marker": gene.upper() if uppercase else gene})

    for prog, source in _iter_flat_programs(program_sets, source="program_set"):
        for marker in list(getattr(prog, "positive_markers", []) or []):
            gene = str(marker).strip()
            if not gene:
                continue
            rows.append({"source": source, "program": str(prog.name), "program_type": "flat_program", "marker": gene.upper() if uppercase else gene})

    for prog, source in _iter_flat_programs(programs, source="program"):
        for marker in list(getattr(prog, "positive_markers", []) or []):
            gene = str(marker).strip()
            if not gene:
                continue
            rows.append({"source": source, "program": str(prog.name), "program_type": "flat_program", "marker": gene.upper() if uppercase else gene})

    fmt = str(return_format).strip().lower()
    if fmt in {"dataframe", "df", "table"}:
        df = pd.DataFrame(rows, columns=["source", "program", "program_type", "marker"])
        if deduplicate and not df.empty:
            df = df.drop_duplicates(subset=["marker"], keep="first").reset_index(drop=True)
        return df
    if fmt not in {"list", "genes", "markers"}:
        raise ValueError("return_format must be 'list' or 'dataframe'")
    markers = [r["marker"] for r in rows]
    return _dedupe_preserve(markers) if deduplicate else markers



def _coerce_flat_program_list(programs):
    """Return flat marker programs as MalignantProgram objects."""
    return [prog for prog, _source in _iter_flat_programs(programs, source="program")]


def _program_role(prog: MalignantProgram) -> str:
    return str((getattr(prog, "metadata", {}) or {}).get("reporting_role", "unknown"))


def _program_group(prog: MalignantProgram) -> str:
    meta = getattr(prog, "metadata", {}) or {}
    return str(meta.get("competition_group") or prog.name)


def _overlap_coefficient(a: set[str], b: set[str]) -> float:
    denom = min(len(a), len(b))
    if denom <= 0:
        return 0.0
    return len(a & b) / float(denom)


def _jaccard(a: set[str], b: set[str]) -> float:
    denom = len(a | b)
    if denom <= 0:
        return 0.0
    return len(a & b) / float(denom)


def _available_marker_filter(markers, gene_set: set[str] | None, *, uppercase: bool = True) -> list[str]:
    out = []
    seen = set()
    for marker in markers or []:
        gene = str(marker).strip()
        if not gene:
            continue
        key = gene.upper() if uppercase else gene
        if gene_set is not None and key not in gene_set:
            continue
        if key not in seen:
            out.append(key if uppercase else gene)
            seen.add(key)
    return out


def validate_marker_programs(
    programs,
    *,
    genes: Iterable[str] | None = None,
    hierarchy=None,
    max_overlap_coefficient: float = 0.75,
    min_markers: int = 3,
    min_unique_markers: int = 2,
    uppercase: bool = True,
):
    """Validate flat tumor/auxiliary marker programs for coverage and overlap.

    This utility is intended for panel design and sanity checks before running
    the flat tumor/auxiliary program track. It does not change any scoring
    behavior. Use :func:`compile_marker_programs_for_panel` when you want to
    filter programs to an input gene panel and optionally merge highly
    overlapping programs.

    Parameters
    ----------
    programs
        Custom :class:`MalignantProgram` objects, dictionaries, or built-in
        marker-program set names such as ``"tumor_colon"``.
    genes
        Optional input gene universe. When provided, coverage metrics are
        computed after filtering to these genes.
    hierarchy
        Optional hierarchy or built-in hierarchy name. When provided, the report
        includes overlap with hierarchy positive markers to help identify tumor
        status programs that mostly recapitulate normal lineage markers.
    max_overlap_coefficient
        Threshold used to flag high-overlap programs within the same reporting
        role and competition group.
    min_markers
        Minimum available positive markers recommended for a program.
    min_unique_markers
        Minimum markers unique within the same role/group recommended for a
        separable competing program.
    uppercase
        Normalize genes/markers to upper case for matching.

    Returns
    -------
    pandas.DataFrame
        One row per program with coverage, uniqueness, overlap, and warning
        fields.
    """
    progs = _coerce_flat_program_list(programs)
    gene_set = None if genes is None else {str(g).strip().upper() if uppercase else str(g).strip() for g in genes if str(g).strip()}
    hierarchy_markers = set()
    if hierarchy is not None:
        hmarkers = collect_positive_marker_genes(hierarchy=hierarchy, uppercase=uppercase, deduplicate=True)
        hierarchy_markers = set(hmarkers)

    pos_sets = {}
    rows = []
    for prog in progs:
        orig_pos = _available_marker_filter(prog.positive_markers, None, uppercase=uppercase)
        avail_pos = _available_marker_filter(prog.positive_markers, gene_set, uppercase=uppercase)
        avail_neg = _available_marker_filter(prog.negative_markers, gene_set, uppercase=uppercase)
        role = _program_role(prog)
        group = _program_group(prog)
        pos_sets[prog.name] = set(avail_pos)
        rows.append({
            "program_name": prog.name,
            "reporting_role": role,
            "competition_group": group,
            "n_positive_markers_original": len(orig_pos),
            "n_positive_markers_available": len(avail_pos),
            "positive_marker_coverage_fraction": (len(avail_pos) / len(orig_pos)) if orig_pos else 0.0,
            "n_negative_markers_available": len(avail_neg),
            "available_positive_markers": ";".join(avail_pos),
            "available_negative_markers": ";".join(avail_neg),
            "n_overlap_with_hierarchy_markers": len(set(avail_pos) & hierarchy_markers),
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    unique_counts = []
    high_overlap_with = []
    max_overlap = []
    max_jaccard = []
    for _, row in df.iterrows():
        name = row["program_name"]
        role = row["reporting_role"]
        group = row["competition_group"]
        cur = pos_sets.get(name, set())
        peers = [r["program_name"] for _, r in df.iterrows() if r["program_name"] != name and r["reporting_role"] == role and r["competition_group"] == group]
        peer_union = set().union(*(pos_sets.get(p, set()) for p in peers)) if peers else set()
        unique_counts.append(len(cur - peer_union))
        hits = []
        ovs = []
        jacs = []
        for peer in peers:
            o = _overlap_coefficient(cur, pos_sets.get(peer, set()))
            j = _jaccard(cur, pos_sets.get(peer, set()))
            ovs.append(o)
            jacs.append(j)
            if o >= max_overlap_coefficient and cur and pos_sets.get(peer, set()):
                hits.append(peer)
        high_overlap_with.append(";".join(hits))
        max_overlap.append(max(ovs) if ovs else 0.0)
        max_jaccard.append(max(jacs) if jacs else 0.0)

    df["n_unique_positive_markers_within_role_group"] = unique_counts
    df["max_overlap_coefficient_within_role_group"] = max_overlap
    df["max_jaccard_within_role_group"] = max_jaccard
    df["high_overlap_with"] = high_overlap_with

    warnings = []
    for _, row in df.iterrows():
        vals = []
        if row["n_positive_markers_available"] < min_markers:
            vals.append("few_available_markers")
        if row["n_unique_positive_markers_within_role_group"] < min_unique_markers:
            vals.append("few_unique_markers_within_role_group")
        if str(row["high_overlap_with"]):
            vals.append("high_overlap_within_role_group")
        if row["n_positive_markers_available"] and row["n_overlap_with_hierarchy_markers"] / row["n_positive_markers_available"] >= 0.75:
            vals.append("mostly_hierarchy_markers")
        warnings.append(";".join(vals))
    df["warnings"] = warnings
    return df


def _merge_program_group(programs: list[MalignantProgram]) -> MalignantProgram:
    if len(programs) == 1:
        return programs[0]
    first = programs[0]
    pos = _dedupe_preserve([g for p in programs for g in p.positive_markers])
    neg = _dedupe_preserve([g for p in programs for g in p.negative_markers])
    meta = dict(getattr(first, "metadata", {}) or {})
    meta["merged_from"] = [p.name for p in programs]
    meta["merged_program_count"] = len(programs)
    label = meta.get("reporting_label") or first.name
    name = str(meta.get("competition_group") or label or first.name)
    return MalignantProgram(
        name=name,
        positive_markers=pos,
        negative_markers=neg,
        description="Merged marker program compiled from: " + ", ".join(p.name for p in programs),
        tags=_dedupe_preserve([t for p in programs for t in getattr(p, "tags", [])]),
        metadata=meta,
    )


def compile_marker_programs_for_panel(
    programs,
    genes: Iterable[str],
    *,
    min_markers: int = 3,
    min_unique_markers: int = 2,
    max_overlap_coefficient: float = 0.75,
    merge_high_overlap: bool = False,
    drop_sparse: bool = True,
    uppercase: bool = True,
    return_report: bool = True,
):
    """Compile flat marker programs to the genes available in a panel.

    The compiler filters positive and negative markers to ``genes`` and can
    optionally merge highly overlapping programs within the same reporting role
    and competition group. Automatic merging is intentionally opt-in so expert
    users can inspect the validation report before changing program definitions.

    Parameters
    ----------
    programs
        Custom programs, dictionaries, or built-in marker-program set names.
    genes
        Input gene universe/panel.
    min_markers
        Drop programs with fewer available positive markers when ``drop_sparse``
        is True.
    min_unique_markers
        Used in the validation report to flag poorly separable programs within
        the same role/group.
    max_overlap_coefficient
        Threshold for high-overlap warnings and optional merging.
    merge_high_overlap
        If True, merge high-overlap programs only when they share the same
        ``reporting_role`` and ``competition_group``.
    drop_sparse
        If True, remove programs with too few available positive markers.
    uppercase
        Normalize returned marker symbols to upper case.
    return_report
        If True, return ``(compiled_programs, report)``. Otherwise return only
        the compiled program list.
    """
    gene_set = {str(g).strip().upper() if uppercase else str(g).strip() for g in genes if str(g).strip()}
    source_programs = _coerce_flat_program_list(programs)
    compiled = []
    for prog in source_programs:
        pos = _available_marker_filter(prog.positive_markers, gene_set, uppercase=uppercase)
        neg = _available_marker_filter(prog.negative_markers, gene_set, uppercase=uppercase)
        if drop_sparse and len(pos) < min_markers:
            continue
        compiled.append(MalignantProgram(
            name=prog.name,
            positive_markers=pos,
            negative_markers=neg,
            description=prog.description,
            tags=list(getattr(prog, "tags", [])),
            metadata=dict(getattr(prog, "metadata", {}) or {}),
        ))

    if merge_high_overlap and len(compiled) > 1:
        # Connected components among same-role/group high-overlap programs.
        n = len(compiled)
        parent = list(range(n))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra
        sets = [set(p.positive_markers) for p in compiled]
        for i in range(n):
            for j in range(i + 1, n):
                if _program_role(compiled[i]) != _program_role(compiled[j]):
                    continue
                if _program_group(compiled[i]) != _program_group(compiled[j]):
                    continue
                if _overlap_coefficient(sets[i], sets[j]) >= max_overlap_coefficient:
                    union(i, j)
        groups = {}
        for i, prog in enumerate(compiled):
            groups.setdefault(find(i), []).append(prog)
        compiled = [_merge_program_group(group) for group in groups.values()]

    report = validate_marker_programs(
        compiled,
        genes=genes,
        max_overlap_coefficient=max_overlap_coefficient,
        min_markers=min_markers,
        min_unique_markers=min_unique_markers,
        uppercase=uppercase,
    )
    if return_report:
        return compiled, report
    return compiled


def compile_builtin_marker_program_set_for_panel(
    program_set: str,
    genes: Iterable[str],
    *,
    include_general: bool | None = None,
    **kwargs,
):
    """Compile a built-in flat marker-program set to an input gene panel.

    This is a convenience wrapper around
    :func:`compile_marker_programs_for_panel` that first resolves
    ``program_set`` with :func:`hierannot.builtin.get_builtin_marker_program_set`.
    """
    from .builtin import get_builtin_marker_program_set
    programs = get_builtin_marker_program_set(program_set, include_general=include_general)
    return compile_marker_programs_for_panel(programs, genes, **kwargs)
