from __future__ import annotations

from typing import List

from ..datamodels import MarkerProgram
from .schema import _meta, _mp


def _immune_core() -> List[MarkerProgram]:
    cd4 = _mp("CD4 T cell", ["IL7R", "LTB", "MAL", "HLA-DRA", "TRAC"], ["NKG7", "PRF1"], metadata=_meta("immune", "pan-tissue", "immune_core", role="major", lineage_module="immune"))
    cd8 = _mp("CD8 T cell", ["NKG7", "CCL5", "PRF1", "GZMB", "TRAC"], ["IL7R"], metadata=_meta("immune", "pan-tissue", "immune_core", role="major", lineage_module="immune"))
    treg = _mp("Treg", ["IL2RA", "FOXP3", "TIGIT", "CTLA4", "LTB"], ["NKG7", "PRF1"], metadata=_meta("immune", "pan-tissue", "immune_core", role="optional", lineage_module="immune"))
    tcell = _mp("T cell", ["CD3D", "CD3E", "TRBC1", "TRAC", "LTB"], ["MS4A1", "LYZ"], [cd4, cd8, treg], "Pan-T lineage markers.", metadata=_meta("immune", "pan-tissue", "immune_core", role="major", lineage_module="immune"))

    plasma = _mp("Plasma cell", ["MZB1", "JCHAIN", "SDC1", "XBP1", "IGKC"], ["CD3D", "LYZ"], metadata=_meta("immune", "pan-tissue", "immune_core", role="major", lineage_module="immune"))
    bcell = _mp("B cell", ["MS4A1", "CD79A", "CD74", "HLA-DRA", "CD79B"], ["CD3D", "LYZ"], [plasma], metadata=_meta("immune", "pan-tissue", "immune_core", role="major", lineage_module="immune"))

    nk = _mp("NK cell", ["NKG7", "KLRD1", "GNLY", "PRF1", "FCGR3A"], ["CD3D", "MS4A1"], metadata=_meta("immune", "pan-tissue", "immune_core", role="major", lineage_module="immune"))
    macro = _mp("Macrophage", ["LYZ", "C1QA", "C1QB", "APOE", "FCER1G"], ["CD3D", "EPCAM"], metadata=_meta("immune", "pan-tissue", "immune_core", role="major", lineage_module="immune"))
    mono = _mp("Monocyte", ["LYZ", "S100A8", "S100A9", "CTSS", "FCN1"], ["CD3D", "EPCAM"], metadata=_meta("immune", "pan-tissue", "immune_core", role="major", lineage_module="immune"))
    dc = _mp("Dendritic cell", ["FCER1A", "CST3", "HLA-DRA", "CD74", "CLEC10A"], ["NKG7", "EPCAM"], metadata=_meta("immune", "pan-tissue", "immune_core", role="major", lineage_module="immune"))
    myeloid = _mp("Myeloid", ["LYZ", "TYMP", "FCER1G", "CTSS", "SAT1"], ["CD3D", "EPCAM"], [macro, mono, dc], metadata=_meta("immune", "pan-tissue", "immune_core", role="major", lineage_module="immune"))

    mast = _mp("Mast cell", ["TPSAB1", "TPSB2", "KIT", "CPA3", "HDC"], ["CD3D", "MS4A1"], metadata=_meta("immune", "pan-tissue", "immune_core", role="major", lineage_module="immune"))
    neut = _mp("Neutrophil", ["S100A8", "S100A9", "FCGR3B", "CXCR2", "CSF3R"], ["CD3D", "MS4A1"], metadata=_meta("immune", "pan-tissue", "immune_core", role="major", lineage_module="immune"))
    return [tcell, bcell, nk, myeloid, mast, neut]



def _parenchymal_fallback_core() -> List[MarkerProgram]:
    basal_ep = _mp(
        "Basal epithelial",
        ["KRT5", "KRT14", "KRT17", "TP63", "DST"],
        ["EPCAM"],
        metadata=_meta("epithelial", "pan-solid-tissue", "solid_tissue_core", role="optional", lineage_module="parenchymal"),
    )
    luminal_ep = _mp(
        "Luminal epithelial",
        ["EPCAM", "KRT8", "KRT18", "KRT19", "MSLN"],
        ["KRT5"],
        metadata=_meta("epithelial", "pan-solid-tissue", "solid_tissue_core", role="optional", lineage_module="parenchymal"),
    )
    epithelial = _mp(
        "Epithelial",
        ["EPCAM", "KRT8", "KRT18", "KRT19", "MSLN"],
        ["COL1A1", "PTPRC", "PECAM1", "VWF"],
        children=[luminal_ep, basal_ep],
        metadata=_meta("solid_tissue", "pan-solid-tissue", "solid_tissue_core", role="major", lineage_module="parenchymal"),
    )
    return [epithelial]


def _stromal_core() -> List[MarkerProgram]:
    fibro = _mp(
        "Fibroblast",
        ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"],
        ["EPCAM", "PTPRC", "PECAM1", "KDR"],
        children=[],
        metadata=_meta("solid_tissue", "pan-solid-tissue", "solid_tissue_core", role="major", lineage_module="stromal"),
    )
    return [fibro]


def _vascular_core() -> List[MarkerProgram]:
    cap_endo = _mp(
        "Capillary endothelial",
        ["PECAM1", "EMCN", "KDR", "RGCC", "CA4"],
        ["ACKR1", "GJA5"],
        metadata=_meta("endothelial", "pan-solid-tissue", "solid_tissue_core", role="optional", lineage_module="vascular"),
    )
    endothelial = _mp(
        "Endothelial",
        ["PECAM1", "VWF", "EMCN", "KDR", "CLDN5"],
        ["COL1A1", "EPCAM", "PTPRC"],
        children=[cap_endo],
        metadata=_meta("solid_tissue", "pan-solid-tissue", "solid_tissue_core", role="major", lineage_module="vascular"),
    )
    return [endothelial]


def _mural_core() -> List[MarkerProgram]:
    pericyte = _mp(
        "Pericyte",
        ["RGS5", "MCAM", "CSPG4", "PDGFRB", "DES"],
        ["MYH11", "CNN1", "EPCAM", "KRT8", "KRT18", "KRT19"],
        metadata=_meta("mural", "pan-solid-tissue", "solid_tissue_core", role="major", lineage_module="mural"),
    )
    mural = _mp(
        "Mural",
        ["RGS5", "MCAM", "CSPG4", "ACTA2", "MYH11"],
        ["EPCAM", "PTPRC"],
        children=[pericyte],
        metadata=_meta("solid_tissue", "pan-solid-tissue", "solid_tissue_core", role="major", lineage_module="mural"),
    )
    return [mural]


def _solid_tissue_core() -> List[MarkerProgram]:
    return _parenchymal_fallback_core() + _stromal_core() + _vascular_core() + _mural_core()

def _brain_core() -> List[MarkerProgram]:
    excit = _mp(
        "Excitatory neuron",
        ["SLC17A7", "CAMK2A", "SATB2", "NRGN", "SYT1"],
        ["GAD1", "AIF1"],
        metadata=_meta("brain", "brain", "brain_core", role="optional", lineage_module="brain_parenchymal"),
    )
    inhib = _mp(
        "Inhibitory neuron",
        ["GAD1", "GAD2", "SLC6A1", "DLX1", "DLX2"],
        ["SLC17A7", "AIF1"],
        metadata=_meta("brain", "brain", "brain_core", role="optional", lineage_module="brain_parenchymal"),
    )
    neuron = _mp(
        "Neuron",
        ["RBFOX3", "SNAP25", "SYT1", "TUBB3", "STMN2"],
        ["AQP4", "AIF1"],
        [excit, inhib],
        metadata=_meta("brain", "brain", "brain_core", role="major", lineage_module="brain_parenchymal"),
    )
    astro = _mp(
        "Astrocyte",
        ["AQP4", "GFAP", "SLC1A3", "ALDH1L1", "GJA1"],
        ["RBFOX3", "MOG"],
        metadata=_meta("brain", "brain", "brain_core", role="major", lineage_module="brain_parenchymal"),
    )
    olig = _mp(
        "Oligodendrocyte",
        ["MOG", "MBP", "PLP1", "MOBP", "MAG"],
        ["AIF1", "AQP4"],
        metadata=_meta("brain", "brain", "brain_core", role="major", lineage_module="brain_parenchymal"),
    )
    opc = _mp(
        "OPC",
        ["PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10"],
        ["MOG", "AIF1"],
        metadata=_meta("brain", "brain", "brain_core", role="major", lineage_module="brain_parenchymal"),
    )
    micro = _mp(
        "Microglia",
        ["AIF1", "P2RY12", "TMEM119", "CX3CR1", "TYROBP"],
        ["RBFOX3", "MOG"],
        metadata=_meta("brain", "brain", "brain_core", role="major", lineage_module="brain_immune"),
    )
    endot = _mp(
        "Endothelial",
        ["PECAM1", "CLDN5", "KDR", "EMCN", "VWF"],
        ["AQP4", "RBFOX3"],
        metadata=_meta("brain", "brain", "brain_core", role="major", lineage_module="vascular"),
    )
    return [neuron, astro, olig, opc, micro, endot]


def _tme_core() -> List[MarkerProgram]:
    roots = _solid_tissue_core()
    roots.append(_mp("Immune", ["PTPRC", "TYROBP", "LST1", "HLA-DRA", "CD53"], ["EPCAM", "COL1A1"], _immune_core(), "Top-level immune branch for mixed tissue microenvironments.", metadata=_meta("tme", "pan-solid-tissue", "tme_core")))
    return roots


def _colon_tme() -> List[MarkerProgram]:
    roots = _tme_core()
    for i, node in enumerate(roots):
        if node.name == "Epithelial":
            roots[i] = _mp(
                "Intestinal epithelial",
                ["KRT20", "CEACAM1", "SLC26A3", "MUC2", "TFF3", "SOX9"],
                ["PTPRC", "COL1A1", "PECAM1", "VWF"],
                children=[
                    _mp("Absorptive-like", ["KRT20", "CA1", "CEACAM1", "SLC26A3", "FABP1"], ["MUC2"], metadata=_meta("epithelial", "colon", "colon_tme", role="major", lineage_module="colon_parenchymal")),
                    _mp("Goblet-like", ["MUC2", "TFF3", "SPINK4", "CLCA1", "AGR2"], ["CA1"], metadata=_meta("epithelial", "colon", "colon_tme", role="major", lineage_module="colon_parenchymal")),
                    _mp("Stem/TA-like", ["LGR5", "OLFM4", "MKI67", "SOX9", "ASCL2"], ["KRT20"], metadata=_meta("epithelial", "colon", "colon_tme", role="optional", lineage_module="colon_parenchymal")),
                ],
                metadata=_meta("solid_tissue", "colon", "colon_tme", role="major", lineage_module="colon_parenchymal"),
            )
    return roots

def _kidney_tme() -> List[MarkerProgram]:
    roots = _tme_core()
    for i, node in enumerate(roots):
        if node.name == "Epithelial":
            roots[i] = _mp(
                "Renal parenchymal",
                ["LRP2", "CUBN", "SLC34A1", "UMOD", "KRT19", "AQP2", "NPHS1"],
                ["PTPRC", "COL1A1", "PECAM1", "VWF"],
                children=[
                    _mp("Proximal tubule", ["LRP2", "CUBN", "SLC34A1", "ALDOB", "AQP1"], ["KRT19"], metadata=_meta("epithelial", "kidney", "kidney_tme", role="major", lineage_module="kidney_parenchymal")),
                    _mp("Distal/TAL-like", ["SLC12A1", "UMOD", "KCNJ1", "FXYD2", "CLDN16"], ["LRP2"], metadata=_meta("epithelial", "kidney", "kidney_tme", role="major", lineage_module="kidney_parenchymal")),
                    _mp("Collecting duct", ["KRT19", "EPCAM", "AQP2", "FXYD4", "SCNN1G"], ["LRP2"], metadata=_meta("epithelial", "kidney", "kidney_tme", role="major", lineage_module="kidney_parenchymal")),
                    _mp("Podocyte", ["NPHS1", "NPHS2", "PODXL", "WT1", "SYNPO"], ["KRT19"], metadata=_meta("epithelial", "kidney", "kidney_tme", role="major", lineage_module="kidney_parenchymal")),
                ],
                metadata=_meta("solid_tissue", "kidney", "kidney_tme", role="major", lineage_module="kidney_parenchymal"),
            )
        elif node.name == "Mural":
            node.children = list(node.children) + [
                _mp("Mesangial/pericyte", ["RGS5", "MCAM", "CSPG4", "PDGFRB", "DES"], ["EPCAM", "PTPRC"], metadata=_meta("mural", "kidney", "kidney_tme", role="optional", lineage_module="kidney_stromal")),
            ]
    return roots

def _tonsil_tme() -> List[MarkerProgram]:
    return _tme_core()


def _breast_tme() -> List[MarkerProgram]:
    roots = _tme_core()
    for i, node in enumerate(roots):
        if node.name == "Epithelial":
            roots[i] = _mp(
                "Mammary epithelial",
                ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14"],
                ["PTPRC", "COL1A1", "PECAM1", "VWF"],
                children=[
                    _mp("Luminal epithelial", ["EPCAM", "KRT8", "KRT18", "ESR1", "PGR"], ["KRT5", "KRT14"], metadata=_meta("epithelial", "breast", "breast_tme", role="major", lineage_module="breast_parenchymal")),
                    _mp("Basal epithelial", ["KRT5", "KRT14", "KRT17", "TP63", "KRT15"], ["ESR1", "PGR"], metadata=_meta("epithelial", "breast", "breast_tme", role="major", lineage_module="breast_parenchymal")),
                    _mp("Secretory/alveolar epithelial", ["EPCAM", "KRT8", "KRT18", "ELF5", "CSN2"], ["KRT5"], metadata=_meta("epithelial", "breast", "breast_tme", role="optional", lineage_module="breast_parenchymal")),
                ],
                metadata=_meta("solid_tissue", "breast", "breast_tme", role="major", lineage_module="breast_parenchymal"),
            )
    return roots

def _pancreas_tme() -> List[MarkerProgram]:
    roots = _tme_core()
    for i, node in enumerate(roots):
        if node.name == "Epithelial":
            roots[i] = _mp(
                "Pancreatic parenchymal",
                ["PRSS1", "CPA1", "CPB1", "REG1A", "KRT19", "EPCAM", "CHGA", "CHGB"],
                ["PTPRC", "TYROBP", "PECAM1", "VWF", "COL1A1", "COL1A2"],
                children=[
                    _mp(
                        "Ductal",
                        ["KRT19", "EPCAM", "MSLN", "KRT8", "KRT18"],
                        ["CPA1", "CPB1", "PRSS1"],
                        metadata=_meta("epithelial", "pancreas", "pancreas_tme", role="major", lineage_module="pancreas_parenchymal"),
                    ),
                    _mp(
                        "Acinar",
                        ["PRSS1", "PRSS3", "CPA1", "CPB1", "REG1A"],
                        ["KRT19", "EPCAM", "PECAM1", "VWF", "PTPRC"],
                        metadata=_meta("epithelial", "pancreas", "pancreas_tme", role="major", lineage_module="pancreas_parenchymal"),
                    ),
                    _mp(
                        "Endocrine",
                        ["CHGA", "CHGB", "ISL1", "PAX6", "PCSK1"],
                        ["KRT19", "PRSS1", "CPA1"],
                        metadata=_meta("epithelial", "pancreas", "pancreas_tme", role="major", lineage_module="pancreas_parenchymal"),
                    ),
                ],
                metadata=_meta("solid_tissue", "pancreas", "pancreas_tme", role="major", lineage_module="pancreas_parenchymal"),
            )
        elif node.name == "Endothelial":
            node.negative_markers = ["COL1A1", "EPCAM", "PTPRC"]
        elif node.name == "Immune":
            node.negative_markers = ["EPCAM", "COL1A1", "PECAM1", "VWF"]
            for child in node.children:
                if child.name == "NK cell":
                    child.negative_markers = list(dict.fromkeys(child.negative_markers + ["PRSS1", "CPA1", "CPB1", "PECAM1", "VWF"]))
                elif child.name in {"Myeloid", "T cell", "B cell"}:
                    child.negative_markers = list(dict.fromkeys(child.negative_markers + ["PECAM1", "VWF"]))
    return roots

def _liver_tme() -> List[MarkerProgram]:
    roots = _tme_core()
    for i, node in enumerate(roots):
        if node.name == "Epithelial":
            roots[i] = _mp(
                "Hepatic parenchymal",
                ["ALB", "APOA1", "TTR", "KRT19", "SOX9"],
                ["PTPRC", "COL1A1", "PECAM1", "VWF"],
                children=[
                    _mp("Hepatocyte-like", ["ALB", "APOA1", "TTR", "HP", "CPS1"], ["KRT19"], metadata=_meta("epithelial", "liver", "liver_tme", role="major", lineage_module="liver_parenchymal")),
                    _mp("Cholangiocyte-like", ["KRT19", "EPCAM", "KRT8", "KRT18", "SOX9"], ["ALB"], metadata=_meta("epithelial", "liver", "liver_tme", role="major", lineage_module="liver_parenchymal")),
                ],
                metadata=_meta("solid_tissue", "liver", "liver_tme", role="major", lineage_module="liver_parenchymal"),
            )
    return roots

def _lung_tme() -> List[MarkerProgram]:
    roots = _tme_core()
    for i, node in enumerate(roots):
        if node.name == "Epithelial":
            roots[i] = _mp(
                "Pulmonary parenchymal",
                ["SFTPA1", "SFTPB", "SFTPC", "SCGB1A1", "KRT19", "KRT5"],
                ["PTPRC", "COL1A1", "PECAM1", "VWF"],
                children=[
                    _mp("Alveolar epithelial", ["SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "NAPSA"], ["KRT5"], metadata=_meta("epithelial", "lung", "lung_tme", role="major", lineage_module="lung_parenchymal")),
                    _mp("Airway epithelial", ["SCGB1A1", "KRT19", "KRT8", "FOXJ1", "PIFO"], ["SFTPC"], metadata=_meta("epithelial", "lung", "lung_tme", role="major", lineage_module="lung_parenchymal")),
                    _mp("Basal epithelial", ["KRT5", "KRT14", "KRT17", "TP63", "KRT15"], ["SFTPC"], metadata=_meta("epithelial", "lung", "lung_tme", role="major", lineage_module="lung_parenchymal")),
                ],
                metadata=_meta("solid_tissue", "lung", "lung_tme", role="major", lineage_module="lung_parenchymal"),
            )
    return roots

def _skin_tme() -> List[MarkerProgram]:
    roots = _tme_core()
    for i, node in enumerate(roots):
        if node.name == "Epithelial":
            roots[i] = _mp(
                "Cutaneous epithelial",
                ["KRT5", "KRT14", "KRT15", "KRT1", "KRT10", "IVL"],
                ["PTPRC", "COL1A1", "PECAM1", "VWF"],
                children=[
                    _mp("Basal keratinocyte", ["KRT5", "KRT14", "KRT15", "TP63", "DST"], ["KRT1", "KRT10"], metadata=_meta("epithelial", "skin", "skin_tme", role="major", lineage_module="skin_parenchymal")),
                    _mp("Differentiated keratinocyte", ["KRT1", "KRT10", "KRTDAP", "IVL", "SPRR1B"], ["KRT14"], metadata=_meta("epithelial", "skin", "skin_tme", role="major", lineage_module="skin_parenchymal")),
                ],
                metadata=_meta("solid_tissue", "skin", "skin_tme", role="major", lineage_module="skin_parenchymal"),
            )
    return roots

def _normalize_query_tokens(text: str) -> List[str]:
    q = str(text).strip().lower()
    q = re.sub(r"[^a-z0-9]+", " ", q)
    q = re.sub(r"\s+", " ", q).strip()
    return [t for t in q.split(" ") if t]


def _canonical_tissue_to_pack(tissue: str) -> str:
    mapping = {
        "breast": "breast_tme",
        "pancreas": "pancreas_tme",
        "kidney": "kidney_tme",
        "colon": "colon_tme",
        "tonsil": "tonsil_tme",
        "brain": "brain_core",
        "liver": "liver_tme",
        "lung": "lung_tme",
        "skin": "skin_tme",
    }
    return mapping.get(tissue, "tme_core")


def _query_specific_tissues(query: str) -> Tuple[List[str], List[str], List[str]]:
    tokens = _normalize_query_tokens(query)
    specific = []
    broad = []
    unknown = []
    for tok in tokens:
        matched = False
        for canon, syns in _SPECIFIC_TISSUE_SYNONYMS.items():
            if tok == canon or tok in syns:
                specific.append(canon)
                matched = True
                break
        if matched:
            continue
        if tok in _BROAD_CONTEXT_TOKENS or tok in _IMMUNE_CONTEXT_TOKENS:
            broad.append(tok)
        else:
            unknown.append(tok)
    return sorted(set(specific)), broad, unknown


def _fuzzy_specific_tissues(tokens: List[str], cutoff: float = 0.75) -> List[str]:
    matches = []
    for tok in tokens:
        if len(tok) < 4:
            continue
        best_tissue = None
        best_score = 0.0
        for canon, syns in _SPECIFIC_TISSUE_SYNONYMS.items():
            for cand in [canon] + list(syns):
                score = SequenceMatcher(None, tok, cand).ratio()
                if score > best_score:
                    best_tissue = canon
                    best_score = score
        if best_tissue is not None and best_score >= cutoff:
            matches.append(best_tissue)
    return sorted(set(matches))


def suggest_builtin_hierarchies(query: str, top_k: int = 5, fuzzy_cutoff: float = 0.75) -> pd.DataFrame:
    """
    Rank built-in hierarchies for a free-text tissue description.

    Matching strategy:
    - exact tissue-token match
    - synonym-expanded tissue match
    - conservative fuzzy tissue match
    - broad context only acts as a weak fallback prior
    """
    specific, broad, unknown = _query_specific_tissues(query)
    fuzzy_specific = _fuzzy_specific_tissues(unknown, cutoff=fuzzy_cutoff) if not specific else []

    rows = []
    for name, spec in BUILTIN_HIERARCHIES.items():
        score = 0.0
        exact = 0.0
        fuzzy = 0.0
        context = 0.0

        pack_scope = str(spec.get("tissue_scope", "")).lower()

        for tissue in specific:
            target_pack = _canonical_tissue_to_pack(tissue)
            if name == target_pack or pack_scope == tissue:
                exact += 3.0

        for tissue in fuzzy_specific:
            target_pack = _canonical_tissue_to_pack(tissue)
            if name == target_pack or pack_scope == tissue:
                fuzzy += 1.0

        if not specific and not fuzzy_specific:
            if "brain" in broad and name == "brain_core":
                context += 1.5
            elif any(tok in _IMMUNE_CONTEXT_TOKENS for tok in broad) and name == "immune_core":
                context += 1.0
            elif any(tok in _BROAD_CONTEXT_TOKENS for tok in broad) and name == "tme_core":
                context += 1.0

        score = exact + fuzzy + context
        rows.append({
            "name": name,
            "score": score,
            "exact_score": exact,
            "fuzzy_score": fuzzy,
            "context_score": context,
            "category": spec["category"],
            "tissue_scope": spec["tissue_scope"],
            "description": spec["description"],
        })

    df = pd.DataFrame(rows).sort_values(
        ["score", "exact_score", "fuzzy_score", "context_score"],
        ascending=False,
    ).reset_index(drop=True)
    if top_k is not None:
        df = df.head(int(top_k))
    return df


def match_builtin_hierarchy(query: str, fallback: str = "tme_core", fuzzy_cutoff: float = 0.75) -> str:
    """
    Return the best-matching built-in hierarchy.

    If no specific tissue match is found, use broad-context fallback:
    - brain-like context -> brain_core
    - immune-only context -> immune_core
    - otherwise -> tme_core
    """
    specific, broad, unknown = _query_specific_tissues(query)
    if specific:
        return _canonical_tissue_to_pack(specific[0])

    fuzzy_specific = _fuzzy_specific_tissues(unknown, cutoff=fuzzy_cutoff)
    if fuzzy_specific:
        return _canonical_tissue_to_pack(fuzzy_specific[0])

    if any(tok in _IMMUNE_CONTEXT_TOKENS for tok in broad):
        return "immune_core"
    if any(tok in _BROAD_CONTEXT_TOKENS for tok in broad):
        return fallback
    return fallback
