from __future__ import annotations

from difflib import SequenceMatcher
from typing import Dict, List, Tuple

import pandas as pd
import re


_SPECIFIC_TISSUE_SYNONYMS: Dict[str, List[str]] = {
    "kidney": ["renal", "rcc", "ccrcc"],
    "liver": ["hepatic", "hcc", "cholangiocarcinoma"],
    "lung": ["pulmonary", "nsclc", "sclc"],
    "brain": ["cerebral", "cns", "glioma", "gbm", "neural", "hippocampus", "cortex"],
    "colon": ["colorectal", "crc", "gastric", "stomach", "gut"],
    "skin": ["dermal", "epidermal", "cutaneous", "melanoma"],
    "breast": ["mammary", "tnbc", "her2", "er_positive", "pr_positive"],
    "pancreas": ["pancreatic", "pdac"],
    "tonsil": ["tonsillar", "oropharynx", "oropharyngeal", "headneck", "hnscc", "tonsillar_squamous", "lymphoid_tonsil"],
}

_BROAD_CONTEXT_TOKENS = {
    "tumor", "tumour", "cancer", "solid", "tme", "microenvironment", "tumor_microenvironment"
}
_IMMUNE_CONTEXT_TOKENS = {"immune", "pbmc", "blood", "hematopoietic", "lymphoid"}


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
    from .registry import BUILTIN_HIERARCHIES

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


_SHORT_IMMUNE_TOKENS = {"dc", "nk", "treg"}
_IMMUNE_EXACT_PHRASES = {
    "immune", "lymphoid", "myeloid", "macrophage", "monocyte", "neutrophil",
    "mast cell", "plasma cell", "b cell", "t cell", "nk cell", "dendritic cell",
    "microglia", "plasmablast", "granulocyte", "lymphocyte",
}


def _normalize_lineage_text(text: str) -> str:
    q = str(text).strip().lower()
    q = re.sub(r"[_\-/]+", " ", q)
    q = re.sub(r"[^a-z0-9 >]+", " ", q)
    q = re.sub(r"\s+", " ", q).strip()
    return q


def _tokenize_lineage_text(text: str) -> List[str]:
    norm = _normalize_lineage_text(text)
    if not norm:
        return []
    return [tok for tok in norm.split(" ") if tok]


def is_immune_like_lineage(*texts: str) -> bool:
    """Return True when any supplied lineage text looks immune-like.

    Strategy:
    - check case-senstive exact text for ``pDC``, ``mDC``, ``MNP``
    - normalize into canonical phrases with word boundaries preserved
    - check exact multiword phrases such as ``b cell`` and ``dendritic cell``
    - check only a restricted set of short exact tokens such as ``dc``/``nk``
    """
    for text in texts:
        if text in {"pDC", "mDC", "MNP"}:
            return True
        
        norm = _normalize_lineage_text(text)
        if not norm:
            continue
        if any(f" {phrase} " in f" {norm} " for phrase in _IMMUNE_EXACT_PHRASES):
            return True
        
        tokens = set(_tokenize_lineage_text(text))
        if any(tok in tokens for tok in _SHORT_IMMUNE_TOKENS):
            return True

    return False
