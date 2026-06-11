from __future__ import annotations

import re
from typing import Any, Dict, Optional


NEW_STYLE_RE = re.compile(
    r"^(?P<base>.+?)_(?P<protomer>p\d+)_(?P<tautomer>t\d+)_(?P<conformer>c\d+)$",
    re.IGNORECASE,
)
LEGACY_POSE_RE = re.compile(r"^(?P<base>.+?)(?:__|_)(?P<pose>pose\d+)$", re.IGNORECASE)


def _to_index(tag: str) -> Optional[int]:
    if not tag:
        return None
    digits = re.search(r"(\d+)$", tag)
    return int(digits.group(1)) if digits else None


def parse_ligand_variant(name: str) -> Dict[str, Any]:
    original = (name or "").strip()
    parsed: Dict[str, Any] = {
        "LigandVariant": original,
        "LigandBase": original,
        "ProtomerTag": "",
        "TautomerTag": "",
        "ConformerTag": "",
        "LegacyPoseTag": "",
        "StateTag": "",
        "ConformerIndex": None,
        "ProtomerIndex": None,
        "TautomerIndex": None,
    }
    if not original:
        return parsed

    match = NEW_STYLE_RE.match(original)
    if match:
        protomer = match.group("protomer")
        tautomer = match.group("tautomer")
        conformer = match.group("conformer")
        parsed.update(
            {
                "LigandBase": match.group("base"),
                "ProtomerTag": protomer,
                "TautomerTag": tautomer,
                "ConformerTag": conformer,
                "StateTag": f"{protomer}_{tautomer}",
                "ProtomerIndex": _to_index(protomer),
                "TautomerIndex": _to_index(tautomer),
                "ConformerIndex": _to_index(conformer),
            }
        )
        return parsed

    match = LEGACY_POSE_RE.match(original)
    if match:
        pose_tag = match.group("pose")
        parsed.update(
            {
                "LigandBase": match.group("base"),
                "LegacyPoseTag": pose_tag,
            }
        )
        return parsed

    return parsed


def ligand_base(name: str) -> str:
    return str(parse_ligand_variant(name).get("LigandBase", name or ""))
