"""Layer 1a: the enzyme universe and its biochemical annotations (UniProt + Rhea).

Serine-hydrolase calls matter because they define what PACE's fluorophosphonate
probes (TAMRA-FP, FP-biotin) can see. We call an enzyme a serine hydrolase when
UniProt annotates a catalytic Ser (active-site residue that is S in the sequence)
on a hydrolase, or when it belongs to a family whose members are FP-reactive.
"""
from __future__ import annotations

import re

import pandas as pd

from .data import RAW

# InterPro families that are FP-reactive serine hydrolase folds
SH_INTERPRO = {
    "IPR029058",  # alpha/beta hydrolase fold
    "IPR002641",  # patatin-like phospholipase (PNPLA)
    "IPR001087",  # GDSL lipase/esterase
    "IPR023631",  # amidase signature (FAAH)
    "IPR013818",  # lipase (pancreatic/lipoprotein lipase)
    "IPR002018",  # carboxylesterase type B (CES)
    "IPR001254",  # trypsin-like serine protease
    "IPR000383",  # xaa-pro dipeptidyl-peptidase-like
}

# Keyword-level lipid pathway tags we surface as columns
PATHWAY_PATTERNS = {
    "sphingolipid": r"sphingo|ceramide|sphingosine|sphinganine",
    "fatty_acid_elongation": r"elongat|3-oxoacyl|very-long-chain",
    "desaturation": r"desaturase|stearoyl-CoA",
    "phospholipid": r"phosphatidyl|glycerophospho|lysophospho|phospholipase",
    "neutral_lipid": r"triacylglycerol|diacylglycerol|monoacylglycerol|cholester(?:yl|ol) ester",
    "eicosanoid": r"prostaglandin|leukotriene|eicosa|arachidon|lipoxygenase",
    "cholesterol": r"cholesterol|sterol|squalene|lanosterol",
}


def _catalytic_ser(active_site: str, seq: str) -> bool:
    # Read each site's own annotation only, so one residue's note never leaks into the next
    for m in re.finditer(r"ACT_SITE (\d+)\b(.*?)(?=ACT_SITE \d|$)", active_site or "", re.S):
        i = int(m.group(1)) - 1
        if 0 <= i < len(seq) and seq[i] == "S":
            note = m.group(2).lower()
            if "nucleophile" in note or "charge relay" in note or "acyl-ester" in note or "note" not in note:
                return True
    return False


def build_catalog() -> pd.DataFrame:
    u = pd.read_csv(RAW / "uniprot" / "human_lipid_hydrolases.tsv", sep="\t").fillna("")
    u = u.rename(columns={
        "Entry": "uniprot", "Gene Names (primary)": "gene", "Protein names": "protein",
        "EC number": "ec", "Keywords": "keywords", "Catalytic activity": "catalytic_activity",
        "ChEMBL": "chembl", "Active site": "active_site", "Sequence": "sequence",
        "InterPro": "interpro", "Subcellular location [CC]": "location",
    })
    u = u[u["gene"] != ""].copy()
    u["gene"] = u["gene"].str.split(";").str[0].str.strip()

    kw = u["keywords"].str.lower()
    ipr = u["interpro"].str.split(";").apply(lambda xs: {x for x in xs if x})
    u["lipid_metabolism"] = kw.str.contains("lipid metabolism|fatty acid|sphingolipid|phospholipid|steroid")
    is_hydrolase = u["ec"].str.contains(r"(?:^|;\s*)3\.") | kw.str.contains("hydrolase")
    cat_ser = [_catalytic_ser(a, s) for a, s in zip(u["active_site"], u["sequence"])]
    cat_ser = pd.Series(cat_ser, index=u.index)
    sh_family = ipr.apply(lambda s: bool(s & SH_INTERPRO))
    patatin = ipr.apply(lambda s: "IPR002641" in s)  # PNPLAs, incl. transacylases like PNPLA1
    u["serine_hydrolase"] = (is_hydrolase & (cat_ser | sh_family)) | (sh_family & cat_ser) | patatin \
        | kw.str.contains("serine esterase")
    # Proteases are FP-reactive but not lipid enzymes; keep them flagged, rank separately
    u["protease"] = kw.str.contains("protease")

    u["rhea_ids"] = u["catalytic_activity"].apply(lambda s: sorted(set(re.findall(r"RHEA:(\d+)", s))))
    u["n_rhea"] = u["rhea_ids"].str.len()
    u["reactions"] = u["catalytic_activity"].apply(
        lambda s: " | ".join(dict.fromkeys(re.findall(r"Reaction=([^;]+);", s)))[:600]
    )
    text = (u["protein"] + " " + u["reactions"] + " " + u["keywords"]).str.lower()
    for tag, pat in PATHWAY_PATTERNS.items():
        u[f"pw_{tag}"] = text.str.contains(pat)
    u["pathways"] = u[[c for c in u if c.startswith("pw_")]].apply(
        lambda r: ",".join(c[3:] for c, v in r.items() if v), axis=1
    )
    u["chembl_target"] = u["chembl"].str.strip(";")
    u["has_chembl_target"] = u["chembl_target"] != ""
    u["location"] = u["location"].str.replace("SUBCELLULAR LOCATION: ", "", regex=False).str[:200]

    keep = ["gene", "uniprot", "protein", "ec", "serine_hydrolase", "protease", "lipid_metabolism",
            "pathways", "n_rhea", "rhea_ids", "reactions", "chembl_target", "has_chembl_target", "location"]
    out = u[keep].drop_duplicates("gene").set_index("gene").sort_index()
    # Restrict to enzymes relevant to lipid biology (lipid-metabolism keyword, or a
    # non-protease serine hydrolase, which may act on lipids without annotation yet)
    return out[out["lipid_metabolism"] | (out["serine_hydrolase"] & ~out["protease"])]
