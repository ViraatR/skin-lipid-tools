"""Public inputs for the in silico PACE pipeline, and how to fetch them.

All sources are public (CC-BY / GEO / UniProt / HPA). Raw files land in
data/raw/ which is git-ignored; run `python -m insilico_pace download` once.
"""
from __future__ import annotations

import gzip
import io
import json
import zipfile
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
RESULTS = ROOT / "results" / "insilico_pace"

UNIPROT_QUERY = (
    "(organism_id:9606) AND (reviewed:true) AND ((keyword:KW-0443) OR (ec:3.1.1.*) "
    "OR (ec:3.1.4.*) OR (ec:3.5.1.*) OR (keyword:KW-0746))"
)
UNIPROT_FIELDS = (
    "accession,gene_primary,protein_name,ec,keyword,cc_catalytic_activity,xref_chembl,"
    "ft_act_site,sequence,xref_interpro,cc_subcellular_location"
)

GEO = {
    # AD + psoriasis + healthy RNA-seq (Tsoi et al. 2019, JID)
    "GSE121212": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121212/suppl/GSE121212_readcount.txt.gz",
    # Primary keratinocytes + IL-4, IL-13, IL-17A, IFN-a, IFN-g, TNF (GPL570)
    "GSE36287": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE36nnn/GSE36287/matrix/GSE36287_series_matrix.txt.gz",
    # Acne lesional vs non-lesional (GPL570), used as a side contrast
    "GSE53795": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53795/matrix/GSE53795_series_matrix.txt.gz",
}
GEO_META = {
    "GSE121212": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121212/matrix/GSE121212_series_matrix.txt.gz",
}
GPL570_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
HPA_SC_CLUSTER = "https://www.proteinatlas.org/download/tsv/rna_single_cell_cluster.tsv.zip"
# Mouse Tmem79 (Spade) dermatitis-model skin lipidomics, RPLC-MS, Keio (MS-DIAL names)
MWB_LIPIDOME = "ST002195"
# Human AD vs control skin biopsies, MALDI direct infusion (sparse; cross-check only)
MWB_HUMAN = "ST003291"


def _get(url: str, dest: Path, timeout: int = 600) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    print(f"  fetching {url}")
    # Stream to a temp file so an interrupted download isn't mistaken for a cached one
    part = dest.with_name(dest.name + ".part")
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with open(part, "wb") as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
    part.replace(dest)
    return dest


def download() -> None:
    print("UniProt enzyme universe")
    dest = RAW / "uniprot" / "human_lipid_hydrolases.tsv"
    if not dest.exists():
        r = requests.get(
            "https://rest.uniprot.org/uniprotkb/stream",
            params={"query": UNIPROT_QUERY, "format": "tsv", "fields": UNIPROT_FIELDS},
            timeout=600,
        )
        r.raise_for_status()
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(r.text)

    print("GEO transcriptomes")
    for acc, url in GEO.items():
        _get(url, RAW / "geo" / Path(url).name)
    for acc, url in GEO_META.items():
        _get(url, RAW / "geo" / Path(url).name)
    _get(GPL570_ANNOT, RAW / "geo" / "GPL570.annot.gz")

    print("HPA skin single-cell clusters (~200 MB zip, skin rows kept)")
    skin = RAW / "hpa" / "skin_clusters_ntpm.tsv"
    if not skin.exists():
        z = _get(HPA_SC_CLUSTER, RAW / "hpa" / "rna_single_cell_cluster.tsv.zip")
        with zipfile.ZipFile(z) as zf, open(skin, "w") as out:
            name = zf.namelist()[0]
            with zf.open(name) as fh:
                for i, line in enumerate(io.TextIOWrapper(fh)):
                    if i == 0 or "\tskin\t" in line:
                        out.write(line)
        z.unlink()

    print("Metabolomics Workbench lipidomes")
    for st in (MWB_LIPIDOME, MWB_HUMAN):
        dest = RAW / "mwb" / f"{st}.json"
        if not dest.exists():
            r = requests.get(f"https://www.metabolomicsworkbench.org/rest/study/study_id/{st}/data", timeout=600)
            r.raise_for_status()
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(r.text)
        fdest = RAW / "mwb" / f"{st}_factors.json"
        if not fdest.exists():
            r = requests.get(f"https://www.metabolomicsworkbench.org/rest/study/study_id/{st}/factors", timeout=600)
            r.raise_for_status()
            fdest.write_text(r.text)
    print("done")


def read_series_matrix(path: Path):
    """Return (sample metadata DataFrame, expression DataFrame) from a GEO series matrix."""
    import pandas as pd

    meta: dict[str, list[str]] = {}
    rows = []
    with gzip.open(path, "rt") as fh:
        in_table = False
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                rows.append(line.rstrip("\n").split("\t"))
            elif line.startswith("!Sample_"):
                parts = line.rstrip("\n").split("\t")
                key = parts[0][8:]
                vals = [p.strip('"') for p in parts[1:]]
                while key in meta:
                    key += "_"
                meta[key] = vals
    m = pd.DataFrame(meta)
    expr = None
    if rows:
        header = [h.strip('"') for h in rows[0]]
        expr = pd.DataFrame([[c.strip('"') for c in r] for r in rows[1:]], columns=header)
        expr = expr.set_index(header[0]).apply(pd.to_numeric, errors="coerce")
    return m, expr


def load_mwb(study: str):
    """Metabolomics Workbench study -> (feature table, sample x feature matrix)."""
    import pandas as pd

    d = json.loads((RAW / "mwb" / f"{study}.json").read_text())
    rows = list(d.values())
    feats = pd.DataFrame(
        [{k: r[k] for k in ("metabolite_id", "metabolite_name", "refmet_name", "analysis_id", "analysis_summary")} for r in rows]
    ).set_index("metabolite_id")
    mat = pd.DataFrame({r["metabolite_id"]: r["DATA"] for r in rows}).apply(pd.to_numeric, errors="coerce")
    return feats, mat
