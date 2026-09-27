import numpy as np
import pandas as pd

from insilico_pace import lipidome, rank
from insilico_pace.transcriptome import bh, contrast
from lipidname import parse


def test_bh_matches_reference():
    p = np.array([0.01, 0.04, 0.03, 0.5])
    # R: p.adjust(c(.01,.04,.03,.5), "BH") -> 0.04 0.0533 0.0533 0.5
    assert np.allclose(bh(p), [0.04, 0.16 / 3, 0.16 / 3, 0.5])


def test_contrast_direction():
    x = pd.DataFrame({"a1": [5, 1], "a2": [5.2, 1], "b1": [1, 1.1], "b2": [1.2, 0.9]}, index=["up", "flat"])
    r = contrast(x, ["a1", "a2"], ["b1", "b2"])
    assert r.loc["up", "lfc"] > 3 and abs(r.loc["flat", "lfc"]) < 0.2


def _toy_lipidome():
    names = ["Cer[NS] 42:1; Cer[NS](d18:1/24:0); [M-H]-", "Cer[NS] 34:1; Cer[NS](d18:1/16:0); [M-H]-",
             "SM 34:1; SM(d18:1/16:0); [M-H]-"]
    feats = pd.DataFrame({"metabolite_name": names, "mode": "neg"}, index=["f1", "f2", "f3"])
    feats["lipid"] = feats["metabolite_name"].map(parse)
    feats["cls"] = feats["lipid"].map(lambda l: l.lipid_class)
    idx = [f"{i}_{g}_{a}" for a in ("4w", "8w") for g in ("WT", "sp") for i in range(3)]
    idx = [f"{n}{s}" for n, s in zip(range(len(idx)), [s[s.index('_'):] for s in idx])]
    rng = np.random.default_rng(0)
    sp = np.array(["_sp_" in s for s in idx])
    x = pd.DataFrame({
        "f1": np.where(sp, 50, 100) * rng.uniform(0.95, 1.05, len(idx)),   # C24 Cer halves in disease
        "f2": np.full(len(idx), 100.0) * rng.uniform(0.95, 1.05, len(idx)),
        "f3": np.full(len(idx), 100.0) * rng.uniform(0.95, 1.05, len(idx)),
    }, index=idx)
    return feats, x


def test_proxy_detects_vlc_loss():
    feats, x = _toy_lipidome()
    design = lipidome.sample_design(x.index)
    eff = lipidome.genotype_effects(lipidome.compute_proxies(feats, x), design, ages=("4w", "8w"))
    assert eff.loc["VLC_ceramide", "lfc_spade"] < -0.8 and eff.loc["VLC_ceramide", "q"] < 0.05
    ev = lipidome.enzyme_lipid_evidence(eff.dropna(subset=["t"]))
    assert ev.loc["CERS2", "activity_t"] < 0 < ev.loc["CERS6", "activity_t"]


def test_linex_names():
    assert lipidome._to_linex(parse("Cer 18:1;O2/24:0")) == "Cer(d42:1)"
    assert lipidome._to_linex(parse("Plasmenyl-PE 40:6; PE(P-18:0/22:6); [M-H]-")) == "PE(O-40:7)"
    assert lipidome._to_linex(parse("Cer[AS] 40:1; Cer[AS](d18:1/22:0); [M-H]-")) is None


def test_rank_prefers_consistent_evidence():
    genes = ["GOOD", "NOISE", "OFF"]
    cat = pd.DataFrame({"serine_hydrolase": [True, False, False]}, index=genes)
    expr = pd.DataFrame({
        "AD_LS_vs_HC_lfc": [-2.0, 0.1, -2.0], "AD_LS_vs_HC_q": [1e-6, 0.9, 1e-6],
        "PSO_LS_vs_HC_lfc": [0, 0, 0], "KC_IL4_lfc": [-1, 0.1, -1], "KC_IL13_lfc": [-1, 0, -1],
        "KC_specificity": [0.8, 0.1, 0.8], "skin_log2cpm_HC": [5, 5, 0], "KC_nCPM": [50, 5, 0],
    }, index=genes)
    ev = pd.DataFrame({"activity_t": [-4.0, 1.0, -4.0]}, index=genes)
    r = rank.rank_enzymes(cat, expr, ev)
    assert r.index[0] == "GOOD"               # consistent across layers and expressed
    assert r.loc["OFF", "score"] < r.loc["GOOD", "score"]  # same evidence, not expressed in skin
