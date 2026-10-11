"""Integrate the layers into one ranked list of candidate AD lipid enzymes.

Score (0-1) = 0.40 AD transcript  + 0.20 Th2 keratinocyte response
            + 0.30 lipidome-inferred activity + 0.10 keratinocyte specificity,
then halved for enzymes barely expressed in skin. Each component is a
percentile so no single dataset's scale dominates. Weights are a judgment
call and are exposed as arguments.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# Known biology the ranking should recover (see README for sources)
POSITIVE_CONTROLS = {
    "ELOVL3": "IL-4/IL-13 suppress via STAT6; down in AD skin",
    "ELOVL6": "IL-4/IL-13 suppress via STAT6; down in AD differentiated KCs (lab scRNA reanalysis)",
    "ASAH1": "beta subunit = SM deacylase; shifts ceramide supply in AD",
    "SCD": "up in AD lesional skin (lab scRNA reanalysis)",
    "CERS3": "ultra-long-chain ceramide synthase; ULC Cer lost in AD",
    "ELOVL4": "ultra-long-chain elongase; ULC Cer lost in AD",
}


def _pct(s: pd.Series) -> pd.Series:
    return s.rank(pct=True).fillna(0)


def rank_enzymes(catalog: pd.DataFrame, expr: pd.DataFrame, lipid_ev: pd.DataFrame,
                 linex: pd.DataFrame | None = None, weights=(0.40, 0.20, 0.30, 0.10)) -> pd.DataFrame:
    weights = tuple(float(w) for w in weights)
    if len(weights) != 4:
        raise ValueError(f"weights must be (ad, th2, lipid, kc); got {len(weights)} values")
    if any(w < 0 for w in weights) or not np.isclose(sum(weights), 1.0):
        raise ValueError(f"weights must be non-negative and sum to 1 so the score stays in 0-1; got {weights}")
    df = catalog.join(expr, how="left").join(lipid_ev, how="left")
    if linex is not None and len(linex):
        df = df.join(linex, how="left")
    df["linex_edges"] = df.get("linex_edges", pd.Series(0, index=df.index)).fillna(0)

    q = df["AD_LS_vs_HC_q"].clip(lower=1e-10)
    df["ad_evidence"] = df["AD_LS_vs_HC_lfc"].abs() * -np.log10(q)
    df["th2_lfc"] = df[["KC_IL4_lfc", "KC_IL13_lfc"]].mean(axis=1)
    agree_th2 = np.sign(df["th2_lfc"]) == np.sign(df["AD_LS_vs_HC_lfc"])
    # LINEX2 membership is supporting evidence (worth at most t = 2), not a driver
    lip = df["activity_t"].abs().fillna(0) + np.minimum(df["linex_edges"], 10) / 5
    # Lipid evidence counts fully when inferred activity moves the same way as the transcript
    agree_lip = (np.sign(df["activity_t"]) == np.sign(df["AD_LS_vs_HC_lfc"])) | df["activity_t"].isna()
    df["lipid_consistent"] = np.where(df["activity_t"].isna(), np.nan, agree_lip.astype(float))

    s_ad = _pct(df["ad_evidence"])
    s_th2 = _pct(df["th2_lfc"].abs()) * np.where(agree_th2, 1.0, 0.5)
    s_lip = (_pct(lip.where(lip > 0)) * np.where(agree_lip, 1.0, 0.5)).fillna(0)
    s_kc = df["KC_specificity"].fillna(0).clip(0, 1)
    w_ad, w_th2, w_lip, w_kc = weights
    score = w_ad * s_ad + w_th2 * s_th2 + w_lip * s_lip + w_kc * s_kc

    expressed = (df["skin_log2cpm_HC"].fillna(0) >= 1) | (df["KC_nCPM"].fillna(0) >= 5)
    df["expressed_in_skin"] = expressed
    df["score"] = np.where(expressed, score, score * 0.5)
    # Enzymes with no AD transcript measurement get no direction rather than "down"
    lfc = df["AD_LS_vs_HC_lfc"]
    df["ad_direction"] = np.select([lfc > 0, lfc < 0], ["up", "down"], default="")
    df["ad_vs_pso_lfc"] = df["AD_LS_vs_HC_lfc"] - df["PSO_LS_vs_HC_lfc"]
    df = df.sort_values("score", ascending=False)
    df.insert(0, "rank", np.arange(1, len(df) + 1))
    df["positive_control"] = df.index.map(POSITIVE_CONTROLS).fillna("")
    return df


def control_recovery(ranked: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """Where do the positive controls land, and is that better than chance?"""
    from scipy.stats import mannwhitneyu

    n = len(ranked)
    ctl = ranked.loc[ranked.index.isin(POSITIVE_CONTROLS)]
    out = ctl[["rank", "score", "AD_LS_vs_HC_lfc", "th2_lfc", "lipid_proxy", "activity_t"]].copy()
    out["percentile"] = 100 * (1 - (out["rank"] - 1) / n)
    rest = ranked.loc[~ranked.index.isin(POSITIVE_CONTROLS), "score"]
    p = mannwhitneyu(ctl["score"], rest, alternative="greater").pvalue
    return out, p


REPORT_COLS = [
    "rank", "score", "serine_hydrolase", "pathways", "ad_direction", "AD_LS_vs_HC_lfc", "AD_LS_vs_HC_q",
    "AD_LS_vs_NL_lfc", "PSO_LS_vs_HC_lfc", "KC_IL4_lfc", "KC_IL13_lfc", "KC_IFNg_lfc", "KC_IL17A_lfc",
    "ACNE_LS_vs_NL_lfc", "lipid_proxy", "activity_t", "lipid_q", "lipid_consistent", "linex_edges",
    "top_skin_cell_type", "KC_nCPM", "KC_specificity", "suprabasal_vs_basal_log2", "has_chembl_target",
    "chembl_target", "uniprot", "protein", "positive_control",
]
