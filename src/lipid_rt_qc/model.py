"""Retention-time model fitted to the dataset's own confident annotations.

RT is column- and gradient-specific, so there is no universal model: the lipids
in each run are the training set. Two learners are averaged:

* per-class robust linear model, RT ~ C + DB (Huber). This is the equivalent-
  carbon-number rule written explicitly, so it extrapolates sensibly;
* gradient boosting on (class, C, DB, OH, ether, base DB), which captures
  curvature and class-specific DB effects.

Predictions are out-of-fold (K-fold), so a mis-annotated feature cannot pull
its own prediction towards itself.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import HuberRegressor
from sklearn.model_selection import KFold

from lipidname import Lipid

FEATURES = ["cls_code", "total_c", "total_db", "total_oh", "ether", "base_db", "n_chains"]


def expected_chains(cls: str) -> int:
    """Chains implied by the class (independent of how detailed the name is)."""
    if cls.startswith("Cer_EO") or cls in ("TG",) or cls.endswith("_EST"):
        return 3 + cls.endswith("_EST")
    if cls in ("CL",):
        return 4
    if cls.startswith(("L", "FA", "CAR", "CE", "MG", "NA", "SPB", "ST", "BA", "VAE")) and not cls.startswith("LacCer"):
        return 1
    return 2


def lipid_features(lipids: pd.Series) -> pd.DataFrame:
    rows = []
    for fid, l in lipids.items():
        if l is None:
            continue
        rows.append({
            "feature_id": fid, "lipid_class": l.lipid_class, "total_c": l.total_c, "total_db": l.total_db,
            "total_oh": l.total_oh, "ether": int(l.lipid_class.endswith(("-O", "-P"))),
            "base_db": l.base[1] if l.base else -1, "n_chains": expected_chains(l.lipid_class),
        })
    return pd.DataFrame(rows).set_index("feature_id")


class _CenteredHuber:
    """Huber regression on centred inputs (converges reliably; raw C ~ 30-60 does not)."""

    def fit(self, X: np.ndarray, y: np.ndarray):
        self.mx, self.my = X.mean(0), np.median(y)
        self.h = HuberRegressor(epsilon=1.35, alpha=1e-4, max_iter=2000).fit(X - self.mx, y - self.my)
        self.coef_ = self.h.coef_
        return self

    def predict(self, X) -> np.ndarray:
        return self.h.predict(np.asarray(X, float) - self.mx) + self.my


class RTModel:
    def __init__(self, min_class_n: int = 4, n_splits: int = 5, seed: int = 0):
        self.min_class_n, self.n_splits, self.seed = min_class_n, n_splits, seed
        self.classes_: list[str] = []

    # ---------------------------------------------------------------- fitting
    def _codes(self, cls: pd.Series) -> np.ndarray:
        lut = {c: i for i, c in enumerate(self.classes_)}
        return cls.map(lambda c: lut.get(c, -1)).to_numpy()

    def _gbm(self):
        # absolute-error loss: a mis-annotated feature must not drag its neighbours
        return HistGradientBoostingRegressor(loss="absolute_error", max_iter=200, learning_rate=0.05, max_leaf_nodes=15,
                                             min_samples_leaf=5, l2_regularization=1.0,
                                             categorical_features=[0], random_state=self.seed)

    def _linear(self, X: pd.DataFrame, y: np.ndarray):
        """Per-class Huber RT ~ C + DB; pooled slopes + class intercept for sparse classes."""
        lin = {}
        for cls, idx in X.groupby("lipid_class").groups.items():
            sub = X.loc[idx]
            if len(sub) >= self.min_class_n and sub["total_c"].nunique() >= 2:
                lin[cls] = _CenteredHuber().fit(sub[["total_c", "total_db"]].to_numpy(float),
                                                y[X.index.get_indexer(idx)])
        pooled = (np.median([h.coef_[0] for h in lin.values()]), np.median([h.coef_[1] for h in lin.values()])) if lin else None
        intercepts = {}
        if pooled:
            for cls, idx in X.groupby("lipid_class").groups.items():
                sub = X.loc[idx]
                intercepts[cls] = np.median(y[X.index.get_indexer(idx)] - pooled[0] * sub["total_c"] - pooled[1] * sub["total_db"])
        return lin, pooled, intercepts

    @staticmethod
    def _lin_predict(lin, pooled, intercepts, X: pd.DataFrame) -> np.ndarray:
        p = np.full(len(X), np.nan)
        for i, (cls, c, d) in enumerate(zip(X["lipid_class"], X["total_c"], X["total_db"])):
            if cls in lin:
                p[i] = lin[cls].predict([[c, d]])[0]
            elif pooled and cls in intercepts:
                p[i] = intercepts[cls] + pooled[0] * c + pooled[1] * d
        return p

    def _fit(self, X: pd.DataFrame, y: np.ndarray):
        """Linear ECN backbone per class + boosted correction of its residuals.

        Boosting on residuals (rather than on RT directly) means that outside a
        class's observed carbon range the prediction falls back to the linear
        trend instead of flat-lining, which is what gradient boosting does."""
        X = X.copy()
        X["cls_code"] = self._codes(X["lipid_class"])
        lin, pooled, intercepts = self._linear(X, y)
        base = self._lin_predict(lin, pooled, intercepts, X)
        have = ~np.isnan(base)
        corr = self._gbm().fit(X.loc[have, FEATURES].to_numpy(float), (y - base)[have]) if have.sum() > 20 else None
        raw = self._gbm().fit(X[FEATURES].to_numpy(float), y)  # for classes with no linear fit
        return lin, pooled, intercepts, corr, raw

    def _predict(self, parts, X: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
        lin, pooled, intercepts, corr, raw = parts
        Xc = X.copy()
        Xc["cls_code"] = self._codes(Xc["lipid_class"])
        F = Xc[FEATURES].to_numpy(float)
        base = self._lin_predict(lin, pooled, intercepts, X)
        pred = base + (corr.predict(F) if corr is not None else 0.0)
        known = Xc["cls_code"].to_numpy() >= 0
        fallback = np.isnan(base) & known
        pred[fallback] = raw.predict(F[fallback]) if fallback.any() else pred[fallback]
        return pred, base

    def fit(self, X: pd.DataFrame, y: pd.Series) -> "RTModel":
        self.classes_ = sorted(X["lipid_class"].unique())
        self.parts_ = self._fit(X, y.to_numpy(float))
        return self

    def predict(self, X: pd.DataFrame) -> np.ndarray:
        return self._predict(self.parts_, X)[0]

    def oof_predict(self, X: pd.DataFrame, y: pd.Series) -> np.ndarray:
        """Out-of-fold predictions for every training row."""
        self.classes_ = sorted(X["lipid_class"].unique())
        yv = y.to_numpy(float)
        pred = np.full(len(X), np.nan)
        k = min(self.n_splits, max(2, len(X) // 10))
        for tr, te in KFold(k, shuffle=True, random_state=self.seed).split(X):
            parts = self._fit(X.iloc[tr], yv[tr])
            pred[te] = self._predict(parts, X.iloc[te])[0]
        self.parts_ = self._fit(X, yv)
        return pred


def robust_z(resid: pd.Series, groups: pd.Series, min_group: int = 8) -> pd.Series:
    """Residual / robust SD. Uses the class's own spread when it has enough members,
    otherwise the global spread; never tighter than the global spread / 2."""
    r = resid.dropna()
    g_mad = 1.4826 * np.median(np.abs(r - np.median(r))) or 1e-6
    z = pd.Series(np.nan, index=resid.index)
    for cls, idx in groups.groupby(groups).groups.items():
        rr = resid.loc[idx]
        ok = rr.dropna()
        if len(ok) >= min_group:
            mad = max(1.4826 * np.median(np.abs(ok - np.median(ok))), g_mad / 2)
            z.loc[idx] = (rr - np.median(ok)) / mad
        else:
            z.loc[idx] = (rr - np.median(r)) / g_mad
    return z
