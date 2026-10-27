import csv

import numpy as np
import pandas as pd
import pytest

from lipid_rt_qc import checks, report
from lipid_rt_qc.cli import main
from lipid_rt_qc.io import read_any, read_msdial, read_skyline
from lipid_rt_qc.masses import enumerate_isobars, formula, mz
from lipid_rt_qc.model import expected_chains
from lipidname import parse

# Reverse-phase-like RT: rises with carbons, falls with double bonds, class offset
OFFSET = {"PC": 0.0, "PE": 0.3, "TG": 5.0, "Cer_NS": 1.0, "Cer_NDS": 1.25, "LPC": -6.0}
NAME = {"PC": "PC {c}:{d}", "PE": "PE {c}:{d}", "TG": "TG {c}:{d}", "LPC": "LPC {c}:{d}",
        "Cer_NS": "Cer 18:1;O2/{a}:{ad}", "Cer_NDS": "Cer 18:0;O2/{a}:{ad}"}
ADD = {"PC": "[M+H]+", "PE": "[M+H]+", "TG": "[M+NH4]+", "LPC": "[M+H]+", "Cer_NS": "[M+H]+", "Cer_NDS": "[M+H]+"}


def true_rt(cls, c, d):
    return 2.0 + OFFSET[cls] + 0.25 * c - 0.45 * d


def synthetic(seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for cls, cs, ds in [("PC", range(30, 41), range(0, 6)), ("PE", range(32, 41), range(0, 6)),
                        ("TG", range(46, 57, 2), range(0, 7)), ("LPC", range(14, 23), range(0, 3)),
                        ("Cer_NS", range(34, 45), range(1, 3)), ("Cer_NDS", range(34, 45), range(0, 2))]:
        for c in cs:
            for d in ds:
                if cls.startswith("Cer"):
                    base_db = 1 if cls == "Cer_NS" else 0
                    name = NAME[cls].format(a=c - 18, ad=d - base_db)
                else:
                    name = NAME[cls].format(c=c, d=d)
                rows.append({"name": name, "cls": cls, "c": c, "d": d, "adduct": ADD[cls],
                             "rt": true_rt(cls, c, d) + rng.normal(0, 0.04),
                             "mz": mz(cls, c, d, ADD[cls])})
    df = pd.DataFrame(rows)
    # plant two errors: a wrong peak and a 2-carbon mislabel
    df.loc[df.index[10], "rt"] += 3.0
    df.loc[df.index[40], "name"] = "PC 38:2"
    df.loc[df.index[40], "c"] = 38
    return df, {df.index[10], df.index[40]}


def write_msdial(df: pd.DataFrame, path, n_samples=4, seed=0):
    rng = np.random.default_rng(seed)
    meta_cols = ["Alignment ID", "Average Rt(min)", "Average Mz", "Metabolite name", "Adduct type",
                 "Ontology", "MS/MS matched"]
    samples = ["Blank_1", "Blank_2"] + [f"S{i}" for i in range(n_samples)] + ["QC_1", "QC_2", "QC_3"]
    types = ["Blank", "Blank"] + ["Sample"] * n_samples + ["QC"] * 3
    lab = len(meta_cols) - 1
    top = []
    for label, vals in (("Class", ["blank"] * 2 + ["group"] * n_samples + ["QC"] * 3),
                        ("File type", types), ("Injection order", list(range(1, len(samples) + 1))),
                        ("Batch ID", [1] * len(samples))):
        top.append([""] * lab + [label] + [str(v) for v in vals] + ["", ""])
    header = meta_cols + samples + ["Average", "Stdev"]
    body = []
    for i, r in df.iterrows():
        inten = rng.lognormal(10, 0.3) * rng.normal(1, 0.05, len(samples))
        inten[:2] = rng.lognormal(5, 0.3, 2)  # blanks low...
        if i == 5:
            inten[:2] = inten[2:4]            # ...except one contaminant
        if i == 7:
            inten[-3:] = [1e4, 5e4, 2e5]      # ...and one imprecise QC feature
        body.append([str(i), f"{r.rt:.3f}", f"{r.mz:.5f}", r["name"], r.adduct, "", "True"]
                    + [f"{v:.1f}" for v in inten] + ["0", "0"])
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerows(top + [header] + body)
    return samples


@pytest.fixture
def msdial_file(tmp_path):
    df, planted = synthetic()
    p = tmp_path / "Height_0_2026.txt"
    write_msdial(df, p)
    return p, df, planted


def test_formula_and_isobar_basics():
    assert str(formula("PC", 34, 1)) == "C42H82NO8P"
    assert abs(mz("PC", 34, 1, "[M+H]+") - 760.5851) < 1e-3
    hits = {(c, t, d, a) for c, t, d, a, _ in enumerate_isobars(mz("PC", 34, 1, "[M+H]+"), {"[M+H]+", "[M+Na]+"}, ppm=3)}
    assert ("PE", 37, 1, "[M+H]+") in hits


# LIPID MAPS reference formulas, one per class rule family
@pytest.mark.parametrize("cls,c,d,expected", [
    ("PE", 34, 1, "C39H76NO8P"),
    ("PS", 34, 1, "C40H76NO10P"),
    ("PG", 34, 1, "C40H77O10P"),
    ("PI", 34, 1, "C43H81O13P"),
    ("PA", 34, 1, "C37H71O8P"),
    ("PE-O", 34, 1, "C39H78NO7P"),  # PE O-16:0/18:1
    ("PE-P", 34, 1, "C39H76NO7P"),  # PE P-16:0/18:1, vinyl-ether DB not counted in y
    ("LPC", 16, 0, "C24H50NO7P"),
    ("TG", 52, 2, "C55H102O6"),
    ("DG", 34, 1, "C37H70O5"),
    ("CE", 18, 1, "C45H78O2"),
    ("FA", 16, 0, "C16H32O2"),
    ("CAR", 16, 0, "C23H45NO4"),
    ("SM", 34, 1, "C39H79N2O6P"),
    ("Cer_NS", 34, 1, "C34H67NO3"),
    ("HexCer_NS", 34, 1, "C40H77NO8"),
])
def test_formula_matches_reference(cls, c, d, expected):
    assert str(formula(cls, c, d)) == expected


def test_formula_unknown_class_or_adduct():
    assert formula("NotAClass", 34, 1) is None
    assert mz("PC", 34, 1, "[M+Xe]+") is None


def test_enumerate_isobars_rejects_unknown_adduct_or_class():
    with pytest.raises(ValueError, match=r"\[M\+Xe\]\+"):
        enumerate_isobars(760.585, {"[M+H]+", "[M+Xe]+"})
    with pytest.raises(ValueError, match="NotAClass"):
        enumerate_isobars(760.585, {"[M+H]+"}, classes=["PC", "NotAClass"])


def test_read_msdial(msdial_file):
    p, df, _ = msdial_file
    ds = read_msdial(p)
    assert len(ds.features) == len(df)
    assert ds.samples["type"].value_counts().to_dict() == {"sample": 4, "qc": 3, "blank": 2}
    assert "Average" not in ds.intensity.columns
    assert read_any(p).source.startswith("MS-DIAL")


def test_read_msdial_rejects_non_msdial_file(tmp_path):
    p = tmp_path / "not_msdial.txt"
    p.write_text("name\trt\tmz\nPC 34:1\t5.1\t760.585\n")
    with pytest.raises(ValueError, match="Alignment ID"):
        read_msdial(p)
    p.write_text("Alignment ID\tMetabolite name\n1\tPC 34:1\n")
    with pytest.raises(ValueError, match="Average Rt"):
        read_msdial(p)


def test_read_generic_csv(tmp_path):
    p = tmp_path / "feats.csv"
    p.write_text("name,rt,mz,QC_1,Blank_1,S1\nPC 34:1,5.1,760.585,10,1,12\n")
    ds = read_any(p)
    assert ds.source.startswith("CSV")
    assert ds.features.iloc[0]["rt"] == 5.1
    assert ds.samples["type"].to_dict() == {"QC_1": "qc", "Blank_1": "blank", "S1": "sample"}
    p.write_text("Compound,RT\nPC 34:1,5.1\n")
    with pytest.raises(ValueError, match="'name'"):
        read_any(p)


def test_read_maf(tmp_path):
    p = tmp_path / "m_lipids.tsv"
    cols = ["metabolite_identification", "mass_to_charge", "retention_time", "modifications",
            "smallmolecule_abundance_sub", "smallmolecule_abundance_std_error_sub", "S1", "QC_1", "Blank_1"]
    rows = [["PC 34:1", "760.585", "306", "[M+H]+", "", "", "100", "90", ""],
            ["TG 52:2", "876.802", "870", "unknown", "", "", "50", "45", "2"]]
    p.write_text("\n".join("\t".join(r) for r in [cols, *rows]) + "\n")
    ds = read_any(p)
    assert ds.source.startswith("MetaboLights MAF")
    # retention times submitted in seconds are converted to minutes
    assert ds.features["rt"].tolist() == pytest.approx([5.1, 14.5])
    # only adduct-shaped modifications are kept
    assert ds.features["adduct"].tolist() == ["[M+H]+", ""]
    # sample columns are the ones after the smallmolecule_abundance block
    assert ds.samples["type"].to_dict() == {"S1": "sample", "QC_1": "qc", "Blank_1": "blank"}
    assert np.isnan(ds.intensity.iloc[0]["Blank_1"])
    p.write_text("metabolite_identification\tmass_to_charge\tS1\nPC 34:1\t760.585\t100\n")
    with pytest.raises(ValueError, match="'retention_time'"):
        read_any(p)


@pytest.mark.parametrize("names", [["foo"], ["PC 34:1", "foo"]])
def test_rt_model_needs_two_annotated_features(tmp_path, names):
    p = tmp_path / "feats.csv"
    p.write_text("name,rt,mz\n" + "".join(f"{n},5.1,760.585\n" for n in names))
    f = checks.annotate(read_any(p))
    with pytest.raises(ValueError, match="at least 2 features"):
        checks.rt_outliers(f)


def test_planted_errors_are_found(msdial_file):
    p, df, planted = msdial_file
    res, _ = checks.run_all(read_msdial(p))
    flagged = set(res.index[res["flag_rt"]].astype(int))
    assert planted <= flagged
    assert len(flagged - planted) <= 3          # few false positives on clean synthetic data
    assert res.loc["5", "flag_blank"] and res.loc["7", "flag_qc_cv"]
    assert res["rt_resid"].abs().median() < 0.1


def test_cer_dhcer_order_rule(msdial_file):
    p, df, _ = msdial_file
    ds = read_msdial(p)
    # push one Cer[NS] after its DHCer partner
    idx = df.index[(df["cls"] == "Cer_NS") & (df["c"] == 40) & (df["d"] == 1)][0]
    ds.features.loc[str(idx), "rt"] += 1.5
    res, _ = checks.run_all(ds)
    assert "DHCer" in res.loc[str(idx), "order_notes"]
    assert res.loc[str(idx), "flag_order"]


def test_bcfa_isomer_call(msdial_file):
    p, df, _ = msdial_file
    ds = read_msdial(p)
    # add an earlier-eluting copy of PC 34:0 (~0.5 C-eq earlier) = mmBCFA-like isomer
    src = ds.features[ds.features["name"] == "PC 34:0"].iloc[0]
    new = src.copy()
    new["rt"] = src["rt"] - 0.5 * 0.25 * 1.0
    ds.features.loc["bcfa"] = new
    res, _ = checks.run_all(ds)
    assert "branched" in res.loc["bcfa", "dup_call"]


def test_skyline_long_format(tmp_path):
    df, _ = synthetic()
    rows = []
    for rep in ["Blank01", "Sample01", "Sample02", "QC01", "QC02"]:
        for _, r in df.iterrows():
            rows.append({"Molecule Name": r["name"], "Precursor Adduct": r.adduct, "Precursor Mz": r.mz,
                         "Replicate Name": rep, "Retention Time": r.rt, "Area": 1000.0})
    p = tmp_path / "skyline.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    ds = read_skyline(p)
    assert len(ds.features) == len(df.drop_duplicates(["name", "adduct"]))
    assert ds.samples["type"].tolist() == ["blank", "qc", "qc", "sample", "sample"]
    res, _ = checks.run_all(ds)
    assert res["flag_rt"].sum() >= 2


def test_cli_and_report(msdial_file, tmp_path, capsys):
    p, _, _ = msdial_file
    out = tmp_path / "qc"
    main(["run", str(p), "--out", str(out)])
    html = (out / "lipid_qc_report.html").read_text()
    assert "Retention-time model" in html and "data:image/png;base64" in html
    flagged = pd.read_csv(out / "lipid_qc_flagged.csv")
    assert {"name", "reasons"} <= set(flagged.columns)



def test_report_uses_cli_thresholds(msdial_file, tmp_path, monkeypatch):
    p, _, _ = msdial_file
    seen = {}
    orig = report.figures

    def spy(res, params=None):
        seen.update(params or {})
        return orig(res, params)

    monkeypatch.setattr(report, "figures", spy)
    out = tmp_path / "qc"
    main(["run", str(p), "--out", str(out), "--z", "4.5", "--max-cv", "20", "--blank-fold", "5"])
    assert seen["max_cv"] == 20 and seen["blank_fold"] == 5
    assert "&gt;4.5 robust SD" in (out / "lipid_qc_report.html").read_text()

@pytest.mark.parametrize("extra", [["--ppm", "-5"], ["--z", "0"], ["--max-cv", "abc"]])
def test_cli_rejects_non_positive_thresholds(msdial_file, tmp_path, extra):
    p, _, _ = msdial_file
    with pytest.raises(SystemExit) as e:
        main(["run", str(p), "--out", str(tmp_path / "qc"), *extra])
    assert e.value.code == 2


@pytest.mark.parametrize("extra", [["--blank", "blank("], ["--qc", "[qc"]])
def test_cli_rejects_invalid_sample_regex(msdial_file, tmp_path, capsys, extra):
    p, _, _ = msdial_file
    with pytest.raises(SystemExit) as e:
        main(["run", str(p), "--out", str(tmp_path / "qc"), *extra])
    assert e.value.code == 2
    assert "invalid regex" in capsys.readouterr().err


def test_cli_reports_missing_input(tmp_path, capsys):
    with pytest.raises(SystemExit) as e:
        main(["run", str(tmp_path / "nope.txt")])
    assert e.value.code == 2
    assert "input file not found" in capsys.readouterr().err


def test_blank_only_feature_is_flagged():
    from lipid_rt_qc.io import Dataset
    feats = pd.DataFrame({"name": ["PC 34:1", "PC 36:2", "PE 38:4"], "rt": [10.0, 10.5, 11.0],
                          "mz": np.nan, "adduct": ""}, index=["a", "b", "c"])
    inten = pd.DataFrame({"Blank_1": [100.0, 5e4, np.nan], "S1": [1e5, np.nan, 1e5], "S2": [1.2e5, np.nan, 9e4]},
                         index=feats.index)
    samples = pd.DataFrame({"type": ["blank", "sample", "sample"]}, index=inten.columns)
    out = checks.blank_and_qc(Dataset(feats, inten, samples))
    # b is only in the blank; c is never in a blank, so it has nothing to compare against
    assert out["flag_blank"].tolist() == [False, True, False]


@pytest.mark.parametrize("name,n", [
    ("TG 18:1_18:2_18:1;O(FA 16:0)", 4), ("DG 16:0_18:1;O(FA 18:0)", 3),
    ("SM 34:1;O2(FA 18:0)", 3), ("TG 52:2", 3), ("PC 34:1", 2), ("LPC 18:1", 1),
])
def test_expected_chains_counts_estolide_fa_once(name, n):
    assert expected_chains(parse(name).lipid_class) == n


@pytest.mark.parametrize("name,n", [
    ("FAHFA 38:5; FAHFA(20:4/18:1); [M-H]-", 2), ("MLCL 52:2", 3), ("CL 72:8", 4),
    ("NAPE 52:4", 3), ("NAE 20:4", 1), ("MGDG 34:3", 2), ("MG 18:1", 1),
])
def test_expected_chains_for_classes_with_misleading_prefixes(name, n):
    assert expected_chains(parse(name).lipid_class) == n


def test_isobar_alternatives_use_class_features_from_the_run():
    X = pd.DataFrame({
        "lipid_class": ["Cer_NP", "Cer_NP", "PC"], "total_c": [42, 40, 34], "total_db": [0, 0, 1],
        "total_oh": [1, 1, 0], "ether": [0, 0, 0], "base_db": [0, 0, -1], "n_chains": [2, 2, 2],
    })
    cand = pd.DataFrame({"lipid_class": ["Cer_NP", "Cer_NS"], "total_c": [44, 44], "total_db": [0, 1],
                         "adduct": ["[M+H]+"] * 2, "ppm": [0.0, 0.0]})
    out = checks._alt_features(cand, X).set_index("lipid_class")
    assert (out.loc["Cer_NP", ["total_oh", "base_db"]] == [1, 0]).all()  # learnt from the run
    assert (out.loc["Cer_NS", ["total_oh", "base_db", "n_chains"]] == [0, 1, 2]).all()  # name-based fallback



def test_isobar_flag_follows_z_threshold():
    # PC 34:1 [M+H]+ is off its RT trend (z = 4); the PE 37:1 isobar fits exactly
    f = pd.DataFrame({"name": ["PC 34:1"], "rt": [10.0], "mz": [mz("PC", 34, 1, "[M+H]+")], "adduct": ["[M+H]+"]},
                     index=["a"])
    f["lipid"] = f["name"].map(parse)
    rt = pd.DataFrame({"rt_resid": [2.0], "rt_z": [4.0]}, index=f.index)
    X = pd.DataFrame({"lipid_class": ["PC", "PE"], "total_c": [34, 36], "total_db": [1, 1], "total_oh": [0, 0],
                      "ether": [0, 0], "base_db": [-1, -1], "n_chains": [2, 2]})

    class AltFitsRT:
        def predict(self, cand):
            return np.full(len(cand), 10.0)

    flag = lambda z: checks.isobars(f, rt, AltFitsRT(), X, z_thresh=z).at["a", "flag_isobar_rt"]
    assert flag(3.0)
    assert not flag(5.0)


@pytest.mark.parametrize("kind", ["carbon_pm2", "db_pm1", "class_isobar", "wrong_peak"])
def test_inject_mutates_only_the_reported_features(msdial_file, kind):
    from lipid_rt_qc.evaluate import inject

    p, _, _ = msdial_file
    ds = read_msdial(p)
    before = ds.features.copy()
    bad, hit = inject(ds, kind, frac=0.05, seed=0)
    pd.testing.assert_frame_equal(ds.features, before)   # input left untouched
    assert len(hit) >= 5
    changed = bad.features.index[(bad.features[["name", "rt"]] != before[["name", "rt"]]).any(axis=1)]
    assert set(changed) == set(hit)
    for fid in hit:
        old, new = parse(before.at[fid, "name"]), parse(bad.features.at[fid, "name"])
        if kind == "wrong_peak":
            assert new == old and abs(bad.features.at[fid, "rt"] - before.at[fid, "rt"]) >= 0.5
        elif kind == "carbon_pm2":
            assert (new.lipid_class, abs(new.total_c - old.total_c), new.total_db) == (old.lipid_class, 2, old.total_db)
        elif kind == "db_pm1":
            assert (new.lipid_class, new.total_c, abs(new.total_db - old.total_db)) == (old.lipid_class, old.total_c, 1)
        else:
            assert {old.lipid_class, new.lipid_class} == {"PC", "PE"} and abs(new.total_c - old.total_c) == 3
    again, hit2 = inject(ds, kind, frac=0.05, seed=0)
    assert list(hit2) == list(hit)                        # same seed, same errors


def test_benchmark_detects_wrong_peaks(msdial_file):
    from lipid_rt_qc.evaluate import benchmark, summarise

    p, _, _ = msdial_file
    b = benchmark(read_msdial(p), kinds=["wrong_peak", "carbon_pm2"], seeds=(0,))
    s = summarise(b)
    assert list(s.index) == ["carbon_pm2", "wrong_peak"]   # reported in KINDS order
    assert s.loc["wrong_peak", "auc"] > 0.9
    assert s.loc["wrong_peak", "recall"] > 0.5
    assert s.loc["wrong_peak", "fpr_clean"] < 0.1
