"""lipid-rt-qc: flag suspicious lipid annotations in MS-DIAL / Skyline / MetaboLights exports.

  lipid-rt-qc run export.txt --out qc/
  lipid-rt-qc run skyline_report.csv --format skyline --blank "blank|PBS" --qc "QC"
  lipid-rt-qc benchmark export.txt --out bench/
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="lipid-rt-qc", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "benchmark"):
        p = sub.add_parser(name)
        p.add_argument("input", type=Path)
        p.add_argument("--out", type=Path, default=Path("lipid_qc"))
        p.add_argument("--format", default="auto", choices=["auto", "msdial", "skyline", "maf", "csv"])
        p.add_argument("--blank", default=r"blank|blk|solvent", help="regex for blank sample names")
        p.add_argument("--qc", default=r"\bqc|pool", help="regex for pooled-QC sample names")
    r = sub.choices["run"]
    r.add_argument("--z", type=float, default=3.0, help="robust-z threshold for RT outliers")
    r.add_argument("--ppm", type=float, default=5.0, help="m/z tolerance for isobar search")
    r.add_argument("--blank-fold", type=float, default=3.0, help="min sample/blank ratio")
    r.add_argument("--max-cv", type=float, default=30.0, help="max QC CV %%")
    a = ap.parse_args(argv)

    from .checks import run_all
    from .io import read_any

    ds = read_any(a.input, a.format, blank_pat=a.blank, qc_pat=a.qc)
    print(f"{ds.source}: {len(ds.features)} features, "
          f"{ds.samples['type'].value_counts().to_dict() if len(ds.samples) else 'no intensity columns'}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if a.cmd == "run":
            from . import report

            params = dict(z=a.z, ppm=a.ppm, blank_fold=a.blank_fold, max_cv=a.max_cv)
            res, _ = run_all(ds, z_thresh=a.z, ppm=a.ppm, blank_fold=a.blank_fold, max_cv=a.max_cv)
            path = report.write(res, a.out, ds.source, params)
            flags = [c for c in res if c.startswith("flag_")]
            print(res[flags].sum().to_string())
            print(f"report: {path}")
        else:
            from .evaluate import benchmark, summarise

            a.out.mkdir(parents=True, exist_ok=True)
            b = benchmark(ds)
            b.to_csv(a.out / "benchmark_runs.csv", index=False)
            s = summarise(b)
            s.to_csv(a.out / "benchmark_summary.csv")
            print(s.round(3).to_string())


if __name__ == "__main__":
    main()
