# skin-lipid-tools

Two computational projects I built for skin lipid biology, both using public data.

| | Project | Question | Status |
|---|---|---|---|
| 1 | [`insilico_pace`](src/insilico_pace) | Using only public data, which skin enzymes should PACE find dysregulated in atopic dermatitis? | Pipeline runs end to end; results are in [`results/insilico_pace`](results/insilico_pace/SUMMARY.md) |
| 2 | [`lipid_rt_qc`](src/lipid_rt_qc) | Which lipid annotations in an MS-DIAL or Skyline export look suspicious? | CLI works; validated on 3 public runs ([results](results/lipid_rt_qc)) |

Both projects share [`lipidname`](src/lipidname/__init__.py), a lipid-name parser. It reads MS-DIAL names
(`Cer[NS] 50:3; Cer[NS](d18:1/32:2); [M+Hac-H]-`), LIPID MAPS shorthand (`Cer 18:1;O2/24:0`) and
RefMet/Skyline names, and turns each into a class, a composition, the chains, the sphingoid base and the adduct.

```bash
python -m venv .venv && .venv/bin/pip install -e ".[linex,dev]"
.venv/bin/insilico-pace all        # download public data (~250 MB) and run Project 1
.venv/bin/pytest
```

---

## Project 1: in silico PACE

PACE (Profiling of Active Cutaneous Enzymes) pairs FP-probe chemoproteomics with lipidomics.
This project is the computational side of that idea. I built a catalog of enzymes that could matter,
scored each one against the public datasets that say something about AD, and ranked them.

### Layers

**Layer 1: catalog and expression** (`catalog.py`, `transcriptome.py`)

| Source | What it adds |
|---|---|
| UniProt (reviewed human; lipid metabolism or hydrolase EC classes) | 764 enzymes, Rhea reactions, ChEMBL targets. **99 FP-reactive serine hydrolases**, called from a catalytic Ser plus serine-hydrolase InterPro folds |
| GSE121212 (Tsoi 2019): AD, psoriasis and healthy skin RNA-seq | AD lesional vs healthy (21 vs 38), AD non-lesional vs healthy, paired AD LS vs NL, psoriasis vs healthy |
| GSE36287: primary keratinocytes + IL-4, IL-13, IL-17A, IFN-α, IFN-γ, TNF | Th2 and other cytokine responses (n=3 per arm, paired by donor) |
| GSE53795: acne lesional vs non-lesional | Side contrast (acne, not AD) |
| Human Protein Atlas skin scRNA-seq | Top skin cell type, keratinocyte expression and specificity, suprabasal vs basal |

**Layer 2: lipidome-inferred enzyme activity** (`lipidome.py`)

The data is Metabolomics Workbench **ST002195**: dorsal skin from *Tmem79* (Spade) dermatitis-model
mice vs WT at P0, 4, 8 and 10 weeks. It is reverse-phase LC-MS with MS-DIAL annotations and 830 parsed lipid species.

- **Activity proxies.** Each proxy is a log2 product/substrate ratio for one enzymatic step, 17 in total.
  Examples: Cer[NS]/SM for sphingomyelinase vs SM synthase, ≥C26 vs ≤C24 acyl Cer[NS] for CERS3/ELOVL4,
  (16:1+18:1)/(16:0+18:0) acyl chains for SCD, and LPC/PC for PLA2 vs LPCAT. Each ratio is fit as
  `ratio ~ genotype + age` over 4–10 weeks, and the result is projected onto the enzymes that raise or lower it.
- **LINEX2.** The reaction-network enrichment runs on the same data as a cross-check. It independently
  selects the PC/PE ↔ LPC/LPE deacylation network (PLA2G4C, ASPG, TAZ).

**Integration** (`rank.py`)

Score = 0.40 × AD transcript + 0.20 × Th2 keratinocyte response + 0.30 × lipidome activity (full weight
only when its direction matches the transcript) + 0.10 × keratinocyte specificity. The score is halved
when the enzyme is barely expressed in skin. Every component is a percentile, and the weights can be changed as arguments.

### Headline results (see [`SUMMARY.md`](results/insilico_pace/SUMMARY.md))

**The method recovers known biology.** The positive controls land at the 78th percentile on average
(one-sided Mann-Whitney vs all other enzymes, p = 0.009):

| Control | Rank / 764 | Why it's a control |
|---|---|---|
| ELOVL3 | 12 | IL-4/IL-13 suppress it via STAT6. It is down 2.1 log2 in AD skin, and the lipidome's ≥C20 FA fraction falls |
| SCD | 60 | Up in AD skin |
| ELOVL6 | 105 | Suppressed by IL-4/IL-13, and down in AD differentiated keratinocytes |
| ASAH1 | 165 | Its β subunit is the SM deacylase |
| CERS3 / ELOVL4 | 306 / 370 | The lipidome shows the strongest ULC-ceramide loss (q < 0.001), but transcripts are flat in bulk AD. A case where activity changes without an mRNA change |

**The lipidome shows the AD chain-shortening signature.** In Spade vs WT skin, VLC and ULC Cer[NS]
fall about 1.9 log2 by 10 weeks. FA ≥C20 falls, SM → Cer drops, and α-hydroxy, phyto and GlcCer
fractions rise. The effects are small at P0 and grow as dermatitis develops ([figure](results/insilico_pace/fig_lipid_proxies.png)).

**Top candidates.** The phospholipase A2 family ranks first. PLA2G2F, PLA2G4D and PLA2G4E are
suprabasal-keratinocyte enzymes that are up in AD lesions, and PC deacylation rises in the lipidome.
Below them come DGAT2, AADACL2, ALOXE3, DEGS2, FA2H and ELOVL3.

**Serine hydrolase predictions.** [`pace_serine_hydrolase_predictions.csv`](results/insilico_pace/pace_serine_hydrolase_predictions.csv)
lists the skin-expressed FP-reactive enzymes in ranked order. The ones I'd expect to change most are
PLA2G4D, PLA2G4E, PLA2G4F, AADACL2, PNLIPRP3, EPHX3, LIPA, PNPLA1 and ABHD12.

### Caveats

- **The lipidome is from mouse.** The obvious human option, the 2026 Zenodo stratum corneum set
  (record 20365907), is polar metabolomics: only 26 of its 195 features are lipids, and it has no ceramide
  species, so LINEX2 cannot use it. The human AD biopsy lipidome ST003291 is 91% missing values. I'll switch to a
  human SC lipidome once a usable one is available. `lipidome.load_spade()` is the only function that has to change.
- **Transcripts are not activity.** The ranking is a prior to test against, not a result.
- **Many-to-one proxies.** LPC/PC can't tell apart the 11 PLA2s it maps to; the transcript and cell-type layers break the tie.
- **Small cytokine arm.** The cytokine dataset has n = 3 per arm. IL-4 → ELOVL3 is −1.3 log2 but p = 0.13.
- **Atlas not used yet.** The Human Skin Cell Atlas objects (Zenodo 21022952; `ife_extended.h5ad` is 3.4 GB)
  aren't downloaded. HPA skin clusters stand in for cell-type localization. The atlas IFE object would add
  keratinocyte layer resolution. See *Next steps*.

### Next steps

1. Compare `pace_serine_hydrolase_predictions.csv` against measured enzyme activity data.
2. Add HSCA `ife_extended.h5ad` pseudobulk by keratinocyte layer (read in backed mode; about 3.4 GB on disk).
3. Replace the mouse lipidome with a human AD stratum corneum lipidome.
4. Add the IL-4/STAT6 differentiating-keratinocyte RNA-seq (GSE58749) once I've confirmed the accession.

---

## Project 2: lipid annotation QC with retention-time prediction

`lipid-rt-qc` reads an MS-DIAL alignment export, a Skyline report or a MetaboLights MAF, and flags
annotations that don't fit the run's own chromatography. It also automates blank filtering and QC-CV checks.

```bash
lipid-rt-qc run Height_0_2026.txt --out qc/            # MS-DIAL alignment export
lipid-rt-qc run skyline_report.csv --blank "blank|PBS"  # Skyline long-format report
lipid-rt-qc benchmark Height_0_2026.txt --out bench/   # how well RT QC works on *your* data
```

The output is `lipid_qc_report.html`, a self-contained page with figures that opens offline. It comes with
`lipid_qc_flagged.csv`, which lists each flagged feature and why, and `lipid_qc_all_features.csv`.

### Checks

| Flag | What it means |
|---|---|
| **RT off-model** | The observed RT is more than 3 robust SD (and more than 0.2 min) from the RT predicted for that class, C and DB. The model is trained on the same run, and each feature is scored out-of-fold, so a bad annotation can't vouch for itself |
| **Elution order** | Reverse-phase ordering is broken: more carbons should elute later and more double bonds earlier. Cer[NS] should elute before its DHCer. Each violating pair is blamed on the member that fits worse |
| **Isobar fits RT better** | Another class or adduct has the same m/z within 5 ppm, and it matches the RT when the current annotation doesn't. Examples: PC x:y [M+H]⁺ vs PC (x−2):(y−3) [M+Na]⁺, and Cer [M+HAc−H]⁻ vs hydroxy-Cer [M−H]⁻ |
| **BCFA isomer candidate** | The same molecular species also appears 0.3–1 carbon-equivalents earlier. That is the elution signature of monomethyl branched-chain (iso/anteiso) acyl isomers. Ether O-/P- pairs, labelled standards and split peaks get their own labels |
| **Blank / QC CV** | Sample/blank ratio below 3, or pooled-QC CV above 30%. Both thresholds are settable |

**RT model.** For each class it fits a robust (Huber) linear trend, RT ~ C + DB, which is the
equivalent-carbon-number rule. Gradient boosting with an absolute-error loss then corrects the residuals.
With boosting applied to residuals rather than to RT directly, predictions stay sensible at the edges of a
class's carbon range, and one mis-annotated feature can't drag its neighbours. Oxidized glycerolipids and
ether TGs get their own trend lines. Classes the run has never seen are reported but not scored.

### Validation on public data

I validated it on three public reverse-phase lipidomics runs, all annotated with MS-DIAL. MTBLS5861 includes skin.
Reports are in [`results/lipid_rt_qc`](results/lipid_rt_qc).

| Run | Samples | Lipids | Classes | Median \|RT error\| (min) | R² | RT off-model | Isobar better | BCFA cand. |
|---|---|---|---|---|---|---|---|---|
| MTBLS5861 | 7 tissues incl. skin, CID+, MS-DIAL 5 | 2285 | 44 | 0.11 | 0.966 | 229 | 50 | 130 |
| MTBLS14841 | mouse liver, pos, MS-DIAL 5 | 895 | 32 | 0.07 | 0.997 | 66 | 13 | 0 |
| MTBLS14692 | C2C12 Scd1/Scd2 KO, alternating, species level | 998 | 28 | 0.22 | 0.931 | 56 | 20 | 5 |

For the benchmark, 5% of annotations in each run were deliberately corrupted, three seeds per error type
([`evaluate.py`](src/lipid_rt_qc/evaluate.py)). The score is |robust z| of the out-of-fold RT residual.

| Injected error | MTBLS5861 AUC / recall | MTBLS14841 AUC / recall | MTBLS14692 AUC / recall |
|---|---|---|---|
| Carbons ±2 | 0.91 / 68% | 0.97 / 91% | 0.83 / 33% |
| Double bond ±1 | 0.83 / 38% | 0.93 / 77% | 0.78 / 18% |
| PC ↔ PE+3C (co-eluting isobar) | 0.64 / 19% | 0.77 / 37% | 0.65 / 2% |
| [M+H]⁺ read as [M+Na]⁺ | 0.86 / 38% | 0.97 / 94% | 0.85 / 20% |
| Wrong peak (RT swapped) | 0.98 / 94% | 0.99 / 98% | 0.95 / 78% |

Flag rate on untouched features at the default threshold: MTBLS5861: 9% · MTBLS14841: 6% · MTBLS14692: 5%.

Power depends on the run. MTBLS14692 annotates only sum compositions, and its RT scatter is twice as wide, so recall there is lower at the same false-positive rate. Use `lipid-rt-qc benchmark` to see what to expect on a given dataset.

RT alone reliably catches wrong-peak picks and chain-length or unsaturation mislabels,
and it catches most sodium-adduct mix-ups. It can't separate PC x:y from PE (x+3):y, because they co-elute;
that's a chemistry limit rather than a model limit, so the report marks those as needing MS/MS head-group fragments. "Clean"
features in public data include the depositors' own mis-annotations, so the false-positive rate is an upper bound.

### Running it on your own data

1. In MS-DIAL, export the alignment result (Height or Area) as a `.txt` file, then run `lipid-rt-qc run export.txt`.
   Blank and QC samples are picked up from MS-DIAL's *File type* row.
2. From Skyline, export a report with Molecule Name, Precursor Adduct, Precursor Mz, Replicate Name,
   Retention Time and Area. Name blanks and QCs so `--blank` / `--qc` can match them.
3. Work through `lipid_qc_flagged.csv` from the top, and check MS/MS in El-MAVEN for anything flagged as RT off-model or isobar.

---

## Data sources and licenses

All inputs are public. Raw files live in `data/raw/`, which is git-ignored and re-created by
`insilico-pace download`. The derived tables in `results/` are small and are committed.

- UniProt REST (CC-BY 4.0) · Rhea (CC-BY 4.0)
- GEO GSE121212, GSE36287, GSE53795
- Human Protein Atlas single-cell clusters (CC-BY-SA 3.0)
- Metabolomics Workbench ST002195 and ST003291 (CC-BY 4.0)
- LINEX2 (Köhler et al., *Brief Bioinform* 2023)
