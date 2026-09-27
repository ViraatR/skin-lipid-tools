# skin-lipid-tools

> **Private.** Keep this repo private until Dr. Kolar has seen it. It uses only public data, but it
> predicts targets for work the lab has not yet published.

This repo holds two computational projects for skin lipid biology. Both are built on public data.

| | Project | Question | Status |
|---|---|---|---|
| 1 | [`insilico_pace`](src/insilico_pace) | Using only public data, which skin enzymes should PACE find dysregulated in atopic dermatitis? | Pipeline runs end to end; results are in [`results/insilico_pace`](results/insilico_pace/SUMMARY.md) |
| 2 | [`lipid_rt_qc`](src/lipid_rt_qc) | Which lipid annotations in an MS-DIAL or Skyline export look suspicious? | See [Project 2](#project-2-lipid-annotation-qc-with-retention-time-prediction) |

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
This project is the computational counterpart. It builds a catalog of the enzymes that could matter,
scores each one against every public dataset that says something about AD, and produces a ranked list.
That list gives the IL-4/IL-13 keratinocyte ABPP experiment something to test against.

### Layers

**Layer 1: catalog and expression** (`catalog.py`, `transcriptome.py`)

| Source | What it adds |
|---|---|
| UniProt (reviewed human; lipid metabolism or hydrolase EC classes) | 764 enzymes, Rhea reactions, ChEMBL targets. **99 FP-reactive serine hydrolases**, called from a catalytic Ser plus serine-hydrolase InterPro folds |
| GSE121212 (Tsoi 2019): AD, psoriasis and healthy skin RNA-seq | AD lesional vs healthy (21 vs 38), AD non-lesional vs healthy, paired AD LS vs NL, psoriasis vs healthy |
| GSE36287: primary keratinocytes + IL-4, IL-13, IL-17A, IFN-α, IFN-γ, TNF | Th2 and other cytokine responses (n=3 per arm, paired by donor) |
| GSE53795: acne lesional vs non-lesional | Side contrast. This is the dataset from lab Slack; it turned out to be acne, not AD |
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
| SCD | 60 | Up in AD skin, matching the lab's scRNA reanalysis |
| ELOVL6 | 105 | Suppressed by IL-4/IL-13, and down in AD differentiated keratinocytes |
| ASAH1 | 165 | Its β subunit is the SM deacylase |
| CERS3 / ELOVL4 | 306 / 370 | The lipidome shows the strongest ULC-ceramide loss (q < 0.001), but transcripts are flat in bulk AD. This is the "activity changes without mRNA change" case that PACE exists to measure |

**The lipidome shows the AD chain-shortening signature.** In Spade vs WT skin, VLC and ULC Cer[NS]
fall about 1.9 log2 by 10 weeks. FA ≥C20 falls, SM → Cer drops, and α-hydroxy, phyto and GlcCer
fractions rise. The effects are small at P0 and grow as dermatitis develops ([figure](results/insilico_pace/fig_lipid_proxies.png)).

**Top candidates.** The phospholipase A2 family ranks first. PLA2G2F, PLA2G4D and PLA2G4E are
suprabasal-keratinocyte enzymes that are up in AD lesions, and PC deacylation rises in the lipidome.
Below them come DGAT2, AADACL2, ALOXE3, DEGS2, FA2H and ELOVL3.

**PACE-visible predictions.** [`pace_serine_hydrolase_predictions.csv`](results/insilico_pace/pace_serine_hydrolase_predictions.csv)
lists the skin-expressed FP-reactive enzymes in ranked order. Those most likely to differ in the
IL-4/IL-13 ABPP gel/MS are PLA2G4D, PLA2G4E, PLA2G4F, AADACL2, PNLIPRP3, EPHX3, LIPA, PNPLA1 and ABHD12.

### Caveats

- **The lipidome is from mouse.** The obvious human option, the 2026 Zenodo stratum corneum set
  (record 20365907), is polar metabolomics: only 26 of its 195 features are lipids, and it has no ceramide
  species, so LINEX2 cannot use it. The human AD biopsy lipidome ST003291 is 91% missing values. Swap in a human
  SC lipidome as soon as one is available, or use the lab's own. `lipidome.load_spade()` is the only function to replace.
- **Transcripts are not activity.** That gap is the reason for PACE. The ranking is a prior to test against, not a result.
- **Many-to-one proxies.** LPC/PC can't tell apart the 11 PLA2s it maps to; the transcript and cell-type layers break the tie.
- **Small cytokine arm.** The cytokine dataset has n = 3 per arm. IL-4 → ELOVL3 is −1.3 log2 but p = 0.13.
- **Atlas not used yet.** The Human Skin Cell Atlas objects (Zenodo 21022952; `ife_extended.h5ad` is 3.4 GB)
  aren't downloaded. HPA skin clusters stand in for cell-type localization. The atlas IFE object would add
  keratinocyte layer resolution. See *Next steps*.

### Next steps

1. When the ABPP hits come back, compare them with `pace_serine_hydrolase_predictions.csv`. That comparison is the "predicted, then measured" figure.
2. Add HSCA `ife_extended.h5ad` pseudobulk by keratinocyte layer (read in backed mode; about 3.4 GB on disk).
3. Replace the mouse lipidome with a human AD stratum corneum lipidome.
4. Add the IL-4/STAT6 differentiating-keratinocyte RNA-seq (GSE58749) once its accession is confirmed.
   The link in the original notes points to E-GEOD-21413.

---

## Project 2: lipid annotation QC with retention-time prediction

See [`src/lipid_rt_qc`](src/lipid_rt_qc) and the section below once built.

---

## Data sources and licenses

All inputs are public. Raw files live in `data/raw/`, which is git-ignored and re-created by
`insilico-pace download`. The derived tables in `results/` are small and are committed.

- UniProt REST (CC-BY 4.0) · Rhea (CC-BY 4.0)
- GEO GSE121212, GSE36287, GSE53795
- Human Protein Atlas single-cell clusters (CC-BY-SA 3.0)
- Metabolomics Workbench ST002195 and ST003291 (CC-BY 4.0)
- LINEX2 (Köhler et al., *Brief Bioinform* 2023)
