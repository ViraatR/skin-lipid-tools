# In silico PACE: results summary

Enzymes scored: 764 (99 FP-reactive serine hydrolases).

## Positive controls
Mean percentile 78; one-sided Mann-Whitney vs all other enzymes p = 0.0089.

| gene   |   rank |   score |   AD_LS_vs_HC_lfc |   th2_lfc | lipid_proxy    |   activity_t |   percentile |
|:-------|-------:|--------:|------------------:|----------:|:---------------|-------------:|-------------:|
| ELOVL3 |     12 |   0.606 |            -2.14  |    -0.973 | FA_VLC         |        -3.27 |         98.6 |
| SCD    |     60 |   0.486 |             0.841 |     0.069 | SCD_index      |        -3.91 |         92.3 |
| ELOVL6 |    105 |   0.435 |            -0.727 |    -0.028 | C18_C16_chains |        -2.71 |         86.4 |
| ASAH1  |    165 |   0.392 |            -0.69  |     0.049 | nan            |       nan    |         78.5 |
| CERS3  |    306 |   0.285 |             0.031 |    -0.171 | ULC_ceramide   |        -5.05 |         60.1 |
| ELOVL4 |    370 |   0.243 |             0.003 |    -0.391 | ULC_ceramide   |        -5.05 |         51.7 |

## Top 25 candidates
| gene     |   rank |   score | serine_hydrolase   | ad_direction   |   AD_LS_vs_HC_lfc |   th2_lfc | lipid_proxy         |   activity_t | top_skin_cell_type         |
|:---------|-------:|--------:|:-------------------|:---------------|------------------:|----------:|:--------------------|-------------:|:---------------------------|
| PLA2G2F  |      1 |   0.688 | False              | up             |             1.74  |    -0.029 | PC_deacylation      |         4.79 | suprabasal keratinocytes   |
| PLA2G4D  |      2 |   0.673 | True               | up             |             1.27  |    -0.121 | PC_deacylation      |         4.79 | suprabasal keratinocytes   |
| DGAT2    |      3 |   0.651 | False              | down           |            -0.983 |    -0.494 | TG_lipolysis        |        -3.32 | cdc                        |
| PLA2G4E  |      4 |   0.64  | True               | up             |             1.11  |   nan     | PC_deacylation      |         4.79 | suprabasal keratinocytes   |
| AADACL2  |      5 |   0.636 | True               | down           |            -0.7   |    -0.695 | nan                 |       nan    | suprabasal keratinocytes   |
| ALOXE3   |      6 |   0.635 | False              | down           |            -0.852 |    -0.506 | Acylceramide        |         3.34 | suprabasal keratinocytes   |
| DEGS2    |      7 |   0.625 | False              | down           |            -0.26  |    -0.519 | Phyto_hydroxylation |        16    | suprabasal keratinocytes   |
| PSAPL1   |      8 |   0.62  | False              | down           |            -0.622 |    -0.765 | nan                 |       nan    | suprabasal keratinocytes   |
| PCSK9    |      9 |   0.617 | True               | up             |             1.44  |     0.505 | nan                 |       nan    | fibroblasts                |
| CYP39A1  |     10 |   0.613 | False              | down           |            -1.49  |    -0.384 | nan                 |       nan    | suprabasal keratinocytes   |
| FA2H     |     11 |   0.61  | False              | down           |            -0.81  |    -0.333 | Alpha_hydroxylation |         7.82 | suprabasal keratinocytes   |
| ELOVL3   |     12 |   0.606 | False              | down           |            -2.14  |    -0.973 | FA_VLC              |        -3.27 | vascular endothelial cells |
| SCD5     |     13 |   0.593 | False              | down           |            -1.62  |     0.194 | SCD_index           |        -3.91 | smooth muscle cells        |
| CYP27B1  |     14 |   0.59  | False              | up             |             0.985 |     0.566 | nan                 |       nan    | basal keratinocytes        |
| ACER1    |     15 |   0.585 | False              | down           |            -0.627 |    -0.344 | nan                 |       nan    | suprabasal keratinocytes   |
| PLA2G4F  |     16 |   0.579 | True               | up             |             0.242 |     0.188 | PC_deacylation      |         4.79 | suprabasal keratinocytes   |
| PNLIPRP3 |     17 |   0.573 | True               | down           |            -1.02  |    -0.626 | nan                 |       nan    | basal keratinocytes        |
| PLA2G3   |     18 |   0.571 | False              | up             |             1.55  |    -0.407 | nan                 |       nan    | suprabasal keratinocytes   |
| HSD11B1  |     19 |   0.568 | False              | down           |            -3.07  |    -0.23  | nan                 |       nan    | fibroblasts                |
| LPCAT3   |     20 |   0.565 | False              | up             |             0.363 |     0.346 | PC_deacylation      |        -4.79 | pericytes                  |
| ALDH3A2  |     21 |   0.563 | False              | down           |            -0.916 |    -0.223 | nan                 |       nan    | suprabasal keratinocytes   |
| CYP7B1   |     22 |   0.56  | False              | up             |             0.586 |     0.52  | nan                 |       nan    | fibroblasts                |
| CYP2W1   |     23 |   0.559 | False              | down           |            -1.66  |    -0.127 | nan                 |       nan    | suprabasal keratinocytes   |
| DEGS1    |     24 |   0.559 | False              | down           |            -0.518 |     0.104 | DHCer_desaturation  |        23.8  | suprabasal keratinocytes   |
| SPTLC3   |     25 |   0.559 | False              | down           |            -0.801 |    -0.241 | nan                 |       nan    | suprabasal keratinocytes   |

## Serine hydrolase predictions
Skin-expressed FP-reactive enzymes, ranked. `th2_lfc` is the mean IL-4/IL-13 transcript change in primary keratinocytes (GSE36287). Transcript changes are a prior, not a measure of enzyme activity.
| gene     |   rank |   score | serine_hydrolase   | ad_direction   |   AD_LS_vs_HC_lfc |   th2_lfc | lipid_proxy    |   activity_t | top_skin_cell_type           |   KC_nCPM |
|:---------|-------:|--------:|:-------------------|:---------------|------------------:|----------:|:---------------|-------------:|:-----------------------------|----------:|
| PLA2G4D  |      2 |   0.673 | True               | up             |             1.27  |    -0.121 | PC_deacylation |         4.79 | suprabasal keratinocytes     |       5.8 |
| PLA2G4E  |      4 |   0.64  | True               | up             |             1.11  |   nan     | PC_deacylation |         4.79 | suprabasal keratinocytes     |     150   |
| AADACL2  |      5 |   0.636 | True               | down           |            -0.7   |    -0.695 | nan            |       nan    | suprabasal keratinocytes     |       3.7 |
| PCSK9    |      9 |   0.617 | True               | up             |             1.44  |     0.505 | nan            |       nan    | fibroblasts                  |       1.1 |
| PLA2G4F  |     16 |   0.579 | True               | up             |             0.242 |     0.188 | PC_deacylation |         4.79 | suprabasal keratinocytes     |      22.6 |
| PNLIPRP3 |     17 |   0.573 | True               | down           |            -1.02  |    -0.626 | nan            |       nan    | basal keratinocytes          |       5.1 |
| LIPA     |     28 |   0.546 | True               | up             |             0.808 |     0.292 | nan            |       nan    | vascular endothelial cells   |      16.7 |
| EPHX3    |     36 |   0.527 | True               | up             |             0.638 |    -0.436 | nan            |       nan    | suprabasal keratinocytes     |      75.3 |
| LIPH     |     45 |   0.506 | True               | down           |            -1.3   |    -0.111 | nan            |       nan    | melanocytes                  |       3.5 |
| NCEH1    |     62 |   0.483 | True               | up             |             0.507 |     0.349 | nan            |       nan    | fibroblasts                  |      19.2 |
| PNPLA5   |     74 |   0.464 | True               | down           |            -1.18  |    -0.182 | nan            |       nan    | basal keratinocytes          |       0   |
| LIPG     |     82 |   0.453 | True               | up             |             1.18  |    -0.164 | nan            |       nan    | suprabasal keratinocytes     |       6.3 |
| PAFAH2   |     84 |   0.45  | True               | down           |            -0.502 |    -0.11  | nan            |       nan    | melanocytes                  |       3.7 |
| EPHX2    |     87 |   0.447 | True               | down           |            -0.704 |     0.283 | nan            |       nan    | vascular smooth muscle cells |      15.8 |
| PNPLA3   |     88 |   0.447 | True               | down           |            -1.17  |     0.168 | nan            |       nan    | fibroblasts                  |       1.3 |
| PNPLA1   |     92 |   0.444 | True               | down           |            -0.227 |    -0.442 | Acylceramide   |         3.34 | suprabasal keratinocytes     |       7.9 |
| ABHD12   |     97 |   0.438 | True               | down           |            -0.34  |    -0.342 | nan            |       nan    | cdc                          |      25.3 |
| PLA2G7   |    110 |   0.427 | True               | down           |            -0.639 |    -0.21  | nan            |       nan    | macrophages                  |       1.6 |
| PAFAH1B3 |    120 |   0.42  | True               | up             |             0.399 |     0.109 | nan            |       nan    | melanocytes                  |      33   |
| FAAH2    |    125 |   0.415 | True               | up             |             0.862 |     0.017 | nan            |       nan    | plasma cells                 |      53.2 |

## Lipidome activity proxies (Spade vs WT, age-adjusted, 4-10 wk)
| proxy               |   lfc_spade |     t |     q |   lfc_P0 |   lfc_4w |   lfc_8w |   lfc_10w | note                                                                               |
|:--------------------|------------:|------:|------:|---------:|---------:|---------:|----------:|:-----------------------------------------------------------------------------------|
| DHCer_desaturation  |       2.18  | 23.8  | 0     |   -0.06  |    2.11  |    2.2   |     2.21  | Cer[NS]/Cer[NDS]: sphinganine -> sphingosine backbone                              |
| Phyto_hydroxylation |       2.44  | 16    | 0     |    0.441 |    1.95  |    2.41  |     2.96  | Cer[NP]/Cer[NDS]: C4-hydroxylation to phytoceramide                                |
| Alpha_hydroxylation |       0.774 |  7.82 | 0     |    0.062 |    0.247 |    0.918 |     1.11  | alpha-OH / non-OH ceramides                                                        |
| SM_hydrolysis       |      -0.613 | -6.71 | 0     |   -0.176 |   -0.341 |   -0.647 |    -0.84  | Cer[NS]/SM: sphingomyelinase vs SM synthase                                        |
| Cer_glucosylation   |       0.737 |  5.75 | 0     |    0.207 |    0.267 |    0.705 |     1.25  | GlcCer/Cer: the SC ceramide pool is released from GlcCer by GBA1                   |
| ULC_ceramide        |      -0.92  | -5.05 | 0     |   -0.27  |   -0.189 |   -0.747 |    -1.88  | acyl >=C26 vs <=C24 in Cer[NS]                                                     |
| PC_deacylation      |       0.392 |  4.79 | 0     |    0.051 |    0.131 |    0.281 |     0.801 | LPC/PC                                                                             |
| VLC_ceramide        |      -0.898 | -4.69 | 0.001 |   -0.105 |   -0.321 |   -0.602 |    -1.87  | C22-24 vs C16-18 Cer[NS]                                                           |
| SCD_index           |      -0.249 | -3.91 | 0.002 |    0.185 |   -0.006 |   -0.181 |    -0.582 | (16:1+18:1)/(16:0+18:0) acyl chains in glycerolipids                               |
| PE_methylation      |       0.202 |  3.85 | 0.002 |    0.141 |   -0.031 |    0.2   |     0.438 | PC/PE                                                                              |
| PS_decarboxylation  |      -0.214 | -3.75 | 0.003 |   -0.125 |    0.02  |   -0.187 |    -0.486 | PE/PS                                                                              |
| PE_deacylation      |       0.546 |  3.46 | 0.005 |    0.04  |   -0.069 |    0.397 |     1.36  | LPE/PE                                                                             |
| TG_lipolysis        |       0.549 |  3.32 | 0.005 |   -0.424 |    0.001 |    0.482 |     1.19  | DG/TG                                                                              |
| Acylceramide        |       0.487 |  3.34 | 0.005 |    0.036 |    0.755 |    0.771 |    -0.16  | Cer[EOS]/ultra-long Cer[NS]: omega-O-acylceramide pathway                          |
| FA_VLC              |      -0.822 | -3.27 | 0.005 |   -0.092 |   -0.114 |   -0.66  |    -1.75  | >=C20 vs C14-18 free fatty acids                                                   |
| C18_C16_chains      |      -0.067 | -2.71 | 0.017 |    0.315 |   -0.058 |   -0.02  |    -0.137 | 18:0+18:1 vs 16:0+16:1 acyl chains in glycerolipids (free C16/C18 FA not reported) |
| Plasmalogen         |       0.075 |  2.55 | 0.021 |    0.237 |    0.019 |    0.026 |     0.197 | PE-P/PE                                                                            |

## LINEX2 enriched subnetwork enzymes
| gene    |   linex_edges |   linex_p | linex_reactions   |
|:--------|--------------:|----------:|:------------------|
| PLA2G4C |            26 |         0 | L_FAdelete        |
| ASPG    |            22 |         0 | L_FAdelete        |
| TAZ     |             4 |         0 | L_FAdelete        |
| GPC1    |             1 |         0 | L_HGmodify        |
