# scATAC mosue digit tip project

# Quality Control (QC) Filtering

Following generation of the initial QC metrics, cells are filtered to remove low-quality nuclei prior to downstream analysis. A cell is retained **only if it satisfies all of the following criteria**:

| Metric | Threshold | Purpose |
|--------|----------:|---------|
| **Number of fragments (`n_fragment`)** | ≥ 1,000 and ≤ 100,000 | Removes empty droplets or poorly sequenced cells with insufficient fragments, while excluding cells with unusually high fragment counts that may represent doublets or other technical artifacts. |
| **Transcription Start Site Enrichment (`TSSe`)** | ≥ 1.5 and ≤ 100 | Removes cells with poor enrichment of reads around transcription start sites, indicative of low-quality chromatin accessibility. The upper threshold of 100 is effectively non-restrictive and primarily serves as a safeguard against extreme outliers. |
| **Duplicate Fragment Fraction (`frac_dup`)** | ≤ 0.5 | Excludes cells with excessive PCR duplication, which generally reflects reduced library complexity and poorer data quality. |

The filtering condition applied to each cell is therefore:

```text
Keep cell if:
    1000 ≤ n_fragment ≤ 100000
AND 1.5 ≤ TSSe ≤ 100
AND frac_dup ≤ 0.5
```

#### SINAA6

<img src="figures/SINAA6_prefilter_n_fragment.png" width="45%" /><img src="figures/SINAA6_postfilter_n_fragment.png" width="45%" />

<img src="figures/SINAA6_prefilter_frac_dup.png" width="45%" /><img src="figures/SINAA6_postfilter_frac_dup.png" width="45%" />

<img src="figures/SINAA6_prefilter_tsse.png" width="45%" /><img src="figures/SINAA6_postfilter_tsse.png" width="45%" />

<img src="figures/SINAA6_prefilter_frac_dup_vs_n_fragment.png" width="45%" /><img src="figures/SINAA6_postfilter_frac_dup_vs_n_fragment.png" width="45%" />

<img src="figures/SINAA6_prefilter_tsse_vs_n_fragment.png" width="45%" /><img src="figures/SINAA6_postfilter_tsse_vs_n_fragment.png" width="45%" />


---

#### SINAB6

<img src="figures/SINAB6_prefilter_n_fragment.png" width="45%" /><img src="figures/SINAB6_postfilter_n_fragment.png" width="45%" />

<img src="figures/SINAB6_prefilter_frac_dup.png" width="45%" /><img src="figures/SINAB6_postfilter_frac_dup.png" width="45%" />

<img src="figures/SINAB6_prefilter_tsse.png" width="45%" /><img src="figures/SINAB6_postfilter_tsse.png" width="45%" />

<img src="figures/SINAB6_prefilter_frac_dup_vs_n_fragment.png" width="45%" /><img src="figures/SINAB6_postfilter_frac_dup_vs_n_fragment.png" width="45%" />

<img src="figures/SINAB6_prefilter_tsse_vs_n_fragment.png" width="45%" /><img src="figures/SINAB6_postfilter_tsse_vs_n_fragment.png" width="45%" />


---

#### SINAA8

<img src="figures/SINAA8_prefilter_n_fragment.png" width="45%" /><img src="figures/SINAA8_postfilter_n_fragment.png" width="45%" />

<img src="figures/SINAA8_prefilter_frac_dup.png" width="45%" /><img src="figures/SINAA8_postfilter_frac_dup.png" width="45%" />

<img src="figures/SINAA8_prefilter_tsse.png" width="45%" /><img src="figures/SINAA8_postfilter_tsse.png" width="45%" />

<img src="figures/SINAA8_prefilter_frac_dup_vs_n_fragment.png" width="45%" /><img src="figures/SINAA8_postfilter_frac_dup_vs_n_fragment.png" width="45%" />

<img src="figures/SINAA8_prefilter_tsse_vs_n_fragment.png" width="45%" /><img src="figures/SINAA8_postfilter_tsse_vs_n_fragment.png" width="45%" />


---

#### SINAH5

<img src="figures/SINAH5_prefilter_n_fragment.png" width="45%" /><img src="figures/SINAH5_postfilter_n_fragment.png" width="45%" />

<img src="figures/SINAH5_prefilter_frac_dup.png" width="45%" /><img src="figures/SINAH5_postfilter_frac_dup.png" width="45%" />

<img src="figures/SINAH5_prefilter_tsse.png" width="45%" /><img src="figures/SINAH5_postfilter_tsse.png" width="45%" />

<img src="figures/SINAH5_prefilter_frac_dup_vs_n_fragment.png" width="45%" /><img src="figures/SINAH5_postfilter_frac_dup_vs_n_fragment.png" width="45%" />

<img src="figures/SINAH5_prefilter_tsse_vs_n_fragment.png" width="45%" /><img src="figures/SINAH5_postfilter_tsse_vs_n_fragment.png" width="45%" />



