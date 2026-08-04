# scATAC mouse digit tip project

# Quality Control (QC) Filtering

Following generation of the initial QC metrics, cells are filtered to remove low-quality nuclei prior to downstream analysis. A cell is retained **only if it satisfies all of the following criteria**:

```
##### `n_fragment`

The total number of fragments assigned to each cell barcode, commonly used as a measure of sequencing depth per cell.

- **Low (<1,000):** Typically indicates empty droplets, debris, or low-quality cells with insufficient data.
- **High (>50,000):** May indicate doublets or multiplets, but can also represent genuine high-quality cells. These cells should be inspected rather than automatically removed.


##### `frac_dup`

The fraction of fragments that are PCR duplicates (identical genomic start/end positions), used as a measure of library complexity.

- **Low (<0.1):** Indicates high library complexity with many unique fragments. Generally considered good quality.
- **High (>0.5):** Indicates low library complexity, where many fragments are duplicated PCR products. These cells contain less unique information and are typically removed.


##### `TSSe`

Transcription Start Site enrichment score: a per-cell measure of signal-to-noise ratio, comparing fragment enrichment around transcription start sites (TSSs) versus surrounding regions.

- **High (>2.0):** Indicates good signal quality, with fragments enriched around promoters. These cells generally have higher-quality chromatin accessibility data.
- **Low (<1.5):** Indicates poor signal quality, with weak TSS enrichment and more background noise. These cells are typically removed.
```


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
#### Stats
| Sample | Cells Before Filtering | Cells After Filtering | Cells Retained (%) |
|--------|-----------------------:|----------------------:|-------------------:|
| **SINAA6** | 83,873 | 14,007 | 16.70% |
| **SINAB6** | 75,844 | 63,199 | 83.33% |
| **SINAA8** | 79,172 | 62,086 | 78.42% |
| **SINAH5** | 83,834 | 49,410 | 58.94% |


#### Note 
##### `frac_dup` vs `n_fragment` Scatter Plot

Shows how PCR duplication changes with sequencing depth per cell.

**Purpose:**

- **Detect low-complexity libraries:** Cells with high `frac_dup` (>0.5), especially at moderate `n_fragment`, contain many duplicated fragments and have reduced unique information.
- **Identify potential doublets:** Cells with very high `n_fragment` (>50k) together with elevated `frac_dup` may represent multiplets, although they require additional confirmation.
- **Assess data quality:** High-quality cells should maintain relatively low `frac_dup` across different fragment depths, indicating a higher proportion of unique fragments.


##### `TSSe` vs `n_fragment` Scatter Plot

Shows how TSS enrichment changes with sequencing depth per cell.

**Purpose:**

- **Detect poor-quality cells:** Cells with high `n_fragment` but low `TSSe` (<1.5) contain many fragments but weak promoter enrichment, indicating high background noise or poor-quality nuclei. These cells are typically removed.
- **Assess sequencing saturation:** Good cells often show increasing `TSSe` with increasing fragment number until reaching a plateau.
- **Identify high-quality cells:** Cells with high `n_fragment` and high `TSSe` (top-right region) generally represent deeply sequenced, high-quality libraries and are usually retained.

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




## References 

SnapATAC2: Zhang, K., Zemke, N. R., Armand, E. J., & Ren, B. (2024). A fast, scalable and versatile tool for analysis of single-cell omics data. Nature methods, 21(2), 217-227.

