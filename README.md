# scATAC mouse digit tip project - using snapATAC2

```text
 __  __                        ____  _       _ _     _____ _       
|  \/  | ___  _   _ ___  ___  |  _ \(_) __ _(_) |_  |_   _(_)_ __  
| |\/| |/ _ \| | | / __|/ _ \ | | | | |/ _` | | __|   | | | | '_ \ 
| |  | | (_) | |_| \__ \  __/ | |_| | | (_| | | |_    | | | | |_) |
|_|  |_|\___/ \__,_|___/\___| |____/|_|\__, |_|\__|   |_| |_| .__/ 
                                       |___/                |_|
``` 

## Quality Control (QC) Filtering

Following generation of the initial QC metrics, cells are filtered to remove low-quality nuclei prior to downstream analysis. A cell is retained **only if it satisfies all of the following criteria**:

> ##### `n_fragment`
>
> The total number of fragments assigned to each cell barcode, commonly used as a measure of sequencing depth per cell.
>
> - **Low (<1,000):** Typically indicates empty droplets, debris, or low-quality cells with insufficient data.
> - **High (>50,000):** May indicate doublets or multiplets, but can also represent genuine high-quality cells. These cells should be inspected rather than automatically removed.
>
> ##### `frac_dup`
>
> The fraction of fragments that are PCR duplicates (identical genomic start/end positions), used as a measure of library complexity.
>
> - **Low (<0.1):** Indicates high library complexity with many unique fragments. Generally considered good quality.
> - **High (>0.5):** Indicates low library complexity, where many fragments are duplicated PCR products. These cells contain less unique information and are typically removed.
>
> ##### `TSSe`
>
> Transcription Start Site enrichment score: a per-cell measure of signal-to-noise ratio, comparing fragment enrichment around transcription start sites (TSSs) versus surrounding regions.
>
> - **High (>2.0):** Indicates good signal quality, with fragments enriched around promoters. These cells generally have higher-quality chromatin accessibility data.
> - **Low (<1.5):** Indicates poor signal quality, with weak TSS enrichment and more background noise. These cells are typically removed.



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


> #### Note
>
> ##### `frac_dup` vs `n_fragment` Scatter Plot
>
> Shows how PCR duplication changes with sequencing depth per cell.
>
> **Purpose:**
>
> - **Detect low-complexity libraries:** Cells with high `frac_dup` (>0.5), especially at moderate `n_fragment`, contain many duplicated fragments and have reduced unique information.
> - **Identify potential doublets:** Cells with very high `n_fragment` (>50k) together with elevated `frac_dup` may represent multiplets, although they require additional confirmation.
> - **Assess data quality:** High-quality cells should maintain relatively low `frac_dup` across different fragment depths, indicating a higher proportion of unique fragments.
>
> ##### `TSSe` vs `n_fragment` Scatter Plot
>
> Shows how TSS enrichment changes with sequencing depth per cell.
>
> **Purpose:**
>
> - **Detect poor-quality cells:** Cells with high `n_fragment` but low `TSSe` (<1.5) contain many fragments but weak promoter enrichment, indicating high background noise or poor-quality nuclei. These cells are typically removed.
> - **Assess sequencing saturation:** Good cells often show increasing `TSSe` with increasing fragment number until reaching a plateau.
> - **Identify high-quality cells:** Cells with high `n_fragment` and high `TSSe` (top-right region) generally represent deeply sequenced, high-quality libraries and are usually retained.

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


### Doublet Removal 

![](figures/noDoublets_before_doublet_removal_doublet_probability.png?v=2) 

![](figures/noDoublets_before_doublet_removal_doublet_probability_vs_n_fragment.png?v=2)

A total of **179,807 cells** were initially loaded. Doublets were identified using a **doublet score cutoff of 0.5**. After filtering, **179,258 cells** remained for downstream analysis.

![](figures/noDoublets_after_doublet_removal_doublet_probability_vs_n_fragment.png?v=1)


## Chromatin Accessibility Profiles and Clustering 

- Fragments are pieces of DNA generated from regions of the genome that were **open and accessible** in a cell.
- Each fragment contains a **cell barcode**, allowing us to assign it back to the cell it came from.
- For each cell, we measure **where fragments are located across the genome**.
- This creates a **chromatin accessibility profile**: a pattern of genomic regions that are accessible in that cell.
- Cells with **similar accessibility profiles** are grouped together computationally.
- Similar accessibility patterns suggest that cells share similar **regulatory programs**, which can reflect the same cell type or biological state.
- Clustering is a method to **identify these groups automatically** based on similarities in chromatin accessibility patterns.

<img src="figures/analysis_umap_clusters.png?v=1" width="45%" /><img src="figures/analysis_umap_sample.png?v=2" width="45%" />

### Check QC per cluster leiden 

![](figures/analysis_qc_violin_by_leiden.png?v=1)

<img src="figures/analysis_umap_n_fragment.png?v=1" width="30%" /><img src="figures/analysis_umap_tsse.png?v=1" width="30%" /><img src="figures/analysis_umap_frac_dup.png?v=1" width="30%" />


### Differential Peak analyis and annotations 

![](figures/3xAmp_vs_1xAmp_volcano_annotated.png?v=1) 

![](figures/5xAmp_vs_Control_volcano_annotated.png?v=1)

![](figures/5xAmp_vs_1xAmp_volcano_annotated.png?v=1)


## References 

SnapATAC2: Zhang, K., Zemke, N. R., Armand, E. J., & Ren, B. (2024). A fast, scalable and versatile tool for analysis of single-cell omics data. Nature methods, 21(2), 217-227.

