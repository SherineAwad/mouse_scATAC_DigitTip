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

Following generation of the initial QC metrics, cells were filtered to remove low-quality nuclei prior to downstream analysis. A cell was retained **only if it satisfied all of the following criteria**:

> ##### `n_fragment`
>
> The total number of fragments assigned to each cell barcode, reflecting the amount of sequencing information available for each cell.
>
> * **Low (<1,000):** Cells with very few fragments contain insufficient information for reliable downstream analysis and were removed.
> * **High (>100,000):** Extremely high fragment counts can indicate potential doublets, multiplets, or other technical abnormalities. Cells above this threshold were excluded during QC. Doublets were subsequently assessed separately using a dedicated doublet-detection step.
>
> ##### `frac_dup`
>
> The fraction of fragments that are PCR duplicates, used as an indicator of library complexity.
>
> * **Low:** A lower duplicate fraction indicates a higher proportion of unique fragments and generally better library complexity.
> * **High (>0.5):** A high duplicate fraction indicates reduced library complexity, with a large proportion of fragments being duplicated. These cells were removed.
>
> ##### `TSSe`
>
> Transcription Start Site enrichment score (TSSe), a measure of the enrichment of accessible chromatin signal around transcription start sites relative to surrounding genomic regions.
>
> * **Low (<1.5):** Indicates weak TSS enrichment and lower signal quality; these cells were removed.
> * **Higher values (>2.0):** Indicate stronger enrichment of accessible chromatin around TSSs and generally better signal quality. **Extremely high values (>100) were excluded because unusually concentrated TSS signal can represent an atypical QC profile rather than simply higher-quality data.**

| Metric                                           | Filtering threshold | Purpose                                                                                                                                                                               |
| ------------------------------------------------ | ------------------: | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Number of fragments (`n_fragment`)**           |   **1,000–100,000** | Removes cells with insufficient sequencing information and excludes cells with extremely high fragment counts that may represent potential doublets or other technical abnormalities. |
| **Transcription Start Site Enrichment (`TSSe`)** |         **1.5–100** | Removes cells with weak TSS enrichment while excluding cells with extremely high, potentially atypical TSS enrichment profiles.                                                       |
| **Duplicate Fragment Fraction (`frac_dup`)**     |            **≤0.5** | Removes cells with excessive PCR duplication and reduced library complexity.                                                                                                          |

> > ##### Note: Library complexity
> >
> > **Library complexity** refers to how many **unique DNA fragments** are present in the sequencing library, rather than repeated copies of the same fragments.
> >
> > - **High library complexity** → many different genomic fragments → more unique information.
> > - **Low library complexity** → many fragments are duplicates of fragments already observed → less unique information.

The filtering condition applied to each cell was therefore:

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


> ##### Note
>
> ##### `frac_dup` vs `n_fragment` Scatter Plot
>
> Shows how the PCR duplication fraction changes relative to sequencing depth per cell.
>
> **Purpose:**
>
> * **Detect low-complexity libraries:** Cells with high `frac_dup` (>0.5), especially at moderate `n_fragment` values, contain an excessively high proportion of duplicated fragments and have reduced library complexity.
> * **Visualize sequencing depth cutoffs:** Identifies cells exceeding the upper fragment limit (>100,000), which are removed to mitigate potential multiplets or other technical abnormalities before downstream doublet detection.
> * **Assess overall library quality:** High-quality single-cell libraries generally show a stable, relatively low `frac_dup` across standard fragment depths, indicating a high proportion of unique genomic fragments.



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

Chromatin accessibility was compared between the two selected groups to identify regions with significantly different accessibility. Peaks with an **absolute log2 fold-change > 0.5** and an **adjusted p-value < 0.05** were considered significant. Both increases and decreases in accessibility were retained.

#### Direction of differential accessibility

For a comparison written as **Group 1 vs Group 2**, a **positive log2 fold-change** indicates higher chromatin accessibility in **Group 1**, while a **negative log2 fold-change** indicates higher chromatin accessibility in **Group 2**. For example, in a **3xAmp vs 1xAmp** comparison, a positive log2 fold-change indicates higher accessibility in **3xAmp**, whereas a negative log2 fold-change indicates higher accessibility in **1xAmp**. Thus, the order of the groups in the comparison determines the interpretation of the direction of differential accessibility.


#### How we calculated Peak-to-TSS annotation and picked the nearest nearby genes

For each significant ATAC peak, the distance to the nearest gene is calculated from the **centre of the peak interval to the gene's transcription start site (TSS)**. We used a **50 kb distance cutoff** for gene assignment, meaning that only genes with a TSS within 50 kb of the peak centre were considered. If multiple genes have exactly the same distance to a peak, the first gene encountered is assigned. Each peak is therefore assigned to **one nearest gene within the 50 kb cutoff**.

This represents our current annotation approach; **alternative peak-to-gene assignment strategies can also be considered if needed**.


##### Top nearby genes: Arid3a, Ltbp4, Gm16031, ENSMUSG00002075506, Ahnak2, ENSMUSG00000121159

![](figures/3xAmp_vs_1xAmp_volcano_annotated.png?v=1) 
[click here for full list of annotated significant peaks 3xAmp vs 1xAmp](https://docs.google.com/spreadsheets/d/1SegpfMxgroVsmVbqVqLy7K9hmB0wNabDq8jmlaRP48Y/edit?usp=sharing) 

##### Top nearby genes: Gtpbp4, Rn18s-rs5, ENSMUSG00002076241, Olfr431-ps1, Limd1, Ndrg1, Lgals3, Mylk2, Gm48776, Khsrp
![](figures/5xAmp_vs_1xAmp_volcano_annotated.png?v=1)

[click here for full list of annotated significant peaks 5xAmp vs 1xAmp](https://docs.google.com/spreadsheets/d/1A455BDEabNjRHyI3AWueVykz6jaaHqeDcAxdOwUgCRU/edit?usp=sharing) 
#### Top nearby genes: Fmnl2, Gm37982, Med13, Dnaja3, Gm45145, Mylk2, Aig1, Gm50077, Gm44210, D630003M21Rik
![](figures/5xAmp_vs_Control_volcano_annotated.png?v=1)

[click here for full list of annotated significant peaks 5xAmp vs Control](https://docs.google.com/spreadsheets/d/1Q0sx-n__Qs15Vokoe4TY2sds9Pd4-b1JngA-nthDp1c/edit?usp=sharing) 

###### More comparisions and analysis are ongoing 

## References 

SnapATAC2: Zhang, K., Zemke, N. R., Armand, E. J., & Ren, B. (2024). A fast, scalable and versatile tool for analysis of single-cell omics data. Nature methods, 21(2), 217-227.

