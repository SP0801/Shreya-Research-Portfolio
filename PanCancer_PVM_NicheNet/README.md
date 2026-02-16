# Perivascular-like macrophages and PDGFRA/PDGFRB fibroblast ligands across cancers (GSE210347)

## Overview
This repository analyzes the pan-cancer scRNA-seq dataset **GSE210347** to identify a perivascular-like macrophage transcriptional program and prioritize stromal ligands that predict this program across cancers. The analysis uses **NicheNet** with **PDGFRA/PDGFRB-expressing fibroblasts** as the sender population and **perivascular-high macrophages** as the receiver population.

The analysis asks:

**Which fibroblast ligands are consistently ranked as candidate upstream regulators of a perivascular-like macrophage program across cancers?**

---

## Dataset
- **Accession:** GSE210347  
- **Modality:** single-cell RNA-seq  
- **Scope:** pan-cancer solid tumors  

---

## Definitions used in this analysis

### Macrophages
Macrophages were enriched from the myeloid compartment using **CD68** as an inclusion marker. CD68 is not macrophage-exclusive, so downstream marker expression is used to interpret macrophage states.

### Perivascular-like macrophages
A perivascular-like score was computed using a composite marker panel. Canonical markers include **LYVE1, FOLR2, MRC1, and CD163**, which have been used to describe vascular-associated or resident-like macrophage populations in tumors and tissue atlases (Ramos et al., 2022; Anstee et al., 2023; Bahri et al., 2023). Additional scavenger and regulatory markers were included to capture broader macrophage programs that co-occur with resident-like states in single-cell profiles.

Cells in the **top quartile** of the score were defined as **perivascular-high**. This cutoff is a pragmatic threshold and does not imply a discrete biological boundary.

### Fibroblast sender population
Fibroblasts were subset based on expression of **PDGFRA** and **PDGFRB**. This choice is intended to capture major CAF compartments described in tumor single-cell studies, including PDGFRA-associated fibroblast programs and PDGFRB-associated myofibroblast/perivascular stromal programs (Öhlund et al., 2017; Elyada et al., 2019). This does not imply that all PDGFRA/PDGFRB fibroblasts are physically perivascular.

---

## Workflow summary
1. Global embedding for dataset overview  
2. Subset myeloid cells  
3. Extract CD68+ macrophage-enriched compartment  
4. Score perivascular-like program and define PV-high cells  
5. Differential expression: PV-high vs PV-low  
6. Functional enrichment (KEGG)  
7. NicheNet per cancer: PDGFRA/PDGFRB fibroblasts (sender) → PV-high macrophages (receiver)  
8. Cross-cancer comparison of ligand rankings  

---

## Results summary 
- Perivascular-high macrophages are detectable across cancer types but vary in abundance by tissue.  
- PV-high macrophages show enrichment for lysosomal and phagocytic programs, consistent with scavenger-like macrophage states.  
- NicheNet prioritizes a small set of fibroblast ligands across cancers, with **CSF1** and **TGFB1** repeatedly ranked among the highest.  
- These results represent computational prioritization and should be interpreted as candidate regulatory signals rather than evidence of causality.  

---

## Figures (key outputs)

### 1) Global UMAP by cell type
![](figures/01_global_umap_tissue.png)

### 2) Myeloid embedding with CD68+ cells highlighted
![](figures/02_Myeloid_CD68.png)

### 3) Marker feature plots on macrophage embedding
![](figures/03_perivascular_markers.png)

### 4) Perivascular module score and PV-high gate
![](figures/04_perivascular_score.png)

### 5) PV score distribution across cancers
![](figures/05_PV%20macrophages%20per%20cancer%20type.png)

### 6) Per-cancer NicheNet ligand activity heatmap
![](figures/06_Per-cancer%20NicheNet%20ligand%20activity%20heatmap.png)

### 7) Conserved ligands summary across cancers
![](figures/07_Conserved%20ligands%20summary.png)

### 8) KEGG enrichment of PV-high macrophages
![](figures/08_pathway%20analysis.png)

---

## References
- Ramos RN et al. *Cell* 2022. Tissue-resident FOLR2+ macrophages associate with perivascular niches in tumor stroma.  
- Anstee JW et al. *Dev Cell* 2023. LYVE-1+ macrophages form CCR5-dependent perivascular niches in tumors.  
- Bahri R et al. *Front Immunol* 2023. Perivascular tumor-associated macrophages and their role in cancer.  
- Öhlund D et al. *J Exp Med* 2017. Distinct fibroblast populations in pancreatic cancer.  
- Elyada E et al. *Cancer Discov* 2019. CAF heterogeneity across species in pancreatic cancer.  
