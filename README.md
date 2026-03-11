# Shreya Prakash | Research Portfolio

I am an immunologist with a growing foundation in computational biology, currently working as a Junior Research Fellow at the Laboratory Oncology Unit, AIIMS New Delhi. My research sits at the intersection of immune cell biology, transcriptomics, and clinical biomarker discovery.

My wet lab background spans primary human immune cell culture, multiparametric flow cytometry, ELISA, and RNA-seq library preparation from low-input clinical material including cell-free RNA. Alongside ongoing research at AIIMS, I have been building computational skills independently, learning by applying tools to biological questions I find genuinely interesting rather than working through tutorials in isolation. This repository documents that process.

The projects here are at different stages and serve different purposes. One reflects active research (pipeline shared without patient data, pending publication). The other is a self-directed analysis on a public dataset, undertaken to develop specific skills while staying anchored to questions that connect to my broader research interests.

---

## Projects

### 1. Pan-cancer perivascular-like macrophages and CAF ligand prioritisation
`PanCancer_PVM_NicheNet/`

A self-directed single-cell RNA-seq analysis using the pan-cancer dataset GSE210347 (~855,000 cells across 11 solid tumour types). I wanted to learn NicheNet and Seurat at scale, and chose a question I found biologically meaningful: whether PDGFRA/PDGFRB-expressing cancer-associated fibroblasts share a conserved set of ligands that predict a perivascular-like transcriptional program in tumour macrophages across cancer types.

The analysis identifies a LYVE1/FOLR2/MRC1/CD163-anchored perivascular module, scores it across the myeloid compartment, and runs NicheNet per cancer type with fibroblasts as the sender population. CSF1 and TGFB1 emerge as the most consistently ranked ligands. KEGG enrichment of perivascular-high macrophages shows lysosomal and phagocytic pathway enrichment and, notably, rheumatoid arthritis and cholesterol metabolism pathways, which connect this tumour-focused analysis to my broader interest in macrophage states in inflammatory disease.

**Tools:** R, Seurat, NicheNet, clusterProfiler, ggplot2, pheatmap  
**Data:** GSE210347 (public)  
**Status:** Complete

---

### 2. AML liquid biopsy: cfRNA and PBMC transcriptomics pipeline
`AML_cfRNA_Pipeline/` *(coming soon)*

This documents the RNA-seq analysis pipeline from my ongoing JRF at AIIMS New Delhi, where I am working on cell-free RNA as a potential liquid biopsy analytic in acute myeloid leukaemia. The biological question is whether plasma-derived cfRNA can reflect disease-relevant immune and tumour-associated transcriptional signatures that mirror what is observed in matched PBMCs, and whether cfRNA-specific signals might provide additional information not accessible from blood cells alone.

The pipeline runs from raw FASTQ files through fastp trimming, STAR alignment to GRCh38/GENCODE v47, and featureCounts quantification (paired-end, reverse-stranded). After CPM filtering and TMM normalisation via edgeR, differential expression is run using limma-voom with batch included as a covariate in the design matrix, so that batch variance is accounted for within the linear model rather than by pre-correcting the data. Three contrasts are tested: AML cfRNA vs healthy controls, AML PBMCs vs healthy controls, and cfRNA vs PBMCs within AML patients. A three-way overlap strategy identifies 1,907 candidate biomarker genes, including 1,416 protein-coding genes and 366 lncRNAs, that are differentially expressed in cfRNA relative to both blood-based comparators. For visualisation only (PCA, heatmaps, z-score plots), `removeBatchEffect` is applied to the voom-transformed expression matrix; this corrected matrix is not used for differential expression.

Patient data are not publicly available pending publication. The pipeline is demonstrated on a simulated count matrix with matching structure and metadata.

**Tools:** fastp, STAR, featureCounts, edgeR, limma-voom, biomaRt, R  
**Data:** Simulated (real data: AIIMS, unpublished)  
**Status:** In preparation

---

## Skills

**Languages:** R (primary), bash  
**Single-cell:** Seurat, NicheNet, UMAP, cell type annotation, module scoring  
**Bulk RNA-seq:** STAR, featureCounts, limma-voom, edgeR, DESeq2, TMM normalisation  
**Functional analysis:** clusterProfiler (KEGG, GO), pathway enrichment, NicheNet ligand activity  
**Wet lab:** PBMC isolation, flow cytometry (multiparametric), ELISA, cfRNA isolation, low-input RNA-seq library preparation, confocal microscopy  

---

## Background

**MSc Immunology and Inflammatory Disease**, University of Glasgow (Merit, 2022-23)  
Dissertation: PIEZO1-mediated macrophage activation in rheumatoid arthritis (Prof. Mariola Kurowska-Stolarska)

**Junior Research Fellow**, Laboratory Oncology Unit, AIIMS New Delhi (2023-present)  
PI: Dr. Deepshi Thakral | Focus: cfRNA liquid biopsy in AML

**Selected presentations:** BSI Congress 2023, BSI SIGNET 2024, EMDCS 2024, Global Immunology Summit 2025  
**BSI Member**
