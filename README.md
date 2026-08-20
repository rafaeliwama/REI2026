# Parasitism, smaller viromes and their consequences on the evolution of the intracellular innate immune system in barnacles
## Supplementary Files and Scripts


This repository contains the supplementary files needed to reproduce the virome
reconstruction analysis and the annotation of innate immune system components in
barnacles. Files are organized in the following directories:

1. `Virome_dir/` — data files, scripts, and notebooks used to run the virome
   analyses and generate the figures.

3. `Immune_annotation/` — the reference sequence sets used as BLAST databases for
   annotating the intracellular innate immune components.

Input data for both analysis is described in the Supplementary Material SXX of the online version of this article and is available on NCBI.


## Detailed directory and file descriptions.
### 1. `Virome_dir/` — virome reconstruction and analysis

Pipeline order: `Virome_inference_REI.sh` → `GetVirome.sh` → `KaijuToDF.py` (+ `ViromeUtils.py`) → `ViromeAnalysisREI2025.ipynb` → figures.

#### Scripts

| File | Description | Requirements |
| :--- | :--- | :--- |
| `Virome_inference_REI.sh` | Master driver. Downloads the reference databases (UniVec, SortMeRNA rRNA databases v4.3.4, Kaiju `kaiju_db_viruses_2024-08-15`), runs `GetVirome.sh` over the accession list, adds Linnean names to the Kaiju output, calls `KaijuToDF.py`, and tallies classified reads per sample. | BBMap; `GetVirome.sh`; `Accession_list_paired.txt` |
| `GetVirome.sh` | Per-sample loop. For each accession: downloads reads, quality/adapter trims against UniVec, removes rRNA reads, re-pairs surviving reads, and classifies them against the Kaiju viral database (`-E 10-5`). Writes one `<accession>.kaiju` per run. | sra-tools, seqyclean, SortMeRNA, BBMap, Kaiju |
| `KaijuToDF.py` | Builds the master virome matrix. Resolves the full Linnean rank path of each host TaxID, counts reads assigned to each viral taxon at the requested rank, and merges both into `final_virome_df.csv`.<br>Usage: `python3 KaijuToDF.py nodes.dmp names.dmp REI_SraRunInfo.csv genus` | pandas; `ViromeUtils.py`; Kaiju `.names` files in the working directory |
| `ViromeUtils.py` | Helper library imported by `KaijuToDF.py`: `getNCBIdics()` loads the NCBI taxonomy dump; `get_NodepPath()` and `getLRanks()` walk a TaxID to the root and return species → domain; `get_readCounts()` converts a Kaiju table into taxon → read-count pairs. | pandas |
| `ViromeAnalysisREI2025.ipynb` | All downstream analyses and figures — see the breakdown below. | pandas, numpy, matplotlib, scikit-bio, scikit-learn, scipy, pingouin |

The notebook downloads the ICTV Virus Metadata Resource (VMR MSL40) to flag bacteriophage genera and genome types, filters the matrix to eukaryotic viruses, computes the percentage of non-rRNA reads assigned to viruses per sample and per barnacle family, calculates Margalef's richness and Shannon's diversity per run, tests Rhizocephala against Thoracica (Levene's test and Welch's ANOVA), compares viral genera shared or unique to each infraclass alongside their genome types, and ordinates samples by non-metric MDS on Bray–Curtis distances.

#### Input and intermediate data

| File | Format | Description |
| :--- | :--- | :--- |
| `Accession_list_paired.txt` | one accession per line | The 43 NCBI SRA run accessions (SRR/DRR/ERR) with paired-end transcriptomic libraries used in the study. Input to `GetVirome.sh`. |
| `REI_SraRunInfo.csv` | CSV, 47 columns | Standard NCBI SRA Run Info export for the 44 runs considered, covering 42 species of thoracican and rhizocephalan barnacles plus outgroups. Supplies host TaxID, library strategy/layout, read length and depth. |
| `Kaiju_total_reads.txt` | headerless CSV | `total_reads,run_accession` — quality-filtered, rRNA-depleted reads submitted to Kaiju per sample. Denominator for normalising viral read counts. |
| `viralCountsReads.csv` | headerless CSV | `viral_reads,run_accession` — reads classified as viral per sample. Numerator for the percentage of viral reads. |
| `final_virome_df.csv` | tab-separated | Main virome matrix, one row per SRA run. Columns 1–20 hold run metadata and host Linnean ranks (`Run`, `BioSample`, `TaxID`, `size_MB`, `avgLength`, `bases`, `LibraryStrategy`, `LibrarySource`, `LibraryLayout`, species → domain); the remaining ~2,000 columns give read counts per viral genus. Primary input to the notebook. |

#### Figures

| File | Produced by | Description |
| :--- | :--- | :--- |
| `richness_diversity.pdf` | notebook | Margalef's richness and Shannon's diversity of viral genera, Rhizocephala vs. Thoracica. |
| `virome_figure_genometype.pdf` | notebook | Viral genera unique to Thoracica, unique to Rhizocephala, or shared, with the genome-type composition of each infraclass-specific fraction. |
| `mds_infraclass.pdf` | notebook | Non-metric MDS of Bray–Curtis distances on viral genus presence/absence, coloured by host infraclass. |

---

### 2. `Immune_annotation/` — reference sets for immune gene annotation

Amino acid FASTA files used as BLAST databases to annotate intracellular innate immune components in the barnacle transcriptomes. Sequences come from UniProtKB (`sp|` / `tr|` headers) and NCBI RefSeq (`XP_` headers), so each header retains its original accession, protein name and source organism.

| File | Target | Sequences | Species | Description |
| :--- | :--- | ---: | ---: | :--- |
| `File_S1.fasta` | cGAS | 752 | ~360 | Cyclic GMP-AMP synthase and cGAS-like proteins, including fragments and isoforms, across vertebrates and invertebrates. |
| `File_S2.fasta` | IRF | 2,593 | ~410 | Interferon regulatory factor family (IRF1–IRF9 and unassigned IRFs), including invertebrate IRF-like proteins. |
| `File_S3.fasta` | NF-κB / Rel | 10 | 10 | Curated NF-κB proteins (p65/RelA, p105/NF-κB1, p100/NF-κB2) from vertebrate and invertebrate models. |
| `File_S4.fasta` | RIG-I / RLR | 30 | 14 | RIG-I (DDX58) and related DExD/H-box RNA helicases. |
| `File_S5.fasta` | STING | 995 | ~690 | STING1/TMEM173 proteins and STING homologues, including the *Drosophila melanogaster* homologue. |
| `File_S6.fasta` | TLR | 12,536 | ~1,090 | Toll-like receptors TLR1–TLR13; the largest reference set in the repository. |
| `File_S7.fasta` | Vago | 1 | 1 | *Drosophila melanogaster* Vago (Q9VZ35), query for Vago-like sequences. |



## Acknowledgements


This README file was formatted with the aid of claude-opus-5. All data, code and scientific content are the responsibility of the authors.
