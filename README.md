# Evaluating a Retrieval-Augmented Method for Curating Plant Gene Sets

This repository contains the full implementation of a retrieval-augmented pipeline for reconstructing plant Gene Ontology (GO) gene sets using large language models (LLMs). Using common wheat (*Triticum aestivum*) as the model organism, the pipeline reconstructs GO gene sets by mining the published literature, extracting candidate genes, verifying each gene-set association against supporting text, and consolidating the results into GMT files for downstream analysis.

The pipeline is demonstrated with **DeepSeek R1 8B**, but the LLM is interchangeable with any model served locally through [Ollama](https://ollama.com).

It includes all code for literature retrieval, gene-set generation, verification, gene identifier conversion, GMT construction, and evaluation.

---

## Installation and Requirements

This project was tested on:

```
Python 3.12.5
Ollama 0.11.11
```

Install Python dependencies:

```
pip install -r requirements.txt
```

Make sure Ollama is installed and running, and that the model you intend to use (e.g. DeepSeek R1 8B) is available locally. We ran Ollama in the background using:

```
nohup ollama serve > ollama.log 2>&1 &
```

This ensures the local LLM server remains active throughout long-running jobs.

### Computing Environment (Cluster)

Our experiments were run on a university research computing cluster, with each run utilizing a single V100-class GPU and 32 GB RAM. Average runtime per GO gene set reconstruction (including retrieval and LLM calls) was approximately 9 minutes, with variation across LLM configurations.

---

## Data Curation

Gene set definitions and gene identifier mappings are derived from Ensembl Plants BioMart exports. These files are too large to include in the repository and must be downloaded manually before running the pipeline.

### Downloading BioMart Exports

**IWGSC RefSeq v2.1 (required)**

1. Open the following URL in your browser:

```
http://plants.ensembl.org/biomart/martview?VIRTUALSCHEMANAME=plants_mart&ATTRIBUTES=tarefseqv2_eg_gene.default.feature_page.ensembl_gene_id|tarefseqv2_eg_gene.default.feature_page.ensembl_transcript_id|tarefseqv2_eg_gene.default.feature_page.external_gene_name|tarefseqv2_eg_gene.default.feature_page.name_1006|tarefseqv2_eg_gene.default.feature_page.definition_1006|tarefseqv2_eg_gene.default.feature_page.namespace_1003|tarefseqv2_eg_gene.default.feature_page.go_linkage_type|tarefseqv2_eg_gene.default.feature_page.go_id&FILTERS=&VISIBLEPANEL=attributepanel
```

2. Click **Results** and export as TSV.
3. Rename the file to `mart_export_2.1.tsv` and place it at `data/mart_exports/mart_export_2.1.tsv`.

**IWGSC RefSeq v1.1 (legacy, required only for v1.1 → v2.1 ID mapping)**

1. Open the following URL in your browser:

```
http://plants.ensembl.org/biomart/martview?VIRTUALSCHEMANAME=plants_mart&ATTRIBUTES=taestivum_eg_gene.default.feature_page.ensembl_gene_id|taestivum_eg_gene.default.feature_page.ensembl_transcript_id|taestivum_eg_gene.default.feature_page.external_gene_name|taestivum_eg_gene.default.feature_page.go_id|taestivum_eg_gene.default.feature_page.name_1006|taestivum_eg_gene.default.feature_page.definition_1006|taestivum_eg_gene.default.feature_page.go_linkage_type|taestivum_eg_gene.default.feature_page.namespace_1003&FILTERS=&VISIBLEPANEL=resultspanel
```

2. Click **Results** and export as TSV.
3. Rename the file to `mart_export.tsv` and place it at `data/mart_exports/mart_export.tsv`.

---

Once both files are in place, run:

- **`biomart_to_gmt.py`** — builds the reference GMT files and `out/go_terms_2.1.json`, which contains metadata for each GO term (gene set name, GO ID, GO domain, definition, etc.).
- **`join_biomart.py`** — builds a combined lookup table from the BioMart exports, enabling translation between identifier schemes (NCBI/Entrez, UniProtKB, WikiGene, protein domain, etc.).

The `v11_to_v21_mapping.json` file maps legacy v1.1 identifiers to v2.1 and is included in the repository.

---

## Pipeline Overview

The pipeline reconstructs gene sets through the following stages:

1. **Literature Retrieval, Gene Maker, and Gene Checker**
2. **Gene Identifier Conversion**
3. **GMT Construction**
4. **Evaluation and Plot Generation**

A direct-prompting baseline is also provided for comparison.

---

## 1. Main Pipeline (Literature Retrieval → Gene Maker → Gene Checker)

Run the main pipeline using:

```
python3 main_llm.py --input_file out/go_terms_2.1.json --llm <model_name>
```

For example:

```
python3 main_llm.py --input_file out/go_terms_2.1.json --llm deepseek-r1:8b
```

This performs:

* PubTator 3.0 literature retrieval and relevance grading
* LLM-based gene extraction ("gene maker", `rag_pipeline_gene_set_maker_llm.py`)
* LLM-based verification ("gene checker", `rag_pipeline_gene_checker_llm.py`)

All prompt templates are defined in `instructs.py`, and shared helper functions are in `utils.py`.

The output consists of extracted and verified gene sets for each GO term, stored under:

```
out/geneset_generations/<model>/
out/geneset_checks/<model>/
```

> **Note:** The maker may produce a `<gene set name>_raw.txt` file when a model response cannot be parsed as JSON. The corresponding gene set still has a valid JSON file; the raw files are ignored during the curation phase.

---

## 2. Gene Identifier Conversion

Genes extracted from the literature are reported using a wide range of naming schemes and must be resolved to stable RefSeq v2.1 TraesCS identifiers before GMT construction.

```
python3 gene_convert.py
```

Conversion uses the BioMart lookup table, the v1.1 → v2.1 mapping, and a FunPlantGenes mapping obtained via `scrape_funplantgenes.py` to assist with identifiers that cannot be resolved through BioMart alone. An optional LLM-assisted normalisation step is available via `gene_convert_llm.py` for identifiers that cannot be resolved by rule-based lookup. Helper conversion functions are provided in `convert_gene_ids.py` and `convert_to_v21.py`.

---

## 3. GMT Construction

To build the reconstructed GMT from the converted and validated genes:

```
python3 build_gmt.py \
    --results  out/conversion/conversion_results.json \
    --checks   out/geneset_checks/<model>/ \
    --original out/wheat_gmt_gene_stable_id_2.1.gmt \
    --out-reconstructed out/reconstructed.gmt \
    --out-subset        out/original_subset.gmt
```

This produces two GMT files:

* `out/reconstructed.gmt` — gene sets rebuilt from converted and validated genes
* `out/original_subset.gmt` — the matching rows from the original GMT, for direct comparison

Use the `--require-full-coverage` flag to include only gene sets where every converted gene has a corresponding check file.

---

## 4. Direct Prompting Baseline

To generate gene sets using direct prompting (no retrieval), run:

```
python3 direct_prompt_baseline.py --input_file out/go_terms_2.1.json --llm <model_name>
```

For example:

```
python3 direct_prompt_baseline.py --input_file out/go_terms_2.1.json --llm deepseek-r1:8b
```

The generated files are written to `out/geneset_generations_direct/<model>/`. To convert these genes and build a curated GMT for the baseline:

```
python3 direct_prompt_convert_and_gmt.py --model <model_name>
```

To restrict the direct-prompting GMT to the same gene sets as the reconstructed set (for a matched comparison):

```
python3 subset_direct_gmt.py \
    --reconstructed out/reconstructed.gmt \
    --direct        out/conversion_direct/<model>/direct_prompt_curated.gmt \
    --out           out/conversion_direct/<model>/direct_prompt_subset.gmt
```

---

## 5. Evaluation

Run evaluation to compute per-gene-set precision, recall, and F1, along with gene addition and loss statistics:

```
python3 evaluation.py --original_gmt <path> --new_gmt <path>
```

This produces, under `out/evaluation/`:

* `gene_set_comparison.csv`
* `gene_set_similarity.csv`
* `per_phenotype_prf.csv` — per-gene-set precision, recall, and F1
* `gene_analysis.txt` — summary statistics (mean loss, mean new genes, overall similarity)

Aggregate reconstruction statistics (mean/max new and lost genes, totals, etc.) are computed with:

```
python3 eval_stats.py
```

Gene set name matching between configurations is handled by `match_gmt.py`.

---

## 6. Plot Generation

The evaluation CSVs can be visualized using:

```
python3 compare_prf_plots.py   # precision, recall, and F1 plots
python3 lost_genes_plot.py     # distribution of % gene loss
python3 new_genes_plot.py      # distribution of new genes added
```

---

## Repository Notes

* The pipeline is model-agnostic: any Ollama-served model can be substituted via the `--llm` flag.
* `out/go_terms_2.1.json` contains the GO term metadata used as pipeline input.
* `data/mart_exports/` contains the BioMart exports for gene identifier conversion (once downloaded per the instructions above).
* All scripts assume that Ollama is installed and running locally.
