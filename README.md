# MGE Scanner

`mge-scanner` is a command-line tool designed to screen genomic neighborhoods around target genes or operons for Horizontal Gene Transfer (HGT) indicators, Mobile Genetic Element (MGE) markers, and local GC content anomalies. It automatically generates detailed text reports, summary tables, and publication-ready GC profile plots.

---

## Features

* **Flexible Operation:** 
  * **Gene Mode:** Directly scans specific target locus tags.
  * **Operon Mode:** Maps multi-gene operon IDs indexed by CSV file to scan neighborhood regions.
* **Flexible Scanning Modes:** Multiple strategy configurations for defining neighborhood boundaries and windows.
* **GC Content Anomaly Detection:** Calculates baseline genomic GC statistics and flags local Z-score drops.
* **Automated Visualization:** Generates detailed neighborhood GC plots highlighting target features and HGT indicators.

---

## Installation

Clone the repository and install the package locally in editable mode using `pip`:

```bash
git clone [https://github.com/your-username/mge_scanner.git
cd mge_scanner
pip install -e .
```


## Usage

Once installed, you can use the mge-scan command directly from your terminal.

```bash
Example Command (Gene Mode)
mge-scan --gbff genome.gbff \
         --mode gene \
         --source-type file \
         --target-list targets.txt \
         --out-dir ./output_results
```

| Argument | Description | Default |
| :--- | :--- | :--- |
| `--gbff` | Path to the GenBank (.gbff/.gbk) annotation file. | *Required* |
| `--mode` | Screening mode: either `gene` or `operon`. | `gene` |
| `--csv` | Path to CSV mapping file (required if using `operon` mode). | None |
| `--source-type` | Target input source: `file` or `folder`. | `file` |
| `--source-list` | Path to a text file containing target locus tags or operons. Required for 'file' mode. | None |
| `--source-folder` | Path to a folder with files containing target names. Requireed for 'folder' mode.  | None |
| `--file-pattern` | Patterns for file names around targets. | `*.0.faa` |
| `--out-dir` | Directory where all reports and plots will be saved. | `.mge_results/` |
| `--plot-dir` | Directory in out-dir where all plots will be saved. | `.meg_results/neighbourhood_plots` |
| `--out-zscore` | Save path for global GC abnormalies file. | `.mge_results/gc_anomolies.txt` |
| `--z-threshold` | Z-score threshold for local GC anomaly detection. | `-2.0` |
| `--z-window-size` | Window size in bp for chromosomal Z-score scan. | `250` |
| `--z-step-size` | Step size in bp for chromosomal Z-score scan. | `50` |
| `--max-genes` | Maximum number of genes up and downstream to scan for annotations. | `15` |
| `--max-bp` | Maximum number of basepairs up and downstream to scan for annotations. | `1500` |
| `--scan-mode` | Neighborhood search strategy (`a`, `b`, or `c`). | `a` |
| `--flank-bp` | Flanking region size in base pairs to scan for GC anomalies. | `1000` |




## Output Files

* **formatted_hits.txt:** A structured, easy-to-read text report organized by target gene/operon. (Often the most useful file for interpretation.)
* **detailed_report.txt:** Line-by-line detailed log of all detected HGT features and GC anomalies.
* **summary_report.tsv:** A machine-readable summary table of all targets and their flagged indicators.
* **gc_anomalies.txt:** All GC abnormalies across all amplicons from your gbff/gbk
* **neigourhood_plots:** A folder containing GC plots for your target and flanking region

License

This project is licensed under the MIT License - see the LICENSE file for details.


source-type 
target-list -> source-list
operons-folder -> source-folder
plot-dir 