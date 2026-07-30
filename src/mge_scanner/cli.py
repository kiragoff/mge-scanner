import os
import argparse
from mge_scanner.scanner import run_pipeline

def main():
    parser = argparse.ArgumentParser(
        description="Screen genomic neighborhoods around target genes or operons for HGT/MGE indicators and GC anomalies."
    )
    
    # Required parameters
    parser.add_argument("--gbff", required=True, help="Path to genome GenBank (.gbff/.gbk) file.")
    
    # Mode and CSV parameters
    parser.add_argument(
        "--mode", 
        choices=["operon", "gene"], 
        default="gene", 
        help="Operational mode: 'operon' (requires CSV mapping) or 'gene' (direct GenBank locus tag search)."
    )
    parser.add_argument("--csv", required=False, help="Path to CSV operon-to-locus mapping file (Required if --mode operon).")
    
    # Input source options
    parser.add_argument("--source-type", choices=["folder", "file"], default="folder", help="Target source type ('folder' or 'file').")
    parser.add_argument("--source-list", help="Path to text file containing target IDs (operons or genes depending on mode).")
    parser.add_argument("--source-folder", help="Path to folder containing target files (if source-type is 'folder').")
    parser.add_argument("--file-pattern", default="*.0.faa", help="Filename pattern/extension for target files in folder mode (default: '*.0.faa').")
    
    # Output file paths
    parser.add_argument("--out-dir", default="mge_results", help="Output directory for results.")
    parser.add_argument("--plot-dir", default="neighbourhood_plots", help="Directory to save local GC neighborhood plots.")
    parser.add_argument("--out-zscore", default="gc_anomalies.txt", help="Path for file with replicon GC anomalies.")
    
    # Threshold options
    parser.add_argument("--z-threshold", type=float, default=-2, help="Z-score threshold for GC anomalies (default: -2).")
    parser.add_argument("--max-genes", type=int, default=15, help="Maximum genes away to scan (default 15).")
    parser.add_argument("--max-bp", type=int, default=15000, help="Maximum base pairs away to scan (default 15000).")
    parser.add_argument("--z-window-size", type=int, default=250, help="Window size in bp for chromosomal Z-score scan (default: 250).")
    parser.add_argument("--z-step-size", type=int, default=50, help="Step size in bp for chromosomal Z-score scan (default: 50).")

    # Z-score scanning options
    parser.add_argument(
        "--scan-mode", 
        choices=["a", "b", "c"], 
        default="a", 
        help=(
            "Scanning behavior mode: "
            "A = Operon/Gene-to-HGTe (and their flanks, or HGTe-HGTe if upstream/downstream found); "
            "B = Strict Z-score anomalies limited to immediate flanks; "
            "C = All Z-score anomalies in the windows around targets."
        )
    )
    parser.add_argument("--flank-bp", type=int, default=1000, help="Maximum base pairs away to scan for local flank GC anomalies.")

    args = parser.parse_args()
    
    # Validate conditional requirements
    if args.mode == "operon" and not args.csv:
        parser.error("--csv mapping file is required when --mode is set to 'operon'.")

    # Ensure the base output directory exists
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Build full paths using out_dir
    detail_path = os.path.join(args.out_dir, "detailed_report.txt")
    summary_path = os.path.join(args.out_dir, "summary_report.txt")
    formatted_path = os.path.join(args.out_dir, "formatted_hits.txt")
    plot_dir = os.path.join(args.out_dir, "neighbourhood_plots")
    zscore_path = os.path.join(args.out_dir, "gc_anomalies.txt")

    run_pipeline(
        gbff_file=args.gbff,
        mode=args.mode,
        csv_file=args.csv,
        input_source_type=args.source_type,
        operons_folder=args.source_folder,
        target_list_file=args.source_list,
        file_pattern=args.file_pattern,
        detail_output_file=detail_path,
        summary_output_file=summary_path,
        formatted_output_file=formatted_path,
        plot_output_dir=plot_dir,
        zscore_output_file=zscore_path,
        max_genes_away=args.max_genes,
        max_bp_away=args.max_bp,
        z_threshold=args.z_threshold,
        scan_mode=args.scan_mode,
        flank_bp=args.flank_bp
    )

if __name__ == "__main__":
    main()