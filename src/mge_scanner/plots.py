import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

def plot_neighbourhood_gc(record, hgt_start_bp, hgt_end_bp, target_start_bp, target_end_bp, target_id, output_dir, mode="gene", window_size=500, step_size=50):
    os.makedirs(output_dir, exist_ok=True)
    
    # Determine replicon type and description from record annotations
    replicon_type = "Replicon"
    if record.annotations.get("topology") == "circular" and len(record.seq) > 1000000:
        replicon_type = "Chromosome"
    elif "plasmid" in record.description.lower() or "plasmid" in record.id.lower():
        replicon_type = "Plasmid"
    else:
        replicon_type = "Replicon"

    replicon_name = record.id
    if record.description and record.description != "<unknown description>":
        replicon_name = f"{record.id} ({record.description})"

    padding = 2000
    plot_start = max(0, hgt_start_bp - padding)
    plot_end = min(len(record.seq), hgt_end_bp + padding)
    
    sub_seq = str(record.seq[plot_start:plot_end])
    gc_vals, positions = [], []
    
    for i in range(0, len(sub_seq) - window_size, step_size):
        window = sub_seq[i:i + window_size]
        g_count = window.upper().count('G')
        c_count = window.upper().count('C')
        gc_vals.append(((g_count + c_count) / window_size) * 100)
        positions.append(plot_start + i + (window_size // 2))
        
    if not gc_vals:
        return

    full_seq = str(record.seq).upper()
    all_window_gc = [
        ((full_seq[i:i+window_size].count('G') + full_seq[i:i+window_size].count('C')) / window_size) * 100
        for i in range(0, len(full_seq) - window_size, step_size)
    ]
    mean_gc = np.mean(all_window_gc) if all_window_gc else 0

    mode_label = "Gene" if mode.lower() == "gene" else "Operon"

    plt.figure(figsize=(10, 3.5))
    plt.axvspan(hgt_start_bp, hgt_end_bp, color='orange', alpha=0.2, label='HGT Neighborhood Span')
    plt.axvspan(target_start_bp, target_end_bp, color='forestgreen', alpha=0.35, label=f'Target {mode_label} {target_id}')
    plt.plot(positions, gc_vals, color='navy', linewidth=1.2, label='Local GC (%)', zorder=3)
    plt.axhline(mean_gc, color='red', linestyle='--', alpha=0.7, label=f'{replicon_type} Mean GC ({mean_gc:.1f}%)')
    
    ax = plt.gca()
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
    
    plt.xlabel(f"{replicon_type} Position (bp)", fontsize=10)
    plt.ylabel("GC Content (%)", fontsize=10)
    plt.title(f"Local GC Neighborhood Profile: Target {mode_label} {target_id}\n[{replicon_name}]", fontsize=10, fontweight='bold')
    plt.legend(loc='upper right', fontsize=8)
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, f"{mode_label.lower()}_{target_id}_gc_neighborhood.png")
    plt.savefig(output_path, dpi=300)
    plt.close()