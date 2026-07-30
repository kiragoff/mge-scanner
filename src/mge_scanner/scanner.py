import os
import glob
import re
import pandas as pd
from Bio import SeqIO
import numpy as np
from mge_scanner.plots import plot_neighbourhood_gc

HGT_KEYWORDS = [
    r"transposase", r"insertion sequence", r"\bis element\b", r"transposon", 
    r"resolvase", r"invertase", r"phage", r"tail fiber", r"capsid", 
    r"portal protein", r"integrase", r"terminase", r"prophage", r"holin", 
    r"endolysin", r"baseplate", r"plasmid", r"mobilization", r"relaxase", 
    r"conjugative", r"\btra protein\b", r"type iv secretion", r"\bt4ss\b", 
    r"partitioning protein", r"genomic island", r"pathogenicity island", 
    r"mobility protein", r"recombinase", r"site-specific recombinase", 
    r"xer[CD]", r"restriction-modification", r"crispr", r"toxin-antitoxin",
    r"integron", r"integron integrase", r"\bintI\b", r"conjugal transfer", 
    r"conjugative transposon", r"type iv coupling protein", r"\bis\d+\b", r"antirepressor"
]
HGT_PATTERN = re.compile("|".join(HGT_KEYWORDS), re.IGNORECASE)

def detect_gc_anomalies(record, window_size=250, step_size=50, z_threshold=-2):
    full_seq = str(record.seq).upper()
    all_gc, all_positions = [], []
    for i in range(0, len(full_seq) - window_size, step_size):
        window = full_seq[i:i + window_size]
        gc_pct = ((window.count('G') + window.count('C')) / window_size) * 100
        all_gc.append(gc_pct)
        all_positions.append(i + (window_size // 2))
        
    if not all_gc:
        return []
        
    mean_gc, std_gc = np.mean(all_gc), np.std(all_gc)
    print(f"[{record.id}] Baseline GC: Mean = {mean_gc:.2f}%, Std Dev = {std_gc:.2f}%")
    
    anomalies = []
    for pos, gc in zip(all_positions, all_gc):
        z_score = (gc - mean_gc) / std_gc if std_gc > 0 else 0
        if z_score <= z_threshold:
            anomalies.append({"position": pos, "gc_pct": gc, "z_score": z_score})
    return anomalies

def run_pipeline(
    gbff_file, csv_file, input_source_type, operons_folder, target_list_file,
    detail_output_file, summary_output_file, formatted_output_file,
    plot_output_dir, zscore_output_file, max_genes_away, max_bp_away, 
    z_threshold, file_pattern="*.0.faa",
    z_window_size=250, z_step_size=50,
    scan_mode="a",
    flank_bp=1000,
    mode="gene"
):
    print("Loading GenBank annotations...")
    seq_records = SeqIO.to_dict(SeqIO.parse(gbff_file, "genbank"))
    
    # --- scan ALL REPLICONS (Chromosomes & Plasmids) for internal GC abnormalities ---
    all_zscore_records = []
    if zscore_output_file and not os.path.exists(zscore_output_file):
        for rec_id, rec in seq_records.items():
            rec_anomalies = detect_gc_anomalies(rec, window_size=z_window_size, step_size=z_step_size, z_threshold=z_threshold)
            for anom in rec_anomalies:
                all_zscore_records.append({
                    "SeqID": rec_id, 
                    "Start_bp": anom["position"] - (z_window_size // 2),
                    "End_bp": anom["position"] + (z_window_size // 2),
                    "GC_Content": anom["gc_pct"], 
                    "Z_Score": anom["z_score"]
                })
        
        with open(zscore_output_file, "w") as z_out:
            z_out.write("SeqID\tStart_bp\tEnd_bp\tGC_Content\tZ_Score\n")
            for z_rec in all_zscore_records:
                z_out.write(f"{z_rec['SeqID']}\t{z_rec['Start_bp']}\t{z_rec['End_bp']}\t{z_rec['GC_Content']:.2f}\t{z_rec['Z_Score']:.2f}\n")
        print(f"Saved global GC anomalies file: {zscore_output_file}")
    elif zscore_output_file and os.path.exists(zscore_output_file):
        print(f"Skipping GC anomalies generation: {zscore_output_file} already exists.")

    locus_to_info = {}
    for contig_id, record in seq_records.items():
        raw_cds_features = [f for f in record.features if f.type == "CDS"]
        deduplicated_cds, seen_tags = [], {}
        
        for feat in raw_cds_features:
            lt_list = feat.qualifiers.get("locus_tag", [])
            lt = lt_list[0] if lt_list else None
            if not lt:
                deduplicated_cds.append(feat)
                continue
            if lt not in seen_tags:
                seen_tags[lt] = feat
                deduplicated_cds.append(feat)
            else:
                existing_feat = seen_tags[lt]
                for key, values in feat.qualifiers.items():
                    if key not in existing_feat.qualifiers:
                        existing_feat.qualifiers[key] = values
                    else:
                        for val in values:
                            if val not in existing_feat.qualifiers[key]:
                                existing_feat.qualifiers[key].append(val)
                                
        for idx, feat in enumerate(deduplicated_cds):
            lt = feat.qualifiers.get("locus_tag", [None])[0]
            if lt:
                locus_to_info[lt] = {
                    "contig": contig_id,
                    "cds_index": idx,
                    "feature": feat,
                    "all_cds": deduplicated_cds,
                    "record": record
                }

    # --- LOAD TARGETS (Operons vs Genes) ---
    target_ids = set()
    if input_source_type == "folder":
        if not operons_folder:
            raise ValueError("--operons-folder must be specified when input-source-type is 'folder'.")
        search_path = os.path.join(operons_folder, file_pattern)
        print(f"Scanning folder for target files: {search_path}")
        target_files = glob.glob(search_path)
        suffix_to_strip = file_pattern.replace("*", "")
        target_ids = {os.path.basename(f).replace(suffix_to_strip, "") for f in target_files}
    elif input_source_type == "file":
        if not os.path.exists(target_list_file):
            raise FileNotFoundError(f"Target list file not found: {target_list_file}")
        with open(target_list_file, "r") as f:
            target_ids = {line.strip() for line in f if line.strip() and not line.startswith("#")}

    # --- SETUP MAPPING & DATA STRUCTURES BASED ON MODE ---
    operon_data = {}
    locus_to_operon = {}

    if mode == "gene":
        print("Running in Gene Mode: Bypassing CSV mapping, searching locus tags directly in GenBank.")
        for gene_id in target_ids:
            if gene_id in locus_to_info:
                operon_data[gene_id] = [gene_id]
            else:
                print(f"Warning: Target locus tag '{gene_id}' not found in GenBank file.")
    else:
        if not csv_file or not os.path.exists(csv_file):
            raise FileNotFoundError(f"CSV mapping file is required for operon mode: {csv_file}")
            
        df = pd.read_csv(csv_file)
        df["operon_id"] = df["operon_id"].astype(str).str.replace(r'\.0$', '', regex=True)
        df["locus_tag"] = df["locus_tag"].astype(str).str.strip()
        locus_to_operon = dict(zip(df["locus_tag"], df["operon_id"]))
        
        filtered_df = df[df["operon_id"].isin(target_ids)]
        operon_data = {op_id: group["locus_tag"].tolist() for op_id, group in filtered_df.groupby("operon_id")}

    detail_report_lines, structured_hits = [], {}
    summary_data = {
        op: {"has_hit": False, "hit_count": 0, "indicators_found": set(), "flank_z_anomaly": False, "flank_anomaly_details": []} 
        for op in target_ids
    }
    
    for operon_id, operon_locus_tags in operon_data.items():
        valid_tags = [lt for lt in operon_locus_tags if lt in locus_to_info]
        if not valid_tags:
            print(f"DEBUG WARNING: Target '{operon_id}' has no valid tags in locus_to_info!")
            continue
            
        contig_id = locus_to_info[valid_tags[0]]["contig"]
        all_cds = locus_to_info[valid_tags[0]]["all_cds"]
        record = locus_to_info[valid_tags[0]]["record"]

        indices = [locus_to_info[lt]["cds_index"] for lt in valid_tags if locus_to_info[lt]["contig"] == contig_id]
        if not indices:
            continue
            
        min_idx, max_idx = min(indices), max(indices)
        operon_start_bp = int(all_cds[min_idx].location.start)
        operon_end_bp = int(all_cds[max_idx].location.end)
        
        # --- SCAN MODES IMPLEMENTATION ---
        neighborhood_start_bp = operon_start_bp
        neighborhood_end_bp = operon_end_bp
        surrounding_neighborhood = []

        if scan_mode == "a":
            temp_up_gene_idx = max(0, min_idx - max_genes_away)
            temp_down_gene_idx = min(len(all_cds), max_idx + 1 + max_genes_away)
            
            hgt_candidates = []
            for idx in range(temp_up_gene_idx, temp_down_gene_idx):
                if idx >= min_idx and idx <= max_idx:
                    continue
                feat = all_cds[idx]
                prod = " ".join(feat.qualifiers.get("product", [""]))
                if HGT_PATTERN.search(prod):
                    stream_type = "upstream" if idx < min_idx else "downstream"
                    hgt_candidates.append({"stream": stream_type, "feature": feat, "index": idx})
            
            if hgt_candidates:
                min_hgt_idx = min(c["index"] for c in hgt_candidates)
                max_hgt_idx = max(c["index"] for c in hgt_candidates)
                core_start_bp = int(all_cds[min(min_idx, min_hgt_idx)].location.start)
                core_end_bp = int(all_cds[max(max_idx, max_hgt_idx)].location.end)
                final_upstream_idx = min(min_idx, min_hgt_idx)
                final_downstream_idx = max(max_idx + 1, max_hgt_idx + 1)
            else:
                core_start_bp = operon_start_bp
                core_end_bp = operon_end_bp
                final_upstream_idx = min_idx
                final_downstream_idx = max_idx + 1

            neighborhood_start_bp = max(0, core_start_bp - flank_bp)
            neighborhood_end_bp = min(len(record.seq), core_end_bp + flank_bp)

            upstream_features = all_cds[final_upstream_idx:min_idx]
            downstream_features = all_cds[max_idx + 1:final_downstream_idx]
            surrounding_neighborhood = [("upstream", f, i) for i, f in enumerate(upstream_features, start=final_upstream_idx)] + \
                                       [("downstream", f, i) for i, f in enumerate(downstream_features, start=max_idx + 1)]

        elif scan_mode == "b":
            final_upstream_idx = max(0, min_idx - max_genes_away)
            final_downstream_idx = min(len(all_cds), max_idx + 1 + max_genes_away)
            neighborhood_start_bp = max(0, operon_start_bp - max_bp_away)
            neighborhood_end_bp = min(len(record.seq), operon_end_bp + max_bp_away)
            surrounding_neighborhood = []

        elif scan_mode == "c":
            final_upstream_idx = max(0, min_idx - max_genes_away)
            final_downstream_idx = min(len(all_cds), max_idx + 1 + max_genes_away)
            neighborhood_start_bp = max(0, operon_start_bp - max_bp_away)
            neighborhood_end_bp = min(len(record.seq), operon_end_bp + max_bp_away)

            upstream_features = all_cds[final_upstream_idx:min_idx]
            downstream_features = all_cds[max_idx + 1:final_downstream_idx]
            surrounding_neighborhood = [("upstream", f, i) for i, f in enumerate(upstream_features, start=final_upstream_idx)] + \
                                       [("downstream", f, i) for i, f in enumerate(downstream_features, start=max_idx + 1)]

        # Flank GC Anomaly Check
        full_seq = str(record.seq).upper()
        all_window_gc = [
            ((full_seq[i:i+z_window_size].count('G') + full_seq[i:i+z_window_size].count('C')) / z_window_size) * 100
            for i in range(0, len(full_seq) - z_window_size, z_step_size)
        ]
        rec_mean = np.mean(all_window_gc) if all_window_gc else 0
        rec_std = np.std(all_window_gc) if all_window_gc else 0
        
        has_neighborhood_anomaly = False
        for w_start in range(neighborhood_start_bp, neighborhood_end_bp - z_window_size, z_step_size):
            w_end = w_start + z_window_size
            window_seq = full_seq[w_start:w_end]
            if len(window_seq) < z_window_size:
                continue
            w_gc = ((window_seq.count('G') + window_seq.count('C')) / z_window_size) * 100
            w_z = (w_gc - rec_mean) / rec_std if rec_std > 0 else 0
            
            if w_z <= z_threshold:
                has_neighborhood_anomaly = True
                already_logged = any(abs(anom["start"] - w_start) < z_window_size for anom in summary_data[operon_id]["flank_anomaly_details"])
                if not already_logged:
                    detail_report_lines.append(f"{operon_id}\tNeighborhood GC Anomaly\tN/A\tGC Drop Window (Z={w_z:.2f})\t{w_start} - {w_end} bp\tGC: {w_gc:.1f}%")
                    summary_data[operon_id]["flank_anomaly_details"].append({"start": w_start, "end": w_end, "z_score": w_z, "gc": w_gc})
                    
        if has_neighborhood_anomaly:
            summary_data[operon_id]["flank_z_anomaly"] = True

        # HGT Keyword Evaluation
        has_operon_hit = False
        hit_min_bp, hit_max_bp = int(operon_start_bp), int(operon_end_bp)

        for direction, feat, feat_idx in surrounding_neighborhood:
            product = " ".join(feat.qualifiers.get("product", ["hypothetical protein"]))
            lt_tag = feat.qualifiers.get("locus_tag", ["unknown"])[0]
            
            if HGT_PATTERN.search(product):
                has_operon_hit = True
                summary_data[operon_id]["has_hit"] = True
                summary_data[operon_id]["hit_count"] += 1
                summary_data[operon_id]["indicators_found"].add(product)
                hit_min_bp = min(hit_min_bp, int(feat.location.start))
                hit_max_bp = max(hit_max_bp, int(feat.location.end))

                dist_bp = int(operon_start_bp - feat.location.end) if direction == "upstream" else int(feat.location.start - operon_end_bp)
                genes_away = min_idx - feat_idx if direction == "upstream" else feat_idx - max_idx

                detail_report_lines.append(f"Operon: {operon_id}\tDirection: {direction}\tHGT Indicator Locus: {lt_tag}\tProduct: {product}\tDistance: {dist_bp} bp ({genes_away} genes away)")

                if operon_id not in structured_hits:
                    structured_hits[operon_id] = []

                collected_items = []
                target_range = range(min_idx, max_idx + 1)
                span_range = range(feat_idx, min_idx) if direction == "upstream" else range(max_idx + 1, feat_idx + 1)
                
                for s_idx in list(span_range) + list(target_range):
                    s_feat = all_cds[s_idx]
                    s_lt = s_feat.qualifiers.get("locus_tag", ["unknown"])[0]
                    collected_items.append({
                        "stream": "upstream" if s_idx < min_idx else ("target" if s_idx <= max_idx else "downstream"),
                        "operon": locus_to_operon.get(s_lt, "NA"),
                        "feature": s_feat,
                        "product": " ".join(s_feat.qualifiers.get("product", ["hypothetical protein"]))
                    })
                
                collected_items.sort(key=lambda x: int(x["feature"].location.start))
                structured_hits[operon_id].extend(collected_items)

        if has_operon_hit or has_neighborhood_anomaly:
            plot_neighbourhood_gc(record, hit_min_bp, hit_max_bp, int(operon_start_bp), int(operon_end_bp), operon_id, output_dir=plot_output_dir, mode=mode)

    mode_label = "Gene" if mode.lower() == "gene" else "Operon"
    col_header = "Gene" if mode.lower() == "gene" else "Operon"

    # --- OUTPUTS ---
    with open(detail_output_file, "w") as out:
        out.write(f"# {mode_label} HGT Neighborhood Detailed Screen Report\n")
        out.write("\n".join(detail_report_lines))
        
    summary_rows = []
    for op_id, data in sorted(summary_data.items(), key=lambda x: str(x[0])):
        summary_rows.append({
            "operon_id": op_id,
            "has_hgt_hit": data["has_hit"],
            "flank_gc_anomaly_z_le_-2": data["flank_z_anomaly"],
            "total_hits": data["hit_count"],
            "indicator_products": "; ".join(sorted(data["indicators_found"])) if data["indicators_found"] else "None"
        })
    
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(summary_output_file, sep="\t", index=False)
    
    print(f"Generating formatted hit list report: {formatted_output_file}")
    with open(formatted_output_file, "w") as out:
        out.write("=" * 115 + "\n")
        out.write(f"DETAILED HGT HIT LIST - SEPARATED BY TARGET {mode_label.upper()}\n")
        out.write("=" * 115 + "\n\n")
        
        all_target_ids = sorted(list(summary_data.keys()))
        
        if not all_target_ids:
            out.write(f"No target {mode_label.lower()}s found to scan.\n")
        else:
            for target_op in all_target_ids:
                out.write(f"Target {mode_label}: {target_op}\n")
                out.write("-" * 115 + "\n")
                out.write(f"{'STREAM':<12} | {col_header:<8} | {'GENE':<15} | {'BP RANGE':<22} | {'PRODUCT':<48}\n")
                out.write("-" * 115 + "\n")
                
                seen_rows = set()
                rows_to_print = []
                
                # UNCONDITIONALLY ADD THE TARGET ITSELF
                if target_op in operon_data:
                    for lt in operon_data[target_op]:
                        if lt in locus_to_info:
                            feat = locus_to_info[lt]["feature"]
                            start_bp = int(feat.location.start)
                            end_bp = int(feat.location.end)
                            bp_str = f"{start_bp:,} - {end_bp:,}"
                            product = " ".join(feat.qualifiers.get("product", ["Unknown product"]))
                            
                            rows_to_print.append({
                                "sort_bp": start_bp,
                                "stream": "target",
                                "operon": "NA",  
                                "gene": lt,
                                "bp": bp_str,
                                "product": product
                            })

                # ADD ANY SURROUNDING HGT HITS 
                if target_op in structured_hits:
                    for hit in structured_hits[target_op]:
                        # Skip adding the target twice if it's already caught in structured_hits
                        feat = hit["feature"]
                        locus_tag = feat.qualifiers.get("locus_tag", ["N/A"])[0]
                        if target_op in operon_data and locus_tag in operon_data[target_op]:
                            continue
                            
                        start_bp = int(feat.location.start)
                        end_bp = int(feat.location.end)
                        bp_str = f"{start_bp:,} - {end_bp:,}"
                        product = " ".join(feat.qualifiers.get("product", ["Unknown product"]))
                        
                        rows_to_print.append({
                            "sort_bp": start_bp,
                            "stream": hit["stream"],
                            "operon": hit["operon"],
                            "gene": locus_tag,
                            "bp": bp_str,
                            "product": product
                        })
                
                # ADD FLANK ANOMALIES (if any)
                if target_op in summary_data and summary_data[target_op].get("flank_z_anomaly"):
                    for anom in summary_data[target_op].get("flank_anomaly_details", []):
                        start_bp = int(anom["start"])
                        end_bp = int(anom["end"])
                        bp_str = f"{start_bp:,} - {end_bp:,}"
                        prod_str = f"Local Flank Anomaly (Z={anom['z_score']:.2f}, GC={anom['gc']:.1f}%)"
                        
                        rows_to_print.append({
                            "sort_bp": start_bp,
                            "stream": "flank",
                            "operon": "NA",
                            "gene": "Flank Region",
                            "bp": bp_str,
                            "product": prod_str
                        })
                
                # Sort everything cleanly by genomic coordinate (bp)
                rows_to_print.sort(key=lambda x: x["sort_bp"])
                
                merged_rows = []
                current_flank_block = None
                
                for r in rows_to_print:
                    if r["stream"] == "flank":
                        parts = r["bp"].replace(",", "").split(" - ")
                        f_start, f_end = int(parts[0]), int(parts[1])
                        
                        if current_flank_block is None:
                            current_flank_block = {
                                "sort_bp": f_start,
                                "stream": "flank",
                                "operon": "NA",
                                "gene": "Flank Region",
                                "start_bp": f_start,
                                "end_bp": f_end,
                                "count": 1,
                                "z_scores": [float(r["product"].split("Z=")[1].split(",")[0])],
                                "gcs": [float(r["product"].split("GC=")[1].replace("%", "").replace(")", "").strip())]
                            }
                        else:
                            if f_start <= current_flank_block["end_bp"]:
                                current_flank_block["end_bp"] = max(current_flank_block["end_bp"], f_end)
                                current_flank_block["count"] += 1
                                z_val = float(r["product"].split("Z=")[1].split(",")[0])
                                gc_val = float(r["product"].split("GC=")[1].replace("%", "").replace(")", "").strip())
                                current_flank_block["z_scores"].append(z_val)
                                current_flank_block["gcs"].append(gc_val)
                            else:
                                avg_z = sum(current_flank_block["z_scores"]) / len(current_flank_block["z_scores"])
                                avg_gc = sum(current_flank_block["gcs"]) / len(current_flank_block["gcs"])
                                c = current_flank_block["count"]
                                merged_rows.append({
                                    "sort_bp": current_flank_block["sort_bp"],
                                    "stream": "flank",
                                    "operon": "NA",
                                    "gene": "Flank Region",
                                    "bp": f"{current_flank_block['start_bp']:,} - {current_flank_block['end_bp']:,}",
                                    "product": f"Local Flank Anomaly (Merged {c}, Avg Z={avg_z:.2f}, Avg GC={avg_gc:.1f}%)"
                                })
                                current_flank_block = {
                                    "sort_bp": f_start,
                                    "stream": "flank",
                                    "operon": "NA",
                                    "gene": "Flank Region",
                                    "start_bp": f_start,
                                    "end_bp": f_end,
                                    "count": 1,
                                    "z_scores": [float(r["product"].split("Z=")[1].split(",")[0])],
                                    "gcs": [float(r["product"].split("GC=")[1].replace("%", "").replace(")", "").strip())]
                                }
                    else:
                        if current_flank_block is not None:
                            avg_z = sum(current_flank_block["z_scores"]) / len(current_flank_block["z_scores"])
                            avg_gc = sum(current_flank_block["gcs"]) / len(current_flank_block["gcs"])
                            c = current_flank_block["count"]
                            merged_rows.append({
                                "sort_bp": current_flank_block["sort_bp"],
                                "stream": "flank",
                                "operon": "NA",
                                "gene": "Flank Region",
                                "bp": f"{current_flank_block['start_bp']:,} - {current_flank_block['end_bp']:,}",
                                "product": f"Local Flank Anomaly (Merged {c}, Avg Z={avg_z:.2f}, Avg GC={avg_gc:.1f}%)"
                            })
                            current_flank_block = None
                        merged_rows.append(r)
                
                if current_flank_block is not None:
                    avg_z = sum(current_flank_block["z_scores"]) / len(current_flank_block["z_scores"])
                    avg_gc = sum(current_flank_block["gcs"]) / len(current_flank_block["gcs"])
                    c = current_flank_block["count"]
                    merged_rows.append({
                        "sort_bp": current_flank_block["sort_bp"],
                        "stream": "flank",
                        "operon": "NA",
                        "gene": "Flank Region",
                        "bp": f"{current_flank_block['start_bp']:,} - {current_flank_block['end_bp']:,}",
                        "product": f"Local Flank Anomaly (Merged {c}, Avg Z={avg_z:.2f}, Avg GC={avg_gc:.1f}%)"
                    })
                
                if not merged_rows:
                    out.write(f"{'flank':<12} | {'NA':<8} | {'None':<15} | {'N/A':<22} | {'No HGT hits or flank anomalies detected':<48}\n")
                else:
                    for r in merged_rows:
                        unique_key = (r["stream"], r["gene"], r["bp"])
                        if unique_key not in seen_rows:
                            seen_rows.add(unique_key)
                            out.write(f"{r['stream']:<12} | {r['operon']:<8} | {r['gene']:<15} | {r['bp']:<22} | {r['product']:<48}\n")
                out.write("\n")

    print(f"Done! All reports saved successfully.")