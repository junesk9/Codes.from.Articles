#!/usr/bin/env python3

"""
Optimized motif analysis script with improved performance:
- Support for flexible k-mer length
- Efficient reverse complement handling
- Optimized mismatches calculation
- Multi-threading for both frequency generation and target analysis
- Progress tracking and caching

Usage:
$python3 optimized_motif_analysis.py 1_target.minusZn.up39.list [--motif_len 8] [--mismatch 1]

April 2025
"""

import sys
import random
import time
import numpy as np
import scipy.stats as stats
import statsmodels.sandbox.stats.multicomp as mc
import pandas as pd
import Levenshtein as lev
import multiprocessing as mp
from functools import partial
import os
import argparse
from collections import defaultdict
import itertools

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Motif analysis with flexible parameters")
    parser.add_argument("target_file", help="File containing target gene IDs")
    parser.add_argument("--motif_len", type=int, default=6, help="Length of motifs to analyze (default: 6)")
    parser.add_argument("--mismatch", type=int, default=0, help="Number of mismatches allowed (default: 0)")
    parser.add_argument("--prom_file", default="IRGSP-1.0_3kb-upstream_2024-07-12.fasta", 
                        help="Promoter sequence file (default: IRGSP-1.0_3kb-upstream_2024-07-12.fasta)")
    parser.add_argument("--prom_range", default="1000-0", help="Promoter range to analyze (default: 1300-0)")
    parser.add_argument("--repeats", type=int, default=1000, help="Number of random samplings (default: 10000)")
    parser.add_argument("--threads", type=int, default=24, 
                        help="Number of threads to use (default: automatic)")
    parser.add_argument("--lsi", action="store_true", help="Enable LSI annotation mode")
    parser.add_argument("--include_gene_ids", action="store_true", 
                        help="Include GeneIDs column in output (can make output file very large)")
    parser.add_argument("--rc", action="store_true", default=True,
                        help="Count reverse complements (default: enabled, use --no-rc to disable)")
    parser.add_argument("--no-rc", action="store_false", dest="rc",
                        help="Do not count reverse complements")
    return parser.parse_args()

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, 'N') for base in reversed(seq.upper()))

def process_kmer_batch(args):
    """Process a batch of kmers for multithreaded precomputation."""
    kmer_batch, all_kmers, mismatch_allowed = args
    
    forward_map_batch = {}
    rc_map_batch = {}
    
    for kmer in kmer_batch:
        rc_kmer = reverse_complement(kmer)
        
        # Maps for this kmer
        forward_matches = set()
        rc_matches = set()
        
        # Find all kmers within allowed mismatch distance
        for other_kmer in all_kmers:
            # Check forward matches
            if lev.distance(kmer, other_kmer) <= mismatch_allowed:
                forward_matches.add(other_kmer)
            
            # Check reverse complement matches
            if lev.distance(rc_kmer, other_kmer) <= mismatch_allowed:
                rc_matches.add(other_kmer)
        
        forward_map_batch[kmer] = forward_matches
        rc_map_batch[kmer] = rc_matches
    
    return forward_map_batch, rc_map_batch

def precompute_kmers_with_mismatches(motif_len, mismatch_allowed, num_threads=None):
    """
    Precompute dictionaries mapping each possible k-mer to:
    1. Its set of k-mers within mismatch_allowed distance
    2. Its set of reverse complements within mismatch_allowed distance
    
    Uses multi-threading to speed up computation for large k-mer lengths.
    """
    if num_threads is None:
        num_threads = max(1, mp.cpu_count() - 1)
    
    bases = ["A", "C", "T", "G"]
    all_kmers = ["".join(combo) for combo in itertools.product(bases, repeat=motif_len)]
    
    # Print information about scale of computation
    total_kmers = len(all_kmers)
    total_comparisons = total_kmers * total_kmers
    print(f"[PROGRESS] Precomputing mismatch maps for {total_kmers} {motif_len}-mers " +
          f"({total_comparisons:,} Levenshtein comparisons)")
    print(f"[PROGRESS] Using {num_threads} threads for precomputation")
    
    # Create maps for fast lookups
    forward_map = {}
    rc_map = {}
    
    # Split kmers into batches for multi-threading
    batch_size = max(1, total_kmers // (num_threads * 10))  # Ensure enough batches
    kmer_batches = [all_kmers[i:i+batch_size] for i in range(0, total_kmers, batch_size)]
    print(f"[PROGRESS] Processing {len(kmer_batches)} batches of ~{batch_size} k-mers each")
    
    # Set up process arguments
    process_args = [(batch, all_kmers, mismatch_allowed) for batch in kmer_batches]
    
    # Create a multiprocessing pool and process batches
    start_time = time.time()
    with mp.Pool(processes=num_threads) as pool:
        results = []
        # Add a tqdm-like progress indicator
        for i, batch_result in enumerate(pool.imap(process_kmer_batch, process_args)):
            results.append(batch_result)
            
            # Report progress periodically
            if (i+1) % max(1, len(kmer_batches)//20) == 0 or i+1 == len(kmer_batches):
                progress = (i+1) / len(kmer_batches)
                elapsed = time.time() - start_time
                remaining = (elapsed / progress) - elapsed if progress > 0 else 0
                
                elapsed_str = time.strftime("%H:%M:%S", time.gmtime(elapsed))
                remaining_str = time.strftime("%H:%M:%S", time.gmtime(remaining))
                
                print(f"[PROGRESS] Precomputation {progress*100:.1f}% complete - " +
                      f"Elapsed: {elapsed_str}, Remaining: {remaining_str}")
    
    # Merge results from all batches
    for forward_batch, rc_batch in results:
        forward_map.update(forward_batch)
        rc_map.update(rc_batch)
    
    total_time = time.time() - start_time
    print(f"[PROGRESS] Precomputation completed in {time.strftime('%H:%M:%S', time.gmtime(total_time))}")
    
    return forward_map, rc_map

def generate_motifs(motif_len):
    """Generate all possible DNA k-mers of given length."""
    bases = ["A", "C", "T", "G"]
    # Use itertools.product for efficient generation
    motifs = ["".join(combo) for combo in itertools.product(bases, repeat=motif_len)]
    return {m: 0 for m in motifs}

def count_motifs_with_mismatch_maps(seq, motif_len, forward_map, rc_map, count_rc=True):
    """
    Count motifs in a sequence using precomputed mismatch maps.
    This is significantly faster than recalculating Levenshtein distances.
    
    Args:
        seq: The DNA sequence to analyze
        motif_len: Length of motifs to count
        forward_map: Precomputed map of forward matches
        rc_map: Precomputed map of reverse complement matches
        count_rc: Whether to count reverse complements (default: True)
    """
    # Initialize counts for both forward and RC
    motif_counts = defaultdict(int)
    gene_found = defaultdict(bool)
    
    # Process each window in the sequence
    for i in range(len(seq) - motif_len + 1):
        window = seq[i:i+motif_len].upper()
        
        # Skip if window contains invalid characters
        if any(c not in "ACGT" for c in window):
            continue
        
        # Check forward matches
        for matching_kmer in forward_map.get(window, set()):
            motif_counts[matching_kmer] += 1
            gene_found[matching_kmer] = True
        
        # Check reverse complement matches if enabled
        if count_rc:
            for matching_kmer in rc_map.get(window, set()):
                motif_counts[matching_kmer] += 1
                gene_found[matching_kmer] = True
    
    return motif_counts, gene_found

def process_target_sequence(args):
    """Process a single target sequence for multithreading."""
    seq, motif_len, forward_map, rc_map, count_rc = args
    return count_motifs_with_mismatch_maps(seq, motif_len, forward_map, rc_map, count_rc)

def process_random_sample(args):
    """Process a single random sample for multithreading."""
    rnum, prom_dict, motif_len, forward_map, rc_map, count_rc = args
    
    # Select random genes
    random_genes = random.sample(list(prom_dict.keys()), rnum)
    seqs = [prom_dict[k] for k in random_genes]
    
    # Initialize counters
    gene_counts = defaultdict(int)
    motif_counts = defaultdict(int)
    
    # Process each sequence
    for seq in seqs:
        seq_motifs, gene_found = count_motifs_with_mismatch_maps(seq, motif_len, forward_map, rc_map, count_rc)
        
        # Update gene counts (each gene counts once per motif)
        for motif in gene_found:
            if gene_found[motif]:
                gene_counts[motif] += 1
            motif_counts[motif] += seq_motifs[motif]
    
    # Normalize counts
    seq_length_sum = sum(len(seq) - motif_len + 1 for seq in seqs)
    motif_freqs = {m: count/seq_length_sum for m, count in motif_counts.items()}
    gene_freqs = {m: count/rnum for m, count in gene_counts.items()}
    
    return (gene_freqs, motif_freqs)

def main():
    # Start timing
    stime = time.time()
    
    # Parse command line arguments
    if len(sys.argv) > 1:
        args = parse_arguments()
    else:
        # For backward compatibility
        class Args:
            pass
        args = Args()
        args.target_file = sys.argv[-1] if len(sys.argv) > 1 else "target.list"
        args.motif_len = 6
        args.mismatch = 0
        args.prom_file = "IRGSP-1.0_3kb-upstream_2024-07-12.fasta"
        args.prom_range = "1300-0"
        args.repeats = 1000
        args.threads = 24
        args.lsi = True
        args.include_gene_ids = False  # Default to not including gene IDs for backward compatibility
        args.rc = True  # Default to counting reverse complements
 
    # Define parameters
    prom_file = args.prom_file
    tgt_file = args.target_file
    prefix = tgt_file.split(".list")[0]
    motif_len = args.motif_len
    mismatch_allowed = args.mismatch
    
    # Read target genes
    with open(tgt_file) as f:
        tgt = [i.strip() for i in f]
    rNum = len(tgt)
    print(f"[PROGRESS] {rNum} genes are input with prefix as {prefix}")
    
    # Parameters for analysis
    prom_rng = args.prom_range
    rpt = args.repeats
    num_threads = args.threads if args.threads > 0 else max(1, mp.cpu_count() - 1)
    print(f"[PROGRESS] Using {num_threads} CPU cores for processing")
    
    # LSI annotation
    goi = {"Os03g0107300":"LSI2", "Os02g0745100":"LSI1", "Os10g0547500":"LSI3", "Os07g0257200":"Nramp5"}
    lsi_mode = args.lsi
    
    # Show settings
    print(f"[SETTINGS] Motif length: {motif_len}")
    print(f"[SETTINGS] Mismatches allowed: {mismatch_allowed}")
    print(f"[SETTINGS] Count reverse complements: {'yes' if args.rc else 'no'}")
    print(f"[SETTINGS] LSI mode: {'enabled' if lsi_mode else 'disabled'}")
    print(f"[SETTINGS] Include gene IDs: {'yes' if args.include_gene_ids else 'no (use --include_gene_ids to show)'}")

    # Update output filenames to indicate when reverse complements are not counted
    # Only add tag when RC is disabled, for backward compatibility with previous versions
    if args.rc:
        # Original filename format for backward compatibility
        wfreq = f"{prefix}.{prom_rng}nt.promoter.{motif_len}mer-{mismatch_allowed}mis.{rNum}n.{rpt}t.freq.txt"
    else:
        # Add 'norc' tag only when reverse complements are not counted
        wfreq = f"{prefix}.{prom_rng}nt.promoter.{motif_len}mer-{mismatch_allowed}mis.norc.{rNum}n.{rpt}t.freq.txt"
    
    outf = ".".join(wfreq.split(".")[:-2]) + ".z-test.txt"
    
    # Read promoter sequences
    dic = {}
    tid = []
    print(f"[PROGRESS] Reading promoter database: {prom_file}")
    with open(prom_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line.split("_")[0][1:] #GID by spliting
                dic[name] = []
                tid.append(name)
            else:
                dic[name].append(line)
    
    # Process promoter range
    prom_rng = [int(i) for i in prom_rng.split("-")]
    prom = {}
    gid_ls = []
    
    for t in tid:
        gid = t.split("-")[0]
        #gid = gid[:4] + "g" + gid[5:]
        if gid not in gid_ls:
            gid_ls.append(gid)
            seq = "".join(dic[t])
            # Subset the promoter sequence by the given range
            seq = seq[::-1][prom_rng[1]:prom_rng[0]]
            seq = seq[::-1]
            prom[gid] = seq
    
    # Get target sequences
    tgt_seqs = [prom[t] for t in tgt if t in prom]
    if len(tgt_seqs) < len(tgt):
        print(f"[WARNING] {len(tgt) - len(tgt_seqs)} target genes not found in promoter database")
    
    # Precompute k-mer mismatch maps with multi-threading
    forward_map, rc_map = precompute_kmers_with_mismatches(motif_len, mismatch_allowed, num_threads)
    
    # Process target data or generate frequency data
    try:
        # Try to use existing frequency file
        print(f"[PROGRESS] Looking for frequency file: {wfreq}")
        with open(wfreq) as f:
            wfreq_data = [i.strip().split("\t") for i in f][1:]
        wfreq_data = {i[0]: i[1:] for i in wfreq_data}
        print(f"[PROGRESS] Frequency file found, starting Z-statistics")
        
        # Process target sequences in parallel
        print(f"[PROGRESS] Processing {len(tgt_seqs)} target sequences")
        
        # Initialize counters for target analysis
        motif_gene_counts = defaultdict(int)
        motif_occurrence_counts = defaultdict(int)
        motif_gene_ids = defaultdict(list)
        
        # Prepare arguments for parallel processing
        process_args = [(seq, motif_len, forward_map, rc_map, args.rc) for seq in tgt_seqs]
        
        # Use multiprocessing for target analysis
        with mp.Pool(processes=num_threads) as pool:
            results = pool.map(process_target_sequence, process_args)
        
        # Collect all results
        for idx, (seq_motifs, gene_found) in enumerate(results):
            gene_id = tgt[idx] if idx < len(tgt) and tgt[idx] in prom else f"Unknown-{idx}"
            
            # Update occurrence counts
            for motif, count in seq_motifs.items():
                motif_occurrence_counts[motif] += count
            
            # Update gene counts
            for motif, found in gene_found.items():
                if found:
                    motif_gene_counts[motif] += 1
                    motif_gene_ids[motif].append(gene_id)
            
            if (idx + 1) % 10 == 0 or idx + 1 == len(results):
                print(f"[PROGRESS] {idx + 1} of {len(results)} target sequences processed")
        
        # Calculate sequence total length for normalization
        total_seq_length = sum(len(seq) - motif_len + 1 for seq in tgt_seqs)
        
        # Prepare columns for results DataFrame
        if args.include_gene_ids:
            columns = ["motif", "GeneC", "GeneIDs", "GeneFC", "GeneZ", "GeneP", 
                      "MotifC", "MotifFC", "MotifZ", "MotifP"]
        else:
            # Skip GeneIDs column to reduce file size
            columns = ["motif", "GeneC", "GeneFC", "GeneZ", "GeneP", 
                      "MotifC", "MotifFC", "MotifZ", "MotifP"]
        
        if lsi_mode:
            columns = columns[:2] + ["Lsi4"] + columns[2:]
        
        print(f"[PROGRESS] Calculating statistics")
        
        # Create results list first (faster than building DataFrame row by row)
        results_data = []
        
        # Calculate statistics for each motif
        for motif in motif_gene_counts:
            # Skip if not in frequency file
            if motif not in wfreq_data:
                continue
                
            # Motif occurrence statistics
            motif_count = motif_occurrence_counts[motif]
            motif_freq = motif_count / total_seq_length if total_seq_length > 0 else 0
            try:
                motif_exp_freq = float(wfreq_data[motif][2])
                motif_std = float(wfreq_data[motif][3])
            except (IndexError, ValueError):
                # Skip motifs with invalid data
                continue
                
            motif_fc = round(motif_freq / motif_exp_freq, 2) if motif_exp_freq > 0 else float('inf')
            motif_z = (motif_freq - motif_exp_freq) / motif_std if motif_std > 0 else 0
            motif_p = 1 - stats.norm.cdf(motif_freq, loc=motif_exp_freq, scale=motif_std) if motif_std > 0 else 1
            
            # Gene count statistics
            gene_count = motif_gene_counts[motif]
            gene_freq = gene_count / len(tgt_seqs) if tgt_seqs else 0
            try:
                gene_exp_freq = float(wfreq_data[motif][0])
                gene_std = float(wfreq_data[motif][1])
            except (IndexError, ValueError):
                # Skip motifs with invalid data
                continue
                
            gene_fc = round(gene_freq / gene_exp_freq, 2) if gene_exp_freq > 0 else float('inf')
            gene_z = (gene_freq - gene_exp_freq) / gene_std if gene_std > 0 else 0
            gene_p = 1 - stats.norm.cdf(gene_freq, loc=gene_exp_freq, scale=gene_std) if gene_std > 0 else 1
            gene_ids = ";".join(sorted(list(set(motif_gene_ids[motif]))))
            
            # Create row data for this motif
            row_data = [motif, gene_count]
            
            # Add LSI annotation if needed
            if lsi_mode:
                goi_set = set(gene_ids.split(";")).intersection(goi.keys())
                goi_ls = ";".join(sorted([goi[k] for k in list(goi_set)])) if goi_set else "NA"
                row_data.append(goi_ls)
            
            # Add gene IDs if requested
            if args.include_gene_ids:
                row_data.append(gene_ids)
            
            # Add remaining statistics
            row_data.extend([gene_fc, gene_z, gene_p, motif_count, motif_fc, motif_z, motif_p])
            
            # Add to results data
            results_data.append(row_data)
        
        # Create DataFrame from results data
        df = pd.DataFrame(results_data, columns=columns)
        
        # Add FDR corrections
        print(f"[PROGRESS] Applying FDR corrections")
        if not df.empty:
            gp = np.array(df.GeneP.tolist())
            gfdr = mc.multipletests(gp, method="fdr_bh")
            df["GeneFDR"] = gfdr[1]
            
            fp = np.array(df.MotifP.tolist())
            ffdr = mc.multipletests(fp, method="fdr_bh")
            df["MotifFDR"] = ffdr[1]
            
            # Reorder columns with FDR
            if args.include_gene_ids:
                columns = ["motif", "GeneC", "GeneIDs", "GeneFC", "GeneZ", "GeneP", "GeneFDR",
                          "MotifC", "MotifFC", "MotifZ", "MotifP", "MotifFDR"]
            else:
                columns = ["motif", "GeneC", "GeneFC", "GeneZ", "GeneP", "GeneFDR",
                          "MotifC", "MotifFC", "MotifZ", "MotifP", "MotifFDR"]
                
            if lsi_mode:
                columns = columns[:2] + ["Lsi4"] + columns[2:]
            
            df = df[columns]
            
            # Sort by gene z-score (descending)
            df = df.sort_values(by="GeneZ", ascending=False)
            
            # Write results
            df.to_csv(outf, sep="\t", index=False)
            print(f"[DONE] Output written to {outf}")
        else:
            print("[WARNING] No motifs found in target sequences that match frequency data")
        
    except FileNotFoundError:
        print(f"[PROGRESS] Frequency file [{wfreq}] not found, calculating frequencies")
        print(f"[PROGRESS] {rNum} promoters will be randomly picked, with {rpt} repeats")
        
        # Initialize data structures for random sampling
        all_motifs = set()
        gene_freq_data = defaultdict(list)
        motif_freq_data = defaultdict(list)
        
        # Set up multiprocessing
        pool = mp.Pool(processes=num_threads)
        
        # Determine batch size
        total_samples = rpt
        batch_size = min(1000, total_samples // (num_threads * 2))  # Smaller batches for better progress tracking
        if batch_size < 10:
            batch_size = 10
            
        # Calculate number of batches
        num_batches = (total_samples + batch_size - 1) // batch_size
        
        print(f"[PROGRESS] Processing {total_samples} samples in {num_batches} batches")
        
        # Process in batches
        completed = 0
        for batch_idx in range(num_batches):
            # Calculate batch size (handle last batch)
            current_batch_size = min(batch_size, total_samples - completed)
            
            print(f"[PROGRESS] Starting batch {batch_idx+1}/{num_batches} ({current_batch_size} samples)")
            
            # Prepare process arguments
            process_args = [(rNum, prom, motif_len, forward_map, rc_map, args.rc) for _ in range(current_batch_size)]
            
            # Process batch and collect results
            batch_results = pool.map(process_random_sample, process_args)
            
            # Update data with batch results
            for gene_freqs, motif_freqs in batch_results:
                # Update all_motifs set
                all_motifs.update(gene_freqs.keys())
                
                # Update frequency data
                for motif, freq in gene_freqs.items():
                    gene_freq_data[motif].append(freq)
                
                for motif, freq in motif_freqs.items():
                    motif_freq_data[motif].append(freq)
            
            # Update counters
            completed += current_batch_size
            
            # Report progress
            ctime = time.time()
            elapsed = ctime - stime
            remaining = (elapsed / completed) * (total_samples - completed) if completed > 0 else 0
            
            elapsed_str = time.strftime("%H:%M:%S", time.gmtime(elapsed))
            remaining_str = time.strftime("%H:%M:%S", time.gmtime(remaining))
            
            print(f"[PROGRESS] Completed {completed}/{total_samples} samples " +
                  f"({completed/total_samples*100:.1f}%) - " +
                  f"Elapsed: {elapsed_str}, Estimated remaining: {remaining_str}")
        
        pool.close()
        pool.join()
        
        # Write frequency file
        print(f"[PROGRESS] Writing frequency file: {wfreq}")
        with open(wfreq, "w") as out:
            # Write header
            print("motif", "gene.freq.avg", "gene.stdev", "motif.freq.avg", "motif.stdev", sep="\t", file=out)
            
            # Calculate and write statistics for each motif
            for motif in sorted(all_motifs):
                gene_freqs = gene_freq_data.get(motif, [0])
                motif_freqs = motif_freq_data.get(motif, [0])
                
                gene_avg = np.mean(gene_freqs)
                gene_std = np.std(gene_freqs)
                motif_avg = np.mean(motif_freqs)
                motif_std = np.std(motif_freqs)
                
                print(motif, gene_avg, gene_std, motif_avg, motif_std, sep="\t", file=out)
        
        print(f"[DONE] Frequency file generated: {wfreq}")
        print(f"[PROGRESS] Re-run the script to perform Z-statistics on your gene set")
    
    # Final timing
    total_time = time.time() - stime
    total_time_str = time.strftime("%H:%M:%S", time.gmtime(total_time))
    print(f"[DONE] Total execution time: {total_time_str}")

if __name__ == "__main__":
    main()
