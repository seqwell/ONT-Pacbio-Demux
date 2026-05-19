#!/usr/bin/env python3
import sys
import matplotlib
matplotlib.use("Agg")  # headless mode for Nextflow
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) != 3:
    print("Usage: leafplot.py <lengths.txt> <output.png>")
    sys.exit(1)

length_file = sys.argv[1]
out_png = sys.argv[2]
out_table = out_png.replace(".png", ".bin_table.txt")
out_weighted_png = out_png.replace(".png", "_weighted.png")

# Load read lengths
with open(length_file) as f:
    lengths = [int(line.strip()) for line in f if line.strip()]

# Define bins
bins = np.linspace(min(lengths), max(lengths), 101)

# Original histogram (for counts)
plt.figure(figsize=(10,5))
n, bin_edges, _ = plt.hist(lengths, bins=bins, color='forestgreen', alpha=0.7)
plt.title("Read Length Distribution")
plt.xlabel("Read Length (bp)")
plt.ylabel("Frequency")
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig(out_png, dpi=300)

# Prepare data for weighted bar plot
bin_centers = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in range(len(n))]
weighted_sums = [((bin_edges[i] + bin_edges[i+1]) / 2) * count for i, count in enumerate(n)]

# Save table with weighted_sum column
with open(out_table, "w") as out:
    out.write("bin_start\tbin_end\tcount\tweighted_sum\n")
    for i in range(len(n)):
        out.write(f"{int(bin_edges[i])}\t{int(bin_edges[i+1])}\t{int(n[i])}\t{weighted_sums[i]:.2f}\n")

# Create bar plot using weighted_sum
plt.figure(figsize=(12,5))
plt.bar(bin_centers, weighted_sums, width=(bin_edges[1]-bin_edges[0])*0.9, color='steelblue', alpha=0.8)
plt.title("Weighted Read Length Contribution")
plt.xlabel("Average Read Length (bp)")
plt.ylabel("Weighted Sum (bp * count)")
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig(out_weighted_png, dpi=300)

