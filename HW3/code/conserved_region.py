import numpy as np
import pandas as pd
import math
from collections import Counter
import matplotlib.pyplot as plt

# Collect summary metrics and block lengths across test cases
summary = []
block_lengths_by_case = {}

for num in range(1, 6):
    fname = f"fast_cpp_result{num}.txt"
    try:
        with open(fname) as f:
            for _ in range(6):
                next(f)
            x = f.readline().strip()
            y = f.readline().strip()
            z = f.readline().strip()
    except FileNotFoundError:
        continue

    L = len(x)
    scores = [1.0 if len({a, b, c}) == 1 else (2/3 if len({a, b, c}) == 2 else 0.0)
              for a, b, c in zip(x, y, z)]

    # identify fully conserved blocks (score==1.0, length>=5)
    blocks = []
    in_block = False
    for i, s in enumerate(scores):
        if s == 1.0 and not in_block:
            start = i; in_block = True
        elif s != 1.0 and in_block:
            end = i - 1
            if end - start + 1 >= 5:
                blocks.append((start, end))
            in_block = False
    if in_block:
        end = L - 1
        if end - start + 1 >= 5:
            blocks.append((start, end))

    lengths = [end - start + 1 for start, end in blocks]
    total_conserved = sum(lengths)
    percent_conserved = total_conserved / L * 100

    summary.append({
        "test_case": num,
        "num_blocks": len(blocks),
        "avg_block_length": np.mean(lengths) if lengths else 0,
        "max_block_length": max(lengths) if lengths else 0,
        "total_conserved": total_conserved,
        "percent_conserved": percent_conserved
    })
    block_lengths_by_case[num] = lengths

# Build summary DataFrame
summary_df = pd.DataFrame(summary)
print("\nOverall Conserved Regions Summary:")
if summary_df.empty:
    print("No data found for any test case.")
else:
    print(summary_df.to_string(index=False))
    summary_df.to_csv("overall_conserved_summary.csv", index=False)
    print("\nSaved summary to 'overall_conserved_summary.csv'.")

    # Prepare plotting data
    cases = summary_df['test_case'].tolist()

    # Plot 1: Number of conserved blocks per test case
    plt.figure()
    plt.bar(cases, summary_df['num_blocks'])
    plt.xlabel("Test Case")
    plt.ylabel("Number of Conserved Blocks")
    plt.title("Conserved Blocks Count by Test Case")
    plt.tight_layout()
    plt.show()

    # Plot 2: Percent conserved positions per test case
    plt.figure()
    plt.bar(cases, summary_df['percent_conserved'])
    plt.xlabel("Test Case")
    plt.ylabel("Percent of Positions Conserved (%)")
    plt.title("Overall Conservation Percentage by Test Case")
    plt.tight_layout()
    plt.show()

    # Plot 3: Distribution of block lengths across test cases
    plt.figure()
    data = [block_lengths_by_case.get(case, []) for case in cases]
    plt.boxplot(data, labels=cases)
    plt.xlabel("Test Case")
    plt.ylabel("Conserved Block Length")
    plt.title("Conserved Block Length Distribution by Test Case")
    plt.tight_layout()
    plt.show()




# 1) Read aligned sequences
with open("fast_cpp_result5.txt") as f:
    for _ in range(6):
        next(f)
    x = f.readline().strip()
    y = f.readline().strip()
    z = f.readline().strip()

# 2) Compute per-column conservation scores (3-way match)
scores = []
for a, b, c in zip(x, y, z):
    unique = len({a, b, c})
    if unique == 1:
        scores.append(1.0)    # all three identical
    elif unique == 2:
        scores.append(2/3)    # two identical, one different
    else:
        scores.append(0.0)    # all different

# 3) Identify fully conserved blocks (score == 1.0) and filter by length >= threshold
min_block_len = 5
blocks = []
in_block = False
for i, s in enumerate(scores):
    if s == 1.0 and not in_block:
        start = i
        in_block = True
    elif s != 1.0 and in_block:
        end = i - 1
        if end - start + 1 >= min_block_len:
            blocks.append((start, end))
        in_block = False
if in_block:
    end = len(scores) - 1
    if end - start + 1 >= min_block_len:
        blocks.append((start, end))

# 4) Build report DataFrame
rows = []
for start, end in blocks:
    seq_block = x[start:end+1]
    length = end - start + 1
    consensus = seq_block  # fully conserved
    avg_score = np.mean(scores[start:end+1])
    # compute average entropy across block
    ent = 0.0
    for pos in range(start, end+1):
        freqs = Counter([x[pos], y[pos], z[pos]])
        total = sum(freqs.values())
        col_ent = 0.0
        for cnt in freqs.values():
            p = cnt / total
            col_ent -= p * math.log2(p)
        ent += col_ent
    entropy = ent / length
    rows.append({
        "start": start,
        "end": end,
        "length": length,
        "consensus": consensus,
        "avg_score": avg_score,
        "avg_entropy": entropy
    })
df = pd.DataFrame(rows)

# 5) Print report to console
print("\nConserved Blocks Report:")
print(df.to_string(index=False))

# Save report to CSV
df.to_csv("conserved_blocks_report.csv", index=False)
print("\nSaved detailed report to 'conserved_blocks_report.csv'.")

# 6) Conservation score plot
plt.figure(figsize=(20, 6))
plt.plot(scores, linewidth=1)
plt.xlabel("Alignment Position", fontsize=16)
plt.xticks(fontsize=14)
plt.ylabel("Conservation Score", fontsize=16)
plt.yticks(fontsize=14)
plt.title("Per-Column Conservation Score (3-way)", fontsize=18)
plt.tight_layout()
plt.show()

# 7) Stacked nucleotide counts per position
letters = ['A', 'C', 'G', 'T']
counts = {l: [] for l in letters}
for pos in range(len(x)):
    col = Counter([x[pos], y[pos], z[pos]])
    for l in letters:
        counts[l].append(col.get(l, 0))

positions = np.arange(len(x))
bottom = np.zeros(len(x))
plt.figure(figsize=(20, 6))
for l in ['A','C','G','T']:
    plt.bar(positions, counts[l], bottom=bottom, label=l)
    bottom += np.array(counts[l])

# only show 20 tick marks along the axis
n_ticks = 20
tick_locs = np.linspace(0, len(x)-1, n_ticks, dtype=int)
plt.xticks(tick_locs, tick_locs, rotation=45, fontsize=12)

plt.xlabel("Position", fontsize=16)
plt.ylabel("Nucleotide Count", fontsize=16)
plt.title("Stacked Nucleotide Counts per Position", fontsize=18)
plt.legend(title="Base")
plt.tight_layout()
plt.show()