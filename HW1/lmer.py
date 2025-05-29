import random
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd

# Function to count l-mers in a sequence
def count_lmers(seq, l):
    return Counter(seq[i:i+l] for i in range(len(seq)-l+1))

# Function to load a genome from a FASTA file
def load_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq.append(line.strip().upper())
    return ''.join(seq)

# Parameters
l = 5
genome_path = 'CP011331.1.fasta'

# Load genome and count l-mers
genome_seq = load_fasta(genome_path)
genome_counts = count_lmers(genome_seq, l)
genome_dist = Counter(genome_counts.values())

# Generate a random sequence of the same length and count l-mers
random_seq = ''.join(random.choice('ACGT') for _ in range(len(genome_seq)))
random_counts = count_lmers(random_seq, l)
random_dist = Counter(random_counts.values())

# Prepare DataFrames for display
df_genome = pd.DataFrame(sorted(genome_dist.items()), columns=['Occurrences', 'Num_lmers'])
df_random = pd.DataFrame(sorted(random_dist.items()), columns=['Occurrences', 'Num_lmers'])

# Plot two distributions of l-mers using subplots
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.scatter(df_genome['Occurrences'], df_genome['Num_lmers'], color='blue', alpha=0.7, label='Genome')
plt.title(f'L-mers Distribution in Genome (l={l})')
plt.xlabel(f'Frequency of {l}-mers')
plt.ylabel(f'Number of distinct {l}-mers')
plt.legend()
plt.grid(True)


plt.subplot(1, 2, 2)
plt.scatter(df_random['Occurrences'], df_random['Num_lmers'], color='orange', alpha=0.7, label='Random')
plt.title(f'L-mers Distribution in Random Sequence (l={l})')
plt.xlabel(f'Frequency of {l}-mers')
plt.ylabel(f'Number of distinct {l}-mers')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('lmer_distribution_both.png')
plt.show()

# plot a combined plot for genome and random sequence
plt.figure(figsize=(10, 6))
plt.scatter(df_genome['Occurrences'], df_genome['Num_lmers'], color='blue', label='Genome')
plt.scatter(df_random['Occurrences'], df_random['Num_lmers'], color='red', label='Random')
plt.title(f'L-mers Distribution (l={l})')
plt.xlabel(f'Frequency of {l}-mers')
plt.ylabel(f'Number of distinct {l}-mers')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('lmer_distribution.png')
plt.show()
