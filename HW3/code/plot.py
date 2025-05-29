import matplotlib.pyplot as plt

# Data from the table
datasets = [1, 2, 3, 4, 5]
opt_score = [5284, 18992, 3427, 39908, 55591]
align_len = [592, 1785, 2793, 4895, 6457]
perfect_match = [425, 1505, 1417, 3371, 4522]
runtime_sec = [2.553, 62.332, 150.935, 1502.2, 6179.08]
peak_mem = [10.5664, 50.3984, 83.2695, 348.867, 610.219]

plt.figure(figsize=(10, 6))

# Plot each column; matplotlib automatically cycles colors
plt.plot(datasets, opt_score, marker='o', label='Optimal Score')
plt.plot(datasets, align_len, marker='o', label='Alignment Length')
plt.plot(datasets, perfect_match, marker='o', label='Perfect Matches')
plt.plot(datasets, runtime_sec, marker='o', label='Runtime (s)')
plt.plot(datasets, peak_mem, marker='o', label='Peak Memory (MiB)')


plt.xticks(datasets, [f'Test Case #{i}' for i in datasets])
plt.yscale('log')  # Log scale to fit all magnitudes nicely
plt.xlabel('Dataset (Test Case #)')
plt.ylabel('Value (log scale)')
plt.title('Alignment Metrics Across Test Cases')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()
plt.savefig('alignment_metrics.png')  # Save the figure