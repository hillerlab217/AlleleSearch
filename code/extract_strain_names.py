import os
import csv

goldenset_dir = "/home/xinruoz/Document/Hiller/data/GoldenSet"
output_csv = "/home/xinruoz/Document/Hiller/strain_names.csv"

strain_names = sorted(
    f[:-6]  # remove ".fasta"
    for f in os.listdir(goldenset_dir)
    if f.endswith(".fasta")
)

with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["strain_name"])
    for name in strain_names:
        writer.writerow([name])

print(f"Written {len(strain_names)} strain names to {output_csv}")
