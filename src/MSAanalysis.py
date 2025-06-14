from Bio import AlignIO
from collections import Counter

# Calculate conservation scores
def get_conservation_scores(alignment):
    scores = []
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        counts = Counter(column.replace("-", ""))
        if counts:
            most_common_freq = counts.most_common(1)[0][1] / len(column.replace("-", ""))
        else:
            most_common_freq = 0
        scores.append(most_common_freq)
    return scores

# Find ranges of conserved regions
def get_conserved_ranges(conserved_positions, min_length=10):
    if not conserved_positions:
        return []

    ranges = []
    start = prev = conserved_positions[0]

    for pos in conserved_positions[1:]:
        if pos == prev + 1:
            prev = pos
        else:
            if (prev - start + 1) >= min_length:
                ranges.append((start, prev))
            start = prev = pos
    # Add the final region if it meets length requirement
    if (prev - start + 1) >= min_length:
        ranges.append((start, prev))

    return ranges

def main():
    # Load MSA
    alignment = AlignIO.read("data/output/MSA.fasta", "fasta")
    scores = get_conservation_scores(alignment)
    threshold = 1
    conserved_positions = [i for i, s in enumerate(scores) if s == threshold]
    ranges = get_conserved_ranges(conserved_positions)
    with open("data/output/conserved_ranges.txt", "w") as f:
        for r in ranges:
            f.write(f"{r[0]}-{r[1]}\n")


if __name__ == "__main__":
    main()