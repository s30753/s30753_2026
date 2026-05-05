import random
import os

def generate_sequence(length: int) -> str:
    nucleotides = ['A', 'C', 'G', 'T']
    return ''.join(random.choices(nucleotides, k=length))


def calculate_stats(sequence: str) -> dict:
    bio_seq = ''.join(ch for ch in sequence if ch.isupper())
    n = len(bio_seq)

    stats = {}
    for nucleotide in ['A', 'C', 'G', 'T']:
        stats[nucleotide] = (bio_seq.count(nucleotide) / n) * 100
    stats['GC'] = stats['G'] + stats['C']

    return stats


def insert_name(sequence: str, name: str) -> str:
    position = random.randint(0, len(sequence))
    return sequence[:position] + name.lower() + sequence[position:]


def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    if description:
        header = f">{seq_id} {description}"
    else:
        header = f">{seq_id}"

    lines = [sequence[i:i + line_width] for i in range(0, len(sequence), line_width)]

    return header + '\n' + '\n'.join(lines) + '\n'


def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    while True:
        user_input = input(prompt)
        try:
            value = int(user_input)
            if min_val <= value <= max_val:
                return value
            else:
                print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")
        except ValueError:
            print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")

def batch_generate(seq_length: int, count: int, base_id: str, description: str) -> str:
    records = []
    for i in range(1, count + 1):
        seq_id = f"{base_id}_{i:03d}"
        sequence = generate_sequence(seq_length)
        records.append(format_fasta(seq_id, description, sequence))
    return ''.join(records)

def get_nucleotide_weights() -> dict:
    while True:
        weights = {}
        print("Enter the percentage for each nucleotide (must sum to 100):")
        valid = True
        for nuc in ['A', 'C', 'G', 'T']:
            raw = input(f"  {nuc}: ")
            try:
                val = float(raw)
                if val < 0:
                    print("Error: percentages must be non-negative")
                    valid = False
                    break
                weights[nuc] = val
            except ValueError:
                print(f"Error: '{raw}' is not a valid number")
                valid = False
                break
        if valid:
            total = sum(weights.values())
            if abs(total - 100.0) < 1e-6:
                return weights
            else:
                print(f"Error: percentages sum to {total:.2f}, must be 100.00")


def generate_sequence_weighted(length: int, weights: dict) -> str:
    nucleotides = ['A', 'C', 'G', 'T']
    w = [weights[n] for n in nucleotides]
    return ''.join(random.choices(nucleotides, weights=w, k=length))


def find_motif(sequence: str, motif: str) -> list:
    bio_seq = ''.join(c for c in sequence if c.isupper())
    motif = motif.upper()
    positions = []
    start = 0
    while start <= len(bio_seq) - len(motif):
        idx = bio_seq.find(motif, start)
        if idx == -1:
            break
        positions.append(idx + 1)
        start = idx + 1
    return positions

def sliding_window_gc(sequence: str, window_size: int) -> list:
    bio_seq = ''.join(ch for ch in sequence if ch.isupper())
    results = []
    for i in range(len(bio_seq) - window_size + 1):
        window = bio_seq[i:i + window_size]
        gc = (window.count('G') + window.count('C')) / window_size * 100
        results.append((i + 1, round(gc, 2)))
    return results


def save_sliding_window_csv(results: list, filename: str) -> None:
    with open(filename, 'w') as f:
        f.write("start_position,gc_content\n")
        for start, gc in results:
            f.write(f"{start},{gc}\n")


def validate_fasta_file(filepath: str) -> list:
    valid_chars = set('ACGTRYSWKMBDHVNacgtryswkmbdhvn')
    errors = []

    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        return [f"Error: file '{filepath}' not found."]

    if not lines:
        return ["Error: file is empty."]

    if not lines[0].startswith('>'):
        errors.append("Line 1: file must start with a header line ('>')")

    previous_was_header = False

    for i, line in enumerate(lines, start=1):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if previous_was_header:
                errors.append(f"Line {i}: consecutive header lines with no sequence between them")
            previous_was_header = True
        else:
            previous_was_header = False
            invalid = set(line) - valid_chars
            if invalid:
                errors.append(f"Line {i}: invalid characters found: {sorted(invalid)}")
            if len(line) > 80:
                errors.append(f"Line {i}: line too long ({len(line)} chars, max 80).")

    return errors