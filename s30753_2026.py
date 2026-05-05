# student number: s30753

# Random DNA sequence generator in FASTA format. The program supports statistics, name embedding,
# and extra features such as batch mode, configurable nucleotide distribution, motif search,
# sliding window, GC analysis, and FASTA file validator

import random
import os

# Returns a random DNA sequence of a length given as argument. Randomly chooses between A, C, G, T and
#joins the whole sequence together
def generate_sequence(length: int) -> str:
    nucleotides = ['A', 'C', 'G', 'T']
    return ''.join(random.choices(nucleotides, k=length))


# calculates the percentage of occurrence of each of the nucleotides as well as the percentage of GC
#it filters out lowercase letters
def calculate_stats(sequence: str) -> dict:
    bio_seq = ''.join(ch for ch in sequence if ch.isupper())
    n = len(bio_seq)

    statistics = {}
    for nucleotide in ['A', 'C', 'G', 'T']:
        statistics[nucleotide] = (bio_seq.count(nucleotide) / n) * 100
    statistics['GC'] = statistics['G'] + statistics['C']

    return statistics

# inserts a name at a random position in the sequence. The name is lowercase so that not to confuse and
#mix with the nucleotides
def insert_name(sequence: str, name: str) -> str:
    position = random.randint(0, len(sequence))
    return sequence[:position] + name.lower() + sequence[position:]

# returns a formatted FASTA record as a string, based on given parameters. Header starting with >
#and including id and potential description. The chunks of the sequence are sliced fixed-size
def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    if description:
        header = f">{seq_id} {description}"
    else:
        header = f">{seq_id}"

    lines = [sequence[i:i + line_width] for i in range(0, len(sequence), line_width)]

    return header + '\n' + '\n'.join(lines) + '\n'

# validates whether the provided number is in the acceptable range
def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    while True:
        user_input = input(prompt)
        try:
            value = int(user_input)
            if min_val <= value <= max_val:
                return value
            else:
                print(f"Error: value must be an integer in the range [{min_val}, {max_val}]")
        except ValueError:
            print(f"Error: value must be an integer in the range [{min_val}, {max_val}]")

# generates multiple DNA sequences, concatenates them and puts them together as a multi-FASTA string
def batch_generate(seq_length: int, count: int, base_id: str, description: str) -> str:
    records = []
    for i in range(1, count + 1):
        seq_id = f"{base_id}_{i:03d}"
        sequence = generate_sequence(seq_length)
        records.append(format_fasta(seq_id, description, sequence))
    return ''.join(records)

# the user is asked to input valid percentages per each nucleotide. The inputs are validated
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

# generates a DNA sequence with nucleotide probabilities defined by the user
def generate_sequence_weighted(length: int, weights: dict) -> str:
    nucleotides = ['A', 'C', 'G', 'T']
    w = [weights[n] for n in nucleotides]
    return ''.join(random.choices(nucleotides, weights=w, k=length))

# searches for all occurrences of a motif in the sequence and returns a list of start positions of motifs.
#Lowercase characters are ignored
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

# calculates GC content in a sliding window along the sequence and returns a list of tuples with
#starting position and gc percentage
def sliding_window_gc(sequence: str, window_size: int) -> list:
    bio_seq = ''.join(ch for ch in sequence if ch.isupper())
    results = []
    for i in range(len(bio_seq) - window_size + 1):
        window = bio_seq[i:i + window_size]
        gc = (window.count('G') + window.count('C')) / window_size * 100
        results.append((i + 1, round(gc, 2)))
    return results

# saves sliding window GC analysis results to a CSV file with columns being the starting position and GC
#content (just like the values in the tuples)
def save_sliding_window_csv(results: list, filename: str) -> None:
    with open(filename, 'w') as f:
        f.write("start_position,gc_content\n")
        for start, gc in results:
            f.write(f"{start},{gc}\n")

# checks if a given file is a valid FASTA file
def validate_fasta_file(filepath: str) -> list:
    #set of valid characters
    valid_chars = set('ACGTRYSWKMBDHVNacgtryswkmbdhvn')
    errors = []

    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        return [f"Error: file '{filepath}' not found"]

    if not lines:
        return ["Error: file is empty"]

    #checks whether the header starts with >
    if not lines[0].startswith('>'):
        errors.append("Line 1: file must start with a header line ('>')")

    previous_was_header = False

    for i, line in enumerate(lines, start=1):
        line = line.rstrip('\n')
        if line.startswith('>'):
            #there cannot be header after header
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

# main() showcases the implemented functionalities. It gets user input, generates a DNA sequence,
#add the user's name, saves the FASTA file, prints statistics of nucleotide occurrences and runs
#provided extra features
def main():

    use_custom = input("Use custom nucleotide distribution? (y/n): ").strip().lower()
    weights = None
    if use_custom.lower() == 'y':
        weights = get_nucleotide_weights()

    length = validate_positive_int("Enter sequence length: ")

    while True:
        seq_id = input("Enter sequence ID: ").strip()
        if not seq_id:
            print("Error: ID cannot be empty.")
        elif any(ch.isspace() for ch in seq_id):
            print("Error: ID cannot contain whitespace.")
        else:
            break

    description = input("Enter a description of the sequence (optional) or click Enter to skip: ").strip()

    name = input("Enter your name: ").strip()

    if weights:
        sequence = generate_sequence_weighted(length, weights)
    else:
        sequence = generate_sequence(length)

    sequence_with_name = insert_name(sequence, name)

    filename = f"{seq_id}.fasta"
    fasta_content = format_fasta(seq_id, description, sequence_with_name)
    with open(filename, 'w') as f:
        f.write(fasta_content)
    print(f"\nSequence saved to {filename}")

    statistics = calculate_stats(sequence_with_name)
    print(f"\nSequence statistics (n={length}):")
    for nuc in ['A', 'C', 'G', 'T']:
        print(f"  {nuc}: {statistics[nuc]:.2f}%")
    print(f"  GC-content: {statistics['GC']:.2f}%")

    do_motif = input("\nSearch for a motif? (y/n): ").strip().lower()
    if do_motif.lower() == 'y':
        motif = input("Enter motif (e.g. ATG): ").strip()
        positions = find_motif(sequence_with_name, motif)
        if positions:
            print(f"Motif '{motif.upper()}' found at positions: {positions}")
        else:
            print(f"Motif '{motif.upper()}' not found in sequence")

    do_window = input("\nRun sliding window GC analysis? (y/n): ").strip().lower()
    if do_window.lower() == 'y':
        window_size = validate_positive_int("Enter window size (nt): ", min_val=1, max_val=length)
        sw_results = sliding_window_gc(sequence_with_name, window_size)
        csv_filename = f"{seq_id}_gc_window.csv"
        save_sliding_window_csv(sw_results, csv_filename)
        print(f"Sliding window GC results saved to: {csv_filename}")

    do_batch = input("\nRun batch mode (generate multiple sequences)? (y/n): ").strip().lower()
    if do_batch.lower() == 'y':
        count = validate_positive_int("How many sequences to generate? ", min_val=1, max_val=100)
        batch_id = input("Enter base ID for batch (e.g. Seq): ").strip()
        batch_desc = input("Enter description for batch sequences: ").strip()
        multi_fasta = batch_generate(length, count, batch_id, batch_desc)
        batch_filename = f"{batch_id}_batch.fasta"
        with open(batch_filename, 'w') as f:
            f.write(multi_fasta)
        print(f"Batch of {count} sequences saved to: {batch_filename}")

    do_validate = input("\nValidate an existing FASTA file? (y/n): ").strip().lower()
    if do_validate.lower() == 'y':
        fasta_path = input("Enter path to FASTA file: ").strip()
        errors = validate_fasta_file(fasta_path)
        if not errors:
            print("File has a valid FASTA format")
        else:
            print(f"Found {len(errors)} error(s):")
            for err in errors:
                print(f"  - {err}")


if __name__ == "__main__":
    main()