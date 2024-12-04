import sys
import os
from os import path
import argparse
import time

# Translate a fasta alignment into a MAPLE format file.

class Sequence:
    '''A sequence parsed from a FASTA file.

    Attributes:
    description: str
        The description of the genome.
    sequence: str
        The genome sequence.
    '''
    def __init__(self, description: str, sequence: str):
        self.name = description
        self.sequence = sequence


def parse_fasta(file_path: str) -> list:
    '''Parse a fasta file and return a list of Sequence objects. '''
    sequences = []
    current_label = None
    current_sequence = ""
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            # Check if the line is a label
            if line.startswith('>'):
                # Flush previous label
                if current_label is not None:
                    sequences.append(Sequence(current_label, current_sequence))

                # Start a new label and reset sequence
                current_label = line[1:]  # Remove the '>'
                current_sequence = ""
            else:
                # Append the sequence lines
                current_sequence += line

        # Save the last label and sequence
        if current_label is not None:
            sequences.append(Sequence(current_label, current_sequence.lower()))
        
        return sequences

def generate_consensus(sequences: list) -> str:
    '''Generate a consensus sequence from a list of sequences. '''
    # Get the length of the sequences
    seq_length = len(sequences[0].sequence)

    for s in sequences:
        if len(s.sequence) != seq_length:
            raise ValueError("All sequences must be the same length.")
    
    consensus = ""
    for i in range(seq_length):
        # Get the nucleotides at position i
        nucleotides = [s.sequence[i] for s in sequences]
        # Count the number of each nucleotide
        counts = {n: nucleotides.count(n) for n in 'acgt-n'}
        # Get the most common nucleotide
        consensus += max(counts, key=counts.get)
    
    print(f"Consensus sequence: {consensus}")
    return consensus

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Translate a fasta alignment into a MAPLE format file.')
    parser.add_argument('--input', type=str, help='The input fasta file.')
    parser.add_argument('--reference', type=str, help='The reference fasta file.')
    parser.add_argument("--output", type=str, help="The output MAPLE file.")
    args = parser.parse_args()

    import time
    start_time = time.time()

    if not path.exists(args.input):
        print(f"ERROR: Input file {args.input} not found.")
        sys.exit(1)
    
    if args.reference:
        if not path.exists(args.reference):
            print(f"ERROR: Reference file {args.reference} not found.")
            sys.exit(1)
    else:
        print("No reference file provided. Consensus sequence will be used as reference.")
    
    reference = None
    sequences = parse_fasta(args.input)
    if args.reference:
        reference = parse_fasta(args.reference)[0].sequence
    else:
        reference = generate_consensus(sequences)
    # Check if the reference sequence length matches the input sequences
    for s in sequences:
        if len(reference) != len(s.sequence):
            print(f"ERROR: Reference sequence length does not match the input sequence {s.name}.")
            sys.exit(1)
    # Write the output MAPLE file
    with open(args.output, 'w') as file:
        file.write(f">reference\n{reference.lower()}\n")

        for s in sequences:
            file.write(f">{s.name}\n")
            diff = ""
            for i in range(len(s.sequence)):
                if s.sequence[i] == 'n' or s.sequence[i] == '-':
                    diff += s.sequence[i]
                elif s.sequence[i] != reference[i]:
                    diff += s.sequence[i]
                else:
                    diff += "0"
            
            
            last_n_or_gap = None
            last_index = -1
            for i in range(len(diff)):
                if diff[i] == 'n':
                    if last_n_or_gap != 'n':
                        # Flush the last gap
                        if last_n_or_gap == '-':
                            file.write(f"-\t{last_index}\t{i - last_index + 1}\n")
                        last_n_or_gap = 'n'
                        last_index = i + 1
                    else:
                        continue
                elif diff[i] == '-':
                    if last_n_or_gap != '-':
                        # Flush the last n
                        if last_n_or_gap == 'n':
                            file.write(f"n\t{last_index}\t{i - last_index + 1}\n")
                        last_n_or_gap = '-'
                        last_index = i + 1
                    else:
                        continue
                elif diff[i] != '0':
                    # Flush the last gap or n
                    if last_n_or_gap == 'n':
                        file.write(f"n\t{last_index}\t{i - last_index + 1}\n")
                    elif last_n_or_gap == '-':
                        file.write(f"n\t{last_index}\t{i - last_index + 1}\n")

                    last_n_or_gap = None
                    last_index = -1

                    file.write(f"{diff[i].lower()}\t{i + 1}\n")
                else:
                    # Flush the last gap or n
                    if last_n_or_gap == 'n':
                        file.write(f"n\t{last_index}\t{i - last_index + 1}\n")
                    elif last_n_or_gap == '-':
                        file.write(f"n\t{last_index}\t{i - last_index + 1}\n")

                    last_n_or_gap = None
                    last_index = -1
                
            # Flush the last gap or n
            if last_n_or_gap == 'n':
                file.write(f"n {last_index} {i - last_index + 1}\n")
            elif last_n_or_gap == '-':
                file.write(f"n {last_index} {i - last_index + 1}\n")

            last_n_or_gap = None
            last_index = -1
    
    print("Time - %s seconds" % (time.time() - start_time))

            

                    
    

        
