from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt

# path to San1 gene of interest
filepath = r"C:\_repos\qPCR-Primer-Designer\SAN1_datasets\ncbi_dataset\data\gene.fna"

# this function reads in a fasta file and returns its records as a list of tuples
def read_fasta(infile: str) -> list[tuple[str, Seq]]:

    # attempts to parse file from given path as .fasta, iteratively appending its records
    try:
        print(f"attempting to pull records from the following path: {infile}")
        records: list[tuple] = []
        for element in SeqIO.parse(infile, 'fasta'):
            ID: str = element.id; sequence: Seq = element.seq
            record: tuple = (ID, sequence)
            records.append(record)
            print(f"record added: {record}")        
        return records

    # exception handling: prints error to terminal and returns empty nested list
    except Exception as error:
        print(f"unable to read file, see exception: {error}")
        return []

# this function generates a list of primer candidates from a parent sequence given minimal and maximal primer lengths
def find_primers(sequence: Seq, min_primer_len: int, max_primer_len: int) -> list[Seq]:
    
    # attempts to iterate through sequence and compile primer candidates as a list
    try:
        primers: list = []
        for i in range(len(sequence)):                                  # iterating through every index in the sequence
            for k in range(min_primer_len, max_primer_len + 1):         # iterating from minimal to maximal primer lengths
                if (i + k <= len(sequence)):                            # ensuring iterator stays within sequence
                    primers.append(sequence[i:i+k])                     # adding primer candidate to the list
        return primers
    
    # exception handling: prints error to terminal and returns empty list
    except Exception as error:
        print(f"unable to determine possible primers, see exception: {error}")
        return []

# this function filters a list of primers by gc_content and returns acceptable candidates as a list
def filter_by_gc(primers: list[Seq], min_gc: float, max_gc: float) -> list[Seq]:
    return [candidate for candidate in primers if (min_gc <= gc_fraction(candidate) <= max_gc)]

# this function filters a list of primers by melting temperature and returns acceptable candidates as a list
def filter_by_mt(primers: list[Seq], min_temp: int, max_temp: int) -> list[Seq]:
    return [candidate for candidate in primers if (min_temp <= mt.Tm_NN(candidate) <= max_temp)]

# this function reads in a list of primers and returns possible primer pairs as a list of tuples
def find_primer_pairs(primers: list[Seq]) -> list[tuple[Seq, Seq]]:
    pairs: list[tuple[Seq, Seq]] = []
    for i, forward_primer in enumerate(primers):
        for reverse_primer in primers[i+1:]:
            pairs.append((forward_primer, reverse_primer))
    return pairs

#this function reads in a list of primer pairs and returns those with compatible melting temperatures
def filter_by_compatible_melting_point(primer_pairs: list[tuple[Seq, Seq]], max_Tm_difference: float) -> list[tuple[Seq, Seq]]:
    return [pair for pair in primer_pairs if abs(mt.Tm_NN(pair[0]) - (mt.Tm_NN(pair[1])) <= max_Tm_difference)]

fasta = read_fasta(filepath)
example_sequence = fasta[0][1]

p1 = find_primers(example_sequence, 23, 23); print(len(p1))
p2 = filter_by_gc(p1, 0.45, 0.55); print(len(p2))
p3 = filter_by_mt(p2, 55, 65); print(len(p3))
p4 = find_primer_pairs(p3); print(len(p4))
p5 = filter_by_compatible_melting_point(p4, .8); print(len(p5))

#for i, seq in enumerate(p4):
#    print(i, seq)