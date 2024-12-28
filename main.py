from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt

# path to San1 gene of interest
filepath = r"C:\_repos\qPCR-Primer-Designer\SAN1_datasets\ncbi_dataset\data\gene.fna"

# this function reads in a fasta file and returns its records as a list of tuples
def read_fasta(infile: str) -> list[tuple]:

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
def find_primers(sequence: str, min: int, max: int) -> list[str]:
    
    # attempts to iterate through sequence and compile primer candidates as a list
    try:
        primers: list = []
        for i in range(len(sequence)):                  # iterating through every index in the sequence
            for k in range(min, max + 1):               # iterating from minimal to maximal primer lengths
                if (i + k <= len(sequence)):            # ensuring iterator stays within sequence
                    primers.append(sequence[i:i+k])     # adding primer candidate to the list
        return primers
    
    # exception handling: prints error to terminal and returns empty list
    except Exception as error:
        print(f"unable to determine possible primers, see exception: {error}")
        return []

# this function filters a list of primers by gc_content and returns acceptable candidates as a list
def filter_by_gc(primers: list[str], min_gc: float, max_gc: float) -> list[str]:
    return [candidate for candidate in primers if (min_gc <= gc_fraction(candidate) <= max_gc)]

# this function filters a list of primers by melting temperature and returns acceptable candidates as a list
def filter_by_mt(primers: list[str], min_temp: int, max_temp: int) -> list[str]:
    return [candidate for candidate in primers if (min_temp <= mt.Tm_NN(candidate) <= max_temp)]

# this function filters a list of primers for those with reverse primers of similar melting temperature
def filter_compatible_pairs(primers: list[str], space: int = 25) -> list[list[Seq]]:
    pairs: list[list[Seq]] = []
    for i, primer in enumerate(primers[0:10]):
        for reverse_primer in primers[i+space:]:
            sequence = Seq(primer)
            reverse_sequence = Seq(primer).reverse_complement()
            if abs(mt.Tm_NN(sequence) - mt.Tm_NN(reverse_sequence) <= 2):
                pairs.append([sequence, reverse_sequence])
    return pairs

fasta = read_fasta(filepath)
example_sequence = fasta[0][1]; print(example_sequence)

p1 = find_primers(example_sequence, 17, 24); print(len(p1))
p2 = filter_by_gc(p1, 0.4, 0.6); print(len(p2))
p3 = filter_by_mt(p2, 55, 65); print(len(p3))
p4 = filter_compatible_pairs(p3); print(len(p4))