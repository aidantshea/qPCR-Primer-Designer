from Bio import SeqIO

# path to San1 gene of interest
filepath = r"C:\_repos\qPCR-Primer-Designer\SAN1_datasets\ncbi_dataset\data\gene.fna"

# this function reads in a fasta file and returns its records as a nested list
def read_fasta(infile: str) -> list[list[str]]:

    # attempts to parse file from given path as .fasta, iteratively appending its records
    try:
        print(f"attempting to pull records from the following path: {infile}")
        records: list = []        
        for element in SeqIO.parse(infile, 'fasta'):
            ID: str = element.id; Seq: str = element.seq
            record: list = [ID, Seq]
            records.append(record)
            print(f"record added: {record}")        
        return records

    # exception handling: prints error to terminal and returns empty nested list
    except Exception as error:
        print(f"unable to read file, see exception: {error}")
        return [[]]

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

sequence = read_fasta(filepath)[0][1]
primers = find_primers(sequence, 25, 27); print(primers[0])