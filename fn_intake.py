#imports
from Bio import SeqIO
from Bio.Seq import Seq



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