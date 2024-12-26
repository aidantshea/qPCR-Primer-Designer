from Bio import SeqIO

# path to San1 gene of interest
filepath = r"C:\_repos\qPCR-Primer-Designer\SAN1_datasets\ncbi_dataset\data\gene.fna"

# this method reads in a fasta file and returns its records as a nested list
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

    # exception handling: prints error to terminal and returns None
    except Exception as error:
        print(f"unable to read file, see exception: {error}")
        return [[]]