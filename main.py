from Bio import SeqIO

# path to San1 gene of interest
filepath = r"C:\_repos\qPCR Primer Algorithm\SAN1_datasets\ncbi_dataset\data\gene.fna"

# this method reads in a fasta file and returns its records as a nested list
def read_fasta(infile):

    # attempts to parse file from given path as .fasta, iteratively appending its records
    try:
        records: list = []
        for element in SeqIO.parse(infile, 'fasta'):
            ID: str = element.id; Seq: str = element.seq
            records.append([ID, Seq])
            print(f"record added: {[ID, Seq]}")
        return records
    
    # exception handling: prints error to terminal and returns None
    except Exception as error:
        print(f"unable to read file, see exception: {error}")
        return None

sequences = read_fasta(filepath)