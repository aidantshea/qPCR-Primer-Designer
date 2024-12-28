from fn_intake import read_fasta
from fn_primer_pipes import find_primers, filter_by_gc, filter_by_mt, find_primer_pairs, filter_by_compatible_melting_point

# path to San1 gene of interest
filepath = r"C:\_repos\qPCR-Primer-Designer\SAN1_datasets\ncbi_dataset\data\gene.fna"

fasta = read_fasta(filepath)
example_sequence = fasta[0][1]

p1 = find_primers(example_sequence, 23, 23)
p2 = filter_by_gc(p1, 0.45, 0.55)
p3 = filter_by_mt(p2, 55, 65)
p4 = find_primer_pairs(p3)
p5 = filter_by_compatible_melting_point(p4, .8)