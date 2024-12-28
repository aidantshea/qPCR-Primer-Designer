# imports
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt



# this function generates a list of primer candidates from a parent sequence given minimal and maximal primer lengths
def find_primers(sequence: Seq, min_primer_len: int, max_primer_len: int) -> list[Seq]:
    
    # attempts to iterate through sequence and compile primer candidates as a list
    try:
        primers: list = []
        for i in range(len(sequence)):                                  # iterating through every index in the sequence
            for k in range(min_primer_len, max_primer_len + 1):         # iterating from minimal to maximal primer lengths
                if (i + k <= len(sequence)):                            # ensuring iterator stays within sequence
                    primers.append(sequence[i:i+k])                     # adding primer candidate to the list
        
        print(f"from the sequence of interest, there are {len(primers)} possible primers with lengths between {min_primer_len} and {max_primer_len}.")
        return primers
    
    # exception handling: prints error to terminal and returns empty list
    except Exception as error:
        print(f"unable to determine possible primers, see exception: {error}")
        return []





# this function filters a list of primers by gc_content and returns acceptable candidates as a list
def filter_by_gc(primers: list[Seq], min_gc: float, max_gc: float) -> list[Seq]:
    
    # attempts to select primers with acceptable gc content through a list comprehension
    try:
        candidates: list = [p for p in primers if (min_gc <= gc_fraction(p) <= max_gc)]
        print(f"from {len(primers)} potential primers, {len(candidates)} contain gc content between {min_gc} and {max_gc}.")
        return candidates
    
    # exception handling: prints error to terminal and returns empty list
    except Exception as error:
        print(f"unable to filter by gc_content, see exception: {error}")
        return []




# this function filters a list of primers by melting temperature and returns acceptable candidates as a list
def filter_by_mt(primers: list[Seq], min_temp: int, max_temp: int) -> list[Seq]:
    
    # attempts to select primers with acceptable melting temperature through a list comprehension
    try:
        candidates: list = [p for p in primers if (min_temp <= mt.Tm_NN(p) <= max_temp)]
        print(f"from {len(primers)} potential primers, {len(candidates)} have a Tm beteween {min_temp} and {max_temp}.")
        return candidates
    
    # exception handling: prints error to terminal and returns empty list
    except Exception as error:
        print(f"unable to filter by melting temperature, see exception: {error}")
        return []





# this function reads in a list of primers and returns possible primer pairs as a list of tuples
def find_primer_pairs(primers: list[Seq]) -> list[tuple[Seq, Seq]]:
    
    # exhaustively lists every primer pair for which the reverse primer is downstream of the forward primer
    try:
        pairs: list[tuple[Seq, Seq]] = []
        for i, forward_primer in enumerate(primers):
            for reverse_primer in primers[i+1:]:
                pairs.append((forward_primer, reverse_primer))
        
        print(f"from {len(primers)} suggested primers, there are {len(pairs)} possible pairings.")
        return pairs

    # exception handling: prints error to terminal and returns empty list
    except Exception as error:
        print(f"unable to generate list of primer pairs, see exception: {error}")
        return []





#this function reads in a list of primer pairs and returns those with compatible melting temperatures
def filter_by_compatible_melting_point(primer_pairs: list[tuple[Seq, Seq]], max_Tm_difference: float) -> list[tuple[Seq, Seq]]:
    
    # attempts to select primer pairs with compatible melting points through a list comprehension
    try:
        pairs: list = [p for p in primer_pairs if abs(mt.Tm_NN(p[0]) - (mt.Tm_NN(p[1])) <= max_Tm_difference)]
        print(f"from {len(primer_pairs)} potential primer pairs, {len(pairs)} have melting temperature within {max_Tm_difference} of each other.")
        return pairs

    # exception handling: prints error to terminal and returns empty list
    except Exception as error:
        print(f"unable to filter primer pairs by Tm compatibility, see exception: {error}")
        return []