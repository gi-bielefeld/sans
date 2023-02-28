"""
This file implements synthetic sequence generation
for testing including the true splits and amino acids
"""
import math
import random
from copy import deepcopy

from .console import log_print
from . import path
from .container import Container


# The DNA bases
BASES = 'A,C,G,T'.split(',')
PAC_BASES = 'A,C,G,T'.split(',')

# Map base : base reverse complements
COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


ORDER = {'A': {'C', 'G', 'T'}, 'C': {'G', 'T'}, 'G': {'T'}, 'T': {}}


# THe amino acid codes
AMINO_CODES = 'A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V'.split(",")

# Map code : [triplets]
TRIPLETS_BY_CODE = {'A': ['GCA', 'GCC', 'GCG', 'GCT'], 'R': ['CGA', 'CGC', 'CGG', 'CGT', 'AGA', 'AGG'],
                    'N': ['AAC', 'AAT'], 'D': ['GAC', 'GAT'], 'C': ['TGC', 'TGT'], 'Q': ['CAA', 'CAG'],
                    'E': ['GAA', 'GAG'], 'G': ['GGA', 'GGC', 'GGG', 'GGT'], 'H': ['CAC', 'CAT'],
                    'I': ['ATA', 'ATC', 'ATT'], 'L': ['CTA', 'CTC', 'CTG', 'CTT'], 'K': ['AAA', 'AAG'],
                    'M': ['ATG'], 'F': ['TTC', 'TTT'], 'P': ['CCA', 'CCC', 'CCG', 'CCT'],
                    'S': ['TCA', 'TCC', 'TCG', 'TCT'], 'T': ['ACA', 'ACC', 'ACG', 'ACT'], 'W': ['TGG'],
                    'Y': ['TAC', 'TAT'], 'V': ['GTA', 'GTC', 'GTG', 'GTT']}

# Map triplet : [code]
CODES_BY_TRIPLETS = dict()
for item in TRIPLETS_BY_CODE.items():
    amino, triplets = item
    for triplet in triplets:
        CODES_BY_TRIPLETS[triplet] = amino

# Amino triplets
TRIPLETS = list(CODES_BY_TRIPLETS.keys())


def __geom2(x, y):
    return round(math.sqrt(x + 1) + math.sqrt(y + 1), 5)


def __minimize_kmer(kmer, k):
    rc_mer = []
    for element in kmer:
        rc_mer.insert(0, COMPLEMENTS[element])
    for i in range(k):
        if kmer[i] in ORDER[rc_mer[i]]:
            return rc_mer
        elif rc_mer[i] in ORDER[kmer[i]]:
            return kmer
    return kmer


def __invert_colors(color, seq_num):
    inverted_color = set()
    for i in range(seq_num):
        if i not in color:
            inverted_color.add(i)
    return inverted_color


def make_synthetic_truth(container: Container):
    seq_path = container.seq_dir.fa_path
    split_path = container.seq_dir.set_path
    seq_num = container.data.seq_num # Number of sequences
    seq_len = container.data.seq_len # Length of sequences

    seqs = []
    amino_seqs = []

    # If the container requires a base truth, the sequences are created accordingly
    log_print("[SEQ_MAKE] Building test data")
    # Make the length a multiple of three for triplet usage
    seq_len = seq_len - (seq_len % 3)       # The length of the sequence as a multiple of three
    amino_len = int(seq_len / 3)            # The length of the amino acid sequence

    # --- [Sequence generation] ---
    # Create a master seq, rc_master seq and amino master seq
    master_seq = []
    rc_master_seq = []
    master_amino_seq = []
    # Create a master sequence
    for base_it in range(0, seq_len, 3):
        triplet = random.choice(TRIPLETS) # Choose a random triplet
        aminoacid = CODES_BY_TRIPLETS[triplet] # Translate
        master_seq.extend(triplet)
        master_amino_seq.append(aminoacid)
    # Compute the reverse complement master seq
    for base in master_seq:
        rc_master_seq.insert(0, COMPLEMENTS[base])

    seqs.append(master_seq)
    seqs.append(rc_master_seq)
    amino_seqs.append(master_amino_seq)
    amino_seqs.append(master_amino_seq)
    for i in range(seq_num - 2):
        seqs.append(deepcopy(master_seq))
        amino_seqs.append(deepcopy(master_amino_seq))

    # --- [Sequence Mutation] ---
    # The distance between two split generating SNPs in the data set
    snp_distance = 3 * container.sans.k
    # The Master seq mutation rate
    m_rate = container.data.m_rate / 100
    # The number of split generating mutations
    mutations = min([int(seq_len * m_rate), int((seq_len - 2 * snp_distance) / snp_distance)])
    true_mutations = 0
    # --- Create mutations and splits ---
    seq_ids = [i for i in range(2, seq_num)]
    for pos in range(snp_distance, seq_len - snp_distance, snp_distance):
        # If the target number of mutations is reached, stop mutation
        if true_mutations == mutations:
            break
        true_mutations += 1

        # The position in the amino sequence
        acid_pos = int((pos / 3) - 1)

        # The number of sequences to mutate
        split_size = random.randint(2, seq_num - 2)
        # The sequences to mutate
        target_ids = random.sample(seq_ids, split_size)

        # The possible bases
        bases = [element for element in BASES if element != master_seq[pos]]

        # Create a new random coding triplet and a mutated base
        new_triplet = ""
        m_base = ""
        while new_triplet not in CODES_BY_TRIPLETS and bases:
            m_base = random.choice(bases)
            new_triplet = "".join(master_seq[pos - 1:pos + 2])
            bases.remove(m_base)
        if new_triplet not in CODES_BY_TRIPLETS:
            continue

        # Apply mutation to all target sequences
        new_code = CODES_BY_TRIPLETS[new_triplet]
        for seq_id in target_ids:
            # Create DNA mutation
            seqs[seq_id][pos - 1] = m_base
            # Create AMINO mutation
            amino_seqs[seq_id][acid_pos] = new_code

    # --- [Write Sequences] ---
    # DNA
    dna_files = [f"{i}.fa" for i in range(seq_num - 1)]
    dna_files.insert(1, f"0r.fa")
    dna_paths = [path.nodes_to_file_path([seq_path, dna_files[i]]) for i in range(seq_num)]
    dna_contents = [f">S{i}\n{''.join(seqs[i])}" for i in range(seq_num)]

    for i in range(seq_num):
        fasta_file = open(dna_paths[i], "w+")
        fasta_file.write(dna_contents[i])
        fasta_file.close()

    # Write lists
    rev_list_content = "\n".join(dna_files)
    rev_list_file = "list_dna.txt"
    rev_list_path = path.nodes_to_file_path([seq_path, rev_list_file])
    rev_list = open(rev_list_path, 'w+')
    rev_list.write(rev_list_content)

    container.seq_dir.dna_list_path = rev_list_path

    n_rev_list_content = "\n".join([file for file in dna_files if file != "0r.fa"])
    n_rev_list_file = "list_dna_nrev.txt"
    n_rev_list_path = path.nodes_to_file_path([seq_path, n_rev_list_file])
    n_rev_list = open(n_rev_list_path,  'w+')
    n_rev_list.write(n_rev_list_content)
    n_rev_list.close()

    container.seq_dir.dna_nrev_list_path = n_rev_list_path


    # AMINO
    amino_files = [f"AMINO_{i}.fa" for i in range(seq_num)]
    amino_paths = [path.nodes_to_file_path([seq_path, amino_files[i]]) for i in range(seq_num)]
    amino_contents = [f">A{i}\n{''.join(amino_seqs[i])}" for i in range(seq_num)]

    for i in range(seq_num):
        fasta_file = open(amino_paths[i], "w+")
        fasta_file.write(amino_contents[i])
        fasta_file.close()

    amino_list_content = "\n".join([file for file in amino_files])
    amino_list_file = "list_amino.txt"
    amino_list_path = path.nodes_to_file_path([seq_path, amino_list_file])
    amino_list = open(amino_list_path, "w+")
    amino_list.write(amino_list_content)
    amino_list.close()

    container.seq_dir.amino_list_path = amino_list_path

    # --- [Split computation] ---
    # DNA
    # COMPUTE KMER MAPS
    nrev_dna_kmer_map = dict()
    rev_dna_kmer_map = dict()
    for i in range(len(seqs)):
        sequence = seqs[i]
        for pos in range(len(sequence) - container.sans.k + 1):
            kmer = sequence[pos:pos+container.sans.k]
            kmer_string = "".join(kmer)

            # NREV
            if kmer_string not in nrev_dna_kmer_map:
                nrev_dna_kmer_map[kmer_string] = {i}
            else:
                nrev_dna_kmer_map[kmer_string].add(i)

            # REV
            canonical = __minimize_kmer(kmer, container.sans.k)
            canonical_string = "".join(canonical)
            if canonical_string not in rev_dna_kmer_map:
                rev_dna_kmer_map[canonical_string] = {i}
            else:
                rev_dna_kmer_map[canonical_string].add(i)

    # Aminio
    amino_kmer_map = dict()
    for i in range(seq_num):
        sequence = amino_seqs[i]
        for pos in range(len(sequence) - container.sans.k + 1):
            kmer = sequence[pos:pos+container.sans.k]
            kmer_string = "".join(kmer)
            if kmer_string not in amino_kmer_map:
                amino_kmer_map[kmer_string] = {i}
            else:
                amino_kmer_map[kmer_string].add(i)

    # COMPUTE COLOR MAPS
    # No reverse
    nrev_dna_color_map = dict()
    for item in nrev_dna_kmer_map.items():
        kmer, colors = item
        inverse_colors = __invert_colors(colors, seq_num)
        ping_pong = 0
        if len(colors) < len(inverse_colors):
            if not colors:
                continue
            color_string = ",".join([str(color) for color in colors])
        else:
            if not inverse_colors:
                continue
            color_string = ",".join([str(color) for color in inverse_colors])
            ping_pong = 1

        if color_string not in nrev_dna_color_map:
            nrev_dna_color_map[color_string] = [0, 0]
        nrev_dna_color_map[color_string][ping_pong] += 1

    # Reverse
    rev_dna_color_map = dict()
    for item in rev_dna_kmer_map.items():
        kmer, colors = item
        inverse_colors = __invert_colors(colors, seq_num)
        ping_pong = 0

        if len(colors) < len(inverse_colors):
            if not colors:
                continue
            color_string = ",".join([str(color) for color in colors])
        else:
            if not inverse_colors:
                continue
            color_string = ",".join([str(color) for color in inverse_colors])
            ping_pong = 1

        if color_string not in rev_dna_color_map:
            rev_dna_color_map[color_string] = [0, 0]
        rev_dna_color_map[color_string][ping_pong] += 1

    # Amino
    amino_color_map = dict()
    for item in amino_kmer_map.items():
        kmer, colors = item
        inverse_colors = __invert_colors(colors, seq_num)
        ping_pong = 0
        if len(colors) < len(inverse_colors):
            if not colors:
                continue
            color_string = ",".join([str(color) for color in colors])
        else:
            if not inverse_colors:
                continue
            color_string = ",".join([str(color) for color in inverse_colors])
            ping_pong = 1

        if color_string not in amino_color_map:
            amino_color_map[color_string] = [0, 0]
        amino_color_map[color_string][ping_pong] += 1

    # COMPUTE SPLITS
    # NREV
    nrev_splits = dict()
    for item in nrev_dna_color_map.items():
        split, counter = item
        nrev_splits[split] = __geom2(*counter)

    # REV
    rev_splits = dict()
    for item in rev_dna_color_map.items():
        split, counter = item
        rev_splits[split] = __geom2(*counter)

    # AMINO
    amino_splits = dict()
    for item in amino_color_map.items():
        split, counter = item
        amino_splits[split] = __geom2(*counter)

    # --- [Write Splits] ---
    # NREV
    nrev_splits_content = []
    for split in nrev_splits.items():
        color_string, weight = split
        colors = [int(color) for color in color_string.split(",")]
        split_line = [str(weight)]
        split_line.extend([dna_files[color] for color in colors])
        split_line = "\t".join(split_line)
        nrev_splits_content.append(split_line)
    nrev_splits_content = "\n".join(nrev_splits_content)

    nrev_splits_path = path.nodes_to_file_path([split_path, "dna_nrev_splits.txt"])
    nrev_splits_file = open(nrev_splits_path, "w+")
    nrev_splits_file.write(nrev_splits_content)
    nrev_splits_file.close()

    container.seq_dir.dna_nrev_splits_path = nrev_splits_path

    # REV
    rev_splits_content = []
    for split in rev_splits.items():
        color_string, weight = split
        colors = [int(color) for color in color_string.split(",")]
        split_line = [str(weight)]
        split_line.extend([dna_files[color] for color in colors])
        split_line = "\t".join(split_line)
        rev_splits_content.append(split_line)
    rev_splits_content = "\n".join(rev_splits_content)

    rev_splits_path = path.nodes_to_file_path([split_path, "dna_splits.txt"])
    rev_splits_file = open(rev_splits_path, "w+")
    rev_splits_file.write(rev_splits_content)
    rev_splits_file.close()

    container.seq_dir.dna_splits_path = rev_splits_path

    # Amino
    amino_splits_content = []
    for split in amino_splits.items():
        color_string, weight = split
        colors = [int(color) for color in color_string.split(",")]
        split_line = [str(weight)]
        split_line.extend([amino_files[color] for color in colors])
        split_line = "\t".join(split_line)
        amino_splits_content.append(split_line)
    amino_splits_content = "\n".join(amino_splits_content)

    amino_splits_path = path.nodes_to_file_path([split_path, "amino_splits.txt"])
    amino_splits_file = open(amino_splits_path, 'w+')
    amino_splits_file.write(amino_splits_content)
    amino_splits_file.close()

    container.seq_dir.amino_splits_path = amino_splits_path

    log_print(f"STATS:")
    log_print(f"\tDNA_SEQ_LEN:\t{seq_len}")
    log_print(f"\tAMINO_SEQ_LEN:\t{amino_len}")
    log_print(f"\tMUTATIONS:\t{true_mutations}/{mutations}")
    log_print(f"\tREV_SPLITS:\t{len(rev_splits)}\n\tNREV_SPLITS:\t{len(nrev_splits)}")
    log_print(f"\tAMINO_SPLITS:\t{len(amino_splits)}")

def make_comp_data(container : Container):
    """
    Simple sequence generation including iupac letters
    """
    seq_path = container.seq_dir.fa_path
    seq_num = container.data.seq_num # Number of sequences
    seq_len = container.data.seq_len # Length of sequences

    seqs = []
    amino_seqs = []

    # If the container requires a base truth, the sequences are created accordingly
    log_print("[SEQ_MAKE] Building COMP test data")
    # Make the length a multiple of three for triplet usage
    seq_len = seq_len - (seq_len % 3)       # The length of the sequence as a multiple of three
    amino_len = int(seq_len / 3)            # The length of the amino acid sequence

    # --- [Sequence generation] ---
    # Create a master seq, rc_master seq and amino master seq
    master_seq = []
    rc_master_seq = []
    master_amino_seq = []
    # Create a master sequence
    for base_it in range(0, seq_len, 3):
        triplet = random.choice(TRIPLETS) # Choose a random triplet
        aminoacid = CODES_BY_TRIPLETS[triplet] # Translate
        master_seq.extend(triplet)
        master_amino_seq.append(aminoacid)
    # Compute the reverse complement master seq
    for base in master_seq:
        rc_master_seq.insert(0, COMPLEMENTS[base])
    
    seqs.append(master_seq)
    seqs.append(rc_master_seq)
    amino_seqs.append(master_amino_seq)
    amino_seqs.append(master_amino_seq)
    for i in range(seq_num - 2):
        seqs.append(deepcopy(master_seq))
        amino_seqs.append(deepcopy(master_amino_seq))
    log_print("[SEQ_MAKE] Created sequences")
    # --- [Sequence Mutation] ---
    # The distance between two split generating SNPs in the data set
    snp_distance = 3 * container.sans.k
    # The Master seq mutation rate
    m_rate = container.data.m_rate / 100
    # The number of split generating mutations
    mutations = min([int(seq_len * m_rate), int((seq_len - 2 * snp_distance) / snp_distance)])
    true_mutations = 0

    # --- Create mutations ---
    seq_ids = [i for i in range(2, seq_num)]
    for pos in range(snp_distance, seq_len - snp_distance, snp_distance):
        # If the target number of mutations is reached, stop mutation
        if true_mutations == mutations:
            break
        true_mutations += 1

        # The position in the amino sequence
        acid_pos = int((pos / 3) - 1)

        # The number of sequences to mutate
        split_size = random.randint(2, seq_num - 2)
        # The sequences to mutate
        target_ids = random.sample(seq_ids, split_size)

        # The possible bases
        bases = [element for element in BASES if element != master_seq[pos]]

        # Create a new random coding triplet and a mutated base
        new_triplet = ""
        m_base = ""
        while new_triplet not in CODES_BY_TRIPLETS and bases:
            m_base = random.choice(bases)
            new_triplet = "".join(master_seq[pos - 1:pos + 2])
            bases.remove(m_base)
        if new_triplet not in CODES_BY_TRIPLETS:
            continue

        # Apply mutation to all target sequences
        new_code = CODES_BY_TRIPLETS[new_triplet]
        for seq_id in target_ids:
            # Create DNA mutation
            seqs[seq_id][pos - 1] = m_base
            # Create AMINO mutation
            amino_seqs[seq_id][acid_pos] = new_code
    log_print(f"[SEQ_MAKE] ADDED MUTATIONS")
    # --- Create IUPAC noise ---
    pacs = 0
    for i in range(int(0.1 * mutations)):
        target_id = random.randint(3, seq_num - 1)
        pos = random.randint(0, seq_len -1)
        seqs[target_id][pos] = 'N'
        pacs += 1
    log_print(f"[SEQ_MAKE] ADDED IUPACS")
    # --- [Write Sequences] ---
    # DNA
    dna_files = [f"{i}.fa" for i in range(seq_num - 1)]
    dna_files.insert(1, f"0r.fa")
    dna_paths = [path.nodes_to_file_path([seq_path, dna_files[i]]) for i in range(seq_num)]
    dna_contents = [f">S{i}\n{''.join(seqs[i])}" for i in range(seq_num)]

    for i in range(seq_num):
        fasta_file = open(dna_paths[i], "w+")
        fasta_file.write(dna_contents[i])
        fasta_file.close()

    # Write lists
    rev_list_content = "\n".join(dna_files)
    rev_list_file = "list_dna.txt"
    rev_list_path = path.nodes_to_file_path([seq_path, rev_list_file])
    rev_list = open(rev_list_path, 'w+')
    rev_list.write(rev_list_content)

    container.seq_dir.dna_list_path = rev_list_path

    n_rev_list_content = "\n".join([file for file in dna_files if file != "0r.fa"])
    n_rev_list_file = "list_dna_nrev.txt"
    n_rev_list_path = path.nodes_to_file_path([seq_path, n_rev_list_file])
    n_rev_list = open(n_rev_list_path,  'w+')
    n_rev_list.write(n_rev_list_content)
    n_rev_list.close()

    container.seq_dir.dna_nrev_list_path = n_rev_list_path


    # AMINO
    amino_files = [f"AMINO_{i}.fa" for i in range(seq_num)]
    amino_paths = [path.nodes_to_file_path([seq_path, amino_files[i]]) for i in range(seq_num)]
    amino_contents = [f">A{i}\n{''.join(amino_seqs[i])}" for i in range(seq_num)]

    for i in range(seq_num):
        fasta_file = open(amino_paths[i], "w+")
        fasta_file.write(amino_contents[i])
        fasta_file.close()

    amino_list_content = "\n".join([file for file in amino_files])
    amino_list_file = "list_amino.txt"
    amino_list_path = path.nodes_to_file_path([seq_path, amino_list_file])
    amino_list = open(amino_list_path, "w+")
    amino_list.write(amino_list_content)
    amino_list.close()

    container.seq_dir.amino_list_path = amino_list_path

    log_print(f"STATS:")
    log_print(f"\tDNA_SEQ_LEN:\t{seq_len}")
    log_print(f"\tAMINO_SEQ_LEN:\t{amino_len}")
    log_print(f"\tMUTATIONS:\t{true_mutations}/{mutations}")
    log_print(f"\tPACS:\t{pacs}")




