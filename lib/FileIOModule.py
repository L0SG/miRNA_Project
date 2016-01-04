__author__ = 'L0SG'


def preprocess_smrna_seq(smrna_file):

    return


def create_ref_seq(ref_file):
    # Create reference sequence from fasta file and Return as List
    ref_name_list = []
    ref_seq_list = []
    tempseq = ""
    # Ignore the first line : must be starting with ">"
    line_split = ref_file.readline().strip().split()
    ref_name_list.append(line_split[0][1:])
    while 1:
        line = ref_file.readline().strip()
        if line == "":
            ref_seq_list.append(tempseq)
            break
        elif line.startswith(">"):
            line_split=line.split()
            ref_name_list.append(line_split[0][1:])
            ref_seq_list.append(tempseq)
            tempseq = ""
            continue
        else:
            tempseq += line
    ref_file.seek(0)
    return ref_name_list, ref_seq_list











