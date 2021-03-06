from lib import FileIOModule
from lib import SeqModule
from collections import Counter, OrderedDict
import subprocess
import re
import multiprocessing
import os
import operator
import argparse
import cPickle
import time
import sys

here = os.path.abspath(os.path.dirname(__file__))
subdir = "result"
path = os.path.join(here, subdir)
if not os.path.exists(path):
    os.mkdir(path)

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--reference", help="reference genome file location(.fa, .fasta)")
parser.add_argument("-i", "--input", help="input small RNA file location(.fa, .fasta)")
parser.add_argument("-o", "--output", help="output file location")
parser.add_argument("--bowtiepath", help="Bowtie location (default: command \"bowtie\")")
parser.add_argument("--RNAfoldpath", help="RNAfold location (default : command \"RNAfold\")")
parser.add_argument("-l", "--minlength", help="min. length of mature miRNA (default : 18)", type=int)
parser.add_argument("-L", "--maxlength", help="max. length of mature miRNA (default : 26)", type=int)
parser.add_argument("--multloci", help="max. multiple loci of miRNA matches to reference genome (default : 20)", type=int)
parser.add_argument("-d", "--distance", help="max. distance between miRNA-miRNA* of precursor hairpin loop (default : 35)", type=int)
parser.add_argument("-a", "--arm", help="length of arms(both ends) of precursor from mature miRNA(miRNA*) (default : 10)", type=int)
parser.add_argument("--mfe", help="min. abs. MFE for valid miRNA precursor (default : 18)", type=int)
parser.add_argument("-m", "--serialmismatch", help="max. serial mismatch of miRNA-miRNA* duplex (default : 2)", type=int)
parser.add_argument("-M", "--multmismatch", help="max. multiple mismatches of miRNA-miRNA* duplex (default : 2)", type=int)
parser.add_argument("-b", "--serialbulge", help="max. serial bulge of miRNA-miRNA* duplex (default : 2)", type=int)
parser.add_argument("-B", "--multbulge", help="max. multiple bulges of miRNA-miRNA* duplex (default : 2)", type=int)
parser.add_argument("-t", "--thread", help="number of CPU threads (default : 1)", type=int)
parser.add_argument("-s", "--step", help="step size of RNAfold precursor extend loop (default : 2)", type=int)
parser.add_argument("-c", "--mincount", help="min. read count of smRNA seq used for mature miRNA finding (default : 2)", type=int)
parser.add_argument("--batch_size", help="size of input data chunk for internal processing (default : number of CPU threads * 3)", type=int)
parser.add_argument("--plot", help="draw simple bar plot of length distribution of genome-aligned RNA sequence, need matplotlib to be installed (default : false)")
parser.add_argument("--annotate", help="automatically annotate miRNA cadidates, blastn required (default : false")
parser.add_argument("--mirbasepath", help="mirbase mature.fa location (default : mature.fa in current loccation)")
parser.add_argument("--blastnpath", help="blastn location (default : command \"blastn\")")
args = parser.parse_args()

# input file list
if args.reference:
    ref_file = open(args.reference, "r")
else:
    ref_file = open("ref.fa", "r")      # Reference genome file
if args.input:
    smrna_file = open(args.input, "r")
else:
    smrna_file = open("smrna.fa", "r")  # smRNA library file

# Bowtie path
if args.bowtiepath:
    bowtie_build_path = os.path.join(args.bowtiepath, 'bowtie-build')
    bowtie_path = os.path.join(args.bowtiepath, 'bowtie')
else:
    bowtie_build_path = "bowtie-build"
    bowtie_path = "bowtie"

# RNAfold path
if args.RNAfoldpath:
    RNAfold_path = os.path.join(args.RNAfoldpath, 'RNAfold')
else:
    RNAfold_path = "RNAfold"

# blastn path
if args.blastnpath:
    blastn_path = os.path.join(args.blastnpath, 'blastn')
else:
    blastn_path = "blastn"

# mirbase path
if args.mirbasepath:
    mirbase_path = args.mirbasepath
else:
    mirbase_path = os.path.join(os.getcwd(), "mature.fa")

# output file list
if args.output:
    path = args.output
    if not os.path.exists(path):
        os.mkdir(path)
if os.path.exists(os.path.join(path, "map")):
    output_map = open(os.path.join(path, "map"), "r")
    output_count_pos = open(os.path.join(path, "count_pos"), "rb")
    output_count_neg = open(os.path.join(path, "count_neg"), "rb")
    ref_count_dump_pos = cPickle.load(output_count_pos)
    ref_count_dump_neg = cPickle.load(output_count_neg)
output_precursor = open(os.path.join(path, "result_precursor.txt"), "w+")
output_precursor_collapsed = open(os.path.join(path, "result_precursor_collapsed.txt"), "w+")
output_mature = open(os.path.join(path, "result_mature.txt"), "w+")
output_distribution = open(os.path.join(path, "result_length_distribution.txt"), "w+")
output_tabular = open(os.path.join(path, "result_tabular_format.txt"), "w+")

# variable list
NUM_THREADS = multiprocessing.cpu_count()
MATURE_MIN_LEN = 18
MATURE_MAX_LEN = 26
MAX_SERIAL_MISMATCH = 2
MAX_MULT_MISMATCH = 2
MAX_SERIAL_BULGE = 2
MAX_MULT_BULGE = 2
MAX_MULTIPLE_LOCI = 20
DISTANCE_THRESHOLD = 35
ARM_EXTEND_THRESHOLD = 10
RNAFOLD_STEP = 2
MIN_ABS_MFE = 18
MIN_READ_COUNT_THRESHOLD = 2
DUPLICATE_FILTER_THRESHOLD = 10
DOMINANT_FACTOR = 0.9
NON_CANONICAL_PREC_FACTOR = 0.01   # smaller value makes it more "canonical"
DISCARD_NO_READ_PREC_FLAG = 1
PLOT_FLAG = 'false'
BATCH_SIZE = NUM_THREADS*3
ANNOTATE_FLAG = 'false'

if args.thread:
    NUM_THREADS = args.thread
if args.minlength:
    MATURE_MIN_LEN = args.minlength
if args.maxlength:
    MATURE_MAX_LEN = args.maxlength
if args.serialmismatch:
    MAX_SERIAL_MISMATCH = args.serialmismatch
if args.multmismatch:
    MAX_MULT_MISMATCH = args.multmismatch
if args.serialbulge:
    MAX_SERIAL_BULGE = args.serialbulge
if args.multbulge:
    MAX_MULT_BULGE = args.multbulge
if args.multloci:
    MAX_MULTIPLE_LOCI = args.multloci
if args.distance:
    DISTANCE_THRESHOLD = args.distance
if args.arm:
    ARM_EXTEND_THRESHOLD = args.arm
if args.step:
    RNAFOLD_STEP = args.step
if args.mfe:
    MIN_ABS_MFE = args.mfe
if args.mincount:
    MIN_READ_COUNT_THRESHOLD = args.mincount
if args.plot == 'true' or args.plot == 'True':
    PLOT_FLAG = args.plot
if args.annotate == 'true' or args.annotate == 'True':
    ANNOTATE_FLAG = args.annotate
if args.batch_size:
    BATCH_SIZE = args.batch_size
# DUPLICATE_FILTER_THRESHOLD, DOMINANT_FACTOR, NON_CANONICAL_PREC_FACTOR are internal variables
# internal variables are not expected to be changed by users, but possible if one wants to experiment

##################################### main script start #####################################

print("\nmiRNA Discovery Project")
print("Program by Lee Sang-Gil")
print("dept. of Applied Biology & Chemistry, Seoul National University")

# Read reference file and generate reference sequence list
print("Loading reference genome file to memory...")
ref_name_list, ref_seq_list = FileIOModule.create_ref_seq(ref_file)

# check whether index file was generated before
if os.path.exists(os.path.join(os.getcwd(), str(ref_file.name)+".1.ebwtl")):
    print("Index file detected, skipping index generation...")
else:
    # generate map file using bowtie

    # generate bowtie index
    # index files are saved to the same path where the reference file is stored
    if args.reference:
        ref_file_path = args.reference
    else:
        ref_file_path = os.path.join(os.getcwd(), "ref.fa")
    print("Generating index file with bowtie_build...")
    bowtie_build = subprocess.Popen([bowtie_build_path,
                               ref_file_path,
                               os.path.join(os.getcwd(), str(ref_file.name)),
                                     '--large-index'],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    bowtie_build.wait()

# map to reference genome and generate bowtie map file
# bowtie map file format is different, need to convert
if os.path.exists(os.path.join(path, "map")):
    print("map file detected, skipping map generation...")
else:
    start = time.time()
    output_map = open(os.path.join(path, "map"), "w+")
    output_count_pos = open(os.path.join(path, "count_pos"), "w+")
    output_count_neg = open(os.path.join(path, "count_neg"), "w+")
    if args.input:
        smrna_file_path = args.input
    else:
        smrna_file_path = os.path.join(os.getcwd(), "smrna.fa")
    print("Mapping smrna-seq to reference genome with bowtie...")
    bowtie = subprocess.Popen([bowtie_path, str(ref_file.name),
                               "-f", smrna_file_path,
                               os.path.join(path, "map_bowtie"),
                               "-v", "0", "-m", str(MAX_MULTIPLE_LOCI), "-a", "-t", "-p", str(NUM_THREADS),
                               '--large-index'],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    bowtie.wait()

    print("Converting bowtie map format to correct map format...")
    # open bowtie-generated map file (read-only, no need to be changed)
    output_bowtie = open(os.path.join(path, "map_bowtie"), "r")

    # convert bowtie map file format to correct form
    SeqModule.convert_bowtie_output(output_bowtie, output_map)
    output_map.seek(0, 0)

    print("Generating count data using map file...")
    # generate count data using map file
    ref_count_dump_pos, ref_count_dump_neg = SeqModule.count_generator(ref_name_list, output_map)
    output_map.seek(0, 0)

    # dump count data file for future usage and skip mapping
    cPickle.dump(ref_count_dump_pos, output_count_pos, -1)
    cPickle.dump(ref_count_dump_neg, output_count_neg, -1)
    output_count_pos.seek(0, 0)
    output_count_neg.seek(0, 0)
    end = time.time()
    print("Elapsed time for mapping : " + str(end - start) + " seconds")
    print("Mapping done")


# generate count list
print("Generating read count information list...")


class count_list(dict):
    def __missing__(self, key):
        return 0
ref_count_list_pos = []
ref_count_list_neg = []
for i in range(0, len(ref_seq_list)):
    ref_count_list_pos.append(count_list(ref_count_dump_pos[i]))
    ref_count_list_neg.append(count_list(ref_count_dump_neg[i]))


# use RNAfold to calculate MFE and select putative precursor
print("\n########## pre-miRNA discovery started ##########")
print("Calculating MFE of putative precursors with RNAfold...")


def precursor_generator(lines):
    output_precursor_infolist = []
    output_precursor_dblist = []
    reads_total_partial = 0
    length_distribution_partial = count_list({})

    for z in range(0, len(lines)):
        line_split = lines[z].split()
        # Rare occasion of improper line data : should skip it
        if len(line_split) != 7:
            continue
        # accumulate raw rna seq read counts for calculation of RPM
        reads_total_partial += int(line_split[1])

        # accumulate length distribution imformation
        # specify 5' end and add to the corresponding index
        seq_dist_check = line_split[6]
        five_prime = seq_dist_check[0]
        seq_length = len(seq_dist_check)
        dict_key = str(seq_length)+str(five_prime)
        length_distribution_partial[dict_key] += int(line_split[1])

        # Screen for Drosha / Dicer cutting sites (Inspired by miREAP)
        qualified_flag = 1
        name_list_index = ref_name_list.index(line_split[2])
        count = 0

        if line_split[5] == "+":
            count = ref_count_list_pos[name_list_index][int(line_split[3])]
            if count < 3:
                continue
            count_region = count
            count_sites = count
            for i in range(1, 20):
                if int(line_split[3])-i < 0 or int(line_split[3])+i >= len(ref_seq_list[name_list_index]):
                    continue
                if ref_count_list_pos[name_list_index][int(line_split[3])-i] > count \
                        or ref_count_list_pos[name_list_index][int(line_split[3])+i] > count:
                    qualified_flag = 0
                    break
                count_region += ref_count_list_pos[name_list_index][int(line_split[3])-i]
                count_region += ref_count_list_pos[name_list_index][int(line_split[3])+i]
                if i < 3:
                    count_sites += ref_count_list_pos[name_list_index][int(line_split[3])-i]
                    count_sites += ref_count_list_pos[name_list_index][int(line_split[3])+i]
            if float(count_sites)/count_region < DOMINANT_FACTOR or float(count)/count_sites < DOMINANT_FACTOR/2.0:
                qualified_flag = 0

        elif line_split[5] == "-":
            count = ref_count_list_neg[name_list_index][int(line_split[4])]
            if count < 3:
                continue
            count_region = count
            count_sites = count
            for i in range(1, 20):
                if int(line_split[4])-i < 0 or int(line_split[4])+i >= len(ref_seq_list[name_list_index]):
                    continue
                if ref_count_list_neg[name_list_index][int(line_split[4])-i] > count \
                        or ref_count_list_neg[name_list_index][int(line_split[4])+i] > count:
                    qualified_flag = 0
                    break
                count_region += ref_count_list_neg[name_list_index][int(line_split[4])-i]
                count_region += ref_count_list_neg[name_list_index][int(line_split[4])+i]
                if i < 3:
                    count_sites += ref_count_list_neg[name_list_index][int(line_split[4])-i]
                    count_sites += ref_count_list_neg[name_list_index][int(line_split[4])+i]
            if float(count_sites)/count_region < DOMINANT_FACTOR or float(count)/count_sites < DOMINANT_FACTOR/2.0:
                qualified_flag = 0

        if qualified_flag == 0:
            continue

        # Precursor Candidate Information Variable List
        pc_seq = ""
        pc_structure = ""
        pc_start = 0
        pc_end = 0
        pc_abs_energy = 0
        pc_norm_abs_energy = 0

        # Find min. MFE of fold structure and save it to result_precursor
        # WARNING : if abs. of calculated free energy is less than 10, output[2] does not contain proper value
        # Skipping this precursor line is proper, since threshold value is at least 18
        for k in range(0, len(ref_seq_list)):  # reference sequence list loop
            # 160205 : No need to search other genomes
            # but fixed, need better implementation (remove ref seq loop)
            if name_list_index != k:
                continue
            # 150907 : No need to loop arm extension? miREAP only uses const FLANK var (10)
            # disabling arm extension loop has no significant difference, but can reduce time complexity
            for i in range(ARM_EXTEND_THRESHOLD, ARM_EXTEND_THRESHOLD+1):  # arm extension loop (disabled)
                for j in range(int(line_split[4])-int(line_split[3]), DISTANCE_THRESHOLD+int(line_split[4])-int(line_split[3]), RNAFOLD_STEP):  # distance loop

                    # Assuming -5p mature sequence
                    start = int(line_split[3])-i
                    end = int(line_split[4])+j+i
                    if start >= 0 and end < len(ref_seq_list[k]):  # continue only if both indices are valid
                        if line_split[5] == "+":
                            rna_fold_seq = ref_seq_list[k][start:end]
                        elif line_split[5] == "-":
                            rna_fold_seq = SeqModule.create_star(ref_seq_list[k][start:end])
                        if "N" in rna_fold_seq:
                            continue
                        rnafold = subprocess.Popen([RNAfold_path, "--noconv", "-d2", "--noPS", "--noLP"],
                                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                        output = rnafold.communicate(rna_fold_seq)[0].split()
                        # Discard non-canonical (i.e. "hard to identify") precursor
                        pc_structure_left = output[1].strip("\n")[0:len(output[1].strip("\n"))/2]
                        pc_structure_right = output[1].strip("\n")[len(output[1].strip("\n"))/2:len(output[1].strip("\n"))]
                        num_open_left = pc_structure_left.count("(")
                        num_close_left = pc_structure_left.count(")")
                        num_open_right = pc_structure_right.count("(")
                        num_close_right = pc_structure_right.count(")")
                        if num_open_left == 0 or num_close_right == 0:
                            continue
                        if float(num_close_left)/num_open_left > NON_CANONICAL_PREC_FACTOR or\
                                                float(num_open_right)/num_close_right > NON_CANONICAL_PREC_FACTOR:
                            continue
                        abs_energy = re.findall(r'\d*\.\d*', str(output[2]))
                        if abs_energy != []:
                            if float(abs_energy[0]) >= MIN_ABS_MFE:
                                norm_abs_energy = float(abs_energy[0])/len(rna_fold_seq)
                                if pc_seq == []: # bad implementation, need to repair
                                    pc_seq = output[0].strip()
                                    pc_structure = output[1].strip("\n")
                                    pc_start = start
                                    pc_end = end
                                    pc_abs_energy = float(abs_energy[0])
                                    pc_norm_abs_energy = norm_abs_energy
                                elif norm_abs_energy > pc_norm_abs_energy:
                                    pc_seq = output[0].strip()
                                    pc_structure = output[1].strip("\n")
                                    pc_start = start
                                    pc_end = end
                                    pc_abs_energy = float(abs_energy[0])
                                    pc_norm_abs_energy = norm_abs_energy

                    # Assuming -3p mature sequence
                    start = int(line_split[3])-j-i
                    end = int(line_split[4])+i
                    if start >= 0 and end < len(ref_seq_list[k]):   # continue only if both indices are valid
                        if line_split[5] == "+":
                            rna_fold_seq = ref_seq_list[k][start:end]
                        elif line_split[5] == "-":
                            rna_fold_seq = SeqModule.create_star(ref_seq_list[k][start:end])
                        if "N" in rna_fold_seq:
                            continue
                        rnafold = subprocess.Popen([RNAfold_path, "--noconv", "-d2", "--noPS", "--noLP"],
                                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                        output = rnafold.communicate(rna_fold_seq)[0].split()
                        # Discard non-canonical (i.e. "hard to identify") precursor
                        pc_structure_left = output[1].strip("\n")[0:len(output[1].strip("\n"))/2]
                        pc_structure_right = output[1].strip("\n")[len(output[1].strip("\n"))/2:len(output[1].strip("\n"))]
                        num_open_left = pc_structure_left.count("(")
                        num_close_left = pc_structure_left.count(")")
                        num_open_right = pc_structure_right.count("(")
                        num_close_right = pc_structure_right.count(")")
                        if num_open_left == 0 or num_close_right == 0:
                            continue
                        if float(num_close_left)/num_open_left > NON_CANONICAL_PREC_FACTOR or\
                                                float(num_open_right)/num_close_right > NON_CANONICAL_PREC_FACTOR:
                            continue
                        abs_energy = re.findall(r'\d*\.\d*', str(output[2]))
                        if abs_energy != []:
                            if float(abs_energy[0]) >= MIN_ABS_MFE:
                                norm_abs_energy = float(abs_energy[0])/len(rna_fold_seq)
                                if pc_seq == []: # bad implementation, need to repair
                                    pc_seq = output[0].strip()
                                    pc_structure = output[1].strip("\n")
                                    pc_start = start
                                    pc_end = end
                                    pc_abs_energy = float(abs_energy[0])
                                    pc_norm_abs_energy = norm_abs_energy
                                elif norm_abs_energy > pc_norm_abs_energy:
                                    pc_seq = output[0].strip()
                                    pc_structure = output[1].strip("\n")
                                    pc_start = start
                                    pc_end = end
                                    pc_abs_energy = float(abs_energy[0])
                                    pc_norm_abs_energy = norm_abs_energy

        if pc_seq != "":
            output_precursor_infolist.append(lines[z].strip()+"\t"+str(pc_abs_energy)+"\t"+str(pc_norm_abs_energy)+"\t"+
                                             str(pc_start)+"\t"+str(pc_end)+"\n")
            output_precursor_dblist.append(pc_seq+"\n"+pc_structure+"\n")
            continue
    # create counters (subclass of dict) to "merge" partial length distribution dicts later
    length_distribution_counter = Counter(length_distribution_partial)
    return output_precursor_infolist, output_precursor_dblist, reads_total_partial, length_distribution_counter


# wrapper function for progress monitoring
def precursor_generator_wrapper(args_):
    input_original, q = args_
    result = precursor_generator(input_original)
    q.put(0)
    return result

# precursor output file header
output_precursor.write("Name\tRead_Count\tChr_Name\tMature_Start\tMature_End\tPos\tMFE\tNorm_MFE\tPrec_Start\tPrec_End\n")

# precursor generator multiprocessing procedure
if __name__ == '__main__':
    start = time.time()
    lines = output_map.readlines()
    manager = multiprocessing.Manager()
    queue = manager.Queue()

    # split original map data to multiple lists for efficient multiprocessing
    numlines = 500
    num_chunk = len(lines)/numlines
    arguments = [(lines[line:line+numlines], queue) for line in range(0, len(lines), numlines)]

    # transform arguments to batches for memory usage suppression
    # optimal batch size is probably dependent on the system
    batch_size = BATCH_SIZE
    arguments = [arguments[i:i+batch_size] for i in range(0, len(arguments), batch_size)]

    # apply multiprocessing, one batch at a time
    # closing pool and extending partial result to original result suppresses memory usage
    result_list = []
    for arguments_partial in arguments:
        pool = multiprocessing.Pool(processes=NUM_THREADS)
        result_list_partial = pool.map_async(precursor_generator_wrapper, arguments_partial)
        while True:
            if result_list_partial.ready():
                break
            size = queue.qsize()
            sys.stdout.write('\r%% of map data processed : %.2f %%' % (float(size)/float(num_chunk)*100))
            time.sleep(0.1)
        result_list.extend(result_list_partial.get())
        pool.close()
        pool.join()
    sys.stdout.write('\r%% of map data processed : %.2f %%' % 100)
    print (' done')

    # merging procedure
    output_precursor_info_result = []
    output_precursor_db_result = []
    reads_total = 0
    length_distribution = Counter({})
    for i in range(0, len(result_list)):
        output_precursor_info_result += result_list[i][0]
        output_precursor_db_result += result_list[i][1]
        reads_total += result_list[i][2]
        length_distribution += result_list[i][3]
    for i in range(0, len(output_precursor_info_result)): # len of info_result, db_result must be same
        output_precursor.write(output_precursor_info_result[i])
        output_precursor.write(output_precursor_db_result[i])
output_precursor.seek(0, 0)
end = time.time()
print("Elapsed time for RNAfold : "+str(end-start)+" seconds")
print("Precursor generation done")

# collapse duplicate precursors
# if either start index or end index of the two precursors are "slightly" different, they're assumed to the duplicate
# same precursor sequences with different genome locations (i.e. different indexes) are NOT assumed to the duplicate
print("Collapsing duplicate precursors...")
output_precursor_collapsed.write("Name\tRead_Count\tChr_Name\tMature_Start\tMature_End\tPos\tMFE\tNorm_MFE\tPrec_Start\tPrec_End\n")
output_precursor.readline()
prec_info_list = []

while 1:
    line = output_precursor.readline()
    if line == "":
        break
    line_split = line.split()
    duplicate_flag = 0
    for i in range(0, len(prec_info_list)):
        if str(line_split[2]) == prec_info_list[i][0]:
            if abs(int(line_split[9])-int(prec_info_list[i][1])) <= DUPLICATE_FILTER_THRESHOLD \
                    or abs(int(line_split[10])-int(prec_info_list[i][2])) <= DUPLICATE_FILTER_THRESHOLD:
                duplicate_flag = 1
                break
    if duplicate_flag == 1:
        output_precursor.readline()
        output_precursor.readline()
        continue
    prec_info_list.append([line_split[2], line_split[9], line_split[10]])
    output_precursor_collapsed.write(line)
    output_precursor_collapsed.write(output_precursor.readline())
    output_precursor_collapsed.write(output_precursor.readline())
output_precursor_collapsed.seek(0, 0)
print("Collapsing done")
print("########## pre-miRNA discovery complete ##########")


print("\n########## mature miRNA-miRNA* duplex calculation started ##########")
# reads_total is calculated before, but for debugging purpose calculate again at this point
reads_total = 0
smrna_file.seek(0, 0)
for i, l in enumerate(smrna_file):
    if l.startswith('>'):
        reads_total += int(l.split()[1])
    pass
smrna_file.seek(0, 0)
if not args.mincount:
    print("min. readcount threshold parameter (-c --mincount) not set, starting automatic read suppression...")
    print("checking smrna file size...")
    smrna_file_length = (i+1)/2
    smrna_file.seek(0, 0)
    print("Total number of unique smrna reads : " + str(smrna_file_length))
    MIN_READ_COUNT_THRESHOLD = smrna_file_length/500000 + 1
    print("min. readcount threshold automatically set to " + str(MIN_READ_COUNT_THRESHOLD))
print("Mapped reads with count lower than " + str(MIN_READ_COUNT_THRESHOLD) +
      " will not be used for further processing")

# load map file to make aligned output
# 2d list, [i] [0]:name, [1]:readcount, [2]:chr, [3]:mature_start, [4]:mature_end, [5]:pos, [6]:seq
# [7]: MFE, [8]:norm_MFE, [9]:prec_start, [10]:prec_end
print("Loading map data for alignment...")
map_data = []
output_map.seek(0, 0)

while 1:
    line_split = output_map.readline().strip().split()
    if line_split == []:
        break
    # if read count is below minimum read count threshold, discard it
    # otherwise there are too many sites to search
    if int(line_split[1]) >= MIN_READ_COUNT_THRESHOLD:
        map_data.append(line_split)
map_data = sorted(map_data, key=operator.itemgetter(3))
print("Loading done")


# Select Precursors which have valid star seq
print("Finding valid star sequence for each putative precursor...")
result_count = 0


def mature_generator(lines):
    global map_data
    # each loop should read exactly 3 lines
    output_list=[]
    iterator = 0
    while 1:
        if iterator == len(lines):
            break
        line_info = lines[iterator].strip()
        if line_info == "":
            break
        line_seq = lines[iterator+1].strip()
        line_db = lines[iterator+2].strip()
        iterator += 3

        # if no read data is matched in putative precursors, discard it
        if DISCARD_NO_READ_PREC_FLAG:
            no_read_prec_flag = SeqModule.check_no_read_prec(line_info, map_data, MIN_READ_COUNT_THRESHOLD)
            if no_read_prec_flag is True:
                continue

        # check conserved sequence with blastn
        # if this line_info is classified as conserved sequence, update line_info
        # no need to find duplex, just mark 5p and 3p index corresponding to matched information
        updated_flag = False
        if ANNOTATE_FLAG == 'true' or ANNOTATE_FLAG == 'True':
            line_info, updated_flag = SeqModule.check_conserved_seq(line_info, line_seq,
                                                                    blastn_path, mirbase_path, ARM_EXTEND_THRESHOLD)
        # if updated_flag is True:
            # start_5p, end_5p, start_3p, end_3p = SeqModule.find_location(line_info, line_seq, line_db)
        # else, do the code below
        ###########################################################

        # Discard non-canonical (i.e. "hard to identify") precursor
        # "Asymmetric" dot-bracket notation precursor : low accuracy, hard to identify star seq, and too many outputs
        # if ")" portion is large in "left side", it's non-canonical
        line_db_left = line_db[0:len(line_db)/2]
        num_open = line_db_left.count("(")
        num_close = line_db_left.count(")")
        if float(num_close)/num_open > NON_CANONICAL_PREC_FACTOR:
            continue

        # find valid star sequence from putative precursors
        start_5p, end_5p, start_3p, end_3p = SeqModule.star_identifier_v2(line_db, MATURE_MIN_LEN, MATURE_MAX_LEN,
                                                                          MAX_SERIAL_MISMATCH, MAX_MULT_MISMATCH,
                                                                          MAX_SERIAL_BULGE, MAX_MULT_BULGE)
        if start_5p == 0 and end_5p == 0 and start_3p == 0 and end_3p == 0:  # star seq not found
            continue

        # write putative precursor to the output file
        output_form = SeqModule.generate_output_form(line_info, line_seq, line_db,
                                                     start_5p, start_3p, end_5p, end_3p,
                                                     map_data, MIN_READ_COUNT_THRESHOLD)
        output_list.append(output_form)
    return output_list


def mature_generator_v2(lines):
    global map_data
    # each loop should read exactly 3 lines
    output_list=[]
    iterator = 0
    while 1:
        if iterator == len(lines):
            break
        line_info = lines[iterator].strip()
        if line_info == "":
            break
        line_seq = lines[iterator+1].strip()
        line_db = lines[iterator+2].strip()
        iterator += 3

        # if no read data is matched in putative precursors, discard it
        if DISCARD_NO_READ_PREC_FLAG:
            no_read_prec_flag = SeqModule.check_no_read_prec(line_info, line_seq, map_data, MIN_READ_COUNT_THRESHOLD)
            if no_read_prec_flag is True:
                continue
        # get alignment form first
        output_form = SeqModule.generate_alignment_form(line_info, line_seq, line_db,
                                                     map_data, MIN_READ_COUNT_THRESHOLD)

        # if mature_prime is bad, discard it
        if output_form[0].split()[12] == 'bad':
            continue

        # check conserved sequence with blastn
        # if this line_info is classified as conserved sequence, update line_info
        # no need to find duplex, just mark 5p and 3p index corresponding to matched information
        updated_flag = False
        if ANNOTATE_FLAG == 'true' or ANNOTATE_FLAG == 'True':
            output_form[0], updated_flag = SeqModule.check_conserved_seq(output_form[0], line_seq,
                                                                    blastn_path, mirbase_path, ARM_EXTEND_THRESHOLD)
            if updated_flag is True:
                start_5p, end_5p, start_3p, end_3p = SeqModule.star_identifier_v2_conserved(output_form[0], line_seq, line_db,
                                                                                            MATURE_MIN_LEN, MATURE_MAX_LEN,
                                                                                            MAX_SERIAL_MISMATCH, MAX_MULT_MISMATCH,
                                                                                            MAX_SERIAL_BULGE, MAX_MULT_BULGE)
            else:
                start_5p, end_5p, start_3p, end_3p = SeqModule.star_identifier_v2(line_db, MATURE_MIN_LEN, MATURE_MAX_LEN,
                                                                                  MAX_SERIAL_MISMATCH, MAX_MULT_MISMATCH,
                                                                                  MAX_SERIAL_BULGE, MAX_MULT_BULGE)
            if start_5p == 0 and end_5p == 0 and start_3p == 0 and end_3p == 0:  # star seq not found
                continue

            # write putative precursor to the output file
            output_form = SeqModule.generate_output_form(output_form[0], line_seq, line_db,
                                                         start_5p, start_3p, end_5p, end_3p,
                                                         map_data, MIN_READ_COUNT_THRESHOLD)
        # else, do the code below
        ###########################################################
        elif ANNOTATE_FLAG == 'false' or ANNOTATE_FLAG == 'False':
            # Discard non-canonical (i.e. "hard to identify") precursor
            # "Asymmetric" dot-bracket notation precursor : low accuracy, hard to identify star seq, and too many outputs
            # if ")" portion is large in "left side", it's non-canonical

            line_db_left = line_db[0:len(line_db)/2]
            num_open = line_db_left.count("(")
            num_close = line_db_left.count(")")
            if float(num_close)/num_open > NON_CANONICAL_PREC_FACTOR:
                continue

            # find valid star sequence from putative precursors
            start_5p, end_5p, start_3p, end_3p = SeqModule.star_identifier_v2(line_db, MATURE_MIN_LEN, MATURE_MAX_LEN,
                                                                              MAX_SERIAL_MISMATCH, MAX_MULT_MISMATCH,
                                                                              MAX_SERIAL_BULGE, MAX_MULT_BULGE)
            if start_5p == 0 and end_5p == 0 and start_3p == 0 and end_3p == 0:  # star seq not found
                continue

            # write putative precursor to the output file
            output_form = SeqModule.generate_output_form(output_form[0], line_seq, line_db,
                                                         start_5p, start_3p, end_5p, end_3p,
                                                         map_data, MIN_READ_COUNT_THRESHOLD)

        output_list.append(output_form)
    return output_list


# wrapper function for progress monitoring
def mature_generator_wrapper(args_):
    input_original, q = args_
    result = mature_generator_v2(input_original)
    q.put(0)
    return result

# mature generator multiprocessing procedure
if __name__ == '__main__':
    start = time.time()
    lines = output_precursor_collapsed.readlines()
    # discard header line
    lines = lines[1:]
    manager = multiprocessing.Manager()
    queue = manager.Queue()
    # numlines MUST be 3
    numlines = 3
    num_chunk = len(lines)/numlines
    arguments = [(lines[line:line+numlines], queue) for line in range(0, len(lines), numlines)]

    # transform arguments to batches for memory usage suppression
    # optimal batch size is probably dependent on the system
    batch_size = BATCH_SIZE*300
    arguments = [arguments[i:i + batch_size] for i in range(0, len(arguments), batch_size)]

    # apply multiprocessing, one batch at a time
    # closing pool and extending partial result to original result suppresses memory usage
    output_list = []
    for arguments_partial in arguments:
        pool = multiprocessing.Pool(processes=NUM_THREADS)
        output_list_partial = pool.map_async(mature_generator_wrapper, arguments_partial)
        while True:
            if output_list_partial.ready():
                break
            size = queue.qsize()
            sys.stdout.write('\r%% of precursor data processed : %.2f %%' % (float(size) / float(num_chunk) * 100))
            time.sleep(0.1)
        output_list.extend(filter(None, output_list_partial.get()))
        pool.close()
        pool.join()
    sys.stdout.write('\r%% of precursor data processed : %.2f %%' % 100)
    print (' done')

    ##############
    # post-process output list using the criteria used for the paper
    output_list_new = []
    for idx in xrange(0, len(output_list)):
        line_info = output_list[idx][0][0]
        line_seq = output_list[idx][0][1]
        line_db = output_list[idx][0][2]
        if len(line_info.split()[6]) < MATURE_MIN_LEN or len(line_info.split()[6]) > MATURE_MAX_LEN:
            continue
        # ad-hoc hack for high-confidence annotation
        if SeqModule.structure_check(line_info, line_seq, line_db, MAX_SERIAL_BULGE, MAX_SERIAL_MISMATCH) is False:
            continue
        output_list_new.append(output_list[idx])
    output_list = output_list_new
    ##############

    # sort output_list to display conserved miRNAs first
    output_list.sort(key=operator.itemgetter(0))
    # result_mature.txt generation

    output_tabular.write(str(("Name\tRead_Count\tRPM\tChr_Name\tMature_Start\tMature_End\tPos\tSeq\tMFE\tNorm_MFE\tPrec_Start\tPrec_End\tMature_len\tMature_prime\n")))
    for i in output_list:
        output_mature.write(str(("Name\tRead_Count\tRPM\tChr_Name\tMature_Start\tMature_End\tPos\tSeq\tMFE\tNorm_MFE\tPrec_Start\tPrec_End\tMature_len\tMature_prime\n")))
        for j in i:
            info = j[0].split()
            RPM = str(format(float(info[1]) / float(reads_total) * 1000000, '.2f'))
            info_new = info[0:2]+[RPM]+info[2:]
            info_new = '\t'.join(info_new)+'\n'
            output_mature.write(info_new)
            output_tabular.write(info_new)
            for k in xrange(1, len(j)):
                output_mature.write(j[k])
        output_mature.write("\n")
    """
    # result_tabular.txt generation
    output_tabular.write("ID\treads\tRead_Count\tRPM\Chromosome\tStart\tEnd\tStrand\tSeq\n")
    for i in output_list:
        info = i[0][0].split()
        ID = str(info[0])
        reads = str(info[1])
        RPM = str(float(info[1])/float(reads_total)*1000000)
        chromosome = str(info[2])
        seq_start = str(info[3])
        seq_end = str(info[4])
        strand = str(info[5])
        seq = str(info[6])
        output_tabular.write(ID+'\t'+reads+'\t'+RPM+'\t'+chromosome+
                             '\t'+seq_start+'\t'+seq_end+'\t'+strand+'\t'+seq+'\n')
    """
    # result_distribution.txt generation
    length_list = length_distribution.items()
    length_list.sort(key=lambda tup: tup[0])
    length_key = [i[0] for i in length_list]
    length_value = [int(i[1]) for i in length_list]
    length_value_RPM = [float(i)/float(reads_total)*1000000 for i in length_value]
    for key in length_key:
        output_distribution.write(str(key)+'\t')
    output_distribution.write('\n')
    for value in length_value:
        output_distribution.write(str(value)+'\t')
    output_distribution.write('\n')
    for value in length_value_RPM:
        output_distribution.write(str(value)+'\t')
    output_distribution.write('\n')
    if PLOT_FLAG == 'true' or PLOT_FLAG == 'True':
        from matplotlib import pyplot
        # simple bar plot generation
        pyplot.bar(xrange(0, len(length_key)), length_value)
        pyplot.xticks(xrange(0, len(length_key)), length_key, rotation='vertical', size='small', ha='left')
        pyplot.autoscale()
        pyplot.savefig(os.path.join(path, "length_distribution.png"), format='png')

end = time.time()
print("Elapsed time for calculating mature miRNA info : "+str(end-start)+" seconds")
print("########## mature miRNA-miRNA* duplex calculation complete ##########")
print("\nDone : "+str(len(output_list))+" miRNA found")
print("Results are generated at : "+str(path))
print("See result_mature.txt for detailed alignment information")
print("See result_tabular_format.txt for simplified tabular format of miRNA candidates")
print("See result_length_distribution.txt for length distribution of genome-mapped RNA-seq reads")


##################################### main script end #####################################
