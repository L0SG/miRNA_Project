from lib import FileIOModule
from lib import SeqModule
import subprocess
import re
import json
import multiprocessing
import os
import operator
import argparse

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
parser.add_argument("-p", "--RNAfoldpath", help="RNAfold location (default : command \"RNAfold\")")
parser.add_argument("-t", "--thread", help="number of CPU threads (default : 1)", type=int)
parser.add_argument("-l", "--minlength", help="min. length of mature miRNA (default : 18)", type=int)
parser.add_argument("-L", "--maxlength", help="max. length of mature miRNA (default : 26)", type=int)
parser.add_argument("-m", "--serialmismatch", help="max. serial mismatch of miRNA-miRNA* duplex (default : 2)", type=int)
parser.add_argument("-M", "--multmismatch", help="max. multiple mismatches of miRNA-miRNA* duplex (default : 2)", type=int)
parser.add_argument("-b", "--serialbulge", help="max. serial bulge of miRNA-miRNA* duplex (default : 2)", type=int)
parser.add_argument("-B", "--multbulge", help="max. multiple bulges of miRNA-miRNA* duplex (default : 2)", type=int)
parser.add_argument("-u", "--multloci", help="max. multiple loci of miRNA matches to reference genome (default : 20)", type=int)
parser.add_argument("-d", "--distance", help="max. distance between miRNA-miRNA* of precursor hairpin loop (default : 35)", type=int)
parser.add_argument("-a", "--arm", help="length of arms(both ends) of precursor from mature miRNA(miRNA*) (default : 10)", type=int)
parser.add_argument("-s", "--step", help="step size of RNAfold precursor extend loop (default : 2)", type=int)
parser.add_argument("-f", "--mfe", help="min. abs. MFE for valid miRNA precursor (default : 18)", type=int)
parser.add_argument("-c", "--mincount", help="min. read count of smRNA seq displayed at results (default : 2)", type=int)
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

# RNAfold path
if args.RNAfoldpath:
    RNAfold_path = args.RNAfoldpath
else:
    RNAfold_path = "RNAfold"

# output file list
if args.output:
    path = args.output
if os.path.exists(os.path.join(path, "map")):
    output_map = open(os.path.join(path, "map"), "r")
    output_count_pos = open(os.path.join(path, "count_pos"), "r")
    output_count_neg = open(os.path.join(path, "count_neg"), "r")
output_precursor = open(os.path.join(path, "result_precursor.txt"), "w+")
output_precursor_collapsed = open(os.path.join(path, "result_precursor_collapsed.txt"), "w+")
output_mature = open(os.path.join(path, "result_mature.txt"), "w+")

# variable list
NUM_THREADS = 12
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
# DUPLICATE_FILTER_THRESHOLD, DOMINANT_FACTOR, NON_CANONICAL_PREC_FACTOR are internal variables

##################################### main script start #####################################

print("miRNA Discovery Project")
print("Program by Lee Sang-Gil")
print("dept. of Applied Biology & Chemistry, Seoul National University")

# Read reference file and generate reference sequence list
ref_name_list, ref_seq_list = FileIOModule.create_ref_seq(ref_file)
ref_count_list_pos = []
ref_count_list_neg = []
for i in range(0, len(ref_seq_list)):
    ref_count_list_pos.append([0]*len(ref_seq_list[i]))
    ref_count_list_neg.append([0]*len(ref_seq_list[i]))


# check whether map file was generated before
if os.path.exists(os.path.join(path, "map")):
    print("map file detected, skipping map generation...")
    # Load previously Dumped count data
    ref_count_list_pos = json.load(output_count_pos)
    ref_count_list_neg = json.load(output_count_neg)

else:
    # read smrna file and generate map file
    output_map = open(os.path.join(path, "map"), "w+")
    output_count_pos = open(os.path.join(path, "count_pos"), "w+")
    output_count_neg = open(os.path.join(path, "count_neg"), "w+")
    print("Generating map file from seq library...")


    def map_generator(lines):
        output_map_linelist = []
        ref_count_list_poslist = []
        ref_count_list_neglist = []

        for i in range(0, len(ref_seq_list)):
            ref_count_list_poslist.append([])
            ref_count_list_neglist.append([])

        for z in range(0, len(lines), 2):
            smrna_info = lines[z].strip().split()
            smrna_seq = lines[z+1].strip()  # seq for pos match
            smrna_star = SeqModule.create_star(smrna_seq)  # seq for neg match
            # decide whether smrna_seq perfectly match genome seq
            # return type : 2D list
            # each matched case returns zero based begin, end index (2-element tuple)
            # character of end index is NOT contained in seq
            # [0][0] : first genome, first matched indices, [0][1] : first genome, second matched indices
            # [1][0] : second genome, first matched indices, [1][1] : second matched indices ...
            pm_index_list_pos = SeqModule.find_perfect_match(smrna_seq, ref_seq_list,
                                                             MATURE_MIN_LEN, MATURE_MAX_LEN, MAX_MULTIPLE_LOCI)
            pm_index_list_neg = SeqModule.find_perfect_match(smrna_star, ref_seq_list,
                                                             MATURE_MIN_LEN, MATURE_MAX_LEN, MAX_MULTIPLE_LOCI)
            for i in range(0, len(pm_index_list_pos)):
                for j in range(0, len(pm_index_list_pos[i])):
                    output_map_linelist.append(smrna_info[0]+"\t"+smrna_info[1]+"\t"+ref_name_list[i]+"\t"+
                                               str(pm_index_list_pos[i][j][0])+"\t"+str(pm_index_list_pos[i][j][1])+
                                               "\t"+"+"+"\t"+smrna_seq+"\n")
                    ref_count_list_poslist[i].append([pm_index_list_pos[i][j][0], int(smrna_info[1])])
            for i in range(0, len(pm_index_list_neg)):
                for j in range(0, len(pm_index_list_neg[i])):
                    output_map_linelist.append(smrna_info[0]+"\t"+smrna_info[1]+"\t"+ref_name_list[i]+"\t"+
                                               str(pm_index_list_neg[i][j][0])+"\t"+str(pm_index_list_neg[i][j][1])+
                                               "\t"+"-"+"\t"+smrna_seq+"\n")
                    ref_count_list_neglist[i].append([pm_index_list_neg[i][j][1], int(smrna_info[1])])
            # count_list_list[i][j][0] : i-th genome, j-th match begin index
            # count_list_list[i][j][1] : i-th genome, j-th match count
        return output_map_linelist, ref_count_list_poslist, ref_count_list_neglist


    if __name__ == '__main__':
        lines = smrna_file.readlines()
        pool = multiprocessing.Pool(processes=NUM_THREADS)
        numlines = 1000
        result_list = pool.map(map_generator, (lines[line:line+numlines] for line in xrange(0, len(lines), numlines)))
        output_map_result=[]
        ref_count_list_pos_result=[]
        ref_count_list_neg_result=[]

        for i in range(0, len(ref_seq_list)):
            ref_count_list_pos_result.append([])
            ref_count_list_neg_result.append([])

        for i in range(0, len(result_list)):
            output_map_result += result_list[i][0]
            for j in range(0, len(ref_seq_list)):
                ref_count_list_pos_result[j].extend(result_list[i][1][j])
                ref_count_list_neg_result[j].extend(result_list[i][2][j])

        for i in range(0, len(output_map_result)):
            output_map.write(str(output_map_result[i]))

        for i in range(0, len(ref_seq_list)):
            for j in range(0, len(ref_count_list_pos_result[i])):
                ref_count_list_pos[i][int(ref_count_list_pos_result[i][j][0])] += ref_count_list_pos_result[i][j][1]
            for j in range(0, len(ref_count_list_neg_result[i])):
                ref_count_list_neg[i][int(ref_count_list_neg_result[i][j][0])] += ref_count_list_neg_result[i][j][1]
    output_map.seek(0, 0)

    # Dump count data file for future usage and skip mapping
    json.dump(ref_count_list_pos, output_count_pos)
    json.dump(ref_count_list_neg, output_count_neg)
    output_count_pos.seek(0, 0)
    output_count_neg.seek(0, 0)
    print("Generating Done")

# use RNAfold to calculate MFE and select putative precursor
print("Calculating MFE of putative precursors with RNAfold...")


def precursor_generator(lines):
    output_precursor_infolist = []
    output_precursor_dblist = []

    for z in range(0, len(lines)):
        line_split = lines[z].split()
        # Rare occasion of improper line data : should skip it
        if len(line_split) != 7:
            continue
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
                if int(line_split[3])-i < 0 or int(line_split[3])+i >= len(ref_count_list_pos[name_list_index]):
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
                if int(line_split[4])-i < 0 or int(line_split[4])+i >= len(ref_count_list_neg[name_list_index]):
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
    return output_precursor_infolist, output_precursor_dblist


output_precursor.write("Name\tRead_Count\tChr_Name\tMature_Start\tMature_End\tPos\tMFE\tNorm_MFE\tPrec_Start\tPrec_End\n")
if __name__ == '__main__':
    lines = output_map.readlines()
    pool = multiprocessing.Pool(processes=NUM_THREADS)
    numlines = 50
    result_list = pool.map(precursor_generator, (lines[line:line+numlines] for line in range(0, len(lines), numlines)))
    output_precursor_info_result = []
    output_precursor_db_result = []
    for i in range(0, len(result_list)):
        output_precursor_info_result += result_list[i][0]
        output_precursor_db_result += result_list[i][1]
    for i in range(0, len(output_precursor_info_result)): # len of info_result, db_result must be same
        output_precursor.write(output_precursor_info_result[i])
        output_precursor.write(output_precursor_db_result[i])
output_precursor.seek(0, 0)
print("Calculating Done")


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
    map_data.append(line_split)
map_data = sorted(map_data, key=operator.itemgetter(3))
print("Loading Done")


# Select Precursors which have valid star seq
print("Finding valid star sequence for each putative precursor...")
result_count = 0
output_precursor_collapsed.readline()  # discard header line

while 1:
    line_info = output_precursor_collapsed.readline().strip()
    if line_info == "":
        break
    line_seq = output_precursor_collapsed.readline().strip()
    line_db = output_precursor_collapsed.readline().strip()

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

    # if no read data is matched in putative precursors, discard it
    if DISCARD_NO_READ_PREC_FLAG:
        count = 0
        if line_info.split()[5] == "+":
            for i in range(0, len(map_data)):
                if map_data[i][2] == line_info.split()[2] and int(map_data[i][3]) >= int(line_info.split()[9])\
                        and int(map_data[i][4]) <= int(line_info.split()[10]) and int(map_data[i][1]) >= MIN_READ_COUNT_THRESHOLD:
                    count += 1
        elif line_info.split()[5] == "-":
            for i in range(0, len(map_data)):
                if map_data[i][2] == line_info.split()[2] and int(map_data[i][3]) >= int(line_info.split()[9])\
                        and int(map_data[i][4]) <= int(line_info.split()[10]) and int(map_data[i][1]) >= MIN_READ_COUNT_THRESHOLD:
                    count += 1
        if count == 0:
            continue

    # write putative precursor to the output file
    result_count += 1
    output_mature.write("Name\tRead_Count\tChr_Name\tMature_Start\tMature_End\tPos\tSeq\tMFE\tNorm_MFE\tPrec_Start\tPrec_End\n")
    output_mature.write(line_info+"\n"+line_seq+"\n"+line_db+"\n")
    output_mature.write('*'*start_5p+line_seq[start_5p:end_5p]+'*'*(len(line_seq)-end_5p)+"\n")
    output_mature.write('*'*start_3p+line_seq[start_3p:end_3p]+'*'*(len(line_seq)-end_3p)+"\n")
    if line_info.split()[5] == "+":
        for i in range(0, len(map_data)):
            if map_data[i][2] == line_info.split()[2] and int(map_data[i][3]) >= int(line_info.split()[9])\
                    and int(map_data[i][4]) <= int(line_info.split()[10]) and int(map_data[i][1]) >= MIN_READ_COUNT_THRESHOLD:
                output_mature.write('-'*(int(map_data[i][3])-int(line_info.split()[9]))+str(map_data[i][6])+
                                    '-'*(int(line_info.split()[10])-int(map_data[i][4]))+'\t'+
                                    str(map_data[i][0])+'\t'+str(map_data[i][1])+'\n')
    elif line_info.split()[5] == "-":
        for i in range(0, len(map_data)):
            if map_data[i][2] == line_info.split()[2] and int(map_data[i][3]) >= int(line_info.split()[9])\
                    and int(map_data[i][4]) <= int(line_info.split()[10]) and int(map_data[i][1]) >= MIN_READ_COUNT_THRESHOLD:
                output_mature.write('-'*(int(line_info.split()[10])-int(map_data[i][4]))+str(map_data[i][6])+
                                    '-'*(int(map_data[i][3])-int(line_info.split()[9]))+'\t'+
                                    str(map_data[i][0])+'\t'+str(map_data[i][1])+'\n')
    output_mature.write("\n")
print("Done : "+str(result_count)+" miRNA found, See result_mature.txt for details")

##################################### main script end #####################################
