# outdated and not used functions
def convert_dump_to_list(ref_count_dump, ref_count_list):
    for i in xrange(0, len(ref_count_dump)):
        for j in ref_count_dump[i]:
            value = ref_count_dump[i][j]
            ref_count_list[i][j] += value


def star_identifier(target_star_seq, search_seq, target_ref_seq):
    if target_star_seq == search_seq:
        return 1
    else:
        return -1


def star_identifier(precursor_db, mature_min_len, mature_max_len, max_serial_mismatch, max_mult_mismatch, max_serial_bulge, max_mult_bulge):
    # find star seq given a dot-bracket seq
    # step 1 : find star seq with min_length, matching all criteria
    # step 2 : extend length of step 1 seq (~max_length) and normalize score
    # step 3 : take the one with best norm.score
    start_5p = 0
    end_5p = 0
    start_3p = 0
    end_3p = 0
    norm_score = 1
    for i in range(0, len(precursor_db)/2 - mature_min_len):  # search potential -5p seq form zero to half the sequence length
        print(i)
        target_mature_init = precursor_db[i:i+mature_min_len]
        if ")" in target_mature_init:
            continue
        for j in range(i+mature_min_len, len(precursor_db)-mature_min_len):  # search potential -3p seq from next to -5p seq to the end of seq
            target_star_init_pre = precursor_db[j:j+mature_min_len]
            target_star_init = target_star_init_pre[::-1]
            if "(" in target_star_init:
                continue
            # consider 3' 2nt overhang and score target seq
            msm, mmm, msb, mmb = score_seq(target_mature_init[0:len(target_mature_init)-2], target_star_init[2:len(target_star_init)])
            if msm <= max_serial_mismatch and mmm <= max_mult_mismatch and msb <= max_serial_bulge and mmb <= max_mult_bulge:
                if norm_score > (msm + mmm + msb + mmb)/len(target_mature_init):
                    start_5p = i
                    end_5p = i+mature_min_len
                    start_3p = j
                    end_3p = j+mature_min_len
                    norm_score = (msm + mmm + msb + mmb)/len(target_mature_init)
                for k in range(0, mature_max_len-mature_min_len):
                    target_mature_extend = precursor_db[i:i+mature_min_len+k]
                    target_star_extend = precursor_db[j:j+mature_min_len+k]
                    msm, mmm, msb, mmb = score_seq(target_mature_extend[0:len(target_mature_init)-2], target_star_extend[2:len(target_star_init)])
                    if norm_score > (msm + mmm + msb + mmb)/len(target_mature_init):
                        start_5p = i
                        end_5p = i+mature_min_len+k
                        start_3p = j
                        end_3p = j+mature_min_len+k
                        norm_score = (msm + mmm + msb + mmb)/len(target_mature_init)
    return start_5p, end_5p, start_3p, end_3p


def find_star(pm_index_list, ref_seq_list, threshold):
    # find star sequence, brute force
    list3d = []
    for i in range(0, len(ref_seq_list)):
        target_ref_seq = ref_seq_list[i]
        list3d.append([])
        for j in range(0, len(pm_index_list[i])):
            start_point = pm_index_list[i][j][0]
            end_point = pm_index_list[i][j][1]
            target_smrna_seq = ref_seq_list[i][start_point:end_point]
            target_star_seq = create_star(target_smrna_seq)
            valid_star_list = []

            # search positive direction, prevent index out of range error
            if end_point + threshold <= len(target_ref_seq):
                for k in range(len(target_smrna_seq), threshold):
                    search_seq_pos = target_ref_seq[start_point+k:end_point+k]
                    if star_identifier(target_star_seq, search_seq_pos, target_ref_seq) == 1:
                        valid_star_list.append((start_point+k, end_point+k))

            # search positive direction, prevent index out of range error
            if start_point - threshold >= 0:
                for k in range(len(target_smrna_seq), threshold):
                    search_seq_pos = target_ref_seq[start_point-k:end_point-k]
                    if star_identifier(target_star_seq, search_seq_pos, target_ref_seq) == 1:
                        valid_star_list.append((start_point-k, end_point-k))
            if(valid_star_list != []):
                print(valid_star_list)
            list3d[i].append(valid_star_list)

    return list3d


def find_precursor(pm_index_list, star_index_list):
    return


"""
print("Generating map file from seq library...")
while 1:
    line = smrna_file.readline().strip()
    if line == "":
        break
    else:
        smrna_info = line.split()
        smrna_seq = smrna_file.readline().strip()  # seq for pos match
        smrna_star = SeqModule.create_star(smrna_seq)  # seq for neg match
        # decide whether smrna_seq perfectly match genome seq
        # return type : 2D list
        # each matched case returns zero based begin, end index (2-element tuple)
        # character of end index is NOT contained in seq
        # [0][0] : first genome, first matched indices, [0][1] : first genome, second matched indices
        # [1][0] : second genome, first matched indices, [1][1] : second matched indices ...
        pm_index_list_pos = SeqModule.find_perfect_match(smrna_seq, ref_seq_list, MATURE_MIN_LEN, MATURE_MAX_LEN, MAX_MULTIPLE_LOCI)
        pm_index_list_neg = SeqModule.find_perfect_match(smrna_star, ref_seq_list, MATURE_MIN_LEN, MATURE_MAX_LEN, MAX_MULTIPLE_LOCI)
        for i in range(0, len(pm_index_list_pos)):
            for j in range(0, len(pm_index_list_pos[i])):
                output_map.write(smrna_info[0]+"\t"+smrna_info[1]+"\t"+ref_name_list[i]+"\t"+str(pm_index_list_pos[i][j][0])+"\t"+str(pm_index_list_pos[i][j][1])+"\t"+"+"+"\t"+smrna_seq+"\n")
                ref_count_list_pos[i][pm_index_list_pos[i][j][0]] += int(smrna_info[1])
        for i in range(0, len(pm_index_list_neg)):
            for j in range(0, len(pm_index_list_neg[i])):
                output_map.write(smrna_info[0]+"\t"+smrna_info[1]+"\t"+ref_name_list[i]+"\t"+str(pm_index_list_neg[i][j][0])+"\t"+str(pm_index_list_neg[i][j][1])+"\t"+"-"+"\t"+smrna_seq+"\n")
                ref_count_list_neg[i][pm_index_list_neg[i][j][1]] += int(smrna_info[1])
output_map.seek(0, 0)
json.dump(ref_count_list_pos, output_count_pos)
json.dump(ref_count_list_neg, output_count_neg)
output_count_pos.seek(0, 0)
output_count_neg.seek(0, 0)
print("Generating Done")
"""
"""
print("Calculating MFE of putative precursors with RNAfold...")
# use RNAfold to calculate MFE and select putative precursor
output_precursor.write("Name\tRead_Count\tChr_Name\tMature_Start\tMature_End\tPos\tMFE\tNorm_MFE\tPrec_Start\tPrec_End\n")
while 1:
    line = output_map.readline().strip()
    line_split = line.split()
    if line == "":
        break
    if len(line_split) !=  7:
        continue

    # Screen for Drosha / Dicer cutting sites (Inspired by miREAP)
    qualified_flag = 1
#    name_list_index = ref_name_list.index(line_split[2])
    name_list_index = 0
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
            if ref_count_list_pos[name_list_index][int(line_split[3])-i] > count or ref_count_list_pos[name_list_index][int(line_split[3])+i] > count:
                qualified_flag = 0
                break
            count_region += ref_count_list_pos[name_list_index][int(line_split[3])-i]
            count_region += ref_count_list_pos[name_list_index][int(line_split[3])+i]
            if i < 3:
                count_sites += ref_count_list_pos[name_list_index][int(line_split[3])-i]
                count_sites += ref_count_list_pos[name_list_index][int(line_split[3])+i]
        if count_sites/count_region < 0.9 or count/count_sites < 0.5:
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
            if ref_count_list_neg[name_list_index][int(line_split[4])-i] > count or ref_count_list_neg[name_list_index][int(line_split[4])+i] > count:
                qualified_flag = 0
                break
            count_region += ref_count_list_neg[name_list_index][int(line_split[4])-i]
            count_region += ref_count_list_neg[name_list_index][int(line_split[4])+i]
            if i < 3:
                count_sites += ref_count_list_neg[name_list_index][int(line_split[4])-i]
                count_sites += ref_count_list_neg[name_list_index][int(line_split[4])+i]
        if count_sites/count_region < 0.9 or count/count_sites < 0.5:
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
    for k in range(0, len(ref_seq_list)):  # reference sequence list loop
        for i in range(0, ARM_EXTEND_THRESHOLD):  # arm extension loop
            for j in range(int(line_split[4])-int(line_split[3]), DISTANCE_THRESHOLD+int(line_split[4])-int(line_split[3])):  # distance loop

                # Assuming -5p mature sequence
                start = int(line_split[3])-i
                end = int(line_split[4])+j
                if start >= 0 and end < len(ref_seq_list[k]):  # continue only if both indices are valid
                    if line_split[5] == "+":
                        rna_fold_seq = ref_seq_list[k][start:end]
                    elif line_split[5] == "-":
                        rna_fold_seq = SeqModule.create_star(ref_seq_list[k][start:end])
                    if "N" in rna_fold_seq:
                        continue
                    rnafold = subprocess.Popen(["RNAfold", "--noconv", "-d2", "--noPS", "--noLP"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                    output = rnafold.communicate(rna_fold_seq)[0].split()
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
                start = int(line_split[3])-j
                end = int(line_split[4])+i
                if start >= 0 and end < len(ref_seq_list[k]):   # continue only if both indices are valid
                    if line_split[5] == "+":
                        rna_fold_seq = ref_seq_list[k][start:end]
                    elif line_split[5] == "-":
                        rna_fold_seq = SeqModule.create_star(ref_seq_list[k][start:end])
                    if "N" in rna_fold_seq:
                        continue
                    rnafold = subprocess.Popen(["RNAfold", "--noconv", "-d2", "--noPS", "--noLP"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                    output = rnafold.communicate(rna_fold_seq)[0].split()
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
        output_precursor.write(line+"\t"+str(pc_abs_energy)+"\t"+str(pc_norm_abs_energy)+"\t"+str(pc_start)+"\t"+str(pc_end)+"\n")
        output_precursor.write(pc_seq+"\n"+pc_structure+"\n")
        print("RNAfold loop")
output_precursor.seek(0, 0)
print("Calculating Done")
"""
"""
        # Perfect match finding results in 4 cases
        # case 1 : pos(+) direction only
        if pm_index_list_pos != [] and pm_index_list_neg == []:
            star_index_list = SeqModule.find_star(pm_index_list_pos, ref_seq_list, FIND_STAR_THRESHOLD)
            precursor_list = SeqModule.find_precursor(pm_index_list_pos, star_index_list)
            FileIOModule.write_mature(output_mature, smrna_info, pm_index_list_pos, star_index_list, "pos")
            FileIOModule.write_precursor(output_precursor, smrna_info, precursor_list, "pos")
            continue
        # case 2 : neg(-) direction only
        elif pm_index_list_pos == [] and pm_index_list_neg != []:
            star_index_list = SeqModule.find_star(pm_index_list_neg, ref_seq_list, FIND_STAR_THRESHOLD)
            FileIOModule.write_mature(output_mature, smrna_info, pm_index_list_neg, star_index_list, "neg")
            FileIOModule.write_precursor(output_precursor, smrna_info, precursor_list, "neg")
            continue
        # case 3 : smrna_seq have both pos(+) and neg(-) matches
        elif pm_index_list_pos != [] and pm_index_list_neg != []:
            print("temp")
            continue
        # case 4 : smrna seq did not match
        else:
            continue
"""