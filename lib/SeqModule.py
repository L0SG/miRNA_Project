__author__ = 'L0SG'
import re

# Pre-Alpha version for testing purpose : does not have proper min, max implementation


def find_perfect_match(seq, ref, maxmult):
    # find perfect match index
    list2d = []
    for i in range(0, len(ref)):
        templist = [temp.span() for temp in re.finditer(seq, str(ref[i]))]
        if templist != [] and len(templist) <= maxmult:
            list2d.append(templist)
        elif templist == []:
            list2d.append(templist)
    return list2d


def convert_bowtie_output(output_bowtie, output_map):
    for line in output_bowtie:
        line_split = line.split()
        if line_split[2] == "-":
            seq = create_star(line_split[5])
        else:
            seq = line_split[5]
        seq_index_end = int(line_split[4])+len(line_split[5])
        output_map.write(line_split[0]+"\t"+line_split[1]+"\t"+line_split[3]+"\t"+
                         line_split[4]+"\t"+str(seq_index_end)+"\t"+line_split[2]+"\t"+seq+"\n")


def create_star(smrna_seq):
    # create reverse complement sequence (i.e. star sequence)
    star_seq = ""
    for i in range(0, len(smrna_seq)):
        if smrna_seq[i] is "A":
            star_seq += "T"
            continue
        elif smrna_seq[i] is "T":
            star_seq += "A"
            continue
        elif smrna_seq[i] is "G":
            star_seq += "C"
            continue
        elif smrna_seq[i] is "C":
            star_seq += "G"
            continue
        else: # N case
            continue
    # reverse the sequence to place 5' in the left
    star_seq_result = star_seq[::-1]
    return star_seq_result


def count_generator(ref_name_list, output_map):
    count_list_pos=[]
    count_list_neg=[]
    for i in xrange(0, len(ref_name_list)):
        count_list_pos.append({})
        count_list_neg.append({})

    for line in output_map:
        line_split = line.split()
        name_list_index = ref_name_list.index(line_split[2])
        count = int(line_split[1])
        index_pos = int(line_split[3])
        index_neg = int(line_split[4])
        sign = line_split[5]
        if sign == "+":
            if not index_pos in count_list_pos[name_list_index]:
                count_list_pos[name_list_index][index_pos] = count
            else:
                count_list_pos[name_list_index][index_pos] += count
        elif sign == "-":
            if not index_neg in count_list_neg[name_list_index]:
                count_list_neg[name_list_index][index_neg] = count
            else:
                count_list_neg[name_list_index][index_neg] += count
    return count_list_pos, count_list_neg


def check_conserved_seq(line_info, line_seq, blastn_path, mirbase_path, arm_length):
    import subprocess
    import os
    from operator import itemgetter

    # run blastn for given line_info
    updated_flag = False
    seq_name = line_info.split()[0]
    seq_seq = line_info.split()[6]
    command = blastn_path
    mirbase = mirbase_path
    command = command + ' -task blastn -db "' + mirbase + '" -outfmt 6 -strand plus -num_alignments 10 -ungapped'
    blastn = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, err = blastn.communicate(seq_seq)
    output = [x.split('\t') for x in output.split('\n')][:-1]
    # if many same matches are detected, annotate it with the matched information
    # get top matches
    top_num = 3
    output_top_check = map(itemgetter(0), output)
    if len(output_top_check) < top_num:
        return line_info, updated_flag
    output_name = map(itemgetter(1), output)[:top_num]
    output_query_start = map(int, map(itemgetter(6), output)[:top_num])
    output_query_end = map(int, map(itemgetter(7), output)[:top_num])
    output_target_start = map(int, map(itemgetter(8), output)[:top_num])
    output_target_end = map(int, map(itemgetter(9), output)[:top_num])

    # check if all value for each lists are same
    if check_same(output_query_start, output_query_end, output_target_start, output_target_end) is False:
        return line_info, updated_flag

    if output_target_start[0] != 1:
        return line_info, updated_flag

    # update line information
    line_split = line_info.split()
    prec_match_start = line_seq.find(seq_seq)
    prec_match_end = prec_match_start + len(seq_seq)

    switch_start = output_target_start[0]-output_query_start[0]
    switch_end = output_target_end[0]-output_query_end[0]
    mature_start_updated = int(line_split[3]) + switch_start
    mature_end_updated = int(line_split[4]) + switch_end
    mature_seq_updated = line_seq[prec_match_start+switch_start:prec_match_end+switch_end]
    # index handling (for RARE case of extremely short arm length specified by user)
    if int(line_split[9]) <= mature_start_updated and int(line_split[10]) >= mature_end_updated:
        line_split[0] = 'xxx-'+output_name[0][4:] # heuristic implementation, need to be improved later
        line_split[3] = str(mature_start_updated)
        line_split[4] = str(mature_end_updated)
        line_split[6] = mature_seq_updated
        line_info = '\t'.join(line_split)
        updated_flag = True

    return line_info, updated_flag


def check_same(*args):
    for index, check_list in enumerate(args):
        if check_list.count(check_list[0]) != len(check_list):
            return False
    return True


def score_seq(target_5p, target_3p):
    i_5p = 0
    i_3p = 0
    max_serial_mismatch = 0
    max_mult_mismatch = 0
    max_serial_bulge = 0
    max_mult_bulge = 0

    while 1:
        if i_5p == len(target_5p) -1 or i_3p == len(target_3p) -1:
            break
        count_5p = 0
        count_3p = 0
        if target_5p[i_5p]=="(" and target_3p[i_3p]==")":
            if i_5p < len(target_5p)-1:
                i_5p += 1
            if i_3p < len(target_3p)-1:
                i_3p += 1
            continue
        else:
            while target_5p[i_5p] == ".":
                if i_5p == len(target_5p)-1:
                    break
                i_5p += 1
                count_5p += 1
            while target_3p[i_3p] == ".":
                if i_3p == len(target_3p)-1:
                    break
                i_3p += 1
                count_3p += 1
            if count_5p == count_3p:  # mismatch case
                if max_serial_mismatch < count_5p:
                    max_serial_mismatch = count_5p
                max_mult_mismatch += 1
                continue
            elif count_5p == 0:  # 3p bulge case
                if max_serial_bulge < count_3p:
                    max_serial_bulge = count_3p
                max_mult_bulge += 1
                continue
            elif count_3p == 0:  # 5p bulge case
                if max_serial_bulge < count_5p:
                    max_serial_bulge = count_5p
                max_mult_bulge += 1
                continue
            elif count_5p > count_3p:  # mismatch and 5p bulge case : consider it as 5p bulge case
                if max_serial_bulge < count_5p:
                    max_serial_bulge = count_5p
                max_mult_bulge += 1
                continue
            elif count_5p < count_3p:  # mismatch and 3p bulge case : consider it as 3p bulge case
                if max_serial_bulge < count_3p:
                    max_serial_bulge = count_3p
                max_mult_bulge += 1
                continue
    return max_serial_mismatch, max_mult_mismatch, max_serial_bulge, max_mult_bulge


def check_no_read_prec(line_info, map_data, MIN_READ_COUNT_THRESHOLD):
    count = 0
    if line_info.split()[5] == "+":
        for i in range(0, len(map_data)):
            if line_info.split()[2] == map_data[i][2] and int(line_info.split()[9]) <= int(map_data[i][3]) \
                    and int(line_info.split()[10]) >= int(map_data[i][4]) and MIN_READ_COUNT_THRESHOLD <= int(map_data[i][1]):
                count += 1
    elif line_info.split()[5] == "-":
        for i in range(0, len(map_data)):
            if line_info.split()[2] == map_data[i][2] and int(line_info.split()[9]) <= int(map_data[i][3]) \
                    and int(line_info.split()[10]) >= int(map_data[i][4]) and MIN_READ_COUNT_THRESHOLD <= int(map_data[i][1]):
                count += 1
    if count == 0:
        return True
    return False


def star_identifier_v2(precursor_db, mature_min_len, mature_max_len, max_serial_mismatch, max_mult_mismatch, max_serial_bulge, max_mult_bulge):
    start_5p = 0
    end_5p = 0
    start_3p = 0
    end_3p = 0
    norm_score = 999
    for i in range(0, len(precursor_db) - mature_min_len):  # search potential -5p seq
        for j in range(mature_min_len, mature_max_len):
            if i+j >= len(precursor_db):
                continue
            target_5p = precursor_db[i:i+j]
            if ")" in target_5p:
                continue
            # consider 2nt overhang
            num_bracket = target_5p[0:len(target_5p)-2].count("(")
            # find target_3p corresponding to db notation of the precursor
            iter_hairpin = i+j-2
            open_count = 0
            close_count = 0
            index_error_flag = 0
            while 1:
                if iter_hairpin == len(precursor_db):
                    index_error_flag = 1
                    break
                if precursor_db[iter_hairpin] == "(":
                    open_count += 1
                    iter_hairpin += 1
                elif precursor_db[iter_hairpin] == ".":
                    iter_hairpin += 1
                elif precursor_db[iter_hairpin] == ")":
                    break
            while open_count != close_count:
                if iter_hairpin == len(precursor_db):
                    index_error_flag = 1
                    break
                if precursor_db[iter_hairpin] == "(":
                    open_count += 1
                    iter_hairpin += 1
                elif precursor_db[iter_hairpin] == ".":
                    iter_hairpin += 1
                elif precursor_db[iter_hairpin] == ")":
                    close_count += 1
                    iter_hairpin += 1
            if index_error_flag == 1:
                continue
            # iter_hairpin is now the candidate start index of target_3p
            # db of this index could be ".", so keep searching until another ")" appears
            # if ")" appears, this index violates db notations of the precursor
            stop_flag = 0
            if iter_hairpin == len(precursor_db):
                continue
            while precursor_db[iter_hairpin] == "." or ")":
                if precursor_db[iter_hairpin] == ")":
                    stop_flag = 1
                for l in range(mature_min_len, mature_max_len):
                    if iter_hairpin+l < len(precursor_db):
                        target_3p_pre = precursor_db[iter_hairpin:iter_hairpin+l]
                        target_3p = target_3p_pre[::-1]
                        if "(" in target_3p:
                            continue
                        # consider 2nt overhang
                        if num_bracket != target_3p[2:len(target_3p)].count(")"):
                            continue
                        # consider 3p 2nt overhang and score target seq
                        msm, mmm, msb, mmb = score_seq(target_5p[0:len(target_5p)-2], target_3p[2:len(target_3p)])
                        if msm <= max_serial_mismatch and mmm <= max_mult_mismatch and msb <= max_serial_bulge and mmb <= max_mult_bulge:
                            if norm_score > msm + mmm + msb + mmb:
                                start_5p = i
                                end_5p = i+j
                                start_3p = iter_hairpin
                                end_3p = iter_hairpin+l
                                norm_score = msm + mmm + msb + mmb
                iter_hairpin += 1
                if stop_flag == 1:
                    break
                if iter_hairpin == len(precursor_db):
                    break
    return start_5p, end_5p, start_3p, end_3p


def generate_output_form(line_info, line_seq, line_db, start_5p, start_3p, end_5p, end_3p, map_data, MIN_READ_COUNT_THRESHOLD):
    output_form = []
    output_form.append(str((line_info + "\n" + line_seq + "\n" + line_db + "\n")))
    output_form.append(str(('*' * start_5p + line_seq[start_5p:end_5p] + '*' * (len(line_seq) - end_5p) + "\n")))
    output_form.append(str(('*' * start_3p + line_seq[start_3p:end_3p] + '*' * (len(line_seq) - end_3p) + "\n")))
    if line_info.split()[5] == "+":
        for i in range(0, len(map_data)):
            if line_info.split()[2] == map_data[i][2]\
                    and int(line_info.split()[9]) <= int(map_data[i][3])\
                    and int(line_info.split()[10]) >= int(map_data[i][4])\
                    and MIN_READ_COUNT_THRESHOLD <= int(map_data[i][1]):
                # filter out reverse complement case
                space_left = int(map_data[i][3]) - int(line_info.split()[9])
                space_right = int(line_info.split()[10]) - int(map_data[i][4])
                seq_align = map_data[i][6]
                seq_ref = line_seq[space_left:space_left + len(seq_align)]
                if seq_align != seq_ref:
                    continue
                output_form.append(str(('-' * space_left + str(map_data[i][6]) +
                                        '-' * space_right + '\t' +
                                        str(map_data[i][0]) + '\t' + str(map_data[i][1]) + '\n')))
    elif line_info.split()[5] == "-":
        for i in range(0, len(map_data)):
            if line_info.split()[2] == map_data[i][2]\
                    and int(line_info.split()[9]) <= int(map_data[i][3])\
                    and int(line_info.split()[10]) >= int(map_data[i][4])\
                    and MIN_READ_COUNT_THRESHOLD <= int(map_data[i][1]):
                # filter out reverse complement case
                space_left = int(line_info.split()[10]) - int(map_data[i][4])
                space_right = int(map_data[i][3]) - int(line_info.split()[9])
                seq_align = map_data[i][6]
                seq_ref = line_seq[space_left:space_left + len(seq_align)]
                if seq_align != seq_ref:
                    continue
                output_form.append(str(('-' * space_left + str(map_data[i][6]) +
                                        '-' * space_right + '\t' +
                                        str(map_data[i][0]) + '\t' + str(map_data[i][1]) + '\n')))
    return output_form


def generate_alignment_form(line_info, line_seq, line_db, map_data, MIN_READ_COUNT_THRESHOLD):
    # this is a variant of original generate_output_form
    # no calculated duplex information, align first
    # get highest reads seq as representative, and search the seq for duplex annotation
    import operator
    output_form = []
    output_form.append(str((line_info + "\n")))
    output_form.append(str((line_seq + "\n")))
    output_form.append(str((line_db + "\n")))
    if line_info.split()[5] == "+":
        for i in range(0, len(map_data)):
            if line_info.split()[2] == map_data[i][2]\
                    and int(line_info.split()[9]) <= int(map_data[i][3])\
                    and int(line_info.split()[10]) >= int(map_data[i][4])\
                    and MIN_READ_COUNT_THRESHOLD <= int(map_data[i][1]):
                # filter out reverse complement case
                space_left = int(map_data[i][3]) - int(line_info.split()[9])
                space_right = int(line_info.split()[10]) - int(map_data[i][4])
                seq_align = map_data[i][6]
                seq_ref = line_seq[space_left:space_left + len(seq_align)]
                if seq_align != seq_ref:
                    continue
                output_form.append(str(('-' * space_left + str(map_data[i][6]) +
                                        '-' * space_right + '\t' +
                                        str(map_data[i][0]) + '\t' + str(map_data[i][1]) + '\n')))
    elif line_info.split()[5] == "-":
        for i in range(0, len(map_data)):
            if line_info.split()[2] == map_data[i][2]\
                    and int(line_info.split()[9]) <= int(map_data[i][3])\
                    and int(line_info.split()[10]) >= int(map_data[i][4])\
                    and MIN_READ_COUNT_THRESHOLD <= int(map_data[i][1]):
                # filter out reverse complement case
                space_left = int(line_info.split()[10]) - int(map_data[i][4])
                space_right = int(map_data[i][3]) - int(line_info.split()[9])
                seq_align = map_data[i][6]
                seq_ref = line_seq[space_left:space_left + len(seq_align)]
                if seq_align != seq_ref:
                    continue
                output_form.append(str(('-' * space_left + str(map_data[i][6]) +
                                        '-' * space_right + '\t' +
                                        str(map_data[i][0]) + '\t' + str(map_data[i][1]) + '\n')))
    aligned_list = output_form[3:]
    for i in xrange(0, len(aligned_list)):
        aligned_list[i] = aligned_list[i].split()
        aligned_list[i][2] = int(aligned_list[i][2])
    aligned_list.sort(key=operator.itemgetter(2), reverse=True)
    seq_rep = aligned_list[0][0].strip('-')
    idx_start = aligned_list[0][0].index(seq_rep)
    idx_end = idx_start + len(seq_rep)
    output_info = output_form[0].split()
    if output_info[5] == '+':
        mature_start = int(output_info[9]) + idx_start
        mature_end = mature_start + len(seq_rep)
    elif output_info[5] == '-':
        mature_end = int(output_info[10]) - idx_start
        mature_start = mature_end - len(seq_rep)

    output_info[3] = str(mature_start)
    output_info[4] = str(mature_end)
    output_info[6] = seq_rep
    output_form[0] = '\t'.join(output_info)+'\n'
    return output_form


def star_identifier_v2_conserved(line_info, line_seq, precursor_db, mature_min_len, mature_max_len, max_serial_mismatch, max_mult_mismatch, max_serial_bulge, max_mult_bulge):
    start_5p = 0
    end_5p = 0
    start_3p = 0
    end_3p = 0
    norm_score = 999
    for i in range(0, len(precursor_db) - mature_min_len):  # search potential -5p seq
        for j in range(mature_min_len, mature_max_len):
            if i+j >= len(precursor_db):
                continue
            target_5p = precursor_db[i:i+j]
            if ")" in target_5p:
                continue
            # consider 2nt overhang
            num_bracket = target_5p[0:len(target_5p)-2].count("(")
            # find target_3p corresponding to db notation of the precursor
            iter_hairpin = i+j-2
            open_count = 0
            close_count = 0
            index_error_flag = 0
            while 1:
                if iter_hairpin == len(precursor_db):
                    index_error_flag = 1
                    break
                if precursor_db[iter_hairpin] == "(":
                    open_count += 1
                    iter_hairpin += 1
                elif precursor_db[iter_hairpin] == ".":
                    iter_hairpin += 1
                elif precursor_db[iter_hairpin] == ")":
                    break
            while open_count != close_count:
                if iter_hairpin == len(precursor_db):
                    index_error_flag = 1
                    break
                if precursor_db[iter_hairpin] == "(":
                    open_count += 1
                    iter_hairpin += 1
                elif precursor_db[iter_hairpin] == ".":
                    iter_hairpin += 1
                elif precursor_db[iter_hairpin] == ")":
                    close_count += 1
                    iter_hairpin += 1
            if index_error_flag == 1:
                continue
            # iter_hairpin is now the candidate start index of target_3p
            # db of this index could be ".", so keep searching until another ")" appears
            # if ")" appears, this index violates db notations of the precursor
            stop_flag = 0
            if iter_hairpin == len(precursor_db):
                continue
            while precursor_db[iter_hairpin] == "." or ")":
                if precursor_db[iter_hairpin] == ")":
                    stop_flag = 1
                for l in range(mature_min_len, mature_max_len):
                    if iter_hairpin+l < len(precursor_db):
                        target_3p_pre = precursor_db[iter_hairpin:iter_hairpin+l]
                        target_3p = target_3p_pre[::-1]
                        if "(" in target_3p:
                            continue
                        # consider 2nt overhang
                        if num_bracket != target_3p[2:len(target_3p)].count(")"):
                            continue
                        # NEW ad-hoc lines for conserved seq
                        line_info_split = line_info.split()
                        if line_info_split[5] == '+':
                            cmp_start = int(line_info_split[3])-int(line_info_split[9])
                            cmp_end = cmp_start + len(line_info_split[6])
                        elif line_info_split[5] == '-':
                            cmp_end = len(precursor_db) - int(line_info_split[3])-int(line_info_split[9])
                            cmp_start = cmp_end - len(line_info_split[6])
                        # 5p match case
                        if cmp_start == i and cmp_end == i+j:
                            start_5p = i
                            end_5p = i + j
                            start_3p = iter_hairpin
                            end_3p = iter_hairpin + l
                            stop_flag = 1
                        # 3p match case
                        if cmp_start == iter_hairpin and cmp_end == iter_hairpin+l:
                            start_5p = i
                            end_5p = i + j
                            start_3p = iter_hairpin
                            end_3p = iter_hairpin + l
                            stop_flag = 1
                        """
                        # consider 3p 2nt overhang and score target seq
                        msm, mmm, msb, mmb = score_seq(target_5p[0:len(target_5p)-2], target_3p[2:len(target_3p)])
                        if msm <= max_serial_mismatch and mmm <= max_mult_mismatch and msb <= max_serial_bulge and mmb <= max_mult_bulge:
                            if norm_score > msm + mmm + msb + mmb:
                                start_5p = i
                                end_5p = i+j
                                start_3p = iter_hairpin
                                end_3p = iter_hairpin+l
                                norm_score = msm + mmm + msb + mmb
                        """
                iter_hairpin += 1
                if stop_flag == 1:
                    break
                if iter_hairpin == len(precursor_db):
                    break
    return start_5p, end_5p, start_3p, end_3p


def builtin_map_generator(smrna_file, ref_name_list, ref_seq_list,
                          output_map, ref_count_list_pos, ref_count_list_neg,
                          NUM_THREADS, MAX_MULTIPLE_LOCI):
    import multiprocessing
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
            smrna_star = create_star(smrna_seq)  # seq for neg match
            # decide whether smrna_seq perfectly match genome seq
            # return type : 2D list
            # each matched case returns zero based begin, end index (2-element tuple)
            # character of end index is NOT contained in seq
            # [0][0] : first genome, first matched indices, [0][1] : first genome, second matched indices
            # [1][0] : second genome, first matched indices, [1][1] : second matched indices ...
            pm_index_list_pos = find_perfect_match(smrna_seq, ref_seq_list, MAX_MULTIPLE_LOCI)
            pm_index_list_neg = find_perfect_match(smrna_star, ref_seq_list, MAX_MULTIPLE_LOCI)
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
        numlines = 5000
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
                if int(ref_count_list_pos_result[i][j][0])>=len(ref_count_list_pos[i]):
                    print int(ref_count_list_pos_result[i][j][0])
                    print "error pos"
                ref_count_list_pos[i][int(ref_count_list_pos_result[i][j][0])] += ref_count_list_pos_result[i][j][1]
            for j in range(0, len(ref_count_list_neg_result[i])):
                if int(ref_count_list_neg_result[i][j][0])>=len(ref_count_list_neg[i]):
                    print int(ref_count_list_neg_result[i][j][0])
                    print "error neg"
                ref_count_list_neg[i][int(ref_count_list_neg_result[i][j][0])] += ref_count_list_neg_result[i][j][1]
    output_map.seek(0, 0)






