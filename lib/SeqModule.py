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


def convert_dump_to_list(ref_count_dump, ref_count_list):
    for i in xrange(0, len(ref_count_dump)):
        for j in ref_count_dump[i]:
            value = ref_count_dump[i][j]
            ref_count_list[i][j] += value

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






