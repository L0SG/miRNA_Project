__author__ = 'L0SG'
import re

# Pre-Alpha version for testing purpose : does not have proper min, max implementation


def find_perfect_match(seq, ref, min, max, maxmult):
    # find perfect match index
    list2d = []
    for i in range(0, len(ref)):
        templist = [temp.span() for temp in re.finditer(seq, str(ref[i]))]
        if templist != [] and len(templist) <= maxmult:
            list2d.append(templist)
    return list2d


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

# Outdated, do not use
"""
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
"""

def star_identifier_v2(precursor_db, mature_min_len, mature_max_len, max_serial_mismatch, max_mult_mismatch, max_serial_bulge, max_mult_bulge):
    start_5p = 0
    end_5p = 0
    start_3p = 0
    end_3p = 0
    norm_score = 999
    for i in range(0, len(precursor_db)/2 - mature_min_len):  # search potential -5p seq form zero to half the sequence length
        for j in range(mature_min_len, mature_max_len):
            target_5p = precursor_db[i:i+j]
            if ")" in target_5p:
                continue
            # consider 2nt overhang
            num_bracket = target_5p[0:len(target_5p)-2].count("(")
            search_3p_seed = 5
            for k in range(len(precursor_db)-(i+j)-search_3p_seed, len(precursor_db)-(i+j)+search_3p_seed):
                for l in range(mature_min_len, mature_max_len):
                    if k+l < len(precursor_db):
                        target_3p_pre = precursor_db[k:k+l]
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
                                start_3p = k
                                end_3p = k+l
                                norm_score = msm + mmm + msb + mmb
    return start_5p, end_5p, start_3p, end_3p


def star_identifier(target_star_seq, search_seq, target_ref_seq):
    if target_star_seq == search_seq:
        return 1
    else:
        return -1


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




