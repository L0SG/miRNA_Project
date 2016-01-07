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
            # find target_3p corresponding to db notation of the precursor
            iter_hairpin = i+j-2
            open_count = 0
            close_count = 0
            while 1:
                if precursor_db[iter_hairpin] == "(":
                    open_count += 1
                    iter_hairpin += 1
                elif precursor_db[iter_hairpin] == ".":
                    iter_hairpin += 1
                elif precursor_db[iter_hairpin] == ")":
                    break
            while open_count != close_count:
                if precursor_db[iter_hairpin] == "(":
                    open_count += 1
                    iter_hairpin += 1
                elif precursor_db[iter_hairpin] == ".":
                    iter_hairpin += 1
                elif precursor_db[iter_hairpin] == ")":
                    close_count += 1
                    iter_hairpin += 1
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







