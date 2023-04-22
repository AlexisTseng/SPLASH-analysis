import sys
import pickle
import os
import math
import numpy as np
from matplotlib import pyplot as plt
import csv
from pycirclize import Circos
#import seaborn as sns #version = 12.2
import pandas as pd
import subprocess
import re

#change window size at windgen, import different picklefile, rename whole_wind directory


def main(genomefasta,pcl1,pcl0,samp_name,Agv,Mgv):
    #, pcl1, pcl0, samp_name
    totreadbp_list = list()
    all_SRI_cd_listlist = list()
    all_intra_cd_listlist = list()
    all_int_cd_listlist = list()
    heat_exp_dir_listlist = list()
    cir_exp_dir_listlist = list()
    exp_name_listlist = list()
    Agv=float(Agv)
    Mgv=float(Mgv)

    #create new directories for storing data
    cd = os.getcwd()
    #print(cd)

#either:
    whole_wind_dir = os.path.join(cd, "results", "window100")
    if not os.path.isdir(whole_wind_dir):
        os.mkdir(whole_wind_dir)
#or:
#    whole_wind_dir = os.path.join(cd, "results", "window200")
#    if not os.path.isdir(whole_wind_dir):
#        os.mkdir(whole_wind_dir)


    samp_wind_dir = os.path.join(whole_wind_dir,samp_name)
    if not os.path.isdir(samp_wind_dir):
        os.mkdir(samp_wind_dir)
    AgMg_dir = os.path.join(samp_wind_dir,f"Ag_{Agv}_Mg_{Mgv}")
    if not os.path.isdir(AgMg_dir):
        os.mkdir(AgMg_dir)
    AgMg_rawdir = os.path.join(samp_wind_dir,f"raw_Ag_{Agv}_Mg_{Mgv}")
    if not os.path.isdir(AgMg_rawdir):
        os.mkdir(AgMg_rawdir)
    ppmdir = os.path.join(samp_wind_dir,"ppm")
    if not os.path.isdir(ppmdir):
        os.mkdir(ppmdir)
    heatmap_dir = os.path.join(whole_wind_dir,"heatmap")
    if not os.path.isdir(heatmap_dir):
        os.mkdir(heatmap_dir)
    circos_dir = os.path.join(whole_wind_dir,"circos_plot")
    if not os.path.isdir(circos_dir):
        os.mkdir(circos_dir)
    
    print("1. Reading Genome File")
    gen_dict = readGenome(genomefasta)

    print("2. Generating Windows")
    segWind_list, all_segname = windgen(gen_dict)
    #print(all_wind_num, all_wind_list)


    print("Produce A count Heatmap")
    polyA_hist(segWind_list, gen_dict,whole_wind_dir)

    exit()
    for pcl in [pcl1,pcl0]:
        print("3. Importing Interaction Picklefile")
        interLRI_list, intraLRI_listdr, SRI_list,exp_name,bam_cnt,junc_cnt,totread_bp, totread_ap, per_inter, per_intra, per_SRI = loadData(pcl)
        totreadbp_list.append(totread_bp)
        #print(SRI_list[0:30])

        exp_dir = os.path.join(heatmap_dir,f"{exp_name}")
        if not os.path.isdir(exp_dir):
            os.mkdir(exp_dir)
        heat_exp_dir_listlist.append(exp_dir)

        cir_exp_dir = os.path.join(circos_dir,f"{exp_name}")
        if not os.path.isdir(cir_exp_dir):
            os.mkdir(cir_exp_dir)
        cir_exp_dir_listlist.append(cir_exp_dir)

        exp_name_listlist.append(exp_name[:-11])

        picklefile = os.path.join(os.path.split(os.path.abspath(pcl))[0], os.path.split(os.path.splitext(os.path.abspath(pcl))[0])[1])

        try:
            ##either:
            SRI_ali_list,intra_ali_list,inter_ali_list = loadData(f"{picklefile}_100window_cnt")
            ##or:
            #SRI_ali_list,intra_ali_list,inter_ali_list = loadData(f"{picklefile}_window_cnt")
        except:
            print("4. Mapping SRI")
            #SRI read all shorter than 200
            #readlen(SRI_list)
            SRI_ali_list = mapSRI(SRI_list,segWind_list)
            
    
            print("5. Mapping intra_LRI")
            #check if all intra read shorter than 200 -> no 
            #readlen_intra(intraLRI_listdr)
            intra_ali_list = mapIntra(intraLRI_listdr,segWind_list)
            
            print("6. Mapping inter_LRI")
            #readlen_inter(interLRI_list)
            inter_ali_list = mapInter(interLRI_list,segWind_list)

            ##either:
            saveData((SRI_ali_list,intra_ali_list,inter_ali_list),f"{picklefile}_100window_cnt")
            ##or:
            #saveData((SRI_ali_list,intra_ali_list,inter_ali_list),f"{picklefile}_window_cnt")
        all_SRI_cd_listlist.append(SRI_ali_list)
        all_int_cd_listlist.append(intra_ali_list)
        all_intra_cd_listlist.append(inter_ali_list)
    SRI_ali_list_1, SRI_ali_list_0 = all_SRI_cd_listlist
    intra_ali_list_1, intra_ali_list_0 = all_intra_cd_listlist
    inter_ali_list_1, inter_ali_list_0 = all_int_cd_listlist
    #print(SRI_ali_list_1)
    all_1 = SRI_ali_list_1 + intra_ali_list_1 + inter_ali_list_1
    all_0 = SRI_ali_list_0 + intra_ali_list_0 + inter_ali_list_0
    #print(type(all_0))
    #print(all_0)
    ### These are InterAlign Objects
    #def __init__(self, name, raw_cd,norm_cd):
    #name (seq,"") or (seqa,seqb), raw count dict, norm count dict
    ####################################################

    exp_dir1, exp_dir0 = heat_exp_dir_listlist
    cir_exp_dir1, cir_exp_dir0 = cir_exp_dir_listlist
    exp_name1, exp_name0 = exp_name_listlist

    print("7. Calculate Normalisation Factor")
    TMM, pair_list,trun_list,all_pr0_SRI_raw,all_pr1_SRI_raw,all_pr1_inter_LRI_raw,all_pr0_inter_LRI_raw,all_pr1_intra_LRI_raw,all_pr0_intra_LRI_raw,all_pr1_raw,all_pr0_raw,all_pr1_ppm,all_pr0_ppm,all_pr1_inter_LRI_ppm,all_pr0_inter_LRI_ppm,all_pr1_intra_LRI_ppm,all_pr0_intra_LRI_ppm,all_pr0_SRI_ppm,all_pr1_SRI_ppm = normf(all_1,all_0,totreadbp_list,Agv,Mgv)
    # pair_list =  class SegInter: (self, name,Mg, Ag, weight,raw1,raw0,norm1,norm0,ppm1,ppm0):#norm1,norm0 read count added later
    #name = seg (pair) name, window (combination) 
    ###########################################################

    #print("8. Normalising")
    #pair_list, all_1,all_0, all_pr1_norm,all_pr0_norm,all_pr1_SRI_norm,all_pr1_LRI_norm,all_pr0_SRI_norm,all_pr0_LRI_norm, all_pr1_inter_LRI_norm,all_pr0_inter_LRI_norm,all_pr1_intra_LRI_norm,all_pr0_intra_LRI_norm = readnorm(TMM, pair_list,all_1,all_0)


    #print("9. Saving summary table")
    #summaryTSV(pair_list,samp_name,whole_wind_dir)
    #trunTSV(trun_list,samp_name,whole_wind_dir)

    #print("10. Producing Fold-change Plots - normalised read")
    #plotLog(all_pr1_SRI_norm,all_pr1_LRI_norm,all_pr0_SRI_norm,all_pr0_LRI_norm,all_pr1_norm,all_pr0_norm,samp_name,samp_wind_dir,Agv,Mgv,AgMg_dir,all_pr1_inter_LRI_norm,all_pr0_inter_LRI_norm,all_pr1_intra_LRI_norm,all_pr0_intra_LRI_norm)

    #print("10. Producing Fold-change Plots - raw read")
    #plotLogRaw(all_pr0_SRI_raw,all_pr1_SRI_raw,all_pr1_inter_LRI_raw,all_pr0_inter_LRI_raw,all_pr1_intra_LRI_raw,all_pr0_intra_LRI_raw,all_pr1_raw,all_pr0_raw,samp_name,samp_wind_dir,Agv,Mgv,AgMg_rawdir)

    #print("10. Producing Fold-change Plots - normalised ppm")
    #plotLogPpm(all_pr1_ppm,all_pr0_ppm,all_pr1_inter_LRI_ppm,all_pr0_inter_LRI_ppm,all_pr1_intra_LRI_ppm,all_pr0_intra_LRI_ppm,all_pr0_SRI_ppm,all_pr1_SRI_ppm,samp_name,samp_wind_dir,ppmdir)

    #print("11. Produce Inter-Seg Interaction Heatmap")
    #interHeatmap1(segWind_list,all_1,exp_dir1)
    #interHeatmap1(segWind_list,all_0,exp_dir0)


    #print("12. Produce test_circos TSV File and generate circos plot")
    #fname1,fname0,fname1w,fname0w = testcircosTSV(pair_list,exp_name1, exp_name0,circos_dir)
    #circosConfig(fname1, circos_dir)
    #ircosConfig(fname0, circos_dir)
    #circosConfig(fname1w, circos_dir)
    #circosConfig(fname0w, circos_dir)
    #doCircos(fname1, circos_dir)
    #doCircos(fname0, circos_dir)
    #doCircos(fname1w, circos_dir)
    #doCircos(fname0w, circos_dir)





def readGenome(genomefasta):
     # load genome fasta file, and return a dictionary
    with open(genomefasta, "r") as genfasta:
        gen_dict = dict()
        seg_name = ""
        for i,line in enumerate(genfasta):
            # this command convert the name ">SA11_S5 WT LC178570.1" into "'>SA11_s5', 'WT', LC..'", i.e. cut at space
            # ">SA11_S3 WT 65433667.7\n"
            line = line.strip()
            # ">SA11_S3 WT 65433667.7"
            if i % 2 == 0:
                # this command removes the > symbol in front of genome segment name (from [1:])
                # ">SA11_S3 WT 65433667.7"
                line_list = line.split()
                # [">SA11_S3" "WT" "65433667.7"]
                seg_name_arrow = line_list[0]
                # ">SA11_S3"
                seg_name = seg_name_arrow[1:]
                # "SA11_S3"
            else:
                gen_dict[seg_name] = line
                #seg_name is a temporary variable which holds name from the previous line
        return gen_dict

def loadData(fname): 
    ## load data with pickle
    with open(f"{fname}.pcl", "r+b") as pcl_in:
        pcl_data = pickle.load(pcl_in)
    return pcl_data

class WindAlign: #segWind_list
    def __init__(self, name, wind_num, wind_list,wind_tick):
        self.name = name
        self.wind_num = wind_num
        self.wind_list = wind_list
        self.wind_tick = wind_tick



def windgen(gen_dict):
    wind_size = 100
    segWind_list = list()
    for key,value in gen_dict.items():
        wind_num = len(value)//wind_size
        rem = len(value) % wind_size
        ext = wind_size - rem
        wind_list_1 = list()
        wind_tick_list = list()
        #wind_list_2 = list()
        start_pos = 0
        end_pos = len(value)
        #print(end_pos)
        i = 0
        k = 0
        if rem < ext: # extend
            pass
        else: # one more window
            wind_num += 1
        while i < wind_num:
            if i == (wind_num-1):
                wind_list_1.append((start_pos,end_pos))
                wind_tick_list.append(start_pos)
                wind_tick_list.append(end_pos)
            else:
                wind_list_1.append((start_pos, start_pos + wind_size))
                wind_tick_list.append(start_pos)
            i += 1
            start_pos += wind_size
        #start_pos1 = 0 + wind_size/2
        #while k < (wind_num -1):
        #    wind_list_2.append((start_pos1, start_pos1 + wind_size))
        #    #print((start_pos + wind_size/2, start_pos + 3*wind_size/2))
        #    k += 1
        #    start_pos1 += wind_size
        #joined_wind_list = wind_list_1 + wind_list_2
        #print(joined_wind_list)
        wind_list_1.sort()
        wind_tick_list.sort()
        #print(wind_tick_list)
        #print(wind_list_1)
        sw = WindAlign(key,wind_num,wind_list_1,wind_tick_list)
        segWind_list.append(sw)
    all_segname = [sw.name for sw in segWind_list]
    all_wind_num = [sw.wind_num for sw in segWind_list]
    #print(all_wind_num) 
    # wind_size = 200, [17, 13, 13, 12, 8, 7, 6, 5, 5, 4, 3]
    # wind_size = 100, [33, 27, 26, 24, 16, 14, 11, 11, 11, 8, 7]

    
    all_wind_list = [sw.wind_list for sw in segWind_list]
    print(all_wind_list) 
    # wind_size = 200, [[(0, 200), (200, 400), (400, 600), (600, 800), (800, 1000), (1000, 1200), (1200, 1400), (1400, 1600), (1600, 1800), (1800, 2000), (2000, 2200), (2200, 2400), (2400, 2600), (2600, 2800), (2800, 3000), (3000, 3200), (3200, 3302)], [(0, 200), (200, 400), (400, 600), (600, 800), (800, 1000), (1000, 1200), (1200, 1400), (1400, 1600), (1600, 1800), (1800, 2000), (2000, 2200), (2200, 2400), (2400, 2693)], [(0, 200), (200, 400), (400, 600), (600, 800), (800, 1000), (1000, 1200), (1200, 1400), (1400, 1600), (1600, 1800), (1800, 2000), (2000, 2200), (2200, 2400), (2400, 2591)], [(0, 200), (200, 400), (400, 600), (600, 800), (800, 1000), (1000, 1200), (1200, 1400), (1400, 1600), (1600, 1800), (1800, 2000), (2000, 2200), (2200, 2362)], [(0, 200), (200, 400), (400, 600), (600, 800), (800, 1000), (1000, 1200), (1200, 1400), (1400, 1610)], [(0, 200), (200, 400), (400, 600), (600, 800), (800, 1000), (1000, 1200), (1200, 1356)], [(0, 200), (200, 400), (400, 600), (600, 800), (800, 1000), (1000, 1105)], [(0, 200), (200, 400), (400, 600), (600, 800), (800, 1059)], [(0, 200), (200, 400), (400, 600), (600, 800), (800, 1063)], [(0, 200), (200, 400), (400, 600), (600, 751)], [(0, 200), (200, 400), (400, 667)]]
    # wind_size = 100, [[(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700), (700, 800), (800, 900), (900, 1000), (1000, 1100), (1100, 1200), (1200, 1300), (1300, 1400), (1400, 1500), (1500, 1600), (1600, 1700), (1700, 1800), (1800, 1900), (1900, 2000), (2000, 2100), (2100, 2200), (2200, 2300), (2300, 2400), (2400, 2500), (2500, 2600), (2600, 2700), (2700, 2800), (2800, 2900), (2900, 3000), (3000, 3100), (3100, 3200), (3200, 3302)], [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700), (700, 800), (800, 900), (900, 1000), (1000, 1100), (1100, 1200), (1200, 1300), (1300, 1400), (1400, 1500), (1500, 1600), (1600, 1700), (1700, 1800), (1800, 1900), (1900, 2000), (2000, 2100), (2100, 2200), (2200, 2300), (2300, 2400), (2400, 2500), (2500, 2600), (2600, 2693)], [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700), (700, 800), (800, 900), (900, 1000), (1000, 1100), (1100, 1200), (1200, 1300), (1300, 1400), (1400, 1500), (1500, 1600), (1600, 1700), (1700, 1800), (1800, 1900), (1900, 2000), (2000, 2100), (2100, 2200), (2200, 2300), (2300, 2400), (2400, 2500), (2500, 2591)], [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700), (700, 800), (800, 900), (900, 1000), (1000, 1100), (1100, 1200), (1200, 1300), (1300, 1400), (1400, 1500), (1500, 1600), (1600, 1700), (1700, 1800), (1800, 1900), (1900, 2000), (2000, 2100), (2100, 2200), (2200, 2300), (2300, 2362)], [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700), (700, 800), (800, 900), (900, 1000), (1000, 1100), (1100, 1200), (1200, 1300), (1300, 1400), (1400, 1500), (1500, 1610)], [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700), (700, 800), (800, 900), (900, 1000), (1000, 1100), (1100, 1200), (1200, 1300), (1300, 1356)], [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700), (700, 800), (800, 900), (900, 1000), (1000, 1105)], [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700), (700, 800), (800, 900), (900, 1000), (1000, 1059)], [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700), (700, 800), (800, 900), (900, 1000), (1000, 1063)], [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700), (700, 751)], [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 667)]]
    
    all_wind_tick = [sw.wind_tick for sw in segWind_list]
    #print(all_wind_tick)
    #wind_size = 200, [[0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3302], [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2693], [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2591], [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2362], [0, 200, 400, 600, 800, 1000, 1200, 1400, 1610], [0, 200, 400, 600, 800, 1000, 1200, 1356], [0, 200, 400, 600, 800, 1000, 1105], [0, 200, 400, 600, 800, 1059], [0, 200, 400, 600, 800, 1063], [0, 200, 400, 600, 751], [0, 200, 400, 667]]
    #wind_size = 100, [[0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3302], [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2693], [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2591], [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2362], [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1610], [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1356], [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1105], [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1059], [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1063], [0, 100, 200, 300, 400, 500, 600, 700, 751], [0, 100, 200, 300, 400, 500, 600, 667]]
    return segWind_list, all_segname

class interAlign: #all_0, all_1
    def __init__(self, name, raw_cd,norm_cd):#name, raw count dict, norm count dict
        self.name = name # (seq,"") or (seqa,seqb)
        self.raw_cd = raw_cd
        self.norm_cd = norm_cd

def readlen(SRI_list):
    for (seq,i,j) in SRI_list:
        if j-i > 200:
            print("long read", (seq,i,j))
        else:
            continue


def mapSRI(SRI_list,segWind_list):
    SRI_ali_list = list()
    for sw in segWind_list:
        print(sw.name)
        SRI_count_dict = dict()
        for (seq,i,j) in SRI_list:
            if seq != sw.name: continue
            for start,end in sw.wind_list:
                #print((start,end))
                if end >= i >= start and start <= j <= end:#1
                    cnt_number = SRI_count_dict.get((start,end),0)
                    cnt_number += 1
                    SRI_count_dict[(start,end)] = cnt_number
                    break
                elif end > i > start and j >end:#2
                    if end-i > (j-i)/2:
                        cnt_number = SRI_count_dict.get((start,end),0)
                        cnt_number += 1
                        SRI_count_dict[(start,end)] = cnt_number
                        break
                    elif end-i == (j-i)/2:
                        cnt_number = SRI_count_dict.get((start,end),0)
                        cnt_number += 1
                        SRI_count_dict[(start,end)] = cnt_number
                        break
                    else:
                        continue
                elif i < start and start < j < end:#3
                    if j - start > (j-i)/2:
                        cnt_number = SRI_count_dict.get((start,end),0)
                        cnt_number += 1
                        SRI_count_dict[(start,end)] = cnt_number
                        break
                    elif j - start == (j-i)/2:
                        continue
                    else:
                        continue
                else: #4
                    continue
        print(dict(sorted(SRI_count_dict.items())))
        sa = interAlign((sw.name,""),dict(sorted(SRI_count_dict.items())),0)
        SRI_ali_list.append(sa)
    return SRI_ali_list

def readlen_intra(intraLRI_listdr):
    for (seq,i,j,k,l) in intraLRI_listdr:
        if j-i > 200 or l-k >200:
            print("long read", (seq,i,j,k,l))
        else:
            continue

def mapIntra(intraLRI_listd,segWind_list):
    intra_ali_list = list()
    for sw in segWind_list:#each segment
        print(sw.name)
        intra_count_dict = dict()
        for (seq,i,j,k,l) in intraLRI_listd: #each read
            if seq != sw.name: continue
            chim_list = ()
            for start1,end1 in sw.wind_list: # alignment 1
                if end1 >= i >= start1 and start1 <= j <= end1:#1
                        chim_list += (start1,end1)
                elif end1 > i > start1 and j >end1:#2
                    if end1-i > (j-i)/2:
                        chim_list += (start1,end1)
                    elif end1-i == (j-i)/2:
                        chim_list += (start1,end1)
                    else:
                        continue
                elif i < start1 and start1 < j < end1:#3
                    if j - start1 > (j-i)/2:
                        chim_list += (start1,end1)
                    elif j - start1 == (j-i)/2:
                        continue
                    else:
                        continue
                else: #4
                    continue
                for start2,end2 in sw.wind_list:
                    if end2 >= k >= start2 and start2 <= l <= end2:#1
                        chim_list += (start2,end2)
                        #print(chim_list)
                        cnt_number = intra_count_dict.get(chim_list,0)
                        cnt_number += 1
                        intra_count_dict[chim_list] = cnt_number
                        break
                    elif end2 > k > start2 and l >end2:#2
                        if end2-k > (l-k)/2:
                            chim_list += (start2,end2)
                            #print(chim_list)
                            cnt_number = intra_count_dict.get(chim_list,0)
                            cnt_number += 1
                            intra_count_dict[chim_list] = cnt_number
                            break
                        elif end2-k == (l-k)/2:
                            chim_list += (start2,end2)
                            #print(chim_list)
                            cnt_number = intra_count_dict.get(chim_list,0)
                            cnt_number += 1
                            intra_count_dict[chim_list] = cnt_number
                            break
                        else:
                            continue
                    elif k < start2 and start2 < l < end2:#3
                        if l - start2 > (j-i)/2:
                            chim_list += (start2,end2)
                            #print(chim_list)
                            cnt_number = intra_count_dict.get(chim_list,0)
                            cnt_number += 1
                            intra_count_dict[chim_list] = cnt_number
                            break
                        elif l - start2 == (l-k)/2:
                            continue
                        else:
                            continue
                    else: #4
                        continue
                    #if len(chim_list) > 4:
                    #    print("error", (seq,i,j,k,l),(start1,end1),(start2,end2),chim_list)
        print(dict(sorted(intra_count_dict.items())))
        intra = interAlign((sw.name,sw.name),dict(sorted(intra_count_dict.items())),0)
        intra_ali_list.append(intra)
    return intra_ali_list

def readlen_inter(interLRI_list):
    for ((seqa,ai,aj),(seqb,bi,bj)) in interLRI_list:
        if aj-ai > 200 or bj-bi >200:
            print("long read", ((seqa,ai,aj),(seqb,bi,bj)))
        else:
            continue

def mapInter(interLRI_list,segWind_list):
    #s = 0
    interLRI_list.sort()
    inter_ali_list = list()
    for sw1 in segWind_list:
        for sw2 in segWind_list:
            if sw1.name >= sw2.name: continue
            #s += 1
            print(sw1.name,sw2.name)
            inter_count_dict = dict()
            for ((seqa,ai,aj),(seqb,bi,bj)) in interLRI_list:
                if seqa != sw1.name or seqb != sw2.name: continue
                for start1,end1 in sw1.wind_list:
                    chim_list = ()
                    if ai < start1 and aj > end1: continue #4
                    if start1 < aj and ai < end1:#1
                        chim_list += (start1,end1)
                    elif end1 > ai > start1 and aj >end1:#2
                        if end1-ai > (aj-ai)/2:
                            chim_list += (start1,end1)
                        elif end1-ai == (aj-ai)/2:
                            chim_list += (sw1.name, start1,end1)
                        else:
                            continue
                    elif ai < start1 and start1 < aj < end1:#3
                        if aj - start1 > (aj-ai)/2:
                            chim_list += (start1,end1)
                        elif aj - start1 == (aj-ai)/2:
                            continue
                        else:
                            continue
                    else: #4
                        continue
                    for start2,end2 in sw2.wind_list:
                        if end2 >= bi >= start2 and start2 <= bj <= end2:#1
                            chim_list += (start2,end2)
                            #print(chim_list)
                            cnt_number = inter_count_dict.get(chim_list,0)
                            cnt_number += 1
                            inter_count_dict[chim_list] = cnt_number
                            break
                        elif end2 > bi > start2 and bj >end2:#2
                            if end2-bi > (bj-bi)/2:
                                chim_list += (start2,end2)
                                #print(chim_list)
                                cnt_number = inter_count_dict.get(chim_list,0)
                                cnt_number += 1
                                inter_count_dict[chim_list] = cnt_number
                                break
                            elif end2-bi == (bj-bi)/2:
                                chim_list += (start2,end2)
                                #print(chim_list)
                                cnt_number = inter_count_dict.get(chim_list,0)
                                cnt_number += 1
                                inter_count_dict[chim_list] = cnt_number
                                break
                            else:
                                continue
                        elif bi < start2 and start2 < bj < end2:#3
                            if bj - start2 > (bj-bi)/2:
                                chim_list += (start2,end2)
                                #print(chim_list)
                                cnt_number = inter_count_dict.get(chim_list,0)
                                cnt_number += 1
                                inter_count_dict[chim_list] = cnt_number
                                break
                            elif bj - start2 == (bj-bi)/2:
                                continue
                            else:
                                continue
                        else: #4
                            continue
            print(dict(sorted(inter_count_dict.items())))
            inter = interAlign((sw1.name,sw2.name),dict(sorted(inter_count_dict.items())),0)
            inter_ali_list.append(inter)
    #print(s), 66 combinations
    return inter_ali_list

class SegInter: #pair_list
    def __init__(self, name,Mg, Ag, weight,raw1,raw0,norm1,norm0,ppm1,ppm0,fc):#norm1,norm0 read count added later
        self.name = name #seg (pair) name, window (combination) 
        self.Mg = Mg
        self.Ag = Ag
        self.weight = weight
        self.raw1 = raw1
        self.raw0 = raw0
        self.norm1 = norm1
        self.norm0 = norm0
        self.ppm1 = ppm1
        self.ppm0 = ppm0
        self.fc = fc

def normf(all_1,all_0,totreadbp_list,Agv,Mgv):
    ## part 1: calculate Mg and Ag
    pair_list = list()
    trun_list = list()
    totr1 = totreadbp_list[0]
    totr0 = totreadbp_list[1]
    print("totr1",totr1)
    print("totr0",totr0)
    for it1,it0 in zip(all_1,all_0):
        #print(it1)
        #print(it1.name)
        #print(it1.raw_cd)
        for k1,v1 in it1.raw_cd.items():
            for k0,v0 in it0.raw_cd.items():
            #print(k1,k0)
                if k1 != k0: continue
                    #print("error")
                    #print(k1,v1,k0,v0,it0.name)
                Mg = math.log((v1/totr1)/(v0/totr0),2)
                #print("Mg",Mg)
                Ag = 0.5*math.log((v1/totr1)*(v0/totr0),2)
                #print("Ag",Ag)
                weight = (totr1-v1)/(totr1*v1) + (totr0 - v0)/(totr0*v0)
                #print("weight",weight)
                #print(it1.name,k1,Mg,Ag,weight)
                ppm1 = 1000000*v1/totr1
                ppm0 = 1000000*v0/totr0
                fc = math.log(ppm1/ppm0,2)
                pair = SegInter((it1.name,k1),Mg,Ag,weight,v1,v0,0,0,ppm1,ppm0,fc)
                pair_list.append(pair)
    # part 2: truncate Mg and Ag value to get G* set
    all_Mg = [pair.Mg for pair in pair_list]
    all_Ag = [pair.Ag for pair in pair_list]
    Mg_uplim = np.quantile(all_Mg,1-Mgv)
    Mg_lowlim = np.quantile(all_Mg,Mgv)
    Ag_uplim = np.quantile(all_Ag,1-Agv)
    Ag_lowlim = np.quantile(all_Ag,Agv)
    for pr in pair_list:
        if Mg_uplim >= pr.Mg >= Mg_lowlim and Ag_uplim >= pr.Ag >= Ag_lowlim:
            trun_list.append(pr)
        else:
            continue
    # part 3: calculate normalisation factor
    w_sum = 0
    WMg_sum = 0
    for tr in trun_list:
        w_sum += tr.weight
        WMg_sum += tr.weight * tr.Mg
    TMM = 2**(WMg_sum/w_sum)
    print(TMM)
    all_pr1_raw = [pr.raw1 for pr in pair_list]
    all_pr0_raw = [pr.raw0 for pr in pair_list]
    all_pr1_inter_LRI_raw = [pr.raw1 for pr in pair_list if pr.name[0][1] and pr.name[0][0] != pr.name[0][1]]
    all_pr0_inter_LRI_raw = [pr.raw0 for pr in pair_list if pr.name[0][1] and pr.name[0][0] != pr.name[0][1]]
    all_pr1_intra_LRI_raw = [pr.raw1 for pr in pair_list if pr.name[0][1] and pr.name[0][0] == pr.name[0][1]]
    all_pr0_intra_LRI_raw = [pr.raw0 for pr in pair_list if pr.name[0][1] and pr.name[0][0] == pr.name[0][1]]
    all_pr0_SRI_raw = [pr.raw0 for pr in pair_list if not pr.name[0][1]]
    all_pr1_SRI_raw = [pr.raw1 for pr in pair_list if not pr.name[0][1]]
    all_pr1_ppm = [pr.ppm1 for pr in pair_list]
    all_pr0_ppm = [pr.ppm0 for pr in pair_list]
    all_pr1_inter_LRI_ppm = [pr.ppm1 for pr in pair_list if pr.name[0][1] and pr.name[0][0] != pr.name[0][1]]
    all_pr0_inter_LRI_ppm = [pr.ppm0 for pr in pair_list if pr.name[0][1] and pr.name[0][0] != pr.name[0][1]]
    all_pr1_intra_LRI_ppm = [pr.ppm1 for pr in pair_list if pr.name[0][1] and pr.name[0][0] == pr.name[0][1]]
    all_pr0_intra_LRI_ppm = [pr.ppm0 for pr in pair_list if pr.name[0][1] and pr.name[0][0] == pr.name[0][1]]
    all_pr0_SRI_ppm = [pr.ppm0 for pr in pair_list if not pr.name[0][1]]
    all_pr1_SRI_ppm = [pr.ppm1 for pr in pair_list if not pr.name[0][1]]
    return TMM, pair_list,trun_list,all_pr0_SRI_raw,all_pr1_SRI_raw,all_pr1_inter_LRI_raw,all_pr0_inter_LRI_raw,all_pr1_intra_LRI_raw,all_pr0_intra_LRI_raw,all_pr1_raw,all_pr0_raw,all_pr1_ppm,all_pr0_ppm,all_pr1_inter_LRI_ppm,all_pr0_inter_LRI_ppm,all_pr1_intra_LRI_ppm,all_pr0_intra_LRI_ppm,all_pr0_SRI_ppm,all_pr1_SRI_ppm

def readnorm(TMM, pair_list,all_1,all_0):
    for pr in pair_list:
        pr.norm1 = pr.raw1* math.sqrt(TMM)
        pr.norm0 = pr.raw0/math.sqrt(TMM)
        pr.fc = math.log(pr.norm1/pr.norm0,2)
        all_pr1_norm = [pr.norm1 for pr in pair_list]
        all_pr0_norm = [pr.norm0 for pr in pair_list]

        all_pr1_LRI_norm = [pr.norm1 for pr in pair_list if pr.name[0][1]]
        all_pr0_LRI_norm = [pr.norm0 for pr in pair_list if pr.name[0][1]]

        all_pr1_inter_LRI_norm = [pr.norm1 for pr in pair_list if pr.name[0][1] and pr.name[0][0] != pr.name[0][1]]
        all_pr0_inter_LRI_norm = [pr.norm0 for pr in pair_list if pr.name[0][1] and pr.name[0][0] != pr.name[0][1]]

        all_pr1_intra_LRI_norm = [pr.norm1 for pr in pair_list if pr.name[0][1] and pr.name[0][0] == pr.name[0][1]]
        all_pr0_intra_LRI_norm = [pr.norm0 for pr in pair_list if pr.name[0][1] and pr.name[0][0] == pr.name[0][1]]

        all_pr0_SRI_norm = [pr.norm0 for pr in pair_list if not pr.name[0][1]]
        all_pr1_SRI_norm = [pr.norm1 for pr in pair_list if not pr.name[0][1]]

    for it1,it0 in zip(all_1,all_0):
        normd_1 = dict()
        normd_0 = dict()
        for pr in pair_list:
            if pr.name[0] != it1.name: continue
            for k1 in it1.raw_cd.keys():
                if pr.name[1] != k1: continue
                normd_1[k1] = pr.raw1* math.sqrt(TMM)
                normd_0[k1] = pr.raw0/math.sqrt(TMM)
        it1.norm_cd = normd_1
        it0.norm_cd = normd_0
    return pair_list, all_1,all_0, all_pr1_norm,all_pr0_norm,all_pr1_SRI_norm,all_pr1_LRI_norm,all_pr0_SRI_norm,all_pr0_LRI_norm, all_pr1_inter_LRI_norm,all_pr0_inter_LRI_norm,all_pr1_intra_LRI_norm,all_pr0_intra_LRI_norm
    #for it1 in all_1:
    #    norm_cd = dict()
    #    for k1,v1 in it1.raw_cd.items():
    #        norm_cd[k1] = v1 * math.sqrt(TMM)
    #    it1.norm_cd = norm_cd
    #all_it1_normcd = [it1.norm_cd for it1 in all_1]
    #for it0 in all_0:
    #    norm_cd = dict()
    #    for k0,v0 in it0.raw_cd.items():
    #        norm_cd[k0] = v0/math.sqrt(TMM)
    #        #sx0.normrd = sx0.read_count/math.sqrt(TMM)
    #    it0.norm_cd = norm_cd
    #all_it0_normcd = [it0.norm_cd for it0 in all_0]
    #return all_it1_normcd, all_it0_normcd


def summaryTSV(pair_list,samp_name,whole_wind_dir):
    with open(f"{whole_wind_dir}/{samp_name}_Window_summaryTable.tsv", "w",newline='') as outtable:
        tsv_writer = csv.writer(outtable,delimiter="\t")
        tsv_writer.writerow(["name","'-NSP2 raw","'+NSP2 raw","'-NSP2 normalised", "'+NSP2 normalised", "'-NSP2 rpm", "+NSP2 rpm", "fold change"])
        for pr in pair_list:
            tsv_writer.writerow([pr.name, pr.raw0, pr.raw1, pr.norm0, pr.norm1, pr.ppm0, pr.ppm1,pr.fc])

#pair = SegInter((it1.name,k1),Mg,Ag,weight,v1,v0,0,0,ppm1,ppm0)
def trunTSV(trun_list,samp_name,whole_wind_dir):
    with open(f"{whole_wind_dir}/trun_{samp_name}_summaryTable.tsv", "w",newline='') as outtable:
        tsv_writer = csv.writer(outtable,delimiter="\t")
        tsv_writer.writerow(["name","'-NSP2 raw","'+NSP2 raw","'-NSP2 rpm", "'+NSP2 rpm", "fold change"])
        for tr in trun_list:
            tr.fc = math.log(tr.ppm1/tr.ppm0,2)
            tsv_writer.writerow([tr.name, tr.raw0,tr.raw1,tr.ppm0, tr.ppm1,tr.fc])



def plotLog(all_pr1_SRI_norm,all_pr1_LRI_norm,all_pr0_SRI_norm,all_pr0_LRI_norm,all_pr1_norm,all_pr0_norm,samp_name,samp_wind_dir,Agv,Mgv,AgMg_dir,all_pr1_inter_LRI_norm,all_pr0_inter_LRI_norm,all_pr1_intra_LRI_norm,all_pr0_intra_LRI_norm):
    s = 1.0
    pz = plt.figure(figsize=(12*s,10*s))
    plt.xscale("log", base = 10)
    plt.yscale("log",base = 10)
    plt.scatter(all_pr1_inter_LRI_norm,all_pr0_inter_LRI_norm, color="blue", label = "inter LRI")
    plt.scatter(all_pr1_intra_LRI_norm,all_pr0_intra_LRI_norm, color="green", label = "intra LRI")
    plt.scatter(all_pr1_SRI_norm,all_pr0_SRI_norm, color="red",edgecolors="black", label = "SRI")
    plt.legend(loc = "upper left",markerscale = 6)
    #ax.axline((0, 0), slope=1)
    if samp_name == "DLP":
        plt.xlabel("Normalised Read Count With NSP2")
        plt.ylabel("Normalised Read Count Without NSP2")
    elif samp_name == "In_Vitro":
        plt.xlabel("Normalised Read Count With NSP2")
        plt.ylabel("Normalised Read Count Without NSP2")
    elif samp_name == "InVitro_VS_DLP_withNSP2":
        plt.xlabel("Normalised read DLP")
        plt.ylabel("Normalised read in vitro")
    else: #samp_name == ""InVitro_VS_DLP_withoutNSP2""
        plt.xlabel("Normalised read DLP")
        plt.ylabel("Normalised read in vitro")
    plt.xlim([1,10**7])
    plt.ylim([1,10**7])
    #sc.SymmetricalLogScale("Normalised Read Count With NSP2",base=10)
    #sc.SymmetricalLogScale("Normalised Read Count Without NSP2",base=10)
    #m,b = np.polyfit(all_pr1_norm,all_pr0_norm, 1)
    #plt.plot()
    #plt.axline(xy1=(0,b),slope = m)
    identity_line = np.linspace(max(min(all_pr1_norm), min(all_pr0_norm)), min(max(all_pr1_norm), max(all_pr0_norm)))
    plt.plot(identity_line, identity_line, color="black", linestyle="dashed", linewidth=3.0)
    #x = range(0,1,0.5)
    #y = range(0,1,0.5)
    #plt.plot(x,y)
    plt.title(f"{samp_name} Ag={Agv} Mg={Mgv}")
    plt.subplots_adjust(bottom=0.2, top=1.2)
    plt.legend()
    plt.show()
    pz.savefig(f"{AgMg_dir}/{samp_name}_window_inter_NSP2.svg", bbox_inches = 'tight', pad_inches = 0.1*s)
    pz.savefig(f"{AgMg_dir}/{samp_name}_window_inter_NSP2.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)


def plotLogRaw(all_pr0_SRI_raw,all_pr1_SRI_raw,all_pr1_inter_LRI_raw,all_pr0_inter_LRI_raw,all_pr1_intra_LRI_raw,all_pr0_intra_LRI_raw,all_pr1_raw,all_pr0_raw,samp_name,samp_wind_dir,Agv,Mgv,AgMg_rawdir):
    s = 1.0
    pz = plt.figure(figsize=(12*s,10*s))
    plt.xscale("log", base = 10)
    plt.yscale("log",base = 10)
    plt.scatter(all_pr1_inter_LRI_raw,all_pr0_inter_LRI_raw, color="blue", label = "inter LRI")
    plt.scatter(all_pr1_intra_LRI_raw,all_pr0_intra_LRI_raw, color="green", label = "intra LRI")
    plt.scatter(all_pr1_SRI_raw,all_pr0_SRI_raw, color="red",edgecolors="black", label = "SRI")
    plt.legend(loc = "upper left")
    #ax.axline((0, 0), slope=1)
    if samp_name == "DLP":
        plt.xlabel("Raw Read Count With NSP2")
        plt.ylabel("Raw Read Count Without NSP2")
    elif samp_name == "In_Vitro":
        plt.xlabel("Raw Read Count With NSP2")
        plt.ylabel("Raw Read Count Without NSP2")
    elif samp_name == "InVitro_VS_DLP_withNSP2":
        plt.xlabel("Raw read DLP")
        plt.ylabel("Raw read in vitro")
    else: #samp_name == ""InVitro_VS_DLP_withoutNSP2""
        plt.xlabel("Raw read DLP")
        plt.ylabel("Raw read in vitro")
    plt.xlim([1,10**7])
    plt.ylim([1,10**7])
    #sc.SymmetricalLogScale("Normalised Read Count With NSP2",base=10)
    #sc.SymmetricalLogScale("Normalised Read Count Without NSP2",base=10)
    #m,b = np.polyfit(all_pr1_norm,all_pr0_norm, 1)
    #plt.plot()
    #plt.axline(xy1=(0,b),slope = m)
    identity_line = np.linspace(max(min(all_pr1_raw), min(all_pr0_raw)), min(max(all_pr1_raw), max(all_pr0_raw)))
    plt.plot(identity_line, identity_line, color="black", linestyle="dashed", linewidth=3.0)
    #x = range(0,1,0.5)
    #y = range(0,1,0.5)
    #plt.plot(x,y)
    plt.title(f"{samp_name} Ag={Agv} Mg={Mgv}")
    plt.subplots_adjust(bottom=0.2, top=1.2)
    plt.show()
    pz.savefig(f"{AgMg_rawdir}/{samp_name}_window_inter_NSP2.svg", bbox_inches = 'tight', pad_inches = 0.1*s)
    pz.savefig(f"{AgMg_rawdir}/{samp_name}_window_inter_NSP2.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)

def plotLogPpm(all_pr1_ppm,all_pr0_ppm,all_pr1_inter_LRI_ppm,all_pr0_inter_LRI_ppm,all_pr1_intra_LRI_ppm,all_pr0_intra_LRI_ppm,all_pr0_SRI_ppm,all_pr1_SRI_ppm,samp_name,samp_wind_dir,ppmdir):
    s = 1.0
    pz = plt.figure(figsize=(12*s,10*s))
    plt.xscale("log", base = 10)
    plt.yscale("log",base = 10)
    
    plt.scatter(all_pr1_inter_LRI_ppm,all_pr0_inter_LRI_ppm, color="blue",label = "Inter LRI")
    #inter_a, inter_b = np.polyfit(math.log(all_pr1_inter_LRI_ppm),math.log(all_pr0_inter_LRI_ppm), 1) #produce gradient and y-int of given data
    ##cor_line = np.linspace(10**(-2),10**6,100) #generate x coordinates for best-fit line
    #cor_line=np.arange(0,math.ceil(max(math.log10(-all_pr1_inter_LRI_ppm),math.log10(all_pr0_inter_LRI_ppm))),step=1)
    #plt.plot(cor_line, inter_a*cor_line + inter_b, color="blue", linewidth=0.8) #plotting the best fit line


    plt.scatter(all_pr1_intra_LRI_ppm,all_pr0_intra_LRI_ppm, color="green",label = "intra LRI")
    #intra_a, intra_b = np.polyfit(math.log(all_pr1_intra_LRI_ppm),math.log(all_pr0_intra_LRI_ppm), 1) #produce gradient and y-int of given data
    ##cor_line = np.linspace(10**(-2),10**6,100) #generate x coordinates for best-fit line
    #cor_line=np.arange(0,math.ceil(max(math.log10(-all_pr1_intra_LRI_ppm),math.log10(all_pr0_intra_LRI_ppm))),step=1)
    #plt.plot(cor_line, intra_a*cor_line + intra_b, color="green", linewidth=0.8) #plotting the best fit line


    plt.scatter(all_pr1_SRI_ppm,all_pr0_SRI_ppm, color="red",edgecolors="black",label = "SRI")
    #sri_a, sri_b = np.polyfit(math.log(all_pr1_SRI_ppm), math.log(all_pr0_SRI_ppm), 1) #produce gradient and y-int of given data
    ##cor_line = np.linspace(10**(-2),10**6,100) #generate x coordinates for best-fit line
    #cor_line=np.arange(0,math.ceil(max(math.log10(-all_pr1_SRI_ppm),math.log10(all_pr0_SRI_ppm))),step=1)
    #plt.plot(cor_line, sri_a*cor_line + sri_b, color="red", linewidth=0.8) #plotting the best fit line

    plt.legend(loc = "upper left",fontsize = 15, markerscale = 6)
    #ax.axline((0, 0), slope=1)
    if samp_name == "DLP":
        plt.xlabel("Read per Million With NSP2",fontsize = 15)
        plt.ylabel("Read per Million Without NSP2",fontsize = 15)
    elif samp_name == "In_Vitro":
        plt.xlabel("Read per Million With NSP2",fontsize = 15)
        plt.ylabel("Read per Million Without NSP2",fontsize = 15)
    elif samp_name == "InVitro_VS_DLP_withNSP2":
        plt.xlabel("Read per Million DLP",fontsize = 15)
        plt.ylabel("Read per Million in vitro",fontsize = 15)
    else: #samp_name == ""InVitro_VS_DLP_withoutNSP2""
        plt.xlabel("Read per Million DLP",fontsize = 15)
        plt.ylabel("Read per Million in vitro",fontsize = 15)
    plt.xlim([10**(-2),10**6])
    plt.ylim([10**(-2),10**6])
    #sc.SymmetricalLogScale("Normalised Read Count With NSP2",base=10)
    #sc.SymmetricalLogScale("Normalised Read Count Without NSP2",base=10)
    #m,b = np.polyfit(all_pr1_norm,all_pr0_norm, 1)
    #plt.plot()
    #plt.axline(xy1=(0,b),slope = m)
    identity_line = np.linspace(max(min(all_pr1_ppm), min(all_pr0_ppm)), min(max(all_pr1_ppm), max(all_pr0_ppm)))
    plt.plot(identity_line, identity_line, color="black", linestyle="dashed", linewidth=3.0)
    #x = range(0,1,0.5)
    #y = range(0,1,0.5)
    #plt.plot(x,y)
    plt.title(f"{samp_name}",fontsize = 20)
    plt.subplots_adjust(bottom=0.2, top=1.2)
    plt.tick_params(axis="both", which="major", labelsize=8)
    plt.show()
    pz.savefig(f"{ppmdir}/{samp_name}_window_inter_NSP2.svg", bbox_inches = 'tight', pad_inches = 0.1*s)
    pz.savefig(f"{ppmdir}/{samp_name}_window_inter_NSP2.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)

    #save bestfit equations
    #with open(f"{ppmdir}/{samp_name}_bestfitline.tsv", "w",newline='') as outtable:
    #    tsv_writer = csv.writer(outtable,delimiter="\t")
    #    tsv_writer.writerow(["inter_grad","inter-yint","intra_grad","intra-yint","SRI_grad","SRI_yint"])
    #    tsv_writer.writerow([inter_a, inter_b,intra_a, intra_b,sri_a, sri_b])


def interHeatmap1(segWind_list,all_1,exp_dir1):
    mat_dict = dict()
    for sw1 in segWind_list:
        for sw2 in segWind_list:
            if sw1.name > sw2.name: continue
            #generate bot intra and inter
            #print(sw1.name,sw2.name)
            print(sw2.name, sw2.wind_num,sw1.name, sw1.wind_num)
            smat = np.zeros(shape = (sw1.wind_num,sw2.wind_num))
            # generate empty matrix of size (x,y), x = number of windows in segA, y = in segB
            #print(smat)
            #exit()
            #pr.name = (seq comb, (s1,e1,s2,e2))
            for it1 in all_1:
                if it1.name != (sw1.name,sw2.name):continue
                #print(f"it1.name {it1.name}")
                for k1,v1 in it1.raw_cd.items():
                    i = sw1.wind_list.index(k1[:2])
                    sw1s = int(sw1.wind_list[i][1])-int(sw1.wind_list[i][0])
                    #print(i)
                    #print("sw1s",sw1s)
                    j = sw2.wind_list.index(k1[2:])
                    sw2s = int(sw2.wind_list[j][1])-int(sw2.wind_list[j][0])
                    #print(j)
                    #print("sw2s",sw2s)
                    #print(f"sw1 = {sw1.name} sw2={sw2.name} sw1s={sw1s} sw2s={sw2s} sw1s*sw2s:{sw1s*sw2s}")
                    if (sw1.name == sw2.name and i >= j) or sw1.name != sw2.name:
                        smat[i,j] += v1/(sw1s*sw2s)
                    else:
                        smat[j,i] += v1/(sw1s*sw2s)
                #print(smat)
                smat[smat == 0] = np.nan
                mat_dict[(sw1.name,sw2.name)] = smat
                s = 1.0
                pz = plt.figure(figsize=(12*s,10*s))
                plt.imshow(smat,cmap = 'viridis')
                plt.colorbar()
                plt.title(f"{sw1.name[5:]} {sw2.name[5:]}",fontsize = 36)
                plt.xlabel(f"{sw2.name[5:]}", fontsize = 25)
                plt.ylabel(f"{sw1.name[5:]}", fontsize = 25)
                xtick_list = list(np.arange(0,sw2.wind_num,step=1)) 
                ytick_list = list(np.arange(0,sw1.wind_num,step=1)) 
                plt.xticks(xtick_list,labels=sw2.wind_tick[1:],rotation=45, fontsize = 13)
                plt.yticks(ytick_list,labels=sw1.wind_tick[1:], fontsize = 13)
                plt.show()
                pz.savefig(f"{exp_dir1}/norm_{sw1.name}_{sw2.name}_heatmap.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)
                pz.savefig(f"{exp_dir1}/norm_{sw1.name}_{sw2.name}_heatmap.svg", bbox_inches = 'tight', pad_inches = 0.1*s)

def polyA_hist(segWind_list, gen_dict,whole_wind_dir):
    CTcountHeat_dir = os.path.join(whole_wind_dir,f"CT_heatmap")
    if not os.path.isdir(CTcountHeat_dir):
        os.mkdir(CTcountHeat_dir)
    mat_dict = dict()
    for sw1 in segWind_list:
        for sw2 in segWind_list:
            if sw1.name < sw2.name: continue
            print(sw1.name, sw1.wind_num,sw2.name, sw2.wind_num)
            smat = np.zeros(shape = (sw1.wind_num,sw2.wind_num))
            for i, (istart, iend) in enumerate(sw1.wind_list):
                Anum1 = len(re.findall("T", gen_dict[sw1.name][istart:iend])) + len(re.findall("C", gen_dict[sw1.name][istart:iend]))
                sw1s = int(iend) - int(istart)
                print(f"{sw1.name} {istart} {iend}, i={i}, Anum1={Anum1}, sw1s={sw1s}")
                for j, (jstart,jend) in enumerate(sw2.wind_list):
                    Anum2 = len(re.findall("T", gen_dict[sw2.name][jstart:jend])) + len(re.findall("C", gen_dict[sw2.name][jstart:jend]))
                    sw2s = int(jend) -int(jstart)
                    print(f"{sw2.name} {jstart} {jend}, j={j}, Anum2={Anum2}, sw2s={sw2s}")
                    if (sw1.name == sw2.name and i >= j) or sw1.name != sw2.name:
                        smat[i,j] = Anum1*Anum2/ (sw1s*sw2s)
                    else:
                        smat[j,i] += Anum1*Anum2/ (sw1s*sw2s)
            smat[smat == 0] = np.nan
            mat_dict[(sw1.name,sw2.name)] = smat
            s = 1.0
            pz = plt.figure(figsize=(12*s,10*s))
            plt.imshow(smat,cmap = 'viridis')
            plt.colorbar()
            plt.title(f"T Count {sw1.name[5:]} {sw2.name[5:]}",fontsize = 36)
            plt.xlabel(f"{sw2.name[5:]}", fontsize = 25)
            plt.ylabel(f"{sw1.name[5:]}", fontsize = 25)
            xtick_list = list(np.arange(0,sw2.wind_num,step=1)) 
            ytick_list = list(np.arange(0,sw1.wind_num,step=1)) 
            plt.xticks(xtick_list,labels=sw2.wind_tick[1:],rotation=45, fontsize = 13)
            plt.yticks(ytick_list,labels=sw1.wind_tick[1:], fontsize = 13)
            plt.show()
            pz.savefig(f"{CTcountHeat_dir}/AT_{sw1.name}_{sw2.name}_heatmap.pdf", bbox_inches = 'tight', pad_inches = 0.1*s)
            pz.savefig(f"{CTcountHeat_dir}/AT_{sw1.name}_{sw2.name}_heatmap.svg", bbox_inches = 'tight', pad_inches = 0.1*s)

                





def testcircosTSV(pair_list,exp_name1, exp_name0,circos_dir):
    ppmcut = 350
    fname1w = f"{exp_name1}_rpmcut{ppmcut}_withIntra"
    fname0w = f"{exp_name0}_rpmcut{ppmcut}_withIntra"
    fname1 = f"{exp_name1}_rpmcut{ppmcut}"
    fname0 = f"{exp_name0}_rpmcut{ppmcut}"
    ############################
    # inter & intra
    #############################
    #with open(f"{cir_exp_dir1}/{fname1w}_circos.tsv", "w",newline='') as outtable1:
    #    tsv_writer = csv.writer(outtable1,delimiter="\t")
    #    for pr in pair_list:
    #        if pr.ppm1 < ppmcut or not pr.name[0][1]: continue
    #        # smaller than cut-off, intra, SRI
    #        #S3 88 125 S4 432 469 color=col5
    #        #print(pr.name) #(('SA11_S1', 'SA11_S3'), (600, 800, 1000, 1200))
    #        #print(pr.name,pr.ppm1)
    #        tsv_writer.writerow([pr.name[0][0][5:], pr.name[1][0], pr.name[1][1], pr.name[0][1], pr.name[1][2], pr.name[1][3]])
    #with open(f"{cir_exp_dir0}/{fname0w}_circos.tsv", "w",newline='') as outtable0:
    #    tsv_writer = csv.writer(outtable0,delimiter="\t")
    #    for pr in pair_list:
    #        #print(pr.name[0][0][5:])
    #        if pr.ppm0 < ppmcut or not pr.name[0][1]:
    #            print("no",pr.name,pr.ppm0)
    #        else:
    #            print("Yes",pr.name,pr.ppm0)
    #        # smaller than cut-off, intra, SRI
    #        #S3 88 125 S4 432 469 color=col5
    #        #print(pr.name) #(('SA11_S1', 'SA11_S3'), (600, 800, 1000, 1200))
    #            tsv_writer.writerow([pr.name[0][0][5:], pr.name[1][0], pr.name[1][1], pr.name[0][1], pr.name[1][2], pr.name[1][3]])
    #########################
    # inter only
    ###########################
    with open(f"{circos_dir}/{fname1}_circos.tsv", "w",newline='') as outtable1:
        tsv_writer = csv.writer(outtable1,delimiter="\t")
        for pr in pair_list:
            if pr.ppm1 < ppmcut or pr.name[0][1] == pr.name[0][0] or not pr.name[0][1]: continue
            # smaller than cut-off, intra, SRI
            #S3 88 125 S4 432 469 color=col5
            #print(pr.name) #(('SA11_S1', 'SA11_S3'), (600, 800, 1000, 1200))
            #print(pr.name,pr.ppm1)
            tsv_writer.writerow([pr.name[0][0][5:], pr.name[1][0], pr.name[1][1], pr.name[0][1][5:], pr.name[1][2], pr.name[1][3]])
    with open(f"{circos_dir}/{fname0}_circos.tsv", "w",newline='') as outtable0:
        tsv_writer = csv.writer(outtable0,delimiter="\t")
        for pr in pair_list:
            #print(pr.name[0][0][5:])
            if pr.ppm0 < ppmcut or pr.name[0][1] == pr.name[0][0] or not pr.name[0][1]: continue
                #print("no",pr.name,pr.ppm0)
                #print("Yes",pr.name,pr.ppm0)
            # smaller than cut-off, intra, SRI
            #S3 88 125 S4 432 469 color=col5
            #print(pr.name) #(('SA11_S1', 'SA11_S3'), (600, 800, 1000, 1200))
            tsv_writer.writerow([pr.name[0][0][5:], pr.name[1][0], pr.name[1][1], pr.name[0][1][5:], pr.name[1][2], pr.name[1][3]])
    return fname1,fname0,fname1w,fname0w

def circosConfig(fname, circos_dir):
    # need: circos_links, circos_histo
    ## tick config
    ticks = f"""
        radius = 1r
        color = black
        thickness = 2p
        multiplier = 1
        format = %d
        <tick>
            spacing = 1u
            size = 10p
        </tick>
        <tick>
            spacing = 10u
            size = 15p
            show_label = yes
            label_size = 20p
            label_offset = 7p
            format = %d
            label_parallel = yes
        </tick>"""
    ## ideogram config
    ideogram = f"""
        <spacing>
            default = 0.015r
        </spacing>
        radius = 0.80r
        thickness = 40p
        fill = yes
        stroke_color = dgrey
        stroke_thickness = 2p
        show_label = yes
        show_bands = yes
        fill_bands = yes
        band_transparency = 0
        band_stroke_color = black
        label_font = default
        label_radius = dims(image,radius) - 60p;
        label_size = 60
        label_parallel = yes"""
    #    label_size = 30
    ## color config
    colors = f"""
        seg1 = 204,222,235
        seg2 = 190,210,226
        seg3 = 176,198,217
        seg4 = 158,176,208
        seg5 = 149,173,199
        seg6 = 158,181,205
        seg7 = 149,173,199
        seg8 = 217,217,217
        seg9 = 222,238,247
        seg10 = 217,234,244
        seg11 = 213,230,241"""
    ## link config
    links = f"""
        <link>
            file = {fname}_circos.tsv
            color = vlgrey
            radius = 0.99r
            bezier_radius = 0.1r
            thickness = 6
        </link>"""
    #        thickness = 4
    #        radius = 0.88r
    ## plot config
    #    <plot>
    #        type = text
    #        color = red
    #        file = circos_bands.tsv
    #        r0 = 1.008r
    #        r1 = 1.3r
    #        label_size = 24
    #        label_font = condensed
    #    </plot>
    if False:
        plots = f"""
            <plot>
                type = histogram
                file = circos_histo.tsv
                r1 = 0.98r
                r0 = 0.88r
                max = 20
                min = 0
                stroke_type = outline
                thickness = 1
                <backgrounds> <background>
                    color = vvlgrey
                </background> </backgrounds>
                <axes> <axis>
                    spacing = 0.2r
                    color = lgrey
                    thickness = 1
                </axis> </axes>
            </plot>"""
    ## main config file
    config = f"""# main circos config file
        show_ticks = yes
        show_tick_labels = yes
        
        <ticks> {ticks}
        </ticks>
        <ideogram> {ideogram}
        </ideogram>
        <colors> {colors}
        </colors>
        
        karyotype = circos_karyo.tsv
        chromosomes_units = 100
        
        <links> {links}
        </links>
        
        # basic settings files
        <image> 
            <<include etc/image.conf>>
            radius* = 600p
        </image>
        <<include etc/colors_fonts_patterns.conf>> 
        <<include etc/housekeeping.conf>> """
    #    <plots> {plots}
    #    </plots>
    ## write circos config file
    circos_conf = os.path.join(f"{circos_dir}/{fname}_circos_config.conf")
    with open(circos_conf, "w") as circonf:
        circonf.write(config)

def doCircos(fname, circos_dir):
    ## start circos and create plots
    circos_conf = os.path.abspath(os.path.join(f"{circos_dir}/{fname}_circos_config.conf"))
    wd = os.getcwd()
    os.chdir(circos_dir)
    C_call = ["circos", "-conf", circos_conf]
    subprocess.call(C_call, shell=False)
    os.rename("circos.png", f"{fname}_circos_plot.png")
    os.rename("circos.svg", f"{fname}_circos_plot.svg")
    os.chdir(wd)



# first gener
#def doCircos(fname, opt):
#    ## start circos and create plots
#    circos_conf = os.path.abspath(os.path.join(opt.pfx, f"{fname}_circos_config.conf"))
#    wd = os.getcwd()
#    os.chdir(opt.pfx)
#    C_call = [opt.cip, "-conf", circos_conf]
#    call(C_call, shell=False)
#    os.rename("circos.png", f"{fname}_circos_plot.png")
#    os.rename("circos.svg", f"{fname}_circos_plot.svg")
#    os.chdir(wd)



def saveData(pcl_data, fname): ## save data with pickle
    with open(f"{fname}.pcl", "w+b") as pcl_out:
        pickle.dump(pcl_data, pcl_out , protocol=4)


if __name__ == "__main__":
    genomefasta = sys.argv[1]
    pcl1 = sys.argv[2]
    pcl0 = sys.argv[3]
    samp_name = sys.argv[4]
    Agv = sys.argv[5]
    Mgv = sys.argv[6]
    main(genomefasta,pcl1,pcl0,samp_name,Agv,Mgv)
    #, pcl1, pcl0, samp_name