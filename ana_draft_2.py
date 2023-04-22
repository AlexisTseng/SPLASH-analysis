
import sys
import os
import pysam
import inspect
import re
import pickle
import gzip

# pysam is for reading sam, bam files.

# in the bash script, 
# data = /mnt/c/data...
# result = /mnt/c/results
# we call python script out with
# python script.py $data/fast.fa $results/fastres.fa
# so $data.. will be the first argument -> sys.argv.[1], assigned to variable infasta


def main(genomefasta,inbam,inJunc):
    # input genome fasta file, return genome dictionary 
    print("Read genome fasta")
    gen_dict = readGenome(genomefasta)
    #print(gen_dict["SA11_S5"]) 
    picklefile = os.path.join(os.path.split(os.path.abspath(inbam))[0], os.path.split(os.path.splitext(os.path.abspath(inbam))[0])[1]) #creating a pcl file in the same directory as inbam, and with the same name as inbam file. First part = path, second path = name
    try: 
        interLRI_list, intraLRI_listdr, SRI_list,exp_name,bam_cnt,junc_cnt,totread_bp, totread_ap, per_inter, per_intra, per_SRI = loadData(picklefile)
    except:
        print("Read bam file")
        SRI_list, intraLRI_listdr,bam_cnt, exp_name = readBam(inbam,gen_dict)
        #print("SRI_list", SRI_list[:20])
        #print("intraLRI_listdr", intraLRI_listdr[:20])
        #print(bam_list[:1][0:])

        print("Read chimeric junction file")
        interLRI_list, intraLRI_listdr,SRI_list,junc_cnt = readChimJunc(inJunc,intraLRI_listdr,SRI_list)

        print("classification summary")
        totread_bp, totread_ap, per_inter, per_intra, per_SRI = summary(bam_cnt, junc_cnt, interLRI_list, intraLRI_listdr, SRI_list)

        saveData((interLRI_list, intraLRI_listdr, SRI_list,exp_name,bam_cnt,junc_cnt,totread_bp, totread_ap, per_inter, per_intra, per_SRI), picklefile)


    print(f"For {exp_name} '\n'total read before processing: {totread_bp} '\n'total interactions after processing: {totread_ap} '\n'% of inter-LRI: {per_inter} '\n'% of intra-LRI: {per_intra}'\n'% of SRI: {per_SRI}")

    #print(f"interlist {interLRI_list[:10]}'\n'intralist {intraLRI_listdr[:10]}'\n' SRI {SRI_list[:10]}")

    #print("seq_name",seq_name[0:5])
    #print("chim_name",chim_name[0:5])
    
    #print("Find equal names")
    #s = findrep(seq_name,chim_name)
    #print("s",s)
    
    #print(type(chim_name),chim_name[0:5],type(seq_name),seq_name[0:5])
    #for counting the number of chimeric CIGAR
    #interLRI_list, intraLRI_listdr, cnt = readChimJunc(inJunc,intraLRI_listdr)
    #print("cnt",cnt)
    #print(chimCIGAR)
    #print("interLRI_list", interLRI_list)
    #print("intraLRI_listdr", intraLRI_listdr)

        
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
            
def readBam(inbam,gen_dict):
    exp_name = os.path.split(os.path.splitext(os.path.abspath(inbam))[0])[1]
    with pysam.AlignmentFile(inbam,"rb") as ib:
        SRI_list = []
        intraLRI_listdr = []
        seq_name = []
        rname = ""
        comp = ()
        CIGAR = ""
        bam_cnt = 0
        for it,line in enumerate(ib):
            #seq_name.append(line.query_name)
            bam_cnt += 1
            rname = line.reference_name
            i = int(line.reference_start) # starting position as integer 0-index
            i,j,k,l = readCIGAR(i,line.cigarstring)
            if k == 0 and l == 0:
                comp = (rname,i,j)
                SRI_list.append(comp)
                #print("bam sequence",line.query_alignment_sequence)
                #print("indexed sequence", comp, gen_dict[rname][i:j])
            elif k-j <20:
                comp = (rname, i,l)
                SRI_list.append(comp)
            elif j-i >= 20 and l-k >= 20:
                comp = (rname, i, j, k, l)
                intraLRI_listdr.append(comp)
            elif j-i >= 20 and l-k < 20:
                comp = (rname, i,j)
                SRI_list.append(comp)
            elif j-i < 20 and l-k >= 20:
                comp = (rname, k,l)
                SRI_list.append(comp)
        #print("bam",bam_cnt)
        return SRI_list, intraLRI_listdr, bam_cnt, exp_name

def readChimJunc(inJunc, intraLRI_listdr,SRI_list):
    with open(inJunc, "r") as ij:
        interLRI_list = []
        #chim_name =[]
        rname1 = ""
        rname2 = ""
        comp = ()
        comp2 = ()
        comp3 = ()
        CIGAR1 = ""
        CIGAR2= ""
        cnt = 0
        junc_cnt = 0
        for it,line in enumerate(ij):
            junc_cnt += 1
            line = line.split()
            if line[0] == "#":
                continue
            rname1 = line[0]
            rname2 = line[3]
            #chim_name.append(line[9])
            rstart1 = int(line[10])
            rstart2 = int(line[12])
            CIGAR1 = line[11]
            CIGAR2 = line[13]
            ai,aj,ak,al = readCIGAR(rstart1,CIGAR1)
            bi,bj,bk,bl = readCIGAR(rstart2,CIGAR2)
            #print (rname1, CIGAR1,rstart1)
            if rname1 == rname2:#intraLRI & SRI
                if (ak,al) != (0,0) and (bk,bl) != (0,0):
                    pass
                elif (ak,al) != (0,0):
                    if ak-aj > 20: #gap >20
                        if aj-ai >= 20 and al-ak >=20: #both match >20
                            comp = (rname1, ai, aj, bi, bj)
                            comp2 = (rname1, ak, al, bi,bj)
                            comp3 = (rname1, ai, aj, ak, al)
                            intraLRI_listdr.append(comp)
                            intraLRI_listdr.append(comp2)
                            intraLRI_listdr.append(comp3)
                        elif aj-ai <20:
                            comp = (rname1, ak,al,bi,bj)
                            intraLRI_listdr.append(comp)
                        elif al-ak <20:
                            comp = (rname1,ai,aj,bi,bj)
                            intraLRI_listdr.append(comp)
                    else:
                        comp = (rname1,ai,al,bi,bj)
                        intraLRI_listdr.append(comp)
                        SRI_list.append(rname1, ai,al)
                elif (bk,bl) != (0,0):
                    if bk-bj > 20:
                        if bj-bi >= 20 and bl-bk >=20:
                            comp = (rname1, bi, bj, ai, aj)
                            comp2 = (rname1, bk, bl, ai,aj)
                            comp3 = (rname1, bi, bj, bk, bl)
                            intraLRI_listdr.append(comp)
                            intraLRI_listdr.append(comp2)
                            intraLRI_listdr.append(comp3)
                        elif bj-bi <20:
                            comp = (rname1, bk,bl,ai,aj)
                            intraLRI_listdr.append(comp)
                        elif bl-bk <20:
                            comp = (rname1,bi,bj,ai,aj)
                            intraLRI_listdr.append(comp)
                    else:
                        comp = (rname1,bi,bl,ai,aj)
                        intraLRI_listdr.append(comp)
                        SRI_list.append(rname1, bi,bl)
                else:
                    intraLRI_listdr.append((rname1,ai,aj,bi,bj))
            else: #interLRI
                if (ak,al) != (0,0) and (bk,bl) != (0,0):
                    pass
                elif (ak,al) != (0,0):
                    comp3 = (rname2, bi, bj)
                    if ak-aj > 20 :
                        if aj-ai >= 20 and al-ak >= 20:
                            comp = (rname1, ai, aj)
                            comp2 = (rname1,ak, al)
                            intraLRI_listdr.append((rname1,ai,aj,ak,al))
                            interLRI_list.append((comp,comp3))
                            interLRI_list.append((comp2,comp3))
                        elif aj-ai <20:
                            comp = (rname1, ak,al)
                            interLRI_list.append((comp,comp3))
                        elif al-ak <20:
                            comp = (rname1,ai,aj)
                            interLRI_list.append((comp,comp3))
                    else: 
                        comp = (rname1,ai,al)
                        interLRI_list.append((comp,comp3))
                        SRI_list.append(comp)
                elif (bk,bl) != (0,0):
                    comp3 = (rname1, ai, aj)
                    if bk-bj > 20 :
                        if bj-bi >= 20 and bl-bk >= 20:
                            comp = (rname2, bi, bj)
                            comp2 = (rname2,bk, bl)
                            intraLRI_listdr.append((rname2,bi,bj,bk,bl))
                            interLRI_list.append((comp,comp3))
                            interLRI_list.append((comp2,comp3))
                        elif bj-bi <20:
                            comp = (rname2, bk,bl)
                            interLRI_list.append((comp,comp3))
                        elif bl-bk <20:
                            comp = (rname2,bi,bj)
                            interLRI_list.append((comp,comp3))
                    else: 
                        comp = (rname2,bi,bl)
                        interLRI_list.append((comp,comp3))
                        SRI_list.append(comp)
                else:
                    comp = (rname1,ai,aj)
                    comp3 = (rname2,bi,bj)
                    interLRI_list.append((comp,comp3))
        #print("ju",junc_cnt)
        return interLRI_list, intraLRI_listdr, SRI_list, junc_cnt

def summary(bam_cnt,junc_cnt, interLRI_list, intraLRI_listdr, SRI_list):
    totread_bp = int(bam_cnt)+int(junc_cnt)
    leninter = len(interLRI_list)
    lenintra = len(intraLRI_listdr)
    lenSRI = len(SRI_list)
    totread_ap = leninter + lenintra+ lenSRI
    per_inter = leninter*100/totread_ap
    per_intra = lenintra*100/totread_ap
    per_SRI = lenSRI*100/totread_ap
    return totread_bp, totread_ap, per_inter, per_intra, per_SRI


def findrep(seq_name,chim_name):
    seq_name.sort()
    chim_name.sort()
    c = 0
    s = 0
    l = len(seq_name)
    for k,name in enumerate(seq_name):
        print(f"{k} of {l}", end="\r")
        for i in range(c,len(chim_name)):
            if name < chim_name[i]:
                continue
            elif name == chim_name[i]:
                s += 1
            else:
                c = i
                break
    print("")
    return s
# no overlap found between bam and chimJunction


def readCIGAR(i,CIGAR): 
    ## read CIGAR string
    # M = match, S = soft clipping, D = deletion, N = intron (longer deletion), I = insertion
    reCig = re.compile(r"([0-9]+)(M|S|D|N|I)")
    j, k, l = i, 0, 0
    for cig in reCig.finditer(CIGAR):
        val,key = int(cig.group(1)), cig.group(2)
        #if   key == "S"  and i == j: i += val; j += val
        if   key in "MD" and k == 0: j += val
        elif key in "MD" and k != 0: l += val
        elif key == "N"  and k == 0: k = j + val; l = j + val
    return i,j,k,l

def loadData(fname): 
    ## load data with pickle
    with open(f"{fname}.pcl", "r+b") as pcl_in:
        pcl_data = pickle.load(pcl_in)
    return pcl_data

def saveData(pcl_data, fname):
    ## save data with pickle
    with open(f"{fname}.pcl", "w+b") as pcl_out:
        pickle.dump(pcl_data, pcl_out , protocol=4)


if __name__ == "__main__":
    genomefasta = sys.argv[1]
    inbam = sys.argv[2]
    inJunc = sys.argv[3]
    #print(infastq)
    #exit()
    main(genomefasta, inbam, inJunc)

exit 
#some_dictionary["new key"] = 700

#line.strip() #for removing characters at the start and end