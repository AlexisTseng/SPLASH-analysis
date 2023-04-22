
import sys
import os
import re
import gzip

########################
## this script is used for removing reads with more than 7 consecutive As in FastQ file
## in SA11, 3 occurance of 7xA, none for 8xA
########################


def main(infastq):
    cd = os.getcwd()
    fastQ_Acut_dir = os.path.join(cd, "results", "FastQ_Acut")
    if not os.path.isdir(fastQ_Acut_dir):
        os.mkdir(fastQ_Acut_dir)
    read_writeFastQ(infastq,fastQ_Acut_dir)


def read_writeFastQ(infastq,fastQ_Acut_dir):
    exp_name = os.path.split(os.path.splitext(os.path.abspath(infastq))[0])[1][:-14]
    with gzip.open(infastq,"rt") as rdfile, gzip.open(f"{fastQ_Acut_dir}/{exp_name}_trimmed-Acut.fastq.gz" ,"wt") as outFastQ:
        for i, line in enumerate(rdfile): #enumerate starts at zero
            line = line.strip()
            #line.strip() #for removing characters at the start and end -> here remove new-line character
            if i % 4 == 0: #id
                id = line
            elif (i-1) % 4 == 0: # seq
                seq = line
                #print(seq)
                if re.search( "AAAAAAAA",seq) == None: pass #if no 8A, keep seq as it is
                else:
                    seq = re.split("AAAAAAAA",seq, 1)[0] # if 8A present, cut seq
            elif (i-2) % 4 == 0: # sign
                sign = line
            elif (i-3) % 4 == 0: #qc
                qc = line[:len(seq)]
                if len(qc) != 0:
                    #print(f"{id}\n{seq}\n{sign}\n{qc}\n")
                    outFastQ.write(f"{id}\n{seq}\n{sign}\n{qc}\n")


if __name__ == "__main__":
    infastq = sys.argv[1]
    main(infastq)

exit 