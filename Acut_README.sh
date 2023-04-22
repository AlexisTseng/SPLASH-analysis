#/bin/bash

#data=/mnt/c/Users/jz531/"OneDrive - University of Cambridge"/partii/data
#results=/mnt/c/Users/jz531/"OneDrive - University of Cambridge"/partii/results
data=/mnt/c/Users/jz531/partii/data
results=/mnt/c/Users/jz531/partii/results
genfasta=$data/SA11_WT.fa
align_Acut=$results/alignment_Acut
fastp=$results/fastp
fastp_QC=$fastp/QC
fasqQ_Acut=$results/FastQ_Acut
resFstqc_Acut=$fasqQ_Acut/QC
AcutNcut=$results/FastQ_AcutNNNcut
AcutNcut_QC=$AcutNcut/QC

## Package Installation
# conda install -c conda-forge biopython
# conda create -n cutadaptenv cutadapt

#conda activate py
#conda install -c bioconda viennarna=2.5.0
#conda install -c conda-forge numpy=1.22.2
#conda install -c lb_arrakistx varna=3.93
#conda install -c conda-forge pandas=1.4.0
#conda install -c conda-forge bokeh=2.4.2
#conda install -c conda-forge matplotlib=3.5.1
#conda install -c bioconda circos=0.69.8



function main(){

    #echo "QC of input adap-trim files"
    #files=("2_SA11_0-NSP2_B"
    #"2_transcribing_DLP_B"
    #"4_SA11_20-NSP2_B" 
    #"5_transcribing_DLP_NSP2")
    #doFastQC1 files


    #echo "Step 1. Clean FastQ file. Remove polyA"
    #files=("2_SA11_0-NSP2_B"
    #"2_transcribing_DLP_B"
    #"4_SA11_20-NSP2_B" 
    #"5_transcribing_DLP_NSP2")
    #cutA files

    #echo "QC of A-cut files"
    #files=("2_SA11_0-NSP2_B"
    #"2_transcribing_DLP_B"
    #"4_SA11_20-NSP2_B" 
    #"5_transcribing_DLP_NSP2")
    #doFastQC2 files

    echo "Remove 5' NNN"
    files=("2_SA11_0-NSP2_B"
    "2_transcribing_DLP_B"
    "4_SA11_20-NSP2_B" 
    "5_transcribing_DLP_NSP2")
    #cutNNN files



    echo "QC of A-cut files"
    files=("2_SA11_0-NSP2_B"
    "2_transcribing_DLP_B"
    "4_SA11_20-NSP2_B" 
    "5_transcribing_DLP_NSP2")
    doFastQC3 files



    #echo "Step 2. Run Alignment"
    #files=(
    ##"2_SA11_0-NSP2_B"
    #"2_transcribing_DLP_B" 
    #"4_SA11_20-NSP2_B" 
    #"5_transcribing_DLP_NSP2")
    #doAlign files
#
    #echo "Step 3. create bam files"
    #files=("2_SA11_0-NSP2_B" 
    #"2_transcribing_DLP_B"
    #"4_SA11_20-NSP2_B" 
    #"5_transcribing_DLP_NSP2")
    #mkBam files

    #echo "Step 5. classifying alignment results"
    ##def main(genomefasta,inbam,inJunc):
    #files=(
    #"2_SA11_0-NSP2_B" 
    #"2_transcribing_DLP_B"
    #"4_SA11_20-NSP2_B" 
    #"5_transcribing_DLP_NSP2")
    #callPython_class files

    
    #echo "Step 6. count and visualise interactions mapped to each segment"
    #NSP_1_file=("5_transcribing_DLP_NSP2" "4_SA11_20-NSP2_B" "5_transcribing_DLP_NSP2" "2_transcribing_DLP_B")
    #NSP_0_file=("2_transcribing_DLP_B" "2_SA11_0-NSP2_B" "4_SA11_20-NSP2_B" "2_SA11_0-NSP2_B")
    #sample_name=("DLP" "In_Vitro" "InVitro_VS_DLP_withNSP2" "InVitro_VS_DLP_withoutNSP2")
    #callPython_count NSP_1_file NSP_0_file sample_name

    #echo "Step 7. Mapping Reads to Interaction Windows"
    #NSP_1_file=("5_transcribing_DLP_NSP2" "4_SA11_20-NSP2_B" "5_transcribing_bDLP_NSP2" "2_transcribing_DLP_B")
    #NSP_0_file=("2_transcribing_DLP_B" "2_SA11_0-NSP2_B" "4_SA11_20-NSP2_B" "2_SA11_0-NSP2_B")
    #sample_name=("DLP" "In_Vitro" "InVitro_VS_DLP_withNSP2" "InVitro_VS_DLP_withoutNSP2")
    #callPython_window genfasta NSP_1_file NSP_0_file sample_name

    #echo " Step 8. Test Different Ag and Mg value For Whole-Seg"
    #NSP_1_file=("5_transcribing_DLP_NSP2" "4_SA11_20-NSP2_B" "5_transcribing_DLP_NSP2" "2_transcribing_DLP_B")
    #NSP_0_file=("2_transcribing_DLP_B" "2_SA11_0-NSP2_B" "4_SA11_20-NSP2_B" "2_SA11_0-NSP2_B")
    #sample_name=("DLP" "In_Vitro" "InVitro_VS_DLP_withNSP2" "InVitro_VS_DLP_withoutNSP2")
    #Agv=("0.05")
    #Mgv=("0.3")
    #callPython_count2 NSP_1_file NSP_0_file sample_name Agv Mgv

    #echo " Step 9. Interaction Windows"
    #NSP_1_file=("5_transcribing_DLP_NSP2" "2_transcribing_DLP_B")
    #NSP_0_file=("4_SA11_20-NSP2_B" "2_SA11_0-NSP2_B")
    #sample_name=("InVitro_VS_DLP_withNSP2" "InVitro_VS_DLP_withoutNSP2")
    #Agv=("0.05")
    #Mgv=("0.3")
    #callPython_window2 genfasta NSP_1_file NSP_0_file sample_name Agv Mgv

    #echo " Step 9. Test Different Ag and Mg value For Interaction Windows"
    #NSP_1_file=("5_transcribing_DLP_NSP2" "4_SA11_20-NSP2_B" "5_transcribing_DLP_NSP2" "2_transcribing_DLP_B")
    #NSP_0_file=("2_transcribing_DLP_B" "2_SA11_0-NSP2_B" "4_SA11_20-NSP2_B" "2_SA11_0-NSP2_B")
    #sample_name=("DLP" "In_Vitro" "InVitro_VS_DLP_withNSP2" "InVitro_VS_DLP_withoutNSP2")
    #Agv=("0.05")
    #Mgv=("0.3")
    #callPython_window2 genfasta NSP_1_file NSP_0_file sample_name Agv Mgv


}

###############################################
## functions


function doFastQC1() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        fastqc -o $fastp_QC --noextract -t 8 $fastp/${one[1]}-trimmed.fastq.gz
    done
}



function cutA() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        python ana_draft_Acut_4.py $fastp/${one[1]}-trimmed.fastq.gz
    done
}


function doFastQC2() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        fastqc -o $resFstqc_Acut --noextract -t 8 $fasqQ_Acut/${one[i]}_trimmed-Acut.fastq.gz
    done
}



function cutNNN() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        cutadapt -m 15 -u 3 $fasqQ_Acut/${one[i]}_trimmed-Acut.fastq.gz > $AcutNcut/${one[i]}_trimmed-Acut-Ncut.fastq.gz
    done
}

function doFastQC3() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        fastqc -o $AcutNcut_QC --noextract -t 8 $AcutNcut/${one[i]}_trimmed-Acut-Ncut.fastq.gz
    done
}




#STAR --runMode genomeGenerate --runThreadN 4 --genomeSAindexNbases 6 --genomeDir "${results}"/SA11_WT_index --genomeFastaFiles "${data}"/SA11_WT.fa 

function doAlign() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        rm -rd ~/STARtemp
        STAR --runMode alignReads --chimSegmentMin 20 --outSAMorder PairedKeepInputOrder  --runThreadN 6 --outFileNamePrefix $align_Acut/${one[i]}/${one[i]} --genomeDir $results/SA11_WT_index --readFilesIn $fasqQ_Acut/${one[i]}_trimmed-Acut.fastq.gz --outFilterMismatchNoverLmax 0.1 --readFilesCommand zcat --outTmpDir ~/STARtemp --chimOutType Junctions --chimOutJunctionFormat 1
    done
}
# minsegmentlength 10 20 40 
# PairedKeepInputOrder -> 
# --outFilterMismatchNoverLmax 0.1 0.3
# STAR wants to create a directory. Remove the directory first so that STAR does not interfere with it

function mkBam() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        samtools view -bh $align_Acut/${one[i]}/${one[i]}Aligned.out.sam | samtools sort -o $align_Acut/${one[i]}/${one[i]}Aligned.out.bam

        samtools index $align_Acut/${one[i]}/${one[i]}Aligned.out.bam $align_Acut/${one[i]}/${one[i]}Aligned.out.bai
    done
}
#conda install -c bioconda samtools=1.15

function callPython_class() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        python ana_draft_2.py $genfasta $align_Acut/${one[i]}/${one[i]}Aligned.out.bam $align_Acut/${one[i]}/${one[i]}*.junction
    done
}


function callPython_count() {
    local -n one=$1
    local -n two=$2
    local -n tre=$3
    for i in "${!one[@]}"
    do
        python ana_compare.py $align_Acut/${one[i]}/${one[i]}Aligned.out $align_Acut/${two[i]}/${two[i]}Aligned.out ${tre[i]}
    done
}

function callPython_count2() {
    local -n one=$1
    local -n two=$2
    local -n tre=$3
    local -n cua=$4
    local -n cin=$5
    for i in "${!one[@]}"
    do
        for k in "${!cin[@]}"
        do
            for j in "${!cua[@]}"
            do
                python ana_compare.py $align/${one[i]}/${one[i]}Aligned.out $align/${two[i]}/${two[i]}Aligned.out ${tre[i]} ${cua[j]} ${cin[k]}
            done
        done
    done
}

function callPython_window() {
    local -n one=$1
    local -n two=$2
    local -n tre=$3
    local -n cua=$4
    for i in "${!two[@]}"
    do
        python ana_window_2.py $genfasta $align/${two[i]}/${two[i]}Aligned.out $align/${tre[i]}/${tre[i]}Aligned.out ${cua[i]}
    done
}

function callPython_window2() {
    local -n one=$1
    local -n two=$2
    local -n tre=$3
    local -n cua=$4
    local -n cin=$5
    local -n sei=$6
    for i in "${!two[@]}"
    do
        for j in "${!cin[@]}"
        do
            python ana_Acut_window.py $genfasta $align_Acut/${two[i]}/${two[i]}Aligned.out $align_Acut/${tre[i]}/${tre[i]}Aligned.out ${cua[i]} ${cin[j]} ${sei[j]}

        done
    done
}

#########################
## call main
main

exit



















