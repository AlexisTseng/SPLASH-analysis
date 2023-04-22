#/bin/bash

#data=/mnt/c/Users/jz531/"OneDrive - University of Cambridge"/partii/data
#results=/mnt/c/Users/jz531/"OneDrive - University of Cambridge"/partii/results
data=/mnt/c/Users/alexi/partii/data
results=/mnt/c/Users/alexi/partii/results
genfasta=$data/SA11_WT.fa
align=$results/alignment
fastp=$results/fastp
FastQC=$results/FastQC

## Step 0. setting up
#conda update conda
#conda update python
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda create -n partii python=3.11.0 
##force the lastest python version into partii

#conda activate partii
#conda install -c bioconda fastqc
#conda install -c bioconda fastp=0.23.2
## in new pc fastp=0.22.0
#conda install -c bioconda samtools=1.15
## in new pc samtools = 1.6
#conda install -c bioconda star=2.7.10a
## in new pc star=2.7.10b
#hisat2 -x ./results/hisat_index/S -U ./results/fastp/2_SA11_0-NSP2_B-trimmed.fastq -S ./results/hisat_alignment/hisat2_NSP.sam -p 4
#conda create -n py python=3.10
#conda activate py
#conda install pysam
## version = 2.35
#conda install -c bioconda samtools=1.15
#samtools view results/alignment/DLP/STAR20DLP*.bam | head



function main(){

    echo "Access Raw data & QC"
    files=("2_SA11_0-NSP2_B" 
    "2_transcribing_DLP_B" 
    "4_SA11_20-NSP2_B" 
    "5_transcribing_DLP_NSP2")
    #doFastQC files

    echo "Step 1. Trimming reads"
    files=("2_SA11_0-NSP2_B" 
    "2_transcribing_DLP_B" 
    "4_SA11_20-NSP2_B" 
    "5_transcribing_DLP_NSP2")
    #Trim files
    
    echo "Access trimmed reads and do QC"
    files=("2_SA11_0-NSP2_B-trimmed" 
    "2_transcribing_DLP_B-trimmed" 
    "4_SA11_20-NSP2_B-trimmed" 
    "5_transcribing_DLP_NSP2-trimmed")
    #doFastQC files


    echo "Step 2. Generate genome index"
    #GenIndex genfasta
    
    echo "Step 3. Run Alignment"
    files=("2_SA11_0-NSP2_B" 
    "2_transcribing_DLP_B" 
    "4_SA11_20-NSP2_B" 
    "5_transcribing_DLP_NSP2")
    #doAlign files

    echo "Step 4. create bam files"
    files=("2_SA11_0-NSP2_B" 
    "2_transcribing_DLP_B" 
    "4_SA11_20-NSP2_B" 
    "5_transcribing_DLP_NSP2")
    #mkBam files

##################################
## default
##################################
    echo "Step 5. classifying alignment results"
    #def main(genomefasta,inbam,inJunc):
    files=("2_SA11_0-NSP2_B" 
    "2_transcribing_DLP_B" 
    "4_SA11_20-NSP2_B" 
    "5_transcribing_DLP_NSP2")
    #callPython_class files
    
    echo "Step 6. count and visualise interactions mapped to each segment"
    NSP_1_file=("5_transcribing_DLP_NSP2" "4_SA11_20-NSP2_B" "5_transcribing_DLP_NSP2" "2_transcribing_DLP_B")
    NSP_0_file=("2_transcribing_DLP_B" "2_SA11_0-NSP2_B" "4_SA11_20-NSP2_B" "2_SA11_0-NSP2_B")
    sample_name=("DLP" "In_Vitro" "InVitro_VS_DLP_withNSP2" "InVitro_VS_DLP_withoutNSP2")
    #callPython_count NSP_1_file NSP_0_file sample_name

    echo "Step 7. Mapping Reads to Interaction Windows"
    NSP_1_file=("5_transcribing_DLP_NSP2" "4_SA11_20-NSP2_B" "5_transcribing_bDLP_NSP2" "2_transcribing_DLP_B")
    NSP_0_file=("2_transcribing_DLP_B" "2_SA11_0-NSP2_B" "4_SA11_20-NSP2_B" "2_SA11_0-NSP2_B")
    sample_name=("DLP" "In_Vitro" "InVitro_VS_DLP_withNSP2" "InVitro_VS_DLP_withoutNSP2")
    #callPython_window genfasta NSP_1_file NSP_0_file sample_name

    echo " Step 8. Test Different Ag and Mg value For Whole-Seg"
    NSP_1_file=("5_transcribing_DLP_NSP2" "4_SA11_20-NSP2_B" "5_transcribing_DLP_NSP2" "2_transcribing_DLP_B")
    NSP_0_file=("2_transcribing_DLP_B" "2_SA11_0-NSP2_B" "4_SA11_20-NSP2_B" "2_SA11_0-NSP2_B")
    sample_name=("DLP" "In_Vitro" "InVitro_VS_DLP_withNSP2" "InVitro_VS_DLP_withoutNSP2")
    Agv=("0.05")
    Mgv=("0.3")
    #callPython_count2 NSP_1_file NSP_0_file sample_name Agv Mgv

    #echo " Step 9. Interaction Windows"
    #NSP_1_file=("5_transcribing_DLP_NSP2" "2_transcribing_DLP_B")
    #NSP_0_file=("4_SA11_20-NSP2_B" "2_SA11_0-NSP2_B")
    #sample_name=("InVitro_VS_DLP_withNSP2" "InVitro_VS_DLP_withoutNSP2")
    #Agv=("0.05")
    #Mgv=("0.3")
    #callPython_window2 genfasta NSP_1_file NSP_0_file sample_name Agv Mgv

    echo " Step 9. Test Different Ag and Mg value For Interaction Windows"
    NSP_1_file=("5_transcribing_DLP_NSP2" "4_SA11_20-NSP2_B" "5_transcribing_DLP_NSP2" "2_transcribing_DLP_B")
    NSP_0_file=("2_transcribing_DLP_B" "2_SA11_0-NSP2_B" "4_SA11_20-NSP2_B" "2_SA11_0-NSP2_B")
    sample_name=("DLP" "In_Vitro" "InVitro_VS_DLP_withNSP2" "InVitro_VS_DLP_withoutNSP2")
    Agv=("0.05")
    Mgv=("0.3")
    callPython_window2 genfasta NSP_1_file NSP_0_file sample_name Agv Mgv

}

###############################################
## functions

function doFastQC() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        fastqc -o $FastQC --noextract -t 8 $data/${one[i]}.fastq.gz
    done
}


function Trim() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        fastp -i $data/${one[i]}.fastq.gz -o $fastp/${one[i]}-trimmed.fastq.gz -D -h $fastp/${one[i]}-trimmed.html -j $fastp/${one[i]}-trimmed.json
    done
}

function GenIndex() {
    local -n one=$1
    STAR --runMode genomeGenerate --runThreadN 4 --genomeSAindexNbases 6 --genomeDir $results/SA11_WT_index --genomeFastaFiles $genfasta
}
#STAR --runMode genomeGenerate --runThreadN 4 --genomeSAindexNbases 6 --genomeDir "${results}"/SA11_WT_index --genomeFastaFiles "${data}"/SA11_WT.fa 

function doAlign() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        rm -rd ~/STARtemp
        STAR --runMode alignReads --chimSegmentMin 20 --outSAMorder PairedKeepInputOrder  --runThreadN 6 --outFileNamePrefix $align/${one[i]}/${one[i]} --genomeDir $results/SA11_WT_index --readFilesIn $fastp/${one[i]}-trimmed.fastq.gz --outFilterMismatchNoverLmax 0.1 --readFilesCommand zcat --outTmpDir ~/STARtemp --chimOutType Junctions --chimOutJunctionFormat 1
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
        samtools view -bh $align/${one[i]}/${one[i]}Aligned.out.sam | samtools sort -o $align/${one[i]}/${one[i]}Aligned.out.bam

        samtools index $align/${one[i]}/${one[i]}Aligned.out.bam $align/${one[i]}/${one[i]}Aligned.out.bai
    done
}
#conda install -c bioconda samtools=1.15

function callPython_class() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        python ana_draft_2.py $genfasta $align/${one[i]}/${one[i]}Aligned.out.bam $align/${one[i]}/${one[i]}*.junction
    done
}

function callPython_classCutA() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        python ana_draft_3.py $genfasta $align/${one[i]}/${one[i]}Aligned.out.bam $align/${one[i]}/${one[i]}*.junction
    done
}

function callPython_count() {
    local -n one=$1
    local -n two=$2
    local -n tre=$3
    for i in "${!one[@]}"
    do
        python ana_compare.py $align/${one[i]}/${one[i]}Aligned.out $align/${two[i]}/${two[i]}Aligned.out ${tre[i]}
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
            python ana_window_3.py $genfasta $align/${two[i]}/${two[i]}Aligned.out $align/${tre[i]}/${tre[i]}Aligned.out ${cua[i]} ${cin[j]} ${sei[j]}

        done
    done
}

#########################
## call main
main

exit



















