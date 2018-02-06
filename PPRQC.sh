#!/bin/bash

downSample() {
    FASTQ=$1
    ORDER=$2
    l=$(wc -l $FASTQ | awk '{print $1}')
    if [ $l -gt 5500000 ]; then
        s=$(($l/44))
        for j in 1 2 3 4 5 6 7 8 9 10
        do
            start=$(($s*4*$j+400000))
            head -$start $FASTQ |tail -400000>>int${ORDER}.fq
        done
    else
        head -4000000 $FASTQ >int${ORDER}.fq
    fi
}

alignSingle() {

    #retrieve the first 30bp of the sequencing reads and align them
    fastx_trimmer -l 30 -i $1 -Q33 -o read1.fq
    bowtie2 -p 8 --reorder -x /home/Genomes/bowtie2_indexes/phiX -U read1.fq -S align.sam 2>&1

    #seperate phiX reas with no indels in the first 30 bp
    MD.pl align.sam T8_tmp.sam

    #seperate forward reads from reverse reads
    head -2 align.sam >sam_header
    cat sam_header T8_tmp.sam >Read1_h.sam
    samtools view -F 0x10 -S Read1_h.sam >phiX_Read1_f.sam
    samtools view -f 0x10 -S Read1_h.sam >phiX_Read1_r.sam

}

alignPair() {

    #retrieve the first 30bp of the sequencing reads and align them
    parallel --link fastx_trimmer -l 30 -i {1} -Q33 -o read{2}.fq ::: $1 $2 ::: 1 2
    bowtie2 -p 8 --reorder -x /home/Genomes/bowtie2_indexes/phiX -1 read1.fq -2 read2.fq -S align.sam 2>&1

    #seperate phiX reas with no indels in the first 30 bp
    MD.pl align.sam T8_tmp.sam

    #separate read1 from read2 and seperate forward reads from reverse reads
    head -2 align.sam >sam_header
    cat sam_header T8_tmp.sam >T8.sam

    samtools view -f 0x40 -S T8.sam >Read1.sam
    cat sam_header Read1.sam >Read1_h.sam
    samtools view -F 0x10 -S Read1_h.sam >phiX_Read1_f.sam
    samtools view -f 0x10 -S Read1_h.sam >phiX_Read1_r.sam

    samtools view -f 0x80 -S T8.sam >Read2.sam
    cat sam_header Read2.sam >Read2_h.sam
    samtools view -F 0x10 -S Read2_h.sam >phiX_Read2_f.sam
    samtools view -f 0x10 -S Read2_h.sam >phiX_Read2_r.sam

}

processSingleMiseqHiseq() {

    alignSingle $1

    #making sequence file of Phix forward/reverse reads
    mk_faf.pl $1 phiX_Read1_f.sam phiX_Read1_f_seq
    mk_far.pl $1 phiX_Read1_r.sam phiX_Read1_r_seq

    #making bed file to retrieve reference sequence
    bed_fs.pl phiX_Read1_f_seq phiX_Read1_f.sam phiX_Read1_f.bed
    bed_rs.pl phiX_Read1_r_seq phiX_Read1_r.sam phiX_Read1_r.bed

    #get reference fasta files from bed
    bedtools getfasta -fi /home/Genomes/bowtie2_indexes/phiX.fa -bed phiX_Read1_f.bed -fo phiX_Read1_f_ref.fa -name
    bedtools getfasta -fi /home/Genomes/bowtie2_indexes/phiX.fa -bed phiX_Read1_r.bed -fo phiX_Read1_r_ref.fa -name

    #format reference fasta file to table file for merging with real sequencing data in next steps
    format.pl phiX_Read1_f_ref.fa phiX_Read1_f_ref
    format.pl phiX_Read1_r_ref.fa phiX_Read1_r_ref

    #merging reference seq table with real sequencing data in next steps
    overlaphash_fs.pl phiX_Read1_f_ref phiX_Read1_f_seq phiX_Read1_f
    overlaphash_rs.pl phiX_Read1_r_ref phiX_Read1_r_seq phiX_Read1_r

    #filter away sequences with indels defined as 3 continuous nt mismatching between ref and seq
    filter_indel.pl phiX_Read1_f phiX_Read1_f_no_indel
    filter_indel.pl phiX_Read1_r phiX_Read1_r_no_indel

    #mismatch counting
    mismatch_count_v4.pl phiX_Read1_f phiX_Read1_r Read1_mm size_r1_c1

    cat $PPPQC_DB_PATH/header Read1_mm > $2
}

processPairMiseqHiseq() {

    alignPair $1 $2

    #making sequence file of forward/reverse reads
    parallel mk_fa1{1}.pl $1 phiX_Read1_{1}.sam phiX_Read1_{1}_seq ::: f r
    parallel mk_fa2{1}.pl $2 phiX_Read2_{1}.sam phiX_Read2_{1}_seq ::: f r

    #making bed file to retrieve reference sequence
    parallel bed_{2}.pl phiX_Read{1}_{2}_seq phiX_Read{1}_{2}.sam phiX_Read{1}_{2}.bed ::: 1 2 ::: f r

    #get reference fasta files from bed
    parallel bedtools getfasta -fi /home/Genomes/bowtie2_indexes/phiX.fa -bed phiX_Read{1}_{2}.bed -fo phiX_Read{1}_{2}_ref.fa -name ::: 1 2 ::: f r

    #format reference fasta file to table file for merging with real sequencing data in next steps
    parallel format.pl phiX_Read{1}_{2}_ref.fa phiX_Read{1}_{2}_ref ::: 1 2 ::: f r

    #merging reference seq table with real sequencing data in next steps
    parallel overlaphash_{2}.pl phiX_Read{1}_{2}_ref phiX_Read{1}_{2}_seq phiX_Read{1}_{2} ::: 1 2 ::: f r

    #filter away sequences with indels defined as 3 continuous nt mismatching between ref and seq
    parallel filter_indel.pl phiX_Read{1}_{2} phiX_Read{1}_{2}_no_indel ::: 1 2 ::: f r

    #mismatch counting
    mismatch_count_v4.pl phiX_Read1_f phiX_Read1_r Read1_mm size_r1
    mismatch_count_v4_r2.pl phiX_Read2_f phiX_Read2_r size_r1 Read2_mm

    LMMH=$3
    cat $PPPQC_DB_PATH/header Read1_mm Read2_mm  > ${LMMH}i
    fix.pl ${LMMH}i ${LMMH}

}

export -f downSample

# Process arguments
while [ $# -ne 0 ]; do
    case $1 in
        -runpath)
            PPPQC_RUN_PATH="$2"
            shift
            ;;
    esac
    shift
    break # We only need to get the -runpath argument
done

set -xv

# Creat the run path if it does not exist and clean it up if it does
if [ ! -d "$PPPQC_RUN_PATH" ]; then
    mkdir -p $PPPQC_RUN_PATH
else
    rm -rf $PPPQC_RUN_PATH/*
fi

cd $PPPQC_RUN_PATH

##########################################
#paired end NextSeq sequencing
if [ $1 = "nextseq_paired_end" ]; then
    shift

    alignPair $1 $2

    #seperate reads by camera
    parallel camera_parsing_seqname.pl phiX_Read{2}_{3}.sam camera{1}_phiX_Read{2}_{3} {1} ::: 1 2 3 4 5 6 ::: 1 2 ::: f r

    #making sequence file of each camera, Read1/2, forward/reverse reads
    parallel mk_fa1{2}.pl $1 camera{1}_phiX_Read1_{2} camera{1}_phiX_Read1_{2}_seq ::: 1 2 3 4 5 6 ::: f r
    parallel mk_fa2{2}.pl $2 camera{1}_phiX_Read2_{2} camera{1}_phiX_Read2_{2}_seq ::: 1 2 3 4 5 6 ::: f r

    #making bed file to retrieve reference sequence
    parallel bed_{3}.pl camera{1}_phiX_Read{2}_{3}_seq phiX_Read{2}_{3}.sam camera{1}_phiX_Read{2}_{3}.bed ::: 1 2 3 4 5 6 ::: 1 2 ::: f r

    #get reference fasta files from bed
    parallel bedtools getfasta -fi /home/Genomes/bowtie2_indexes/phiX.fa -bed camera{1}_phiX_Read{2}_{3}.bed -fo camera{1}_phiX_Read{2}_{3}_ref.fa -name ::: 1 2 3 4 5 6 ::: 1 2 ::: f r

    #format reference fasta file to table file for merging with real sequencing data in next steps
    parallel format.pl camera{1}_phiX_Read{2}_{3}_ref.fa camera{1}_phiX_Read{2}_{3}_ref ::: 1 2 3 4 5 6 ::: 1 2 ::: f r

    #merging reference seq table with real sequencing data in next steps
    parallel overlaphash_{3}.pl camera{1}_phiX_Read{2}_{3}_ref camera{1}_phiX_Read{2}_{3}_seq camera{1}_phiX_Read{2}_{3} ::: 1 2 3 4 5 6 ::: 1 2 ::: f r

    #filter away sequences with indels defined as 3 continuous nt mismatching between ref and seq
    parallel filter_indel.pl camera{1}_phiX_Read{2}_{3} camera{1}_phiX_Read{2}_{3}_no_indel ::: 1 2 3 4 5 6 ::: 1 2 ::: f r

    #mismatch counting
    parallel mismatch_count_v4.pl    camera{1}_phiX_Read1_f camera{1}_phiX_Read1_r camera{1}_Read1_mm size_r1_c{1} ::: 1 2 3 4 5 6
    parallel mismatch_count_v4_r2.pl camera{1}_phiX_Read2_f camera{1}_phiX_Read2_r size_r1_c{1} camera{1}_Read2_mm ::: 1 2 3 4 5 6

    parallel "cat $PPPQC_DB_PATH/header camera{1}_Read1_mm camera{1}_Read2_mm  >camera{1}_mmhi" ::: 1 2 3 4 5 6

    parallel fix.pl camera{1}_mmhi camera{1}_mmh ::: 1 2 3 4 5 6

    #plotting
    R --slave --vanilla <$PPPQC_DB_PATH/plot.r

    #################################################
    #single end NextSeq sequencing
elif [ $1 = "nextseq_single_end" ]; then
    shift

    alignSingle $1

    #seperate reads by camera
    parallel camera_parsing_seqname.pl phiX_Read1_{2}.sam camera{1}_phiX_Read1_{2} {1} ::: 1 2 3 4 5 6 ::: f r

    #making sequence file of each camera, forward/reverse reads
    parallel mk_fa{2}.pl $1 camera{1}_phiX_Read1_{2} camera{1}_phiX_Read1_{2}_seq ::: 1 2 3 4 5 6 ::: f r

    #making bed file to retrieve reference sequence
    parallel bed_{2}s.pl camera{1}_phiX_Read1_{2}_seq phiX_Read1_{2}.sam camera{1}_phiX_Read1_{2}.bed ::: 1 2 3 4 5 6 ::: f r 

    #get reference fasta files from bed
    parallel bedtools getfasta -fi /home/Genomes/bowtie2_indexes/phiX.fa -bed camera{1}_phiX_Read1_{2}.bed -fo camera{1}_phiX_Read1_{2}_ref.fa -name ::: 1 2 3 4 5 6 ::: f r

    #format reference fasta file to table file for merging with real sequencing data in next steps
    parallel format.pl camera{1}_phiX_Read1_{2}_ref.fa camera{1}_phiX_Read1_{2}_ref ::: 1 2 3 4 5 6 ::: f r

    #merging reference seq table with real sequencing data in next steps
    parallel overlaphash_{2}s.pl camera{1}_phiX_Read1_{2}_ref camera{1}_phiX_Read1_{2}_seq camera{1}_phiX_Read1_{2} ::: 1 2 3 4 5 6 ::: f r

    #filter away sequences with indels defined as 3 continuous nt mismatching between ref and seq
    parallel filter_indel.pl camera{1}_phiX_Read1_{2} camera{1}_phiX_Read1_{2}_no_indel ::: 1 2 3 4 5 6 ::: f r

    #mismatch counting
    parallel mismatch_count_v4.pl camera{1}_phiX_Read1_f camera{1}_phiX_Read1_r camera{1}_Read1_mm size_r1_c{1} ::: 1 2 3 4 5 6

    parallel "cat $PPPQC_DB_PATH/header camera{1}_Read1_mm >camera{1}_mmh" ::: 1 2 3 4 5 6

    #plotting
    R --slave --vanilla <$PPPQC_DB_PATH/plot.r

    ############################################################
    #paired end MiSeq sequencing
elif [ $1 = "miseq_paired_end" ]; then
    shift

    processPairMiseqHiseq $1 $2 mmh

    #plotting
    R --slave --vanilla <$PPPQC_DB_PATH/plot_miseq.r

    ############################################################
    #single end MiSeq sequencing
elif [ $1 = "miseq_single_end" ] ; then
    shift

    processSingleMiseqHiseq $1 mmh

    #plotting
    R --slave --vanilla <$PPPQC_DB_PATH/plot_miseq.r

    ############################################################
    #paired end HiSeq sequencing (essentially the same code as miseq_paired_end with looping on lane fastq files)
elif [ $1 = "hiseq_paired_end" ]; then
    shift
    n="$#"
    # downsample
    parallel --link downSample {1} {2} ::: "$@" ::: $(seq 1 $n)

    #loop through each lane and perform the related tasks
    i=1
    until [ $i -gt $n ];
    do
        j=$(($i+1))
        k=$(($j/2))
        processPairMiseqHiseq int$i.fq int$j.fq L${k}_mmh
        if [ $(wc -l L${k}_mmh | awk '{print $1}') -eq "1" ]; then
            rm L${k}_mmh
        fi
        i=`expr $i + 2`
    done

    #plotting
    Rscript $PPPQC_DB_PATH/plot_hiseq.r `ls *mmh|sort`

    ############################################################
    #single end HiSeq sequencing (essentially the same as miseq_single_end with looping on lane fastq files)
elif [ $1 = "hiseq_single_end" ]; then
    shift
    n="$#"
    # downsample
    parallel --link downSample {1} {2} ::: "$@" ::: $(seq 1 $n)

    # We need loop here to avoid different lane processing steps into each other if using parallel
    i=1
    until [ $i -gt $n ];
    do
        processSingleMiseqHiseq int${i}.fq L${i}_mmh
        if [ $(wc -l L${i}_mmh | awk '{print $1}') -eq "1" ]; then
            rm L${i}_mmh
        fi
        i=`expr $i + 1`
    done

    #plotting
    Rscript $PPPQC_DB_PATH/plot_hiseq.r `ls *mmh|sort`

fi
rm *sam *fq *seq *_r *_f *bed *mm *mmh size* *indel *fa *ref sam_header
