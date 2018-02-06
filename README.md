# BioMicroCenter-PPR_Program Version5
Guaranteeing high quality next-generation sequencing (NGS) data in a rapidly changing environment is an ongoing challenge. The depreciation of specific metrics from Illumina's Sequencing Analysis Viewer (SAV) have made it more difficult to directly determine the baseline error rate of sequencing runs. We have created an open-source tool to construct the Percent Perfect Reads (PPR) plot previously provided by the Illumina sequencers. The PPR program is compatible with HiSeq2000/2500, MiSeq, and NextSeq500 instruments, and provides an alternative to Illumina's Q scores for determining run quality.

In version5, we upgraded dependencies to newer versions including:
bwa/0.7.16a
samtools/1.5
r/3.4.0
fastxtoolkit/0.0.13
bowtie2/2.3.2
bedtools/2.26.0
GNU Parallel 20170222

To run the tool, you need to install all dependencies to you UNIX/Lunix machine. Download PPR source code and uncompress it to a folder. Then include the folder to path.
Example: 
PPPQC_DB_PATH=/home/PPRv5SourceCode
export PATH=$PPPQC_DB_PATH:$PATH

Commands:
./PPPQC.sh -runpath path_of_outputs NGS_RunType path/fastq_file(s)

Explanation of commands:
./PPRQC.sh: run the shell script PPRQC.sh
-runpath: first argument (mandatory)
path_of_outputs: second argument (mandatory). The path of the output folder, such as /home/user/PPRQC/run1/output
NGS_RunType: third argument (mandatory). The options are 
      nextseq_paired_end
      nextseq_single_end
      miseq_paired_end
      miseq_single_end
      hiseq_paired_end
      hiseq_single_end
path/fastq_file(s):forth to the last arguments (mandatory). It is the paths and file names of the fastq input files. See Examples of commands section

Examples of commands:
Paired end NextSeq sequencing:
./PPPQC.sh -runpath /home/user/PPRQC/PairedEndNextSeq/Run1 nextseq_paired_end /home/sequencing/SequencingRead1.fq /home/sequencing/SequencingRead2.fq

Single end NextSeq sequencing:
./PPPQC.sh -runpath /home/user/PPRQC/SingleEndNextSeq/Run1 nextseq_single_end /home/sequencing/SequencingRead1.fq

Paired end MiSeq sequencing:
./PPPQC.sh -runpath /home/user/PPRQC/PairedEndMiSeq/Run1 miseq_paired_end /home/sequencing/SequencingRead1.fq /home/sequencing/SequencingRead2.fq

Single end MiSeq sequencing:
./PPPQC.sh -runpath /home/user/PPRQC/SingleEndMiSeq/Run1 miseq_single_end /home/sequencing/SequencingRead1.fq

Paired end HiSeq sequencing:
./PPPQC.sh -runpath /home/user/PPRQC/PairedEndMiSeq/Run1 hiseq_paired_end /home/sequencing/Lane1_1.fq /home/sequencing/Lane1_2.fq /home/sequencing/Lane2_1.fq /home/sequencing/Lane2_2.fq /home/sequencing/Lane3_1.fq /home/sequencing/Lane3_2.fq /home/sequencing/Lane4_1.fq /home/sequencing/Lane4_2.fq /home/sequencing/Lane5_1.fq /home/sequencing/Lane5_2.fq /home/sequencing/Lane6_1.fq /home/sequencing/Lane6_2.fq /home/sequencing/Lane7_1.fq /home/sequencing/Lane7_2.fq /home/sequencing/Lane8_1.fq /home/sequencing/Lane8_2.fq

Single end HiSeq sequencing:
./PPPQC.sh -runpath /home/user/PPRQC/PairedEndMiSeq/Run1 hiseq_single_end /home/sequencing/Lane1_1.fq /home/sequencing/Lane2_1.fq /home/sequencing/Lane3_1.fq /home/sequencing/Lane4_1.fq /home/sequencing/Lane5_1.fq /home/sequencing/Lane6_1.fq /home/sequencing/Lane7_1.fq /home/sequencing/Lane8_1.fq

Note:
The input fastq files cannot be zipped or tarred. Otherwise, the script will not work.

Results:
See a jpg file named by the time of job submission
