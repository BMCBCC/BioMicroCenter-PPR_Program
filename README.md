# BioMicroCenter-PPR_Program Version5
Guaranteeing high quality next-generation sequencing (NGS) data in a rapidly changing environment is an ongoing challenge. The depreciation of specific metrics from Illumina's Sequencing Analysis Viewer (SAV) has made it more difficult to directly determine the baseline error rate of sequencing runs. We have created an open-source tool to construct the Percent Perfect Reads (PPR) plot previously provided by the Illumina sequencers. The PPR program is compatible with HiSeq2000/2500, MiSeq, and NextSeq500 instruments, and provides an alternative to Illumina's Q scores for determining run quality.

How to run BioMicroCenter-PPR_Program:
Please refer Manual.txt 

Note:
The input fastq files cannot be zipped or tarred. Otherwise, the script will not work.

Results:
See a png file

Test input fastq file and expected output:
Example.fq.gz is an example Hiseq single end test fastq input file. This example file is compressed. Make sure to unzip it before your test run. Otherwise, the tool will not work. The expected output basing on the above test file should look like ExampleL1.jpg.If you see the same image, you have all the dependencies installed correctly.
