FROM ubuntu:xenial
MAINTAINER Ivo Buchhalter @ DKFZ

RUN \
    umask 000 && \
    apt-get -y update && \
    apt-get -y upgrade && \
    apt-get -y install apt-utils && \
    apt-get -y install apt-utils apt-transport-https && \
    apt-get -y install autoconf make && \
    apt-get -y install build-essential && \
    apt-get -y install zlibc zlib1g zlib1g-dev && \
    apt-get -y install libncurses5-dev && \
    apt-get -y install sudo && \
    apt-get -y install wget && \
    apt-get -y install git && \
    apt-get -y install unzip && \
    apt-get -y install vim 

RUN \
    apt-get -y install python && \
    apt-get -y install python-pysam && \
    apt-get -y install python-numpy && \
    apt-get -y install python-scipy && \
    apt-get -y install python-matplotlib

RUN \
    mkdir -p /home/pcawg/scripts && \
    mkdir -p /home/pcawg/results

ADD scripts/ /home/pcawg/scripts

ENTRYPOINT python /home/pcawg/scripts/biasFilter.py -q --mapq=1 --baseq=1 --tempFolder=/home/pcawg/ /home/pcawg/input.vcf /home/pcawg/tumor.bam /home/pcawg/hs37d5.fa /home/pcawg/filtered.vcf && mv /home/pcawg/filtered.vcf /home/pcawg/results/ && mv /home/pcawg/filtered_qcSummary /home/pcawg/results/

# docker run -v /ibios/co01/buchhalt/gits/mixed_projects/PanCanQC/database_files/:/home/pcawg/data/database_files/ -v /ibios/co01/buchhalt/temp/tumor_MD-207_merged.mdup.bam:/home/pcawg/data/bams/tumor.bam -v /ibios/co01/buchhalt/temp/tumor_MD-207_merged.mdup.bam:/home/pcawg/data/bams/control.bam -v /ibios/co01/buchhalt/gits/mixed_projects/PanCanQC/results/:/home/pcawg/data/results/ pcawg-qc
