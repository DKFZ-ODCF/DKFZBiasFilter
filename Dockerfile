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
    mkdir -p /home/pcawg/results

ADD scripts/ /usr/local/bin

CMD /usr/local/bin/run_biasfilter.sh -q /home/pcawg/input.vcf /home/pcawg/tumor.bam /home/pcawg/hs37d5.fa
