FROM ubuntu:18.04 

COPY ./pyscripts/ /data/pyscripts
COPY ./nfNeST_ver03.nf/ /data/nfNeST_ver03.nf
COPY ./6genes_ver03/ /data/6genes_ver03
COPY ./testrun/ /data/testrun
WORKDIR /data/

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git
RUN apt-get install -y wget
RUN apt-get update && apt-get install -y \
  curl \
  unzip \
  perl \
  openjdk-11-jre-headless
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip
RUN unzip gatk-4.2.0.0.zip
RUN curl -s https://get.nextflow.io | bash
#RUN apt install -y picard-tools
RUN chmod +x /data/nextflow
RUN chmod 777 /data/nextflow
RUN mv nextflow /bin/
RUN wget https://github.com/broadinstitute/picard/releases/download/2.25.6/picard.jar
RUN chmod 777 picard.jar
#RUN unzip picard-tools-1.124.zip
RUN apt-get update && apt-get install -y fastqc
RUN apt-get update && apt-get install -y bowtie2
RUN apt-get update && apt-get install -y freebayes
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.7 \
    python3-pip \
    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN pip3 install -U pip setuptools 
RUN pip3 install pandas
RUN pip3 install matplotlib
RUN pip3 install pysam
RUN pip3 install xlrd==1.2.0
ENV BBMAP_VERSION 38.20
ENV BBMAP_DIR /working/bbmap/
RUN wget -O bbmap.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VERSION}.tar.gz/download \
    && tar -zxf bbmap.tar.gz \
    && rm bbmap.tar.gz
RUN apt-get update && apt-get -y upgrade && \
	apt-get install -y build-essential wget \
		libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
	tar jxf samtools-1.13.tar.bz2 && \
	rm samtools-1.13.tar.bz2 && \
	cd samtools-1.13 && \
	./configure --prefix $(pwd) && \
	make

RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
	tar jxf bcftools-1.9.tar.bz2 && \
	rm bcftools-1.9.tar.bz2 && \
	cd bcftools-1.9 && \
	make install

ENV package_version 4_3t
ADD https://downloads.sourceforge.net/project/snpeff/snpEff_v${package_version}_core.zip ./tmp/
RUN cd ./tmp/ && unzip snpEff_v${package_version}_core.zip \
         && rm snpEff_v${package_version}_core.zip
RUN apt-get update && apt-get install -y locales
RUN locale-gen "en_US.UTF-8"
RUN update-locale LC_ALL="en_US.UTF-8"
ENV LANG="en_US.UTF-8"
ENV PATH=/data/samtools-1.13/:$PATH
ENV PATH=$PATH:/data/nextflow