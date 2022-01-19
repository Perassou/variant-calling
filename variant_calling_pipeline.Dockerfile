FROM python:3.8-slim-buster
MAINTAINER Marie Peras

RUN apt-get update && apt-get install -y build-essential samtools wget libz-dev freebayes

# get bwa
RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2 \
    && bunzip2 bwa-0.7.17.tar.bz2 \
    && tar xvf bwa-0.7.17.tar \
    && cd bwa-0.7.17 && make \
    && cp bwa /usr/bin/bwa

#  build the bwa index
ARG genomes_dir
ADD $genomes_dir /genomes/
RUN mv /genomes/*.fasta /genomes/reference.fasta
RUN bwa index -p genomes_index /genomes/reference.fasta

RUN mkdir /data
WORKDIR "/data"
ENV PYTHONPATH '/'

COPY alignment/variant_calling.py /variant_calling.py
RUN chmod +x /variant_calling.py

# we want to manually specify the existing index we have; 
# the whole point of this docker container is to have the index ready to go
ENTRYPOINT ["/variant_calling.py", "-i", "/genomes_index", "-r", "/genomes/reference.fasta"]