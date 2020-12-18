FROM ubuntu:20.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu


RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y vim git wget build-essential python3 python3-dev python3-pip python3-matplotlib

RUN pip3 install whatshap==1.0 numpy matplotlib

COPY chr_lengths /root/chr_lengths
COPY merge_results.py /root/merge_results.py

WORKDIR /data
