# Applied Bioinformatics Group
# YLoc Docker Image
#
# Philipp Thiel

FROM ubuntu:18.04

# Update package repository
RUN apt-get update
RUN apt-get -y upgrade


# ----------------------------------------------------------
# Install Useful stuff
# ----------------------------------------------------------
RUN apt-get install -y -q dirmngr software-properties-common \
                          perl python python3 vim wget


# ----------------------------------------------------------
# Install and BLAST
# ----------------------------------------------------------
RUN apt-get install -y ncbi-blast+


# ----------------------------------------------------------
# Install YLoc
# ----------------------------------------------------------
RUN apt-get install -y default-jre

ADD YLoc /YLoc
WORKDIR /YLoc

# ----------------------------------------------------------
# Test YLoc
# ----------------------------------------------------------

RUN python yloc.py test.fasta "YLoc-LowRes* Animals"
