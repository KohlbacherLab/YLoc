# Applied Bioinformatics Group
# YLoc Docker Image
#
# Philipp Thiel

FROM ubuntu:18.04

# Update package repository
RUN apt-get update


# ----------------------------------------------------------
# Install Useful stuff
# ----------------------------------------------------------
RUN apt-get install -y -q dirmngr software-properties-common vim wget


# ----------------------------------------------------------
# Install YLoc and dependencies
# ----------------------------------------------------------
RUN  apt-get install -y ncbi-blast+ pftools default-jre
COPY YLoc/ YLoc/
RUN  chown -R www-data:www-data /YLoc
RUN  chmod -R 775 /YLoc

RUN mkdir /yljobs
RUN chmod 777 /yljobs

# ----------------------------------------------------------
# Test YLoc
# ----------------------------------------------------------
RUN python /YLoc/yloc.py /YLoc/test.fasta "YLoc-LowRes* Animals"


# ----------------------------------------------------------
# Setup YLoc Webservice
# ----------------------------------------------------------
RUN apt-get -y install mysql-server python-mysql.connector apache2
RUN a2enmod cgid
ADD webservice/apache2.conf         /etc/apache2/apache2.conf
ADD webservice/serve-cgi-bin.conf   /etc/apache2/conf-available/serve-cgi-bin.conf

COPY webservice/webloc.cgi  /var/www/html/cgi-bin/
COPY webservice/downloads/  /var/www/html/cgi-bin/downloads/
COPY webservice/images/     /var/www/html/cgi-bin/images/

RUN mkdir /webservice
ADD webservice/job_cleanup.sh      /webservice/job_cleanup.sh
ADD webservice/ylsetup.py          /webservice/ylsetup.py
ADD webservice/yloc_entrypoint.sh  /webservice/yloc_entrypoint.sh
ADD webservice/ylocdb.sql          /webservice/ylocdb.sql


EXPOSE 80

CMD ["sh", "/webservice/yloc_entrypoint.sh"]
