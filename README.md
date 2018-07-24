## YLoc
### An interpretable web server for predicting subcellular localization

Travis CI  
[![Build Status](https://travis-ci.org/KohlbacherLab/YLoc.svg?branch=master)](https://travis-ci.org/KohlbacherLab/YLoc)  

This repository mainly provides a containerized installation of YLoc using Docker.  
In order to install the software outside of a container please read the instructions   
in this document and the Dockerfile.

**Citing & Further Information**  

If you use YLoc please cite the following publications:

Briesemeister S., Rahnenführer J., and Kohlbacher O. (2010)  
[YLoc—an interpretable web server for predicting subcellular localization.](https://doi.org/10.1093/bioinformatics/btq115)  
Bioinformatics, 26. 1232-8

Briesemeister S., Rahnenführer J., and Kohlbacher O. (2010)  
[YLoc—an interpretable web server for predicting subcellular localization.](https://dx.doi.org/10.1093%2Fnar%2Fgkq477)  
Nucleic Acids Res., 38, W497–W502  
  
  
**Installation Requirements**  

- Linux OS
- Docker
- Python >= 2.7
- [BLAST (no legacy BLAST)](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [PrositeScan](https://prosite.expasy.org/scanprosite/)


**Installation**

The easiest option is to build the Docker image from this repository using the following steps:  
` $ git clone https://github.com/KohlbacherLab/YLoc.git`  
` $ docker build --no-cache -t <your_image_name> YLoc/`  

**YLoc Usage**  

YLoc general usage:  
` $ python yloc.py <fasta_file> <model_name> <prediction_id(optional)> <print_result(y/n)(optional)>`  

You can print the usage description and available models by executing  
` $ python yloc.py`  

**Running YLoc**  

You can either start your container interactively and run YLoc  
` $ docker run --rm -it <your_image_name> /bin/bash`  
`root@<some_hash>:/YLoc# python yloc.py test.fasta "some_model"`  


