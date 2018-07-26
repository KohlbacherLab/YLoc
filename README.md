## YLoc
### An interpretable web server for predicting subcellular localization

Travis CI  
[![Build Status](https://travis-ci.org/KohlbacherLab/YLoc.svg?branch=master)](https://travis-ci.org/KohlbacherLab/YLoc)  

This repository mainly provides a containerized installation of YLoc using Docker.  
In order to install the software elsewhere please read the instructions in the  
[README](YLoc/README.txt) file of the YLoc source folder and have a look in the Dockerfile.  


**Citing & Further Information**  

If you use YLoc please cite the following publications:

Briesemeister S., Rahnenführer J., and Kohlbacher O. (2010)  
[YLoc—an interpretable web server for predicting subcellular localization.](https://doi.org/10.1093/bioinformatics/btq115)  
Bioinformatics, 26. 1232-8

Briesemeister S., Rahnenführer J., and Kohlbacher O. (2010)  
[YLoc—an interpretable web server for predicting subcellular localization.](https://dx.doi.org/10.1093%2Fnar%2Fgkq477)  
Nucleic Acids Res., 38, W497–W502  
  
  
**Requirements**  

- Linux OS
- Docker


**Installation**

The easiest option is to build the Docker image from this repository using the following steps:  
` $ git clone https://github.com/KohlbacherLab/YLoc.git`  
` $ docker build --no-cache -t <your_image_name> YLoc/`  

**YLoc Usage (Interactively in container)**  

YLoc general usage:  
` $ python yloc.py <fasta_file> <model_name> <prediction_id(optional)> <print_result(y/n)(optional)>`  

You can print the usage description and available models by executing  
` $ python yloc.py`  

**Running YLoc Using Docker run**  

You can either start your container interactively and run YLoc  
` $ docker run --rm -it <your_image_name> /bin/bash`  
`root@<some_hash>:/YLoc# python yloc.py test.fasta "some_model"`  


