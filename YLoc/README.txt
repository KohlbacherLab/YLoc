YLoc - Interpretable Subcellular Localization Prediction
*stand-alone version, packed 02/01/2011, copyright Sebastian Briesemeister*
-----------------------------------------------------------------------------

This is the first attempt of a stand-alone version of YLoc. It basically equals
the web service version of YLoc but does not connect to a MySQL database to save
the prediction outcomes. Therefore, some lines were set into comments. Moreover,
the source code is not well documentated and consequently hard to understand.
If you want to make changes and get into troubles, please contact 
briese@informatik.uni-tuebingen.de. Please be also aware of the fact, that
the file sortsignals.py does contain a very large set of procedures to
calculate possible sorting signals (mostly from the literature). Most of them
are not used in the final models but were only used during the evaluation and
feature selection.

!!! If you re-use part of the codes for your own project, please contact
	the authors to get permission !!!

If you use YLoc please cite:

Sebastian Briesemeister, Jörg Rahnenführer, and Oliver Kohlbacher, (2010). 
Going from where to why - interpretable prediction of protein subcellular 
localization, Bioinformatics, 26(9):1232-1238.

Sebastian Briesemeister, Jörg Rahnenführer, and Oliver Kohlbacher, (2010).
YLoc - an interpretable web server for predicting subcellular localization,
Nucleic Acids Research, 38:W497-W502.


To set-up this stand-alone version make sure of the following things:
1.) Enter the path of this directory (the installation directory)
    in the first line of config.py.
2.) This versions of YLoc works with a 32-bit BLAST version and a
    32-bit PrositeScan version, both compiled for Linux. The paths
	in config.py are thus set to the directory Blast and Prosite.
	If use a different system or want to use different versions of
	this software, indicate the paths to them in the config.py file.
3.) Make sure that you have all the neccessary Python packages.

If you still have problems to setup YLoc, please contact
briese@informatik.uni-tuebingen.de.

This code is not ment to be distributed!!!
