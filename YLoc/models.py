from aasequences import *;
from featurecreator import *;
import copy;
import config;
import string;

webloc_abbreviations = { "cytoplasm" : "Cy", "ER": "ER", "extracellular space" : "Ex", "Golgi apparatus" : "Go", "lysosome" : "Ly",
        "mitochondrion" : "Mi", "nucleus" : "Nu", "peroxisome" : "Pe", "plasma membrane" : "Pm" , "secreted pathway" : "SP", "chloroplast" : "Ch", "vacuole" : "Va", "transmembrane" : "TM", "secreted vesicle compartment" : "Ve"};

class Models(object):

	PATH_MODELS = config.PATH_MODELS;
	PATH_DATA = config.PATH_DATA;


	__locations = { "cyt" : "cytoplasm", "mTP" : "mitochondrion", "nuc" : "nucleus", "SR" : "secreted pathway",
	"ER" : "ER", "ex" : "extracellular space", "gol" : "Golgi apparatus", "lys" : "lysosome", "per" : "peroxisome",
	"mem" : "plasma membrane", "vac" : "vacuole", "cTP" : "chloroplast", "cyt_nuc" : "cytoplasm & nucleus",
	"ex_mem" : "extracellular space & plasma membrane", "cyt_mem" : "cytoplasm & plasma membrane",
	"cyt_mTP" : "cytoplasm & mitochondrion", "nuc_mTP" : "nucleus & mitochondrion", "ER_ex" : "ER & extracellular space",
	"ex_nuc" : "nucleus & extracellular space", "tmem" : "transmembrane", "mTP_cTP" : "mitochondrion & chloroplast", "ves" : "secreted vesicle compartment"};


	__model_names = [ "bacello_animals_p+p" ,
					"bacello_animals" ,
					"bacello_fungi",
					"bacello_fungi_p+p",
					"bacello_plants" ,
					"bacello_plants_p+p" ,
					"multiloc_animals",
					"multiloc_fungi" ,
					"multiloc_plants",
					"multiloc_animals_p+p",
					"multiloc_fungi_p+p" ,
					"multiloc_plants_p+p",
					"multiloc_animals_dbm",
					"multiloc_fungi_dbm" ,
					"multiloc_plants_dbm",
					"multiloc_animals_dbm_p+p",
					"multiloc_fungi_dbm_p+p" ,
					"multiloc_plants_dbm_p+p"];
	__external_model_names = ["YLoc-LowRes* Animals","YLoc-LowRes Animals","YLoc-LowRes Fungi","YLoc-LowRes* Fungi","YLoc-LowRes Plants","YLoc-LowRes* Plants","YLoc-HighRes Animals","YLoc-HighRes Fungi","YLoc-HighRes Plants","YLoc-HighRes* Animals","YLoc-HighRes* Fungi","YLoc-HighRes* Plants","YLoc+ Animals","YLoc+ Fungi","YLoc+ Plants","YLoc+* Animals","YLoc+* Fungi","YLoc+* Plants"];

	"Class that saves models fasta files and feature lists"

	def __init__(self, id = "1"):
		self.process_id = id;

	def __createFeaturesForSequences(self, model_name,sequences, response = '?',return_most_similar=False):
		features = FeatureCreator(copy.copy(sequences), id = self.process_id);
		classes = [];
		most_similar = [];

		if model_name == "bacello_animals_p+p":
			features.addSize();
			features.addMonoNLS26SortingSignal();
			features.addMIP13SortingSignal();
			features.addAutoCorrelation(1,20,3,4,features.PROPERTIES[1]);
			features.addAutoCorrelation(1,30,5,6,features.PROPERTIES[2]);
			features.addAutoCorrelation(1,-1,1,2,features.PROPERTIES[2]);
			features.addCountLetterNormed(1, 20, "W" , features.ALPHABETS[0]); #neu
			features.addPseudoComposition(1,30,3,4,"leucine");
			features.addCountLetter(1, 30, "C" , features.ALPHABETS[0]);
			features.addMaxPseudoCompositionNormed(1,60,"alanine");
			features.addPseudoComposition(1,-1,2,3,"cysteine");
			features.addMaxPseudoCompositionNormed(1,90,"hydrophob_very");
			features.addCountLetterNormed(1, 20, "1" , features.ALPHABETS[3]);
			features.addMinPseudoComposition(-80,-1, "charge_none");
			features.addMaxPseudoComposition(-100,-1,"structure_none");
			features.addMaxPseudoComposition(1,-1,"structure_sheet");
			features.addPrositeMotifs(["clust_cyt","clust_SR","clust_mTP","clust_nuc"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "bacello_animals":
			features.addSize();
			features.addMonoNLS26SortingSignal();
			features.addMIP13SortingSignal();
			features.addAutoCorrelation(1,20,3,4,features.PROPERTIES[1]);
			features.addAutoCorrelation(1,30,5,6,features.PROPERTIES[2]);
			features.addAutoCorrelation(1,-1,1,2,features.PROPERTIES[2]);
			features.addPseudoComposition(1,30,3,4,"leucine");
			features.addCountLetter(1, 30, "C" , features.ALPHABETS[0]);
			features.addMaxPseudoCompositionNormed(-30,-1,"cysteine");
			features.addMaxPseudoCompositionNormed(1,90,"hydrophob_very");
			features.addCountLetterNormed(1, 20, "1" , features.ALPHABETS[3]);  #not zero as in file
			features.addMaxPseudoComposition(-100,-1,"structure_none"); #neu
			features.addMaxPseudoComposition(1,-1,"structure_sheet");
			features.addPrositeMotifs(["clust_cyt","clust_mem"]);
			most_similar = features.addGoTerms(["0005739","0005634","0016740","0005576","go_one"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "bacello_fungi":
			features.addSize();
			features.addMonoNLS27SortingSignal();
			features.addVacuole9SortingSignal();
			features.addAutoCorrelation(1,20,3,4,features.PROPERTIES[2]);
			features.addMaxAutoCorrelation(1,110, features.PROPERTIES[2]);
			features.addCountLetter(1, 20, "L" , features.ALPHABETS[0]);
			features.addCountLetter(1, 180, "A" , features.ALPHABETS[0]);  #neu
			features.addMaxPseudoComposition(1,70,"hydrophob_slightly");
			features.addPseudoComposition(-100,-1,4,5,"hydrophob_slightly");
			features.addPseudoComposition(1,20,5,6,"unpolar");
			features.addMaxPseudoCompositionNormed(1,-1,"polar_polar");
			features.addCountLetterNormed(1, 20, "1" , "charge_negative");   # not 0 = no charge as in file
			features.addMaxPseudoCompositionNormed(1,10,"aliphatic_helix_other");
			features.addCountLetter(1, 50, "1" , "aliphatic_helix_uncharged");
			features.addPseudoComposition(1,-1,4,5,"aliphatic_helix_other");
			features.addPrositeMotifs(["clust_cyt","clust_nuc"]);
			most_similar = features.addGoTerms(["0005739","0005576","go_one"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "bacello_fungi_p+p":
			features.addSize();
			features.addMonoNLS27SortingSignal();
			features.addVacuole9SortingSignal(); #neu
			features.addVacuole31SortingSignal(); #neu
			features.addAutoCorrelation(1,20,3,4,features.PROPERTIES[2]);
			features.addMaxAutoCorrelation(1,110, features.PROPERTIES[2]);
			features.addCountLetter(1, 20, "L" , features.ALPHABETS[0]);
			features.addCountLetter(1, 40, "A" , features.ALPHABETS[0]);
			features.addCountLetter(1, 180, "G" , features.ALPHABETS[0]);
			features.addCountLetter(-90, -1, "G" , features.ALPHABETS[0]);
			features.addMaxPseudoComposition(1,70,"hydrophob_slightly");
			features.addPseudoComposition(1,20,5,6,"unpolar");
			features.addMaxPseudoCompositionNormed(1,-1,"polar_polar");
			features.addCountLetterNormed(1, 20, "1" , "charge_negative");   # not 0 = no charge as in file
			features.addPseudoComposition(1,-1,4,5,"size_large");
			features.addMaxPseudoCompositionNormed(1,10,"aliphatic_helix_other");
			features.addCountLetter(1, 50, "1" , "aliphatic_helix_uncharged");
			features.addPseudoComposition(1,-1,4,5,"aliphatic_helix_other");
			features.addPrositeMotifs(["clust_cyt","clust_nuc"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "bacello_plants":
			features.addSize();
			features.addMonoNLS26SortingSignal();
			features.addVacuole37SortingSignal();
			features.addMIP2SortingSignal();
			features.addMIP13SortingSignal();
			features.addAliphaticHelixNormedSortingSignal();
			features.addAutoCorrelation(1,60,4,5,features.PROPERTIES[2]);
			features.addCountLetterNormed(1, 10, "0" , features.ALPHABETS[1]); #not 2 as in file  #neu
			features.addCountLetter(1, 50, "1" , features.ALPHABETS[3]);  #not 0 as written in file
			features.addCountLetter(1, 170, "1" , features.ALPHABETS[3]); #not 0 as written in file
			features.addMaxPseudoCompositionNormed(-100,-1,"charge_none");
			features.addPseudoComposition(1,-1,4,5,"charge_positive");
			features.addMaxPseudoComposition(1,160,"structure_sheet");  #neu
			features.addPseudoComposition(1,140,4,5,"hydroxylated_none");
			features.addPseudoComposition(1,-1,4,5,"aliphatic_helix_other");
			features.addPrositeMotifs(["clust_nuc"]);
			most_similar = features.addGoTerms(["0005739","0005576","0009507","go_one"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "bacello_plants_p+p":
			features.addSize();
			features.addMonoNLS26SortingSignal();
			features.addVacuole37SortingSignal();
			features.addMIP2SortingSignal();
			features.addMIP13SortingSignal();
			features.addAliphaticHelixNormedSortingSignal();
			features.addAutoCorrelation(1,60,4,5,features.PROPERTIES[2]);
			features.addCountLetter(-70, -1, "C" , features.ALPHABETS[0]);
			features.addCountLetterNormed(1, 10, "0" , features.ALPHABETS[1]); # not 2 as written in file
			features.addCountLetter(1, 50, "1" , features.ALPHABETS[3]);  #not 0 as written in file
			features.addCountLetter(1, 170, "1" , features.ALPHABETS[3]); #not 0 as written in file
			features.addMaxPseudoCompositionNormed(-100,-1,"charge_none");
			features.addPseudoComposition(1,-1,4,5,"charge_positive");
			features.addMaxPseudoComposition(1,160,"structure_sheet"); #neu
			features.addPseudoComposition(1,-1,6,7,"size_large");
			features.addPseudoComposition(1,140,4,5,"hydroxylated_none");
			features.addMaxPseudoComposition(-90,-1,"aliphatic_helix_basic");  #neu
			features.addPseudoComposition(1,-1,4,5,"aliphatic_helix_other");
			features.addPrositeMotifs(["clust_cyt"]); #neu
			features.addPrositeMotifs(["clust_nuc"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "multiloc_animals":
			features.addSize();
			features.addMonoNLS26SortingSignal();
			features.addER5SortingSignal();
			features.addVacuole31SortingSignal();
			features.addGolgi4SortingSignal();
			features.addTransmembrane2SortingSignal();
			features.addAutoCorrelation(1,20,1,2,features.PROPERTIES[1]);
			features.addAutoCorrelation(1,70,6,7,features.PROPERTIES[1]);
			features.addProperties(1,30, features.PROPERTIES[2]);
			features.addAutoCorrelation(1,-1,2,3,features.PROPERTIES[2]);
			features.addCountLetter(1, 20, "K" , features.ALPHABETS[0]);
			features.addMaxPseudoCompositionNormed(1,50,"leucine");
			features.addMaxPseudoComposition(1,170,"cysteine");
			features.addPseudoComposition(1,50,2,3,"hydrophob_very");
			features.addPseudoComposition(1,-1,2,3,"hydrophob_very");
			features.addCountLetter(1, 20, "1" , features.ALPHABETS[3]);  ## nicht 0 wie in datei vermutet
			features.addPseudoComposition(1,-1,2,3,"charge_none");
			features.addPseudoComposition(1,-1,6,7,"size_small");
			features.addMaxPseudoComposition(1,-1,"aromatic_none"); #neu
			features.addMinPseudoComposition(50,-20, "aliphatic_helix_uncharged");
			features.addPrositeMotifs(["clust_cyt","clust_nuc","clust_mem"]);
			most_similar = features.addGoTerms(["0005737","0005783","0005739","0005576","0005777","0005764","go_one"]); #neu5764
			classes = self.__getModelClasses(model_name);

		if model_name == "multiloc_animals_p+p":
			features.addSize();
			features.addMonoNLS26SortingSignal();
			features.addER5SortingSignal();
			features.addPeroxi3SortingSignal(); #neu
			features.addVacuole31SortingSignal();
			features.addGolgi4SortingSignal();
			features.addTransmembrane2SortingSignal();
			features.addAutoCorrelation(1,20,1,2,features.PROPERTIES[1]);
			features.addAutoCorrelation(1,70,6,7,features.PROPERTIES[1]);
			features.addProperties(1,30, features.PROPERTIES[2]);
			features.addAutoCorrelation(1,-1,2,3,features.PROPERTIES[2]);
			features.addCountLetter(1, 20, "K" , features.ALPHABETS[0]);
			features.addMaxPseudoComposition(1,30,"asparagine");    #neu
			features.addMaxPseudoCompositionNormed(1,50,"leucine");
			features.addCountLetter(1, 90, "W" , features.ALPHABETS[0]);
			features.addMaxPseudoComposition(1,170,"cysteine");
		        features.addPseudoComposition(1,50,2,3,"hydrophob_very");
			features.addMaxPseudoComposition(-40,-1,"hydrophob_very");
 			features.addPseudoComposition(1,-1,2,3,"hydrophob_very");
  			features.addCountLetter(1, 20, "1" , features.ALPHABETS[3]);
  			features.addPseudoComposition(1,-1,2,3,"charge_none");
			features.addPseudoComposition(1,-1,6,7,"size_small");
			features.addMaxPseudoComposition(1,-1,"aromatic_none");
			features.addMaxPseudoCompositionNormed(1,-1,"hydroxylated_hydroxy");
 			features.addMinPseudoComposition(50,-20, "aliphatic_helix_uncharged");
			features.addPrositeMotifs(["clust_cyt","clust_SR","clust_mTP","clust_nuc","clust_mem"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "multiloc_fungi":
			features.addSize();
			features.addMonoNLS26SortingSignal();
			features.addER5SortingSignal();
			features.addVacuole31SortingSignal();
			features.addGolgi4SortingSignal();
			features.addTransmembrane2SortingSignal();
			features.addMIP13SortingSignal();
			features.addAutoCorrelation(1,20,2,3,features.PROPERTIES[1]);
			features.addAutoCorrelation(1,90,1,2,features.PROPERTIES[1]);
			features.addProperties(1,30, features.PROPERTIES[2]);
			features.addAutoCorrelation(1,-1,2,3,features.PROPERTIES[2]);
			features.addMaxPseudoCompositionNormed(1,40,"leucine");
			features.addMaxPseudoComposition(1,40,"lysine"); #neu
			features.addMaxPseudoComposition(1,170,"cysteine");
			features.addMaxPseudoComposition(1,-1,"serine"); #neu
			features.addPseudoComposition(1,60,2,3,"hydrophob_very");
			features.addMaxPseudoCompositionNormed(1,-1,"hydrophob_very");
			features.addCountLetter(1, 20, "1" , features.ALPHABETS[3]);  ## nicht 0 wie in datei vermutet
			features.addPseudoComposition(1,-1,2,3,"charge_none");
			features.addMinPseudoComposition(1,-1, "structure_sheet");
			features.addMaxPseudoCompositionNormed(1,170, "aliphatic_helix_other");
			features.addPrositeMotifs(["clust_cyt","clust_nuc","clust_mem"]);
			most_similar = features.addGoTerms(["0005737","0005783","0005739","0005576","0005777","go_one"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "multiloc_fungi_p+p":
			features.addSize();
			features.addMonoNLS26SortingSignal();
			features.addER5SortingSignal();
			features.addPeroxi3SortingSignal(); #neu
			features.addVacuole31SortingSignal();
			features.addGolgi4SortingSignal();
			features.addTransmembrane2SortingSignal();
			features.addMIP13SortingSignal();
			features.addAutoCorrelation(1,20,2,3,features.PROPERTIES[1]);
			features.addAutoCorrelation(1,90,1,2,features.PROPERTIES[1]);
			features.addAutoCorrelation(-30,-1,3,4,features.PROPERTIES[1]); #neu
			features.addProperties(1,30, features.PROPERTIES[2]);
			features.addAutoCorrelation(1,-1,2,3,features.PROPERTIES[2]);
			features.addMaxPseudoCompositionNormed(1,40,"leucine");
			features.addMaxPseudoComposition(1,40, "lysine");
			features.addMaxPseudoComposition(1,170,"cysteine");
			features.addMaxPseudoComposition(1,-1, "serine");
			features.addPseudoComposition(1,60,2,3,"hydrophob_very");
			features.addMaxPseudoCompositionNormed(1,-1,"hydrophob_very");
			features.addCountLetter(1, 20, "1" , features.ALPHABETS[3]);  ## nicht 0 wie in datei vermutet
			features.addPseudoComposition(1,-1,2,3,"charge_none");
			features.addMinPseudoComposition(1,-1, "structure_sheet");
			features.addMaxPseudoCompositionNormed(1,-1,"aromatic_aromatic");
			features.addCountLetter(1,110,"0", features.ALPHABETS[7]);	## nicht 1 wie in datei vermutet
			features.addMaxPseudoCompositionNormed(1,170, "aliphatic_helix_other");
			features.addPrositeMotifs(["clust_cyt","clust_SR","clust_mTP","clust_nuc","clust_mem"]);
			classes = self.__getModelClasses(model_name);


		if model_name == "multiloc_plants":
			features.addSize();
			features.addMonoNLS26SortingSignal();
			features.addER5SortingSignal();
			features.addVacuole31SortingSignal();
			features.addGolgi4SortingSignal();
			features.addGolgi6SortingSignal();
			features.addTransmembrane2SortingSignal();
			features.addMIP2SortingSignal();
			features.addProperties(1,10, features.PROPERTIES[0]);
			features.addAutoCorrelation(1,20,2,3,features.PROPERTIES[1]);
			features.addAutoCorrelation(1,70,6,7,features.PROPERTIES[1]);
			features.addAutoCorrelation(1,60,6,7,features.PROPERTIES[2]);
			features.addMaxPseudoCompositionNormed(1,40,"leucine");
			features.addMaxPseudoComposition(1,180,"cysteine");
			features.addMaxPseudoComposition(1,40,"hydrophob_slightly"); #neu
			features.addPseudoComposition(1,50,2,3,"hydrophob_very");
			features.addPseudoComposition(1,-1,2,3,"hydrophob_very");
			features.addCountLetter(1, 20, "1" , features.ALPHABETS[3]);  #nich wie 0 in datei
			features.addPseudoComposition(1,-1,2,3,"charge_none");
			features.addCountLetter(1, 30, "2" , "aliphatic_helix_other");  # nich wie 0(basic) in datei
			features.addMaxPseudoCompositionNormed(1,-1, "aliphatic_helix_other");
			features.addPrositeMotifs(["clust_cyt","clust_nuc","clust_mem"]);
			most_similar = features.addGoTerms(["0005783","0005739","0005576","0005777","0009507","go_one"]);#neu 0005783
			classes = self.__getModelClasses(model_name);

		if model_name == "multiloc_plants_p+p":
			features.addSize();
			features.addMonoNLS26SortingSignal();
			features.addER5SortingSignal();
			features.addVacuole31SortingSignal();
			features.addGolgi4SortingSignal();
			features.addGolgi6SortingSignal();
			features.addTransmembrane2SortingSignal();
			features.addMIP2SortingSignal();
			features.addProperties(1,10, features.PROPERTIES[0]);
			features.addAutoCorrelation(1,20,2,3,features.PROPERTIES[1]);
			features.addAutoCorrelation(1,70,6,7,features.PROPERTIES[1]);
			features.addAutoCorrelation(1,60,6,7,features.PROPERTIES[2]);
			features.addCountLetter(1, 20, "K" , features.ALPHABETS[0]);
			features.addMaxPseudoCompositionNormed(1,40,"leucine");
			features.addMaxPseudoComposition(1,180,"cysteine");
			features.addMaxPseudoComposition(1,40,"hydrophob_slightly");
			features.addPseudoComposition(1,50,2,3,"hydrophob_very");
			features.addPseudoComposition(1,-1,2,3,"hydrophob_very");
			features.addCountLetter(1, 20, "1" , features.ALPHABETS[3]);  #nich wie 0 in datei
			features.addPseudoComposition(1,-1,2,3,"charge_none");
			features.addCountLetter(1, 160, "2" , features.ALPHABETS[4]); #nicht 0 wie in datei
			features.addMaxPseudoComposition(1,-1,"aromatic_none"); #neu
			features.addCountLetter(1, 110, "1" , features.ALPHABETS[7]); #neu, nich 0 wie in file
			features.addCountLetter(1, 30, "2" , "aliphatic_helix_other");  # nich wie 0(basic) in datei
			features.addMaxPseudoCompositionNormed(1,-1, "aliphatic_helix_other");
			features.addPrositeMotifs(["clust_cyt","clust_SR","clust_mTP","clust_nuc","clust_mem"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "multiloc_animals_dbm":
			features.addSize();
			features.addER5SortingSignal();
			features.addPeroxi1SortingSignal();
			features.addVacuole21SortingSignal();
			features.addVacuole31SortingSignal();
			features.addTransmembrane2SortingSignal();
			features.addMIP12SortingSignal();
			features.addAutoCorrelation(1,20,2,3,features.PROPERTIES[1]);
			features.addMaxAutoCorrelation(1,30, features.PROPERTIES[2]);
			features.addCountLetter(1, 70, "M" , features.ALPHABETS[0]); #neu
			features.addCountLetter(1, 70, "N" , features.ALPHABETS[0]); #neu
			features.addMaxPseudoComposition(1,120,"cysteine");
			features.addMaxPseudoComposition(1,120,"lysine");
			features.addCountLetter(1, 120, "W" , features.ALPHABETS[0]);
			features.addPseudoComposition(1,130,2,3,"hydrophob_very");
			features.addMaxPseudoCompositionNormed(-40,-1,"hydrophob_very");
			features.addMaxPseudoComposition(1,30,"charge_negative");
			features.addMinPseudoComposition(1,-1,"structure_helix");
			features.addCountLetter(1, 20, "1" , "size_small");  # not 0 (tiny) as in file
			features.addMaxPseudoCompositionNormed(1,-1,"aromatic_aromatic");
			features.addPseudoComposition(-100,-1,3,4,"hydroxylated_hydroxy"); #neu
			features.addPseudoComposition(1,-1,2,3,"aliphatic_helix_uncharged");
			features.addPrositeMotifs(["clust_mem","clust_nuc"]);
			most_similar = features.addGoTerms(["0005783","0005739","0005576","0042025","0005778","go_one"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "multiloc_animals_dbm_p+p":
			features.addSize();
			features.addER5SortingSignal();
			features.addPeroxi1SortingSignal();
			features.addVacuole21SortingSignal();
			features.addVacuole31SortingSignal();
			features.addTransmembrane2SortingSignal();
			features.addMIP12SortingSignal();
			features.addAutoCorrelation(1,20,2,3,features.PROPERTIES[1]);
			features.addMaxAutoCorrelation(1,30, features.PROPERTIES[2]);
			features.addCountLetter(1, 70, "N" , features.ALPHABETS[0]); #neu
			features.addPseudoComposition(1,90,2,3,"serine");
			features.addMaxPseudoComposition(1,120,"cysteine");
			features.addMaxPseudoComposition(1,120,"lysine");
			features.addCountLetter(1, 120, "M" , features.ALPHABETS[0]);
			features.addCountLetter(1, 120, "W" , features.ALPHABETS[0]);
			features.addPseudoComposition(1,130,2,3,"hydrophob_very");
			features.addMaxPseudoCompositionNormed(-40,-1,"hydrophob_very");
			features.addMaxPseudoComposition(1,30,"charge_negative");
			features.addMinPseudoComposition(1,-1,"structure_helix");
			features.addCountLetter(1, 20, "1" , "size_small");  # not 0 (tiny) as in file
			features.addMaxPseudoCompositionNormed(1,-1,"aromatic_aromatic");
			features.addPseudoComposition(-100,-1,3,4,"hydroxylated_hydroxy"); #neu
			features.addPseudoComposition(1,-1,2,3,"aliphatic_helix_uncharged");
			features.addPrositeMotifs(["PS00639","PS50041","clust_SR","clust_mem","clust_nuc","clust_mTP","clust_cyt"]); #PS00639 neu
			classes = self.__getModelClasses(model_name);

		if model_name == "multiloc_fungi_dbm":
			features.addSize();
			features.addER5SortingSignal();
			features.addPeroxi1SortingSignal();
			features.addVacuole21SortingSignal();
			features.addVacuole31SortingSignal();
			features.addTransmembrane2SortingSignal();
			features.addMIP12SortingSignal();
			features.addAutoCorrelation(1,20,2,3,features.PROPERTIES[1]);
			features.addMaxAutoCorrelation(1,30, features.PROPERTIES[2]);
			features.addCountLetter(1, 70, "N" , features.ALPHABETS[0]); #neu
			features.addMaxPseudoComposition(1,120,"lysine");
			features.addCountLetter(1, 120, "W" , features.ALPHABETS[0]);
			features.addMaxPseudoCompositionNormed(1,130,"phenylalanine");
			features.addMaxPseudoComposition(1,130,"cysteine");
			features.addPseudoComposition(1,130,2,3,"hydrophob_very");
			features.addMaxPseudoCompositionNormed(-40,-1,"hydrophob_very");
			features.addMaxPseudoComposition(1,30,"charge_negative");
			features.addMinPseudoComposition(1,-1,"structure_helix");
			features.addCountLetter(1, 20, "1" , "size_small");  # not 0 (tiny) as in file
			features.addPseudoComposition(-100,-1,3,4,"hydroxylated_hydroxy"); #neu
			features.addPseudoComposition(1,-1,2,3,"aliphatic_helix_uncharged");
			features.addPrositeMotifs(["clust_mem","clust_nuc"]);
			most_similar = features.addGoTerms(["0005783","0005739","0005576","0042025","0005773","0005778","go_one"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "multiloc_fungi_dbm_p+p":
			features.addSize();
			features.addER5SortingSignal();
			features.addPeroxi1SortingSignal();
			features.addVacuole21SortingSignal();
			features.addVacuole31SortingSignal();
			features.addTransmembrane2SortingSignal();
			features.addMIP12SortingSignal();
			features.addAutoCorrelation(1,20,2,3,features.PROPERTIES[1]);
			features.addMaxAutoCorrelation(1,30, features.PROPERTIES[2]);
			features.addCountLetter(1, 70, "N" , features.ALPHABETS[0]); #neu
			features.addMaxPseudoCompositionNormed(1,120,"tyrosine"); #neu
			features.addMaxPseudoComposition(1,120,"lysine");
			features.addCountLetter(1, 120, "M" , features.ALPHABETS[0]);
			features.addCountLetter(1, 120, "W" , features.ALPHABETS[0]);
			features.addMaxPseudoCompositionNormed(1,130,"phenylalanine");
			features.addMaxPseudoComposition(1,130,"cysteine");
			features.addPseudoComposition(1,130,2,3,"hydrophob_very");
			features.addMaxPseudoCompositionNormed(-40,-1,"hydrophob_very");
			features.addMaxPseudoComposition(1,30,"charge_negative");
			features.addMinPseudoComposition(1,-1,"structure_helix");
			features.addCountLetter(1, 20, "1" , "size_small");  # not 0 (tiny) as in file
			features.addMaxPseudoComposition(-100,-1,"aromatic_aromatic");
			features.addPseudoComposition(-100,-1,3,4,"hydroxylated_hydroxy"); #neu
			features.addPseudoComposition(1,-1,2,3,"aliphatic_helix_uncharged");
			features.addPrositeMotifs(["PS50041","clust_SR","clust_mem","clust_nuc","clust_mTP","clust_cyt"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "multiloc_plants_dbm":
			features.addSize();
			features.addER5SortingSignal();
			features.addVacuole31SortingSignal();
			features.addTransmembrane2SortingSignal();
			features.addMIP12SortingSignal();
			features.addAutoCorrelation(1,20,1,2,features.PROPERTIES[1]);
			features.addProperties(1,10, features.PROPERTIES[2]);
			features.addAutoCorrelation(1,50,6,7,features.PROPERTIES[2]);
			features.addCountLetter(1, 50, "S" , features.ALPHABETS[0]);
			features.addMaxPseudoCompositionNormed(1,80,"leucine");
			features.addMaxPseudoComposition(1,120,"cysteine");
			features.addMaxPseudoComposition(1,120,"lysine");
			features.addCountLetter(1, 120, "F" , features.ALPHABETS[0]); #neu
			features.addCountLetter(1, 120, "W" , features.ALPHABETS[0]);
			features.addPseudoComposition(1,100,2,3,"hydrophob_very");
			features.addMaxPseudoCompositionNormed(-40,-1,"hydrophob_very");
			features.addMaxPseudoComposition(1,30,"charge_negative");
			features.addPseudoComposition(1,20,3,4,"structure_helix");
			features.addMinPseudoComposition(1,-1,"structure_helix");
			features.addCountLetter(1, 20, "2" , "size_large");
			features.addMaxPseudoComposition(1,-1,"hydroxylated_hydroxy"); #neu
			features.addPseudoComposition(1,-1,2,3,"aliphatic_helix_uncharged");
			features.addPrositeMotifs(["clust_mem","clust_nuc"]);
			most_similar = features.addGoTerms(["0005783","0005739","0005576","0042025","0009507","go_one"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "multiloc_plants_dbm_p+p":
			features.addSize();
			features.addER5SortingSignal();
			features.addPeroxi1SortingSignal(); #neu
			features.addVacuole31SortingSignal();
			features.addTransmembrane2SortingSignal();
			features.addMIP12SortingSignal();
			features.addChloro2SortingSignal();
			features.addAutoCorrelation(1,20,1,2,features.PROPERTIES[1]);
			features.addProperties(1,10, features.PROPERTIES[2]);
			features.addAutoCorrelation(1,50,6,7,features.PROPERTIES[2]);
			features.addCountLetter(1, 50, "S" , features.ALPHABETS[0]);
			features.addMaxPseudoCompositionNormed(1,80,"leucine");
			features.addMaxPseudoComposition(1,120,"cysteine");
			features.addMaxPseudoComposition(1,120,"lysine");
			features.addCountLetter(1, 120, "F" , features.ALPHABETS[0]);
			features.addCountLetter(1, 120, "W" , features.ALPHABETS[0]);
			features.addPseudoComposition(1,100,2,3,"hydrophob_very");
			features.addMaxPseudoCompositionNormed(-40,-1,"hydrophob_very");
			features.addMaxPseudoComposition(1,30,"charge_negative");
			features.addPseudoComposition(1,20,3,4,"structure_helix");
			features.addMinPseudoComposition(1,-1,"structure_helix");
			features.addCountLetter(1, 20, "2" , "size_large");
			features.addMaxPseudoComposition(1,-1,"hydroxylated_hydroxy"); #neu
			features.addPseudoComposition(1,-1,2,3,"aliphatic_helix_uncharged");
			features.addPrositeMotifs(["clust_SR","clust_mem","clust_nuc","clust_mTP","clust_cTP","clust_cyt"]);
			classes = self.__getModelClasses(model_name);

		if model_name == "suba_plus":
			features.addSize();
			features.addAutoCorrelation(1,20,1,2,features.PROPERTIES[1]);
			features.addAutoCorrelation(1,-1,1,2,features.PROPERTIES[1]);
			features.addAutoCorrelation(1,120,1,2,features.PROPERTIES[2]);
			features.addCountLetter(1, 10, "A" , features.ALPHABETS[0]);
			features.addCountLetter(1, 10, "K" , features.ALPHABETS[0]);
			features.addCountLetter(1, 20, "G" , features.ALPHABETS[0]);
			features.addCountLetter(1, 70, "Y" , features.ALPHABETS[0]);
			features.addMaxPseudoCompositionNormed(1,80,"tryptophane");
			features.addMaxPseudoComposition(1,80,"glutamine");
			features.addMaxPseudoCompositionNormed(1,190,"serine");
			features.addCountLetter(1, 190, "C" , features.ALPHABETS[0]);
			features.addCountLetter(-80,-1, "S" , features.ALPHABETS[0]);
			features.addPseudoComposition(-70,-1,4,5,"proline");
			features.addMaxPseudoComposition(-60,-1,"leucine");
			features.addCountLetter(-30, -1, "I" , features.ALPHABETS[0]);
			features.addMaxPseudoCompositionNormed(-20,-1,"glutamate");
			features.addCountLetter(1,40, "1" , "charge_negative");
			features.addPseudoComposition(1,-1,3,4,"structure_helix");
			features.addMaxPseudoComposition(1,40,"size_small");
			features.addPseudoComposition(1,-1,3,4,"size_large");
			features.addPhylogeneticProfiles([14]);
			features.addPhylogeneticProfiles([77]);
			features.addPrositeMotifs(["PS50191","PS50920","clust_nuc"]);
			most_similar = features.addGoTerms(["0005789","0005576","0004298","go_one"]);
			classes = self.__getModelClasses(model_name);

		if response == '?' :
			features.addResponse("?");
			features.setResponseClasses(classes);
		else:
			features.addResponse(response);

		if return_most_similar:
			return (features, most_similar);
		else:
			return features;

	def createFeaturesForSequences(self, sequences, model_name):
		return self.__createFeaturesForSequences(model_name,sequences,return_most_similar=True)

	def createFeaturesForDataSet(self, data_set, model_name):
		sequences = AASequences();
		allfeatures = FeatureCreator();

		if data_set == "decoys":
			ff_name = self.PATH_DATA + "Decoys/decoys.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			allfeatures = self.__createFeaturesForSequences(model_name,sequences, "X");

			return allfeatures;

		if data_set == "suba_plus":

			ff_name = self.PATH_DATA + "SUBA_Plus/nucleus.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			allfeatures = self.__createFeaturesForSequences(model_name,sequences, "nuc");

			ff_name = self.PATH_DATA + "SUBA_Plus/plastid.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "cTP");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "SUBA_Plus/plasma_membrane.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "mem");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "SUBA_Plus/transmembrane.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "tmem");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "SUBA_Plus/mitochondrion.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "mTP");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "SUBA_Plus/extracellular.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "ex");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "SUBA_Plus/cytoplasm.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "cyt");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "SUBA_Plus/peroxisome.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "per");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "SUBA_Plus/endoplasmatic_reticulum.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "ves");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "SUBA_Plus/vacuole.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "ves");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "SUBA_Plus/golgi.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "ves");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "SUBA_Plus/mitochondrion_plastid.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "mTP_cTP");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "SUBA_Plus/cytoplasm_nucleus.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "cyt_nuc");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "SUBA_Plus/cytoplasm_plasma_membrane.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "cyt_mem");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "SUBA_Plus/extracellular_plasma_membrane.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "ex_mem");
			allfeatures = allfeatures +  features;

			return allfeatures;

		if data_set == "bacello_animals" or data_set == "bacello_animals_p+p":

			ff_name = self.PATH_DATA + "BaCelLo/Animals/Cytoplasm.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			allfeatures = self.__createFeaturesForSequences(model_name,sequences, "cyt");

			ff_name = self.PATH_DATA + "BaCelLo/Animals/Mitochondion.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "mTP");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "BaCelLo/Animals/Nucleus.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name, sequences,  "nuc");
			allfeatures = allfeatures + features;

			ff_name = self.PATH_DATA + "BaCelLo/Animals/Secretory.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name, sequences, "SR");
			allfeatures = allfeatures + features;

			return allfeatures;

		if data_set == "bacello_fungi" or data_set == "bacello_fungi_p+p":

			ff_name = self.PATH_DATA + "BaCelLo/Fungi/Cytoplasm.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			allfeatures = self.__createFeaturesForSequences(model_name,sequences, "cyt");

			ff_name = self.PATH_DATA + "BaCelLo/Fungi/Mitochondion.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "mTP");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "BaCelLo/Fungi/Nucleus.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name, sequences,  "nuc");
			allfeatures = allfeatures + features;

			ff_name = self.PATH_DATA + "BaCelLo/Fungi/Secretory.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name, sequences, "SR");
			allfeatures = allfeatures + features;

			return allfeatures;

		if data_set == "bacello_plants" or data_set == "bacello_plants_p+p":

			ff_name = self.PATH_DATA + "BaCelLo/Plants/Cytoplasm.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			allfeatures = self.__createFeaturesForSequences(model_name,sequences, "cyt");

			ff_name = self.PATH_DATA + "BaCelLo/Plants/Mitochondion.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "mTP");
			allfeatures = allfeatures +  features;

			ff_name = self.PATH_DATA + "BaCelLo/Plants/Nucleus.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name, sequences,  "nuc");
			allfeatures = allfeatures + features;

			ff_name = self.PATH_DATA + "BaCelLo/Plants/Secretory.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name, sequences, "SR");
			allfeatures = allfeatures + features;

			ff_name = self.PATH_DATA + "BaCelLo/Plants/Chloroplast.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name, sequences, "cTP");
			allfeatures = allfeatures + features;

			return allfeatures;

		split = data_set.split('_');

		if split[0] == "multiloc":

			ff_name = self.PATH_DATA + "MultiLocDataset/cytoplasmic.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			allfeatures = self.__createFeaturesForSequences(model_name,sequences, "cyt");

			if split[1] == "plants":
				ff_name = self.PATH_DATA + "MultiLocDataset/chloroplast.fasta";
				sequences.clear();
				sequences.read_fasta_file(ff_name);
				features = self.__createFeaturesForSequences(model_name,sequences, "cTP");
				allfeatures = allfeatures + features;

			ff_name = self.PATH_DATA + "MultiLocDataset/ER.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "ER");
			allfeatures = allfeatures + features;

			ff_name = self.PATH_DATA + "MultiLocDataset/extracellular.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "ex");
			allfeatures = allfeatures + features;

			ff_name = self.PATH_DATA + "MultiLocDataset/Golgi.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "gol");
			allfeatures = allfeatures + features;

			if split[1] == "animals":
				ff_name = self.PATH_DATA + "MultiLocDataset/lysosomal.fasta";
				sequences.clear();
				sequences.read_fasta_file(ff_name);
				features = self.__createFeaturesForSequences(model_name,sequences, "lys");
				allfeatures = allfeatures + features;

			ff_name = self.PATH_DATA + "MultiLocDataset/mitochondrial.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "mTP");
			allfeatures = allfeatures + features;

			ff_name = self.PATH_DATA + "MultiLocDataset/nuclear.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "nuc");
			allfeatures = allfeatures + features;

			ff_name = self.PATH_DATA + "MultiLocDataset/peroxisomal.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "per");
			allfeatures = allfeatures + features;

			ff_name = self.PATH_DATA + "MultiLocDataset/plasma_membrane.fasta";
			sequences.clear();
			sequences.read_fasta_file(ff_name);
			features = self.__createFeaturesForSequences(model_name,sequences, "mem");
			allfeatures = allfeatures + features;

			if split[1] != "animals":
				ff_name = self.PATH_DATA + "MultiLocDataset/vacuolar.fasta";
				sequences.clear();
				sequences.read_fasta_file(ff_name);
				features = self.__createFeaturesForSequences(model_name,sequences, "vac");
				allfeatures = allfeatures + features;


			if len(split) >= 3 and split[2] == "dbm":

				ff_name = self.PATH_DATA + "DBMLocReduced80/cytoplasm_nucleus_.fasta";
				sequences.clear();
				sequences.read_fasta_file(ff_name);
				features = self.__createFeaturesForSequences(model_name,sequences, "cyt_nuc");
				allfeatures =  allfeatures + features;

				ff_name = self.PATH_DATA + "DBMLocReduced80/secreted_plasma_membrane_.fasta";
				sequences.clear();
				sequences.read_fasta_file(ff_name);
				features = self.__createFeaturesForSequences(model_name,sequences, "ex_mem");
				allfeatures =  allfeatures + features;

				ff_name = self.PATH_DATA + "DBMLocReduced80/cytoplasm_plasma_membrane_.fasta";
				sequences.clear();
				sequences.read_fasta_file(ff_name);
				features = self.__createFeaturesForSequences(model_name,sequences, "cyt_mem");
				allfeatures =  allfeatures + features;

				ff_name = self.PATH_DATA + "DBMLocReduced80/cytoplasm_mitochondrion_.fasta";
				sequences.clear();
				sequences.read_fasta_file(ff_name);
				features = self.__createFeaturesForSequences(model_name,sequences, "cyt_mTP");
				allfeatures =  allfeatures + features;

				ff_name = self.PATH_DATA + "DBMLocReduced80/nucleus_mitochondrion_.fasta";
				sequences.clear();
				sequences.read_fasta_file(ff_name);
				features =self.__createFeaturesForSequences(model_name,sequences, "nuc_mTP");
				allfeatures =  allfeatures + features;

				ff_name = self.PATH_DATA + "DBMLocReduced80/secreted_ER_.fasta";
				sequences.clear();
				sequences.read_fasta_file(ff_name);
				features = self.__createFeaturesForSequences(model_name,sequences, "ER_ex");
				allfeatures =  allfeatures + features;

				ff_name = self.PATH_DATA + "DBMLocReduced80/secreted_nucleus_.fasta";
				sequences.clear();
				sequences.read_fasta_file(ff_name);
				features = self.__createFeaturesForSequences(model_name,sequences, "ex_nuc");
				allfeatures =  allfeatures + features;

			return allfeatures;

	def getModelFilePath(self, model_name):
		if model_name == "bacello_animals_p+p":
			return self.PATH_MODELS + "bacello_animals_p+p.model"
		if model_name == "bacello_animals":
			return self.PATH_MODELS + "bacello_animals.model"
		if model_name == "bacello_fungi":
			return self.PATH_MODELS + "bacello_fungi.model"
		if model_name == "bacello_fungi_p+p":
			return self.PATH_MODELS + "bacello_fungi_p+p.model"
		if model_name == "bacello_plants":
			return self.PATH_MODELS + "bacello_plants.model"
		if model_name == "bacello_plants_p+p":
			return self.PATH_MODELS + "bacello_plants_p+p.model"
		if model_name == "multiloc_animals":
			return self.PATH_MODELS + "multiloc_animals.model"
		if model_name == "multiloc_animals_p+p":
			return self.PATH_MODELS + "multiloc_animals_p+p.model"
		if model_name == "multiloc_fungi":
			return self.PATH_MODELS + "multiloc_fungi.model"
		if model_name == "multiloc_fungi_p+p":
			return self.PATH_MODELS + "multiloc_fungi_p+p.model"
		if model_name == "multiloc_plants":
			return self.PATH_MODELS + "multiloc_plants.model"
		if model_name == "multiloc_plants_p+p":
			return self.PATH_MODELS + "multiloc_plants_p+p.model"
		if model_name == "multiloc_animals_dbm":
			return self.PATH_MODELS + "multiloc_animals_dbm.model"
		if model_name == "multiloc_animals_dbm_p+p":
			return self.PATH_MODELS + "multiloc_animals_dbm_p+p.model"
		if model_name == "multiloc_fungi_dbm":
			return self.PATH_MODELS + "multiloc_fungi_dbm.model"
		if model_name == "multiloc_fungi_dbm_p+p":
			return self.PATH_MODELS + "multiloc_fungi_dbm_p+p.model"
		if model_name == "multiloc_plants_dbm":
			return self.PATH_MODELS + "multiloc_plants_dbm.model"
		if model_name == "multiloc_plants_dbm_p+p":
			return self.PATH_MODELS + "multiloc_plants_dbm_p+p.model"
		if model_name == "suba_plus":
			return self.PATH_MODELS + "suba_plus.model"

	def __getModelClasses(self, model_name):
		if model_name == "bacello_animals_p+p":
			return ["cyt","mTP","nuc","SR"];
		if model_name == "bacello_animals":
			return ["cyt","mTP","nuc","SR"];
		if model_name == "bacello_fungi":
			return ["cyt","mTP","nuc","SR"];
		if model_name == "bacello_fungi_p+p":
			return ["cyt","mTP","nuc","SR"];
		if model_name == "bacello_plants":
			return ["cyt","mTP","nuc","SR","cTP"];
		if model_name == "bacello_plants_p+p":
			return ["cyt","mTP","nuc","SR","cTP"];
		if model_name == "multiloc_animals" or model_name == "multiloc_animals_p+p":
			return ["cyt","ER","ex","gol","lys","mTP","nuc","per","mem"];
		if model_name == "multiloc_fungi" or model_name == "multiloc_fungi_p+p":
			return ["cyt","ER","ex","gol","mTP","nuc","per","mem","vac"];
		if model_name == "multiloc_plants" or model_name == "multiloc_plants_p+p":
			return ["cyt","cTP","ER","ex","gol","mTP","nuc","per","mem","vac"];
		if model_name == "multiloc_animals_dbm" or model_name == "multiloc_animals_dbm_p+p":
			return ["cyt","ER","ex","gol","lys","mTP","nuc","per","mem","cyt_nuc","ex_mem","cyt_mem","cyt_mTP","nuc_mTP","ER_ex","ex_nuc"];
		if model_name == "multiloc_fungi_dbm" or model_name == "multiloc_fungi_dbm_p+p":
			return ["cyt","ER","ex","gol","mTP","nuc","per","mem","vac","cyt_nuc","ex_mem","cyt_mem","cyt_mTP","nuc_mTP","ER_ex","ex_nuc"];
		if model_name == "multiloc_plants_dbm" or model_name == "multiloc_plants_dbm_p+p":
			return ["cyt","cTP","ER","ex","gol","mTP","nuc","per","mem","vac","cyt_nuc","ex_mem","cyt_mem","cyt_mTP","nuc_mTP","ER_ex","ex_nuc"];
		if model_name == "suba_plus":
			return ["nuc","cTP","mem","tmem","mTP","ex","cyt","per","ves","mTP_cTP","cyt_nuc","cyt_mem","ex_mem"];

	def __getModelClassesForPrediction(self, model_name):
		if model_name == "bacello_animals_p+p":
			return ["cyt","mTP","nuc","SR"];
		if model_name == "bacello_animals":
			return ["cyt","mTP","nuc","SR"];
		if model_name == "bacello_fungi":
			return ["cyt","mTP","nuc","SR"];
		if model_name == "bacello_fungi_p+p":
			return ["cyt","mTP","nuc","SR"];
		if model_name == "bacello_plants":
			return ["cyt","mTP","nuc","SR","cTP"];
		if model_name == "bacello_plants_p+p":
			return ["cyt","mTP","nuc","SR","cTP"];
		if model_name == "multiloc_animals" or model_name == "multiloc_animals_p+p":
			return ["cyt","ER","ex","gol","lys","mTP","nuc","per","mem"];
		if model_name == "multiloc_fungi" or model_name == "multiloc_fungi_p+p":
			return ["cyt","ER","ex","gol","mTP","nuc","per","mem","vac"];
		if model_name == "multiloc_plants" or model_name == "multiloc_plants_p+p":
			return ["cyt","cTP","ER","ex","gol","mTP","nuc","per","mem","vac"];
		if model_name == "multiloc_animals_dbm" or model_name == "multiloc_animals_dbm_p+p":
			return ["cyt","ER","ex","gol","lys","mTP","nuc","per","mem"];
		if model_name == "multiloc_fungi_dbm" or model_name == "multiloc_fungi_dbm_p+p":
			return ["cyt","ER","ex","gol","mTP","nuc","per","mem","vac"];
		if model_name == "multiloc_plants_dbm" or model_name == "multiloc_plants_dbm_p+p":
			return ["cyt","cTP","ER","ex","gol","mTP","nuc","per","mem","vac"];
		if model_name == "suba_plus":
			return ["nuc","cTP","mem","tmem","mTP","ex","cyt","per","ves"];

	def getModelClasses(self, model_name):
		classes = [];
		classes = self.__getModelClassesForPrediction(model_name);
		for i in range(len(classes)):
			classes[i] = self.__locations[classes[i]];
		return classes;

	def getAvailableModels(self):
		return ["bacello_animals_p+p","bacello_animals","bacello_fungi","bacello_plants","bacello_plants_p+p","multiloc_animals","multiloc_fungi","multiloc_plants","multiloc_animals_p+p","multiloc_fungi_p+p","multiloc_plants_p+p","multiloc_animals_dbm","multiloc_fungi_dbm","multiloc_plants_dbm","multiloc_animals_dbm_p+p", "multiloc_fungi_dbm_p+p" , 	"multiloc_plants_dbm_p+p"];

	def getAvailableModelsOutside(self):
		return self.__external_model_names;

	def getAvailableModelCombinations(self):
		return [["YLoc-LowRes","YLoc-HighRes","YLoc+"],["Animals","Fungi","Plants"],["Yes","No"]];

	def getModelFromCombinations(self, sel1, sel2, sel3):
		model_name = "bacello";
		if (len(sel1) == 5 or sel1[5] == "H"):
			model_name = "multiloc";
		model_name += "_"+sel2.lower();
		if (len(sel1) == 5):
			model_name += "_dbm";
		if (sel3 == "No"):
			model_name += "_p+p";
		return model_name;

	def getGOIntervals(self, model_name, short=True):
		if model_name == "bacello_animals" or model_name == "YLoc-LowRes Animals":
			if short:
				return ["no","cy","ex","mi","nu","pm"];
			else:
				return ["no","cytoplasm","extracellular","mitochondrion","nuclear","plasma membrane"];
		elif model_name == "bacello_fungi" or model_name == "YLoc-LowRes Fungi":
			if short:
				return ["no","cy","ex","mi","nu","pm"];
			else:
				return ["no","cytoplasm","extracellular","mitochondrion","nuclear","plasma membrane"];
		elif model_name == "bacello_plants" or model_name == "YLoc-LowRes Plants":
			if short:
				return ["no","ch","cy","ER","ex","mi","nu","pm/va"];
			else:
				return ["no","chloroplast","cytoplasm","endoplasmatic reticulum","extracellular","mitochondrion","nuclear","plasma membrane or vacuolar"];
		elif model_name == "multiloc_animals" or model_name == "YLoc-HighRes Animals":
			if short:
				return ['none','cy','ER','ex','go','ly','mi','nu','pe','pm'];
			else:
				return ["no","cytoplasm","endoplasmatic reticulum","extracellular","Golgi appartus","lysosome","mitochondrion","nuclear","peroxisomal","plasma membrane"];
		elif model_name == "multiloc_fungi" or model_name == "YLoc-HighRes Fungi":
			if short:
				return ['none','cy','ER','ex','go','mi','nu','pe','pm','va']
			else:
				return ["no","cytoplasm","endoplasmatic reticulum","extracellular","Golgi appartus","mitochondrion","nuclear","peroxisomal","plasma membrane","vacuolar"];
		elif model_name == "multiloc_plants" or model_name == "YLoc-HighRes Plants":
			if short:
				return ['none','ch','cy','ER','ex','go','mi','nu','pe','pm','va']
			else:
				return ["no","chloroplast","cytoplasm","endoplasmatic reticulum","extracellular","Golgi appartus","mitochondrion","nuclear","peroxisomal","plasma membrane","vacuolar"];
		elif model_name == "multiloc_animals_dbm" or model_name == "YLoc+ Animals":
			if short:
				return ['none','ch','cy','ER','ex','go','ly','mi','nu','pe','pm']
			else:
				return ["no","chloroplast","cytoplasm","endoplasmatic reticulum","extracellular","Golgi appartus","lysosome","mitochondrion","nuclear","peroxisomal","plasma membrane"];
		else:
			if short:
				return ['none','ch','cy','ER','ex','go','ly','mi','nu','pe','pm','va']
			else:
				return ["no","chloroplast","cytoplasm","endoplasmatic reticulum","extracellular","Golgi appartus","lysosome","mitochondrion","nuclear","peroxisomal","plasma membrane","vacuolar"];


	def getModelName(self,model_name):
		if model_name in self.__model_names:
			pos = self.__model_names.index(model_name);
			return self.__external_model_names[pos];
		else:
			return "unknown model";

	def getModelNameIntern(self,model_name):
		if model_name in self.__external_model_names:
			pos = self.__external_model_names.index(model_name);
			return self.__model_names[pos];
		else:
			return "unknown model";

	def getDualOptions(self, model_name):
		if model_name == "multiloc_animals_dbm" or model_name == "multiloc_animals_dbm_p+p":
			return [9,[0,6,-1,2,8,-1,0,8,-1,0,5,-1,5,6,-1,1,2,-1,2,6,-1]];
		if model_name == "multiloc_fungi_dbm" or model_name == "multiloc_fungi_dbm_p+p":
			return [9,[0,5,-1,2,7,-1,0,7,-1,0,4,-1,4,5,-1,1,2,-1,2,5,-1]];
		if model_name == "multiloc_plants_dbm" or model_name == "multiloc_plants_dbm_p+p":
			return [10,[0,6,-1,3,8,-1,0,8,-1,0,5,-1,5,6,-1,2,3,-1,3,6,-1]];
		if model_name == "suba_plus":
			return [9,[1,4,-1,0,6,-1,2,6,-1,2,5,-1]];
		return [];

	def getModelFeatureNames(self, model_name):
		names = [];
		if model_name == "bacello_animals_p+p":
			names.append("size");
			names.append("Mono_NLS_Signal26");
			names.append("MIP_signal13");
			names.append("AutoCorrelatuion|d=3|1-20|hydrophob");
			names.append("AutoCorrelatuion|d=5|1-30|charge");
			names.append("AutoCorrelatuion|d=1|1--1|charge");
			names.append("hist_normed|1-20|aminoacid|W");
			names.append("Paacount|d=3|1-30|aminoacid|L");
			names.append("hist|1-30|aminoacid|C");
			names.append("maxPaacountNormed|1-60|aminoacid|A");
			names.append("Paacount|d=2|1--1|aminoacid|C");
			names.append("maxPaacountNormed|1-90|hydrophob|2");
			names.append("hist_normed|1-20|charged|0");
			names.append("minPaacount|-80--1|charged|0");
			names.append("maxPaacount|-100--1|structure|0");
			names.append("maxPaacount|1--1|structure|2");
			names.append("prosite_clust_conf_80_cyt");
			names.append("prosite_clust_conf_80_SR");
			names.append("prosite_clust_conf_80_mTP");
			names.append("prosite_clust_conf_80_nuc");

		if model_name == "bacello_animals":
			names.append("size");
			names.append("Mono_NLS_Signal26");
			names.append("MIP_signal13");
			names.append("AutoCorrelatuion|d=3|1-20|hydrophob");
			names.append("AutoCorrelatuion|d=5|1-30|charge");
			names.append("AutoCorrelatuion|d=1|1--1|charge");
			names.append("Paacount|d=3|1-30|aminoacid|L");
			names.append("hist|1-30|aminoacid|C");
			names.append("maxPaacountNormed|-30--1|aminoacid|C");
			names.append("maxPaacountNormed|1-90|hydrophob|2");
			names.append("hist_normed|1-20|charged|1");
			names.append("maxPaacount|-100--1|structure|0");
			names.append("maxPaacount|-1--1|structure|2");
			names.append("prosite_clust_cyt");
			names.append("prosite_clust_mem");
			names.append("GO:0005739");
			names.append("GO:0005634");
			names.append("GO:0016740");
			names.append("GO:0005576");
			names.append("go_one");

		if model_name == "bacello_fungi":
			names.append("size");
			names.append("Mono_NLS_Signal27");
			names.append("Vacuole_Signal9");
			names.append("AutoCorrelatuion|d=3|1-20|charge");
			names.append("maxAutoCorrelation|1-110|charge");
			names.append("hist|1-20|aminoacid|L");
			names.append("hist|1-180|aminoacid|A");
			names.append("maxPaacount|1-70|hydrophob|1");
			names.append("Paacount|d=4|-100--1|hydrophob|1");
			names.append("Paacount|d=5|1-20|polar|0");
			names.append("maxPaacountNormed|1--1|polar|1");
			names.append("hist_normed|1-20|charged|1");
			names.append("maxPaacountNormed|1-10|aliphatic_helix|2");
			names.append("hist|1-50|aliphatic_helix|1");
			names.append("Paacount|d=4|1--1|aliphatic_helix|2");
			names.append("prosite_clust_cyt");
			names.append("prosite_clust_nuc");
			names.append("GO:0005739");
			names.append("GO:0005576");
			names.append("go_one");

		if model_name == "bacello_fungi_p+p":
			names.append("size");
			names.append("Mono_NLS_Signal27");
			names.append("Vacuole_Signal9");
			names.append("Vacuole_Signal31");
			names.append("AutoCorrelatuion|d=3|1-20|charge");
			names.append("maxAutoCorrelation|1-110|charge");
			names.append("hist|1-20|aminoacid|L");
			names.append("hist|1-40|aminoacid|A");
			names.append("hist|1-180|aminoacid|G");
			names.append("hist|-90--1|aminoacid|G");
			names.append("maxPaacount|1-70|hydrophob|1");
			names.append("Paacount|d=5|1-20|polar|0");
			names.append("maxPaacountNormed|1--1|polar|1");
			names.append("hist_normed|1-20|charged|1");	#not 0 as in file
			names.append("Paacount|d=4|1--1|size|2");
			names.append("maxPaacountNormed|1-10|aliphatic_helix|2");
			names.append("hist|1-50|aliphatic_helix|1");
			names.append("Paacount|d=4|1--1|aliphatic_helix|2");
			names.append("prosite_clust_cyt");
			names.append("prosite_clust_nuc");


		if model_name == "bacello_plants":
			names.append("size");
			names.append("Mono_NLS_Signal26");
			names.append("Vacuole_Signal37");
			names.append("MIP_signal2");
			names.append("MIP_signal13");
			names.append("aliphatic_helix_normed");
			names.append("AutoCorrelatuion|d=4|1-60|charge");
			names.append("hist_normed|1-10|hydrophob|0")
			names.append("hist|1-50|charged|1");
			names.append("hist|1-170|charged|1");
			names.append("maxPaacountNormed|-100--1|charged|0");
			names.append("Paacount|d=4|1--1|charged|2");
			names.append("maxPaacount|1-160|structure|2");
			names.append("Paacount|d=4|1-140|hydroxylated|0");
			names.append("Paacount|d=4|1--1|aliphatic_helix|2");
			names.append("prosite_clust_nuc");
			names.append("GO:0005739");
			names.append("GO:0005576");
			names.append("GO:0009507");
			names.append("go_one");

		if model_name == "bacello_plants_p+p":
			names.append("size");
			names.append("Mono_NLS_Signal26");
			names.append("Vacuole_Signal37");
			names.append("MIP_signal2");
			names.append("MIP_signal13");
			names.append("aliphatic_helix_normed");
			names.append("AutoCorrelatuion|d=4|1-60|charge");
			names.append("hist|-70--1|aminoacid|C");
 			names.append("hist_normed|1-10|hydrophob|0");
			names.append("hist|1-50|charged|1");
			names.append("hist|1-170|charged|1");
			names.append("maxPaacountNormed|-100--1|charged|0");
			names.append("Paacount|d=4|1--1|charged|2");
			names.append("maxPaacount|1-160|structure|2");
			names.append("Paacount|d=6|1--1|size|2");
			names.append("Paacount|d=4|1-140|hydroxylated|0");
			names.append("maxPaacount|-90--1|aliphatic_helix|0");
			names.append("Paacount|d=4|1--1|aliphatic_helix|2");
			names.append("prosite_clust_cyt");
			names.append("prosite_clust_nuc");

		if model_name == "multiloc_animals":
			names.append("size");
			names.append("Mono_NLS_signal26");
			names.append("ER_signal5");
			names.append("Vacuole_signal31");
			names.append("Golgi_signal4");
			names.append("Transmembrane_signal2");
			names.append("AutoCorrelatuion|d=1|1-20|hydrophob");
			names.append("AutoCorrelatuion|d=6|1-70|hydrophob");
			names.append("property_sum|1-30|charge");
			names.append("AutoCorrelatuion|d=2|1--1|charge");
			names.append("letterCount|1-20|K|aminoacid");
			names.append("maxPaacountNormed|1-50|leucine|L");
			names.append("maxPaacount|1-170|cysteine|C");
			names.append("Paacount|d=2|1-50|hydrophob_very|2");
			names.append("Paacount|d=2|1--1|hydrophob_very|2");
			names.append("letterCount|1-20|1|charged");
			names.append("Paacount|d=2|1--1|charge_none|0");
			names.append("Paacount|d=6|1--1|size_small|1");
			names.append("maxPaacount|1--1|aromatic|1");
			names.append("minPaacount|50--20|aliphatic_helix_uncharged|1");
			names.append("prosite_clust_cyt");
			names.append("prosite_clust_nuc");
			names.append("prosite_clust_mem");
			names.append("GO:0005737");
			names.append("GO:0005783");
			names.append("GO:0005739");
			names.append("GO:0005576");
			names.append("GO:0005777");
			names.append("GO:0005764");
			names.append("go_one");

		if model_name == "multiloc_animals_p+p":
			names.append("size");
			names.append("Mono_NLS_signal26");
			names.append("ER_signal5");
			names.append("Peroxi_Signal3");
			names.append("Vacuole_signal31");
			names.append("Golgi_signal4");
			names.append("Transmembrane_signal2");
			names.append("AutoCorrelatuion|d=1|1-20|hydrophob");
			names.append("AutoCorrelatuion|d=6|1-70|hydrophob");
			names.append("property_sum|1-30|charge");
			names.append("AutoCorrelatuion|d=2|1--1|charge");
			names.append("letterCount|1-20|K|aminoacid");
			names.append("maxPaacount|1-30|aminoacid|N");
			names.append("maxPaacountNormed|1-50|leucine|L");
			names.append("letterCount|1-90|W|aminoacid");
			names.append("maxPaacount|1-170|cysteine|C");
			names.append("Paacount|d=2|1-50|hydrophob_very|2");
			names.append("maxPaacount|-40--1|hydrophob_very|2");
			names.append("Paacount|d=2|1--1|hydrophob_very|2");
			names.append("letterCount|1-20|1|charged");
			names.append("Paacount|d=2|1--1|charge_none|0");
			names.append("Paacount|d=6|1--1|size_small|1");
 			names.append("maxPaacount|1--1|aromatic_none|1");
			names.append("maxPaacountNormed|1--1|hydroxylated_hydroxy|1");
			names.append("minPaacount|50--20|aliphatic_helix_uncharged|1");
			names.append("prosite_clust_cyt");
			names.append("prosite_clust_SR");
			names.append("prosite_clust_mTP");
			names.append("prosite_clust_nuc");
			names.append("prosite_clust_mem");

		if model_name == "multiloc_fungi":
			names.append("size");
			names.append("Mono_NLS_signal26");
			names.append("ER_signal5");
			names.append("Vacuole_signal31");
			names.append("Golgi_signal4");
			names.append("Transmembrane_signal2");
			names.append("MIP_signal13");
			names.append("AutoCorrelatuion|d=2|1-20|hydrophob");
			names.append("AutoCorrelatuion|d=1|1-90|hydrophob");
			names.append("property_sum|1-30|charge");
			names.append("AutoCorrelatuion|d=2|1--1|charge");
			names.append("maxPaacountNormed|1-40|leucine|L");
			names.append("maxPaacount|1-40|aminoacid|K");
			names.append("maxPaacount|1-170|cysteine|C");
			names.append("maxPaacount|1--1|aminoacid|S");
			names.append("Paacount|d=2|1-60|hydrophob_very|2");
			names.append("maxPaacountNormed|1--1|hydrophob|2");
			names.append("hist|1-20|charged|1");
			names.append("Paacount|d=2|1--1|charge_none|0");
			names.append("minPaacount|1--1|structure|2");
			names.append("maxPaacountNormed|1-170|aliphatic_helix|2");
			names.append("prosite_clust_cyt");
			names.append("prosite_clust_nuc");
			names.append("prosite_clust_mem");
			names.append("GO:0005737");
			names.append("GO:0005783");
			names.append("GO:0005739");
			names.append("GO:0005576");
			names.append("GO:0005777");
			names.append("go_one");

		if model_name == "multiloc_fungi_p+p":
			names.append("size");
			names.append("Mono_NLS_signal26");
			names.append("ER_signal5");
			names.append("Peroxi_Signal3");
			names.append("Vacuole_signal31");
			names.append("Golgi_signal4");
			names.append("Transmembrane_signal2");
			names.append("MIP_signal13");
			names.append("AutoCorrelatuion|d=2|1-20|hydrophob");
			names.append("AutoCorrelatuion|d=1|1-90|hydrophob");
			names.append("AutoCorrelatuion|d=3|-30--1|hydrophob");
			names.append("property_sum|1-30|charge");
			names.append("AutoCorrelatuion|d=2|1--1|charge");
			names.append("maxPaacountNormed|1-40|leucine|L");
			names.append("maxPaacount|1-40|lysine|K");
			names.append("maxPaacount|1-170|cysteine|C");
			names.append("maxPaacount|1--1|serine|S");
			names.append("Paacount|d=2|1-60|hydrophob_very|2");
			names.append("maxPaacountNormed|1--1|hydrophob|2");
			names.append("hist|1-20|charged|1");
			names.append("Paacount|d=2|1--1|charge_none|0");
			names.append("minPaacount|1--1|structure|2");
			names.append("maxPaacountNormed|1--1|aromatic_aromatic|0");
			names.append("hist|1-110|hydroxylated_hydroxy|0");
			names.append("maxPaacountNormed|1-170|aliphatic_helix|2");
			names.append("prosite_clust_cyt");
			names.append("prosite_clust_SR");
			names.append("prosite_clust_mTP");
			names.append("prosite_clust_nuc");
			names.append("prosite_clust_mem");

		if model_name == "multiloc_plants":
			names.append("size");
			names.append("Mono_NLS_signal26");
			names.append("ER_signal5");
			names.append("Vacuole_signal31");
			names.append("Golgi_signal4");
			names.append("Golgi_signal6");
			names.append("Transmembrane_signal2");
			names.append("MIP_signal2");
			names.append("property_sum|1-10|volume");
			names.append("AutoCorrelatuion|d=2|1-20|hydrophob");
			names.append("AutoCorrelatuion|d=6|1-70|hydrophob");
			names.append("AutoCorrelatuion|d=6|1-60|charge");
			names.append("maxPaacountNormed|1-40|aminoacid|L");
			names.append("maxPaacount|1-180|aminoacid|C");
			names.append("maxPaacount|1-40|hydrophob|1");
			names.append("Paacount|d=2|1-50|hydrophob|2");
			names.append("Paacount|d=2|1--1|hydrophob|2");
			names.append("hist|1-20|charged|1");
			names.append("Paacount|d=2|1--1|charged|0");
			names.append("hist|1-30|aliphatic_helix|2");
			names.append("maxPaacountNormed|1--1|aliphatic_helix|2");
			names.append("prosite_clust_cyt");
			names.append("prosite_clust_nuc");
			names.append("prosite_clust_mem");
			names.append("GO:0005783");
			names.append("GO:0005739");
			names.append("GO:0005576");
			names.append("GO:0005777");
			names.append("GO:0009507");
			names.append("go_one");

		if model_name == "multiloc_plants_p+p":
			names.append("size");
			names.append("Mono_NLS_signal26");
			names.append("ER_signal5");
			names.append("Vacuole_signal31");
			names.append("Golgi_signal4");
			names.append("Golgi_signal6");
			names.append("Transmembrane_signal2");
			names.append("MIP_signal2");
			names.append("property_sum|1-10|volume");
			names.append("AutoCorrelatuion|d=2|1-20|hydrophob");
			names.append("AutoCorrelatuion|d=6|1-70|hydrophob");
			names.append("AutoCorrelatuion|d=6|1-60|charge");
			names.append("hist|1-20|aminoacid|K");
			names.append("maxPaacountNormed|1-40|aminoacid|L");
			names.append("maxPaacount|1-180|aminoacid|C");
			names.append("maxPaacount|1-40|hydrophob|1");
			names.append("Paacount|d=2|1-50|hydrophob|2");
			names.append("Paacount|d=2|1--1|hydrophob|2");
			names.append("hist|1-20|charged|1");
			names.append("Paacount|d=2|1--1|charged|0");
			names.append("hist|1-160|structure|2");
			names.append("maxPaacount|1--1|aromatic|1");
			names.append("hist|1-110|hydroxylated|1");
			names.append("hist|1-30|aliphatic_helix|2");
			names.append("maxPaacountNormed|1--1|aliphatic_helix|2");
			names.append("prosite_clust_cyt");
			names.append("prosite_clust_SR");
			names.append("prosite_clust_mTP");
			names.append("prosite_clust_nuc");
			names.append("prosite_clust_mem");

		if model_name == "multiloc_animals_dbm":
			names.append("size");
			names.append("ER_signale5");
			names.append("Peroxi_Signal1");
			names.append("Vacuole_signal21");
			names.append("Vacuole_signal31");
			names.append("Transmembrane_signal2");
			names.append("MIP_signal12");
			names.append("AutoCorrelatuion|d=2|1-20|hydrophob");
			names.append("maxAutoCorrelation|1-30|charge");
			names.append("hist|1-70|aminoacid|M");
			names.append("hist|1-70|aminoacid|N");
			names.append("maxPaacount|1-120|aminoacid|C");
			names.append("maxPaacount|1-120|aminoacid|K");
			names.append("hist|1-120|aminoacid|W");
			names.append("Paacount|d=2|1-130|hydrophob|2");
			names.append("maxPaacountNormed|-40--1|hydrophob|2");
			names.append("maxPaacount|1-30|charged|1");
			names.append("minPaacount|1--1|structure|1");
			names.append("hist|1-20|size|1");
			names.append("maxPaacountNormed|1--1|aromatic|0");
			names.append("Paacount|d=3|-100--1|hydroxylated|1");
			names.append("Paacount|d=2|1--1|aliphatic_helix|1");
			names.append("prosite_clust_mem");
			names.append("prosite_clust_nuc");
			names.append("GO:0005783");
			names.append("GO:0005739");
			names.append("GO:0005576");
			names.append("GO:0042025");
			names.append("GO:0005778");
			names.append("go_one");

		if model_name == "multiloc_animals_dbm_p+p":
			names.append("size");
			names.append("ER_signale5");
			names.append("Peroxi_Signal1");
			names.append("Vacuole_signal21");
			names.append("Vacuole_signal31");
			names.append("Transmembrane_signal2");
			names.append("MIP_signal12");
			names.append("AutoCorrelatuion|d=2|1-20|hydrophob");
			names.append("maxAutoCorrelation|1-30|charge");
			names.append("hist|1-70|aminoacid|N");
			names.append("Paacount|d=2|1-90|aminoacid|S");
			names.append("maxPaacount|1-120|aminoacid|C");
			names.append("maxPaacount|1-120|aminoacid|K");
			names.append("hist|1-120|aminoacid|M");
			names.append("hist|1-120|aminoacid|W");
			names.append("Paacount|d=2|1-130|hydrophob|2");
			names.append("maxPaacountNormed|-40--1|hydrophob|2");
			names.append("maxPaacount|1-30|charged|1");
			names.append("minPaacount|1--1|structure|1");
			names.append("hist|1-20|size|1");
			names.append("maxPaacountNormed|1--1|aromatic|0");
			names.append("Paacount|d=3|-100--1|hydroxylated|1");
			names.append("Paacount|d=2|1--1|aliphatic_helix|1");
			names.append("phylogenetic_profile_species_9606");
			names.append("phylogenetic_profile_species_3702");
			names.append("PS00639");
			names.append("PS50041");
			names.append("prosite_clust_SR");
			names.append("prosite_clust_mem");
			names.append("prosite_clust_nuc");
			names.append("prosite_clust_mTP");
			names.append("prosite_clust_cyt");

		if model_name == "multiloc_fungi_dbm":
			names.append("size");
			names.append("ER_signale5");
			names.append("Peroxi_Signal1");
			names.append("Vacuole_signal21");
			names.append("Vacuole_signal31");
			names.append("Transmembrane_signal2");
			names.append("MIP_signal12");
			names.append("AutoCorrelatuion|d=2|1-20|hydrophob");
			names.append("maxAutoCorrelation|1-30|charge");
			names.append("hist|1-70|aminoacid|N");
			names.append("maxPaacount|1-120|aminoacid|K");
			names.append("hist|1-120|aminoacid|W");
			names.append("maxPaacountNormed|1-130|aminoacid|F");
			names.append("maxPaacount|1-130|aminoacid|C");
			names.append("Paacount|d=2|1-130|hydrophob|2");
			names.append("maxPaacountNormed|-40--1|hydrophob|2");
			names.append("maxPaacount|1-30|charged|1");
			names.append("minPaacount|1--1|structure|1");
			names.append("hist|1-20|size|1");
			names.append("Paacount|d=3|-100--1|hydroxylated|1");
			names.append("Paacount|d=2|1--1|aliphatic_helix|1");
			names.append("prosite_clust_mem");
			names.append("prosite_clust_nuc");
			names.append("GO:0005783");
			names.append("GO:0005739");
			names.append("GO:0005576");
			names.append("GO:0042025");
			names.append("GO:0005773");
			names.append("GO:0005778");
			names.append("go_one");

		if model_name == "multiloc_fungi_dbm_p+p":
			names.append("size");
			names.append("ER_signale5");
			names.append("Peroxi_Signal1");
			names.append("Vacuole_signal21");
			names.append("Vacuole_signal31");
			names.append("Transmembrane_signal2");
			names.append("MIP_signal12");
			names.append("AutoCorrelatuion|d=2|1-20|hydrophob");
			names.append("maxAutoCorrelation|1-30|charge");
			names.append("hist|1-70|aminoacid|N");
			names.append("maxPaacountNormed|1-120|aminoacid|Y");
			names.append("maxPaacount|1-120|aminoacid|K");
			names.append("hist|1-120|aminoacid|M");
			names.append("hist|1-120|aminoacid|W");
			names.append("maxPaacountNormed|1-130|aminoacid|F");
        		names.append("maxPaacount|1-130|aminoacid|C");
			names.append("Paacount|d=2|1-130|hydrophob|2");
			names.append("maxPaacountNormed|-40--1|hydrophob|2");
			names.append("maxPaacount|1-30|charged|1");
			names.append("minPaacount|1--1|structure|1");
			names.append("hist|1-20|size|1");
			names.append("maxPaacount|-100--1|aromatic|0");
			names.append("Paacount|d=3|-100--1|hydroxylated|1");
			names.append("Paacount|d=2|1--1|aliphatic_helix|1");
			names.append("PS50041");
			names.append("prosite_clust_SR");
			names.append("prosite_clust_mem");
			names.append("prosite_clust_nuc");
			names.append("prosite_clust_mTP");
			names.append("prosite_clust_cyt");

		if model_name == "multiloc_plants_dbm":
			names.append("size");
			names.append("ER_signale5");
			names.append("Vacuole_signal31");
			names.append("Transmembrane_signal2");
			names.append("MIP_signal12");
			names.append("AutoCorrelatuion|d=1|1-20|hydrophob");
			names.append("property_sum|1-10|charge");
			names.append("AutoCorrelatuion|d=6|1-50|charge");
			names.append("hist|1-50|aminoacid|S");
			names.append("maxPaacountNormed|1-80|aminoacid|L");
			names.append("maxPaacount|1-120|aminoacid|C");
			names.append("maxPaacount|1-120|aminoacid|K");
			names.append("hist|1-120|aminoacid|F");
			names.append("hist|1-120|aminoacid|W");
			names.append("Paacount|d=2|1-100|hydrophob|2");
			names.append("maxPaacountNormed|-40--1|hydrophob|2");
			names.append("maxPaacount|1-30|charged|1");
			names.append("Paacount|d=3|1-20|structure|1");
			names.append("minPaacount|1--1|structure|1");
			names.append("hist|1-20|size|2");
			names.append("maxPaacount|1--1|hydroxylated|1");
			names.append("Paacount|d=2|1--1|aliphatic_helix|1");
			names.append("prosite_clust_conf_80_mem");
			names.append("prosite_clust_conf_80_nuc");
			names.append("GO:0005783");
			names.append("GO:0005739");
			names.append("GO:0005576");
			names.append("GO:0042025");
			names.append("GO:0009507");
			names.append("go_one");

		if model_name == "multiloc_plants_dbm_p+p":
			names.append("size");
			names.append("ER_signale5");
			names.append("Peroxi_Signal1");
			names.append("Vacuole_signal31");
			names.append("Transmembrane_signal2");
			names.append("MIP_signal12");
			names.append("Chloro_signal2");
			names.append("AutoCorrelatuion|d=1|1-20|hydrophob");
			names.append("property_sum|1-10|charge");
			names.append("AutoCorrelatuion|d=6|1-50|charge");
			names.append("hist|1-50|aminoacid|S");
			names.append("maxPaacountNormed|1-80|aminoacid|L");
			names.append("maxPaacount|1-120|aminoacid|C");
			names.append("maxPaacount|1-120|aminoacid|K");
			names.append("hist|1-120|aminoacid|F");
			names.append("hist|1-120|aminoacid|W");
			names.append("Paacount|d=2|1-100|hydrophob|2");
			names.append("maxPaacountNormed|-40--1|hydrophob|2");
			names.append("maxPaacount|1-30|charged|1");
			names.append("Paacount|d=3|1-20|structure|1");
			names.append("minPaacount|1--1|structure|1");
			names.append("hist|1-20|size|2");
			names.append("maxPaacount|1--1|hydroxylated|1");
			names.append("Paacount|d=2|1--1|aliphatic_helix|1");
			names.append("prosite_clust_conf_80_SR");
			names.append("prosite_clust_conf_80_mem");
			names.append("prosite_clust_conf_80_nuc");
			names.append("prosite_clust_conf_80_mTP");
			names.append("prosite_clust_conf_80_cTP");
			names.append("prosite_clust_conf_80_cyt");

		if model_name == "suba_plus":
			names.append("size");
			names.append("AutoCorrelatuion|d=1|1-20|hydrophob");
			names.append("AutoCorrelatuion|d=1|1--1|hydrophob");
			names.append("AutoCorrelation|d=1|1-120|charge");
			names.append("hist|1-10|aminoacid|A");
			names.append("hist|1-10|aminoacid|K");
			names.append("hist|1-20|aminoacid|G");
			names.append("hist|1-70|aminoacid|Y");
			names.append("maxPaacountNormed|1-80|aminoacid|W");
			names.append("maxPaacount|1-80|aminoacid|Q");
			names.append("maxPaacountNormed|1-190|aminoacid|S");
			names.append("hist|1-190|aminoacid|C");
			names.append("hist|-80--1|aminoacid|S");
			names.append("Paacount|d=4|-70--1|aminoacid|P");
			names.append("maxPaacount|-60--1|aminoacid|L");
			names.append("hist|-30--1|aminoacid|I");
			names.append("maxPaacountNormed|-20--1|aminoacid|E");
			names.append("hist|1-40|charged|1");
			names.append("Paacount|d=3|1--1|structure|1");
			names.append("maxPaacount|1-40|size|1");
			names.append("Paacount|d=3|1--1|size|2");
			names.append("phylogenetic_profile_species_240292");
			names.append("phylogenetic_profile_species_39947");
			names.append("PS50191");
			names.append("PS50920");
			names.append("prosite_clust_conf_80_nuc");
			names.append("GO:0005789");
			names.append("GO:0005576");
			names.append("GO:0004298");
			names.append("go_one");

		return names;

	def getModelFeatureDescription(self, model_name,features, nr):
		description = [];
		# each desription includes [ understandable feature name, detailed description of feature ,  sequence property, biological property, functional property]
		# plus code 1 = low/medium/high       2=small/average/large         3=present/no present    4=barely/medium/very
		# 5 = typical/unusal   6=no present/present
		if model_name == "bacello_animals_p+p":
			description.append([[2,"number of amino acids"],[2,"protein size"]]);
			description.append([[2,"number non-hydrophobic amino acids in N-terminus and mono NLS signals"],[4,"non-hydrophobic N-terminus and mono NLS signals"],[7,"mono NLS sorting signal"]]);
			description.append([[2,"weighted sum of typical amino acids for mitochondrial and secreted proteins in N-terminus"],[7,"putative mitochondrial or secretory pathway sorting signal"]]);
			description.append([[1,"autocorrelation of every third hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[13,"secretory pathway sorting signal"]]);
			description.append([[1,"autocorrelation of every fifth charged amino acid within the first 30 amino acids in the N-terminus"],[4,"charged N-terminus"],[12,"putative mitochondrial sorting signal"]]);

			description.append([[1,"overall autocorrelation of charged amino acid"],[4,"charged protein"]]);
			description.append([[2,"number of tryptophane within the first 20 amino acids in the N-terminus"],[2,"number of tryptophane in the N-terminus"]]);
			description.append([[1,"pseudo amino acid count of leucine residues in a distant of three within the first 30 amino acids in the N-terminus"],[2,"number of LxxL patterns in the N-terminus"]]);
			description.append([[2,"number of cysteine within the first 30 amino acids in the N-terminus"],[2,"number of cysteine in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of alanine within the the first 60 amino acids in the N-terminus"],[2,"number of alanine in the N-terminus"]]);

			description.append([[1,"overall pseudo amino acid count of cysteine residues in a distance of two"],[2,"number of CxC patterns"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of very hydrophobic residues within the the first 90 amino acids in the N-terminus"],[2,"number of very hydrophobic amino acids in the N-terminus"],[4,"hydrophobic N-terminus"]]);
			description.append([[2,"normalized number of negatively charged residues within the first 20 amino acids in the N-terminus"],[4,"negatively charged N-terminus"]]);
			description.append([[1,"minimal pseudo amino acid count of uncharged residues within the the last 80 amino acids in the C-terminus"],[4,"uncharged C-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of turn favored amino acids within the the last 100 amino acids in the C-terminus"],[4,"loop rich C-terminus"]]);

			description.append([[1,"overall maximal pseudo amino acid count of hydrophobic amino acids that often occur in beta strands [CITVWY]"],[4,"hydrophobic protein"]]);
			descr = (features.getDescriptions())[nr];
			p = ""
			if 15 in descr.keys() and len(descr[15])>1:
				p = "(" + descr[15] + ") ";
			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);
			p = ""
			if 16 in descr.keys() and len(descr[16])>1:
				p = "(" + descr[16] + ") ";
			description.append([[6,"typical secreted pathway prosite pattern "+p],[6,"typical secreted pathway prosite pattern"]]);
			p = ""
			if 17 in descr.keys() and len(descr[17])>1:
				p = "(" + descr[17] + ") ";
			description.append([[6,"typical mitochondrial prosite pattern "+p],[6,"typical mitochondrial prosite pattern"]]);
			p = ""
			if 18 in descr.keys() and len(descr[18])>1:
				p = "(" + descr[18] + ") ";
			description.append([[6,"prosite pattern "+p+", typical for nuclear proteins"],[6,"prosite pattern "+p]]);
		#features.addPrositeMotifs(["clust_cyt","clust_SR","clust_mTP","clust_nuc"]);

		# plus code 1 = low/medium/high       2=small/average/large         3=present/no present    4=barely/medium/very
		# 5 = typical/unusal   6=no present/present  7=no/weak/strong
		# 8=negative/neutral/positive
		# 9=go clusters
		if model_name == "bacello_animals":
			description.append([[2,"number of amino acids"],[2,"protein size"]]);
			description.append([[2,"number non-hydrophobic amino acids in N-terminus and mono NLS signals"],[7,"mono NLS sorting signal"]]);
			description.append([[2,"weighted sum of typical amino acids for mitochondrial and secreted proteins in N-terminus"],[7,"putative mitochondrial or secretory pathway sorting signal"]]);
			description.append([[1,"autocorrelation of every third hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[13,"secretory pathway sorting signal"]]);
			description.append([[1,"autocorrelation of every fifth charged amino acid within the first 30 amino acids in the N-terminus"],[4,"charged N-terminus"],[12,"putative mitochondrial sorting signal"]]);
			description.append([[1,"overall autocorrelation of charged amino acid"],[4,"charged protein"]]);
			description.append([[1,"pseudo amino acid count of Lysine residues in a distant of three within the first 20 amino acids in the N-terminus"],[2,"number of LxxL patterns in the N-terminus"]]);
			description.append([[2,"number of Cysteine within the first 30 amino acids in the N-terminus"],[2,"number of Cysteine in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of Cysteine within the the last 30 amino acids in the C-terminus"],[2,"number of Cysteine in the C-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of very hydrophobic residues within the the first 90 amino acids in the N-terminus"],[2,"number of very hydrophobic amino acids in the N-terminus"],[4,"hydrophobic N-terminus"]]);

			description.append([[2,"normalized number of negatively charged residues within the first 20 amino acids in the N-terminus"],[4,"negatively charged N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of turn favored amino acids within the the last 100 amino acids in the C-terminus"],[4,"loop rich C-terminus"]]);
			description.append([[1,"overall maximal pseudo amino acid count of hydrophobic amino acids that often occur in beta strands [CITVWY]"],[4,"hydrophobic protein"]]);
			descr = (features.getDescriptions())[nr];
			p = ""
			if 13 in descr.keys() and len(descr[13])>1:
				p = "(" + descr[13] + ") ";
			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);
			p = ""
			if 14 in descr.keys() and len(descr[14])>1:
				p = "(" + descr[14] + ") ";
			description.append([[6,"typical plasma membrane prosite pattern "+p],[6,"typical plasma membrane prosite pattern"]]);
			p = ""
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a> typical for mitochondrial proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005634>GO:0005634 (nucleus)</a>  typical for nuclear proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005634>GO:0005634 (nucleus)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (transferase activity)</a> untypical for secreted proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (transferase activity)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a> typical for extracellular proteins"],[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a>"]]);
			p = ""
			if 19 in descr.keys() and len(descr[19])>1:
				if descr[19] != "":
					p = "(GO terms: " + str(descr[19]) + ")";
			description.append([[9,"typical GO term "+p],[9,"typical GO term"]]);


		if model_name == "bacello_fungi":
			description.append([[2,"number of amino acids"],[2,"protein size"]]);
			description.append([[1,"number of NLS signal sequences from NLSDB."],[6,"NLS sorting signal"]]);
			description.append([[1,"autocorrelation of hydrophilic amino acids within the first 10 amino acids in the N-terminus"],[4,"hydrophilic N-terminus"],[6,"alternative NLS sorting signal"]]);
			description.append([[1,"autocorrelation of every third charged amino acid within the first 20 amino acids in the N-terminus"],[4,"positive charged amphiphilic helix in the N-terminus"],[11,"mitochondrial targeting peptide"]]);
			description.append([[1,"maximal autocorrelation of charged amino acids within the first 110 amino acids in the N-terminus"],[4,"positively charged N-terminus"]]);

			description.append([[2,"number of leucine within the first 20 amino acids in the N-terminus"],[2,"number of leucine in the N-terminus"]]);
			description.append([[2,"number of alanine within the first 180 amino acids in the N-terminus"],[2,"number of alanine in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of slighly hydrophobic residues within the the first 70 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"]]);
			description.append([[1,"pseudo amino acid count of slightly hydrophobic residues in a distant of four within the last 100 amino acids in the C-terminus"],[4,"hydrophobic C-terminus"]]);
			description.append([[1,"pseudo amino acid count of unpolar residues in a distant of five within the first 20 amino acids in the N-terminus"],[4,"unpolar N-terminus"],[6,"secretory pathway sorting signal"]]);
			description.append([[1,"maximal normalized overall pseudo amino acid count of polar residues"],[2,"number of polar amino acids"]]);

			description.append([[2,"normalized number of negatively charged residues within the first 20 amino acids in the N-terminus"],[4,"negatively charged N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of [DEHNPQS] within the the first 10 amino acids in the N-terminus"],[2,"number of polar but not positively charged residues in N-terminus"]]);
			description.append([[2,"number of hydrophic and uncharged residues within the first 50 amino acids in the N-terminus"],[4,"hydrophobic and uncharged N-terminus"]]);
			description.append([[1,"Overall pseudo amino acid count of [DEHNPQS] in a distance of four"],[4,"unpolar N-terminus"],[4,"polar but not positively charged protein"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 15 in descr.keys() and len(descr[15])>1:
				p = "(" + descr[15] + ") ";
			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);
			p = ""
			if 16 in descr.keys() and len(descr[16])>1:
				p = "(" + descr[16] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a> typical for mitochondrial proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a> typical for extracellular proteins"],[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a>"]]);
			p = ""
			if 19 in descr.keys() and len(descr[19])>1:
				if descr[19] != "":
					p = "(GO terms: " + str(descr[19]) + ")";
			description.append([[9,"typical GO term "+p],[9,"typical GO term"]]);

		if model_name == "bacello_fungi_p+p":
			description.append([[2,"number of amino acids"],[2,"protein size"]]);
			description.append([[1,"number of NLS signal sequences from NLSDB."],[6,"NLS sorting signal"]]);
			description.append([[1,"autocorrelation of hydrophilic amino acids within the first 10 amino acids in the N-terminus"],[4,"hydrophilic N-terminus"],[6,"alternative NLS sorting signal"]]);
			description.append([[2,"sum of hydrophobicity of 10 amino acids before the N-terminal cleavage site."],[4,"hydrophobic region before the cleavage site"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[1,"autocorrelation of every third charged amino acid within the first 20 amino acids in the N-terminus"],[4,"positive charged amphiphilic helix in the N-terminus"],[11,"mitochondrial targeting peptide"]]);
			description.append([[1,"maximal autocorrelation of charged amino acids within the first 110 amino acids in the N-terminus"],[4,"positively charged N-terminus"]]);
			description.append([[2,"number of leucine within the first 20 amino acids in the N-terminus"],[2,"number of leucine in the N-terminus"]]);

			description.append([[2,"number of alanine within the first 40 amino acids in the N-terminus"],[2,"number of alanine in the N-terminus"]]);
			description.append([[2,"number of glycine within the first 180 amino acids in the N-terminus"],[2,"number of glycine in the N-terminus"]]);
			description.append([[2,"number of glycine within the last 90 amino acids in the C-terminus"],[2,"number of glycine in the C-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of slighly hydrophobic residues within the the first 70 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"]]);
			description.append([[1,"pseudo amino acid count of unpolar residues in a distant of five within the first 20 amino acids in the N-terminus"],[4,"unpolar N-terminus"],[6,"secretory pathway sorting signal"]]);

			description.append([[1,"maximal normalized overall pseudo amino acid count of polar residues"],[2,"number of polar amino acids"]]);
			description.append([[2,"normalized number of negatively charged residues within the first 20 amino acids in the N-terminus"],[4,"negatively charged N-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of large residues in a distance of four"],[2,"number of large amino acids"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of [DEHNPQS] within the the first 10 amino acids in the N-terminus"],[2,"number of polar but not positively charged residues in N-terminus"]]);
			description.append([[2,"number of hydrophic and uncharged residues within the first 50 amino acids in the N-terminus"],[4,"hydrophobic and uncharged N-terminus"]]);

			description.append([[1,"overall pseudo amino acid count of [DEHNPQS] in a distance of four"],[4,"unpolar N-terminus"],[4,"polar but not positively charged protein"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 18 in descr.keys() and len(descr[18])>1:
				p = "(" + descr[18] + ") ";
			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);
			p = ""
			if 19 in descr.keys() and len(descr[19])>1:
				p = "(" + descr[19] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);


		if model_name == "bacello_plants":
			description.append([[2,"number of amino acids"],[2,"protein size"],[]]);
			description.append([[2,"number non-hydrophobic amino acids in N-terminus and mono NLS signals"],[6,"mono NLS sorting signal"]]);
			description.append([[2,"sum of hydrophobicity of the first 20 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[6,"secretory pathway sorting signal"]]);
			description.append([[2,"number of positively charged amino acid"],[6,"alternative mitochondrial targeting peptide"]]);
			description.append([[2,"weighted sum of typical amino acids for mitochondrial and secreted proteins in N-terminus"],[7,"putative mitochondrial or secretory pathway sorting signal"]]);
			description.append([[4,"maximum hydrophobic moment in a window of 18 residues"],[6,"amphiphilic helix"]]);
			description.append([[1,"autocorrelation of every fourth charged amino acid within the first 60 amino acids in the N-terminus"],[4,"positive charged amphiphilic helix in the N-terminus"],[6,"mitochondrial or chloroplast targeting peptide"]]);
			description.append([[2,"normalized number of hydrophilic residues within the first 10 amino acids in the N-terminus"],[4,"hydrophilic N-terminus"]]);
			description.append([[2,"number of negatively charged residues within the first 50 amino acids in the N-terminus"],[4,"negatively charged very N-terminus"]]);
			description.append([[2,"number of negatively charged residues within the first 170 amino acids in the N-terminus"],[4,"negatively charged N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of uncharged residues within the the last 100 amino acids in the C-terminus"],[4,"uncharged C-terminus"]]);

			description.append([[1,"overall pseudo amino acid count of positively charged residues in a distant of four"],[4,"positve charge of protein"]]);
			description.append([[1,"maximal pseudo amino acid count of hydrophobic amino acids that often occur in beta strands [CITVWY] within the first 160 amino acids"],[4,"hydrophobic N-terminus"]]);
			description.append([[1,"pseudo amino acid count of potentially hydroxylated residues in a distant of four within the first 140 amino acids in the N-terminus"],[2,"number of potentially hydroxylated residues in the N-terminus"]]);
			description.append([[1,"Overall pseudo amino acid count of [DEHNPQS] in a distance of four"],[4,"unpolar N-terminus"],[4,"polar but not positively charged protein"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 15 in descr.keys() and len(descr[15])>1:
				p = "(" + descr[15] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a> typical for mitochondrial proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a> typical for extracellular proteins"],[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0009507>GO:0009507 (chloroplast)</a> typical for chloroplast proteins"],[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0009507>GO:0009507 (chloroplast)</a>"]]);
			p = ""
			if 19 in descr.keys() and len(descr[19])>1:
				if descr[19] != "":
					p = "(GO terms: " + str(descr[19]) + ")";
			description.append([[9,"typical GO term "+p],[9,"typical GO term"]]);

		if model_name == "bacello_plants_p+p":
			description.append([[2,"number of amino acids"],[2,"protein size"],[]]);
			description.append([[2,"number non-hydrophobic amino acids in N-terminus and mono NLS signals"],[6,"mono NLS sorting signal"]]);
			description.append([[2,"sum of hydrophobicity of the first 20 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[6,"secretory pathway sorting signal"]]);
			description.append([[2,"number of positively charged amino acid"],[6,"alternative mitochondrial targeting peptide"]]);
			description.append([[2,"weighted sum of typical amino acids for mitochondrial and secreted proteins in N-terminus"],[7,"putative mitochondrial or secretory pathway sorting signal"]]);

			description.append([[4,"maximum hydrophobic moment in a window of 18 residues"],[6,"amphiphilic helix"]]);
			description.append([[1,"autocorrelation of every fourth charged amino acid within the first 60 amino acids in the N-terminus"],[4,"positive charged amphiphilic helix in the N-terminus"],[6,"mitochondrial or chloroplast targeting peptide"]]);
			description.append([[2,"number of cysteine within the last 70 amino acids in the C-terminus"],[2,"number of cysteine in the C-terminus"]]);
			description.append([[2,"normalized number of hydrophilic residues within the first 10 amino acids in the N-terminus"],[4,"hydrophilic N-terminus"]]);
			description.append([[2,"number of negatively charged residues within the first 50 amino acids in the N-terminus"],[4,"negatively charged very N-terminus"]]);

			description.append([[2,"number of negatively charged residues within the first 170 amino acids in the N-terminus"],[4,"negatively charged N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of uncharged residues within the the last 100 amino acids in the C-terminus"],[4,"uncharged C-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of positively charged residues in a distant of four"],[4,"positve charge of protein"]]);
			description.append([[1,"maximal pseudo amino acid count of hydrophobic amino acids that often occur in beta strands [CITVWY] within the first 160 amino acids"],[4,"hydrophobic N-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of large residues in a distance of six"],[2,"number of large amino acids"]]);
			description.append([[1,"pseudo amino acid count of potentially hydroxylated residues in a distant of four within the first 140 amino acids in the N-terminus"],[2,"number of potentially hydroxylated residues in the N-terminus"]]);

			description.append([[1,"maximal pseudo amino acid count of basic residues within the last 90 amino acids in the C-terminus"],[2,"number of basic residues in the C-terminus"]]);
			description.append([[1,"Overall pseudo amino acid count of [DEHNPQS] in a distance of four"],[4,"unpolar N-terminus"],[4,"polar but not positively charged protein"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 18 in descr.keys() and len(descr[18])>1:
				p = "(" + descr[18] + ") ";
			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);
			p = ""
			if 19 in descr.keys() and len(descr[19])>1:
				p = "(" + descr[19] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);

		if model_name == "multiloc_animals":
			description.append([[2,"number of amino acids"],[2,"protein size"],[]]);
			description.append([[1,"number of NLS signal sequences from NLSDB."],[6,"NLS sorting signal"]]);
			description.append([[1,"number of weak ER retention signals (number of K,D,E,L in the very C-terminus) at C-terminus times a factor that depends on the fact that no transmembrane helix and no glycolization signal is present. If a strong ER retention signal (KDEL, [KRHQSA][DENQ]EL, DEL, EL) at C-terminus is present it is 100."],[3,"KDEL or derivate in C-terminus"],[7,"ER retention signal"]]);
			description.append([[2,"sum of hydrophobicity of 10 amino acids before the N-terminal cleavage site."],[4,"hydrophobic region before the cleavage site"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[2,"maximal number of hydrophobic unpolar residues in a window of size 24 within the last 150 residues in the C-terminus"],[4,"hydrophobic and unpolar C-terminus"]]);
			description.append([[2,"maximal length of very hydrophobic region where the iterative sum of hydrophobicity normed to mean zero drops below zero in at most two cases"],[2,"length of longest very hydrophobic region"]]);
			description.append([[1,"autocorrelation of every hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic very N-terminus"],[7,"secretory pathway sorting signal"]]);
			description.append([[1,"autocorrelation of every sixth hydrophobic amino acid within the first 70 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[4,"hydrophobic N-terminus"]]);
			description.append([[1,"sum of charge of first 30 amino acids in N-terminus"],[7,"mitochondrial targeting peptide"]]);
			description.append([[1,"autocorrelation of every second charged amino acid"],[4,"positive charged amphiphilic helix in the N-terminus"],[4,"charged protein"]]);

			description.append([[2,"number of lysine within the first 10 amino acids in the N-terminus"],[2,"number of lysine in the N-terminus"],[6,"weak NLS signal sequence"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of leucine within the first 50 amino acids in the N-terminus"],[1,"number of leucine in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of cysteine within the first 170 amino acids in the N-terminus"],[1,"number of cysteine in the N-terminus"]]);
			description.append([[1,"pseudo amino acid count of very hydrophobic residues in a distance of two within the first 50 amino acids in the N-terminus"],[2,"number of very hydrophobic residues in the N-terminus"],[7,"alternative secretory pathway sorting signal"]]);

			description.append([[1,"overall pseudo amino acid count of very hydrophobic residues in a distance of two"],[2,"number of very hydrophobic residues"],[4,"hydrophobic protein"]]);
			description.append([[2,"number of negatively charged residues within the first 20 amino acids in the N-terminus"],[4,"negatively charged  N-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of uncharged residues in a distance of two"],[2,"number of uncharged residues in the N-terminus"],[4,"uncharged protein"]]);
			description.append([[1,"overall pseudo amino acid count of small residues in a distance of six"],[2,"number of small residues"]]);
			description.append([[1,"maximal overall pseudo amino acid count of non-aromatic residues"],[1,"number of non-aromatic residues"]]);
			description.append([[1,"minimal pseudo amino acid count of hydrophic uncharged residues [ILVMFYWCTAG] from the amino acid at position 50 to the 20th last amino acid"],[2,"number of hydrophobic uncharged residues in the middle part of the protein"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 20 in descr.keys() and len(descr[20])>1:
				p = "(" + descr[20] + ") ";
			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);
			p = ""
			if 21 in descr.keys() and len(descr[21])>1:
				p = "(" + descr[21] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			if 22 in descr.keys() and len(descr[22])>1:
				p = "(" + descr[22] + ") ";
			description.append([[6,"typical plasma membrane prosite pattern "+p],[6,"typical plasma membrane prosite pattern"]]);
			p = ""
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005737>GO:0005737 (cytoplasm)</a> typical for cytoplasmatic proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005737>GO:0005737 (cytoplasm)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005783>GO:0005783 (endoplasmic reticulum)</a> typical for ER proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005783>GO:0005783 (endoplasmic reticulum)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a> typical for mitochondrial proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a> typical for extracellular proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005777>GO:0005777 (peroxisome)</a> typical for peroxisomal proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005777>GO:0005777 (peroxisome)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005764>GO:0005764 (lysosome)</a> typical for lysosomal proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005764>GO:0005764 (lysosome)</a>"]]);

			p = ""
			if 29 in descr.keys() and len(descr[29])>1:
				if descr[29] != "":
					p = "(GO terms: " + str(descr[29]) + ")";
			description.append([[9,"typical GO term "+p],[9,"typical GO term"]]);

		if model_name == "multiloc_animals_p+p":
			description.append([[2,"number of amino acids"],[2,"protein size"],[]]);
			description.append([[1,"number of NLS signal sequences from NLSDB."],[6,"NLS sorting signal"]])
			description.append([[1,"number of weak ER retention signals (number of K,D,E,L in the very C-terminus) at C-terminus times a factor that depends on the fact that no transmembrane helix and no glycolization signal is present. If a strong ER retention signal (KDEL, [KRHQSA][DENQ]EL, DEL, EL) at C-terminus is present it is 100."],[3,"KDEL or derivate in C-terminus"],[7,"ER retention signal"]]);
			description.append([[6,"C-terminal PTS1 signal sequence or N-terminal PTS2 signal sequence"],[6,"peroxisomal targeting sequence (PTS)"],[]]);
			description.append([[2,"sum of hydrophobicity of 10 amino acids before the N-terminal cleavage site."],[4,"hydrophobic region before the cleavage site"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[2,"maximal number of hydrophobic unpolar residues in a window of size 24 within the last 150 residues in the C-terminus"],[4,"hydrophobic and unpolar C-terminus"]]);

			description.append([[2,"maximal length of very hydrophobic region where the iterative sum of hydrophobicity normed to mean zero drops below zero in at most two cases"],[2,"length of longest very hydrophobic region"]]);
			description.append([[1,"autocorrelation of every hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic very N-terminus"],[7,"secretory pathway sorting signal"]]);
			description.append([[1,"autocorrelation of every sixth hydrophobic amino acid within the first 70 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[4,"hydrophobic N-terminus"]]);
			description.append([[1,"sum of charge of first 30 amino acids in N-terminus"],[7,"mitochondrial targeting peptide"]]);
			description.append([[1,"autocorrelation of every second charged amino acid"],[4,"positive charged amphiphilic helix in the N-terminus"],[4,"charged protein"]]);

			description.append([[2,"number of lysine within the first 20 amino acids in the N-terminus"],[2,"number of lysine in the N-terminus"],[6,"weak NLS signal sequence"]]);
			description.append([[1,"maximal pseudo amino acid count of asparagine within the first 30 amino acids in the N-terminus"],[1,"number of asparagine in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of leucine within the first 50 amino acids in the N-terminus"],[1,"number of leucine in the N-terminus"]]);
			description.append([[2,"number of tryptophane within the first 90 amino acids in the N-terminus"],[2,"number of tryptophane in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of cysteine within the first 170 amino acids in the N-terminus"],[1,"number of cysteine in the N-terminus"]]);
			description.append([[1,"pseudo amino acid count of very hydrophobic residues in a distance of two within the first 50 amino acids in the N-terminus"],[2,"number of very hydrophobic residues in the N-terminus"],[7,"alternative secretory pathway sorting signal"]]);

			description.append([[1,"maximal pseudo amino acid count of very hydrophobic residues within the last 40 amino acids in the C-terminus"],[1,"number of very hydrophobic residues in the C-terminus"],[4,"hydrophobic C-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of very hydrophobic residues in a distance of two"],[2,"number of very hydrophobic residues"],[4,"hydrophobic protein"]]);
			description.append([[2,"number of negatively charged residues within the first 20 amino acids in the N-terminus"],[4,"negatively charged  N-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of uncharged residues in a distance of two"],[2,"number of uncharged residues in the N-terminus"],[4,"uncharged protein"]]);
			description.append([[1,"overall pseudo amino acid count of small residues in a distance of six"],[2,"number of small residues"]]);

			description.append([[1,"maximal overall pseudo amino acid count of non-aromatic residues"],[1,"number of non-aromatic residues"]]);
			description.append([[1,"maximal normalized overall pseudo amino acid count of potentially hydroxylated residues"],[1,"number of potentially hydroxylated residues"]]);
			description.append([[1,"minimal pseudo amino acid count of hydrophic uncharged residues [ILVMFYWCTAG] from the amino acid at position 50 to the 20th last amino acid"],[2,"number of hydrophobic uncharged residues in the middle part of the protein"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 25 in descr.keys() and len(descr[25])>1:
				p = "(" + descr[25] + ") ";
			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);
			p = ""
			if 26 in descr.keys() and len(descr[26])>1:
				p = "(" + descr[26] + ") ";
			description.append([[6,"typical extracellular prosite pattern "+p],[6,"typical extracellular prosite pattern"]]);
			p = ""
			if 27 in descr.keys() and len(descr[27])>1:
				p = "(" + descr[27] + ") ";
			description.append([[6,"typical mitochondrial prosite pattern "+p],[6,"typical mitochondrial prosite pattern"]]);
			p = ""
			if 28 in descr.keys() and len(descr[28])>1:
				p = "(" + descr[28] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			if 29 in descr.keys() and len(descr[29])>1:
				p = "(" + descr[29] + ") ";
			description.append([[6,"typical plasma membrane prosite pattern "+p],[6,"typical plasma membrane prosite pattern"]]);

		if model_name == "multiloc_fungi":
			description.append([[2,"number of amino acids"],[2,"protein size"]]);
			description.append([[1,"number of NLS signal sequences from NLSDB."],[6,"NLS sorting signal"]]);
			description.append([[1,"number of weak ER retention signals (number of K,D,E,L in the very C-terminus) at C-terminus times a factor that depends on the fact that no transmembrane helix and no glycolization signal is present. If a strong ER retention signal (KDEL, [KRHQSA][DENQ]EL, DEL, EL) at C-terminus is present it is 100."],[3,"KDEL or derivate in C-terminus"],[7,"ER retention signal"]]);
			description.append([[2,"sum of hydrophobicity of 10 amino acids before the N-terminal cleavage site."],[4,"hydrophobic region before the cleavage site"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[2,"maximal number of hydrophobic unpolar residues in a window of size 24 within the last 150 residues in the C-terminus"],[4,"hydrophobic and unpolar C-terminus"]]);

			description.append([[2,"maximal length of very hydrophobic region where the iterative sum of hydrophobicity normed to mean zero drops below zero in at most two cases"],[2,"length of longest very hydrophobic region"]]);
			description.append([[2,"weighted sum of typical amino acids for mitochondrial and secreted proteins in N-terminus"],[7,"putative mitochondrial or secretory pathway sorting signal"]]);
			description.append([[1,"autocorrelation of every second hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic very N-terminus"],[7,"secretory pathway sorting signal"]]);
			description.append([[1,"autocorrelation of every hydrophobic amino acid within the first 90 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[1,"sum of charge of first 30 amino acids in N-terminus"],[7,"mitochondrial targeting peptide"]]);

			description.append([[1,"autocorrelation of every second charged amino acid"],[4,"positive charged amphiphilic helix in the N-terminus"],[4,"charged protein"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of leucine within the first 40 amino acids in the N-terminus"],[1,"number of leucine in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of lysine within the first 40 amino acids in the N-terminus"],[1,"number of lysine in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of cysteine within the first 170 amino acids in the N-terminus"],[1,"number of cysteine in the N-terminus"]]);
			description.append([[1,"overall maximal pseudo amino acid count of serine"],[1,"overall number of serine"]]);
			description.append([[1,"pseudo amino acid count of very hydrophobic residues in a distance of two within the first 60 amino acids in the N-terminus"],[2,"number of very hydrophobic residues in the N-terminus"],[7,"alternative secretory pathway sorting signal"]]);

			description.append([[1,"maximal overall pseudo amino acid count of very hydrophobic residues"],[2,"number of very hydrophobic residues"],[4,"hydrophobic protein"]]);
			description.append([[2,"number of negatively charged residues within the first 20 amino acids in the N-terminus"],[4,"negatively charged  N-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of uncharged residues in a distance of two"],[2,"number of uncharged residues in the N-terminus"],[4,"uncharged protein"]]);
			description.append([[1,"minimal overall pseudo amino acid count of hydrophobic residues often occur in beta strands [CITVWY]"],[2,"number of hydrophobic residues"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of hydrophilic, non-positively charged residues [DEHNPQS] within the first 170 amino acids in the N-terminus"],[1,"number of hydrophilic, non-positively charged residues [DEHNPQS] in the N-terminus"],[4,"hydrophilic and non-positively charged N-terminus"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 21 in descr.keys() and len(descr[21])>1:
				p = "(" + descr[21] + ") ";
			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);
			p = ""
			if 22 in descr.keys() and len(descr[22])>1:
				p = "(" + descr[22] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			if 23 in descr.keys() and len(descr[23])>1:
				p = "(" + descr[23] + ") ";
			description.append([[6,"typical plasma membrane prosite pattern "+p],[6,"typical plasma membrane prosite pattern"]]);
			p = ""
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005737>GO:0005737 (cytoplasm)</a> typical for cytoplasmatic proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005737>GO:0005737 (cytoplasm)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005783>GO:0005783 (endoplasmic reticulum)</a> typical for ER proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005783>GO:0005783 (endoplasmic reticulum)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a> typical for mitochondrial proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a> typical for extracellular proteins"],[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005777>GO:0005777 (peroxisome)</a> typical for peroxisomal proteins"],[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005777>GO:0005777 (peroxisome)</a>"]]);
			p = ""
			if 29 in descr.keys() and len(descr[29])>1:
				if descr[29] != "":
					p = "(GO terms: " + str(descr[29]) + ")";
			description.append([[9,"typical GO term "+p],[9,"typical GO term"]]);

		if model_name == "multiloc_fungi_p+p":
			description.append([[2,"number of amino acids"],[2,"protein size"]]);
			description.append([[1,"number of NLS signal sequences from NLSDB."],[6,"NLS sorting signal"]]);
			description.append([[1,"number of weak ER retention signals (number of K,D,E,L in the very C-terminus) at C-terminus times a factor that depends on the fact that no transmembrane helix and no glycolization signal is present. If a strong ER retention signal (KDEL, [KRHQSA][DENQ]EL, DEL, EL) at C-terminus is present it is 100."],[3,"KDEL or derivate in C-terminus"],[7,"ER retention signal"]]);
			description.append([[6,"C-terminal PTS1 signal sequence or N-terminal PTS2 signal sequence"],[6,"peroxisomal targeting sequence (PTS)"],[]]);
			description.append([[2,"sum of hydrophobicity of 10 amino acids before the N-terminal cleavage site."],[4,"hydrophobic region before the cleavage site"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[2,"maximal number of hydrophobic unpolar residues in a window of size 24 within the last 150 residues in the C-terminus"],[4,"hydrophobic and unpolar C-terminus"]]);

			description.append([[2,"maximal length of very hydrophobic region where the iterative sum of hydrophobicity normed to mean zero drops below zero in at most two cases"],[2,"length of longest very hydrophobic region"]]);
			description.append([[2,"weighted sum of typical amino acids for mitochondrial and secreted proteins in N-terminus"],[7,"putative mitochondrial or secretory pathway sorting signal"]]);
			description.append([[1,"autocorrelation of every second hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic very N-terminus"],[7,"secretory pathway sorting signal"]]);
			description.append([[1,"autocorrelation of every hydrophobic amino acid within the first 90 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[1,"autocorrelation of every third hydrophobic amino acid within the last 30 amino acids in the C-terminus"],[4,"hydrophobic C-terminus"]]);
			description.append([[1,"sum of charge of first 30 amino acids in N-terminus"],[7,"mitochondrial targeting peptide"]]);

			description.append([[1,"autocorrelation of every second charged amino acid"],[4,"positive charged amphiphilic helix in the N-terminus"],[4,"charged protein"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of leucine within the first 40 amino acids in the N-terminus"],[1,"number of leucine in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of lysine within the first 40 amino acids in the N-terminus"],[1,"number of lysine in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of cysteine within the first 170 amino acids in the N-terminus"],[1,"number of cysteine in the N-terminus"]]);
			description.append([[1,"maximal overall pseudo amino acid count of serine"],[1,"number of serines"]]);

			description.append([[1,"pseudo amino acid count of very hydrophobic residues in a distance of two within the first 60 amino acids in the N-terminus"],[2,"number of very hydrophobic residues in the N-terminus"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[1,"maximal overall pseudo amino acid count of very hydrophobic residues"],[2,"number of very hydrophobic residues"],[4,"hydrophobic protein"]]);
			description.append([[2,"number of negatively charged residues within the first 20 amino acids in the N-terminus"],[4,"negatively charged  N-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of uncharged residues in a distance of two"],[2,"number of uncharged residues in the N-terminus"],[4,"uncharged protein"]]);
			description.append([[1,"minimal overall pseudo amino acid count of very hydrophobic residues [CITVWY]"],[2,"number of beta-sheet preferred residues [CITVWY]"]]);

			description.append([[1,"maximal normalized overall pseudo amino acid count of aromatic residues"],[1,"number of aromatic residues"]]);
			description.append([[2,"number of potentially hydroxylated residues within the first 110 amino acids in the N-terminus"],[1,"number of potentially hydroxylated residues in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of hydrophilic, non-positively charged residues [DEHNPQS] within the first 170 amino acids in the N-terminus"],[1,"number of hydrophilic, non-positively charged residues [DEHNPQS] in the N-terminus"],[4,"hydrophilic and non-positively charged N-terminus"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 25 in descr.keys() and len(descr[25])>1:
				p = "(" + descr[25] + ") ";
			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);
			p = ""
			if 26 in descr.keys() and len(descr[26])>1:
				p = "(" + descr[26] + ") ";
			description.append([[6,"typical extracellular prosite pattern "+p],[6,"typical extracellular prosite pattern"]]);
			p = ""
			if 27 in descr.keys() and len(descr[27])>1:
				p = "(" + descr[27] + ") ";
			description.append([[6,"typical mitochondrial prosite pattern "+p],[6,"typical mitochondrial prosite pattern"]]);
			p = ""
			if 28 in descr.keys() and len(descr[28])>1:
				p = "(" + descr[28] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			if 29 in descr.keys() and len(descr[29])>1:
				p = "(" + descr[29] + ") ";
			description.append([[6,"typical plasma membrane prosite pattern "+p],[6,"typical plasma membrane prosite pattern"]]);


		if model_name == "multiloc_plants":
			description.append([[2,"number of amino acids"],[2,"protein size"],]);
			description.append([[1,"number of NLS signal sequences from NLSDB."],[6,"NLS sorting signal"]]);
			description.append([[1,"number of weak ER retention signals (number of K,D,E,L in the very C-terminus) at C-terminus times a factor that depends on the fact that no transmembrane helix and no glycolization signal is present. If a strong ER retention signal (KDEL, [KRHQSA][DENQ]EL, DEL, EL) at C-terminus is present it is 100."],[3,"KDEL or derivate in C-terminus"],[7,"ER retention signal"]]);
			description.append([[2,"sum of hydrophobicity of 10 amino acids before the N-terminal cleavage site."],[4,"hydrophobic region before the cleavage site"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[2,"maximal number of hydrophobic unpolar residues in a window of size 24 within the last 150 residues in the C-terminus"],[4,"hydrophobic and unpolar C-terminus"]]);

			description.append([[2,"maximal normalized number of unpolar residues in a window of size 6,8,10 or 12 within the first 50 residues in the N-terminus"],[4,"unpolar cluster in the N-terminus"]]);
			description.append([[2,"maximal length of very hydrophobic region where the iterative sum of hydrophobicity normed to mean zero drops below zero in at most two cases"],[2,"length of longest very hydrophobic region"]]);
			description.append([[2,"number of positively charged amino acid"],[6,"alternative mitochondrial targeting peptide"]]);
			description.append([[1,"sum of volume of first 10 amino acids in N-terminus"],[2,"number of large amino acids in the N-terminus"]]);
			description.append([[1,"autocorrelation of every second hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic very N-terminus"],[7,"secretory pathway sorting signal"]]);
			description.append([[1,"autocorrelation of every sixth hydrophobic amino acid within the first 70 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[4,"hydrophobic N-terminus"]]);
			description.append([[1,"autocorrelation of every sixth charged amino acid within the first 60 amino acids in the N-terminus"],[8,"charged N-terminus"]]);

			description.append([[1,"maximal normalized pseudo amino acid count of leucine within the first 40 amino acids in the N-terminus"],[1,"number of leucine in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of cysteine within the first 180 amino acids in the N-terminus"],[1,"number of cysteine in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of weakly hydrophic [AGHPSTY] within the first 40 amino acids in the N-terminus"],[1,"number of weakly hydrophobic residues [AGHPSTY] in the N-terminus"],[7,"alternative chloroplast targeting signal"]]);

			description.append([[1,"pseudo amino acid count of very hydrophobic residues in a distance of two within the first 50 amino acids in the N-terminus"],[2,"number of very hydrophobic residues in the N-terminus"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[1,"overall pseudo amino acid count of very hydrophobic residues in a distance of two"],[2,"number of very hydrophobic residues"],[4,"hydrophobic protein"]]);
			description.append([[2,"number of negatively charged residues within the first 20 amino acids in the N-terminus"],[4,"negatively charged  N-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of uncharged residues in a distance of two"],[2,"number of uncharged residues in the N-terminus"],[4,"uncharged protein"]]);
			description.append([[2,"number of hydrophilic, non-positively charged residues [DEHNPQS] within the first 30 amino acids in the N-terminus"],[4,"hydrophilic and non-positively charged N-terminus"]]);
			description.append([[1,"maximal overall normalized pseudo amino acid count of hydrophilic, non-positively charged residues [DEHNPQS]"],[1,"number of hydrophilic, non-positively charged residues [DEHNPQS]"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 21 in descr.keys() and len(descr[21])>1:
				p = "(" + descr[21] + ") ";
			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);
			p = ""
			if 22 in descr.keys() and len(descr[22])>1:
				p = "(" + descr[22] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			if 23 in descr.keys() and len(descr[23])>1:
				p = "(" + descr[23] + ") ";
			description.append([[6,"typical plasma membrane prosite pattern "+p],[6,"typical plasma membrane prosite pattern"]]);
			p = ""
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005783>GO:0005783 (endoplasmic reticulum)</a> typical for endoplasmic reticulum proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005783>GO:0005783 (endoplasmic reticulum)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a> typical for mitochondrial proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a> typical for extracellular proteins"],[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005777>GO:0005777 (peroxisome)</a> typical for peroxisomal proteins"],[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005777>GO:0005777 (peroxisome)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0009507>GO:0009507 (chloroplast)</a> typical for chloroplast proteins"],[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0009507>GO:0009507 (chloroplast)</a>"]]);
			p = ""
			if 29 in descr.keys() and len(descr[29])>1:
				if descr[29] != "":
					p = "(GO terms: " + str(descr[29]) + ")";
			description.append([[9,"typical GO term "+p],[9,"typical GO term"]]);

		if model_name == "multiloc_plants_p+p":
			description.append([[2,"number of amino acids"],[2,"protein size"],]);
			description.append([[1,"number of NLS signal sequences from NLSDB."],[6,"NLS sorting signal"]]);
			description.append([[1,"number of weak ER retention signals (number of K,D,E,L in the very C-terminus) at C-terminus times a factor that depends on the fact that no transmembrane helix and no glycolization signal is present. If a strong ER retention signal (KDEL, [KRHQSA][DENQ]EL, DEL, EL) at C-terminus is present it is 100."],[3,"KDEL or derivate in C-terminus"],[7,"ER retention signal"]]);
			description.append([[2,"sum of hydrophobicity of 10 amino acids before the N-terminal cleavage site."],[4,"hydrophobic region before the cleavage site"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[2,"maximal number of hydrophobic unpolar residues in a window of size 24 within the last 150 residues in the C-terminus"],[4,"hydrophobic and unpolar C-terminus"]]);

			description.append([[2,"maximal normalized number of unpolar residues in a window of size 6,8,10 or 12 within the first 50 residues in the N-terminus"],[4,"unpolar cluster in the N-terminus"]]);
			description.append([[2,"maximal length of very hydrophobic region where the iterative sum of hydrophobicity normed to mean zero drops below zero in at most two cases"],[2,"length of longest very hydrophobic region"]]);
			description.append([[2,"number of positively charged amino acid"],[6,"alternative mitochondrial targeting peptide"]]);
			description.append([[1,"sum of volume of first 10 amino acids in N-terminus"],[2,"number of large amino acids in the N-terminus"]]);
			description.append([[1,"autocorrelation of every second hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic very N-terminus"],[7,"secretory pathway sorting signal"]]);

			description.append([[1,"autocorrelation of every sixth hydrophobic amino acid within the first 70 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[4,"hydrophobic N-terminus"]]);
			description.append([[1,"autocorrelation of every sixth charged amino acid within the first 60 amino acids in the N-terminus"],[8,"charged N-terminus"]]);
			description.append([[2,"number of lysine within the first 20 amino acids in the N-terminus"],[2,"number of lysine in the N-terminus"],[6,"weak NLS signal sequence"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of leucine within the first 40 amino acids in the N-terminus"],[1,"number of leucine in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of cysteine within the first 180 amino acids in the N-terminus"],[1,"number of cysteine in the N-terminus"]]);

			description.append([[1,"maximal pseudo amino acid count of slightly hydrophobic residues within the first 40 amino acids in the N-terminus"],[1,"number of slightly hydrophobic residues in the N-terminus"]]);
			description.append([[1,"pseudo amino acid count of very hydrophobic residues in a distance of two within the first 50 amino acids in the N-terminus"],[2,"number of very hydrophobic residues in the N-terminus"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[1,"overall pseudo amino acid count of very hydrophobic residues in a distance of two"],[2,"number of very hydrophobic residues"],[4,"hydrophobic protein"]]);
			description.append([[2,"number of negatively charged residues within the first 20 amino acids in the N-terminus"],[4,"negatively charged  N-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of uncharged residues in a distance of two"],[2,"number of uncharged residues in the N-terminus"],[4,"uncharged protein"]]);

			description.append([[2,"number of very hydrophobic residues [CITVWY] within the first 160 amino acids in the N-terminus"],[2,"number of beta-sheet preferred residues [CITVWY] in the N-terminus"]]);
			description.append([[1,"maximal overall pseudo amino acid count of non-aromatic residues"],[1,"number of non-aromatic residues"]]);
			description.append([[2,"number of potentially hydroxylated residues within the first 110 amino acids in the N-terminus"],[2,"number of potentially hydroxylated residues in the N-terminus"]]);
			description.append([[2,"number of hydrophilic, non-positively charged residues [DEHNPQS] within the first 30 amino acids in the N-terminus"],[4,"hydrophilic and non-positively charged N-terminus"]]);
			description.append([[1,"maximal overall normalized pseudo amino acid count of hydrophilic, non-positively charged residues [DEHNPQS]"],[1,"number of hydrophilic, non-positively charged residues [DEHNPQS]"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 25 in descr.keys() and len(descr[25])>1:
				p = "(" + descr[25] + ") ";
			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);
			p = ""
			if 26 in descr.keys() and len(descr[26])>1:
				p = "(" + descr[26] + ") ";
			description.append([[6,"typical extracellular prosite pattern "+p],[6,"typical extracellular prosite pattern"]]);
			p = ""
			if 27 in descr.keys() and len(descr[27])>1:
				p = "(" + descr[27] + ") ";
			description.append([[6,"typical mitochondrial prosite pattern "+p],[6,"typical mitochondrial prosite pattern"]]);
			p = ""
			if 28 in descr.keys() and len(descr[28])>1:
				p = "(" + descr[28] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			if 29 in descr.keys() and len(descr[29])>1:
				p = "(" + descr[29] + ") ";
			description.append([[6,"typical plasma membrane prosite pattern "+p],[6,"typical plasma membrane prosite pattern"]]);


		if model_name == "multiloc_animals_dbm":
		# plus code 1 = low/medium/high       2=small/average/large         3=present/no present    4=barely/medium/very
		# 5 = typical/unusal   6=no present/present  7=no/weak/strong
		# 8=negative/neutral/positive
		# 9=go clusters
			description.append([[2,"number of amino acids"],[2,"protein size"]]);
			description.append([[1,"number of weak ER retention signals (number of K,D,E,L in the very C-terminus) at C-terminus times a factor that depends on the fact that no transmembrane helix and no glycolization signal is present. If a strong ER retention signal (KDEL, [KRHQSA][DENQ]EL, DEL, EL) at C-terminus is present it is 100."],[3,"KDEL or derivate in C-terminus"],[7,"ER retention signal"]]);
			description.append([[1,"weighted PTS sum = SKL*0.83+SKF*0.5+[SAGCN][RKH][LIVMAF]*0.25 (presented by Nakai et al.)"],[7,"peroxisomal targeting signal (PTS)"],[7,"peroxisomal targeting signal (PTS)"]]);
			description.append([[1,"product of scores for leucine clusters (1.5 for LL, 2.5 for LLL, 5 for LLLL, 12 for LLLLL, 20 for LLLLLL) in the first 50 amino acids in the N-terminus"],[2,"number of leucine and leucine clusters in the N-terminus"]]);
			description.append([[2,"sum of hydrophobicity of 10 amino acids before the N-terminal cleavage site."],[4,"hydrophobic region before the cleavage site"],[7,"alternative secretory pathway sorting signal"]]);

			description.append([[2,"maximal length of very hydrophobic region where the iterative sum of hydrophobicity normed to mean zero drops below zero in at most two cases"],[2,"length of longest very hydrophobic region"]]);
			description.append([[1,"weighted sum of typical amino acids for mitochondrial proteins in N-terminus"],[5,"amino acids for mitochondrial protein in N-terminus"],[6,"putative mitochondrial sorting signal"]]);
			description.append([[1,"autocorrelation of every second hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[6,"secretory pathway sorting signal"]]);
			description.append([[1,"maximal autocorrelation of every sixth charged amino acid within the first 30 amino acids in the N-terminus"],[4,"positive charged amphiphilic helix in the N-terminus"],[7,"mitochondrial targeting peptide"]]);
			description.append([[2,"number of methionine within the first 70 amino acids in the N-terminus"],[2,"number of methionine in the N-terminus"]]);
			description.append([[2,"number of asparagine within the first 70 amino acids in the N-terminus"],[2,"number of asparagine in the N-terminus"]]);
			description.append([[2,"maximal pseudo amino acid count of cysteine within the first 120 amino acids in the N-terminus"],[2,"number of cysteine in the N-terminus"]]);

			description.append([[1,"maximal normalized pseudo amino acid count of lysine within the the first 120 amino acids in the N-terminus"],[2,"number of lysine in the N-terminus"]]);
			description.append([[2,"number of tryptophane within the first 120 amino acids in the N-terminus"],[2,"number of tryptophane in the N-terminus"]]);
			description.append([[1,"pseudo amino acid count of very hydrophobic residues in a distance of two within the first 130 amino acids in the N-terminus"],[2,"number of very hydrophobic residues in the N-terminus"],[4,"hydrophobic N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of very hydrophobic residues in the last 40 amino acids in the C-terminus"],[2,"number of very hydrophobic residues in the C-terminus"],[4,"hydrophobic C-terminus"]]);

			description.append([[1,"maximal pseudo amino acid count of negatively charged residues within the the first 30 amino acids in the N-terminus"],[4,"negatively charged N-terminus"]]);
			description.append([[1,"minimal overall pseudo amino acid count of alpha-helix preferred residues [AFGHKLMNR]"],[2,"number of alpha-helix preferred residues [AFGHKLMNR]"]]);
			description.append([[2,"number of small amino acids within the first 20 amino acids in the N-terminus"],[2,"number of small amino acids in the N-terminus"]]);
			description.append([[1,"maximal normalized overall pseudo amino acid count of aromatic residues"],[1,"number of aromatic residues"]]);
			description.append([[1,"pseudo amino acid count of potentially hydroxylated residues in a distance of three within the last 100 amino acids in the C-terminus"],[2,"number of potentially hydroxylated residues in the C-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of uncharged hydrophobic residues (I,L,V,M,F,Y,W,C,T,A,G) in a distant of two"],[4,"uncharged and hydrophobic"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 22 in descr.keys() and len(descr[22])>1:
				if descr[22] != "":
					p = "(" + descr[22] + ") ";
			description.append([[6,"typical plasma membrane prosite pattern "+p],[6,"typical plasma membrane prosite pattern"]]);
			p = ""
			if 23 in descr.keys() and len(descr[23])>1:
				if descr[23] != "":
					p = "(" + descr[23] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			description.append([[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005783>GO term GO:0005783 (endoplasmic reticulum)</a> typical for ER proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005783>GO:0005783 (endoplasmic reticulum)</a>"]]);

			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a> typical for mitochondrial proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a> typical for extracellular proteins"],[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0042025>GO:0042025 (host cell nucleus)</a> typical for viral proteins present in extracellular space and nucleus"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0042025>GO:0042025 (host cell nucleus)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005778>GO:0005778 (peroxisomal membrane)</a> typical for peroxisomal proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005778>GO:0005778 (peroxisomal membrane)</a>"]]);
			p = ""
			if 29 in descr.keys() and len(descr[29])>1:
				if descr[29] != "":
					p = "(GO terms: " + str(descr[29]) + ")";
			description.append([[9,"typical GO term "+p],[9,"typical GO term"]]);

		if model_name == "multiloc_animals_dbm_p+p":
		# plus code 1 = low/medium/high       2=small/average/large         3=present/no present    4=barely/medium/very
		# 5 = typical/unusal   6=no present/present  7=no/weak/strong
		# 8=negative/neutral/positive
		# 9=go clusters
			description.append([[2,"number of amino acids"],[2,"protein size"]]);
			description.append([[1,"number of weak ER retention signals (number of K,D,E,L in the very C-terminus) at C-terminus times a factor that depends on the fact that no transmembrane helix and no glycolization signal is present. If a strong ER retention signal (KDEL, [KRHQSA][DENQ]EL, DEL, EL) at C-terminus is present it is 100."],[3,"KDEL or derivate in C-terminus"],[7,"ER retention signal"]]);
			description.append([[1,"weighted PTS sum = SKL*0.83+SKF*0.5+[SAGCN][RKH][LIVMAF]*0.25 (presented by Nakai et al.)"],[7,"peroxisomal targeting signal (PTS)"],[7,"peroxisomal targeting signal (PTS)"]]);
			description.append([[1,"product of scores for leucine clusters (1.5 for LL, 2.5 for LLL, 5 for LLLL, 12 for LLLLL, 20 for LLLLLL) in the first 50 amino acids in the N-terminus"],[2,"number of leucine and leucine clusters in the N-terminus"]]);
			description.append([[2,"sum of hydrophobicity of 10 amino acids before the N-terminal cleavage site."],[4,"hydrophobic region before the cleavage site"],[7,"alternative secretory pathway sorting signal"]]);

			description.append([[2,"maximal length of very hydrophobic region where the iterative sum of hydrophobicity normed to mean zero drops below zero in at most two cases"],[2,"length of longest very hydrophobic region"]]);
			description.append([[1,"weighted sum of typical amino acids for mitochondrial proteins in N-terminus"],[5,"amino acids for mitochondrial protein in N-terminus"],[6,"putative mitochondrial sorting signal"]]);
			description.append([[1,"autocorrelation of every second hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[6,"secretory pathway sorting signal"]]);
			description.append([[1,"maximal autocorrelation of every sixth charged amino acid within the first 30 amino acids in the N-terminus"],[4,"positive charged amphiphilic helix in the N-terminus"],[7,"mitochondrial targeting peptide"]]);
			description.append([[2,"number of asparagine within the first 70 amino acids in the N-terminus"],[2,"number of asparagine in the N-terminus"]]);
			description.append([[1,"pseudo amino acid count of serines in a distance of two within the first 90 amino acids in the N-terminus"],[2,"number of serines in the N-terminus"]]);

			description.append([[2,"maximal pseudo amino acid count of cysteine within the first 120 amino acids in the N-terminus"],[2,"number of cysteine in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of lysine within the the first 120 amino acids in the N-terminus"],[2,"number of lysine in the N-terminus"]]);
			description.append([[2,"number of methionine within the first 120 amino acids in the N-terminus"],[2,"number of methionine in the N-terminus"]]);
			description.append([[2,"number of tryptophane within the first 120 amino acids in the N-terminus"],[2,"number of tryptophane in the N-terminus"]]);

			description.append([[1,"pseudo amino acid count of very hydrophobic residues in a distance of two within the first 130 amino acids in the N-terminus"],[2,"number of very hydrophobic residues in the N-terminus"],[4,"hydrophobic N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of very hydrophobic residues in the last 40 amino acids in the C-terminus"],[2,"number of very hydrophobic residues in the C-terminus"],[4,"hydrophobic C-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of negatively charged residues within the the first 30 amino acids in the N-terminus"],[4,"negatively charged N-terminus"]]);
			description.append([[1,"minimal overall pseudo amino acid count of alpha-helix preferred residues [AFGHKLMNR]"],[2,"number of alpha-helix preferred residues [AFGHKLMNR]"]]);
			description.append([[2,"number of small amino acids within the first 20 amino acids in the N-terminus"],[2,"number of small amino acids in the N-terminus"]]);

			description.append([[1,"maximal normalized overall pseudo amino acid count of aromatic residues"],[1,"number of aromatic residues"]]);
			description.append([[1,"pseudo amino acid count of potentially hydroxylated residues in a distance of three within the last 100 amino acids in the C-terminus"],[2,"number of potentially hydroxylated residues in the C-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of uncharged hydrophobic residues (I,L,V,M,F,Y,W,C,T,A,G) in a distant of two"],[4,"uncharged and hydrophobic"]]);
			description.append([[6,"Prosite pattern <a href=http://www.expasy.ch/prosite/PS00639>PS00639</a> (cysteine proteases) often in lysosomal proteins"],[6,"Prosite pattern <a href=http://www.expasy.ch/prosite/PS00639>PS00639</a> (cysteine proteases)"]]);
			description.append([[6,"Prosite pattern <a href=http://www.expasy.ch/prosite/PS50041>PS50041</a> (C-type lectin domain) present in plasma membrane and extracellular proteins"],[6,"Prosite pattern <a href=http://www.expasy.ch/prosite/PS50041>PS50041</a> (C-type lectin domain)"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 25 in descr.keys() and len(descr[25])>1:
				p = "(" + descr[25] + ") ";
			description.append([[6,"typical extracellular prosite pattern "+p],[6,"typical extracellular prosite pattern"]]);
			p = ""
			if 26 in descr.keys() and len(descr[26])>1:
				p = "(" + descr[26] + ") ";
			description.append([[6,"typical plasma membrane prosite pattern "+p],[6,"typical plasma membrane prosite pattern"]]);
			p = ""
			if 27 in descr.keys() and len(descr[27])>1:
				p = "(" + descr[27] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			if 28 in descr.keys() and len(descr[28])>1:
				p = "(" + descr[28] + ") ";
			description.append([[6,"typical mitochondrial prosite pattern "+p],[6,"typical mitochondrial prosite pattern"]]);
			p = ""
			if 29 in descr.keys() and len(descr[29])>1:
				p = "(" + descr[29] + ") ";

			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);

		if model_name == "multiloc_fungi_dbm":
			description.append([[2,"number of amino acids"],[2,"protein size"]]);
			description.append([[1,"number of weak ER retention signals (number of K,D,E,L in the very C-terminus) at C-terminus times a factor that depends on the fact that no transmembrane helix and no glycolization signal is present. If a strong ER retention signal (KDEL, [KRHQSA][DENQ]EL, DEL, EL) at C-terminus is present it is 100."],[3,"KDEL or derivate in C-terminus"],[7,"ER retention signal"]]);
			description.append([[1,"weighted PTS sum = SKL*0.83+SKF*0.5+[SAGCN][RKH][LIVMAF]*0.25 (presented by Nakai et al.)"],[7,"peroxisomal targeting signal (PTS)"],[7,"peroxisomal targeting signal (PTS)"]]);
			description.append([[1,"product of scores for leucine clusters (1.5 for LL, 2.5 for LLL, 5 for LLLL, 12 for LLLLL, 20 for LLLLLL) in the first 50 amino acids in the N-terminus"],[2,"number of leucine and leucine clusters in the N-terminus"]]);
			description.append([[2,"sum of hydrophobicity of 10 amino acids before the N-terminal cleavage site."],[4,"hydrophobic region before the cleavage site"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[2,"maximal length of very hydrophobic region where the iterative sum of hydrophobicity normed to mean zero drops below zero in at most two cases"],[2,"length of longest very hydrophobic region"]]);
			description.append([[1,"weighted sum of typical amino acids for mitochondrial proteins in N-terminus"],[5,"amino acids for mitochondrial protein in N-terminus"],[6,"putative mitochondrial sorting signal"]]);
			description.append([[1,"autocorrelation of every second hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[7,"secretory pathway sorting signal"]]);
			description.append([[1,"maximal autocorrelation of every sixth charged amino acid within the first 30 amino acids in the N-terminus"],[4,"positive charged amphiphilic helix in the N-terminus"],[7,"mitochondrial targeting peptide"]]);
			description.append([[2,"number of asparagine within the first 70 amino acids in the N-terminus"],[2,"number of asparagine in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of Lysine within the the first 120 amino acids in the N-terminus"],[2,"number of Lysine in the N-terminus"]]);

			description.append([[2,"number of tryptophanee within the first 120 amino acids in the N-terminus"],[2,"number of tryptophane in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of phenylalanine within the the first 130 amino acids in the N-terminus"],[2,"number of phenylalanine in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of cysteine within the the first 130 amino acids in the N-terminus"],[2,"number of cysteine in the N-terminus"]]);
			description.append([[1,"pseudo amino acid count of very hydrophobic residues in a distance of two within the first 130 amino acids in the N-terminus"],[2,"number of very hydrophobic residues in the N-terminus"],[4,"hydrophobic N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of very hydrophobic residues in the last 40 amino acids in the C-terminus"],[2,"number of very hydrophobic residues in the C-terminus"],[4,"hydrophobic C-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of negatively charged residues within the the first 30 amino acids in the N-terminus"],[4,"negatively charged N-terminus"]]);
			description.append([[1,"minimal overall pseudo amino acid count of alpha-helix preferred residues [AFGHKLMNR]"],[2,"number of alpha-helix preferred residues [AFGHKLMNR]"]]);
			description.append([[2,"number of small amino acids within the first 20 amino acids in the N-terminus"],[2,"number of small amino acids in the N-terminus"]]);
			description.append([[1,"pseudo amino acid count of potentially hydroxylated residues in a distance of three within the last 100 amino acids in the C-terminus"],[2,"number of potentially hydroxylated residues in the C-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of uncharged hydrophobic residues (I,L,V,M,F,Y,W,C,T,A,G) in a distant of two"],[4,"uncharged and hydrophobic"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 21 in descr.keys() and len(descr[21])>1:
				if descr[21] != "":
					p = "(" + descr[21] + ") ";
			description.append([[6,"typical plasma membrane prosite pattern "+p],[6,"typical plasma membrane prosite pattern"]]);
			p = ""
			if 22 in descr.keys() and len(descr[22])>1:
				if descr[22] != "":
					p = "(" + descr[22] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			description.append([[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005783>GO term GO:0005783 (endoplasmic reticulum)</a> typical for ER proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005783>GO:0005783 (endoplasmic reticulum)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a> typical for mitochondrial proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a> typical for extracellular proteins"],[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0042025>GO:0042025 (host cell nucleus)</a> typical for viral proteins present in extracellular space and nucleus"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0042025>GO:0042025 (host cell nucleus)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005773>GO:0005773 (vacuole)</a> typical for vacuolar proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005773>GO:0005773 (vacuole)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005778>GO:0005778 (peroxisomal membrane)</a> typical for peroxisomal proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005778>GO:0005778 (peroxisomal membrane)</a>"]]);
			p = ""
			if 29 in descr.keys() and len(descr[29])>1:
				if descr[29] != "":
					p = "(GO terms: " + str(descr[29]) + ")";
			description.append([[9,"typical GO term "+p],[9,"typical GO term"]]);

		if model_name == "multiloc_fungi_dbm_p+p":
			description.append([[2,"number of amino acids"],[2,"protein size"]]);
			description.append([[1,"number of weak ER retention signals (number of K,D,E,L in the very C-terminus) at C-terminus times a factor that depends on the fact that no transmembrane helix and no glycolization signal is present. If a strong ER retention signal (KDEL, [KRHQSA][DENQ]EL, DEL, EL) at C-terminus is present it is 100."],[3,"KDEL or derivate in C-terminus"],[7,"ER retention signal"]]);
			description.append([[1,"weighted PTS sum = SKL*0.83+SKF*0.5+[SAGCN][RKH][LIVMAF]*0.25 (presented by Nakai et al.)"],[7,"peroxisomal targeting signal (PTS)"],[7,"peroxisomal targeting signal (PTS)"]]);
			description.append([[1,"product of scores for leucine clusters (1.5 for LL, 2.5 for LLL, 5 for LLLL, 12 for LLLLL, 20 for LLLLLL) in the first 50 amino acids in the N-terminus"],[2,"number of leucine and leucine clusters in the N-terminus"]]);
			description.append([[2,"sum of hydrophobicity of 10 amino acids before the N-terminal cleavage site."],[4,"hydrophobic region before the cleavage site"],[7,"alternative secretory pathway sorting signal"]]);

			description.append([[2,"maximal length of very hydrophobic region where the iterative sum of hydrophobicity normed to mean zero drops below zero in at most two cases"],[2,"length of longest very hydrophobic region"]]);
			description.append([[1,"weighted sum of typical amino acids for mitochondrial proteins in N-terminus"],[5,"amino acids for mitochondrial protein in N-terminus"],[6,"putative mitochondrial sorting signal"]]);
			description.append([[1,"autocorrelation of every second hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[7,"secretory pathway sorting signal"]]);
			description.append([[1,"maximal autocorrelation of every sixth charged amino acid within the first 30 amino acids in the N-terminus"],[4,"positive charged amphiphilic helix in the N-terminus"],[7,"mitochondrial targeting peptide"]]);
			description.append([[2,"number of asparagine within the first 70 amino acids in the N-terminus"],[2,"number of asparagine in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of tryptophane within the the first 120 amino acids in the N-terminus"],[2,"number of tryptophane in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of lysine within the the first 120 amino acids in the N-terminus"],[2,"number of lysine in the N-terminus"]]);

			description.append([[2,"number of methionine within the first 120 amino acids in the N-terminus"],[2,"number of methionine in the N-terminus"]]);
			description.append([[2,"number of tryptophanee within the first 120 amino acids in the N-terminus"],[2,"number of tryptophane in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of phenylalanine within the the first 130 amino acids in the N-terminus"],[2,"number of phenylalanine in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of cysteine within the the first 130 amino acids in the N-terminus"],[2,"number of cysteine in the N-terminus"]]);

			description.append([[1,"pseudo amino acid count of very hydrophobic residues in a distance of two within the first 130 amino acids in the N-terminus"],[2,"number of very hydrophobic residues in the N-terminus"],[4,"hydrophobic N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of very hydrophobic residues in the last 40 amino acids in the C-terminus"],[2,"number of very hydrophobic residues in the C-terminus"],[4,"hydrophobic C-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of negatively charged residues within the the first 30 amino acids in the N-terminus"],[4,"negatively charged N-terminus"]]);
			description.append([[1,"minimal overall pseudo amino acid count of alpha-helix preferred residues [AFGHKLMNR]"],[2,"number of alpha-helix preferred residues [AFGHKLMNR]"]]);
			description.append([[2,"number of small amino acids within the first 20 amino acids in the N-terminus"],[2,"number of small amino acids in the N-terminus"]]);

			description.append([[1,"maximal pseudo amino acid count of aromatic residues within the last 100 residues in the C-terminus"],[1,"number of aromatic residues in the C-terminus"]]);
			description.append([[1,"pseudo amino acid count of potentially hydroxylated residues in a distance of three within the last 100 amino acids in the C-terminus"],[2,"number of potentially hydroxylated residues in the C-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of uncharged hydrophobic residues (I,L,V,M,F,Y,W,C,T,A,G) in a distant of two"],[4,"uncharged and hydrophobic"]]);

			description.append([[6,"Prosite pattern <a href=http://www.expasy.ch/prosite/PS50041>PS50041</a> (C-type lectin domain) present in plasma membrane and extracellular proteins"],[6,"Prosite pattern <a href=http://www.expasy.ch/prosite/PS50041>PS50041</a> (C-type lectin domain)"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 25 in descr.keys() and len(descr[25])>1:
				p = "(" + descr[25] + ") ";
			description.append([[6,"typical extracellular prosite pattern "+p],[6,"typical extracellular prosite pattern"]]);
			p = ""
			if 26 in descr.keys() and len(descr[26])>1:
				p = "(" + descr[26] + ") ";
			description.append([[6,"typical plasma membrane prosite pattern "+p],[6,"typical plasma membrane prosite pattern"]]);
			p = ""
			if 27 in descr.keys() and len(descr[27])>1:
				p = "(" + descr[27] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			if 28 in descr.keys() and len(descr[28])>1:
				p = "(" + descr[28] + ") ";
			description.append([[6,"typical mitochondrial prosite pattern "+p],[6,"typical mitochondrial prosite pattern"]]);
			p = ""
			if 29 in descr.keys() and len(descr[29])>1:
				p = "(" + descr[29] + ") ";

			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);

		if model_name == "multiloc_plants_dbm":
			# plus code 1 = low/medium/high       2=small/average/large         3=present/no present    4=barely/medium/very
			# 5 = typical/unusal   6=no present/present  7=no/weak/strong
			# 8=negative/neutral/positive
			# 9=go clusters
			description.append([[2,"number of amino acids"],[2,"protein size"]]);
			description.append([[1,"number of weak ER retention signals (number of K,D,E,L in the very C-terminus) at C-terminus times a factor that depends on the fact that no transmembrane helix and no glycolization signal is present. If a strong ER retention signal (KDEL, [KRHQSA][DENQ]EL, DEL, EL) at C-terminus is present it is 100."],[3,"KDEL or derivate in C-terminus"],[7,"ER retention signal"]]);
			description.append([[2,"sum of hydrophobicity of 10 amino acids before the N-terminal cleavage site."],[4,"hydrophobic region before the cleavage site"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[2,"maximal length of very hydrophobic region where the iterative sum of hydrophobicity normed to mean zero drops below zero in at most two cases"],[2,"length of longest very hydrophobic region"]]);
			description.append([[1,"weighted sum of typical amino acids for mitochondrial proteins in N-terminus"],[5,"amino acids for mitochondrial protein in N-terminus"],[6,"putative mitochondrial sorting signal"]]);
			description.append([[1,"autocorrelation of every hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[6,"secretory pathway sorting signal"]]);
			description.append([[1,"sum of charge of first 10 amino acids in N-terminus"],[7,"mitochondrial targeting peptide"]]);
			description.append([[1,"autocorrelation of every sixth charged amino acid within the first 50 amino acids in the N-terminus"],[6,"charged N-terminus"],[6,"putative mitochondrial or chloroplast targeting signal"]]);
			description.append([[2,"number of Serine within the first 50 amino acids in the N-terminus"],[2,"number of Serine in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of leucine within the first 80 amino acids in the N-terminus"],[1,"number of leucine in the N-terminus"]]);

			description.append([[1,"maximal pseudo amino acid count of cysteine within the the first 120 amino acids in the N-terminus"],[2,"number of cysteine in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of lysine within the the first 120 amino acids in the N-terminus"],[2,"number of lysine in the N-terminus"]]);
			description.append([[2,"number of phenylalanine within the first 120 amino acids in the N-terminus"],[2,"number of phenylalanine in the N-terminus"]]);
			description.append([[2,"number of tryptophane within the first 120 amino acids in the N-terminus"],[2,"number of tryptophane in the N-terminus"]]);
			description.append([[1,"pseudo amino acid count of very hydrophobic residues in a distance of two within the first 100 amino acids in the N-terminus"],[2,"number of very hydrophobic residues in the N-terminus"],[7,"alternative secretory pathway sorting signal"]]);

			description.append([[1,"maximal normalized pseudo amino acid count of very hydrophobic residues in the last 40 amino acids in the C-terminus"],[2,"number of very hydrophobic residues in the C-terminus"],[4,"hydrophobic C-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of negatively charged residues within the the first 30 amino acids in the N-terminus"],[4,"negatively charged N-terminus"]]);
			description.append([[1,"pseudo amino acid count of alpha-helix preferred residues [AFGHKLMNR] in a distance of three in the first 20 amino acids"],[2,"number of alpha-helix preferred residues [AFGHKLMNR] in the N-terminus"]]);
			description.append([[1,"minimal overall pseudo amino acid count of alpha-helix preferred residues [AFGHKLMNR]"],[2,"number of alpha-helix preferred residues [AFGHKLMNR]"]]);
			description.append([[2,"number of large amino acids within the first 20 amino acids in the N-terminus"],[2,"number of large amino acids in the N-terminus"]]);
			description.append([[1,"maximal overall pseudo amino acid count of potentially hydroxylated residues"],[1,"number of potentially hydroxylated residues"]]);
			description.append([[1,"overall pseudo amino acid count of uncharged hydrophobic residues (I,L,V,M,F,Y,W,C,T,A,G) in a distant of two"],[4,"uncharged and hydrophobic"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 22 in descr.keys() and len(descr[22])>1:
				if descr[22] != "":
					p = "(" + descr[22] + ") ";
			description.append([[6,"typical for plasma membrane prosite pattern "+p],[6,"typical for plasma membrane prosite pattern"]]);
			p = ""
			if 23 in descr.keys() and len(descr[23])>1:
				if descr[23] != "":
					p = "(" + descr[23] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			description.append([[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005783>GO term GO:0005783 (endoplasmic reticulum)</a> typical for ER proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005783>GO:0005783 (endoplasmic reticulum)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a> typical for mitochondrial proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005739>GO:0005739 (mitochondrion)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a> typical for extracellular proteins"],[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0042025>GO:0042025 (host cell nucleus)</a> typical for viral proteins present in extracellular space and nucleus"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0042025>GO:0042025 (host cell nucleus)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0009507>GO:0009507 (chloroplast)</a> typical for chloroplast proteins"],[6,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0009507>GO:0009507 (chloroplast)</a>"]]);
			p = ""
			if 29 in descr.keys() and len(descr[29])>1:
				if descr[29] != "":
					p = "(GO terms: " + str(descr[29]) + ")";
			description.append([[9,"typical GO term "+p],[9,"typical GO term"]]);

		if model_name == "multiloc_plants_dbm_p+p":
			# plus code 1 = low/medium/high       2=small/average/large         3=present/no present    4=barely/medium/very
			# 5 = typical/unusal   6=no present/present  7=no/weak/strong
			# 8=negative/neutral/positive
			# 9=go clusters
			description.append([[2,"number of amino acids"],[2,"protein size"]]);
			description.append([[1,"number of weak ER retention signals (number of K,D,E,L in the very C-terminus) at C-terminus times a factor that depends on the fact that no transmembrane helix and no glycolization signal is present. If a strong ER retention signal (KDEL, [KRHQSA][DENQ]EL, DEL, EL) at C-terminus is present it is 100."],[3,"KDEL or derivate in C-terminus"],[7,"ER retention signal"]]);
			description.append([[1,"weighted PTS sum = SKL*0.83+SKF*0.5+[SAGCN][RKH][LIVMAF]*0.25 (presented by Nakai et al.)"],[7,"peroxisomal targeting signal (PTS)"],[7,"peroxisomal targeting signal (PTS)"]]);
			description.append([[2,"sum of hydrophobicity of 10 amino acids before the N-terminal cleavage site."],[4,"hydrophobic region before the cleavage site"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[2,"maximal length of very hydrophobic region where the iterative sum of hydrophobicity normed to mean zero drops below zero in at most two cases"],[2,"length of longest very hydrophobic region"]]);
			description.append([[1,"weighted sum of typical amino acids for mitochondrial proteins in N-terminus"],[5,"amino acids for mitochondrial protein in N-terminus"],[6,"putative mitochondrial sorting signal"]]);

			description.append([[6,"chloroplast motif [MA] at the beginning of the protein"],[6,"typical chloroplast protein N-terminus (MA motif)"]]);
			description.append([[1,"autocorrelation of every hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic N-terminus"],[6,"secretory pathway sorting signal"]]);
			description.append([[1,"sum of charge of first 10 amino acids in N-terminus"],[7,"mitochondrial targeting peptide"]]);
			description.append([[1,"autocorrelation of every sixth charged amino acid within the first 50 amino acids in the N-terminus"],[6,"charged N-terminus"],[6,"putative mitochondrial or chloroplast targeting signal"]]);
			description.append([[2,"number of Serine within the first 50 amino acids in the N-terminus"],[2,"number of Serine in the N-terminus"]]);

			description.append([[1,"maximal normalized pseudo amino acid count of leucine within the first 80 amino acids in the N-terminus"],[1,"number of leucine in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of cysteine within the the first 120 amino acids in the N-terminus"],[2,"number of cysteine in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of lysine within the the first 120 amino acids in the N-terminus"],[2,"number of lysine in the N-terminus"]]);
			description.append([[2,"number of phenylalanine within the first 120 amino acids in the N-terminus"],[2,"number of phenylalanine in the N-terminus"]]);
			description.append([[2,"number of tryptophane within the first 120 amino acids in the N-terminus"],[2,"number of tryptophane in the N-terminus"]]);

			description.append([[1,"pseudo amino acid count of very hydrophobic residues in a distance of two within the first 100 amino acids in the N-terminus"],[2,"number of very hydrophobic residues in the N-terminus"],[7,"alternative secretory pathway sorting signal"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of very hydrophobic residues in the last 40 amino acids in the C-terminus"],[2,"number of very hydrophobic residues in the C-terminus"],[4,"hydrophobic C-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of negatively charged residues within the the first 30 amino acids in the N-terminus"],[4,"negatively charged N-terminus"]]);
			description.append([[1,"pseudo amino acid count of alpha-helix preferred residues [AFGHKLMNR] in a distance of three in the first 20 amino acids"],[2,"number of alpha-helix preferred residues [AFGHKLMNR] in the N-terminus"]]);
			description.append([[1,"minimal overall pseudo amino acid count of alpha-helix preferred residues [AFGHKLMNR]"],[2,"number of alpha-helix preferred residues [AFGHKLMNR]"]]);

			description.append([[2,"number of large amino acids within the first 20 amino acids in the N-terminus"],[2,"number of large amino acids in the N-terminus"]]);
			description.append([[1,"maximal overall pseudo amino acid count of potentially hydroxylated residues"],[1,"number of potentially hydroxylated residues"]]);
			description.append([[1,"overall pseudo amino acid count of uncharged hydrophobic residues (I,L,V,M,F,Y,W,C,T,A,G) in a distant of two"],[4,"uncharged and hydrophobic"]]);

			descr = (features.getDescriptions())[nr];
			p = ""
			if 24 in descr.keys() and len(descr[24])>1:
				p = "(" + descr[24] + ") ";
			description.append([[6,"typical extracellular prosite pattern "+p],[6,"typical extracellular prosite pattern"]]);

			p = ""
			if 25 in descr.keys() and len(descr[25])>1:
				p = "(" + descr[25] + ") ";
			description.append([[6,"typical plasma membrane prosite pattern "+p],[6,"typical plasma membrane prosite pattern"]]);
			p = ""
			if 26 in descr.keys() and len(descr[26])>1:
				p = "(" + descr[26] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nuclear prosite pattern"]]);
			p = ""
			if 27 in descr.keys() and len(descr[27])>1:
				p = "(" + descr[27] + ") ";
			description.append([[6,"typical mitochondrial prosite pattern "+p],[6,"typical mitochondrial prosite pattern"]]);
			p = ""
			if 28 in descr.keys() and len(descr[28])>1:
				p = "(" + descr[28] + ") ";
			description.append([[6,"typical chloroplast prosite pattern "+p],[6,"typical chloroplast prosite pattern"]]);
			p = ""
			if 29 in descr.keys() and len(descr[29])>1:
				p = "(" + descr[29] + ") ";
			description.append([[6,"typical cytoplasmic prosite pattern "+p],[6,"typical cytoplasmic prosite pattern"]]);


		if model_name == "suba_plus":
		# plus code 1 = low/medium/high       2=small/average/large         3=present/no present    4=barely/medium/very
		# 5 = typical/unusal   6=no present/present  7=no/weak/strong
		# 8=negative/neutral/positive 10=no/weak/strong
		# 9=go clusters
			description.append([[2,"number of amino acids"],[2,"protein size"]]);
			description.append([[1,"autocorrelation of every hydrophobic amino acid within the first 20 amino acids in the N-terminus"],[4,"hydrophobic very N-terminus"],[7,"secretory pathway sorting signal"]]);
			description.append([[1,"overall autocorrelation of every hydrophobic amino acid"],[4,"hydrophobic protein"],[2,"hydrophobic regions (possibly transmembrane regions)"]]);
			description.append([[1,"autocorrelation of every charged amino acid within the first 120 amino acids in the N-terminus"],[4,"positive charged N-terminus"],[6,"chloroplast or mitochondrial targeting peptide"]]);
			description.append([[2,"number of Alanine within the first 10 amino acids in the N-terminus"],[2,"number of Alanine in the N-terminus"]]);
			description.append([[2,"number of Lysine within the first 10 amino acids in the N-terminus"],[2,"number of Lysine in the N-terminus"]]);
			description.append([[2,"number of Glycine within the first 20 amino acids in the N-terminus"],[2,"number of Glycine in the N-terminus"]]);
			description.append([[2,"number of Tyrosine within the first 70 amino acids in the N-terminus"],[2,"number of Tyrosine in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of Tryptophane within the the first 80 amino acids in the N-terminus"],[2,"number of Tryptophane in the N-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of Glutamine within the the first 80 amino acids in the N-terminus"],[2,"number of Glutamine in the N-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of Serine within the the first 190 amino acids in the N-terminus"],[2,"number of Serine in the N-terminus"]]);

			description.append([[2,"number of Cysteine within the first 190 amino acids in the N-terminus"],[2,"number of Cysteine in the N-terminus"]]);
			description.append([[2,"number of Serine within the last 80 amino acids in the C-terminus"],[2,"number of Serine in the C-terminus"]]);
			description.append([[1,"pseudo amino acid count of Proline in a distance of four within the last 70 amino acids in the C-terminus"],[2,"number of Proline in the C-terminus"]]);
			description.append([[1,"maximal pseudo amino acid count of Leucine within the the last 60 amino acids in the C-terminus"],[2,"number of Leucine in the C-terminus"]]);
			description.append([[2,"number of Isoleucine within the last 30 amino acids in the C-terminus"],[2,"number of Isoleucine in the C-terminus"]]);
			description.append([[1,"maximal normalized pseudo amino acid count of Glutamate within the the last 20 amino acids in the C-terminus"],[2,"number of Glutamate in the C-terminus"]]);
			description.append([[2,"number of negatively charged residues within the first 40 amino acids in the N-terminus"],[4,"negatively charged  N-terminus"]]);
			description.append([[1,"overall pseudo amino acid count of alpha helix preferred amino acids in a distance of three"],[2,"number of alpha helix preferred acids"]]);
			description.append([[1,"maximal pseudo amino acid count of small residues within the first 40 amino acids in the N-terminus"],[2,"number of small amino acids in the N-terminus"]]);

			description.append([[1,"overall pseudo amino acid count of large amino acids in a distance of three"],[2,"number of large amino acids"]]);
			description.append([[1,"sequence similarity to a protein in the Anabaena variabilis genome"],[6,"homolog in Anabaena variabilis"]]);
			description.append([[1,"sequence similarity to a protein in the Oryza sativa Japonica genome"],[6,"homolog in Oryza sativa Japonica"]]);
			descr = (features.getDescriptions())[nr];
			p = ""
			description.append([[6,"CRAL-TRIO lipid binding domain profile prosite pattern: constitute a hydrophobic lipid binding pocket (typical for cytoplasmic and plasma membrane proteins)"],[6,"prosite pattern PS50191"]]);
			description.append([[6,"Solute carrier (Solcar) repeat profile prosite pattern: present in substrate carrier proteins involved in energy transfer in the inner mitochondrial membrane (typical for mitochondrial proteins)"],[6,"prosite pattern PS50920"]]);
			p = ""
			if 25 in descr.keys() and len(descr[25])>1:
				if descr[25] != "":
					p = "(" + descr[25] + ") ";
			description.append([[6,"typical nuclear prosite pattern "+p],[6,"typical nucler prosite pattern"]]);
			p = ""
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005789>GO:0005789 (endoplasmic reticulum membranen)</a> typical for endoplasmatic and transmembrane proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005789>GO:0005789 (endoplasmic reticulum membrane)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a> typical for extracellular proteins"],[14,"<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0005576>GO:0005576 (extracellular region)</a>"]]);
			description.append([[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0004298>GO:0004298 (threonine-type endopeptidase activity)</a> typical for cytoplasmic and nuclear proteins"],[14,"GO term <a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:0004298>GO:0004298 (threonine-type endopeptidase activity)</a>"]]);
			p = ""
			if 29 in descr.keys() and len(descr[29])>1:
				if descr[29] != "":
					p = "(GO terms: " + str(descr[29]) + ")";
			description.append([[9,"typical GO term "+p],[9,"typical GO term"]]);


		return description;

	def getQuantityDescriptors(self,model_name):
		quantity = [];
		quantity.append([]);
		quantity.append(["low", "medium", "high"]);
		quantity.append(["small", "average", "large"]);
		quantity.append(["present", "no"]);
		quantity.append(["barely", "medium", "very"]);
		quantity.append(["typical", "unusal"]);
		quantity.append(["no", "present"]);
		quantity.append(["no", "weak", "strong"]);
		quantity.append(["negative", "neutral", "positive"]);
		quantity.append(self.getGOIntervals(model_name, False));
		quantity.append(["strong", "weak", "no"]);
		#special case where discretization is somehow weird
		quantity.append(["no", "no","present","no"]);
		quantity.append(["no", "no","weak","strong","no"]); #12
		quantity.append(["no", "no","weak","strong"]); #13
		quantity.append(["absent", "present"]); #14

		return quantity;

	# save raw model and return attribute distributions
	def loadRawModel(self, model_name):
		f = open(config.PATH_MODELS+model_name+".model_raw", "r");
		line = f.readline();
		model_distributions = [];

		while line:
			if line[0] == "*":
				line = f.readline();
				attribute = [];
				while line and line[0] != "*":
					line_s = line.strip(" \n\r\t").split(" ");
					attribute.append(line_s);
					line = f.readline();
				model_distributions.append(attribute);
			else:
				line = f.readline();

		return model_distributions;

	def getAttributeDistributions(self, external_model_name, a_index):
		if external_model_name in self.__external_model_names:
			pos = self.__external_model_names.index(external_model_name);
			model_name = self.__model_names[pos];
			model_distribution = self.loadRawModel(model_name);
			return model_distribution[a_index];
		else:
			return "ERROR!";
