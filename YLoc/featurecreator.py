from aasequences import *;
import types;
from sets import Set;
import random;
import re;
from math import sqrt, exp;

from sortsignals import *;
from go2 import *;
from prosite import *;

#####################################
######### Utils #####################
#####################################

def uniqify(seq): 
	# order preserving
	def idfun(x): return x
	seen = {}
	result = []
	for item in seq:
		marker = idfun(item)
		if marker in seen: continue	
		seen[marker] = 1
		result.append(item)
	return result
	


class FeatureCreator(object):

	"Class is initialized with a set of sequences and its methods create features."

#####################################
############ alphabet names #########
#####################################

	ALPHABETS = [ "aminoacid", "hydrophob", "polar", "charged", "structure" , "size", "aromatic", "hydroxylated" ,"aliphatic_helix"];
	PROPERTIES = [ "volume" , "hydrophob" , "charge" ];

#####################################
######### constructor ###############
#####################################

	def __init__(self, aa_sequences = None, feature_vectors=None, response = None, labels = None, id = "1" , description = None):
		if aa_sequences == None:
			aa_sequences = AASequences();
		self.__aa_sequences = aa_sequences;
		self.__feature_vectors = feature_vectors;
		self.__response = response;
		self.__description = description;
		if response != None:
			self.__response_classes = uniqify(response);
		else:
			self.__response_classes = None;
		self.__labels = labels;
		if feature_vectors is None:
			self.__feature_vectors = [];
			self.__labels = [];
			for sequence_nr in range(aa_sequences.size()):
				self.__feature_vectors.append([]);
		if description is None:
			self.__description = [];
			for sequence_nr in range(aa_sequences.size()):
				self.__description.append({});
		self.process_id = id;


#####################################
######### Operator ##################
#####################################

	def __add__(self, x):
		if self.__response == None:
			new_response = None;
		else:
			new_response = self.__response + x.__response;
		return FeatureCreator(self.__aa_sequences + x.__aa_sequences, self.__feature_vectors + x.__feature_vectors, new_response , self.__labels, self.process_id, self.__description + x.__description);
        
	def size(self):
		return self.__aa_sequences.size();
	
	def getDescriptions(self):
		return self.__description;
	
	# use with caution!
	def getFeature(self, index):
		f_list = [];
		for i in range(self.__aa_sequences.size()):
			f_list.append(self.__feature_vectors[i][index]);
		return f_list;
		
	def setSequences(self, aa_sequences):
		self.__aa_sequences = aa_sequences;
		
		
#####################################
######## alphabet access ############
#####################################

	def __get_properties(self, property_name):
		if (property_name == "volume"):
			properties = {"A" : 88.6, "C" : 108.5, "D" : 111.1, "E" : 138.4, "F" : 189.9,
				"G" : 60.1, "H" : 153.2, "I" : 166.7, "K" : 168.6, "L" : 166.7,
				"M" : 162.9, "N" : 114.1, "P" : 112.7, "Q" : 143.8, "R" : 173.4,
				"S" : 89.0, "T" : 116.1, "V" : 140.0, "W" : 227.8, "Y" : 193.6};
			return properties;
		if (property_name == "charge"):
			properties = {"A" : 6.107, "C" : 5.02, "D" : 2.98, "E" : 3.08, "F" : 5.91,
				"G" : 6.064, "H" : 7.64, "I" : 6.038, "K" : 9.47, "L" : 6.038,
				"M" : 5.74, "N" : 5.41, "P" : 6.3, "Q" : 5.65, "R" : 10.76,
				"S" : 5.68, "T" : 5.64, "V" : 6.002, "W" : 5.88, "Y" : 5.63};
			return properties;
		if (property_name == "hydrophob"):
			properties = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880};
			return properties;


		raise "FeatureCreator: get_properties: unknown property name";

	def __get_alphabet(self, alphabet_name):
		if (alphabet_name=="aminoacid"):
			# 2 very hydrophobic, 1 slighly hydrophobic
			alphabet = {"A" : "A", "C" : "C", "D" : "D", "E" : "E", "F" : "F",
				"G" : "G", "H" : "H", "I" : "I", "K" : "K", "L" : "L",
				"M" : "M", "N" : "N", "P" : "P", "Q" : "Q", "R" : "R",
				"S" : "S", "T" : "T", "V" : "V", "W" : "W", "Y" : "Y"};
			return alphabet;
		if (alphabet_name=="hydrophob"):
			# 2 very hydrophobic, 1 slighly hydrophobic
			alphabet = {"A" : "1", "C" : "2", "D" : "0", "E" : "0", "F" : "2",
				"G" : "1", "H" : "1", "I" : "2", "K" : "0", "L" : "2",
				"M" : "2", "N" : "0", "P" : "1", "Q" : "0", "R" : "0",
				"S" : "1", "T" : "1", "V" : "2", "W" : "2", "Y" : "1"};
			return alphabet;
		if (alphabet_name=="polar"):
			# 1 polar residues
			alphabet = {"A" : "0", "C" : "0", "D" : "1", "E" : "1", "F" : "0",
				"G" : "0", "H" : "1", "I" : "0", "K" : "1", "L" : "0",
				"M" : "0", "N" : "1", "P" : "0", "Q" : "1", "R" : "1",
				"S" : "0", "T" : "0", "V" : "0", "W" : "0", "Y" : "0"};
			return alphabet;
		if (alphabet_name=="charged"):
			# 1 negative, 2 positive
			alphabet = {"A" : "0", "C" : "0", "D" : "1", "E" : "1", "F" : "0",
				"G" : "0", "H" : "2", "I" : "0", "K" : "2", "L" : "0",
				"M" : "0", "N" : "0", "P" : "0", "Q" : "0", "R" : "2",
				"S" : "0", "T" : "0", "V" : "0", "W" : "0", "Y" : "0"};
			return alphabet;
		if (alphabet_name=="structure"):
			# 1 helix, 2 barrel 
			alphabet = {"A" : "1", "C" : "2", "D" : "0", "E" : "0", "F" : "1",
				"G" : "1", "H" : "1", "I" : "2", "K" : "1", "L" : "1",
				"M" : "1", "N" : "1", "P" : "0", "Q" : "0", "R" : "1",
				"S" : "0", "T" : "2", "V" : "2", "W" : "2", "Y" : "2"};
			return alphabet;
		if (alphabet_name=="size"):
			# 0 tiny, 1 small, 2 other 
			alphabet = {"A" : "0", "C" : "1", "D" : "1", "E" : "2", "F" : "2",
				"G" : "0", "H" : "2", "I" : "2", "K" : "2", "L" : "2",
				"M" : "2", "N" : "1", "P" : "1", "Q" : "2", "R" : "2",
				"S" : "0", "T" : "1", "V" : "1", "W" : "2", "Y" : "2"};
			return alphabet;
		if (alphabet_name=="aromatic"):
			# 0 aromatic, 1 other
			alphabet = {"A" : "1", "C" : "1", "D" : "1", "E" : "1", "F" : "0",
				"G" : "1", "H" : "0", "I" : "1", "K" : "1", "L" : "1",
				"M" : "1", "N" : "1", "P" : "1", "Q" : "1", "R" : "1",
				"S" : "1", "T" : "1", "V" : "1", "W" : "0", "Y" : "0"};
			return alphabet;
		if (alphabet_name=="hydroxylated"):
			# 0 hydroxy, 1 other
			alphabet = {"A" : "0", "C" : "0", "D" : "0", "E" : "0", "F" : "0",
				"G" : "0", "H" : "0", "I" : "0", "K" : "0", "L" : "1",
				"M" : "0", "N" : "0", "P" : "1", "Q" : "0", "R" : "1",
				"S" : "1", "T" : "0", "V" : "0", "W" : "0", "Y" : "1"};
			return alphabet;
		if (alphabet_name=="aliphatic_helix"):
			# 0 basic [RK}
			# 1 hydrophobic uncharged [ILVMFYWCTAG]
			# 2 others
			alphabet = {"A" : "1", "C" : "1", "D" : "2", "E" : "2", "F" : "1",
				"G" : "1", "H" : "2", "I" : "1", "K" : "0", "L" : "1",
				"M" : "1", "N" : "2", "P" : "2", "Q" : "2", "R" : "0",
				"S" : "2", "T" : "1", "V" : "1", "W" : "1", "Y" : "1"};
			return alphabet;

		if (alphabet_name=="glutamate"):
			alphabet = {"E" : "E"};
			return alphabet;
			
		if (alphabet_name=="leucine"):
			alphabet = {"L" : "L"};
			return alphabet;
		
		if (alphabet_name=="isoleucine"):
			alphabet = {"I" : "I"};
			return alphabet;
			
		if (alphabet_name=="alanine"):
			alphabet = {"A" : "A"};
			return alphabet;
		
		if (alphabet_name=="aspartate"):
			alphabet = {"D" : "D"};
			return alphabet;
			
		if (alphabet_name=="arginine"):
			alphabet = {"R" : "R"};
			return alphabet;
		
		if (alphabet_name=="tyrosine"):
			alphabet = {"Y" : "Y"};
			return alphabet;
		
		if (alphabet_name=="lysine"):
			alphabet = {"K" : "K"};
			return alphabet;
		
		if (alphabet_name=="phenylalanine"):
			alphabet = {"F" : "F"};
			return alphabet;
		
		if (alphabet_name=="methionine"):
			alphabet = {"M" : "M"};
			return alphabet;
		
		if (alphabet_name=="asparagine"):
			alphabet = {"N" : "N"};
			return alphabet;
	
		if (alphabet_name=="cysteine"):
			alphabet = {"C" : "C"};
			return alphabet;
	
		if (alphabet_name=="histidine"):
			alphabet = {"H" : "H"};
			return alphabet;
	
		if (alphabet_name=="glutamine"):
			alphabet = {"Q" : "Q"};
			return alphabet;
	
		if (alphabet_name=="proline"):
			alphabet = {"P" : "P"};
			return alphabet;
			
		if (alphabet_name=="serine"):
			alphabet = {"S" : "S"};
			return alphabet;
			
		if (alphabet_name=="valine"):
			alphabet = {"V" : "V"};
			return alphabet;
		
		if (alphabet_name=="threonine"):
			alphabet = {"T" : "T"};
			return alphabet;
			
		if (alphabet_name=="tryptophane"):
			alphabet = {"W" : "W"};
			return alphabet;
	
		if (alphabet_name=="structure_none"):
			alphabet = {"D" : "0", "E" : "0", "P" : "0", "Q" : "0", "S" : "0"};
			return alphabet
			
		if (alphabet_name=="structure_helix"):
			alphabet = {"A" : "1",  "F" : "1","G" : "1", "H" : "1", "K" : "1", "L" : "1",
				"M" : "1", "N" : "1", "R" : "1"};
			return alphabet
		
		if (alphabet_name=="structure_sheet"):
			# 1 helix, 2 barrel
			alphabet = {"C" : "2", "I" : "2", "T" : "2", "V" : "2", "W" : "2", "Y" : "2"};
			return alphabet;
	
		if (alphabet_name=="unpolar"):
			# 1 polar residues
			alphabet = {"A" : "0", "C" : "0", "F" : "0",
				"G" : "0", "I" : "0", "L" : "0",
				"M" : "0", "P" : "0", 
				"S" : "0", "T" : "0", "V" : "0", "W" : "0", "Y" : "0"};
			return alphabet
		
		if (alphabet_name=="polar_polar"):
			# 1 polar residues
			alphabet = {"D" : "1", "E" : "1", "H" : "1",  "K" : "1", 
				"N" : "1", "Q" : "1", "R" : "1"};
			return alphabet
	
		if (alphabet_name=="charge_none"):
			# 1 negative, 2 positive
			alphabet = {"A" : "0", "C" : "0",  "F" : "0",
				"G" : "0", "I" : "0", "L" : "0",
				"M" : "0", "N" : "0", "P" : "0", "Q" : "0", 
				"S" : "0", "T" : "0", "V" : "0", "W" : "0", "Y" : "0"};
			return alphabet
		
		if (alphabet_name=="charge_negative"):
			# 1 negative, 2 positive
			alphabet = {"D" : "1", "E" : "1"};
			return alphabet;
		
		if (alphabet_name=="charge_positive"):
			# 1 negative, 2 positive
			alphabet = {"H" : "2", "K" : "2", "R" : "2" };
			return alphabet
		
		if (alphabet_name=="size_tiny"):
			alphabet = {"A" : "0", 
				"G" : "0",
				"S" : "0"};
			return alphabet;
		
		if (alphabet_name=="size_small"):
			# 0 tiny, 1 small, 2 other
			alphabet = {"C" : "1", "D" : "1", 
				 "N" : "1", "P" : "1", 
				 "T" : "1", "V" : "1"};
			return alphabet;
		
		if (alphabet_name=="size_large"):
			# 0 tiny, 1 small, 2 other
			alphabet = {"E" : "2", "F" : "2",
				"H" : "2", "I" : "2", "K" : "2", "L" : "2",
				"M" : "2", "Q" : "2", "R" : "2",
				"W" : "2", "Y" : "2"};
			return alphabet;
		
		if (alphabet_name=="hydrophob_very"):
			alphabet = {"C" : "2", "F" : "2",
				"I" : "2", "L" : "2",
				"M" : "2", "V" : "2", "W" : "2"};
			return alphabet;
		
		if (alphabet_name=="hydrophob_slightly"):
			alphabet = {"A" : "1", 
				"G" : "1", "H" : "1", 
				"P" : "1", 
				"S" : "1", "T" : "1", "Y" : "1"};
			return alphabet;
		
		if (alphabet_name=="hydrophob_none"):
			alphabet = {"D" : "0", "E" : "0", "K" : "0", "N" : "0", "Q" : "0", 
			"R" : "0"};
			return alphabet;
			
		if (alphabet_name=="hydroxylated_hydroxy"):
			# 0 hydroxy, 1 other
			alphabet = {"L" : "1", "P" : "1", "R" : "1","S" : "1", "Y" : "1"};
			return alphabet;
			
		if (alphabet_name=="hydroxylated_none"):
			# 0 hydroxy, 1 other
			alphabet = {"A" : "0", "C" : "0", "D" : "0", "E" : "0", "F" : "0",
				"G" : "0", "H" : "0", "I" : "0", "K" : "0",
				"M" : "0", "N" : "0", "Q" : "0",
				"T" : "0", "V" : "0", "W" : "0"};
			return alphabet;
			
		if (alphabet_name=="aromatic_none"):
			# 0 aromatic, 1 other
			alphabet = {"A" : "1", "C" : "1", "D" : "1", "E" : "1", 
				"G" : "1", "I" : "1", "K" : "1", "L" : "1",
				"M" : "1", "N" : "1", "P" : "1", "Q" : "1", "R" : "1",
				"S" : "1", "T" : "1", "V" : "1"};
			return alphabet;
		
		if (alphabet_name=="aromatic_aromatic"):
			alphabet = {"F" : "0", "H" : "0", "W" : "0", "Y" : "0"};
			return alphabet;
	
		if (alphabet_name=="aliphatic_helix_basic"):
			# 0 basic [RK}
			# 1 hydrophobic uncharged [ILVMFYWCTAG]
			# 2 others
			alphabet = {"K" : "0", "R" : "0"};
			return alphabet;
			
		if (alphabet_name=="aliphatic_helix_uncharged"):
			# 0 basic [RK}
			# 1 hydrophobic uncharged [ILVMFYWCTAG]
			# 2 others
			alphabet = {"A" : "1", "C" : "1",  "F" : "1",
				"G" : "1", "I" : "1", "L" : "1",
				"M" : "1",
				"T" : "1", "V" : "1", "W" : "1", "Y" : "1"};
			return alphabet;
		
		if (alphabet_name=="aliphatic_helix_other"):
			alphabet = {"D" : "2", "E" : "2", "H" : "2", "N" : "2",
				"P" : "2", "Q" : "2", 
				"S" : "2"};
			return alphabet;
			
		raise "FeatureCreator: get_alphabet: unknown alphabet name";

#####################################
######### feature functions #########
#####################################

	def getFeatures(self, i):
		return self.__feature_vectors[i];

	def addNames(self):
		for sequence_nr in range(self.__aa_sequences.size()):
			name = self.__aa_sequences.get(sequence_nr)[0][:6];
			self.__feature_vectors[sequence_nr].append(name);
		self.__labels.append("name");

	def addResponse(self, feature):
		if self.__response == None:
			self.__response = [];
			for sequence_nr in range(self.__aa_sequences.size()):
				self.__response.append(feature);
		else:
			for sequence_nr in range(self.__aa_sequences.size()):
				self.__response[sequence_nr] = feature;
		self.__response_classes = uniqify(self.__response);
		
	def setResponseClasses(self, response_vector):
		self.__response_classes = response_vector;

	def __addDictToFeatures(self,sequence_nr, hist, alphabet_values = []):
		if len(alphabet_values) == 0:
			for key in hist.keys():
				self.__feature_vectors[sequence_nr].append(hist[key]);
		else:
			for key in alphabet_values:
				self.__feature_vectors[sequence_nr].append(hist[key]);

	def __translateSequence(self,sequence, alphabet):
		new_sequence = "";
		for elem in sequence:
			if elem in alphabet.keys():
				new_sequence = new_sequence + alphabet[elem];
			else:
				new_sequence = new_sequence + elem;
		return new_sequence;

	#### simple 
	
	def addSize(self):
		for sequence_nr in range(self.__aa_sequences.size()):		
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			self.__feature_vectors[sequence_nr].append(len(sequence));
		self.__labels.append("size");

	#### count letters

	def __countLetter(self, translated_sequence, letter):
		count = 0;
		for k in range(len(translated_sequence)-len(letter)+1):
			if (translated_sequence[k:(k+len(letter))] == letter):
				count += 1;
		return count;

	def addCountLetter(self, start_pos, end_pos,letter, alphabet_name):
		translation_alphabet = self.__get_alphabet(alphabet_name); 
		for sequence_nr in range(self.__aa_sequences.size()):	
			sequence = self.__translateSequence(self.__aa_sequences.get(sequence_nr)[1][start_pos:end_pos], translation_alphabet);
			counts = self.__countLetter(sequence, letter);
			self.__feature_vectors[sequence_nr].append(counts);
		self.__labels.append("letterCount|"+str(start_pos)+"-"+str(end_pos)+"|"+letter+"|"+alphabet_name);
		
	def addCountLetterNormed(self, start_pos, end_pos,letter, alphabet_name):
		translation_alphabet = self.__get_alphabet(alphabet_name); 
		for sequence_nr in range(self.__aa_sequences.size()):		
			sequence = self.__translateSequence(self.__aa_sequences.get(sequence_nr)[1], translation_alphabet);
			counts = self.__countLetter(sequence[start_pos:end_pos], letter);
			wholecounts = self.__countLetter(sequence, letter);
			self.__feature_vectors[sequence_nr].append(float(counts) / (wholecounts+1));
		self.__labels.append("letterCountNormed|"+str(start_pos)+"-"+str(end_pos)+"|"+letter+"|"+alphabet_name);
		
		
	def addMaxWindowCountLetter(self, start_pos, end_pos, wide, letter, alphabet_name):
		translation_alphabet = self.__get_alphabet(alphabet_name); 
		for sequence_nr in range(self.__aa_sequences.size()):		
			sequence = self.__translateSequence(self.__aa_sequences.get(sequence_nr)[1][start_pos:end_pos], translation_alphabet);
			max_count = 0;
			for pos in range(1, len(sequence)-(wide*2-1)):
				counts = self.__countLetter(sequence[pos-wide:pos+wide], letter);
				if counts > max_count:
					max_count = counts;
			self.__feature_vectors[sequence_nr].append(max_count);
		self.__labels.append("maxLetterCountWin|"+str(start_pos)+"-"+str(end_pos)+"|"+str(wide)+"|"+letter+"|"+alphabet_name);
	

	#### count alphabet  = hist

	def __getHist(self, translated_sequence, alphabet):
		# create alphabet values
		alphabet_values = uniqify(alphabet.values());
		# create empty histogram dictionary from translation alphabet 
		hist = {}; 
		for key in alphabet_values:
			hist[key] = 0;
		# fill histogram by previous translation of sequences
		for aa in translated_sequence:
			if aa in alphabet_values:
				hist[aa] += 1;
        	return(hist);                

	def addHist(self, start_pos, end_pos, alphabet_name = "aminoacid"):
		# get alphabet
		translation_alphabet = self.__get_alphabet(alphabet_name); 
		alphabet_values = uniqify(translation_alphabet.values());
		# get histogramm for every sequence
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start_pos:end_pos];
			sequence = self.__translateSequence(sequence, translation_alphabet);
			hist = self.__getHist(sequence, translation_alphabet);
			self.__addDictToFeatures(sequence_nr, hist, alphabet_values);
		# add labels to label vector

		for key in alphabet_values:
			self.__labels.append("hist|"+str(start_pos)+"-"+str(end_pos)+"|"+alphabet_name+"|"+key);
			
	def addNormalizedHist(self, start_pos, end_pos, alphabet_name = "aminoacid"):
		# get alphabet
		translation_alphabet = self.__get_alphabet(alphabet_name); 
		alphabet_values = uniqify(translation_alphabet.values());
		# get histogramm for every sequence
		for sequence_nr in range(self.__aa_sequences.size()):
			wholesequence = self.__aa_sequences.get(sequence_nr)[1];
			wholesequence = self.__translateSequence(wholesequence, translation_alphabet);
			sequence = wholesequence[start_pos:end_pos];
			hist = self.__getHist(sequence, translation_alphabet);
			histwhole = self.__getHist(wholesequence, translation_alphabet);
			for key in hist.keys():
				hist[key] = float(hist[key])/(histwhole[key]+1);
			self.__addDictToFeatures(sequence_nr, hist, alphabet_values);
		# add labels to label vector

		for key in alphabet_values:
			self.__labels.append("hist_normed|"+str(start_pos)+"-"+str(end_pos)+"|"+alphabet_name+"|"+key);


	#### count alphabet in window

	def addWindowHist(self, start_pos, end_pos, window_size, alphabet_name = "aminoacid"):
		# get alphabet
		translation_alphabet = self.__get_alphabet(alphabet_name); 
		alphabet_values = uniqify(translation_alphabet.values());
		# get histogramm for every sequence and every window
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			sequence = self.__translateSequence(sequence, translation_alphabet);
			for aa_pos in range(start_pos, end_pos):
				window_left = start_pos - window_size;
				window_right = end_pos + window_size;
				if (window_left < 1 & window_right > len(sequence)):
					raise "FeatureCreator: addWindowHist: Start or end position with window is out of sequence bounds";
				hist = self.__getHist(sequence[window_left:window_right], translation_alphabet);
				self.__addDictToFeatures(sequence_nr, hist, alphabet_values);
		# add labels to label vector
		for aa_pos in range(start_pos, end_pos):
			window_left = start_pos - window_size;
			window_right = end_pos + window_size;
			for key in alphabet_values:
				self.__labels.append("hist|win="+str(window_size)+"|"+str(aa_pos)+"|"+alphabet_name+"|"+key);

	#### count alphabet kmers

	def __getAlphabetKmers(self,alphabet, kmer):
		hist ={}
		alphabet_values = uniqify(alphabet.values());
		# create 1mer histogram
		for elem in alphabet_values:
			hist[elem] = 0;
		# create kmer histogram
		for i in range(1,kmer):
			new_hist = {};
			for k in hist.keys():
				for elem in alphabet_values:
					new_hist[(k+elem)] = 0;
			hist = new_hist;
		# return final empty histogram
		return hist;


	def addHistKmer(self, start_pos, end_pos, kmer, alphabet_name = "aminoacid"):
		# create empty hist
		translation_alphabet = self.__get_alphabet(alphabet_name); 
		hist_blank = self.__getAlphabetKmers(translation_alphabet, kmer);
		# get histogramm for every sequence
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start_pos:end_pos];
			sequence = self.__translateSequence(sequence, translation_alphabet);
			# create hist
			hist = hist_blank.copy();
			for i in range(1,len(sequence)-kmer):
				if sequence[i:(i+kmer)] in hist.keys():
					hist[ sequence[i:(i+kmer)]] = hist[ sequence[i:(i+kmer)]] + 1;
			#hist = self.__getHist(sequence, translation_alphabet);
			# add hist to features
			self.__addDictToFeatures(sequence_nr, hist);
		# add labels to label vector
		for key in hist_blank.keys():
			self.__labels.append("histkmer|"+str(start_pos)+"-"+str(end_pos)+"|"+str(kmer)+"|"+alphabet_name+"|"+key);

	#### add properties
	
	def __getProperties(self, sequence, properties):
		# sum over property
		sum_property = 0;
		for aa in sequence:
			if aa in properties.keys():
				sum_property += properties[aa];
        	return(sum_property);          

	def  addProperties(self, start_pos, end_pos, property_name = "volume" ):
		# get property table
		properties = self.__get_properties(property_name);
		# iterate over sequences
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start_pos:end_pos];
			value = self.__getProperties(sequence, properties);
			self.__feature_vectors[sequence_nr].append(value);	
		self.__labels.append("property_sum|"+str(start_pos)+"-"+str(end_pos)+"|"+property_name);
	
	def addWindowProperties(self, start_pos, end_pos, window_size, property_name = "volume" ):
		# get property table
		properties = self.__get_properties(property_name);
		# get histogramm for every sequence and every window
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			for aa_pos in range(start_pos, end_pos):
				window_left = start_pos - window_size;
				window_right = end_pos + window_size;
				if (window_left < 1 & window_right > len(sequence)):
					raise "FeatureCreator: addWindowProperties: Start or end position with window is out of sequence bounds";
				value = self.__getProperties(sequence[window_left:window_right], properties);
				self.__feature_vectors[sequence_nr].append(value);	
		# add labels to label vector
		for aa_pos in range(start_pos, end_pos):
			window_left = start_pos - window_size;
			window_right = end_pos + window_size;
			self.__labels.append("property_sum|win="+str(window_size)+"|"+str(aa_pos)+"|"+property_name);

	def addMaxProperties(self, start_pos, end_pos, wide, property_name = "volume" ):
		# get property table
		properties = self.__get_properties(property_name);
		# iterate over sequences
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start_pos:end_pos];
			max = -100000;
			for aa_pos in range(1,len(sequence)-wide):
				sum = self.__getProperties(sequence[aa_pos:(aa_pos+wide)], properties);
				if (sum > max):
					max = sum;
			self.__feature_vectors[sequence_nr].append(max);	
		self.__labels.append("property_max|"+str(start_pos)+"-"+str(end_pos)+"|"+str(wide)+"|"+property_name);

	def addMinProperties(self, start_pos, end_pos, wide, property_name = "volume" ):
		# get property table
		properties = self.__get_properties(property_name);
		# iterate over sequences
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start_pos:end_pos];
			min = 100000;
			for aa_pos in range(1,len(sequence)-wide):
				sum = self.__getProperties(sequence[aa_pos:(aa_pos+wide)], properties);
				if (sum < min):
					min = sum;
			self.__feature_vectors[sequence_nr].append(min);	
		self.__labels.append("property_min|"+str(start_pos)+"-"+str(end_pos)+"|"+str(wide)+"|"+property_name);


	#### pseudo aa composition 

	def __getPseudoComposition(self, sequence, dist, letter):
		count = 0;
		if len(sequence) <= dist:
			return count;
		for pos in range(1,len(sequence) - dist):
			if ((sequence[pos] == letter) & (sequence[pos] == sequence[pos+dist])):
				count += 1;
		count = float(count) / (len(sequence)-dist);
		return count;
				
	def addMaxPseudoComposition(self, start_pos, end_pos, alphabet_name):
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start_pos:end_pos];
			#get translation alphabet and translate
			alphabet = self.__get_alphabet(alphabet_name); 
			alphabet_values = uniqify(alphabet.values());
			sequence = self.__translateSequence(sequence, alphabet)
			for letter in alphabet_values:
				max = 0;
				#for dist in range(1, end_pos - start_pos):  !!!!
				for dist in range(1, 20):
					count = -1;
					if (len(sequence)>dist):
						count = self.__getPseudoComposition(sequence, dist, letter);
					if count > max:
						max = count;
				self.__feature_vectors[sequence_nr].append(max);
		# add labels to label vector
		translation_alphabet = self.__get_alphabet(alphabet_name); 
		alphabet_values = uniqify(translation_alphabet.values());
		for key in alphabet_values:
			self.__labels.append("maxPaacount|"+str(start_pos)+"-"+str(end_pos)+"|"+alphabet_name+"|"+key);

	def addMaxPseudoCompositionNormed(self, start_pos, end_pos, alphabet_name):
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start_pos:end_pos];
			#get translation alphabet and translate
			alphabet = self.__get_alphabet(alphabet_name); 
			alphabet_values = uniqify(alphabet.values());
			sequence = self.__translateSequence(sequence, alphabet)
			for letter in alphabet_values:
				max = 0;
				for dist in range(1, 20): 
				#for dist in range(1, end_pos - start_pos): !!!!!!!!!!!!
					count = -1;
					if (len(sequence)>dist):
						count = self.__getPseudoComposition(sequence, dist, letter) / dist;
					if count > max:
						max = count;
				self.__feature_vectors[sequence_nr].append(max);
		# add labels to label vector
		translation_alphabet = self.__get_alphabet(alphabet_name); 
		alphabet_values = uniqify(translation_alphabet.values());
		for key in alphabet_values:
			self.__labels.append("maxPaacountNormed|"+str(start_pos)+"-"+str(end_pos)+"|"+alphabet_name+"|"+key);

	def addMinPseudoComposition(self, start_pos, end_pos, alphabet_name):
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start_pos:end_pos];
			#get translation alphabet and translate
			alphabet = self.__get_alphabet(alphabet_name); 
			alphabet_values = uniqify(alphabet.values());
			sequence = self.__translateSequence(sequence, alphabet)
			for letter in alphabet_values:
				min = 100000;
				#for dist in range(1, end_pos - start_pos): !!!!!!!!!!
				for dist in range(10, 20): 
					count = -1;
					if (len(sequence)>dist):
						count = self.__getPseudoComposition(sequence, dist, letter) / dist;
					if count < min:
						min = count;
				self.__feature_vectors[sequence_nr].append(min);
		# add labels to label vector
		translation_alphabet = self.__get_alphabet(alphabet_name); 
		alphabet_values = uniqify(translation_alphabet.values());
		for key in alphabet_values:
			self.__labels.append("minPaacount|"+str(start_pos)+"-"+str(end_pos)+"|"+alphabet_name+"|"+key);

				
	def addPseudoComposition(self, start_pos, end_pos, dist_start, dist_end, alphabet_name):
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start_pos:end_pos];
			#get translation alphabet and translate
			alphabet = self.__get_alphabet(alphabet_name); 
			alphabet_values = uniqify(alphabet.values());
			sequence = self.__translateSequence(sequence, alphabet)
			for dist in range(dist_start, dist_end):
				for letter in alphabet_values:
					count = self.__getPseudoComposition(sequence, dist, letter);
					self.__feature_vectors[sequence_nr].append(count);
		# add labels to label vector
		translation_alphabet = self.__get_alphabet(alphabet_name); 
		alphabet_values = uniqify(translation_alphabet.values());
		for dist in range(dist_start, dist_end):
			for key in alphabet_values:
				self.__labels.append("Paacount|d="+str(dist)+"|"+str(start_pos)+"-"+str(end_pos)+"|"+alphabet_name+"|"+key);


#### pseudo aa composition 
	
	def __getAutoCorrelation(self, sequence, dist):
		corr = 0;
		for pos in range(1,len(sequence) - dist):
			corr += sequence[pos] * sequence[pos+dist];
		corr = corr / (len(sequence)-dist);
		return corr;

	def __seq2List(self, sequence, properties):
		seq = []
		for aa in sequence:
			if aa in properties.keys():
				seq.append(properties[aa]);
			else:
				seq.append(0.0);
		return seq;

	def addMaxAutoCorrelation(self, start_pos, end_pos, property_name):
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start_pos:end_pos];
			#get translation alphabet and translate
			properties = self.__get_properties(property_name); 
			seq = self.__seq2List(sequence, properties);
			max = 0;
			for dist in range(1, 20): 
			#for dist in range(1, end_pos - start_pos): !!!!!!!!
				count = -1;
				if (len(sequence)>dist):
					count = self.__getAutoCorrelation(seq, dist);
				if count > max:
					max = count;
			self.__feature_vectors[sequence_nr].append(max);
		# add label to label vector
		self.__labels.append("maxAutoCorrelation|"+str(start_pos)+"-"+str(end_pos)+"|"+property_name);

	def addAutoCorrelation(self, start_pos, end_pos, dist_start, dist_end, property_name):
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start_pos:end_pos];
			#get translation alphabet and translate
			properties = self.__get_properties(property_name); 
			seq = self.__seq2List(sequence, properties);
			for dist in range(dist_start, dist_end):
				count = 0;
				# !!!!!!!!
				if (len(sequence)>dist):
					count = self.__getAutoCorrelation(seq, dist);
				self.__feature_vectors[sequence_nr].append(count);
		# add labels to label vector
		for dist in range(dist_start, dist_end):
			self.__labels.append("AutoCorrelatuion|d="+str(dist)+"|"+str(start_pos)+"-"+str(end_pos)+"|"+property_name);

		
	#### GO terms
	
	def addGoTerms(self, select = []):
		go = GoFilter(self.process_id);
		(labels, features, description, most_similar_list) = go.getGoFeatures(self.__aa_sequences, select);
		for sequence_nr in range(self.__aa_sequences.size()):
			self.__feature_vectors[sequence_nr] += features[sequence_nr];
		self.__labels += labels;
		# save found Go terms
		if len(description) > 0:
			offset = len(self.__labels) - len(description[0]);
			for sequence_nr in range(self.__aa_sequences.size()):
				for feature_nr in range(len(description[sequence_nr])):
					self.__description[sequence_nr][feature_nr+offset] = description[sequence_nr][feature_nr];
		return most_similar_list;
	
	#### Prosite Motifs
	
	def addPrositeMotifs(self, cluster = ""):
		prosite = PrositeFilter(self.process_id);
		(labels, features,description) = prosite.getPrositeFeatures(self.__aa_sequences, cluster);
		for sequence_nr in range(self.__aa_sequences.size()):
			self.__feature_vectors[sequence_nr] += features[sequence_nr];
		self.__labels += labels;
		# save found prosite patterns
		if len(description) > 0:
			offset = len(self.__labels) - len(description[0]);
			for sequence_nr in range(self.__aa_sequences.size()):
				for feature_nr in range(len(description[sequence_nr])):
					self.__description[sequence_nr][feature_nr+offset] = description[sequence_nr][feature_nr];
		
	
	#### Sorting Signals
	
	def addSortingSignals(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			features = signals.getAllSignals(sequence);
			self.__feature_vectors[sequence_nr] += features;
		self.__labels += signals.getAllLabels();

	def addCleavageSite(self):
		signals = SortingSignals();
                for sequence_nr in range(self.__aa_sequences.size()):
                        sequence = self.__aa_sequences.get(sequence_nr)[1];
                        feature = signals.getCleavageSite(sequence);
                        self.__feature_vectors[sequence_nr].append(feature);
                self.__labels.append("Cleavage_site");

	def addNLSDBSortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getNLSDBSignals(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("NLSDBB_signal");

	def addMonoNLSSortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
                        sequence = self.__aa_sequences.get(sequence_nr)[1];
                        feature = signals.getMonoNLSSignal(sequence);
                        self.__feature_vectors[sequence_nr].append(feature);
                self.__labels.append("Mono_NLS_signal");

	def addMonoNLS1SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal1(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal1");
	
	def addMonoNLS3SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal3(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal3");
	
	def addMonoNLS4SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal4(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal4");
	
	def addMonoNLS6SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal6(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal6");
	
	def addMonoNLS11SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal11(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal11");
	
	def addMonoNLS12SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal12(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal12");
	
	def addMonoNLS15SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal15(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal15");
	
	def addMonoNLS18SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal18(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal18");
	
	def addMonoNLS19SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal19(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal19");
	
	def addMonoNLS21SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal21(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal21");
	
	def addMonoNLS22SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal22(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal22");
	
	def addMonoNLS25SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal25(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal25");
	
	def addMonoNLS26SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal26(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal26");
	
	def addMonoNLS27SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMonoNLSSignal27(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mono_NLS_signal27");
	
	def addBipartiteNLSSortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getBipartiteNLSSignal(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Bipartite_NLS_signal");
	
	def addERSortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getERSignal(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("ER_signal");
	
	def addER4SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getERSignal4(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("ER_signal4");
	
	def addER5SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getERSignal5(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("ER_signal5");
	
	def addPeroxi1SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getPeroxiSignal(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Peroxi_signal1");

        def addPeroxi3SortingSignal(self):
                signals = SortingSignals();
                for sequence_nr in range(self.__aa_sequences.size()):
                        sequence = self.__aa_sequences.get(sequence_nr)[1];
                        feature = signals.getPeroxiSignal3(sequence);
                        self.__feature_vectors[sequence_nr].append(feature);
                self.__labels.append("Peroxi_signal3");
	
	def addPeroxi4SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getPeroxiSignal4(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Peroxi_signal4");
	
	def addPeroxi8SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getPeroxiSignal8(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Peroxi_signal8");
	
	def addPeroxi9SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getPeroxiSignal9(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Peroxi_signal9");
		
 	def addGolgi1SortingSignal(self):
                signals = SortingSignals();
                for sequence_nr in range(self.__aa_sequences.size()):
                        sequence = self.__aa_sequences.get(sequence_nr)[1];
                        feature = signals.getGolgiSignal(sequence);
                        self.__feature_vectors[sequence_nr].append(feature);
                self.__labels.append("Golgi_signal1");

	def addGolgi3SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getGolgiSignal3(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Golgi_signal3");
	
	def addGolgi4SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getGolgiSignal4(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Golgi_signal4");
	
	def addGolgi5SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getGolgiSignal5(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Golgi_signal5");
	
	def addGolgi6SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getGolgiSignal6(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Golgi_signal6");	
	
	def addGolgi7SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getGolgiSignal7(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Golgi_signal7");
	
	def addGolgi9SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getGolgiSignal9(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Golgi_signal9");
	
	def addGolgi11SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getGolgiSignal11(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Golgi_signal11");
	
	def addGolgi12SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getGolgiSignal12(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Golgi_signal12");
	
	def addVacuole6SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal6(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal6");
		
	def addVacuole7SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal7(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal7");
	
	def addVacuole8SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal8(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal8");
		
	def addVacuole9SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal9(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal9");
	
	def addVacuole18SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal18(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal18");
	
	def addVacuole21SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal21(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal21");
	
	def addVacuole26SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal26(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal26");
	
	def addVacuole27SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal27(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal27");
	
	def addVacuole30SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal30(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal30");
	
	def addVacuole31SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal31(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal31");
	
	def addVacuole33SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal33(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal33");
	
	def addVacuole34SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal34(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal34");
		
	def addVacuole35SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal35(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal35");
		
	def addVacuole36SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal36(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal36");
	
	def addVacuole37SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getVacuoleSignal37(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Vacuole_signal37");
	
	def addLysosomal2SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getLysosomalSignal2(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Lysosomal_signal2");
		
	def addLysosomal4SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getLysosomalSignal4(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Lysosomal_signal4");
	
	def addLysosomal6SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getLysosomalSignal6(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Lysosomal_signal6");
		
	def addTransmembrane2SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getTransMembraneSignal2(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Transmembrane_signal2");
	
	def addTransmembrane3SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getTransMembraneSignal3(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Transmembrane_signal3");
	
	def addGlycositeSortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getGlycoSiteSignal(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Glycosite_signal");
	
	def addNESSortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getNESSignal(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("NES_signal");
	

	def addMIP2SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMIPSignal2(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mip_signal2");
	
	def addMIP3SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMIPSignal3(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mip_signal3");
	
	def addMIP5SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMIPSignal5(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mip_signal5");

	def addMIP9SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMIPSignal9(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mip_signal9");
		
	def addMIP10SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMIPSignal10(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mip_signal10");	
		
	def addMIP12SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMIPSignal12(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mip_signal12");
		
	def addMIP13SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getMIPSignal13(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Mip_signal13");
		
	def addChloro2SortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getChlorSignal2(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("Chloro_signal2");
		
	def addAliphaticHelixNormedSortingSignal(self):
		signals = SortingSignals();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			feature = signals.getAliphaticHelixNormed(sequence);
			self.__feature_vectors[sequence_nr].append(feature);
		self.__labels.append("aliphatic_helix_normed");
	
	
	#### Motifs
	
	def addMotifs(self):
		motifs = MotifFilter();
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			features = motifs.getAllMotifs(sequence);
			self.__feature_vectors[sequence_nr] += features;
		self.__labels += motifs.getAllLabels();
		
	#### add Isoeletric Point
	def addIsoelectricPoint(self):
		emboss = {"C" : 8.5, "D" : 3.9, "E" : 4.1, "H" : 6.5, 
				"K" : 10.8, "R" : 12.5, "Y" : 10.1};
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
		
			nr_c = len(re.findall("C",sequence));
			nr_d = len(re.findall("D",sequence));
			nr_e = len(re.findall("E",sequence));
			nr_h = len(re.findall("H",sequence));
			nr_k = len(re.findall("K",sequence));
			nr_y = len(re.findall("Y",sequence));
			nr_r = len(re.findall("R",sequence));
			
			ph = 0.0;
			netto_charge = 1;
			
			while netto_charge > 0:
				qn1 = -1/ (1+10**(3.65-ph));
				qn2 = -nr_d / (1+10**(emboss["D"] - ph));
				qn3 = -nr_e / (1+10**(emboss['E'] - ph));
				qn4 = -nr_c / (1+10**(emboss['C'] - ph));
				qn5 = -nr_y / (1+10**(emboss['Y'] - ph));
				qp1 = 1 / (1+10**(ph-8.2));
				qp2 = nr_k / (1+10**(ph-emboss['K']));
				qp3 = nr_h / (1+10**(ph-emboss['H']));
				qp4 = nr_r / (1+10**(ph-emboss['R']));
		
				netto_charge = qn1+qn2+qn3+qn4+qn5+qp1+qp2+qp3+qp4;
				
				ph += 0.01
			
			self.__feature_vectors[sequence_nr].append(ph);
			
		self.__labels.append("isoelectric_point");
	
	
	def addExtinctionCoefficient(self):
		coeff = {"W" : 5690, "C" : 120, "Y" : 1280};
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
		
			nr_c = len(re.findall("C",sequence));
			nr_w = len(re.findall("W",sequence));
			nr_y = len(re.findall("Y",sequence));
			
			value = nr_y * coeff['Y'] + nr_c * coeff['C']+nr_w*coeff['W'];
			
			self.__feature_vectors[sequence_nr].append(value);
			
		self.__labels.append("extinction_coefficient");	
		
	def addMass(self):
		mass = {"A" : 71.0788, "C" : 103.1388, "D" : 115.0866, "E" : 129.1155, "F" : 147.1766,
				"G" : 57.0519, "H" : 137.1411, "I" : 113.1594, "K" : 128.1741, "L" : 113.1594,
				"M" : 131.1926, "N" : 114.1038, "P" : 97.1167, "Q" : 128.1307, "R" : 156.18575,
				"S" : 87.0782, "T" : 101.1051, "V" : 99.1326, "W" : 186.07931, "Y" : 163.1760};
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			score = 0;
			for aa in sequence:
				if aa in mass.keys():
					score += mass[aa];
			self.__feature_vectors[sequence_nr].append(score);
			
		self.__labels.append("peptide_mass");	
		
	def addSP(self):
		charge = self.__get_properties("charge");
		hydr = self.__get_properties("hydrophob");
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			
			c = self.__getProperties(sequence[1:4], charge);
			h = self.__getProperties(sequence[5:15], hydr);
			cleav = len(re.findall("[GAS].[GAS]",sequence[15:30]));
			value = (c+h)*cleav;
			self.__feature_vectors[sequence_nr].append(value);
			
		self.__labels.append("SP");	
	
	def addSP2(self):
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			
			value = len(re.findall("^M.{0,4}[KR].{0,10}LL.{0,10}[GAS].[GAS]",sequence));
			self.__feature_vectors[sequence_nr].append(value);
			
		self.__labels.append("SP");	
	
	def addSolventaccesibility(self):
		#http://curser.science.ru.nl/content-e/pub_NWI/Bioinformatics%20Summerschool%202006%20light%2003092006/documents/do_4623.htm
		sea30 = {"S" : 0.7, "T" : 0.71, "A" : 0.48, "G" : 0.51, "P" : 0.78, "C" : 0.32,
				"D" : 0.81, "E" : 0.93, "Q" : 0.81, "N" : 0.82, "L" : 0.41, "I" : 0.39,
				"V" : 0.40, "M" : 0.44, "F" : 0.42, "Y" : 0.67, "W" : 0.49, "K" : 0.93,
				"R" : 0.84,  "H" : 0.66};
		sea20 = {"S" : 0.1, "T" : 0.13, "A" : 0.17, "G" : 0.13, "P" : 0.09, "C" : 0.14,
				"D" : 0.10, "E" : 0.03, "Q" : 0.09, "N" : 0.08, "L" : 0.10, "I" : 0.14,
				"V" : 0.10, "M" : 0.36, "F" : 0.16, "Y" : 0.13, "W" : 0.07, "K" : 0.05,
				"R" : 0.11,  "H" : 0.15};
		sea10 = {"S" : 0.2, "T" : 0.16, "A" : 0.35, "G" : 0.36, "P" : 0.13, "C" : 0.54,
				"D" : 0.09, "E" : 0.04, "Q" : 0.10, "N" : 0.10, "L" : 0.49, "I" : 0.47,
				"V" : 0.50, "M" : 0.20, "F" : 0.42, "Y" : 0.20, "W" : 0.44, "K" : 0.02,
				"R" : 0.05,  "H" : 0.19};
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			value = 0;
			for aa in sequence:
				if aa in sea30.keys():
					value += sea30[aa]*30+sea20[aa]*20+sea10[aa]*10;
			self.__feature_vectors[sequence_nr].append(value);
			
		self.__labels.append("Solvent_accesibility");
	
	def addSurfaceCharge(self):
		#http://curser.science.ru.nl/content-e/pub_NWI/Bioinformatics%20Summerschool%202006%20light%2003092006/documents/do_4623.htm
		sea30 = {"S" : 0.7, "T" : 0.71, "A" : 0.48, "G" : 0.51, "P" : 0.78, "C" : 0.32,
				"D" : 0.81, "E" : 0.93, "Q" : 0.81, "N" : 0.82, "L" : 0.41, "I" : 0.39,
				"V" : 0.40, "M" : 0.44, "F" : 0.42, "Y" : 0.67, "W" : 0.49, "K" : 0.93,
				"R" : 0.84,  "H" : 0.66};
		sea20 = {"S" : 0.1, "T" : 0.13, "A" : 0.17, "G" : 0.13, "P" : 0.09, "C" : 0.14,
				"D" : 0.10, "E" : 0.03, "Q" : 0.09, "N" : 0.08, "L" : 0.10, "I" : 0.14,
				"V" : 0.10, "M" : 0.36, "F" : 0.16, "Y" : 0.13, "W" : 0.07, "K" : 0.05,
				"R" : 0.11,  "H" : 0.15};
		sea10 = {"S" : 0.2, "T" : 0.16, "A" : 0.35, "G" : 0.36, "P" : 0.13, "C" : 0.54,
				"D" : 0.09, "E" : 0.04, "Q" : 0.10, "N" : 0.10, "L" : 0.49, "I" : 0.47,
				"V" : 0.50, "M" : 0.20, "F" : 0.42, "Y" : 0.20, "W" : 0.44, "K" : 0.02,
				"R" : 0.05,  "H" : 0.19};
		charge = self.__get_properties("charge");
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			value = 0;
			for aa in sequence:
				if aa in sea30.keys():
					value += charge[aa]*sea30[aa]*30+sea20[aa]*20+sea10[aa]*10;
			self.__feature_vectors[sequence_nr].append(value);
			
		self.__labels.append("Surface_charge");
	
	def countAverageRe(self,rex,start,stop):
		count = 0.0;
		min = 10000;
		max = 0;
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start:stop];
			c = len(re.findall(rex,sequence)) ;
			if c < min:
				min = c
			if c > max:
				max = c;
			if c>0:
				count += 1
		count = count / self.__aa_sequences.size();
		print "min"+str(min);
		print "max"+str(max);
		return count;
	
	def test(self,start,stop):
		empty = {'A': 0, 'C':0,'D' :0, 'E' :0, 'F':0, 'G':0,'H':0,'I':0,'K':0,'L':0,
			'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':0};
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start:stop];
			hist = empty.copy();
			for aa in sequence:
				if aa in hist.keys():
					hist[aa] += 1;
			mean = 0;
			for aa in hist.keys():
				mean += hist[aa];
			mean = mean / 20;
			sd = 0;
			for aa in hist.keys():
				sd += (hist[aa]-mean)**2;
			sd = sqrt(sd / 20);
			self.__feature_vectors[sequence_nr].append(sd);
		self.__labels.append("test");
	
			
	def test2(self,start,stop,w):
		signals = SortingSignals();
		#signals.statNLS();
		#empty = {'A': 0, 'C':0,'D' :0, 'E' :0, 'F':0, 'G':2,'H':0,'I':0,'K':10,'L':2,
		#	'M':0,'N':0,'P':2,'Q':2,'R':10,'S':2,'T':0,'V':0,'W':0,'Y':0};
		empty = {'A': 0.03, 'C':0.008,'D' :0.02, 'E' :0.03, 'F':0.09, 'G':0.03,'H':0.009,'I':0.016,'K':0.264,'L':0.05,
			'M':0.02,'N':0.02,'P':0.06,'Q':0.04,'R':0.266,'S':0.03,'T':0.03,'V':0.02,'W':0.004,'Y':0.008};
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start:stop];
			hist = empty;
			max = 0;
			for window in range(len(sequence)-w):
				value = 0;
				for aa in sequence[window:window+w]:
					if aa in hist.keys():
						value += hist[aa];
				if value > max:
					max = value
			self.__feature_vectors[sequence_nr].append(max);
		self.__labels.append("test2|"+str(start)+"-"+str(stop)+"|"+str(w));

	def test3(self,start,stop,w1,w2):
		signals = SortingSignals();
		helix = {"A" : "1",  "F" : "1","G" : "1", "H" : "1", "K" : "1", "L" : "1",
				"M" : "1", "N" : "1", "R" : "1"};
		empty = {'A': 0.03, 'C':0.008,'D' :0.02, 'E' :0.03, 'F':0.09, 'G':0.03,'H':0.009,'I':0.016,'K':0.264,'L':0.05,
			'M':0.02,'N':0.02,'P':0.06,'Q':0.04,'R':0.266,'S':0.03,'T':0.03,'V':0.02,'W':0.004,'Y':0.008};
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1][start:stop];
			hist = empty;
			max = 0;
			for window in range(len(sequence)-(w1-w2)):
				value1 = w1;
				for aa in sequence[window:window+w1]:
					if aa in helix.keys():
						value1 -= 1;
				value2 = 0;
				for aa in sequence[window+w1:window+w1+w2]:
					if aa in hist.keys():
						value2 += hist[aa];
						
				value = value1 + value2*w1*8;
				if value > max:
					max = value
			self.__feature_vectors[sequence_nr].append(max);
		self.__labels.append("test3|"+str(start)+"-"+str(stop)+"|"+str(w1)+"|"+str(w2));

	def test4(self,w):
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			max = 0;
			for window in range(len(sequence)-w):
				value = 0;
				acid = 0;
				for aa in sequence[window:window+w]:
					if aa == 'L':
						value += 1;
					if aa == 'D' or aa == 'E':
						acid += 1;
				value += acid;
				if value > max:
					max = value
			self.__feature_vectors[sequence_nr].append(max);
		self.__labels.append("test4|"+str(w));

	
	def evalNormDistribution(self, mean, sd, value):
		#print("m="+str(mean)+" var="+str(sd)+" value="+str(value));
		#print sd;
		if sd != 0.0:
			#print "h";
			res = 1/(sd * sqrt(2*3.14)) * exp( -1* (value-mean)**2 / (2*sd**2));
		else:
			#print "KKKKKKKKK";
			res = 1;
		#print(res);
		return res;
	
	def classify(self):
		response = uniqify(self.__response);
		zero_vector = [];
		for elem in response:
			zero_vector.append(0.0);
		prior = zero_vector;
		
		count_matrix = [];
		for i in range(len(response)):
			count_matrix.append([]);
			for k in range(len(self.__feature_vectors[0])):
				count_matrix[i].append([]);
		
		# list of lists
		# each list for a attribute with values of each class
		for i in range(len(self.__feature_vectors[0])):
			for sequence_nr in range(self.__aa_sequences.size()):
				r = response.index(self.__response[sequence_nr]);
				count_matrix[r][i].append(self.__feature_vectors[sequence_nr][i]);
		for sequence_nr in range(self.__aa_sequences.size()):
			r = response.index(self.__response[sequence_nr]);
			prior[r] += 1.0;
		for i in range(len(response)):
			prior[i] = prior[i] / self.__aa_sequences.size();
		
		# calculate distributions
		mean_matrix = []
		sd_matrix = [];
		for i in range(len(response)):
			mean_matrix.append([]);
			sd_matrix.append([]);
			for k in range(len(self.__feature_vectors[0])):
				mean_matrix[i].append(0.0);
				sd_matrix[i].append(0.0);
				for elem in count_matrix[i][k]:
					mean_matrix[i][k] += float(elem);
				mean_matrix[i][k] = mean_matrix[i][k] / len(count_matrix[i][k]);
				for elem in count_matrix[i][k]:
					sd_matrix[i][k] += (float(elem) - mean_matrix[i][k])**2;
				sd_matrix[i][k]  = sqrt(sd_matrix[i][k]  / len(count_matrix[i][k]));
		
		#print(mean_matrix);
		
		# classify
		error = 0;
		#ig = [0.215, 0.186,0.423,0.127,0.582,0.585,0.399,0.655,0.181,0.405,0.644,0.426,0.401,0.157,0.228];
		for sequence_nr in range(self.__aa_sequences.size()):
			max = 0.0;
			max_class = 0;
			#print("sequence "+str(sequence_nr));
			probs = []
			sum = 0;
			for i in range(len(response)):
				prob = prior[i];
				#prob = 1.0;
				for k in range(len(self.__feature_vectors[sequence_nr])):
					prob *= self.evalNormDistribution(mean_matrix[i][k],sd_matrix[i][k],float(self.__feature_vectors[sequence_nr][k]))  ;
					#print("p="+str(prob));
						
				#print "final"
				#print(prob);
				#print("class "+str(i)+" with prob=" +str(prob));
				probs.append(prob);
				sum += prob;
			
			#print(probs);
			#print(sum);
			
			for k in range(len(response)):
				if sum != 0:
					probs[k] = probs[k] / sum;
				else:
					probs[k] = 0;
				if probs[k] > max:
					max_class = k;
					max = probs[k];
			#print(probs);
			#print("real="+self.__response[sequence_nr]+" predicted="+response[max_class]);
			if response[max_class] != self.__response[sequence_nr]:
				error +=1;
		print(1.0 - float(error) / self.__aa_sequences.size());
		
			
	
	def delme(self):
		count = 0;
		comp = [ "SP", "OT",  "CTP", "OT", "MTP"];
		#prior = [0.0105, 0.1415, 0.1405, 0.026, 0.2078,  0.025, 0.075,
		#0.2368, 0.03323, 0.0856];
		u = [0,0,0,0,0];
		unsure = 0;
		print( Set(self.__response));
		for sequence_nr in range(self.__aa_sequences.size()):
			sequence = self.__aa_sequences.get(sequence_nr)[1];
			prob = [];
			prob.append(1-float(self.__feature_vectors[sequence_nr][1]));
			prob.append(float(self.__feature_vectors[sequence_nr][2]));
			prob.append(1-float(self.__feature_vectors[sequence_nr][3]));
			prob.append(1-float(self.__feature_vectors[sequence_nr][4]));
			prob.append(1-float(self.__feature_vectors[sequence_nr][5]));
			print(prob);
						
			prob[4] *= 1-float(self.__feature_vectors[sequence_nr][0]);
			prob[2] *= float(self.__feature_vectors[sequence_nr][0]);
			
			
			print(str(sequence_nr)+" from "+self.__response[sequence_nr]+ " with probs: ");
			print(prob);
			#for i in range(len(prob)):
			#	prob[i] *= (prior[i]**0.1);
			max_i = 0;
			max = 0;
			second = 0;
			
			for i in range(len(prob)):
				if prob[i] > max:
					second = max;
					max = prob[i];
					max_i = i;
			print("-->"+ comp[max_i]);
			if (max-0.1 < second or max < 0.3):
				unsure += 1;
			if self.__response[sequence_nr] == comp[max_i]:
				#print(str(sequence_nr)+" was right");
				count +=1;
		print(float(count) / float(self.__aa_sequences.size()) );
		print unsure;
		print(float(count) / float(self.__aa_sequences.size() - unsure) );
		
				
	
	#### Sub Predictor
	
	def addSubPredictor(self, model_name):
		predictor = Predictor(self.__aa_sequences);
		#print("use sub predictor");
		result = predictor.predict(model_name, sort=False);
		for sequence_nr in range(self.__aa_sequences.size()):
			if len(result[0]) > 2:
				class_prob = [];
				for elem in result[sequence_nr]:
					class_prob.append(elem[1]);
				self.__feature_vectors[sequence_nr] += class_prob;
			else:
				self.__feature_vectors[sequence_nr].append(result[sequence_nr][0][1]);
		if len(result[0]) > 2:
			new_labels = [];
			for elem_nr in range(len(result[0])):
				new_labels.append("Subpredictor:"+model_name+"_Prob_"+result[0][elem_nr][0].strip(" "));
			self.__labels += new_labels;
		else:
			self.__labels.append("Subpredictor:"+model_name+"_Prob_"+result[0][0][0].strip(" "));

		
#####################################
##### read from file functions ######
#####################################

        def readARFFFeatureFile(self, file_name):
		ff_file = open(file_name,"r");			
		line = ff_file.readline();
		data=False
		while line:
			line2 = line.upper();
			if(line2[0]=='@'):
				if(line2[0:10]=="@ATTRIBUTE"):
					if line.find("{") == -1:
						if '\'' in line:
							self.__labels.append('\''+(line[11:-1].split('\''))[1]+'\'');
						else:
							self.__labels.append((line[11:-1].split(' '))[0]);
					else:
						self.__response = [];
				if (line2[0:5]=="@DATA"):
					data=True;
			else:
				line = line.strip(' \n');
				if ((len(line)>2) & (data==True)):
					features = line.split(',');
					self.__aa_sequences.append(("A","A"));
					if self.__response == None:
						self.__feature_vectors.append(features);
					else:
						self.__feature_vectors.append(features[0:-1]);
						self.__response.append(features[-1]);
			line = ff_file.readline();
		self.__response_classes = uniqify(self.__response);
		ff_file.close();


#####################################
##### write to file functions #######
#####################################

	def writeARFFFeatureFile(self, file_name,shuffle=False):
	
		ff_file = open(file_name,"w");	
			
		ff_file.write("@RELATION "+str(file_name)+ "\n");
		for label_nr in range(0,len(self.__labels)):
			ff_file.write("@ATTRIBUTE " + str(self.__labels[label_nr]) + " numeric \n");
	
		if self.__response_classes == None:
			if self.__response != None:
				#get all class labels
				classes = uniqify(self.__response);
				ff_file.write("@ATTRIBUTE response {");
				for nr in range(0,len(classes)-1):
					ff_file.write(str(classes[nr])+",");
				ff_file.write(str(classes[-1])+"} \n \n");
		else:
			classes = self.__response_classes;
			ff_file.write("@ATTRIBUTE response {");
			for nr in range(0,len(classes)-1):
				ff_file.write(str(classes[nr])+",");
			ff_file.write(str(classes[-1])+"} \n \n");
			
		ff_file.write("@DATA \n");
		
		if (shuffle):
			random.shuffle(self.__feature_vectors);
		
		for sequence_nr in range(self.__aa_sequences.size()):
			for feature in self.__feature_vectors[sequence_nr]:
				ff_file.write(str(feature) + ",");
			if self.__response != None:
				ff_file.write(str(self.__response[sequence_nr]));
			ff_file.write(" \n");
		ff_file.close();


	def writeScaleFeatureFile(self, file_name):
		ff_file = open(file_name,"w");					
		
		for sequence_nr in range(self.__aa_sequences.size()):
			ff_file.write(str(self.__response[sequence_nr]) + " ");
			for feature_nr in range(0,len(self.__feature_vectors[sequence_nr])):
				feature = self.__feature_vectors[sequence_nr][feature_nr];
				label = self.__labels[feature_nr];
				ff_file.write(str(feature_nr)+":"+str(feature) + " ");
			ff_file.write("\n");
		ff_file.close();
		

	def writePlainFeatureFile(self, file_name, shuffle=False, labels=True):
		ff_file = open(file_name,"w");
		
		if (labels):
			#write label string
			label_string = "";
			for label in self.__labels:
				label_string += label + "\t";
			label_string += "\n";
			ff_file.write(label_string);
			
		if (shuffle):
			random.shuffle(self.__feature_vectors);

		self.__writeFeaturesToFile(ff_file);


	def __writeFeaturesToFile(self, ff_file):
		for sequence_nr in range(self.__aa_sequences.size()):
			out_string = "";
			for feature in self.__feature_vectors[sequence_nr]:
				out_string += str(feature) + "\t";
			if self.__response != None:
				out_string += self.__response[sequence_nr];
			out_string += "\n";
			ff_file.write(out_string);
		ff_file.close();


