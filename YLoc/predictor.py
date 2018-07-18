from weka import *;
from models import *;
import copy;
from prediction import *;
import time;

class Predictor(object):

	def __init__(self, aa_sequences, id = "1"):
		self.__aa_sequences = aa_sequences;
		self.__loc_names = { "SP": "Secreted", "SR": "Secreted", "Nuc" : "Nucleus    ", "nuc" : "Nucleus    ", "cTP" : "Chloroplast" , "mTP" : "Mitochondria", "Cyt" : "Cytosol    " , "OT" : "Nuclues_or_Cytosol" , "ER" : "Endoplasmatic_Reticulum", "X" : "Other", "cyt" : "Cytosol    ", "gol" : "Golgi     ", "mem" : "Plasma_membrane", "per" : "Peroxisomen", "vac" : "Vacuole    ", "lys" : "Lysosomen"};
		self.process_id = id;

	def __createFeaturesForModel(self, model_name):
		m = Models(self.process_id);
		(features, most_similar) = m.createFeaturesForSequences(copy.copy(self.__aa_sequences),model_name);
		model_file = m.getModelFilePath(model_name);
		return( (model_file, features, most_similar) );

	def __getDescription(self, class_name):
		return self.__loc_names[class_name];
	
	def __compare_classes(self,list1, list2):
		if list1[1] < list2[1]:
			return 1;
		elif list1[1] > list2[1]:
			return -1;
		else:
			return 0;
	
	def __compare_attributes(self,tuple1, tuple2):
		(name1, prob1) = tuple1;
		(name2, prob2) = tuple2;
		if prob1 < prob2:
			return 1;
		elif prob1 > prob2:
			return -1;
		else:
			return 0;

	def useRightQuantity(self,quantities, nr, prediction, attr_index):
		#print(str(attr_index)+"<BR>");
		possible_list = quantities[nr];
		#print(str(possible_list)+"<BR>");
		all = len(prediction.discretization_intervals[attr_index]);
		#print(str(all)+"<BR>");
		this = prediction.matched_intervals[attr_index] + 1;
		#print(str(this)+"<BR>");
		index = int(float(this)/float(all) * (len(possible_list)-1));  # -1 ????????
		#print(str(index)+" done<BR>");
		return possible_list[index];
	
	def fillToDigets(self,nr, digets):
		x = str(nr);
		while len(x) < digets:
			x = "0"+x;
		return x;
	
	def predict(self, model_name, sort = True):
		# get features for this model and given sequences
		(model_file, features, most_similar_list) = self.__createFeaturesForModel(model_name);
		
		# load weka interface and predict
		weka = WekaInterface(self.process_id);
		prediction_list = weka.predictModel( model_file, features);
		
		# fill prediction object with class names and feature names
		m = Models(self.process_id);
		classes = m.getModelClasses(model_name);
		feature_names = m.getModelFeatureNames(model_name);
		for elem in prediction_list:
			elem.class_names = classes;
			elem.feature_names = feature_names;
		
			
		# fill prediction object with feature descriptions
		quantities = m.getQuantityDescriptors(model_name);
		for nr in range(len(prediction_list)):
			elem = prediction_list[nr];
			if len(most_similar_list)>0:
				elem.most_similar = most_similar_list[nr];
			descriptions = m.getModelFeatureDescription(model_name, features, nr);
			elem.feature_description = [];
			# calc for every attribute in every prediction a sentence
			#print(elem.attr_nr);
			for i in range(elem.attr_nr):
				sentences = [];
				sentence_list = descriptions[i];
				#print(sentence_list);
				# create all three type of sentences
				for sentence in sentence_list:
					# skip if no sentence is known
					if len(sentence) > 0:
						s = "";
						if sentence[0] != 0:
							sentences.append( self.useRightQuantity(quantities, sentence[0], elem,i) );
						else:
							sentences.append("");
						sentences.append(sentence[1]);
				#print(sentences);
				elem.feature_description.append(sentences);
		
		# add other information
		current_time = time.localtime();
		time_string = str(current_time[0])+self.fillToDigets(current_time[1],2)+self.fillToDigets(current_time[2],2)+self.fillToDigets(current_time[3],2)+self.fillToDigets(current_time[4],2)+self.fillToDigets(current_time[5],2);
		for i in range(len(prediction_list)):
			prediction_list[i].id = self.process_id;
			prediction_list[i].date = time_string;
			prediction_list[i].model = m.getModelName(model_name);
			(name, seq) = self.__aa_sequences.get(i);
			prediction_list[i].sequence = seq;
			prediction_list[i].sequence_name = name;
		
		return prediction_list;





