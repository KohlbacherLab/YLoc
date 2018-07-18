import os, sys, re;
import copy;
from math import *;
class Prediction(object):

	def __init__(self):
		self.locations = [];
		self.putative_locations = [];
		self.probability_distribution = [];
		self.attr_nr = 0;
		self.probability_quantity = [];
		self.probability_quality = [];
		self.discrimination_scores = [];
		self.discretization_intervals = [];
		self.matched_intervals = [];
		self.merged_true_classes = [];
		self.true_class_attribute_distribution = [];
		self.pos_classes = [];
		self.pos_class_attribute_distribution = [];
		self.neg_classes = [];
		self.neg_class_attribute_distribution = [];
		
		self.class_names = [];
		self.feature_names = [];
		
		self.feature_description = [];
		
		self.id = 0;
		self.date = "";
		self.ip = "";
		self.model = "";
		self.sequence = "";
		self.sequence_name = "";
		self.sequence_id = 0;
		self.null_prob = 0.0;
		self.most_similar = [];
	
	
	
	def getObjectAsList(self):
		o_list = [];
		o_list.append(self.locations);
		o_list.append(self.putative_locations);
		o_list.append(self.probability_distribution);
		o_list.append(self.attr_nr);
		o_list.append(self.probability_quantity);
		o_list.append(self.probability_quality);
		o_list.append(self.discrimination_scores);
		o_list.append(self.discretization_intervals);
		o_list.append(self.matched_intervals);
		o_list.append(self.merged_true_classes);
		o_list.append(self.true_class_attribute_distribution);
		o_list.append(self.pos_classes);
		o_list.append(self.pos_class_attribute_distribution);
		o_list.append(self.neg_classes);
		o_list.append(self.neg_class_attribute_distribution);
		o_list.append(self.class_names);
		o_list.append(self.feature_names);
		o_list.append(self.feature_description);
		o_list.append(self.id);
		o_list.append(self.date);
		o_list.append(self.model);
		o_list.append(self.sequence);
		o_list.append(self.sequence_name);
		o_list.append(self.sequence_id);
		o_list.append(self.null_prob);
		o_list.append(self.most_similar);
		return o_list;		

	def fillObjEctFromStringOfList(self, s):
		o_list = eval(s);
		self.locations = o_list[0];
		self.putative_locations = o_list[1];
		self.probability_distribution = o_list[2];
		self.attr_nr = o_list[3];
		self.probability_quantity = o_list[4];
		self.probability_quality = o_list[5];
		self.discrimination_scores = o_list[6];
		self.discretization_intervals = o_list[7];
		self.matched_intervals = o_list[8];
		self.merged_true_classes = o_list[9];
		self.true_class_attribute_distribution = o_list[10];
		self.pos_classes = o_list[11];
		self.pos_class_attribute_distribution = o_list[12];
		self.neg_classes = o_list[13];
		self.neg_class_attribute_distribution = o_list[14];
		self.class_names = o_list[15];
		self.feature_names = o_list[16];
		self.feature_description = o_list[17];
		self.id = o_list[18];
		self.date = o_list[19];
		self.model = o_list[20];
		self.sequence = o_list[21];
		self.sequence_name = o_list[22];
		self.sequence_id = o_list[23];
		self.null_prob = o_list[24];
		self.most_similar = o_list[25];

	def getSortedAttributeList(self):
		help = copy.copy(self.discrimination_scores);
		help2 = copy.copy(self.discrimination_scores);
		for i in range(len(help)):
			help[i] = abs(help[i]);
			help2[i] = abs(help2[i]);
		help.sort();
		help = help[::-1];
		indices = [];
		for elem in help:
			indices.append(help2.index(elem));
		return indices;
		
	
	def getSortedLocationList(self):
		help = copy.copy(self.probability_distribution);
		help.sort();
		help = help[::-1];
		indices = [];
		for elem in help:
			indices.append(self.probability_distribution.index(elem));
		return indices;
		
	def getClassNameString(self, list, capitalize = False):
		s = "";
		if len(list) > 0:
			if capitalize:
				s += self.class_names[list[0]][:1].upper()+self.class_names[list[0]][1:];
			else:
				s += str(self.class_names[list[0]]);
			index = 1;
			while index < len(list)-1:
				s += ", " + str(self.class_names[list[index]]);
				index += 1;
			if len(list) > 2:
				s += ", and "+str(self.class_names[list[-1]]);
			elif len(list) == 2:
				s += " and "+str(self.class_names[list[-1]]);
		return s;
		
		
	### help function to extract links
	def __extract_links(self,s):
		pos = s.find("<a");
		pos_end = s[pos:].find(">")+pos;
		pos2 = s.find("</a>");
		i = 1;
		while pos != -1:
			s = s[:pos]+s[pos_end+1:pos2]+s[pos2+4:];
			pos = s.find("<a");
			pos_end = s[pos:].find(">")+pos;
			pos2 = s.find("</a>");
			i+=1;
		return s;	
	
	def getShortReasoning(self,detailed=False):
		
		output_string = "";

		#calc prob
		prob = 0;
		for i in range(len(self.locations)):
			prob += self.probability_distribution[self.locations[i]];

		## calc trust
		qualitative_word = "";
		if self.null_prob >= 0.95:
			qualitative_word = "very strong"
		elif self.null_prob >= 0.8:
			qualitative_word = "strong";
		elif self.null_prob >=0.3:
			qualitative_word = "normal";
		else:
			qualitative_word = "small";
		
		## more than one
		plural="";
		if len(self.locations)>1:
			plural = "s";

		output_string += "<p>";
		output_string +=  self.model+" predicted that protein sequence <i>"+self.sequence_name+"</i> is located in the "+self.getClassNameString(self.locations)+" with a ";
		#print "<a class='helptext' href=#G>";
		output_string +=  "probability of %.1f" % (prob*100) +"%. ";
		#print "<span style='position:absolute;'><i>"+ "bla" +"</i></a> </span>";
		output_string +=  "YLoc has a "+str(qualitative_word)+ " confidence (%.2f" % (self.null_prob) +") that this prediction is reliable.<BR>";
		output_string +=  "</p>";

		ordered_attributes = self.getSortedAttributeList();
		most_important_attr = -1;
		second_important_attr = -1;
		i=0;
		while i<self.attr_nr and second_important_attr==-1:
			if self.discrimination_scores[ordered_attributes[i]]>0 and len(self.pos_classes[ordered_attributes[i]]) > 0:
				if most_important_attr == -1:
					most_important_attr = ordered_attributes[i];	
				else:
					second_important_attr = ordered_attributes[i];	
			i += 1;
		

		output_string +=  "<p class=underl>";
		
		binding_word = "the protein has a";
		if self.feature_description[second_important_attr][-1][-7:] == "protein":
			binding_word = "it is a";		
		
		output_string += "The most important reason for making this prediction is the "+str(self.feature_description[most_important_attr][-2])+" ";
		
		if detailed:
			output_string += str(self.feature_description[most_important_attr][-1])+".";
		else:
			output_string += str(self.__extract_links(self.feature_description[most_important_attr][-1]))+".";

		output_string +=   (" %.0f" % (self.true_class_attribute_distribution[most_important_attr][self.matched_intervals[most_important_attr ]]*100) ) +"% of the proteins from the "+ self.getClassNameString(self.locations) +" have a similar attribute, whereas only about "+( "%.0f" % (self.pos_class_attribute_distribution[most_important_attr][self.matched_intervals[most_important_attr]]*100) )+"% of the proteins form the "+str(self.getClassNameString(self.pos_classes[most_important_attr]))+" show this property. ";
		
		output_string +=  "Moreover, "+str(binding_word)+" "+str(self.feature_description[second_important_attr][-2])+" ";
		if detailed:
			output_string += str(self.feature_description[second_important_attr][-1])+".";
		else:
			output_string += str(self.__extract_links(self.feature_description[second_important_attr][-1]))+".";
		output_string +=   (" %.0f" % (self.true_class_attribute_distribution[second_important_attr][self.matched_intervals[second_important_attr ]]*100) ) +"% of the proteins from the "+ self.getClassNameString(self.locations) +" have a similar attribute, whereas only about "+( "%.0f" % (self.pos_class_attribute_distribution[second_important_attr][self.matched_intervals[second_important_attr]]*100) )+"% of the proteins form the "+str(self.getClassNameString(self.pos_classes[second_important_attr]))+" show this property. There are more properties that support the predicted location"+plural+".";
		output_string +=  "</p>";	

		if detailed:
			output_string +=  "<p class=underl>";
			if len(self.most_similar)>3:
				(id, e_value, identity,go_string) = self.most_similar;
				output_string +=  "The most similar protein from Swiss-Prot 42.0 to the query is <a href='http://www.expasy.org/uniprot/"+str(id)+"' onclick=\"openNewWindow(this.href); return false;\">"+str(id)+"</a> with a local sequence similarity of "+str(identity)+"% (E-value="+str(e_value)+") and is associated with the following GO terms: "+go_string+".";
			else:
				output_string +=  "There exist no very similar protein in Swiss-Prot 42.0.";
			output_string += "</p>";


		return output_string;
