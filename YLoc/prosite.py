from math import log;
import os, sys, re;
from aasequences import *;
import copy;
import config;

class PrositeFilter(object):
	PATH = config.PATH_PROSITE;
	PATH_PS = config.PATH_PSSCAN; 
	PATH_TMP = config.PATH_TMP; 
	
	"Class creates features from sequences by extracting prosite domains"
	
	#process_id = "1";
	
	def __init__(self, id = "4"):
		self.process_id = id;
		
	def __getPrositeMatches(self, prosite_ids, sequences,cluster=""):
		# include psscan path
		current_path = os.environ['PATH'];
		if current_path.find(self.PATH_PS[:-1]) == -1:
			os.environ['PATH'] = self.PATH_PS[:-1] + ":" + current_path
		
		#replace names in sequences
		for i in range(sequences.size()):
			(name,seq) = sequences.get(i);
			name = "Sq%05d" % i;
			sequences.set( (name,seq), i);
		
		sequences.write_fasta_file(self.PATH_TMP+"tmp_prosite_"+self.process_id+".fasta");
		prosite_fname = self.PATH + "prosite.dat";
		if len(cluster) > 0:
			prosite_fname += "_"+cluster;
		#print(prosite_fname);
		os.system("perl %sps_scan.pl -d %s %stmp_prosite_%s.fasta > %s" % (self.PATH_PS,prosite_fname, self.PATH_TMP,self.process_id, self.PATH_TMP+"prositeout_"+self.process_id));
		# create prosite term list	
		prosite_list = [];
		# current position and entry name in list
		current = -1;
		current_name = "";
		# read-in output of blastp
		f = open(self.PATH_TMP+"prositeout_"+self.process_id,"r");
		line = f.readline();
		while line:
			if line[0]==">":
				new_name = line[1:8];
				#print line;
				#print new_name;
				if current < 0:
					current += 1;
					prosite_list.append([]);

				while sequences.get(current)[0][:7] != new_name:
					current += 1;
					pos = current+1;
					prosite_list.append([]);
					
				current_name = new_name;
				prosite_id = re.findall("PS[0123456789]{5} ",line)[0].strip(" \t\n");
				#print prosite_id;
				if prosite_id in prosite_ids:
					prosite_list[current].append(prosite_id);
					#print("add to "+str(current));
			line = f.readline();
		f.close();
		
		# in case one of the last is missing...fill up prosite list
		#print("fill");
		for i in range(sequences.size()-current):
			#print(i);
			prosite_list.append([]);
		
		# delete files
		os.system("rm %s" % (self.PATH_TMP+"tmp_prosite_"+self.process_id+".fasta"));
		os.system("rm %s" % (self.PATH_TMP+"prositeout_"+self.process_id));
		
		return prosite_list;
	
	def __uniqify(self, seq): 
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
	
	def list2String(self, list):
		s = "";
		for elem in list:
			s += "<a href=http://www.expasy.ch/prosite/"+str(elem)+">"+str(elem)+"</a> ";
		return s;
	
	def __getPrositeList(self, file_name):
		prosite_ids = [];
		f = open(file_name,"r");
		line = f.readline();
		while line:
			line = line.strip(" \n\r");
			prosite_ids.append(line);
			line = f.readline();		
		f.close();
		return prosite_ids;	
	
	def getPrositeFeatures(self,sequences,cluster= []):
		# create list with all prosite terms
		prosite_ids = [];
		file_name = self.PATH+"prosite_ids.txt";
		cluster_dict = {};
		one_cluster_only = [];
		if len(cluster) > 0:
			for cluster_name in cluster:
				if "_" in cluster_name:
					file_name = self.PATH+"prosite_ids.txt_"+cluster_name;
					new_ids = self.__getPrositeList(file_name);
					for elem in new_ids:
						if elem not in cluster_dict.keys():
							cluster_dict[elem] = [];
						help = cluster_dict[elem];
						help.append(cluster_name);
						cluster_dict[elem] = help;
					prosite_ids += new_ids;
				else:
					prosite_ids.append(cluster_name);
			if len(cluster) ==1:
				if "_" in cluster[0]:
					one_cluster_only = cluster[0];
		else:
			prosite_ids = self.__getPrositeList(file_name);
		
		#print(sequences.size());
		#get prosite ids for sequences
		(name,seq) = sequences.get(0)
		prosite_list = self.__getPrositeMatches(prosite_ids, copy.deepcopy(sequences),one_cluster_only);
		(name,seq) = sequences.get(0)
		#print(prosite_list);
		#print(len(prosite_list));
		#print(sequences.size());
		
		if len(cluster) == 0:
			#create features
			empty = [];
			for elem in prosite_ids:
				empty.append(0);
			feature_list = [];
			for i in range(sequences.size()):
				feature_list.append(copy.copy(empty));
				for elem in prosite_list[i]:
					k = prosite_ids.index(elem);
					feature_list[-1][k] += 1;

			return (prosite_ids,feature_list,[]);
		else:
			feature_list = [];
			description_list = [];
			empty = [];
			description_empty = [];
			for elem in cluster:
				empty.append(0);
				description_empty.append("");
			for i in range(sequences.size()):
				feature_list.append(copy.copy(empty));
				description_list.append(copy.copy(description_empty));
				for elem in prosite_list[i]:
					if elem in cluster:
						index = cluster.index(elem);
						feature_list[-1][index] = 1;
						description_list[-1][index] = "<a href=http://www.expasy.org/prosite/"+str(elem)+">"+str(elem)+"</a>";
					else:
						#print elem;
						if elem in cluster_dict.keys():
							found_cluster = cluster_dict[elem];
							for cl in found_cluster:
								index = cluster.index(cl);
								feature_list[-1][index] = 1;
								if len(description_list[-1][index])== 0:
									description_list[-1][index] += "<a href=http://www.expasy.org/prosite/"+str(elem)+">"+str(elem)+"</a>";
								elif not elem in description_list[-1][index]:
									description_list[-1][index] += ", <a href=http://www.expasy.org/prosite/"+str(elem)+">"+str(elem)+"</a>";
								
			
			#print feature_list;
			#print description_list;
			#x = 8/0;
			for i in range(len(cluster)):
				if "_" in cluster[i]:
					cluster[i] = "prosite_"+cluster[i];
			prosite_ids = cluster;
			return (prosite_ids,feature_list,description_list);
