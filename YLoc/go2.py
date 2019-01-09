import config;
from math import log;
import os, sys;
from aasequences import *;
import copy;

class GoFilter(object):
	PATH = config.PATH_GO;
	PATH_TMP = config.PATH_TMP;
	CLUSTER_FILE = "multi_goterms_homolog42_3_clustering_0.95.txt";

	"Class creates features from sequences by extracting Go terms of homologeus"

	def __init__(self,id = "1"):
		self.process_id = id;

	def __getGoTerms(self,gomap,sequences):
		sequences.write_fasta_file(self.PATH_TMP+"tmp_go_"+self.process_id+".fasta");
		#os.system("%sblastall -p blastp -d %sswissprot/swissprot -i %stmp_go_%s.fasta >& %s" % (self.PATH_BLAST,self.PATH, self.PATH_TMP, self.process_id,self.PATH_TMP+"blastout_"+self.process_id));

		os.system("blastp -db %sswissprot/swissprot -query %stmp_go_%s.fasta > %s" % (self.PATH, self.PATH_TMP, self.process_id, self.PATH_TMP+"blastout_"+self.process_id));

		# create lis of go terms
		go_list = [];
		most_similar = [];
		# read-in output of blastp
		f = open(self.PATH_TMP+"blastout_"+self.process_id,"r");
		line = f.readline();
		while line:
			if line[:7]=="Query= ":
				go_list.append([]);
				most_similar.append(());
			if line[:19]=="[blastall] WARNING:":
				print "FOUND WARNING";
				line = f.readline();
				go_list.append([]);
				most_similar.append(());
			if line[:42]=="Sequences producing significant alignments":
				line = f.readline();
				line = f.readline();
				stop = False;
				count = 0;
				while (line and not stop):

					if len(line)>0 and line[0] == ">":
						id = line[1:7];
						line = f.readline();
						line = f.readline();
						line = f.readline();
						# get e-value
						line_s = line.split(",");
						e_value = line_s[1][10:].strip(" \s\t\r\n,");
						if e_value[0] == "e":
							e_value = "1"+e_value;
						e_value = float(e_value);
						# get identity
						line = f.readline();
						line_s = line.split("(");
						line_s = line_s[1].split("%");
						identity = int(line_s[0]);

						#add info to most_similar_list
						if count == 0:
							go_string = "";
							if id in gomap.keys():
								help_list = gomap[id];
								help_list = self.__uniqify(help_list);
								for elem in help_list:
									go_string += "<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:"+str(elem)+" >GO:"+str(elem)+"</a> ";
							most_similar[-1] = (id,e_value,identity, go_string);
						count += 1;
						if (not e_value == 0):
							#print(e_value);
							e_value = -1*log(e_value)
							if (e_value < 5000 and e_value >= 10 and identity >= 30):
							#if (e_value > 10):
								#print("take");
								#print(e_value);
								#print(id);
								if id in gomap.keys():
									#print("insert");
									#print(gomap[id]);
									#print((e_value/70.0));
									for elem in gomap[id]:
										go_list[-1].append((elem,1.0));
										#go_list[-1].append((elem,(e_value/70.0)));
										#go_list[-1].append((elem,(e_value)));
							elif (e_value < 10 or identity < 30):
							#e < 0.001 & i < 0
								stop = True;
							#if (count > 5):
							#	stop = True;
					else:
						if len(line)>5 and line[0:6] == "BLASTP":
							stop = True;
					line = f.readline();

				#print go_list[-1];
				#print "LENGTH:  "+str(len(go_list))+ " #####";

			line = f.readline();
		f.close();

		# delete files
		os.system("rm %s" % (self.PATH_TMP+"tmp_go_"+self.process_id+".fasta"));
		os.system("rm %s" % (self.PATH_TMP+"blastout_"+self.process_id));

		return (go_list, most_similar);

	# reads in clusters as dictionary GO-Id = key, entry=nr of cluster
	# stops at line "attr_class....", so no anti-cluster will be read
	def readClusterFile(self):
		f = open(self.PATH+self.CLUSTER_FILE);
		line = f.readline();
		clusters = {};
		while line:
			if line[0:11] == "clust_class":
				nr = int(line[12:14]);
				line_s = line.split("]");
				line_s = line_s[0].split("[");
				line_s = line_s[1].split(",");
				for i in range(len(line_s)):
					line_s[i] = line_s[i].strip(" ");
					line_s[i] = line_s[i].strip("X");
					clusters[line_s[i]] = nr;
			if line[0:10] == "attr_class":
				break;
			line = f.readline();
		#print clusters;
		return clusters;

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

	def __getTopScoring(self, found_go_clusters):
		go_clusters = {};
		for elem in found_go_clusters:
			if not elem in go_clusters.keys():
				go_clusters[elem] = 1;
			else:
				go_clusters[elem] = go_clusters[elem] + 1;

		top_hit = 0;
		top_score = 0;
		for key in go_clusters.keys():
			if go_clusters[key] > top_score:
				top_score = go_clusters[key];
				top_hit = key;
		return top_hit;

	def getGoFeatures(self,sequences, select = [] ):
		# read-in go map
		gomap = {};
		goterms = [];

		f = open(self.PATH+"uniprot2go.map","r");
		line = f.readline();
		while line:
			line_s = line.split(" ");
			gomap[line_s[0]] = line_s[1:-1];
			for elem in line_s[1:-1]:
				goterms.append(elem);
			line = f.readline();
		goterms = self.__uniqify(goterms);
		f.close();

		#get go terms
		(go_list, most_similar) = self.__getGoTerms(gomap, sequences);

		description_list = [];
		#create features
		if len(select) == 0:
			empty = [];
			description_empty = [];
			for elem in goterms:
				empty.append(0);
				description_empty.append("");
			feature_list = [];
			#print(sequences.size());
			#print(len(go_list));
			for i in range(sequences.size()):
				description_list.append(copy.copy(description_empty));
				feature_list.append(copy.copy(empty));
				for elem in go_list[i]:
					(go,score) = elem;
					k = goterms.index(go);
					feature_list[-1][k] = score;

		else:
			# read in go cluster file
			clusters = self.readClusterFile();
			goterms = select;
			empty = [];
			description_empty = [];
			for elem in goterms:
				empty.append(0);
				description_empty.append("");
			feature_list = [];
			#print(sequences.size());
			#print(len(go_list));
			#print(go_list);
			for i in range(sequences.size()):
				description_list.append(copy.copy(description_empty));
				feature_list.append(copy.copy(empty));
				adds = [];
				for elem in go_list[i]:
					(go,score) = elem;
					if go in goterms:
						k = goterms.index(go);
						feature_list[-1][k] = score;
						description_list[-1][k] = "GO:"+str(go);
						###check go clusters
					if go in clusters.keys():
						adds.append(clusters[go]);
						k = goterms.index("go_one");
						if not str(go) in description_list[-1][k]:
							description_list[-1][k] += "<a href=http://www.ebi.ac.uk/ego/GTerm?id=GO:"+str(go)+">GO:"+str(go)+"</a> ";
				# uniquify adds and add cluster
				if "go_one" in goterms:
					k = goterms.index("go_one");
					#print(adds);
					top_hit = self.__getTopScoring(adds);
					adds = self.__uniqify(adds);
					#print(adds);
					sum = 0;
					for elem in adds:
						sum += elem;
					#print(sum)
					#feature_list[-1][k] = sum;
					feature_list[-1][k] = top_hit
					#description_list[-1][k] = "GO???";

		return (goterms, feature_list, description_list, most_similar);
