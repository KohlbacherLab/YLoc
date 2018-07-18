import os, sys, re;
from math import sqrt;

from featurecreator import *;
from aasequences import *;
import config;
from prediction import *;

class WekaInterface(object):

	"Class that communicates with the weka toolbox"
	
	__PATH_TMP = config.PATH_TMP;
	__PATH_WEKA = config.PATH_WEKA_LEARN;
	__PATH_UNBIAS_WEKA = config.PATH_WEKA_LEARN;
	__PATH_WEKA_PREDICT = config.PATH_WEKA_PREDICT;
	__PATH_WEKA_MCC_OPTIMIZER = config.PATH_WEKA_LEARN;
	
	def __init__(self, id = "1"):
		self.process_id = id;
				

	def normalizeFeatures(self, features):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");
		old_first_element = features.getFeatures(0);
		os.system("java -classpath %s -Xmx1500m weka.filters.unsupervised.attribute.Normalize  -i %s -o %s" % (self.__PATH_WEKA,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",self.__PATH_TMP+"tmp_norm_weka_"+self.process_id+".arff"));
		new_features = FeatureCreator();
		new_features.readARFFFeatureFile(self.__PATH_TMP+"tmp_norm_weka_"+self.process_id+".arff");
		new_first_element = new_features.getFeatures(0);
		norm_coefficients = [];
		for i in range(len(old_first_element)):
			norm_coefficients.append( float(new_first_element[i]) / float(old_first_element[i]) );
		return new_features;
	
	def runFeatureSelection(self, features):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");		
		os.system("java -classpath %s -Xmx3000m weka.filters.supervised.attribute.AttributeSelection -E \"weka.attributeSelection.CfsSubsetEval -M\" -S \"weka.attributeSelection.BestFirst -D 2 -N 5 -S 10 \" -i %s -o %s" % (self.__PATH_WEKA,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",self.__PATH_TMP+"tmp_sel_weka_"+self.process_id+".arff"));
		new_features = FeatureCreator();
		new_features.readARFFFeatureFile(self.__PATH_TMP+"tmp_sel_weka_"+self.process_id+".arff");
		return new_features;
		
	def runNBFeatureSelection(self, features,steps, mcc_opt=False, unbias = False):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");		
		weka_path = self.__PATH_WEKA;
		if unbias:
			weka_path = self.__PATH_UNBIAS_WEKA;
		if mcc_opt:
			weka_path = self.__PATH_WEKA_MCC_OPTIMIZER;
		os.system("java -classpath %s -Xmx2000m weka.filters.supervised.attribute.AttributeSelection -E \"weka.attributeSelection.WrapperSubsetEval -B weka.classifiers.bayes.NaiveBayes -F 5 -T 0.01 -R 1 --\" -S \"weka.attributeSelection.GeneticSearch -Z 20 -G %s -C 0.6 -M 0.033 -R 20 -S 1\" -i %s -o %s" % (weka_path,steps,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",self.__PATH_TMP+"tmp_selnb_weka_"+self.process_id+".arff"));
		new_features = FeatureCreator();
		new_features.readARFFFeatureFile(self.__PATH_TMP+"tmp_selnb_weka_"+self.process_id+".arff");
		return new_features;

	def runNBSimpleFeatureSelection(self, features,steps, mcc_opt=False, unbias = False):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");		
		weka_path = self.__PATH_WEKA;
		if unbias:
			weka_path = self.__PATH_UNBIAS_WEKA;
		if mcc_opt:
			weka_path = self.__PATH_WEKA_MCC_OPTIMIZER;
		os.system("java -classpath %s -Xmx5000m weka.filters.supervised.attribute.AttributeSelection -E \"weka.attributeSelection.WrapperSubsetEval -B weka.classifiers.bayes.NaiveBayesSimple -F 5 -T 0.01 -R 1 --\" -S \"weka.attributeSelection.GeneticSearch -Z 20 -G %s -C 0.6 -M 0.033 -R 20 -S 1\" -i %s -o %s" % (weka_path,steps,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",self.__PATH_TMP+"tmp_selnbsimple_weka_"+self.process_id+".arff"));
		new_features = FeatureCreator();
		new_features.readARFFFeatureFile(self.__PATH_TMP+"tmp_selnbsimple_weka_"+self.process_id+".arff");
		return new_features;
		
	def runLogisticFeatureSelection(self, features,steps):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");		
		os.system("java -classpath %s -Xmx3500m weka.filters.supervised.attribute.AttributeSelection -E \"weka.attributeSelection.WrapperSubsetEval -B weka.classifiers.functions.Logistic  -F 5 -T 0.01 -R 1 -- -R 1.0E-8 -M -1\" -S \"weka.attributeSelection.GeneticSearch -Z 20 -G %s -C 0.6 -M 0.033 -R 20 -S 1\" -i %s -o %s" % (self.__PATH_WEKA,steps,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",self.__PATH_TMP+"tmp_sellog_weka_"+self.process_id+".arff"));
		new_features = FeatureCreator();
		new_features.readARFFFeatureFile(self.__PATH_TMP+"tmp_sellog_weka_"+self.process_id+".arff");
		return new_features;
	
	def runBNFeatureSelection(self, features,steps):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");		
		os.system("java -classpath %s -Xmx1500m weka.filters.supervised.attribute.AttributeSelection -E \"weka.attributeSelection.WrapperSubsetEval -B weka.classifiers.bayes.BayesNet  -F 5 -T 0.01 -R 1 -- -D -Q weka.classifiers.bayes.net.search.local.K2 -- -P 1 -S BAYES -E weka.classifiers.bayes.net.estimate.SimpleEstimator -- -A 0.5\" -S \"weka.attributeSelection.GeneticSearch -Z 20 -G %s -C 0.6 -M 0.033 -R 20 -S 1\" -i %s -o %s" % (self.__PATH_WEKA,steps,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",self.__PATH_TMP+"tmp_selbn_weka_"+self.process_id+".arff"));
		new_features = FeatureCreator();
		new_features.readARFFFeatureFile(self.__PATH_TMP+"tmp_selbn_weka_"+self.process_id+".arff");
		return new_features;
	
	def __calculateAverageResult(self, result_list):
		times = len(result_list);
		mean_v = [];
		sd_v = [];
		for k in range(0,len(result_list[0])):
				mean = 0;
				for i in range(0,times):
					mean += result_list[i][k];
				mean = mean / times;
				sd = 0;
				for i in range(0,times):
					sd += (mean-result_list[i][k])*(mean-result_list[i][k]);
				sd = sqrt(sd);
				mean_v.append(mean);
				sd_v.append(sd);
		return [mean_v,sd_v];
	
	def crossValidSeriesFeaturesNB(self, features, fold):
		times = 5;
		res = [self.crossValidFeaturesNB(features,fold)];
		for i in range(2,times+1):
			res.append(self.crossValidFeaturesNB(features,fold,i));

		return self.__calculateAverageResult(res);
	
	def crossValidSeriesFeaturesDiscretizedNB(self, features, fold):
		times = 5;
		res = [self.crossValidFeaturesDiscretizedNB(features,fold)];
		for i in range(2,times+1):
			res.append(self.crossValidFeaturesDiscretizedNB(features,fold,i));

		return self.__calculateAverageResult(res);
	
	def crossValidSeriesFeaturesNBSimple(self, features, fold):
		times = 5;
		res = [self.crossValidFeaturesNBSimple(features,fold)];
		for i in range(2,times+1):
			res.append(self.crossValidFeaturesNB(features,fold,i));

		return self.__calculateAverageResult(res);
	
	# gets name of weka classifier output file and parses this file
	# returns list with [accuracy, average mcc, mcc_class_1,...,mcc_class_k]
	def __getStatisticFromClassifierOutput(self, output_file_name):
		f = open(output_file_name,"r");
		line = f.readline();
		
		acc = 0.0;
		inst = 0;
		
		while(line):
			if (line[0:30]=="Correctly Classified Instances"):
				line  = line[57:-1];
				line_s = line.split(" ");
				acc = float(line_s[0].strip(" "));
			line = f.readline();
			if (line[0:25]=="Total Number of Instances"):
				line = line[36:-1]
				inst = int(line);
			if (line[0:24]=="=== Confusion Matrix ==="):
				line = f.readline();				
				line = f.readline();
				line_s = line.split("<");
				line = line_s[0].replace(" ","");
				classes = len(line);
				confusion_matrix = [];
				line = f.readline();
				for i in range(0,classes):
					confusion_line = [];
					re_result = re.findall("[1234567890]{1,}\s",line);
					for elem in re_result:
						confusion_line.append(float(elem));
					confusion_matrix.append(confusion_line);
					re_result2 = re.findall("[1234567890]$",line);
					line = f.readline();

				mcc = [];
				sum_mcc = 0;
				for i in range(0,classes):
					tp = confusion_matrix[i][i];
					fn = -1* tp;
					for k in confusion_matrix[i]:
						fn += k;
					fp = -1*tp;
					for k in range(0,classes):
						fp += confusion_matrix[k][i];
					tn = inst - tp - fn - fp;
					mcc_class = (tp*tn - fp*fn) / sqrt( (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) );
					sum_mcc += mcc_class;
					mcc.append(mcc_class);
				sum_mcc = sum_mcc / classes;
		
		return [acc,sum_mcc]+mcc;
	

	def crossValidFeaturesNB(self, features, fold,seed = 1):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");		
		os.system("java -classpath %s -Xmx1500m weka.classifiers.bayes.NaiveBayes -v -s %s -t %s -x %s > %s" % (self.__PATH_WEKA,seed,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",fold,self.__PATH_TMP+"weka_"+self.process_id+".out"));	
		return self.__getStatisticFromClassifierOutput(self.__PATH_TMP+"weka_"+self.process_id+".out");

	def crossValidFeaturesDiscretizedNB(self, features, fold,seed = 1):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");		
		os.system("java -classpath %s -Xmx1500m weka.classifiers.meta.FilteredClassifier -F \"weka.filters.supervised.attribute.Discretize -R first-last\" -W \"weka.classifiers.bayes.NaiveBayes\" -v -s %s -t %s -x %s > %s" % (self.__PATH_WEKA_PREDICT,seed,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",fold,self.__PATH_TMP+"weka_"+self.process_id+".out"));	
		return self.__getStatisticFromClassifierOutput(self.__PATH_TMP+"weka_"+self.process_id+".out");

		
	def crossValidFeaturesNBSimple(self, features, fold,seed = 1):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");		
		os.system("java -classpath %s -Xmx1500m weka.classifiers.bayes.NaiveBayesSimple -v -s %s -t %s -x %s > %s" % (self.__PATH_WEKA,seed,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",fold,self.__PATH_TMP+"weka_"+self.process_id+".out"));	
		return self.__getStatisticFromClassifierOutput(self.__PATH_TMP+"weka_"+self.process_id+".out");
		
	def createNBModel(self, features, model_file):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");	
		print("wrote mp");	
		os.system("java -classpath %s -Xmx1500m weka.classifiers.bayes.NaiveBayes -v -t %s -d %s > %s" % (self.__PATH_WEKA,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",model_file,self.__PATH_TMP+"weka_"+self.process_id+".out"));	

	def createNBSimpleModel(self, features, model_file):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");	
		print("wrote mp");	
		os.system("java -classpath %s -Xmx1500m weka.classifiers.bayes.NaiveBayesSimple -v -t %s -d %s > %s" % (self.__PATH_WEKA,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",model_file,self.__PATH_TMP+"weka_"+self.process_id+".out"));	

	def createDiscretizedNBModel(self, features, model_file):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");	
		os.system("java -classpath %s -Xmx1500m weka.classifiers.meta.FilteredClassifier -F \"weka.filters.unsupervised.attribute.Normalize\" -W \"weka.classifiers.bayes.NaiveBayesStdOut\" -v -t %s -d %s > %s" % (self.__PATH_WEKA_PREDICT,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",model_file,self.__PATH_TMP+"weka_"+self.process_id+".out"));
	
	def createNormalizedNBSimpleModel(self, features, model_file):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");	
		os.system("java -classpath %s -Xmx1500m weka.classifiers.meta.FilteredClassifier -F \"weka.filters.unsupervised.attribute.Normalize\" -W \"weka.classifiers.bayes.NaiveBayesSimple\" -v -t %s -d %s > %s" % (self.__PATH_WEKA,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",model_file,self.__PATH_TMP+"weka_"+self.process_id+".out"));
		
	def createDiscretizeNBModel(self, features, model_file):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");	
		os.system("java -classpath %s -Xmx1500m weka.classifiers.meta.FilteredClassifier -F \"weka.filters.supervised.attribute.Discretize -R first-last \" -W \"weka.classifiers.bayes.NaiveBayesStdOut\" -v -t %s -d %s > %s" % (self.__PATH_WEKA_PREDICT,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",model_file,self.__PATH_TMP+"weka_"+self.process_id+".out"));
		
	def createDiscretizeMLNBModel(self, features, model_file, first_dual, class_keys):
		class_keys_string = "";
		for elem in class_keys:
			class_keys_string += str(elem) +",";
		class_keys_string = class_keys_string[:-1];
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");	
		os.system("java -classpath %s -Xmx1500m weka.classifiers.meta.FilteredClassifier -v -t %s -d %s -F \"weka.filters.supervised.attribute.Discretize -R first-last \" -W weka.classifiers.bayes.MultiLabelNaiveBayesStdOut -- -F %s -C \"%s\"   > %s" % (self.__PATH_WEKA_PREDICT, self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",model_file, first_dual, class_keys_string,self.__PATH_TMP+"weka_"+self.process_id+".out"));	
	
	
	def __makeNumberCompact(self, s):
		if s != "inf" and s != "-inf":
			p = s.find(".");
			if p != -1 and len(s)>p+3:
				return s[:p+3];
				
		return s;
	
	def __makeDiscrCompact(self, discr_string):
		
		div=discr_string.find("-");
		if div == 0:
			div=discr_string[1:].find("-")+1;
		first = discr_string[0:div];
		second = discr_string[div+1:];
		new = self.__makeNumberCompact(first) + " to " + self.__makeNumberCompact(second);
		return new;
			
	def __parseSebsProbabilities(self, weka_out_file, first_line):
		p = Prediction();
		## get probs from first line
		probs = [];
		line_s = first_line[9:].strip("\r\n ").split(" ");
		for elem in line_s:
			probs.append( float(elem));
		p.probability_distribution = probs;
		## get null prob
		line = weka_out_file.readline();
		np = float(line[13:]);
		p.null_prob = np;
		#print line[13:];
		#print np;
		## read true classes
		line = weka_out_file.readline();
		line_s = line[16:].strip("\r\n ").split(" ");
		locs = [];
		for elem in line_s:
			locs.append( int(elem) );
		p.locations = locs;
		## read putative classes
		line = weka_out_file.readline();
		line_s = line[21:].strip("\r\n ").split(" ");
		if line_s[0] == '':
			line_s = [];
		locs = [];
		for elem in line_s:
			locs.append( int(elem) );
		p.putative_locations = locs;
		## number of attributes
		line = weka_out_file.readline();
		p.attr_nr = int(line[25:]);
		### get attr probability table and signed results
		for i in range(p.attr_nr):
			line = weka_out_file.readline();
			line = weka_out_file.readline();
			line_s = line.strip(" \r\n").split(" ");
			probs = [];
			for elem in line_s:
				probs.append( float(elem));
			p.probability_quantity.append(probs);
			line = weka_out_file.readline();
			line_s = line.strip(" \r\n").split(" ");
			p.probability_quality.append(line_s);
		### read attr discr scores
		line = weka_out_file.readline();
		line_s = line[35:].strip(" \n\r").split(" ");
		attrDiscrScore = [];
		for elem in line_s:
			attrDiscrScore.append( float (elem));
		p.discrimination_scores = attrDiscrScore;
		### read discritization intervals
		discretizationI = [];
		for i in range(len(attrDiscrScore)):
			line = weka_out_file.readline();
			line_s = line[32:].strip(" \n\r").split(" ");
			for k in range(len(line_s)):
				line_s[k] = self.__makeDiscrCompact(line_s[k][2:-2]);
			discretizationI.append(line_s);
		p.discretization_intervals = discretizationI;
		### prepare each attribute by reading true interval and distributions
		line = weka_out_file.readline();
		attrNr = -1;
		while line:
			if line[0:16] == "***True Interval":
				attrNr += 1;
				# get true interval
				interval = int( line[17:] );
				p.matched_intervals.append(interval);
				# get true merged classes
				line = weka_out_file.readline();
				line_s = line[24:].strip("\n\r ").split(" ");
				merged_classes = [];
				for elem in line_s:
					merged_classes.append( int(elem) );
				p.merged_true_classes.append(merged_classes);
				# get true distribution
				line = weka_out_file.readline();
				line_s = line[22:].strip("\n\r ").split(" ");
				true_distr = [];
				for elem in line_s:
					true_distr.append( float(elem) );
				p.true_class_attribute_distribution.append(true_distr);
				# check if pos or neg distribution
				line = weka_out_file.readline();
				if line[0:15] == "***Pos Classes:":
					# add pos classes
					line_s = line[16:].strip(" \n\r").split(" ");
					pos_classes = [];
					for elem in line_s:
						pos_classes.append( int(elem));
					p.pos_classes.append(pos_classes);
					# add pos distribution
					line = weka_out_file.readline();
					line_s = line[21:].strip(" \n\r").split(" ");
					pos_distr = [];
					for elem in line_s:
						pos_distr.append( float(elem) );
					p.pos_class_attribute_distribution.append(pos_distr);
					line = weka_out_file.readline();
				else:
					p.pos_classes.append([]);
					p.pos_class_attribute_distribution.append([]);
				if line[0:15] == "***Neg Classes:":
					# add pos classes
					line_s = line[16:].strip(" \n\r").split(" ");
					neg_classes = [];
					for elem in line_s:
						neg_classes.append( int( elem) );
					p.neg_classes.append(neg_classes);
					# add pos distribution
					line = weka_out_file.readline();
					line_s = line[21:].strip(" \n\r").split(" ");
					neg_distr = [];
					for elem in line_s:
						neg_distr.append( float(elem) );
					p.neg_class_attribute_distribution.append(neg_distr);
					line = weka_out_file.readline();
				else:
					 p.neg_classes.append([]);
					 p.neg_class_attribute_distribution.append([]);
			if line[0:3] == "###":
				line = False;
			else:		
				line = weka_out_file.readline();
				
		return p;

	def __parsePredictionResults(self, weka_output_name):
		f = open(weka_output_name,"r");
		line = f.readline();

		all_results = [];
		class_names = [];

		c = 0;
		while line:
			# get probabilities
			if line[0:8]=="###Probs":
				predict_result = self.__parseSebsProbabilities(f, line);
				predict_result.sequence_id = c;
				all_results.append(predict_result);
				c += 1;
			line = f.readline();

		return all_results;		
			
	def predictModel(self, model_file, features):
		features.writeARFFFeatureFile(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff");
		s = "java -classpath %s -Xmx1500m weka.classifiers.meta.FilteredClassifier -p -l %s -T %s > %s" % (self.__PATH_WEKA_PREDICT,model_file,self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff",self.__PATH_TMP+"weka_"+self.process_id+".out");
		os.system(s);	
		
		results = self.__parsePredictionResults(self.__PATH_TMP+"weka_"+self.process_id+".out");
		#results = self.__parsePredictionResults(self.__PATH_TMP+"weka.out");
		
		os.system( "rm %s"  %(self.__PATH_TMP+"weka_"+self.process_id+".out"));
		os.system( "rm %s"  %(self.__PATH_TMP+"tmp_weka_"+self.process_id+".arff"));
		
		return results;


		
#java -classpath weka/weka.jar -Xmx4500m weka.filters.supervised.attribute.AttributeSelection -E "weka.attributeSelection.CfsSubsetEval -M" -S "weka.attributeSelection.BestFirst -D 1 -N 5" -i test_large3.arff -o test_large3_sel.arff &
#weka.filters.supervised.attribute.AttributeSelection -E "weka.attributeSelection.WrapperSubsetEval -B weka.classifiers.meta.FilteredClassifier -F 5 -T 0.01 -R 1 -- -F \"weka.filters.supervised.attribute.Discretize -E -K -R first-last\" -W weka.classifiers.bayes.NaiveBayes --" -S "weka.attributeSelection.GeneticSearch -Z 20 -G 50 -C 0.6 -M 0.033 -R 20 -S 1"
