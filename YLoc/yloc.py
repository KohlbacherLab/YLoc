from aasequences import *;
import config;
from predictor import *;
from featurecreator import *;
from prediction import *;

if os.getenv("YLOC_WS"):
    from dbsupport import *;

import models;

import copy;
import md5;
import os;
import random;
import time;


os.environ['PYTHON_EGG_CACHE'] = config.PATH_TMP;


class YLoc(object):

    def __init__(self, aa_sequences, id = "1"):
        if id == "1":
            id = str(random.random());
            md5object = md5.new(id);
            id = md5object.hexdigest();

        self.process_id = id;
        for i in range(aa_sequences.size()):
            (name, seq) = aa_sequences.get(i);
            name = name.replace("'"," ");
            name = name.replace("\""," ");
            aa_sequences.set((name,seq),i);
        self.__aa_sequences = aa_sequences;

        self.__predictor = Predictor(aa_sequences, self.process_id);
        self.__result = None;


    def getIP(self):
        sql_string = "SELECT ip FROM pending WHERE id='"+str(self.process_id)+"';";
        res=db_query(sql_string);

        if len(res) > 0:
            return str(res[0][0]);

        return "";


    def readFastaFile(self,fasta_name):
        self.__aa_sequences.clear();
        self.__aa_sequences.read_fasta_file(fasta_name);


    def predictLocation(self, model_name, ip=""):
        result = self.__predictor.predict(model_name);
        c = 0;

        if os.getenv("YLOC_WS"):
            ip = self.getIP();

            for elem in result:
                sql_string = "Insert Into queries (id,prediction,ip) VALUES ('"+str(self.process_id)+"',"+"\""+str(elem.getObjectAsList())+"\",'"+str(ip)+"');";
                res=db_query(sql_string);

                c += 1;
                alterJobStatus(self.process_id, 1, False, c);

            alterJobStatus(self.process_id, 1, False);

        self.__result = result;
        return result;


    def __printHeader(self):
        print("YLoc Prediction Result \n");
        print("========================== \n");
        print("PROCESS-ID "+self.process_id+ "\n");


    def printPredictionResultSimple(self):
        result = self.__result;
        self.__printHeader();
        for nr in range(0, len(result)):
            print("=== Prediction for sequence: "+str(nr)+"\t( "+self.__aa_sequences.get(nr)[0]+" ) === ");
            print("Location:\t"+str(result[nr].class_names[result[nr].location])+"\t"+str(result[nr].probability_distribution[result[nr].location]*100)+"%");
            print("\n");


    def __printDetailedPrediction(self, one_result):
        top = 10;
        print("\tLocation\t\tProbability\t  \n");
        loc_list = one_result.getSortedLocationList();
        line_string = "";
        for elem in loc_list:
            line_string += "\t"+str(one_result.class_names[elem])+"\t"+str(one_result.probability_distribution[elem]*100)+"%\t\n";
        print(line_string);
        # trust
        print("\tTrust value: "+str(one_result.null_prob));
        # print attribute prob table
        s = "\tAttribute Nr";
        for k in range(len(one_result.class_names)):
            s += "\t" + one_result.class_names[k];
        print(s);
        for i in range(one_result.attr_nr):
            #s = "\t" + str(one_result.feature_names[i]);
            s = "\t" + str(i);
            for k in range(len(one_result.class_names)):
                s += "\t" + one_result.probability_quality[i][k];
            print(s);
        print;
        # print attributes orderd by importance
        print("\tAttribute List\t\t Value \t\tEvidence\n");
        attr_list = one_result.getSortedAttributeList();

        all_e = 0;
        for elem in one_result.discrimination_scores:
            all_e += abs(elem);

        evidence = 0.0;
        c = 0;
        while evidence < 0.8*all_e:
            e = one_result.discrimination_scores[attr_list[c]];
            evidence += abs(e);
            true_class_string = one_result.getClassNameString(one_result.merged_true_classes[attr_list[c]]);
            #print one_result.feature_names;
            #print one_result.discretization_intervals;
            print("\t"+str(one_result.feature_names[attr_list[c]])+"\t\t"+str(one_result.discretization_intervals[attr_list[c]][one_result.matched_intervals[attr_list[c]]])+"\t\tdiscr score: "+str(e));
            print("\tExplanation: "+str(one_result.feature_description[attr_list[c]][-2])+" "+str(one_result.feature_description[attr_list[c]][-1]));
            if len(one_result.pos_classes[attr_list[c]]) > 0:
                print("\tTypical for "+ true_class_string + "("+str(one_result.true_class_attribute_distribution[attr_list[c]][one_result.matched_intervals[attr_list[c]]]*100)+"%) compared to "+str(one_result.getClassNameString(one_result.pos_classes[attr_list[c]])))+ "("+str(one_result.pos_class_attribute_distribution[attr_list[c]][one_result.matched_intervals[attr_list[c]]]*100)+"%)";
            if len(one_result.neg_classes[attr_list[c]]) > 0:
                print("\tUntypical for "+ true_class_string + "("+str(one_result.true_class_attribute_distribution[attr_list[c]][one_result.matched_intervals[attr_list[c]]]*100)+"%) compared to "+str(one_result.getClassNameString(one_result.neg_classes[attr_list[c]])))+ "("+str(one_result.neg_class_attribute_distribution[attr_list[c]][one_result.matched_intervals[attr_list[c]]]*100)+"%)";
            c += 1;
            print;

        """
        for attr_nr in range(2,top+2):
            (name,prob) = one_result[class_nr][attr_nr];
            line_string += name+" ["+str(prob*100)+"%]\n \t\t\t\t\t";
        """


    def printPredictionResult(self):
        result = self.__result;
        self.__printHeader();
        for nr in range(len(result)):
            print("=== Prediction for sequence: "+str(nr)+"\t( "+self.__aa_sequences.get(nr)[0]+" ) === ");
            s = "Locations:\t";
            for elem in result[nr].locations:
                s += str(result[nr].class_names[elem])+"("+str(result[nr].probability_distribution[elem]*100)+"%)\t";
            print(s);
            print("Details:");
            self.__printDetailedPrediction(result[nr]);
            print("\n");


if os.getenv("YLOC_WS"):
    sequence_file = sys.argv[1];
    model_name = sys.argv[2];
    id = sys.argv[3];

    aa = AASequences();
    aa.read_fasta_file(sequence_file);

    m = YLoc(aa, id);
    m.predictLocation(model_name);

    m.printPredictionResult();
else:
    if len(sys.argv)>2:
        model = models.Models();
        exit = False;
        print_results = False;

        sequence_file = sys.argv[1];
        model_name = sys.argv[2];
        if not model_name in model.getAvailableModelsOutside():
            exit = True;

        if len(sys.argv)>3:
            id = sys.argv[3];
            if len(sys.argv)>4:
                print_results = True;
        else:
            print_results = True;
            id = "1";

        if exit:
            print "Usage YLoc: python yloc.py <fasta_file> <model_name> <prediction_id(optional)> <print_result(y/n)(optional)>"
            print
            m = models.Models();
            print "Available models:"
            models = m.getAvailableModelsOutside();
            for elem in models:
                print elem;
        else:
            aa = AASequences();
            aa.read_fasta_file(sequence_file);

            m = YLoc(aa, id);

            name = model.getModelNameIntern(model_name);
            print name;
            m.predictLocation(name);

            if print_results:
                m.printPredictionResult();
    else:
        print "Usage YLoc: python yloc.py <fasta_file> <model_name> <prediction_id(optional)> <print_result(y/n)(optional)>"
        print
        m = models.Models();
        print "Available models:"
        models = m.getAvailableModelsOutside();
        for elem in models:
            print elem;
