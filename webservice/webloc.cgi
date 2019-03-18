#!/usr/bin/python

import sys, os, re, cgi, time, traceback
from subprocess import PIPE, Popen

sys.stderr = sys.stdout

try:
  sys.path.append('/YLoc');

  os.environ['PYTHON_EGG_CACHE'] = '/var/www/html/cgi-bin/';

  os.environ['YLOC_WS'] = '1';

  python_path = "/usr/bin/python";

  img_path = "images/";
  download_path = "downloads/";

  from ylconfig import *

  import cgitb;
  from dbsupport import *;
  from aasequences import *;
  from prediction import *;
  from math import *;
  from models import *;
  import config;
  import random;
  import md5;
  import re;
  cgitb.enable();
except:
  print "\n\n<PRE>"
  traceback.print_exc()


def __print_header(refresh=0,id=1):
  print "Content-type: text/html\n"
  print "<HTML><HEAD>";

  if refresh >0:
    print "<meta http-equiv='refresh' content='"+str(refresh)+"; URL=webloc.cgi?id="+str(id)+"'>"

  print "<TITLE>YLoc</TITLE>"
  print "<style type='text/css'>";
  print " h1 { color:#734D38; font-family:Verdana; font-variant:small-caps; font-size:30pt; font-weight:bold; }  "
  print " h2 { color:#000000; font-family:Verdana; font-size:13pt; font-weight:bold; }  "
  print " h3 { color:#000000; font-family:Verdana; font-size:10pt; font-weight:bold; }  "
  print " #navi { margin: 0px 10px; height:32px; background-color:#FFFFFF; border:0px; border-top:1px; border-style:solid; border-color:#000000;}"
  print " #navi a { color:#FFFFFF; font-size:10pt; background-color: #003380; padding:0px 18px; padding-bottom:3px; margin:0pt; border:1px; border-color:#000000; border-style:solid; text-decoration:none;margin-right:3pt;}"
  print " #navi a:hover {  background-color: #aaccff; color:#000000; border-color:#000000;  }"
  print " #navi .current { background-color: #cbdbff; color:#000000; border-color:#0F0F0F;}"
  print " #subnavi { padding-bottom:5px; padding-top:0px; text-align:right; margin-left:0px;width:190px;background-color:#FFFFFF; border:0px; border-top:2px; border-bottom:2px; border-color:#B38D68; border-style:solid;} ";
  print " #subnavi ul { margin:0px;padding:0px; text-align:right;list-style:none; width:140px; }";
  print " #subnavi li { text-align:right; }";
  print " #subnavi li a {  font-size:10pt; text-align:left; padding:2px; margin-left:20px; height:15px; width:140px; background:#EEEEEE; display:block; margin-top:5px; text-decoration:none; border:0px; border-left:4px; border-color:#CCCCCC; border-style:solid;}"
  print " #subnavi li a:hover { border-color:#777777;  }"
  print " #subnavi .current { background-color: #DDDDDD; color:#444444; }";
  print " a { font-family:Verdana; color:#000000; } ";
  print " img { border:none; }";
  print " hr { border:0px solid #000000; height:1px; background-color:#000000;}";
  print " a.helptext { color:#000000; } ";
  print " p.underl a { color:#000000; text-decoration:underline;} ";
  print " th.underl a {text-decoration:underline;};"

  print " p { font-family:Verdana; vertical-align:top; font-size:10pt; font-weight:normal; color:#000000; text-decoration:none; padding:3px; border:0px; text-align:left; } ";
  print " div { font-family:Verdana; vertical-align:top; font-size:10pt; font-weight:normal; color:#000000; text-decoration:none; padding:0px; border:0px; text-align:left; } ";
  print " .input { font-family:Verdana; vertical-align:top; font-size:10pt; font-weight:normal; color:#000000; text-decoration:none; padding:3px; text-align:left; } ";
  print " .button {  background-color:#cbdbff; } ";
  print " .button2 { cursor:pointer; background-color:#dee7ec; font-family:Verdana; vertical-align:top; font-size:10pt; font-weight:normal; color:#000000; text-decoration:none;  border:0px solid #000000; padding:2px; margin2px; text-align:center;} ";
  print " td { font-family:Verdana; vertical-align:top; font-size:10pt; font-weight:normal; color:#000000; text-decoration:none; padding:3px; border:0px; text-align:left; } ";
  print " th { font-family:Verdana; vertical-align:top; font-size:10pt; font-weight:normal; color:#000000; text-decoration:none; padding:20px; border:0px; text-align:left; } ";
  print " th a {  font-size:10pt; color:#000000; text-decoration:none; } ";
  print " .headerprediction { vertical-align:center; height:10pt; font-size:8pt; color:#000000; text-decoration:none; font-weight:normal; text-align:left;} ";
  print " .predict1 {  vertical-align:center; background-color:#CCE6FF; height:12pt; font-size:10pt; color:#000000; text-decoration:none; font-weight:normal; text-align:left;} ";
  print " .attr1 {   background-color:#EEEEFF; font-size:12pt; color:#000000; text-decoration:none; font-weight:normal; text-align:left;} ";
  print " .predict2 { vertical-align:center;  background-color:#DDF6FF; height:12pt;  font-size:10pt; color:#222222; text-decoration:none; font-weight:normal; text-align:left;} ";
  print " .attr2 {  background-color:#EEEEFF; font-size:10pt; color:#222222; text-decoration:none; font-weight:normal; text-align:left;} ";
  print " .tablematrix { font-family:Arial; color:#000000; background-color:#FFFFFF; vertical-align:top;  font-size:10pt; padding: 1px 1px; border:0px; }";
  print " .amatrix1 { vertical-align:center;  background-color:#EEEEFE; height:12pt;  font-size:10pt; color:#222222; text-decoration:none; font-weight:normal; text-align:left; padding:1px;} ";
  print " .amatrix2 { vertical-align:center;  background-color:#E2E2F2; height:12pt;  font-size:10pt; color:#222222; text-decoration:none; font-weight:normal; text-align:left; padding:1px;} ";
  print " .amatrix3 { vertical-align:center;  background-color:#FFFFFF; height:12pt;  font-size:10pt; color:#222222; text-decoration:none; font-weight:normal; text-align:left; padding:1px; border:3px solid #D2D2D2;} ";
  print " .amatrixside { vertical-align:center;  background-color:#cbdbff;    font-weight:bold; font-size:10pt; color:#222222; text-decoration:none;  text-align:left; } ";
  print " .empty {  background-color:#FFFFFF; font-size:10pt; color:#444444; text-decoration:none; font-weight:normal; text-align:left;} ";
  print " .thpredict {  background-color:#EEEEEE; font-size:12pt; color:#00000; text-decoration:none; font-weight:bold; padding-bottom: 2px; padding-top:2px; text-align:left;} ";
  print " .thheading {  background-color:#FFFFFF; font-size:11pt; color:#00000; text-decoration:none; font-weight:bold; padding-bottom: 2px; padding-top:2px; text-align:left;} ";
  #print " th.thheading span  { display:none;  position:absolute; left:30em;  } ";
  #print " th.thheading:hover span { display:block; background-color:#FFFFC5; color:#000000; font-size:10pt; cursor:help; text-decoration:none; font-weight:normal; cursor:help; width:400pt; border:1px solid #000000; padding:3px; } ";
  #print " th.thheading:hover {   cursor:help; } ";
  print " .thheading2 {  background-color:#FFFFFF; font-size:11pt; color:#00000; text-decoration:none; font-weight:bold; padding-bottom: 2px; padding-top:2px; text-align:left;} ";
  print " .thheading span {  color:#000000; font-size:10pt; cursor:help; text-decoration:none; font-weight:normal; cursor:help; padding:3px; } ";
  print " a.thheading2 span  { display:none;  position:absolute; left:30em;  } ";
  print " a.thheading2:hover span { display:block; background-color:#FFFFC5; color:#000000; font-size:10pt; cursor:help; text-decoration:none; font-weight:normal; cursor:help; width:400pt; border:1px solid #000000; padding:3px; } ";
  print " a.thheading2:hover {   cursor:help; } ";

  print " a.helptext span  { display:none;  position:absolute; } ";
  print " a.helptext:hover span { display:block; background-color:#FFFFC5; color:#000000; font-size:10pt; cursor:help; text-decoration:none; font-weight:normal; cursor:help; width:200pt; border:1px solid #000000; padding:3px; } ";
  print " a.helptext:hover {   cursor:help; } ";

  print " a.helptextsmall span  { display:none;  position:absolute; } ";
  print " a.helptextsmall:hover span { display:block; background-color:#FFFFC5; color:#000000; font-size:10pt; cursor:help; text-decoration:none; font-weight:normal; cursor:help; width:70pt; border:1px solid #000000; padding:3px; } ";
  print " a.helptextsmall:hover {   cursor:help; } ";

  print " td.probtable span  { display:none;  position:absolute; top:30em; right:6em;  } "
  print " td.probtable:hover span  { display:block; background-color:#FFFFFF; color:#000000; font-size:10pt; cursor:help; text-decoration:none; font-weight:normal; cursor:help; width:330pt;padding:3px; border:1px solid #cbdbff;} ";
  print " td.probtable:hover {   cursor:help; } ";
  print " a.probtable span  { display:none;  position:absolute; top:30em; right:6em;  } "
  print " a.probtable:hover span  { display:block; background-color:#FFFFFF; color:#000000; font-size:10pt; cursor:help; text-decoration:none; font-weight:normal; cursor:help; width:330pt;padding:3px; border:1px solid #cbdbff;} ";
  print " a.probtable:hover {   cursor:help; } ";

  print " td.attinfo span  { display:none;  position:absolute; left:5em; } "
  print " td.attinfo:hover span  { display:block; background-color:#FFFFC5; color:#000000; font-size:10pt; cursor:help; text-decoration:none; font-weight:normal; cursor:help; width:600pt;padding:3px; border:1px solid #000000; padding:3px;} ";
  print " td.attinfo:hover {   cursor:help; } ";
  print " a.attinfo span  { display:none;  position:absolute; left:5em;  } "
  print " a.attinfo:hover span  { display:block; background-color:#FFFFC5; color:#000000; font-size:10pt; cursor:help; text-decoration:none; font-weight:normal; cursor:help; width:600pt;padding:3px; border:1px solid #000000; padding:3px; } ";
  print " a.attinfo:hover {   cursor:help; } ";

  print " a.amatrix3 span  { display:none;  position:absolute; left:30em; } ";
  print " a.amatrix3:hover span { display:block; background-color:#FFFFC5; color:#000000; font-size:10pt; cursor:help; text-decoration:none; font-weight:normal; cursor:help; width:400pt; border:1px solid #000000; padding:3px; } ";
  print " a.amatrix3:hover {   cursor:help; } ";

  print " .amatrixlarge { font-size:12pt; padding:8px; padding-left:1px;}";
  print " .amatrixmedium {font-size:10pt; padding:3px; padding-left:1px;}";
  print " .amatrixsmall {font-size:8pt; padding:1px;}";


  print " .tablepredict { font-family:Arial; color:#000000; background-color:#FFFFFF; width:100%; vertical-align:top;  font-size:10pt; padding: 1px 1px; padding-left:0pt; border:0px;text-decoration:none;  }";
  print " .popup { color:#000000; cursor:help; text-decoration:none;}";


  print "</style>";
  print "<SCRIPT LANGUAGE='javascript'>";
  print "<!--";
  print "function openNewWindow(address) {";
  print "   x=window.open (address, 'new', config='height=650, width=800, toolbar=no, menubar=no, scrollbars=yes, resizable=yes, location=no, directories=no, status=no');"
  print "  x.focus();";
  print " };"
  print "function openHelpWindow(address) {";
  print "   x=window.open (address, 'help', config='height=600, width=800, toolbar=no, menubar=no, scrollbars=yes, resizable=yes, location=no, directories=no, status=no');"
  print "  x.focus();";
  print " };"
  print "function openExternalWindow(address) {";
  print "   x=window.open (address, 'external', config='height=600, width=800, toolbar=no, menubar=no, scrollbars=yes, resizable=yes, location=no, directories=no, status=no');"
  print "  x.focus();";
  print " };"
  print "function openNewPage(address) {";
  print "   this.document.open(address);"
  print " };"
  print "-->";
  print "</SCRIPT>";

  print "</HEAD></HTML>";

def __print_head():
  print "<BODY bgcolor=#DDDDDD>";
  print "<center>";
  print "<TABLE bgcolor=#FFFFFF  border=0 width = 95% height=95%>";
  print "<TR height=80pt><TD width=150 style='padding-left:20pt'><img src='"+str(img_path)+"yloc_1.png'></TD>";
  print "<TD width=100%  style='vertical-align:middle; padding:0px; padding-left:20pt;'><h2>Interpretable Subcellular Localization Prediction</h2></TD></TR>";
  print "<TR><TH colspan=2 width=100% height=40pt>";
  print "<div id='navi'>"

  class_list = ["","",""];

  if not form.has_key("page"):
    class_list[0] = "class='current'";
  elif form["page"].value=="help":
    class_list[1] = "class='current'";
  elif form["page"].value=="info":
    class_list[2] = "class='current'";
  else:
    class_list[0] = "class='current'";

  print "<a "+class_list[0]+" href='webloc.cgi'>Predict with YLoc</a>";
  print "<a "+class_list[1]+" href='webloc.cgi?page=help'>Tutorial</a>";
  print "<a "+class_list[2]+" href='webloc.cgi?page=info'>Information</a>";
  print "</div></TH></TR>"
  print "<TR><TH colspan=2 height=80% style='padding-left:30pt;'>";

def __print_foot():
  print "</TH></TR>";
  print "<TR><TH colspan=2 bgcolor=#FFFFFF>"
  print "<hr>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;<a style='font-size:8pt;color:blue' href=mailto:" + contact_email + "?subject=YLoc%20Webservice>Mail to YLoc Admin</a></p>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;<a href=" + imprint_url + " style='font-size:8pt;color:blue' target='_blank'>Imprint (Impressum)</a></p>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;<a href=" + gdpr_url + " style='font-size:8pt;color:blue' target='_blank'>GDPR Declaration (Datenschutzerklaerung)</a></p>";
  print "</TH></TR>"
  print "</TABLE>";
  print "</BODY></HTML>";

def __print_start_screen(error_msg = ""):
  __print_header();
  __print_head();

  print "<p>YLoc is a interpretable prediction system for protein subcellular localization prediction. In addition to the predicted location, YLoc gives a reasoning why this prediction was made and which biological properties of the protein sequence lead to this prediction. Moreover, a confidence estimate helps users to rate predictions as trustworthy. YLoc+ is able to predict the location of multiple-targeted proteins with high accuracy.</p>";
  print "<BR><BR>";

  m = Models();
  model_combinations = m.getAvailableModelCombinations();

  print "<form name=yloc enctype='multipart/form-data' action='webloc.cgi' method='post'>";

  if error_msg != "":
    print "<p style='color:#FF0000;'><b>An error occured: "+error_msg+"</b></p>";

  #example_Q805D4 = ">sp|Q805D4|CNP3 (Animals)\\nMSLNLPGYALFFILLVASSGAKPAPDLQILEPPLSSLEEQEEMQEEVQEKVQEQQEEVQE\\nKVQEQQEEVQEQQEEVQEQQEEQQEEVQERGRGTGDVLLRAQLDSSTWALQKDDVLMRLF\\nKDLLRTSKRSRSRYKKGGLRSCFGVRLARIGSFSGLGC";

  example_Q75WG7 = ">sp|Q75WG7|U13-HTXT (Animals)\\nMKLSALVFVASVMLVAASPVKDVEEPVETHLAADLKTIEELAKYEEAAVQKRSCIVGSKN\\nIGETCVASCQCCGATVRCIGEGTKGICNNYQTNNILGQILLYAKDTVVNTAGLLVCAQDL\\nSEYE";

  example_P00724 = ">sp|P00724|SUC2 (Fungi)\\nMLLQAFLFLLAGFAAKISASMTNETSDRPLVHFTPNKGWMNDPNGLWYDEKDAKWHLYFQ\\nYNPNDTVWGTPLFWGHATSDDLTNWEDQPIAIAPKRNDSGAFSGSMVVDYNNTSGFFNDT\\nIDPRQRCVAIWTYNTPESEEQYISYSLDGGYTFTEYQKNPVLAANSTQFRDPKVFWYEPS\\nQKWIMTAAKSQDYKIEIYSSDDLKSWKLESAFANEGFLGYQYECPGLIEVPTEQDPSKSY\\nWVMFISINPGAPAGGSFNQYFVGSFNGTHFEAFDNQSRVVDFGKDYYALQTFFNTDPTYG\\nSALGIAWASNWEYSAFVPTNPWRSSMSLVRKFSLNTEYQANPETELINLKAEPILNISNA\\nGPWSRFATNTTLTKANSYNVDLSNSTGTLEFELVYAVNTTQTISKSVFADLSLWFKGLED\\nPEEYLRMGFEVSASSFFLDRGNSKVKFVKENPYFTNRMSVNNQPFKSENDLSYYKVYGLL\\nDQNILELYFNDGDVVSTNTYFMTTGNALGSVNMTTGVDNLFYIDKFQVREVK";

  example_P21549 = ">sp|P21549|AGT1 (Animals)\\nMFQALAKASAAPGSRAAGWVRTMASHKLLVTPPKALLKPLSIPNQLLLGPGPSNLPPRIMAAGGLQMIGSMSKDMYQIMDEI\\nKEGIQYVFQTRNPLTLVISGSGHCALEAALVNVLEPGDSFLVGANGIWGQRAVDIGERIG\\nARVHPMTKDPGGHYTLQEVEEGLAQHKPVLLFLTHGESSTGVLQPLDGFGELCHRYKCLL\\nLVDSVASLGGTPLYMDRQGIDILYSGSQKALNAPPGTSLISFSDKAKKKMYSRKTKPFSF\\nYLDIKWLANFWGCDDQPRMYHHTIPVISLYSLRESLALIAEQGLENSWRQHREAAAYLHG\\nRLQALGLQLFVKDPALRLPTVTTVAVPAGYDWRDIVSYVIDHFDIEIMGGLGPSTGKVLR\\nIGLLGCNATRENVDRVTEALRAALQHCPKKKL";

  example_P41921 = ">sp|P41921|GLR1 (Fungi)\\nMLSATKQTFRSLQIRTMSTNTKHYDYLVIGGGSGGVASARRAASYGAKTLLVEAKALGGT\\nCVNVGCVPKKVMWYASDLATRVSHANEYGLYQNLPLDKEHLTFNWPEFKQKRDAYVHRLN\\nGIYQKNLEKEKVDVVFGWARFNKDGNVEVQKRDNTTEVYSANHILVATGGKAIFPENIPG\\nFELGTDSDGFFRLEEQPKKVVVVGAGYIGIELAGVFHGLGSETHLVIRGETVLRKFDECI\\nQNTITDHYVKEGINVHKLSKIVKVEKNVETDKLKIHMNDSKSIDDVDELIWTIGRKSHLG\\nMGSENVGIKLNSHDQIIADEYQNTNVPNIYSLGDVVGKVELTPVAIAAGRKLSNRLFGPE\\nKFRNDKLDYENVPSVIFSHPEAGSIGISEKEAIEKYGKENIKVYNSKFTAMYYAMLSEKS\\nPTRYKIVCAGPNEKVVGLHIVGDSSAEILQGFGVAIKMGATKADFDNCVAIHPTSAEELV\\nTMR";

  example_P07954 = ">sp|P07954|Fumerate hydratase (FH) (Animals)\\nMYRALRLLARSRPLVRAPAAALASAPGLGGAAVPSFWPPNAARMASQNSFRIEYDTFGEL\\nKVPNDKYYGAQTVRSTMNFKIGGVTERMPTPVIKAFGILKRAAAEVNQDYGLDPKIANAI\\nMKAADEVAEGKLNDHFPLVVWQTGSGTQTNMNVNEVISNRAIEMLGGELGSKIPVHPNDH\\nVNKSQSSNDTFPTAMHIAAAIEVHEVLLPGLQKLHDALDAKSKEFAQIIKIGRTHTQDAV\\nPLTLGQEFSGYVQQVKYAMTRIKAAMPRIYELAAGGTAVGTGLNTRIGFAEKVAAKVAAL\\nTGLPFVTAPNKFEALAAHDALVELSGAMNTTACSLMKIANDIRFLGSGPRSGLGELILPE\\nNEPGSSIMPGKVNPTQCEAMTMVAAQVMGNHVAVTVGGSNGHFELNVFKPMMIKNVLHSA\\nRLLGDASVSFTENCVVGIQANTERINKLMNESLMLVTALNPHIGYDKAAKIAKTAHKNGS\\nTLKETAIELGYLTAEQFDEWVKPKDMLGPK"

  print "<table style='border:1px solid #cbdbff;width:100%;'><tr><td>"
  print " Please copy your protein sequence(s) in one letter code in the box below. A single sequence can be either raw or in FASTA format, whereas multiple sequences need to be in FASTA format. As an alternative, you can upload a FASTA file.<BR>";
  print "<BR>Click here for example sequences: ";
  print "<a target='_self' style='text-decoration:underline;' href=\"javascript:void(document.forms['yloc'].plain_sequence.value='"+str(example_Q75WG7)+"');\">U13-HTXT(Q75WG7)</a>&nbsp;&nbsp;";
  print "<a style='text-decoration:underline;'  target=_self href=\"javascript:void(document.forms['yloc'].plain_sequence.value='"+str(example_P00724)+"')\">SUC2(P00724)</a>&nbsp;&nbsp;";
  print "<a style='text-decoration:underline;'  target=_self href=\"javascript:void(document.forms['yloc'].plain_sequence.value='"+str(example_P21549)+"')\">AGT1(P21549)</a>&nbsp;&nbsp;";
  print "<a style='text-decoration:underline;'  target=_self href=\"javascript:void(document.forms['yloc'].plain_sequence.value='"+str(example_P41921)+"')\">GLR1(P41921)</a>&nbsp;&nbsp;";
  print "<a style='text-decoration:underline;'  target=_self href=\"javascript:void(document.forms['yloc'].plain_sequence.value='"+str(example_P07954)+"')\">FH(P07954)</a>";

  print "</td><td style='vertical-align:bottom;width:50pt;'><a href='webloc.cgi?page=help#start' onclick=\"openHelpWindow(this.href); return false;\"><img src='"+str(img_path)+"helpBlackexp.png'></a></td></tr></table><BR>";

  print "<textarea class=input name='plain_sequence' cols=100 rows=10></textarea><BR>";
  print "OR upload a FASTA file <input type='file' name='fastafile' size='60'><BR><BR>"
  print "Select prediction model: ";
  print "<select name='model' size=1>";

  for elem in model_combinations[0]:
    print "<option value='"+elem+"'>"+str(elem)+"</option>";

  print "</select>&nbsp;&nbsp;&nbsp;Version:";
  print "<select name='origin' size=1>";

  for elem in model_combinations[1]:
    print "<option value='"+elem+"'>"+str(elem)+"</option>";

  print "</select>&nbsp;&nbsp;Use GO term features:";
  print "<select name='goterms' size=1>";

  for elem in model_combinations[2]:
    print "<option value='"+elem+"'>"+str(elem)+"</option>";

  print "</select><BR>";

  print "<input type='hidden' name='page' value='predict'>";
  print "<input class=button type='submit' name='Submit' value='Predict'>";
  print "</form><BR>";

  print "<form enctype='multipart/form-data' action='webloc.cgi' method='post'>";
  print "<table style='border:1px solid #cbdbff;width:100%;'><tr><td>"
  print "To view the results of a previous prediction enter the query ID.";
  print "</td><td style='vertical-align:bottom;width:50pt;'><a href='webloc.cgi?page=help#view' onclick=\"openHelpWindow(this.href); return false;\"><img src='"+str(img_path)+"helpBlackexp.png'></a></td></tr></table><BR>";
  print "Query ID: <input name=id type=text size=40></input>";
  print "<input type='hidden' name='page' value='predicted'>";
  print "<input  class=button  type='submit' name='Submit' value='View'>";
  print "</form>";

  __print_foot();

def __createSequenceObject(input_field, input_file = ""):
  if input_file != "" and not input_file.filename == "":
    input_field = input_file.file.read();

  input_field = input_field.upper();
  lines = input_field.splitlines();
  aasequences = AASequences();
  alphabet = "ACDERTKSLQWYPGHVNMFI ";

  if len(lines) > 0 and len(lines[0]) > 0 and lines[0][0] == '>':
    name = "";
    sequence = "";
    for line in lines:
      # clean line
      line = line.lstrip(" \t");
      line = line.replace(' ','');
      line = line.replace('-','');
      if len(line) > 0 and line[0] == '>':
        if sequence != "":
          aasequences.append((name, sequence));
        name = line[1:];
        sequence = "";
      else:
        for elem in line:
          if not elem in alphabet:
            return None;
        sequence += line;
    if sequence != "":
      aasequences.append((name, sequence));
  else:
    name = "unknown_sequence";
    sequence = "";
    for line in lines:
      # clean line
      line = line.lstrip(" \t");
      line = line.replace(' ','');
      line = line.replace('-','');
      for elem in line:
          if not elem in alphabet:
            return None;

      sequence += line;
    aasequences.append((name, sequence));

  if aasequences.size() == 0:
    aasequences = None;

  return aasequences;

"""
Intermediate Screen and Start of ML3 via independent thread
"""

def __print_intermediate_screen(nr,id,finished=0):
  __print_header(nr*10, id);
  __print_head();

  print "<h2>Prediction in progress</h2>";
  print "Query ID: "+str(id)+"<BR>";
  print "Number of Sequences: "+str(nr)+"<BR>";

  if finished == 0:
    pred_time = "< 1 minute";

    if nr > 5:
      pred_time = "< 2 minutes";
    elif nr > 10:
      pred_time = "< 3 minutes";
    elif nr > 20:
      pred_time = "< "+str(nr/50+1)+"0 minutes";
    #print "Approximated time for prediction: "+pred_time+"<BR>";
    print "<b>Please wait while YLoc is predicting!</b><BR>"
    s="";
    if nr > 1:
      s="s";
    print "Your request contained "+str(nr)+" protein sequence"+s+". The prediction will take approximately "+pred_time+".<BR><BR>"
  else:
    print "<b>YLoc prediction is finished!</b><BR>";
    print "For "+str(nr-finished)+" of "+str(nr)+" protein sequences YLoc needs to create images.<BR><BR>";

  print "You will be redirected to the prediction summary page as soon as the prediction process is finished.<BR>"
  print "Use the query ID to come back later or simply <b>bookmark this page</b><BR>"

  __print_foot();

  sys.stdout.flush()


def __createID():
  id = str(random.random());
  md5object = md5.new(id);
  id = md5object.hexdigest();
  return id;


def __startML(sequences, model, ip, id):
  try:
    sequences.write_fasta_file(config.PATH_TMP+"tmp_sequences_"+str(id)+".fasta");
    os.system("echo " + python_path + " " + config.PATH_MAIN + "yloc.py " + config.PATH_TMP + "tmp_sequences_" + str(id) + ".fasta " + model + " " + str(id) + " >> " + config.PATH_TMP + "ml_out_" + str(id));
    os.system(python_path + " " + config.PATH_MAIN + "yloc.py " + config.PATH_TMP + "tmp_sequences_" + str(id) + ".fasta " + model + " " + str(id) + " >> " + config.PATH_TMP + "ml_out_" + str(id) + " & ");
  except:
    print "\n\n<PRE>"
    traceback.print_exc()


def __getPrediction(sequences, model):
  id = __createID();
  ip = os.environ["REMOTE_ADDR"];
  #ip = str(cgi.escape(os.environ["REMOTE_ADDR"]));
  #ip = gethostbyname(gethostname())
  __startML(sequences,model,ip,id);

  # insert ID to pending table in DB
  alterJobStatus(id,sequences.size(),True,-1,ip);
  # print intermediate screen with refresh
  __print_intermediate_screen(sequences.size(),id);


"""
Decision whether waiting screen is printed or result screen is shown.
"""
def __printResultOrWaitingScreen(id):
  info = getJobInfo(id);

  if len(info) == 1:
    __printResultTable(id);
    os.system("rm "+config.PATH_TMP+"tmp_sequences_"+str(id)+".fasta");
    os.system("rm "+config.PATH_TMP+"ml_out_"+str(id));
  elif len(info) == 3:
    __print_intermediate_screen(info[1],info[0],info[2]);
  else:
    __print_start_screen("Could not find a query with ID "+str(id));


def __id2Prediction(result_id):
  sql_string = "SELECT prediction FROM queries where id='"+str(result_id)+"';";
  res=db_query(sql_string);
  predictions = [];

  for elem in res:
    p = Prediction();
    p.fillObjEctFromStringOfList(elem[0]);
    predictions.append(p);

  return predictions;


##########################################
##### Attribute matrix with +/- ##########
##########################################

def __printAttributeMatrix(p, result_id, sequence_id):
  #name2abbreviation = models.webloc_abbreviations;
  name2abbreviation = webloc_abbreviations;

  help_text = " The table below lists the most important attributes for this particular YLoc prediction beginning with attributes that influenced the prediction at most. A (double) plus indicates that the attribute (strongly) supports the decision for this localization, whereas a (double) minus indicates that the attribute (strongly) supports a decision against this localization. Place the mouse cursor over a field to see the ratio of proteins from this localization having the same attribute value as the query protein.<BR>Click on <i>Attribute Details</i> to see details about the protein concerning this particular attribute.";


  print "<tr><th colspan=2 class=thheading>Attribute Influence<BR>";
  print "<BR><table style='border:1px solid #cbdbff;width:100%;'><tr><td>";
  print help_text;

  #<span><a href='webloc.cgi?page=help#attributes' onclick=\"openHelpWindow(this.href); return false;\"><div style='border:1px solid #cbdbff;'><img height=18 src='"+str(img_path)+"helpBlack18.png'><BR>"+help_text+"<BR>";
  abbreviation_legend = "<BR>Used abbreviations:";
  for i in range(len(p.class_names)):
    if p.class_names[i] in name2abbreviation.keys():
      abbreviation_legend += " " + name2abbreviation[p.class_names[i]] + " = '" + p.class_names[i][:1].upper()+ p.class_names[i][1:] + "',";
  print abbreviation_legend[:-1];
  print "</td><td style='vertical-align:bottom;width:50pt;'><a href='webloc.cgi?page=help#attributes' onclick=\"openHelpWindow(this.href); return false;\"><img src='"+str(img_path)+"helpBlackexp.png'></a></td></tr></table>";
  print "</th></tr>";

  print "<tr><th colspan=3><table class=tablematrix>";
  # get list of attribute order

  help_attribute = "Property and name of attribute";
  help_discr_score = "Discrimination score indicates how well a attributes helps to discriminate against other locations. In detail it is the maximum log-ratio of the probability for the predicted location against another location.";
  help_details = "Click on the <i>Attribute Details</i> button to see details about this attribute.";


  ordered_attributes = p.getSortedAttributeList();
  print "<tr><th colspan=2 class=amatrixside style='vertical-align:top;padding:3pt;'><a class='helptextsmall' href=#G>Attribute<span>"+help_attribute+"</span></a></th>";
  print "<td class=amatrixside><a class='helptext' href='webloc.cgi?page=help#discr' onclick=\"openHelpWindow(this.href); return false;\">Discrimination Score<span>"+help_discr_score+"</span></a></td>";
  for i in range(len(p.class_names)):
    class_name = p.class_names[i];
    class_abv = class_name
    if class_name in name2abbreviation.keys():
      class_abv = name2abbreviation[class_name];
    if i in p.locations:
      print "<td class=amatrixside style='width:20pt;background-color:#FFDDDD; font-weight:bold;'><a class='helptextsmall' href=#G>" + class_abv + "<span>" + class_name+"</span></a></td>";
    else:
      print "<td class=amatrixside style='width:20pt;'><a class='helptextsmall' href=#G>" + class_abv + "<span>"+class_name+"</span></a></td>";
  print "<td class=amatrixside><a class='helptext' href=#G>Detailed Attribute Information<span>"+help_details+"</span></a></td></tr>";
  # sum all discrimations scores
  all_e = 0;
  for elem in p.discrimination_scores:
    all_e += abs(elem);
  # print attributes in order
  summed_e = 0;
  crossed = False;

  for j in range(p.attr_nr):
    i = ordered_attributes[j];
    e = abs(p.discrimination_scores[ordered_attributes[j]]);
    summed_e += e;
    # set row size
    fs = 13 - int(all_e/e)*0.1;
    pd = 6 - int(all_e/e)*0.1;
    if summed_e - e > 0.8* all_e :
      if crossed == False:
        print "<tr><th colspan=4 style='height:8pt; padding:0pt;'></th></tr>";
      crossed = True;
      fs = 7;
      pd = 0;
    style_status = "font-size:"+str(fs)+"pt; padding:"+str(pd)+"pt; padding-left:1pt;";

    # print style for even and even rows
    css_class='amatrix2';
    if j % 2 == 0:
      css_class='amatrix1';


    print "<tr><td class='"+css_class+"' style='"+style_status+"'>"+str(p.feature_description[i][-2])+"</td>";


    print "<td class='"+css_class+"' style='"+style_status+"'>";
    #check for detailed description that can be shown in popup
    #if "<a href" in str(p.feature_description[i][-1]):
    #  descr = str(p.feature_description[i][-1]).strip();
    #  descr = string.replace(descr,"<a href","<a style='"+style_status+"text-decoration:underline;' href");
    #  print "<a class='helptext' style='"+style_status+"' href=#G>"+ descr;
    #else:
    print "<a class='helptext' style='"+style_status+"' href=#G>"+ __extract_links(p.feature_description[i][-1]).strip();
    print "<span><i>"+ str(p.feature_description[i][0]) +"</i> " + __extract_links(p.feature_description[i][1]).strip() +  "</span>";
    print "</a></td>";
    # print discrimination score
    print "<td class='"+css_class+"' style='"+style_status+"'>";
    print "%.2f" % p.discrimination_scores[ordered_attributes[j]] +"</td>";

    # print quantitative supports
    for k in range(len(p.class_names)):
      if k in p.locations:
        print "<td class='amatrix1' style='background-color:#FFDDDD; "+style_status+" '>";
      else:
        print "<td class='"+css_class+"' style='"+style_status+"'>";
      print "<a class='helptext' href=#G style='"+style_status+"'>" + str(p.probability_quality[i][k]) +  "<span>"+("%.0f" % (p.probability_quantity[i][k]*100))+"% of the proteins from the "+str(p.class_names[k])+" have the same property as the query sequence.</span></a></td>";
    # print more infos
    #attribute_info = __getAttributeInfo(p, i);
    print "<td class='"+css_class+"' style='"+style_status+"'><a class='helptextsmall'  target=new href=webloc.cgi?id="+str(result_id)+"&detailed="+str(sequence_id)+"&attribute="+str(i)+" onclick=\"openNewWindow(this.href); return false;\"><img  src='"+str(img_path)+"helpBlackdetails.png'><span>Click the help button to see attribute details.</span></a></td>";

    print "</tr>";

  print "</table></th></tr>";


def __getGOInterval(model):
  m = Models();
  l = m.getAvailableModels();
  for i in range(len(l)):
    l[i] = m.getModelName(l[i]);
  pos = l.index(model);
  if pos != -1:
    return m.getGOIntervals(model);
  return ['?','?','?','?','?','?','?','?','?','?','?','?','?'];

###################################################
##### Very detailed feature description page ######
###################################################

def __getAttributeInfo(p, nr):
  #get attribute distribution
  m = Models();

  distr = m.getAttributeDistributions(p.model,nr);


  #set level of discrimination
  e=p.discrimination_scores[nr];
  if e > 0:
    qualitative = "+";
    if len(p.neg_classes[nr]) > 0:
      qualitative += "/-";
  if e > 2:
    qualitative = "++";
    if len(p.neg_classes[nr]) > 0:
      qualitative += "/--";
  if e < 0:
    qualitative = "-";
    if len(p.pos_classes[nr]) > 0:
      qualitative += "/+";
  if e < -2:
    qualitative = "--";
    if len(p.pos_classes[nr]) > 0:
      qualitative += "/++";
  s = "<table class=tablematrix>";
  true_class_string = p.getClassNameString(p.merged_true_classes[nr]);
  s += "<tr><td class=headerprediction>Attribute</td><td class=headerprediction>Support for "+true_class_string+"</td></tr>";
  s += "<tr><td class=amatrixside style='height:35pt;'><i>" + str(p.feature_description[nr][-2]) + "</i> "+ str(p.feature_description[nr][-1]) + "</td><td class=amatrixside style='font-size:9pt;'><b style='font-size:13pt;'>"+ str(qualitative) + "</b> (discrimination score: "+("%.2f" % (e))+")</td></tr>";
  s += "<tr><th colspan=2 class=amatrix1>";
  if len(p.feature_description[nr])>5:
    s += "<b>Alternative attribute name:</b> <i>"+str(p.feature_description[nr][2])+"</i> "+str(p.feature_description[nr][3])+"<BR>";
  s += "<b>Detailed attribute:</b> <i>"+str(p.feature_description[nr][0]).strip()+"</i> " + str(p.feature_description[nr][1]).strip()+"<BR>";
  s += "</th></tr>";


  if len(p.pos_classes[nr]) > 0:
    s += "<tr><th colspan=2 class=amatrix2>Typical for "+ true_class_string + ( " (%.0f" % (p.true_class_attribute_distribution[nr][p.matched_intervals[nr]]*100) ) +"%) compared to "+str(p.getClassNameString(p.pos_classes[nr]))+ ( " (%.0f" % (p.pos_class_attribute_distribution[nr][p.matched_intervals[nr]]*100) )+"%)</th></tr>";
  if len(p.neg_classes[nr]) > 0:
    s += "<tr><th colspan=2 class=amatrix2>Untypical for "+ true_class_string +  ( " (%.0f" % (p.true_class_attribute_distribution[nr][p.matched_intervals[nr]]*100) ) +"%) compared to "+str(p.getClassNameString(p.neg_classes[nr]))+ ( " (%.0f" % (p.neg_class_attribute_distribution[nr][p.matched_intervals[nr]]*100) )+"%)</th></tr>";


  ## help box text
  help_text = "The plot above displays the distribution of proteins for this attribute.  The x-axis shows the different discretization intervals for this attribute.  In this case, the feature was discretized into "+str(len(p.discretization_intervals[nr]))+" intervals.  The height of the bars indicates the ratio of proteins from a specific localization having an attribute value in this interval. The attribute value of the query protein is located in the highlighted interval, interval "+str(p.matched_intervals[nr]+1)+".";
  if len(p.pos_classes[nr]) > 0:
    help_text += "We observe that " + ("%.0f" % (p.true_class_attribute_distribution[nr][p.matched_intervals[nr]]*100) ) +"% of the proteins from "+ true_class_string +" have a similar attribute value, whereas only "+( "%.0f" % (p.pos_class_attribute_distribution[nr][p.matched_intervals[nr]]*100) )+"% of the proteins form "+str(p.getClassNameString(p.pos_classes[nr]))+" have an attribute value in this interval. Hence, this attribute discriminates well against the localization "+str(p.getClassNameString(p.pos_classes[nr]))+".";
  if len(p.neg_classes[nr]) > 0:
    help_text += "We observe that only " + ("%.0f" % (p.true_class_attribute_distribution[nr][p.matched_intervals[nr]]*100) ) +"% of the proteins from "+ true_class_string +" have a similar attribute value, whereas "+( "%.0f" % (p.neg_class_attribute_distribution[nr][p.matched_intervals[nr]]*100) )+"% of the proteins form "+str(p.getClassNameString(p.neg_classes[nr]))+" have an attribute value in this interval. Hence, this attribute discriminates well against the localization "+str(p.getClassNameString(p.neg_classes[nr]))+".";
  s += "<tr><th colspan=2 class=amatrix3>";

  #include java applet
  s += "<APPLET CODEBASE='"+str(img_path)+"' CODE='Histogram.class' WIDTH='620' HEIGHT='280' ALIGN='BOTTOM' HSPACE='40' VSPACE='10'>\n";
  s += "<param name='classes' value='"+str(len(p.class_names))+"'>\n";
  s += "<param name='intervals' value='"+str(len(distr[0]))+"'>\n";
  for i in range(len(distr)):
    s += "<param name='values"+str(i+1)+"' value='";
    for k in range(len(distr[i])):
      s += str(distr[i][k])+" ";
    s += "'>\n";
  s += "<param name='intervalnames' value='";
  values = [];
  if p.feature_names[nr] == "go_one":
    values = __getGOInterval(p.model);
  else:
    values = p.discretization_intervals[nr]
  for i in range(len(values)):
    s += values[i] + "&";
  s += "'>\n";
  s += "<param name='classnames' value='";
  for i in range(len(p.class_names)):
    #s += models.webloc_abbreviations[p.class_names[i]] + "&";
    s += webloc_abbreviations[p.class_names[i]] + "&";
  s += "'>\n";
  s += "<param name='active' value='";
  for i in range(len(p.class_names)):
    if (i in p.merged_true_classes[nr]) or (i in p.pos_classes[nr]) or (i in p.neg_classes[nr]):
      s += "1 ";
    else:
      s += "0 ";
  s += "'>\n";
  s += "<param name='interval' value='"+str(p.matched_intervals[nr]+1)+"'>\n";
  s += "<param name='predicted' value='";
  for i in range(len(p.class_names)):
    if i in p.merged_true_classes[nr]:
      s += "1 ";
    else:
      s += "0 ";
  s += "'>\n";
  s += "<param name='positive' value='";
  for i in range(len(p.class_names)):
    if i in p.pos_classes[nr]:
      s += "1 ";
    else:
      s += "0 ";
  s += "'>\n";
  s += "<param name='negative' value='";
  for i in range(len(p.class_names)):
    if i in p.neg_classes[nr]:
      s += "1 ";
    else:
      s += "0 ";
  s += "'>\n";
  s += "You need a Java-enabled browser to view this."
  s += "</APPLET>\n";

  #s += "<img border=0 src='"+config.WEBPATH_PLOTS+"plot_"+str(p.id)+"_seq_"+str(p.sequence_id)+"_attr_"+str(nr)+".jpeg'>";
  s += "</th></tr>";
  s += "<tr><th colspan=2 class=amatrix2>"+help_text+"</th></tr>";
  s += "</table>";
  return s;

### return probability table as string
def __getProbabilityTable(p):
  s = "<table>";
  s += "<td class=amatrixside>Location</td><td class=amatrixside>Probability of Location</td></tr>";
  loc_list = p.getSortedLocationList();
  locs = len(p.locations);
  for i in range(locs):
    cname = p.class_names[p.locations[i]]
    s += "<tr><td class=amatrix1 style='height:25pt; font-size:10pt; background-color:#FFDDDD;'>"+cname[:1].upper()+ cname[1:]+"</td><td class=amatrix1 style='height:25pt; font-size:10pt; background-color:#FFDDDD;text-align:center;'> %.1f " % (p.probability_distribution[p.locations[i]]*100);
    s += "%</td></tr>";

  for i in range(locs,len(loc_list)):
    elem = loc_list[i];
    cname = p.class_names[elem];
    if i % 2 == 0:
      s += "<tr><td class=amatrix1>"+cname[:1].upper()+ cname[1:]+"</td><td class=amatrix1 style='text-align:center;'> %.1f " % (p.probability_distribution[elem]*100);
    else:
      s += "<tr><td class=amatrix2>"+cname[:1].upper()+ cname[1:]+"</td><td class=amatrix2 style='text-align:center;'> %.1f " % (p.probability_distribution[elem]*100);
    s += "%</td></tr>";

  s += "</table>";

  return s;

### help function to extract links
def __extract_links(s):
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


#########################################
######## Print short reasoning ##########
#########################################

def __print_short_reasoning(p,detailed=False):
  return p.getShortReasoning();
  """
  output_string = "";

  #calc prob
  prob = 0;
  for i in range(len(p.locations)):
    prob += p.probability_distribution[p.locations[i]];

  ## calc trust
  qualitative_word = "";
  if p.null_prob >= 0.95:
    qualitative_word = "very strong"
  elif p.null_prob >= 0.8:
    qualitative_word = "strong";
  elif p.null_prob >=0.3:
    qualitative_word = "normal";
  else:
    qualitative_word = "small";

  ## more than one
  plural="";
  if len(p.locations)>1:
    plural = "s";

  output_string += "<p>";
  output_string +=  p.model+" predicted that protein sequence <i>"+p.sequence_name+"</i> is located in the "+p.getClassNameString(p.locations)+" with a ";
  #print "<a class='helptext' href=#G>";
  output_string +=  "probability of %.1f" % (prob*100) +"%. ";
  #print "<span style='position:absolute;'><i>"+ "bla" +"</i></a> </span>";
  output_string +=  "YLoc has a "+str(qualitative_word)+ " confidence (%.2f" % (p.null_prob) +") that this prediction is reliable.<BR>";
  output_string +=  "</p>";

  ordered_attributes = p.getSortedAttributeList();
  most_important_attr = -1;
  second_important_attr = -1;
  i=0;
  while i<p.attr_nr and second_important_attr==-1:
    if p.discrimination_scores[ordered_attributes[i]]>0 and len(p.pos_classes[ordered_attributes[i]]) > 0:
      if most_important_attr == -1:
        most_important_attr = ordered_attributes[i];
      else:
        second_important_attr = ordered_attributes[i];
    i += 1;


  output_string +=  "<p class=underl>";
  output_string += "The most important reason for making this prediction is the "+str(p.feature_description[most_important_attr][-2])+" ";

  if detailed:
    output_string += str(p.feature_description[most_important_attr][-1])+".";
  else:
    output_string += str(__extract_links(p.feature_description[most_important_attr][-1]))+".";

  output_string +=   (" %.0f" % (p.true_class_attribute_distribution[most_important_attr][p.matched_intervals[most_important_attr ]]*100) ) +"% of the proteins from the "+ p.getClassNameString(p.locations) +" have a similar attribute, whereas only about "+( "%.0f" % (p.pos_class_attribute_distribution[most_important_attr][p.matched_intervals[most_important_attr]]*100) )+"% of the proteins form the "+str(p.getClassNameString(p.pos_classes[most_important_attr]))+" show this property. ";

  output_string +=  "Moreover, the protein has a "+str(p.feature_description[second_important_attr][-2])+" ";
  if detailed:
    output_string += str(p.feature_description[second_important_attr][-1])+".";
  else:
    output_string += str(__extract_links(p.feature_description[second_important_attr][-1]))+".";
  output_string +=   (" %.0f" % (p.true_class_attribute_distribution[second_important_attr][p.matched_intervals[second_important_attr ]]*100) ) +"% of the proteins from the "+ p.getClassNameString(p.locations) +" have a similar attribute, whereas only about "+( "%.0f" % (p.pos_class_attribute_distribution[second_important_attr][p.matched_intervals[second_important_attr]]*100) )+"% of the proteins form the "+str(p.getClassNameString(p.pos_classes[second_important_attr]))+" show this property. There are more properties that support the predicted location"+plural+".";
  output_string +=  "</p>";

  if detailed:
    output_string +=  "<p class=underl>";
    if len(p.most_similar)>3:
      (id, e_value, identity,go_string) = p.most_similar;
      output_string +=  "The most similar protein from Swiss-Prot 42.0 to the query is <a href='http://www.expasy.org/uniprot/"+str(id)+"' onclick=\"openNewWindow(this.href); return false;\">"+str(id)+"</a> with a local sequence similarity of "+str(identity)+"% (E-value="+str(e_value)+") and is associated with the following GO terms: "+go_string+".";
    else:
      output_string +=  "There exist no very similar protein in Swiss-Prot 42.0.";
    output_string += "</p>";


  return output_string;
  """

## help function to get qualitative word for confidence
def __getQualitativeConfidence(null_prob):

  qualitative_word = "";
  if null_prob >= 0.95:
    qualitative_word = "very strong"
  elif null_prob >= 0.8:
    qualitative_word = "strong";
  elif null_prob >=0.3:
    qualitative_word = "normal";
  else:
    qualitative_word = "small";

  return qualitative_word;


###################################################
##### Detailed prediction description #############
###################################################
def __printDetailedResultTable(predictions,result_id, sequence_id):
  p = predictions[sequence_id];

  print "<BR><BR>";
  print "<table class=tablepredict style='border:1px solid #000000; padding:3px;'>";
  print "<tr><td>";
  print "<h2>Why was this prediction made?</h2>";

  print __print_short_reasoning(p, True);

  print  "</td></tr>";
  print  "</table>";

  print "<BR><p>More details concerning the probability distribution of locations and the protein attributes and their influence on the prediction are shown below.</p>";
  print "<BR><BR>";

  ## print header for details
  print "<h2>Prediction details</h2>";
  print "<table class=tablepredict>";
  print "<form enctype='multipart/form-data' action='webloc.cgi' method='post'>";
  print "<input type='hidden' name='id' value='"+str(predictions[0].id)+"'>";

  ## print probability distribution


  table_text = __getProbabilityTable(p);
  help_text = "The table below displays the probability output of YLoc. The subcellular locations are ordered by their probability beginning with the most probable one. The most probable location or location combination is highlighted in red.";

  print "<tr><th colspan=2 class=thheading><BR>Probability Distribution of Subcellular Locations<BR>";

  print "<BR><table style='border:1px solid #cbdbff;width:100%;'><tr><td>";
  print help_text;
  print "</td><td style='vertical-align:bottom;width:50pt;'><a href='webloc.cgi?page=help#distribution' onclick=\"openHelpWindow(this.href); return false;\"><img src='"+str(img_path)+"helpBlackexp.png'></a></td></tr></table>";
  print "</th></tr>";

  #<span><a href='webloc.cgi?page=help#distribution' onclick=\"openHelpWindow(this.href); return false;\"><div style='border:1px solid #cbdbff;'><img height=18 src='"+str(img_path)+"helpBlack18.png'><BR>"+help_text+"<BR></div></a></span></th></tr>";

  print "<tr><th colspan=2 class=thheading><BR>"+table_text+"</th></tr>";

  ## calc trust
  qualitative_word = __getQualitativeConfidence(p.null_prob);

  ## print trust
  print "<tr><th colspan=2>YLoc has a "+str(qualitative_word)+ " confidence (%.2f" % (p.null_prob) +") that this prediction is reliable.</th></tr>";

  ## print most similar

  if not "*" in p.model:
    print "<tr><th colspan=2 class=thheading><BR>Most similar protein<BR>";

    help_text="The most similar protein is the protein from Swiss-Prot 42.0 with the highest local sequence identity to the query protein. Among others, the GO terms associated with these protein were used for prediction. However, since it is the most similar protein we list it here again.";

    print "<BR><table style='border:1px solid #cbdbff;width:100%;'><tr><td>";
    print help_text;
    print "</td><td style='vertical-align:bottom;width:50pt;'><a href='webloc.cgi?page=help#mostsimilar' onclick=\"openHelpWindow(this.href); return false;\"><img src='"+str(img_path)+"helpBlackexp.png'></a></td></tr></table>";
    print "</th></tr>";

    if len(p.most_similar)>3:
      (id, e_value, identity,go_string) = p.most_similar;
      print "<tr><th colspan=2 class=underl>The most similar protein from Swiss-Prot 42.0 to the query is <a href='http://www.expasy.org/uniprot/"+str(id)+"' onclick=\"openNewWindow(this.href); return false;\">"+str(id)+"</a> with a local sequence similarity of "+str(identity)+"% (E-value="+str(e_value)+") and is associated with the following GO terms: "+go_string+".</th></tr>";
    else:
      print "<tr><th colspan=2 class=underl>There exist no very similar protein in Swiss-Prot 42.0.</th></tr>";
    print "<BR>";

  __printAttributeMatrix(p,result_id,sequence_id);

def __printAttributeTable(predictions,result_id, sequence_id, attribute_id):
  p = predictions[sequence_id];

  __print_header();

  print "<BODY bgcolor=#FFFFFF>";
  print "<table style='width:100%;'><tr><td>";
  print "<center>";
  print "<h3>YLoc Prediction Result for sequence <i>"+predictions[sequence_id].sequence_name+"</i><BR>-Attribute Details-</h3>";
  print "</td><td style='vertical-align:bottom;width:50pt;'><input class=button style='font-size:8pt;' type='button' onClick='javascript:window.close();' value='close window'></td></tr>";
  print "</table>";

  print "<BR><table style='border:1px solid #cbdbff;width:100%;'><tr><td>";
  print "The table below summarizes an attribute. The header shows the name and the discrimination score of this attribute. The first row describes the property of the sequence in more detail. In addition, it is shown whether the attribute value is typical for the predicted location(s). Subcellular locations that have the most differently properties concerning this attribute are named and displayed in the barplot below.";
  print "</td><td style='vertical-align:bottom;width:50pt;'><a href='webloc.cgi?page=help#details' onclick=\"openHelpWindow(this.href); return false;\"><img src='"+str(img_path)+"helpBlackexp.png'></a></td></tr></table><BR>";


  #print "<p><span><a href='webloc.cgi?page=help#details' onclick=\"openHelpWindow(this.href); return false;\" style='text-decoration:none;'><div style='border:1px solid #cbdbff;border:bottom:0pt;'><img height=18 src='"+str(img_path)+"helpBlack18.png'><BR><BR></div></a></span></p>";
  print "<center>";

  attribute_info = __getAttributeInfo(p, attribute_id);
  print attribute_info;

  print "<BR><table style='width:100%;'><tr><td>";
  print "<center>";
  print "</td><td style='vertical-align:bottom;width:50pt;'><input class=button style='font-size:8pt;' type='button' onClick='javascript:window.close();' value='close window'></td></tr>";
  print "</table>";


#######################################
##### Prediction overview table #######
#######################################

def __printOverviewTable(predictions, result_id):
  print "<table class=tablepredict>";
  print "<form enctype='multipart/form-data' action='webloc.cgi' method='post'>";
  print "<input type='hidden' name='id' value='"+str(predictions[0].id)+"'>";
  print "<tr><th colspan=2 class=thheading>";
  help_text = "The prediction of YLoc is displayed in the table below. Place the mouse cursor over the table to see more details. Click on the <i>Elucidate</i> button to obtain an explanation why this prediction has been made.";
  print "<BR><table style='border:1px solid #cbdbff;width:100%;'><tr><td>";
  print help_text;
  print "</td><td style='vertical-align:bottom;width:50pt;'><a href='webloc.cgi?page=help#overview' onclick=\"openHelpWindow(this.href); return false;\"><img src='"+str(img_path)+"helpBlackexp.png'></a></td></tr></table>";
  print "</th></tr>";

  print "<tr><th colspan=2><table>";
  help_query_sequence = "Name of the sequence extracted from the given fasta format.";
  help_predicted_location = "The most probable subcellular localization(s) predicted by YLoc.";
  help_probability = "The probability of the predicted subcellular localizations(s). ";
  help_confidence = "The confidence score that the predicted subcellular localizations(s) are correct. Proteins which are typical for YLoc can be predicted with a higher reliability and therefore are assigned with a higher confidence score.";
  help_details = "Click on the <i>Elucidate</i> button to see a more detailed prediction including a reasoning. ";

  print "<tr><td class=amatrixside><a class='helptext' href=#G>Query Sequence<span>"+help_query_sequence+"</span></a></td><td class=amatrixside><a class='helptext' href=#G>Predicted Location<span>"+help_predicted_location+"</span></a></td><td class=amatrixside><a class='helptext' href=#G>Probability<span>"+help_probability+"</span></a></td><td class=amatrixside><a class='helptext' href=#G>Confidence<span>"+help_confidence+"</span></a></td><td class=amatrixside><a class='helptext' href=#G>Detailed Info<span>"+help_details+"</span></a></td><td width=0></td></tr>";

  even = False;
  for j in range(len(predictions)):
    p = predictions[j]
    table_text = "<h3>Probability Distribution for "+str(p.sequence_name)+"</h3>"+ __getProbabilityTable(p);
    table_text = "<h3>Why?</h3>"+ __print_short_reasoning(p, False)+"<p><b>For more details click the <i>Elucidate</i> button.</b></p>";

    # change view from line to line
    even = not even;
    viewtype="amatrix1";

    # print query sequence name
    print "<tr><td class='"+viewtype+" probtable'>"
    print "<a class='"+viewtype+" probtable' href=webloc.cgi?id="+str(result_id)+"&detailed="+str(j)+">"
    print str(p.sequence_name)
    print "</a></td>";

    # print predicted location
    true_class_string = p.getClassNameString(p.locations, True);
    putative_class_string = p.getClassNameString(p.putative_locations);
    if len(putative_class_string) > 2:
      print "<td class='"+viewtype+" probtable'>";
      print "<a class='"+viewtype+" probtable' href=webloc.cgi?id="+str(result_id)+"&detailed="+str(j)+">";
      print true_class_string+"</b> (also possible: "+putative_class_string+")";
      print "</a></td>";
    else:
      print "<td class='"+viewtype+" probtable'>";
      print "<a class='"+viewtype+" probtable' href=webloc.cgi?id="+str(result_id)+"&detailed="+str(j)+">";
      print true_class_string;
      print "</a></td>";

    # print confidence
    conf = 0;
    for i in range(len(p.locations)):
      conf += p.probability_distribution[p.locations[i]];
    conf = conf/len(p.locations);
    print "<td class='"+viewtype+" probtable'>";
    print "<a class='"+viewtype+" probtable' href=webloc.cgi?id="+str(result_id)+"&detailed="+str(j)+">"
    print "%.2f " % (conf*100) + "% ";
    print "</a></td>";

    # print confidence
    print "<td class='"+viewtype+" probtable'>";
    print "<a class='"+viewtype+" probtable' href=webloc.cgi?id="+str(result_id)+"&detailed="+str(j)+">"
    print __getQualitativeConfidence(p.null_prob)+" (%.2f) " % (p.null_prob);
    print "</a></td>";

    # print link to detailed description
    print "<td class="+viewtype+" probtable'><a class='probtable' href=webloc.cgi?id="+str(result_id)+"&detailed="+str(j)+"><img  src='"+str(img_path)+"helpBlackelucidate.png'>"
    #<span>"+table_text+"</a></span></td>";
    print "</a></td>";

    #print description WHY
    print "</tr><tr><th colspan=5 class=amatrix3>"
    print table_text
    print "</th></tr>";

    print "<tr><td></td></tr>";


  print "</table>";
  nr_brs=10-int(len(predictions)*1.3);
  for i in range(nr_brs):
    print "<BR>"

###################################################
##### General results screen which redirects ######
##### to detailed view or overview           ######
###################################################

def __printResultTable(result_id):

  #get data from DB
  predictions = __id2Prediction(result_id);

  # check display mode
  detailed = False;
  attribute = False;
  sequence_id = 0;
  attribute_id = 0;
  if form.has_key("detailed"):
    sequence_id = int(form["detailed"].value);
    if not sequence_id in range(len(predictions)):
      __print_error("There exist no sequence with sequence_id "+str(sequence_id)+" for this query!");
    else:
      detailed = True;

  if form.has_key("attribute") and detailed:
    attribute_id = int(form["attribute"].value);
    if not attribute_id in range(predictions[sequence_id].attr_nr):
      __print_error("There exist no attribute "+str(attribute_id) +" for sequence with sequence_id "+str(sequence_id)+" for this query!");
    else:
      attribute = True;

  if attribute:
    __printAttributeTable(predictions, result_id, sequence_id, attribute_id);

  else:
    __print_header();
    __print_head();

    # print heading
    if not detailed:
      print "<h2>Prediction Summary</h2>";

    date = predictions[0].date[6:8] + "/" + predictions[0].date[4:6] + "/" + predictions[0].date[0:4]+ " " + predictions[0].date[8:10] + ":" + predictions[0].date[10:12] + ":" +predictions[0].date[12:14];
    print "<table class=tablepredict><tr><td>";
    print "<p style='padding-left:20pt; font-size:8pt;'>Query ID: " + str(predictions[0].id) + " Query Date: " + str(date) +"<BR>";
    print "Prediction based on model <i> "+str(predictions[0].model)+"</i>.</p>";
    print "</td><td style='text-align:right;'>";
    if detailed:
      print "<a href=webloc.cgi?id="+str(result_id)+"><input class=button type='button' onClick=\"window.location.href='webloc.cgi?id="+str(result_id)+"';\" value='back to prediction summary'></a>";
    print "</td></tr></table>";


    # get real fillings for page
    if not detailed:
      __printOverviewTable(predictions, result_id);

    else:
      __printDetailedResultTable(predictions, result_id, sequence_id);
      print "</table>";
      print "<BR><a href=webloc.cgi?id="+str(result_id)+"><input class=button type='button' onClick=\"window.location.href='webloc.cgi?id="+str(result_id)+"';\" value='back to prediction summary'></a>";



    print "</th></tr>";
    print "</form>";

    __print_foot();


def __print_prediction_page():
  if form.has_key("plain_sequence"):
    if form.has_key("fastafile"):
      aasequences = __createSequenceObject(form["plain_sequence"].value, form["fastafile"]);
    else:
      aasequences = __createSequenceObject(form["plain_sequence"].value);

    if aasequences == None:
      __print_start_screen("Given protein sequence(s) have wrong format. You have been redirected to the start page.");
      return;
    too_small = False;
    for i in range(aasequences.size()):
      (name,seq) = aasequences.get(i);
      if len(seq) < 20:
        too_small = True;
    if too_small:
      __print_start_screen("Given protein sequence too small (less than 20 amino acids). You have been redirected to the start page.");
    elif aasequences.size() > 20:
      __print_start_screen("Please use at most 20 protein sequences in the web service. If you need to predict more proteins, please contact us or use the YLoc SOAP interface.");
    elif not form.has_key("model"):
      __print_start_screen("No model given for prediction. You have been redirected to the start page.");
    elif not form.has_key("origin"):
      __print_start_screen("No origin given for prediction. You have been redirected to the start page.");
    elif not form.has_key("goterms"):
      __print_start_screen("No GO terms selection given for prediction. You have been redirected to the start page.");
    else:
      m = Models();
      model = m.getModelFromCombinations(form["model"].value, form["origin"].value, form["goterms"].value);
      __getPrediction(aasequences, model);
      #__printResultTable(result_id);
  else:
    __print_start_screen("Missing value of parameter 'plain_sequence'. You have been redirected to the start page.");

def __print_error(error_msg):
  print "<p><b><font color=#FF0000>An error occured: "+error_msg+"</b></p>";

def __print_help_page():
  __print_header();
  __print_head();

  print "<h2>YLoc Tutorial</h2>";
  print "<table>";
  print "<tr><td class=thheading style='background-color:#cbdbff;' id='start'>Starting Predictions</td></tr>";
  print "<tr><td>To start a new YLoc prediction, copy your protein sequences in one letter code into the provided textbox. For a single protein sequence this could be either simple one letter code or FASTA format. For multiple sequences FASTA format is required. As an alternative, you may provide your sequences as FASTA file. <BR>Then, select a model for your prediction. The following models are available:";
  print "<center><table width=80%>";
  print "<tr><td width=100pt>YLoc-LowRes</td><td>predicts into 4 locations (nucleus, cytoplasm, mitochodrion, secretory pathway for the animal and fungi version) or 5 locations (in addition chloroplast for the plant version), respectively.</td></tr>";
  print "<tr><td>YLoc-HighRes</td><td>predicts into 9 or 10 locations, respectively. These are nucleus, cytoplasm, mitochodrion, plasma membrane, extracellular space, endoplasmic reticulum, peroxisome, and Golgi apparatus for all models. In addition, lysosome for the animal model, vacuole for the fungi model, and vacuole and chloroplast for the plant model.</td></tr>";
  print "<tr><td>YLoc+</td><td>predicts into 9 or 10 locations, as described above. In addition, it allows to predict multiple locations. It was trained, in addition to the 11 main eukaryotic location classes, on 7 multi-location classes.";
  print "</td></tr>";
  print "</table></center><BR>";
  print "Every model is available in a version specialized on animal, fungal, and plant proteins. Moreover, the use of GO terms transfered from close homologous proteins can be switched off.<BR><BR>";
  print "Click the <i>Predict</i> button to start the prediction. During the prediction process a waiting screen in displayed. If you wish to skip the waiting and instead come back later to view the prediction results, write down the query ID."
  print "<BR><table bgcolor=#F5F5F5 width=100%><tr><td>";
  print "<h3>Example Input:</h3>";
  print "<p>>sp|Q75WG7|U13-HTXT (Animals)<BR>";
  print "MKLSALVFVASVMLVAASPVKDVEEPVETHLAADLKTIEELAKYEEAAVQKRSCIVGSKN<BR>";
  print "IGETCVASCQCCGATVRCIGEGTKGICNNYQTNNILGQILLYAKDTVVNTAGLLVCAQDL<BR>"
  print "SEYE</p>";
  print "<p><i>To follow this example, copy the protein sequence into the text box on the <a href=http://www-abi.informatik.uni-tuebingen.de/Services/YLoc/webloc.cgi>YLoc startpage</a> and click the 'predict' button.</i></p>"
  print "<h3>Example Output: </h3>"
  print "<img src="+str(img_path)+"example1.png><BR>";
  print "<BR>The corresponding waiting page: <BR>";
  print "<img src="+str(img_path)+"example12.png>";
  print "</td></tr></table>";
  print "<BR><input class=button style='font-size:8pt;' type='button' onClick='javascript:window.close();' value='close window'><BR><BR></td></tr>";
  print "<tr><td class=thheading style='background-color:#cbdbff;' id='view'>View Previous Prediction Results</td></tr>";
  print "<tr><td>Input the query ID from a previous predictions into the query ID field of the start page. By clicking the <i>View</i> button you will be redirected to the result page of this particular prediction. Note, in the case that your prediction process hasn't finish yet you will be redirected to the waiting page.";
  print "<BR><table bgcolor=#F5F5F5 width=100%><tr><td>";
  print "<h3>Example:</h3>";
  print "<p>a typical query ID: <a href=http://www-abi.informatik.uni-tuebingen.de/Services/YLoc/webloc.cgi?id=42847fd6e7a81248451be0cf325b3b43>42847fd6e7a81248451be0cf325b3b43</a></p>";
  print "<img src="+str(img_path)+"example2.png>";
  print "</td></tr></table>";
  print "<BR><input class=button style='font-size:8pt;' type='button' onClick='javascript:window.close();' value='close window'><BR><BR></td></tr>";
  print "<tr><td class=thheading style='background-color:#cbdbff;' id='overview'>Prediction Summary</td></tr>";
  print "<tr><td>The result overview page shows you a short overview over all predictions that have been made for your query. The predictions are displayed in a table, one row for every prediction. For every prediction the name of the protein sequence is displayed, the predicted location, the probability of this location, and the confidence of YLoc that the predicted location is correct. If the sequence was entered in raw sequence format, the name of the sequence will be <i>unknown sequence</i>. Below a short reasoning is displayed in reader friendly format. The short reasining refers to the two most important properties that lead to this particular prediction. By clicking on the 'Elucidate' button or any other cell in this row, you will get additional information.<BR>At the top of the page you see additional information concerning your prediction like the date and time the query was submitted, the query ID, and the model used for predicting."
  print "<BR><table bgcolor=#F5F5F5 width=100%><tr><td>";
  print "<h3>Example:</h3>";
  print "<p>Prediction summary for query <a href=http://www-abi.informatik.uni-tuebingen.de/Services/YLoc/webloc.cgi?id=42847fd6e7a81248451be0cf325b3b43>42847fd6e7a81248451be0cf325b3b43</a>:</p>";
  print "<img src="+str(img_path)+"example3.png>";
  print "</td></tr></table>";
  print "<BR><input class=button style='font-size:8pt;' type='button' onClick='javascript:window.close();' value='close window'><BR><BR></td></tr>";

  print "<tr><td class=thheading style='background-color:#cbdbff;' id='distribution'>Probability Distribution and Confidence Score</td></tr>";
  print "<tr><td>The YLoc webservice returns for a query protein a probability distribution of subcellular locations. The probability for every location is an estimate of YLoc how likely this protein is present in this particular location. The locations are ordered in a table according to their probability, beginning with the most probable one. The predicted location or location combination (for YLoc+) is highlighted with red background. For YLoc-LowRes and YLoc-HighRes this is always the most probable location. YLoc+ is able to predict multiple locations and, thus, highlights all predicted locations. <BR> In addition, YLoc returns a confidence score that lies between 0 and 1. The larger the confidence score, the higher the confidence that this prediction is correct."
  print "<BR><table bgcolor=#F5F5F5 width=100%><tr><td>";
  print "<h3>Example:</h3>";
  print "<p>Probability Distribution of Locations for protein <i>Q75WG7</i> from query <a href=http://www-abi.informatik.uni-tuebingen.de/Services/YLoc/webloc.cgi?id=42847fd6e7a81248451be0cf325b3b43&detailed=0>42847fd6e7a81248451be0cf325b3b43</a>:</p>";
  print "<img src="+str(img_path)+"example4.png>";
  print "</td></tr></table>";
  print "<BR><input class=button style='font-size:8pt;' type='button' onClick='javascript:window.close();' value='close window'><BR><BR></td></tr>";

  print "<tr><td class=thheading style='background-color:#cbdbff;' id='mostsimilar'>Most Similar Protein</td></tr>";
  print "<tr><td>As an additional information, YLoc displays the most similar protein in the training dataset of YLoc. That is, it will only show proteins from Swiss-Prot 42.0. Moreover, it shows the significance of the BLAST hit and shows which GO terms are associated to this protein according to the YLoc Swiss-Prot to GO mapping."
  print "<BR><table bgcolor=#F5F5F5 width=100%><tr><td>";
  print "<h3>Example:</h3>";
  print "<p>Most similar protein of <i>Q75WG7</i> from query <a href=http://www-abi.informatik.uni-tuebingen.de/Services/YLoc/webloc.cgi?id=42847fd6e7a81248451be0cf325b3b43&detailed=0>42847fd6e7a81248451be0cf325b3b43</a>:</p>";
  print "<img src="+str(img_path)+"example5.png>";
  print "</td></tr></table>";
  print "<BR><input class=button style='font-size:8pt;' type='button' onClick='javascript:window.close();' value='close window'><BR><BR></td></tr>";

  print "<tr><td class=thheading style='background-color:#cbdbff;' id='attributes'>Attribute Influence Table</td></tr>";
  print "<tr><td>The attribute influence table is a summary which shows how attributes influence the final prediction. Every attribute is displayed in one row of the table. The attributes are ordered according to their influence to the prediction, beginning with most influencing attributes. To visualize this sorting more clearly, the rows have different height and font size. An extra gap is included to separate the attributes that have 80&#37; of the influence on the prediction from the other not so important attributes. In every row the attribute name and value, the discrimination score, and the influence on every subcellular location is displayed.  The discrimination score measures how strongly the attribute influences the prediction. A large positive value indicates a strong support for the predicted location. That is, the observed attribute is more likely for proteins from the predicted location than for proteins of other locations. In contrast, a negative discrimination score opposes the predicted location. That is, it is more likely to observe the given attribute in proteins from a location other than the predicted location. By moving the mouse cursor of the attribute name, the user can gain additional information about this attribute.The influence of the attribute is displayed with a (double) plus if it (stongly) supports this location and displayed with a (double) minus if it (strongly) opposes this location, respectively. When the mouse cursor is moved over the table cell, the percentage of proteins from this location having the same property as the query protein is displayed. A detailed description, graphical view, and some more details about its discrimination ability is available for every protein by clicking on the <i>Attribute Details</i> button in the very right column."
  print "<BR><table bgcolor=#F5F5F5 width=100%><tr><td>";
  print "<h3>Example:</h3>";
  print "<p>Attribute influence table for protein <i>Q75WG7</i> from query <a href=http://www-abi.informatik.uni-tuebingen.de/Services/YLoc/webloc.cgi?id=42847fd6e7a81248451be0cf325b3b43&detailed=0>42847fd6e7a81248451be0cf325b3b43</a>:</p>";
  print "<img src="+str(img_path)+"example6.png>";
  print "</td></tr></table>";
  print "<BR><input class=button style='font-size:8pt;' type='button' onClick='javascript:window.close();' value='close window'><BR><BR></td></tr>";
  print "<tr><td class=thheading style='background-color:#cbdbff;' id='discr'>Discrimination Score</td></tr>";
  print "<tr><td>The discrimination score is a measure how strongly the attribute discriminates between the predicted location and the other locations. A positive discrimination score shows that the attribute supports this location, whereas a negative discrimination score shows that the attribute opposes this location.  The discrimination score is based on the probability of observering this particular attribute value in the different locations.  If the attribute value is more likely to be observed in the predicted location than in another location, the feature discriminates well and supports the predicted location the prediction model. Thus, the attribute gets a high discrimination score. For more details concerning the calculation of the discrimination score see the manuscript or the example below."
  print "<BR><table bgcolor=#F5F5F5 width=100%><tr><td>";
  print "<h3>Example:</h3>";
  print "<p>Discrimination score of a <i>strong secretory pathway sorting signal</i> of protein <i>Q75WG7</i> from query <a href=http://www-abi.informatik.uni-tuebingen.de/Services/YLoc/webloc.cgi?id=42847fd6e7a81248451be0cf325b3b43&detailed=0>42847fd6e7a81248451be0cf325b3b43</a>:</p>";
  print "<img src="+str(img_path)+"example7.png>";
  print "<p>In this example 69&#37; of the secreted proteins contain a similarly strong secretory pathway sorting signal. In contrast, only 0&#37;, 1&#37;, and 2&#37; of the cytoplasmic, mitochondrial, and nuclear proteins have such an attribute, respectively. Hence, the attribute is a good discriminator and discriminates particularly well against cytoplasmic proteins. The support score is calculated by ln(0.69/0.00226) = 5.72. For the two other locations, the same calculation would result in the values 4.2 and 3.5. However, only the highest discrimination score is displayed in this example. In the detailed attribute page (see below), the distributions of the attribute for the different locations is plotted.</p>"
  print "</td></tr></table>";
  print "<BR><input class=button style='font-size:8pt;' type='button' onClick='javascript:window.close();' value='close window'><BR><BR></td></tr>";
  print "<tr><td class=thheading style='background-color:#cbdbff;' id='details'>Attribute Details</td></tr>";
  print "<tr><td>The attribute details page summarizes various informations about an attribute beginning with the name of the attribute together with its discrimination score. In addition, alternative descriptions of the attribute are displayed below. Above the plot, a short sentence summarizes whether the attribute is typical for the predicted location or not and which locations show opposite behaviour. The subcelluar locations that significantly differ in this attribute value displayed in a barplot. The barplot shows the attribute with its discretization intervals on the x-axis. The y-axis shows the ratio of proteins from a locations having an attribute value within a discretization interval. The attribute value of the query protein is located in the highlighted discretization interval. Here, it can be observed how many proteins (ratio) from the predicted class have this particular attribute values. In additon, the user can see how many proteins of the other location(s) have these attribute values. The difference between the bar heights shows whether the attribute value is typical or not for the predicted location(s). The user can include other locations in the barplot and exclude already shown locations. "
  print "<BR><input class=button style='font-size:8pt;' type='button' onClick='javascript:window.close();' value='close window'><BR><BR></td></tr>";
  print "<BR><table bgcolor=#F5F5F5 width=100%><tr><td>";
  print "<h3>Example:</h3>";
  print "<p>Attribute details of a <i>strong secretory pathway sorting signal</i> of protein <i>Q75WG7</i> from query <a href=http://www-abi.informatik.uni-tuebingen.de/Services/YLoc/webloc.cgi?id=42847fd6e7a81248451be0cf325b3b43&detailed=0&attribute=3>42847fd6e7a81248451be0cf325b3b43</a>:</p>";
  print "<img src="+str(img_path)+"example8.png>";
  print "<p>The displayed attribute strongly supports the secretory pathway, indicated by the double plus in the header of the table. The strong secretory pathway signal corresponds to a very hydrophobic N-terminus and was calculated using the autocorrelation of every third hydrophobic residue within the first 20 N-terminal amino acids.  In the plot, we observe that a far more secreted proteins share such a property, whereas nuclear and cytoplasmic proteins only rarely have this property.  </p>"
  print "</td></tr></table>";

  __print_foot();

def __uniqify(seq):
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


def __print_admin_page():
  __print_header();
  __print_head();

  if form.has_key("delete"):
    id_list = form.getlist("delete");
    id_list = __uniqify(id_list);
    print "Start deletion...<BR>";

    for elem in id_list:
      sql_string = "DELETE FROM queries where id='"+str(elem)+"';";
      res=db_query(sql_string);
      sql_string = "DELETE FROM pending where id='"+str(elem)+"';";
      res=db_query(sql_string);
      #os.system("rm "+config.PATH_PLOTS+"plot_"+str(elem)+"* ");
      print "deleted "+str(elem)+"<BR>";

  else:

    sql_string = "SELECT id,date,ip FROM queries";
    res=db_query(sql_string);

    t0 = time.time();
    tm = time.localtime(t0);
    (year, month, day, hour, minute, second,a,b,c) = tm;
    easytime = (int(year)-2008)*365+int(month)*31+int(day);

    print "<h2>YLoc Admin Area</h2>";
    print "<p>Displayed below are all saved predictions from YLoc. Predictions older than 20 days are marked for deletion.</p>";
    print "<table>";
    print "<form enctype='multipart/form-data' action='webloc.cgi' method='post'>";
    print "<input type='hidden' name='page' value='MuLESB82'>";

    print "<tr class=amatrixside><td>Delete option</td><td>query ID</td><td>Prediction date</td><td>IP</td><td>View?</td></tr>";
    c = 1;
    for elem in res:
      (id, date,ip) = elem;
      easytime2 = (int(date.year)-2008)*365+int(date.month)*31+int(date.day);
      value="";
      if (easytime - easytime2 > 20):
        value="checked";

      css = "amatrix1"
      if c % 2:
        css = "amatrix2"
      print "<tr class="+css+"><td><input type='checkbox' name='delete' value='"+str(id)+"' "+value+">delete?</td><td>"+str(id)+"</td><td>"+str(date)+"</td><td>"+str(ip)+"</td><td><a class='helptextsmall' href=webloc.cgi?id="+str(id)+"><img height=20 src='"+str(img_path)+"helpBlack20.png'></a></td></tr>";



      c += 1
    print "</table><BR><BR>"

    print "<input class=button type='submit' name='Submit' value='Delete marked predictions.'>";
    print "</form>";


    sql_string = "SELECT id,date FROM pending";
    res=db_query(sql_string);

    t0 = time.time();
    tm = time.localtime(t0);
    (year, month, day, hour, minute, second,a,b,c) = tm;
    easytime = (int(year)-2008)*365+int(month)*31+int(day);

    print "<p>Displayed below are all pending predictions from YLoc.</p>";
    print "<table>";
    print "<form enctype='multipart/form-data' action='webloc.cgi' method='post'>";
    print "<input type='hidden' name='page' value='MuLESB82'>";

    print "<tr class=amatrixside><td>Delete option</td><td>query ID</td><td>Prediction date</td></tr>";
    c = 1;
    for elem in res:
      (id, date) = elem;
      easytime2 = (int(date.year)-2008)*365+int(date.month)*31+int(date.day);
      value="";
      if (easytime - easytime2 > 20):
        value="checked";

      css = "amatrix1"
      if c % 2:
        css = "amatrix2"
      print "<tr class="+css+"><td><input type='checkbox' name='delete' value='"+str(id)+"' "+value+">delete?</td><td>"+str(id)+"</td><td>"+str(date)+"</td></tr>";
      c += 1
    print "</table><BR><BR>"

    print "<input class=button type='submit' name='Submit' value='Delete marked pending predictions.'>";
    print "</form>";

  __print_foot();


def __blastp_version_info():
  try:
    p = Popen("blastp -version", shell=True, stdout=PIPE, stderr=PIPE);
    blastp_v, stderr  = p.communicate();

    print "<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" + blastp_v.split("\n")[0].strip() + "</p>";
    print "<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" + blastp_v.split("\n")[1].strip() + "</p>";
  except:
    if os.getenv("YLOC_DEBUG"):
      print "\n\n<PRE>"
      traceback.print_exc()
    else:
      print "<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Cannot detect BLASTp version</p>";


def __pfscan_version_info():
  try:
    p = Popen("pfscan", shell=True, stdout=PIPE, stderr=PIPE)
    stdout, pfscan_v = p.communicate()

    print "<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" + pfscan_v.split("\n")[1].strip() + "</p>";
  except:
    if os.getenv("YLOC_DEBUG"):
      print "\n\n<PRE>"
      traceback.print_exc()
    else:
      print "<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Cannot detect pfscan version</p>";


def __print_info_page():
  __print_header();
  __print_head();

  print "<h2>Terms of Service</h2>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;<a href=" + imprint_url + " style='color:blue' target='_blank'>Imprint (Impressum)</a></p>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;<a href=" + gdpr_url + " style='color:blue' target='_blank'>GDPR Declaration (Datenschutzerklaerung)</a></p>";
  print "<BR><BR>";

  print "<h2>Contact</h2>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;Please write an eMail to  <a style='color:blue' href=mailto:" + contact_email + "?subject=YLoc%20Webservice>YLoc Admin</a></p>";
  print "<BR><BR>";

  print "<h2>How to Cite</h2>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;Sebastian Briesemeister, J&ouml;rg Rahnenf&uuml;hrer, and Oliver Kohlbacher, (2010): Going from where to why - interpretable prediction of protein subcellular localization.</p>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;<a href='https://dx.doi.org/10.1093/bioinformatics/btq115' style='color:blue' target='_blank'>Bioinformatics, 26(9):1232-1238.</a></p>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;Sebastian Briesemeister, J&ouml;rg Rahnenf&uuml;hrer, and Oliver Kohlbacher, (2010): YLoc - an interpretable web server for predicting subcellular localization.</p>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;<a href='https://dx.doi.org/10.1093%2Fnar%2Fgkq477' style='color:blue' target='_blank'>Nucleic Acids Research, 38:W497-W502.</a></p>";
  print "<BR><BR>";

  print "<h2>License</h2>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;YLoc is distributed under the GNU General Public License (GPL).</p>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;YLoc is using LIBSVM, BLAST, and InterProScan. These software tools have their own license terms.</p>";
  print "<BR><BR>";

  print "<h2>Third Party Software Tools</h2>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;YLoc is implemented in Python 2.7 and uses BLAST and pfscan.</p>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;BLASTp:</p>";
  __blastp_version_info();
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;pfscan:  </p>";
  __pfscan_version_info();
  print "<BR><BR>";

  print "<h2>Run YLoc Locally</h2>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;YLoc is available as a repository to build a docker image, which also includes this websever: <a href='https://github.com/KohlbacherLab/YLoc' style='color:blue' target='_blank'>YLoc at KohlbacherLab GitHub</a></p>";
  print "<BR>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;Please note:</p>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;Using other versions of third party software tools, especially of BLAST and pfscan, in your local installation may result in slightly different prediction scores compared to those calculated by the online service.</p>";
  print "<BR><BR>";

  print "<h3>Downloads</h3>";
  print "<p>&nbsp;&nbsp;&nbsp;&nbsp;<a href='http://gpcr.biocomp.unibo.it/bacello/dataset.htm' style='color:blue' target='_blank'>BaCelLo datasets</a><BR><BR>";
  print "&nbsp;&nbsp;&nbsp;&nbsp;<a href='" + download_path + "multiloc2_datasets.tar.bz2' style='color:blue' download>H&ouml;glund Datasets (Multiloc / Multiloc2) (BZ2-archive, 1,7 MB)</a><BR><BR>";
  print "&nbsp;&nbsp;&nbsp;&nbsp;<a href='" + download_path + "DBMLocDataset.zip' style='color:blue' download>DBMLoc dataset (ZIP-archive, 883 KB)</a><BR><BR>";
  print "&nbsp;&nbsp;&nbsp;&nbsp;<a href='" + download_path + "YLocFeatureFiles.zip' style='color:blue' download>YLoc feature files (Arff format in ZIP-archive, 6.78 MB)</a><BR><BR>";
  print "&nbsp;&nbsp;&nbsp;&nbsp;<a href='" + download_path + "swissprot2go.zip' style='color:blue' download>SwissProt to GO-term map (ZIP-archive, 4.82 MB)</a><BR><BR>";
  print "&nbsp;&nbsp;&nbsp;&nbsp;<a href='ftp://ftp.expasy.org/databases/swiss-prot/sw_old_releases/' style='color:blue' target='_blank'>Link to SwissProt release 42.0</a><BR><BR>";
  print "&nbsp;&nbsp;&nbsp;&nbsp;<a href='ftp://ftp.expasy.org/databases/prosite/old_releases/' style='color:blue' target='_blank'>Link to PROSITE release 20.33</a></p>";
  print "<BR><BR><BR>";

  __print_foot();


# ### Main ###

form=cgi.FieldStorage()

try:
  if not (form.has_key("page") or form.has_key("id")):
    __print_start_screen();
  elif form.has_key("id"):
    __printResultOrWaitingScreen(form["id"].value);
  elif form.has_key("page"):
    if form["page"].value == "predict":
      __print_prediction_page();
    elif form["page"].value == "help":
      __print_help_page();
    elif form["page"].value == "info":
      __print_info_page();
    else:
      __print_start_screen("Unkown value of parameter 'page'. You have been redirected to the start page.");
except:
  print "\n\n<PRE>"
  traceback.print_exc()
