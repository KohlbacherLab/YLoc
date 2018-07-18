import re, os, sys;
from math import log,exp,sin,cos,sqrt,ceil;

class SortingSignals(object):

	"Class creates features from known sorting signals"
	
	# neccessary for presto...not used anymore
	PATH_PRESTO = "python /share/usr/briese/Pythonscrips/presto/presto_V5.py";
	PATH_TMP = "/share/usr/briese/Pythonscrips/tmp/";
	PATH_MOTIFS = "/share/usr/briese/Pythonscrips/motifs/";
	
	def __init__(self):
		# add NLS DB regular expression terms
		self.__nlsdb=[]
		self.__nlsdb.append("APKRKSGVSKC")
		self.__nlsdb.append("APTKRKGS")
		self.__nlsdb.append("CKRKTTNADRRKA")
		self.__nlsdb.append("CYGSKNTGAKKRKIDDA")
		self.__nlsdb.append("DK[QL]KK[QL]")
		self.__nlsdb.append("DR[MN]KKKKE")
		self.__nlsdb.append("D[KR].{0,1}[QL][RK]{2,3}R")
		self.__nlsdb.append("EEDGPQKKKRRL")
		self.__nlsdb.append("EYLSRKGKLEL")
		self.__nlsdb.append("GGG.{3}KNRR.{6}RGGRN")
		self.__nlsdb.append("GKKKYKLKH")
		self.__nlsdb.append("GKKRSKA")
		self.__nlsdb.append("GRKRKKRT")
		self.__nlsdb.append("GR[RK]{2,4}..[RK][QL]")
		self.__nlsdb.append("G{2,4}[RK].{1,3}G{3}")
		self.__nlsdb.append("HKKKKIRTSPTFTTPKTLRLRRQPKYPRKSAPRRNKLDHY")
		self.__nlsdb.append("HRIEEKRKRTYETFKSI")
		self.__nlsdb.append("HRKYEAPRH.{6}PRKR")
		self.__nlsdb.append("IKYFKKFPKD")
		self.__nlsdb.append("KAKRQR")
		self.__nlsdb.append("KDCVINKHHRNRCQYCRLQR")
		self.__nlsdb.append("KHLKGR")
		self.__nlsdb.append("KHRKHPG")
		self.__nlsdb.append("KKEKKKSKK")
		self.__nlsdb.append("KKKKKEEEGEGKKK")
		self.__nlsdb.append("KKKKK.{3,6}KK")
		self.__nlsdb.append("KKKKRKREK")
		self.__nlsdb.append("KKKKR[KR]")
		self.__nlsdb.append("KKKRERLD")
		self.__nlsdb.append("KKKRRSREK")
		self.__nlsdb.append("KKKR[KR][VPL]")
		self.__nlsdb.append("KKKYKLK")
		self.__nlsdb.append("KKP.{6,9}K.{1,3}RK")
		self.__nlsdb.append("KKQTTLAFKPIKKGKKR")
		self.__nlsdb.append("KKRKRT")
		self.__nlsdb.append("KKRKR[KR]")
		self.__nlsdb.append("KKRKR[ST]")
		self.__nlsdb.append("KKRRK")
		self.__nlsdb.append("KKRR[DE]K")
		self.__nlsdb.append("KKRR.K")
		self.__nlsdb.append("KKSKKGRQEALERLKKA")
		self.__nlsdb.append("KK[MNQSTC]R[MNQSTC]K[MNQSTC]")
		self.__nlsdb.append("KK.R.{3,5}R[PVL]K")
		self.__nlsdb.append("KK.{1,7}K[PL][PLIV]KK")
		self.__nlsdb.append("KK.{15}KKRK")
		self.__nlsdb.append("KRAAEDDEDDDVDTKKQK")
		self.__nlsdb.append("KRGRGRPRK")
		self.__nlsdb.append("KRKKEMANKSAPEAKKKK")
		self.__nlsdb.append("KRKRRP")
		self.__nlsdb.append("KRK.{0,8}KR[PL]K")
		self.__nlsdb.append("KRK.{10,14}[KR]{3,}?.[KR]K")
		self.__nlsdb.append("KRK.{11}KKKSKK")
		self.__nlsdb.append("KRK.{2,4}DRRK")
		self.__nlsdb.append("KRK.{22}KELQKQITK")
		self.__nlsdb.append("KRK.{5,10}KK[PL]K")
		self.__nlsdb.append("KRMRNRIAASKCRKRKL")
		self.__nlsdb.append("KRPAATKKAGQAKKKK")
		self.__nlsdb.append("KRPACTLKPECVQQLLVCSQEAKK")
		self.__nlsdb.append("KRPAEDMEEEQAFKRSR")
		self.__nlsdb.append("KRPMNAFIVWSRDQRRK")
		self.__nlsdb.append("KRPMNAFMVWAQAARRK")
		self.__nlsdb.append("KRPRP")
		self.__nlsdb.append("KRQR.{20}KKSKK")
		self.__nlsdb.append("KRSAEGGNPPKPLKKLR")
		self.__nlsdb.append("KR[GPL]R[GPL]R[GLP]RK")
		self.__nlsdb.append("KR[MNSQ]R[MNSQ]R")
		self.__nlsdb.append("KR[PLV][GA]KRK[PL]")
		self.__nlsdb.append("KR[RK][RK].{2,4}[RK].{0,2}R.{3,5}[RK].{0,2}[RK].{0,2}[RK][RK]K")
		self.__nlsdb.append("KR[ST]R..R{2,4}[QL]K")
		self.__nlsdb.append("KR.R.R.{2,6}RKRK")
		self.__nlsdb.append("KR.R..RRLK")
		self.__nlsdb.append("KR.[DE][KR][KR].K")
		self.__nlsdb.append("KR..KK.K[DE]")
		self.__nlsdb.append("KR.{1,3}H.{3,5}R[LQ]RR")
		self.__nlsdb.append("KR.{10}KKKL")
		self.__nlsdb.append("KR.{7,9}PQPKKKP")
		self.__nlsdb.append("KR.{9}KTKK")
		self.__nlsdb.append("KR{2,4}.{3,6}[RK]{2,4}.{0,2}KR")
		self.__nlsdb.append("KR{3,}?[LVI]")
		self.__nlsdb.append("KSKKKAQ")
		self.__nlsdb.append("KTRKHRG")
		self.__nlsdb.append("KVNSRKRRKEVPGPNGATEED")
		self.__nlsdb.append("KVTKRKHDNEGSGSKRPK")
		self.__nlsdb.append("K[GA]K[AG]KK[AG]")
		self.__nlsdb.append("K[IVQM]RR[VI][STK]L")
		self.__nlsdb.append("K[KR][KR]RR[KR]")
		self.__nlsdb.append("K[KR][QMN][RK]R[QMN]R")
		self.__nlsdb.append("K[MNQ]RR[PLVI]K[PL]")
		self.__nlsdb.append("K[PLMN]RRK[MNQ]")
		self.__nlsdb.append("K[PL]K{2,3}.{1,3}[RK]{2,4}.{6,9}K[KR]")
		self.__nlsdb.append("K[PL]K{3,}?.KK")
		self.__nlsdb.append("K[RK]{2,4}[ST]H")
		self.__nlsdb.append("K[RK]{2,}?[QL].{3,8}R{3}")
		self.__nlsdb.append("K[RK]{3,5}.{11,18}[RK]K.{2,3}K")
		self.__nlsdb.append("K.KRQR")
		self.__nlsdb.append("K.K.K.....RKK")
		self.__nlsdb.append("K.[PLV][RK][RK]RK")
		self.__nlsdb.append("K..K.K.K.....RKK")
		self.__nlsdb.append("K{1,}?R{2,}?[QM]R{2,}")
		self.__nlsdb.append("K{3,4}R{2,3}")
		self.__nlsdb.append("LEKKVKKKFDWCA")
		self.__nlsdb.append("LKDVRKRKLGPGH")
		self.__nlsdb.append("LKKIKQ")
		self.__nlsdb.append("LKRKLQR")
		self.__nlsdb.append("LKRPRSPSS")
		self.__nlsdb.append("MAPSAKATAAKKAVVKGTNGKKALKVRTSATFRLPKTLKLAR")
		self.__nlsdb.append("MNKIPIKDLLNPG")
		self.__nlsdb.append("MPKTRRRPRRSQRKRPPT")
		self.__nlsdb.append("MPTEERVRKRKESNRESARRSRYRKAAHLK")
		self.__nlsdb.append("NQSSNFGPMKGGNFGGRSSGPYGGGGQYFAKPRNQGGY")
		self.__nlsdb.append("N[QR]RQ[RK][EG]KR[IVLS]")
		self.__nlsdb.append("PAAKRVKLD")
		self.__nlsdb.append("PAKRARRGYK")
		self.__nlsdb.append("PKKARED")
		self.__nlsdb.append("PKKKRKV")
		self.__nlsdb.append("PKKK.RK")
		self.__nlsdb.append("PKKNRLRRP")
		self.__nlsdb.append("PKRPRDRHDGELGGRKRARG")
		self.__nlsdb.append("PKRP.{5,8}L.{2,4}R.K.K")
		self.__nlsdb.append("PK[KR][KRP][RAK][KT][VSE]")
		self.__nlsdb.append("PLLKKIKQ")
		self.__nlsdb.append("PNKKKRK")
		self.__nlsdb.append("PPQKKIKS")
		self.__nlsdb.append("PPRIYPQLPSAPT")
		self.__nlsdb.append("PPRKKRTVV")
		self.__nlsdb.append("PPVKRERTS")
		self.__nlsdb.append("PQPKKKP")
		self.__nlsdb.append("PQSRKKLR")
		self.__nlsdb.append("PRGRRQPIPKARQP")
		self.__nlsdb.append("PRPRKIPR")
		self.__nlsdb.append("PRRGPR")
		self.__nlsdb.append("PRRRK")
		self.__nlsdb.append("PYLNKRKGKP")
		self.__nlsdb.append("P.[PQLVMN][KR]{2,3}.KQ")
		self.__nlsdb.append("QNRR.K.[RK][RK][DQE]")
		self.__nlsdb.append("QRKRQK")
		self.__nlsdb.append("Q[RK][HRK][RK].RR")
		self.__nlsdb.append("RAIKRRPGLDFDDDGEGNSKFLR")
		self.__nlsdb.append("REKKEKEQKEKCA")
		self.__nlsdb.append("RER[MNQ]K.{4,8}R[MNQ]RR")
		self.__nlsdb.append("RGRRRRQR")
		self.__nlsdb.append("RH[RK]H.{2,4}[RK]{2,4}[PL]R")
		self.__nlsdb.append("RIRKKLR")
		self.__nlsdb.append("RKCLQAGMNLEARKTKK")
		self.__nlsdb.append("RKEWLTNFMEDRRQRKL")
		self.__nlsdb.append("RKKRKR")
		self.__nlsdb.append("RKKRK.{9}KAKKSK")
		self.__nlsdb.append("RKKRRQRRR")
		self.__nlsdb.append("RKKRR.R")
		self.__nlsdb.append("RKRAFHGDDPFGEGPPDKK")
		self.__nlsdb.append("RKRIREDRKATTAQKVQQMKQRLNENERKRKR")
		self.__nlsdb.append("RKRKK")
		self.__nlsdb.append("RKRKKMPASQRSKRRKT")
		self.__nlsdb.append("RKRKR[KR]")
		self.__nlsdb.append("RKRRR")
		self.__nlsdb.append("RKR[PLQMN]R[PLQMN]R")
		self.__nlsdb.append("RKR.{12,16}RRKK")
		self.__nlsdb.append("RKR{3,5}[ST]")
		self.__nlsdb.append("RK[IVE]W[ML][TQR]N[HF]")
		self.__nlsdb.append("RK[PL][PLV]KK[RKH]")
		self.__nlsdb.append("RK[RK][QML][RK].R")
		self.__nlsdb.append("RK]{2,4}[PL][RK].{7,11}[RK][QL]KH")
		self.__nlsdb.append("RK.{7,12}RK[STMNQ]KK")
		self.__nlsdb.append("RLKKLKCSK.{19}KTKR")
		self.__nlsdb.append("RPRRK")
		self.__nlsdb.append("RQARRNRRRRWR")
		self.__nlsdb.append("RRERNKMAAAKCRNRRR")
		self.__nlsdb.append("RRER[MNQ]K.{4,8}R[MNQ]RRR")
		self.__nlsdb.append("RRER.{4}RPRKIPR")
		self.__nlsdb.append("RRKGKEK")
		self.__nlsdb.append("RRK.{3,5}R[DE]R{3,}?[PLV]")
		self.__nlsdb.append("RRK.{5,7}RRR")
		self.__nlsdb.append("RRMKWKK")
		self.__nlsdb.append("RRPS.{22}RRKRQ")
		self.__nlsdb.append("RRRKKR")
		self.__nlsdb.append("RRRK[STC]K")
		self.__nlsdb.append("RRRRRR")
		self.__nlsdb.append("RRRRRR.{0,2}R")
		self.__nlsdb.append("RRRRR.RR")
		self.__nlsdb.append("RRR[LP]..R[PLQ]")
		self.__nlsdb.append("RRR[PL]RK")
		self.__nlsdb.append("RRR.{11}KRRK")
		self.__nlsdb.append("RRR{3,5}T")
		self.__nlsdb.append("RRSMKRK")
		self.__nlsdb.append("RR[PLIV]RK.K")
		self.__nlsdb.append("RR[PLQMN].RRRR")
		self.__nlsdb.append("RR[TS].[QK][KR][KNS]")
		self.__nlsdb.append("RR[TS].[QK][KR][KN]")
		self.__nlsdb.append("RR.KR.K[PLV]")
		self.__nlsdb.append("RR.RRRRR")
		self.__nlsdb.append("RR.R[PVL]RK")
		self.__nlsdb.append("RR.R.K.KQ")
		self.__nlsdb.append("RR.R.RKQ")
		self.__nlsdb.append("RR..KRK")
		self.__nlsdb.append("RR.{0,1}RRRRR")
		self.__nlsdb.append("RVHPYQR")
		self.__nlsdb.append("R[GA][IVLP]KRR")
		self.__nlsdb.append("R[GA].{0,2}[GA]R[GA].[GA]R[GA]")
		self.__nlsdb.append("R[GVLIP]RRRR.R")
		self.__nlsdb.append("R[IVLP][IVLP]KRR")
		self.__nlsdb.append("R[KR]RRRR.R")
		self.__nlsdb.append("R[KR][RK].{0,2}[RK].{0,2}[RK].{3,5}[RK].{0,2}[RK][RK][RK][RK][PMQL]")
		self.__nlsdb.append("R[KR]{3,4}K[DE]")
		self.__nlsdb.append("R[MNQ]RRRR.R")
		self.__nlsdb.append("R[MNQ].{4,8}R[MNQ]RR")
		self.__nlsdb.append("R[PL].G.[KR][KR].K")
		self.__nlsdb.append("R[PL]..[KR]{2,}?..[KR]V")
		self.__nlsdb.append("R[QMPL]RR[DE]R")
		self.__nlsdb.append("R[RK]K[RK]KR")
		self.__nlsdb.append("R[RK].[KR].[RK]{2,}?[DE]")
		self.__nlsdb.append("R[RK].{4,6}[RK][RK].[RK].{1,3}[RK][RK][PLQ]")
		self.__nlsdb.append("R[RK]{2,4}[PL][RK][MNQ]R")
		self.__nlsdb.append("R[RK]{2,4}.{15,19}[RK]{2,4}[QLM]K")
		self.__nlsdb.append("R[RK]{3,}?[DE]K")
		self.__nlsdb.append("R[STCMNQ]R[STCMNQ]KR")
		self.__nlsdb.append("R.KKKK[DE]")
		self.__nlsdb.append("R.RR.{4,6}RKK")
		self.__nlsdb.append("R.RSRS.{0,1}R.R")
		self.__nlsdb.append("R.R.R.R.R.R")
		self.__nlsdb.append("R.R.R.R.R.R.K")
		self.__nlsdb.append("R.R{2,}?[QL].[ST]R")
		self.__nlsdb.append("R.[KR][KR]K[PLQM]R")
		self.__nlsdb.append("R.[KR][KR][KR]..RKKR")
		self.__nlsdb.append("R.{2,3}H.{3,5}RRRR")
		self.__nlsdb.append("R.{2,3}RRRRRR")
		self.__nlsdb.append("R{2,3}K{3,4}[PLRKE]")
		self.__nlsdb.append("R{2,3}.K{2,3}R[ST]")
		self.__nlsdb.append("R{2,}?PR{3,}?")
		self.__nlsdb.append("R{2,}?[QMN]R{3,}?")
		self.__nlsdb.append("SANKVTKNKSNSSPYLNKRGKPGPDS")
		self.__nlsdb.append("SDKKVRSRLIECA")
		self.__nlsdb.append("SKRVAKRKL")
		self.__nlsdb.append("S.GTKRSY..M")
		self.__nlsdb.append("TEKK[QG]KSILYDCA")
		self.__nlsdb.append("TKRS...M")
		self.__nlsdb.append("T[PLV]KRC")
		self.__nlsdb.append("VNEAFETLKRC")
		self.__nlsdb.append("VSRKRPR")
		self.__nlsdb.append("WKQ[KR]RKF")
		self.__nlsdb.append("YKRPCKRSFIRFI")
		self.__nlsdb.append("YLTQETNKVETYKEQPLKTPGKKKKGKP")
		self.__nlsdb.append("YNNQSSNFGPMKGGN")
		self.__nlsdb.append("[DE]KK[PL][GL]K[GL]")
		self.__nlsdb.append("[DE]KR[MQN]R[MQN]R")
		self.__nlsdb.append("[DE]K[NIF]RR[DEK][STMNQ]")
		self.__nlsdb.append("[DE]K.RRK[MNQ]")
		self.__nlsdb.append("[DE]RKRR[DEPLQ]")
		self.__nlsdb.append("[DE]R.KKKK")
		self.__nlsdb.append("[DE]R{2,4}.RK[PL]")
		self.__nlsdb.append("[DE][KR]RR[KR][FYW]")
		self.__nlsdb.append("[DE][RK]{2,4}[GA]R[PL][GA]")
		self.__nlsdb.append("[DE][RK]{3,}?.[KR]{2,}?[PL]")
		self.__nlsdb.append("[DE][ST][PL]KR[STC]")
		self.__nlsdb.append("[ED]R{4,}?[ED]")
		self.__nlsdb.append("[GAPLV]RKRKKR")
		self.__nlsdb.append("[GA]K.KKK[MNQ]")
		self.__nlsdb.append("[GA]R.[RK].[RK][RK].[QM]")
		self.__nlsdb.append("[GA][KR]KR.[KR][GA]")
		self.__nlsdb.append("[KAR]TPIQKHWRPTVLTEGPPVKIRIETGEWE[KA]")
		self.__nlsdb.append("[KR]G{2,}?..G{3,}?[RK]")
		self.__nlsdb.append("[KR]KRKK")
		self.__nlsdb.append("[KR]XXKNKX{6,8}K[KR]")
		self.__nlsdb.append("[KR][DE][KR][DE]..[KR]{4,}?")
		self.__nlsdb.append("[KR][KR][KR][KR][KR][KR][KR]")
		self.__nlsdb.append("[KR][KR][QMN]R[RK][QMN]R")
		self.__nlsdb.append("[KR][KR].[KR][KR][KR].[KR][KR]")
		self.__nlsdb.append("[KR]{2,3}..KR[KR][QLM]")
		self.__nlsdb.append("[KR]{2,}?[PL].{1,4}[KR]{2,}?.{1,5}K{3,}?")
		self.__nlsdb.append("[KR]{2}.{0,1}[KR]{2,4}.{25,34}K{2,4}.{1,2}K")
		self.__nlsdb.append("[KR]{4}.{20,24}K{1,4}.K")
		self.__nlsdb.append("[LF][STK][VIQM][KR]R[QMVI][STK]L")
		self.__nlsdb.append("[MI]VWSRD[HEQ]RRK")
		self.__nlsdb.append("[PLQMKR]R[KR][QM][KR]R.K")
		self.__nlsdb.append("[PLQMNKR]K[KR][KR]R.K[PLQMNKR]")
		self.__nlsdb.append("[PLQ]K[RK].{1,2}[RK].{3,6}[RK][RK].{1,2}[RK].{1,2}[RK][RK]")
		self.__nlsdb.append("[PLQ][KR].{3,4}KKRK")
		self.__nlsdb.append("[PLV]K[RK].[QMN][RK]R")
		self.__nlsdb.append("[PLV]K[RK].[RK][RK][RK][PL]")
		self.__nlsdb.append("[PLV]RK[ST]R[DE]K")
		self.__nlsdb.append("[PL]K..KRR")
		self.__nlsdb.append("[PL]RKRK[PL]")
		self.__nlsdb.append("[PL]R[DE]K[DE]R")
		self.__nlsdb.append("[PL][KR]{5,7}[PL]")
		self.__nlsdb.append("[PL][PL].[KR]R[DE][KR][QST]")
		self.__nlsdb.append("[PL][RK][RK][DEP]R[RK][FYW]")
		self.__nlsdb.append("[PL][RK][RK][KR][GAPL][RK][STQM]")
		self.__nlsdb.append("[PL][RK]{2,3}K[PLI][RK].[PLI].K")
		self.__nlsdb.append("[PL]..KR[IV]K[PL][DE]")
		self.__nlsdb.append("[PVLI][RK][RK][RK][RK][RK][QMN]K")
		self.__nlsdb.append("[QL]K{2,4}.{8,12}[RK][QL][RK][QL]KR")
		self.__nlsdb.append("[QL].KR.K.KK")
		self.__nlsdb.append("[QMN]R[RK].K.[RK][RK]")
		self.__nlsdb.append("[RK]H[RK]...[RK]{2,4}.R")
		self.__nlsdb.append("[RK]K{2,4}.[RK][QL][RK][PL]")
		self.__nlsdb.append("[RK]R[MS]K.K[KR]")
		self.__nlsdb.append("[RK][PLIV][KR][RK]{2,4}[PLVI]R")
		self.__nlsdb.append("[RK].[RK].[KR].{4,6}RKK")
		self.__nlsdb.append("[RK]{2,4}.{1,2}[RK].{0,2}[RK].{3,5}[RK].{0,2}[RK][RK]{2,4}[PL]")
		self.__nlsdb.append("[RK]{2,4}.{2,4}[QLM][RK].{2,3}[RK]KR")
		self.__nlsdb.append("[RK]{3,}?.[RK].[RK].{4,9}[RK]{3,}?")
		self.__nlsdb.append("[RK]{3,}?.{8,16}[RK]{4,}?")
		self.__nlsdb.append("[RK]{4,}?[QMNPL][RK].{3,4}[RK]{2}")
		self.__nlsdb.append("[STQM]RKRK[STQM]")
		self.__nlsdb.append("[STQM]RKRR[STQM]")
		self.__nlsdb.append("[STQM]RRRK[STQM]")
		self.__nlsdb.append("[ST]G.{1,3}G{3,}?.{1,2}G{3,}?[ST]")
		self.__nlsdb.append("[TS][RK]KK[VLI]R[PL]")
		self.__nlsdb.append("[YFW]RRRR[PL]")

	def statNLS(self):
		empty = {'A': 0, 'C':0,'D' :0, 'E' :0, 'F':0, 'G':0,'H':0,'I':0,'K':0,'L':0,
			'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':0};
		count = 0;
		for elem in self.__nlsdb:
			for aa in elem:
				if aa in empty.keys():
					empty[aa] += 1;
					count += 1;
		for elem in empty.keys():
			print(elem+" => "+str(float(empty[elem])/count));
		
	def writeToFile(self):
		f = open("nlsdb.motifs","w");
		it = 0;
		for motif in self.__nlsdb:
			f.write(str(it)+" "+motif+" # NLS DB Motif Nr. "+str(it)+"\n");
			it += 1;
		f.close();

	def getAllSignals(self, sequence):
		features = [];

		features.append(self.getCleavageSite(sequence));
		#features.append(self.getPrestoNLSDBSignals(sequence));
		features.append(self.getNLSDBSignals(sequence)); #good
		features.append(self.getMonoNLSSignal(sequence));
		features.append(self.getMonoNLSSignal1(sequence)); #good
		features.append(self.getMonoNLSSignal2(sequence)); #good
		features.append(self.getMonoNLSSignal3(sequence)); #good
		
		features.append(self.getMonoNLSSignal4(sequence));
		features.append(self.getMonoNLSSignal5(sequence)); #good
		features.append(self.getMonoNLSSignal6(sequence)); 
		features.append(self.getMonoNLSSignal7(sequence)); 
		features.append(self.getMonoNLSSignal8(sequence)); 
		features.append(self.getMonoNLSSignal9(sequence)); 
		features.append(self.getMonoNLSSignal10(sequence)); 
		features.append(self.getMonoNLSSignal11(sequence));
		features.append(self.getMonoNLSSignal12(sequence));
		features.append(self.getMonoNLSSignal13(sequence));
		features.append(self.getMonoNLSSignal14(sequence));
		features.append(self.getMonoNLSSignal15(sequence));
		features.append(self.getMonoNLSSignal16(sequence));
		features.append(self.getMonoNLSSignal17(sequence));
		features.append(self.getMonoNLSSignal18(sequence));
		features.append(self.getMonoNLSSignal19(sequence));
		features.append(self.getMonoNLSSignal20(sequence));
		features.append(self.getMonoNLSSignal21(sequence));
		features.append(self.getMonoNLSSignal22(sequence));
		features.append(self.getMonoNLSSignal23(sequence));
		features.append(self.getMonoNLSSignal24(sequence));
		features.append(self.getMonoNLSSignal25(sequence));
		features.append(self.getMonoNLSSignal26(sequence));
		features.append(self.getMonoNLSSignal27(sequence));
		features.append(self.getBipartiteNLSSignal(sequence));
		features.append(self.getBipartiteNLSSignal1(sequence));
		features.append(self.getBipartiteNLSSignal2(sequence));
		
		features.append(self.getBipartiteNLSSignal3(sequence));
		features.append(self.getERSignal(sequence));  #good
		features.append(self.getERSignal2(sequence));
		features.append(self.getERSignal3(sequence));
		features.append(self.getERSignal4(sequence));
		features.append(self.getERSignal5(sequence));
		features.append(self.getPeroxiSignal(sequence));
		features.append(self.getPeroxiSignal2(sequence));
		features.append(self.getPeroxiSignal3(sequence)); #good
		features.append(self.getPeroxiSignal4(sequence));
		features.append(self.getPeroxiSignal5(sequence));
		features.append(self.getPeroxiSignal6(sequence));
		features.append(self.getPeroxiSignal7(sequence));
		features.append(self.getPeroxiSignal8(sequence));
		features.append(self.getPeroxiSignal9(sequence));
		
		features.append(self.getPeroxiSignalPlant(sequence)); #good
		features.append(self.getVacuoleSignal(sequence));
		features.append(self.getVacuoleSignal2(sequence));
		features.append(self.getVacuoleSignal3(sequence));
		features.append(self.getVacuoleSignal4(sequence));
		features.append(self.getVacuoleSignal5(sequence));
		features.append(self.getVacuoleSignal6(sequence)); #good
		features.append(self.getVacuoleSignal7(sequence)); #some
		features.append(self.getVacuoleSignal8(sequence));
		features.append(self.getVacuoleSignal9(sequence));
		features.append(self.getVacuoleSignal10(sequence));
		features.append(self.getVacuoleSignal11(sequence));
		features.append(self.getVacuoleSignal12(sequence));
		features.append(self.getVacuoleSignal13(sequence));
		features.append(self.getVacuoleSignal14(sequence));
		features.append(self.getVacuoleSignal15(sequence));
		features.append(self.getVacuoleSignal16(sequence));
		features.append(self.getVacuoleSignal17(sequence));
		features.append(self.getVacuoleSignal18(sequence)); #some
		features.append(self.getVacuoleSignal19(sequence)); #some
		features.append(self.getVacuoleSignal20(sequence)); #some
		features.append(self.getVacuoleSignal21(sequence));
		features.append(self.getVacuoleSignal22(sequence));
		features.append(self.getVacuoleSignal23(sequence)); #good
		features.append(self.getVacuoleSignal24(sequence)); #some
		features.append(self.getVacuoleSignal25(sequence)); 
		features.append(self.getVacuoleSignal26(sequence)); #good
		features.append(self.getVacuoleSignal27(sequence));
		features.append(self.getVacuoleSignal28(sequence));
		features.append(self.getVacuoleSignal29(sequence)); #good
		features.append(self.getVacuoleSignal30(sequence));
		features.append(self.getVacuoleSignal31(sequence));
		features.append(self.getVacuoleSignal32(sequence));
		features.append(self.getVacuoleSignal33(sequence));
		features.append(self.getVacuoleSignal34(sequence));
		features.append(self.getVacuoleSignal35(sequence));
		features.append(self.getVacuoleSignal36(sequence));
		features.append(self.getVacuoleSignal37(sequence));
		
		features.append(self.getLysosomalSignal1(sequence));
		features.append(self.getLysosomalSignal2(sequence));
		features.append(self.getLysosomalSignal3(sequence));
		features.append(self.getLysosomalSignal4(sequence));
		features.append(self.getLysosomalSignal5(sequence));
		features.append(self.getLysosomalSignal6(sequence));
		features.append(self.getLysosomalSignal7(sequence));
		
		#6+23+26 42
		#6+18+24+26 38
	
		
		features.append(self.getTransGolgiSignal(sequence));
		features.append(self.getGolgiSignal(sequence));
		features.append(self.getGolgiSignal2(sequence));
		features.append(self.getGolgiSignal3(sequence));
		features.append(self.getGolgiSignal4(sequence));
		features.append(self.getGolgiSignal5(sequence));
		features.append(self.getGolgiSignal6(sequence));
		features.append(self.getGolgiSignal7(sequence));
		features.append(self.getGolgiSignal8(sequence));
		features.append(self.getGolgiSignal9(sequence));
		features.append(self.getGolgiSignal10(sequence));
		features.append(self.getGolgiSignal11(sequence));
		features.append(self.getGolgiSignal12(sequence));
		features.append(self.getTransMembraneSignal(sequence));
		features.append(self.getTransMembraneSignal2(sequence));
		features.append(self.getTransMembraneSignal3(sequence));
				
		features.append(self.getEndocytosisSignal(sequence));
		features.append(self.getGlycoSiteSignal(sequence));
		features.append(self.getAnchorSignal(sequence)); #good
		features.append(self.getNESSignal(sequence)); #good
		
		features.append(self.getMIPSignal(sequence));
		features.append(self.getMIPSignal2(sequence));
		features.append(self.getMIPSignal3(sequence));
		features.append(self.getMIPSignal4(sequence)); #good!!!
		features.append(self.getMIPSignal5(sequence));
		features.append(self.getMIPSignal6(sequence));
		features.append(self.getMIPSignal7(sequence));
		features.append(self.getMIPSignal8(sequence));
		features.append(self.getMIPSignal9(sequence));
		features.append(self.getMIPSignal10(sequence));
		features.append(self.getMIPSignal11(sequence)); #ok
		features.append(self.getMIPSignal12(sequence));
		features.append(self.getMIPSignal13(sequence));
		features.append(self.getAliphaticHelix(sequence)); #good
		features.append(self.getAliphaticHelixNormed(sequence)); #good
		features.append(self.getAliphaticSheet(sequence)); 
		features.append(self.getAliphaticSheetNormed(sequence)); 
		features.append(self.getChlorSignal(sequence));
		features.append(self.getChlorSignal2(sequence));
		features.append(self.getVesicularSignal(sequence));
		features.append(self.getSecretedSignal1(sequence));

			
		return features;

	def getAllLabels(self):
	  #"Presto_NLSDB_prob"
		return [ "cleavage_site", "NLSDB_Signal", "Mono_NLS_Signal","Mono_NLS_Signal1","Mono_NLS_Signal2","Mono_NLS_Signal3","Mono_NLS_Signal4","Mono_NLS_Signal5","Mono_NLS_Signal6","Mono_NLS_Signal7","Mono_NLS_Signal8","Mono_NLS_Signal9","Mono_NLS_Signal10","Mono_NLS_Signal11","Mono_NLS_Signal12","Mono_NLS_Signal13","Mono_NLS_Signal14","Mono_NLS_Signal15","Mono_NLS_Signal16","Mono_NLS_Signal17","Mono_NLS_Signal18","Mono_NLS_Signal19","Mono_NLS_Signal20","Mono_NLS_Signal21","Mono_NLS_Signal22","Mono_NLS_Signal23","Mono_NLS_Signal24","Mono_NLS_Signal25", "Mono_NLS_Signal26", "Mono_NLS_Signal27","Bipartite_NLS_Signal","Bipartite_NLS_Signal1","Bipartite_NLS_Signal2","Bipartite_NLS_Signal3", "ER_Signal","ER_Signal2","ER_Signal3", "ER_Signal4","ER_Signal5","Peroxi_Signal1", "Peroxi_Signal2", "Peroxi_Signal3", "Peroxi_Signal4", "Peroxi_Signal5","Peroxi_Signal6","Peroxi_Signal7","Peroxi_Signal8","Peroxi_Signal9","Peroxi_Signal_Plant", "Vacuole_Signal", "Vacuole_Signal2","Vacuole_Signal3","Vacuole_Signal4","Vacuole_Signal5","Vacuole_Signal6","Vacuole_Signal7","Vacuole_Signal8","Vacuole_Signal9","Vacuole_Signal10","Vacuole_Signal11","Vacuole_Signal12","Vacuole_Signal13","Vacuole_Signal14","Vacuole_Signal15","Vacuole_Signal16","Vacuole_Signal17","Vacuole_Signal18","Vacuole_Signal19","Vacuole_Signal20","Vacuole_Signal21","Vacuole_Signal22","Vacuole_Signal23","Vacuole_Signal24","Vacuole_Signal25","Vacuole_Signal26","Vacuole_Signal27","Vacuole_Signal28","Vacuole_Signal29","Vacuole_Signal30","Vacuole_Signal31","Vacuole_Signal32","Vacuole_Signal33","Vacuole_Signal34","Vacuole_Signal35","Vacuole_Signal36","Vacuole_Signal37","Lysosomal_Signal1","Lysosomal_Signal2","Lysosomal_Signal3","Lysosomal_Signal4","Lysosomal_Signal5","Lysosomal_Signal6","Lysosomal_Signal7","Trans_Golgi_Signal", "Golgi_signal", "Golgi_signal2","Golgi_signal3","Golgi_signal4","Golgi_signal5","Golgi_signal6","Golgi_signal7", "Golgi_signal8","Golgi_signal9","Golgi_signal10","Golgi_signal11","Golgi_signal12","Transmembrane_signal","Transmembrane_signal2","Transmembrane_signal3", "Endocytosis_Signal", "Glycosite_Signal", "Anchor_Signal","NES_signal" ,"MIP_signal","MIP_signal2","MIP_signal3","MIP_signal4","MIP_signal5","MIP_signal6","MIP_signal7","MIP_signal8","MIP_signal9","MIP_signal10","MIP_signal11","MIP_signal12","MIP_signal13","aliphatic_helix","aliphatic_helix_normed","aliphatic_sheet","aliphatic_sheet_normed","chloro_signal","chloro_signal2","vesicular_signal","secreted_signal"];


	def getPrestoNLSDBSignals(self, sequence):
		fasta_file = open(self.PATH_TMP+"tmp.fasta","w");
		fasta_file.write(">NLSDBSignal_Sequence, delete me\n");
		fasta_file.write(sequence);
		fasta_file.close();
		
		print(".");
		os.system("%s searchdb=%stmp.fasta motifs=%snls.motifs mismatch=1,0 i=-1 outfile=%stmp.tdt > %spresto.out" % (self.PATH_PRESTO,self.PATH_TMP,self.PATH_MOTIFS,self.PATH_TMP,self.PATH_TMP));
			
		result_file = open(self.PATH_TMP+"tmp.tdt.presto.tdt","r");
		line = result_file.readline();
		line = result_file.readline();
		overallprob = 1;
		while line:
			line_s = line.split("\t");
			e_value = float(line_s[15].strip(" "));
			prob = 1/(1+exp(-1*log(0.00001/e_value)));
			print(prob);
			overallprob = overallprob * (1-prob);
			line = result_file.readline();
		print(1-overallprob);
		return (1-overallprob);
		
	def getNLSDBSignals(self, sequence):
		nr_signals = 0;
		for nls in self.__nlsdb:
			if len(re.findall(nls,sequence))>0:
				nr_signals += 1;
		return nr_signals;

	def getMonoNLSSignal(self, sequence):
		nls_mono = len(re.findall("K[KR].[KR]",sequence)); # multiLoc1
		nls_mono += len(re.findall("[KR]{4}[PGRHKED]",sequence)); # Ema02
		nls_mono += len(re.findall("[KR]{4,7}",sequence)); # Seb
		nls_mono += len(re.findall("[KR]{2}.[KR]{2,4}",sequence)); # Seb
		nls_mono += len(re.findall("[KR]{4,7}",sequence)); # Ema02
		#own
		return nls_mono;
		
	def getMonoNLSSignal1(self, sequence):
		nls_mono = len(re.findall("K[KR].[KR]",sequence)); # multiLoc1
		return nls_mono;
		
	def getMonoNLSSignal2(self, sequence):
		nls_mono = len(re.findall("[KR]{4}[PGRHKED]",sequence)); # Ema02
		return nls_mono;
		
	def getMonoNLSSignal3(self, sequence):
		nls_mono = len(re.findall("[KR]{4,7}",sequence)); # Seb
		return nls_mono;
		
	def getMonoNLSSignal4(self, sequence):
		nls_mono = len(re.findall("[KR]{2}.[KR]{2,4}",sequence)); # Seb
		return nls_mono;
		
	def getMonoNLSSignal5(self, sequence):
		nls_mono = len(re.findall("K[KR]{4,5}[KR]",sequence[1:100])); 
		return nls_mono;
	
	def getMonoNLSSignal6(self, sequence):
		nls_mono = len(re.findall("K[KR]",sequence[1:100])); 
		return nls_mono;
		
	def getMonoNLSSignal7(self, sequence):
		nls_mono = len(re.findall("K[KR][KR]",sequence[1:100])); 
		return nls_mono;
		
	def getMonoNLSSignal8(self, sequence):
		nls_mono = len(re.findall("K[KR].{2,4}[KR]",sequence[-100:])); 
		return nls_mono;
	
	def getMonoNLSSignal9(self, sequence):
		nls_mono = len(re.findall("K.{2,4}[KR][KR]",sequence[-100:])); 
		return nls_mono;
		
	def getMonoNLSSignal10(self, sequence):
	
		nls_mono = self.getMonoNLSSignal6(sequence)+self.getMonoNLSSignal7(sequence)+self.getMonoNLSSignal8(sequence)
		nls_mono += self.getBipartiteNLSSignal1(sequence);
		return nls_mono;
  
	def getMonoNLSSignal11(self, sequence):
		nls_mono = len(re.findall("[KR]{4}[PRHK]",sequence)); 
		return nls_mono;
		
	def getMonoNLSSignal12(self, sequence):
		nls_mono = len(re.findall("[KR]{5}",sequence)); 
		return nls_mono;
		
	def getMonoNLSSignal13(self, sequence):
		nls_mono = len(re.findall("[KR]{4}",sequence)); 
		return nls_mono;
	
	def getMonoNLSSignal14(self, sequence):
		nls_mono = len(re.findall("[KR]{6}",sequence)); 
		return nls_mono;
	
	def getMonoNLSSignal15(self, sequence):
		nls_mono = len(re.findall("[KR]{3}[HP]",sequence)); 
		return nls_mono;

  
	def getMonoNLSSignal16(self, sequence):
		nls_mono = len(re.findall("P.K{3}",sequence)); 
		return nls_mono;
	
	def getMonoNLSSignal17(self, sequence):
		nls_mono = len(re.findall("P..K{3}",sequence)); 
		return nls_mono;
  
	def getMonoNLSSignal18(self, sequence):
		nls_mono = 0;
		#base distribution
		score = self.getMonoNLSSignal1(sequence);
		# multiply by 1 if 20% occurence in nuclear proteins
		nls_mono += self.getMonoNLSSignal6(sequence)*0.3;
		nls_mono += self.getMonoNLSSignal7(sequence);
		nls_mono += self.getMonoNLSSignal11(sequence);
		nls_mono += self.getMonoNLSSignal12(sequence)*3;
		nls_mono += self.getMonoNLSSignal13(sequence)*2;
		nls_mono += self.getMonoNLSSignal14(sequence)*4;
		nls_mono += self.getMonoNLSSignal15(sequence);
		nls_mono += self.getMonoNLSSignal16(sequence)*4;
		nls_mono += self.getMonoNLSSignal17(sequence)*2;
		
		score += nls_mono * 0.2;
		
		return score;

	def getMonoNLSSignal19(self, sequence):
		nls_mono = 0;
		#base distribution
		score = self.getMonoNLSSignal1(sequence);
		# multiply by 1 if 20% occurence in nuclear proteins
		nls_mono += self.getMonoNLSSignal6(sequence)*0.2;
		nls_mono += self.getMonoNLSSignal7(sequence)*0.26;
		nls_mono += self.getMonoNLSSignal7(sequence)*0.17;
		nls_mono += self.getMonoNLSSignal13(sequence)*0.83;
		nls_mono += self.getMonoNLSSignal14(sequence)*-0.76;
		nls_mono += self.getMonoNLSSignal15(sequence)*0.96;
		nls_mono += self.getMonoNLSSignal16(sequence)*0.7;
		
		score += nls_mono * 5;
		
		return score;		
  
  	def getMonoNLSSignal20(self, sequence):
		nls_mono = 0;

		# multiply by 1 if 20% occurence in nuclear proteins
		nls_mono += self.getMonoNLSSignal6(sequence)*0.2;
		nls_mono += self.getMonoNLSSignal7(sequence)*0.26;
		nls_mono += self.getMonoNLSSignal7(sequence)*0.17;
		nls_mono += self.getMonoNLSSignal13(sequence)*0.83;
		nls_mono += self.getMonoNLSSignal14(sequence)*-0.76;
		nls_mono += self.getMonoNLSSignal15(sequence)*0.96;
		nls_mono += self.getMonoNLSSignal16(sequence)*0.7;
		
		score = nls_mono * 5;
		
		return score;		
  
  	def getMonoNLSSignal21(self, sequence):
		nls_mono = 0;
		#base distribution
		score = self.getMonoNLSSignal1(sequence);
		# multiply by 1 if 20% occurence in nuclear proteins
		nls_mono += self.getMonoNLSSignal6(sequence)*0.3;
		nls_mono += self.getMonoNLSSignal7(sequence);
		nls_mono += self.getMonoNLSSignal11(sequence);
		nls_mono += self.getMonoNLSSignal12(sequence)*3;
		nls_mono += self.getMonoNLSSignal13(sequence)*2;
		nls_mono += self.getMonoNLSSignal14(sequence)*4;
		nls_mono += self.getMonoNLSSignal15(sequence);
		nls_mono += self.getMonoNLSSignal16(sequence)*4;
		nls_mono += self.getMonoNLSSignal17(sequence)*2;
		nls_mono += self.getNLSDBSignals(sequence)*2;
		
		score += nls_mono * 0.2;
		
		return score;
		
	def getMonoNLSSignal22(self, sequence):
		nls_mono = 0;

		# multiply by 1 if 20% occurence in nuclear proteins
		nls_mono += self.getMonoNLSSignal6(sequence)*0.2;
		nls_mono += self.getMonoNLSSignal7(sequence)*0.26;
		nls_mono += self.getMonoNLSSignal7(sequence)*0.17;
		nls_mono += self.getMonoNLSSignal13(sequence)*0.83;
		nls_mono += self.getMonoNLSSignal14(sequence)*-0.76;
		nls_mono += self.getMonoNLSSignal15(sequence)*0.96;
		nls_mono += self.getMonoNLSSignal16(sequence)*0.7;
		nls_mono += self.getNLSDBSignals(sequence)*0.55;
		
		score = nls_mono * 5;
		
		return score;
		
	def getMonoNLSSignal23(self, sequence):
		nls = len(re.findall("[DEHKNQR][DEHKNQR]",sequence[1:190])); 
		return nls;
	
	def getMonoNLSSignal24(self, sequence):
		nls = len(re.findall("[KRH][KRH]",sequence[1:190]))*2; 
		nls += len(re.findall("[DE][DE]",sequence[1:190]))*-2; 
		nls += len(re.findall("[NQ][NQ]",sequence[1:190]))*1; 
		nls += len(re.findall("[DE][NQ]",sequence[1:190]))*-1; 
		nls += len(re.findall("[NQ][DE]",sequence[1:190]))*-1; 
		nls += len(re.findall("[KRH][NQ]",sequence[1:190]))*1; 
		nls += len(re.findall("[NQ][KRH]",sequence[1:190]))*1; 
		return nls;
		
	def getMonoNLSSignal25(self, sequence):
		nls = len(re.findall("[DEHKR][DEHKR]",sequence[1:190])); 
		return nls;
		
	# a good NLS distribution
	def getMonoNLSSignal26(self, sequence):	
		# find polar tuples with first 190 residues
		nls = len(re.findall("[DEHKNQR][DEHKNQR]",sequence[1:190])); 
		
		# non-hydrophobic n-terminus
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880, "X" : 0.5};
		sum = 0;
		for aa in sequence[1:10]:
			if aa in hydr.keys():
				sum += 1-hydr[aa];
		nls += sum*4;
		
		# occurences of NLS Signal (weights learned via linear regression)
		nls_mono = len(re.findall("K[KR]",sequence[1:100]))*0.2;
		nls_mono += len(re.findall("K[KR][KR]",sequence[1:100]))*0.43;
		nls_mono += len(re.findall("[KR]{4}",sequence))*0.83;
		nls_mono += len(re.findall("[KR]{6}",sequence))*-0.76;
		nls_mono += len(re.findall("[KR]{3}[HP]",sequence))*0.96;
		nls_mono += len(re.findall("P.K{3}",sequence))*0.7;
		nls_mono += self.getNLSDBSignals(sequence)*0.55;
			
		nls += nls_mono*15;
		return nls;
	
	# a hard NLS Bound
	def getMonoNLSSignal27(self, sequence):
		nls = self.getNLSDBSignals(sequence);
		nls += self.getMonoNLSSignal2(sequence);
		nls += self.getMonoNLSSignal3(sequence);
		nls += self.getMonoNLSSignal11(sequence);
		nls += self.getMonoNLSSignal12(sequence);
		nls15 = self.getMonoNLSSignal15(sequence);
		if nls15 > 1:
			nls += nls15 - 1;
		nls += self.getMonoNLSSignal16(sequence);
		
		return nls;
  
	def getBipartiteNLSSignal(self, sequence):
		nls_bi = 0;
		nls_bi += len(re.findall("[KR]{2,}.{8,14}[RK]{2,}",sequence)); # Ema02
		nls_bi += len(re.findall("G{20,40}",sequence)); #Nakai
		count = 0;
		for i in range(0,9):
			if sequence[i]=='R':
				count += 1;
		if count>=6:
			nls_bi += 1 #Matsuda
		return nls_bi;
			
	def getBipartiteNLSSignal1(self, sequence):
		nls_bi = 0;
		nls_bi += len(re.findall("[KR]{2,}.{8,14}[RK]{2,}",sequence)); # Ema02
		return nls_bi;
	
	def getBipartiteNLSSignal2(self, sequence):
		nls_bi = 0;
		nls_bi += len(re.findall("G{20,40}",sequence)); #Nakai
		return nls_bi;
		
	def getBipartiteNLSSignal3(self, sequence):
		nls_bi = 0;
		count = 0;
		for i in range(0,9):
			if sequence[i]=='R':
				count += 1;
		if count>=6:
			nls_bi += 1 #Matsuda
		return nls_bi;

	def getNESSignal(self, sequence):
		nes = re.findall("L[^DERKH]{2,3}[HWYFM][^DERKH]{2,3}L[^DERKH][LI]",sequence); #Nakai
		return len(nes);
		
	def getERSignal(self, sequence):
		# KDEL
		er_target = len(re.findall("KDEL$",sequence))*20;
		er_target += len(re.findall("[KRHQSA][DENQ]EL$",sequence))*2;
		er_target += len(re.findall("DEL$",sequence))*4;
		er_target += len(re.findall("EL$",sequence))*0.5;
		# Dilysine motif
		#er_target += len(re.findall("KK..$",sequence));
		#er_target += len(re.findall("K.K..$",sequence));
		# other N-terminal sequence
		#er_target += len(re.findall("^.{1,3}[KR]{2}",sequence));
		#er_target += len(re.findall("^.{1,2}[KR].[KR]",sequence));
		return er_target;
	
	def getERSignal2(self, sequence):
		# Dilysine motif
		er_target = len(re.findall("KK..$",sequence));
		er_target += len(re.findall("K.K..$",sequence));
		# other N-terminal sequence
		er_target += len(re.findall("^.{1,3}[KR]{2}",sequence));
		er_target += len(re.findall("^.{1,2}[KR].[KR]",sequence));
		return er_target;
	
	def getERSignal3(self, sequence):
		score = 0;
		kdel = [ "K", "D", "E" ,"L"];
		for aa in sequence[-10:]:
			if aa in kdel:
				score +=1;
		return score;
		
	def getERSignal4(self, sequence): 
		score = 1;
		#  weak ER signal
		er = self.getERSignal3(sequence)
		score *= er;
		# definitly  strong ER signal
		er_strong = self.getERSignal(sequence)
		if (er_strong > 0):
			score = 11;
		
		return score;
	
	def getERSignal5(self, sequence):
		score = 1;
		# glycostie  -> should be less then 12
		glyco = self.getGlycoSiteSignal(sequence);
		if glyco <= 10:
			score += 3;
		# transmembrane -> should be less than 9
		trans = self.getTransMembraneSignal2(sequence);
		if trans < 8.5:
			score += 5;
		elif trans < 7:
			score += 2;
		# weak ER signal
		er = self.getERSignal3(sequence)
		score *= er;
		# definitly  strong ER signal
		er_strong = self.getERSignal(sequence)
		if (er_strong > 0):
			score = 100;
		return score;

	
	def getPeroxiSignal(self, sequence):
		# scoring by Nakai
		score = 0;
		score += len(re.findall("SKL$",sequence)) * 0.83333;
		score += len(re.findall("SKF$",sequence)) * 0.5;
		score += len(re.findall("[SAGCN][RKH][LIVMAF]$",sequence)) * 0.25;
		return score; 

	def getPeroxiSignal2(self, sequence):
		score = 0;
		score += len(re.findall("SKL$",sequence))*6;
		score += len(re.findall("SKF$",sequence))*4;
		score += len(re.findall("[SAGCN][RKH][LIVMAF]$",sequence));
		#PTS1
		score += len(re.findall("[SA][QN]L$$",sequence))*2;
		score += len(re.findall("[SAC][KHR]L$$",sequence));
		#PTS2
		score += len(re.findall("[RK][LVI].{4,5}[HQ][LA]",sequence[1:50]))*3;
		
		return score; 

	def getPeroxiSignal3(self, sequence):
		score = 0;
		score += len(re.findall("[SAGCN][RKH][LIVMAF]$",sequence));
		#PTS1
		score += len(re.findall("[SA][QN]L$$",sequence));
		score += len(re.findall("[SAC][KHR]L$$",sequence));
		#PTS2
		score += len(re.findall("[RK][LVI].{5}[HQ][LA]",sequence[1:50]));
		return score; 
	
	def getPeroxiSignal4(self, sequence):
		score = 0;
		# Last is not E,C,W
		last = [ "E" , "C", "W" ];
		if sequence[-1] in last:
			score += -10;
		
		# pattern P 1
		pos1 = { "L" : 30, "M" : 6, "I" : 6};
		pos2 = { "K" : 20, "R" : 10, "A" : 8};
		pos3 = { "S" : 10, "A" : 6, "F" : 3};
		
		if sequence[-1] in pos1.keys():
			score += pos1[sequence[-1]];
		if sequence[-2] in pos2.keys():
			score += pos2[sequence[-2]];
		if sequence[-3] in pos3.keys():
			score += pos3[sequence[-3]];
				
		return score; 

	def getPeroxiSignal5(self, sequence):
		score = 0;
		posneg = [ "ER", "RE" ,"KE","EK", "QR","RQ","KQ","QK","DK","KD","DR", "RD"]
		
		for elem in posneg:
			if elem in sequence[-10:]:
				score += 10;				
		return score; 

	def getPeroxiSignal6(self, sequence):
		score = 0;
		posneg = [ "ER", "RE" ,"KE","EK", "QR","RQ","KQ","QK","DK","KD","DR", "RD"]
		
		for elem in posneg:
			if elem in sequence[-20:]:
				score += 10;				
		return score; 

	def getPeroxiSignal7(self, sequence):
		score = 0;
		posneg = [ "ER", "RE" ,"KE","EK", "QR","RQ","KQ","QK","DK","KD","DR", "RD"]
		
		for elem in posneg:
			if elem in sequence:
				score += 10;				
		return score; 
		
	def getPeroxiSignal8(self, sequence):
		score = 0;
		coeff = { 'F' : 1.384, 'W': 1.905, 'N' : -0.462, 'E' : -0.177, 'C' : -0.585, 'Y': 0.457,
				'H' : 0.521, 'G' : 0.351, 'T' : 0.226, 'P' : -0.106 };
		
		for aa in sequence:
			if aa in coeff.keys():
				score += coeff[aa];				
		return score; 

	def getPeroxiSignal9(self, sequence):
		score = len(re.findall("[SA][KRH]L$",sequence));
		score -= len(re.findall("SKL.",sequence));
		return score
		
		
	def getPeroxiSignalPlant(self, sequence):
		# McNew & Goodman, "The targeting and assembly of peroxisomal proteins: some old rules do not apply.") Trends Biochem Sci., 21(2):54-8, 1996). 
		score = 0;
		score += len(re.findall("[RK][LI].{5}[HQ]L",sequence));
		return score; 	
		
	def getVacuoleSignal(self, sequence):
		score = 0;
		# PSORT II Signal
		score += len(re.findall("[TIK]LP[NKI]",sequence))*5;
		# QRPL from yeast (Valls et al 1994)
		score += len(re.findall("QRPL",sequence[1:200]))*6;
		# weak QRPL signal from Voorst - Mutational Analysis of the Vacuolar Sorting Signal of Procarboxypeptidase Y in Yeast Shows a Low Requirement for Sequence Conservation
		score += len(re.findall("Q[RLVKQCSNGE][PFVLATSQ][LFVI]",sequence[1:200]));
		#Dileucline related
		score += len(re.findall("[DE]E.{3}L[ILV]",sequence))*2;
		# non-hydrophobic tail
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880, "X" : 0.5};
		sum = 0;
		for aa in sequence[1:10]:
			if aa in hydr.keys():
				sum += hydr[aa];
		score += sum/2
		
		return score; 

	def getVacuoleSignal2(self, sequence):
		score = 0;
		#Dileucline related
		score += len(re.findall("[DE]E.{3}L[ILV]",sequence))*2;		
		return score; 
		
	def getVacuoleSignal3(self, sequence):
		score = 0;
		# PSORT II Signal
		score += len(re.findall("[TIK]LP[NKI]",sequence));
		return score; 
		
	def getVacuoleSignal4(self, sequence):
		score = 0;
		score += len(re.findall("QRPL",sequence[1:200]))*6;
		return score; 

	def getVacuoleSignal5(self, sequence):
		score = 0;
		# weak QRPL signal from Voorst - Mutational Analysis of the Vacuolar Sorting Signal of Procarboxypeptidase Y in Yeast Shows a Low Requirement for Sequence Conservation
		score += len(re.findall("Q[RLVKQCSNGE][PFVLATSQ][LFVI]",sequence[1:200]));
		return score; 
		
	def getVacuoleSignal6(self, sequence):
		score = 0;
		# non-hydrophobic tail
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880, "X" : 0.5};
		sum = 0;
		for aa in sequence[1:10]:
			if aa in hydr.keys():
				sum += 1-hydr[aa];
		score += sum;
		return score;
		
	def getVacuoleSignal7(self, sequence):
		score = 0;
		# non-hydrophobic tail
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880, "X" : 0.5};
		sum = 0;
		for aa in sequence[2:5]:
			if aa in hydr.keys():
				sum += 1-hydr[aa];
		score += sum;
		return score;
		
	def getVacuoleSignal8(self, sequence):
		score = 0;
		# non-hydrophobic tail
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880, "X" : 0.5};
		sum = 0;
		for aa in sequence[2:8]:
			if aa in hydr.keys():
				sum += 1-hydr[aa];
		score += sum;
		return score;
		
	def getVacuoleSignal9(self, sequence):
		score = 0;
		# non-hydrophobic tail
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880};
		sum = 1;
		for aa in sequence[2:10]:
			if aa in hydr.keys():
				sum *= (1-hydr[aa]);
		score += sum;
		return score;

	def getVacuoleSignal10(self, sequence):
		score = 0;
		score += len(re.findall("[QN][KR]P[LI]",sequence[1:200]));
		return score; 
		
	def getVacuoleSignal11(self, sequence):
		score = 0;
		score += len(re.findall("[QN][KR]P[LI]",sequence));
		return score;
		 
	def getVacuoleSignal12(self, sequence):
		score = 0;
		score += len(re.findall("QRPL",sequence));
		return score; 
		
	def getVacuoleSignal13(self, sequence):
		score = 0;
		score += len(re.findall("Q.PL",sequence[1:200]));
		score += len(re.findall("QR.L",sequence[1:200]));
		score += len(re.findall("QRP",sequence[1:200]));
		score += len(re.findall("RPL",sequence[1:200]));
		return score; 

	def getVacuoleSignal14(self, sequence):
		score = 0;
		# weak QRPL signal from Voorst - Mutational Analysis of the Vacuolar Sorting Signal of Procarboxypeptidase Y in Yeast Shows a Low Requirement for Sequence Conservation
		score += len(re.findall("Q[RLVKQCSNGE][PFVLATSQ][LFVI]",sequence[1:100]));
		return score; 

	def getVacuoleSignal15(self, sequence):
		score = 0;
		# weak QRPL signal from Voorst - Mutational Analysis of the Vacuolar Sorting Signal of Procarboxypeptidase Y in Yeast Shows a Low Requirement for Sequence Conservation
		score += len(re.findall("Q[RLVKQCSNGE][PFVLATSQ][LFVI]",sequence[1:50]));
		return score; 
	
	def getVacuoleSignal16(self, sequence):
		score = 0;
		# weak QRPL signal from Voorst - Mutational Analysis of the Vacuolar Sorting Signal of Procarboxypeptidase Y in Yeast Shows a Low Requirement for Sequence Conservation
		score += len(re.findall("Q[RLVKQCSNGE][PFVLATSQ][LFVI]",sequence[1:20]));
		return score; 

	def getVacuoleSignal17(self, sequence):
		score = 0;
		# weak QRPL signal from Voorst - Mutational Analysis of the Vacuolar Sorting Signal of Procarboxypeptidase Y in Yeast Shows a Low Requirement for Sequence Conservation
		score += len(re.findall("Q[RLVKQCSNGE][PFVLATSQ][LFVI]",sequence[1:10]));
		return score; 

	def getVacuoleSignal18(self, sequence):
		score = 0;
		# non-hydrophobic tail
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880, "X" : 0.5};
		sum = 0;
		for aa in sequence[2:5]:
			if aa in hydr.keys():
				sum += 1-hydr[aa];
		score += sum;
		
		# QRPL from yeast (Valls et al 1994)
		score = (len(re.findall("QRPL",sequence[1:200]))+1) * score;
		return score;
		
	def getVacuoleSignal19(self, sequence):
		score = 0;
		# Sebs score
		signals = { "QR" : 2, "LL" : 1.5, "LLL" : 2.5, "LLLL" : 5, "LLLLL" : 12, "LLLLLL" : 20};
		
		for key in signals.keys():
			nr = len(re.findall(key,sequence[1:50]));
			score += nr * signals[key];	
		
		return score;

	def getVacuoleSignal20(self, sequence):
		score = 0;
		# Sebs score
		signals = { "QQQ" : 2, "LL" : 1.5, "LLL" : 2.5, "LLLL" : 5, "LLLLL" : 12, "LLLLLL" : 20};
		
		for key in signals.keys():
			nr = len(re.findall(key,sequence[1:50]));
			score += nr * signals[key];	
		
		return score;

	def getVacuoleSignal21(self, sequence):
		score = 1;
		# Sebs score
		signals = { "QR" : 2, "LL" : 1.5, "LLL" : 2.5, "LLLL" : 5, "LLLLL" : 12, "LLLLLL" : 20};
		
		for key in signals.keys():
			nr = len(re.findall(key,sequence[1:50]))+1;
			score *= nr * signals[key];	
		
		return score;

	def getVacuoleSignal22(self, sequence):
		score = 0;
		# Sebs score
		signals = { "QQQ" : 2, "QQQQ" : 4, "LL" : 1.5, "LLL" : 2.5, "LLLL" : 5, "LLLLL" : 12, "LLLLLL" : 20};
		
		for key in signals.keys():
			nr = len(re.findall(key,sequence[1:50]));
			score += nr * signals[key];	
		
		return score;
		
	def getVacuoleSignal23(self, sequence):
		score = 0;
		# Sebs score
		signals = { "EEE" : 2, "EEEE" : 4, "EEEEE" : 7.5, "LL" : 1.5, "LLL" : 2.5, "LLLL" : 5, "LLLLL" : 12, "LLLLLL" : 20};
		
		for key in signals.keys():
			nr = len(re.findall(key,sequence[1:50]));
			score += nr * signals[key];	
		
		return score;

	def getVacuoleSignal24(self, sequence):
		score = 0;
		# Sebs score
		signals = {  "EEE" : 2, "EEEE" : 4, "EEEEE" : 7.5, "QQQ" : 2, "QQQQ" : 4, "LL" : 1.5, "LLL" : 2.5, "LLLL" : 5, "LLLLL" : 12, "LLLLLL" : 20};
		
		for key in signals.keys():
			nr = len(re.findall(key,sequence[1:50]));
			score += nr * signals[key];	
		
		return score;

	def getVacuoleSignal25(self, sequence):
		score = 0;
		# Sebs score
		signals = {  "EEE" : 2, "EEEE" : 4, "EEEEE" : 7.5, "QQQ" : 2, "QQQQ" : 4, "LL" : 1.5, "LLL" : 2.5, "LLLL" : 5, "LLLLL" : 12, "LLLLLL" : 20};
		
		for key in signals.keys():
			nr = len(re.findall(key,sequence[1:50]));
			score += nr * signals[key];	
			nr = len(re.findall(key,sequence[1:100]));
			score += nr * signals[key] * 0.5;
			nr = len(re.findall(key,sequence));
			score += nr * signals[key] * 0.1;		
		
		return score;

	def getVacuoleSignal26(self, sequence): 
		score = 0;
		score += len(re.findall("L",sequence[1:10]));
		return score;

	def getVacuoleSignal27(self, sequence): 
		score = 0;
		score += len(re.findall("LL",sequence[1:10]))*2;
		score += len(re.findall("LL",sequence[1:10]));
		return score;

	def getVacuoleSignal28(self, sequence): 
		score = 0;
		score += len(re.findall("LL",sequence[1:10]));
		return score;

	# sebs final?
	def getVacuoleSignal29(self, sequence): 
		score = 0;
		# cleavage site within 20-40 (less frequent in 0-20 or 40-60)
		cs = self.getCleavageSite(sequence);
		if (cs >=20 & cs <=40):
			score += 30;
		if (cs < 20 or (cs>40 & cs <=60)):
			score += 10;
		# no weak ER signal
		er = self.getERSignal3(sequence)
		score *= (10-er);
		# definitly no strong ER signal
		er_strong = self.getERSignal(sequence)
		if (er_strong > 0):
			score = 0;
		
		return score;

	def getVacuoleSignal30(self, sequence): 
		score = 0;
		# no plasma & extracellular (they have only few QQ and QQQs)
		qqq = len(re.findall("QQQ",sequence))*20;
		qq = len(re.findall("QQ",sequence))*4;
		score += qqq+qq;
		return score;

	def getVacuoleSignal31(self, sequence):
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880};
		cs = self.getCleavageSite(sequence);
		score = 0;
		for aa in sequence[cs-10:cs]:
			if aa in hydr.keys():
				score += hydr[aa];
		return score;
		
	def getVacuoleSignal32(self, sequence):
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880};
		cs = self.getCleavageSite(sequence);
		score = 0;
		for aa in sequence[cs:cs+10]:
			if aa in hydr.keys():
				score += hydr[aa];
		return score;
		
	def getVacuoleSignal33(self, sequence):
		score = 0;
		score = len(re.findall("[CW]",sequence[1:20]));
		score += len(re.findall("C",sequence[1:100]));
		return score;
		
	def getVacuoleSignal34(self, sequence):
		score = 0;
		score = len(re.findall("^M[MK]",sequence));
		return score;
		
	def getVacuoleSignal35(self, sequence):
		score = 0;
		before = '';
		stretch = 0;
		sum = 0;
		for aa in sequence:
			if aa == before:
				stretch += 1;
			else:
				sum += stretch;
				stretch = 0;
			before = aa;
		score = float(sum) / len(sequence);
			
		return score;

	def getVacuoleSignal36(self, sequence):
		score = 0;
		score = len(re.findall("^M[MK]",sequence))*2;
		c = len(re.findall("[CW]",sequence[1:20]));
		c += len(re.findall("C",sequence[1:100]));
		if c <= 2:
			score +=2
		score -= len(re.findall("VD",sequence[1:100])) *2;
		score -= len(re.findall("[ECW]$",sequence))*3;
			
		return score;
	
	def getVacuoleSignal37(self, sequence):
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880};
		cs = self.getCleavageSite(sequence);
		score = 0;
		for aa in sequence[1:20]:
			if aa in hydr.keys():
				score += hydr[aa];
		return score;
	
	def getLysosomalSignal1(self, sequence):
		score = len(re.findall("N.[ST]",sequence));
		if score > 2:
			score = 1;
		else:
			score = 0;
		return score;
		
	def getLysosomalSignal2(self, sequence):
		score = len(re.findall("[GF]L",sequence[1:20]));
		return score;
		
	def getLysosomalSignal3(self, sequence):
		score = len(re.findall("P",sequence));
		return score;
		
	def getLysosomalSignal4(self, sequence):
		score = len(re.findall("GY",sequence[-100:]));
		return score;
		
	def getLysosomalSignal5(self, sequence):
		score = len(re.findall("[^ECW]$",sequence));
		return score;
	
	def getLysosomalSignal6(self,sequence):
		score = self.getGlycoSiteSignal(sequence);
		if score > 1.5:
			score = 2;
		else:
			score = 0;
			
		return score;
	
	def getLysosomalSignal7(self,sequence):
		score = 0;
		score  += self.getLysosomalSignal1(sequence)
		p = self.getLysosomalSignal3(sequence)
		if p>=8 and p<=40:
			score += 3;
		if p<=16:
			score -=1;
		score += self.getLysosomalSignal5(sequence) 
		score += self.getLysosomalSignal6(sequence)
		return score;
		
	def getTransGolgiSignal(self, sequence):
		# Humphrey et al., "Localization of TGN38 to the trans-Golgi network: involvement of a cytoplasmic tyrosine-containing sequence", The Journal of Cell Biology, 120:1123-35, 1993
		score = 0;
		score += len(re.findall("YQRL",sequence));
		return score; 		
		
	def getGolgiSignal(self, sequence):
		# Sun - Sterol-regulated transport of SREBPs from
		score = len(re.findall("MELADL",sequence));
		# GRIP domain obtained from SMART
		score += len(re.findall("[ETRQAD]Y[LATI][KR][NHK][VI][IVLMF]",sequence));
		# part 2
		score += len(re.findall("[LIAVT].[VAL][VIL]..[IVLAM][LA]",sequence));
		return score;
		
	def getGolgiSignal2(self, sequence):
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880};
		flank=2;
		window=20;
		max = -10000;
		sum = 0;
		for pos in range(1,len(sequence)-2*flank-window):
					sum = 0;
					for aa in sequence[pos:pos+flank]:
						if aa in hydr.keys():
							sum += -1* (hydr[aa]-0.5);
					for aa in sequence[pos+flank:pos+flank+window]:
						if aa in hydr.keys():
							sum += hydr[aa]-0.5;
					for aa in sequence[pos+flank+window:pos+2*flank+window]:
						if aa in hydr.keys():
							sum += -1* (hydr[aa]-0.5);
					if sum > max:
						max = sum;

		score = sum;
		return score;

	def getGolgiSignal3(self, sequence):
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880};
		flank=2;
		window=20;
		max = -10000;
		for pos in range(len(sequence)-150,len(sequence)-2*flank-window):
					sum = 0;
					for aa in sequence[pos:pos+flank]:
						if aa in hydr.keys():
							sum += -1* (hydr[aa]-0.5);
					for aa in sequence[pos+flank:pos+flank+window]:
						if aa in hydr.keys():
							sum += hydr[aa]-0.5;
					for aa in sequence[pos+flank+window:pos+2*flank+window]:
						if aa in hydr.keys():
							sum += -1* (hydr[aa]-0.5);
					if sum > max:
						max = sum;

		score = max;
		return score;

	def getGolgiSignal4(self, sequence):
		hydrnopolar = [ "I", "L", "V","G","A","M","F"];
		basic = [ "H", "K" ,"R"];
		flank=2;
		window=20;
		max = 0;
		for pos in range(len(sequence)-150,len(sequence)-2*flank-window):
					sum = 0;
					for aa in sequence[pos:pos+flank]:
						if aa in basic:
							sum += 1
					for aa in sequence[pos+flank:pos+flank+window]:
						if aa in hydrnopolar:
							sum += 1;
					for aa in sequence[pos+flank+window:pos+2*flank+window]:
						if aa in basic:
							sum += 1;
					if sum > max:
						max = sum;

		score = max;
		return score;

	def getGolgiSignal5(self, sequence):
		nopolar = [ "I", "L", "V","G","A","M","F","P"];
		windows = [ 6,8,10,12];
		max = 0;
		for w in windows:
			max_w = 0;
			for pos in range(0,50):
				sum = 0;
				for aa in sequence[pos:pos+w]:
					if aa in nopolar:
						sum += 1
				if sum > max_w:
					max_w = sum;
			if float(max_w)/float(w) > max:
				max = float(max_w)/float(w);
		
				

		score = max		
		return score;
	
	def getGolgiSignal6(self, sequence):
		return self.getGolgiSignal4(sequence)/20 + self.getGolgiSignal5(sequence);

	def getGolgiSignal7(self, sequence):
		return self.getGolgiSignal4(sequence)/20 * self.getGolgiSignal5(sequence);
		
	def getGolgiSignal8(self, sequence):
		score = 0;
		# glycostie  -> should be less then 12
		glyco = self.getGlycoSiteSignal(sequence);
		if glyco <= 12:
			score += 30;
		# transmembrane -> should be less than 9
		trans = self.getTransMembraneSignal2(sequence);
		if trans < 9:
			score += 30;
		elif trans < 7.5:
			score += 15;
		# cleavage site within 20-40 (less frequent in 0-20 or 40-60)
		cs = self.getCleavageSite(sequence);
		if (cs >=20 & cs <=50):
			score += 30;
		# no weak ER signal
		er = self.getERSignal3(sequence)
		score *= (10-er);
		# definitly no strong ER signal
		er_strong = self.getERSignal(sequence)
		if (er_strong > 0):
			score = 0;
		return score;
	
	def getGolgiSignal9(self, sequence):
		score = 0;
		
		if len(re.findall("C",sequence[1:40])) <= 2:
			score += 1;
		score += len(re.findall("W.W",sequence));
		
		g4 = self.getGolgiSignal4(sequence)
		if g4 <=11 and g4 >= 15:
			score += 2;
			
		p8 = self.getPeroxiSignal8(sequence)
		if p8 >=26 and p8 <=71:
			score +=3;
		
		# no weak ER signal
		er = self.getERSignal3(sequence)
		if er >= 5:
			score -= er;
		# definitly no strong ER signal
		er_strong = self.getERSignal(sequence)
		if (er_strong > 0):
			score -= 10;
		
		g3 = self.getGolgiSignal3(sequence)
		if g3 >= 5:
			score -= 3;
		
		score -= len(re.findall("[ECW]$",sequence));
		
		p = len(re.findall("P",sequence));
		if p <= 8 and p>=40:
			score +=2;
		
		return score;
		
	def getGolgiSignal10(self, sequence):
		score = 0;
		
		if len(re.findall("C",sequence[1:40])) <= 2:
			score += 1;
		
		g4 = self.getGolgiSignal4(sequence)
		if g4 <=11 and g4 >= 15:
			score += 2;
			
		p8 = self.getPeroxiSignal8(sequence)
		if p8 >=26 and p8 <=71:
			score +=3;
		
		# no weak ER signal
		er = self.getERSignal3(sequence)
		if er >= 5:
			score -= er;
		# definitly no strong ER signal
		er_strong = self.getERSignal(sequence)
		if (er_strong > 0):
			score -= 10;
		
		g3 = self.getGolgiSignal3(sequence)
		if g3 >= 5:
			score -= 3;
		
		score -= len(re.findall("[ECW]$",sequence));
		
		p = len(re.findall("P",sequence));
		if p <= 8 and p>=40:
			score +=2;
		
		return score;

	def getGolgiSignal11(self, sequence):
		score = 0;
		
		if len(re.findall("C",sequence[1:40])) <= 2:
			score += 1;
		score += len(re.findall("W.W",sequence));
		
		glyco = self.getGlycoSiteSignal(sequence);
		if glyco <= 12 and glyco > 1.5:
			score += 3;
		
		g4 = self.getGolgiSignal4(sequence)
		if g4 <11 or g4 > 15:
			score -= 2;
			
		p8 = self.getPeroxiSignal8(sequence)
		if p8 <26 or p8 >71:
			score -=3;
		
		# no weak ER signal
		er = self.getERSignal3(sequence)
		if er >= 5:
			score -= er;
		# definitly no strong ER signal
		er_strong = self.getERSignal(sequence)
		if (er_strong > 0):
			score -= 10;
		
		g3 = self.getGolgiSignal3(sequence)
		if g3 >= 5:
			score -= 3;
		
		score -= len(re.findall("[ECW]$",sequence));
		
		p = len(re.findall("P",sequence));
		if p <= 8 and p>=40:
			score +=2;
		
		return score;

	def getGolgiSignal12(self, sequence):
		score = 0;
		
		if len(re.findall("C",sequence[1:40])) <= 2:
			score += 1;
		score += len(re.findall("W.W",sequence));		
		# no weak ER signal
		er = self.getERSignal3(sequence)
		if er >= 5:
			score -= er;
		# definitly no strong ER signal
		er_strong = self.getERSignal(sequence)
		if (er_strong > 0):
			score -= 10;
	
		score -= len(re.findall("[ECW]$",sequence));
		
		p = len(re.findall("P",sequence));
		if p <= 8 and p>=40:
			score +=2;
		
		return score;
		
	def getTransMembraneSignal(self, sequence):
		score = 0;
		# Di-leucine motifs
		score += len(re.findall("DD.{2}LL",sequence)) * 5;
		# Di-leucine motifs weaker
		score += len(re.findall("DD.{2}[LMV]{2}",sequence));
		return score;
	
	def getTransMembraneSignal2(self, sequence):
		score = 0;
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880};
		# find longest hydrophobic part
		max = 0;
		pos = 0;
		while (pos < len(sequence)):
			sum = 0;
			length = 0;
			wrong = 0;
			while length < len(sequence)-pos :
				if sequence[pos+length] in hydr.keys():
					value = hydr[sequence[pos+length]] - 0.5;
					sum += value;
					if value > 0:
						wrong = 0;
					else:
						wrong += 1;
				
				length += 1;
				if sum > max:
					max = sum;
				if sum < 0 or wrong >= 2:
					break;
					
			pos += length;
					
		score = max;
		return score;
	
	def getTransMembraneSignal3(self, sequence):
		score = 0;
		g = self.getGolgiSignal(sequence);
		if (g>=2):
			score += 5;
		g = self.getGolgiSignal3(sequence);
		if (g>=5):
			score += 5;
		g = self.getGolgiSignal6(sequence);
		if (g>=1.5):
			score += 5;
		g = self.getGolgiSignal7(sequence);
		if (g>=0.5):
			score += 5;
		if len(sequence) > 800:
			score += 5;
		g = self.getGlycoSiteSignal(sequence)
		if g >= 1.5:
			score += 5;
		#glycosite signal
		g = len(re.findall("N.[ST]",sequence));
		if g > 2:
			score +=5;
		score += self.getTransMembraneSignal2(sequence)*2;
		
		return score;
		
	
	def getEndocytosisSignal(self, sequence):
		# need tm helix number!!
		score = 0;
		score += len(re.findall("LL",sequence));
		# Di-leucine motifs
		score += len(re.findall("DD.{2}LL",sequence));
		# all from Bonifacino - IGNALS FOR SORTING OF TRANSMEMBRANE PROTEINS TO ENDOSOMES AND LYSOSOMES
		score += len(re.findall("[DE].{3}L[LI]",sequence))*2;
		score += len(re.findall("NP.[LY]",sequence))*2;
		score += len(re.findall("NFF.{1,2}D",sequence))*2;
		score += len(re.findall("Y.{2}[IFVLML]",sequence));
		score += len(re.findall("DP[FW]",sequence));
		score += len(re.findall("F.D.F",sequence));
		score += len(re.findall("NPF",sequence));
		score += len(re.findall("L[IFVLML].[IFVLML][DE]",sequence))*3;
		score += len(re.findall("LLDLL",sequence))*3;
		score += len(re.findall("PWDLW",sequence))*3;
		return score; 		
	
	def getGlycoSiteSignal(self, sequence):
		glyco_sites = 1000*(float(len(re.findall("N[^P][S|T][^P]",sequence))) /float( len(sequence)))
		return glyco_sites;
		
	def getAnchorSignal(self, sequence):
		score = 0;
    # Myristoylation
		score += len(re.findall("^MG[^DERKH]{3}[ACGSTPN]",sequence));
		# Palmitoylation
		score += len(re.findall("[IFVLML]LCC[^DERKH][RK][RK]",sequence));
		score += len(re.findall("IPCCPV",sequence));
		# Prenylation
		score += len(re.findall("C[ILV][ILV][LSMCAQ]$",sequence)); #CaaX motif
		score += len(re.findall("[LSMCAQ][LSMCAQ]CC$",sequence));
		score += len(re.findall("[LSMCAQ]CC[LSMCAQ]$",sequence));
		score += len(re.findall("CC[LSMCAQ][LSMCAQ]$",sequence));
		score += len(re.findall("[LSMCAQ]C[LSMCAQ]C$",sequence));
		
		return score; 

	def getMIPSignal(self, sequence):
		score = 0;
		# MIP signal [Ema02]
		score += len(re.findall("R.[FLI].{2}[TSG].{4}",sequence));
		# Ema02 aa histogram, lots R,A,S  but rare D,E
		pos = ['R','A','S'];
		neg = ['D','E'];
		signal = 0;
		for aa in sequence[1:70]:
			if aa in pos:
				signal += 1;
			if aa in neg:
				signal -= 1;
		signal = (score+1) * signal;
		return score;
	
	def getMIPSignal2(self, sequence):
		score = 0;
		charge = {"A" : 6.107, "C" : 5.02, "D" : 2.98, "E" : 3.08, "F" : 5.91,
				"G" : 6.064, "H" : 7.64, "I" : 6.038, "K" : 9.47, "L" : 6.038,
				"M" : 5.74, "N" : 5.41, "P" : 6.3, "Q" : 5.65, "R" : 10.76,
				"S" : 5.68, "T" : 5.64, "V" : 6.002, "W" : 5.88, "Y" : 5.63};
		charge_sum = 0;
		for aa in sequence[1:20]:
			if aa in charge.keys():
				charge_sum += charge[aa]
		if charge_sum > 120:
			score += 10;
		if charge_sum > 125:
			score += 30;
		
		s_sum = 0;
		for aa in sequence[1:50]:
			if aa == 'S':
				s_sum += 1;
		if s_sum < 7:
			score += 5;
		if s_sum < 4:
			score += 15;
			
		return score;

	def getMIPSignal3(self, sequence):
		score = 0;
		charge = {"A" : 6.107, "C" : 5.02, "D" : 2.98, "E" : 3.08, "F" : 5.91,
				"G" : 6.064, "H" : 7.64, "I" : 6.038, "K" : 9.47, "L" : 6.038,
				"M" : 5.74, "N" : 5.41, "P" : 6.3, "Q" : 5.65, "R" : 10.76,
				"S" : 5.68, "T" : 5.64, "V" : 6.002, "W" : 5.88, "Y" : 5.63};
		charge_sum = 0;
		for aa in sequence[1:20]:
			if aa in charge.keys():
				charge_sum += charge[aa]
		score += (charge_sum-120);
		score += (charge_sum-125)*3;
		
		s_sum = 0;
		for aa in sequence[1:50]:
			if aa == 'S':
				s_sum += 1;
				
		score += (7-s_sum)*3;
		score += (4-s_sum)*5;
			
		return score;

	def getMIPSignal4(self, sequence):
		score = 0;
		charge = {"A" : 6.107, "C" : 5.02, "D" : 2.98, "E" : 3.08, "F" : 5.91,
				"G" : 6.064, "H" : 7.64, "I" : 6.038, "K" : 9.47, "L" : 6.038,
				"M" : 5.74, "N" : 5.41, "P" : 6.3, "Q" : 5.65, "R" : 10.76,
				"S" : 5.68, "T" : 5.64, "V" : 6.002, "W" : 5.88, "Y" : 5.63};
		charge_sum = 0;
		for aa in sequence[1:20]:
			if aa in charge.keys():
				charge_sum += charge[aa]
		score += (charge_sum-120);
		score += (charge_sum-125)*3;
		
		s_sum = 0;
		for aa in sequence[1:50]:
			if aa == 'S':
				s_sum += 1;
				
		score += (7-s_sum)*3;
		score += (4-s_sum)*5;
			
		score -= len(re.findall("SSS",sequence))*10;
			
		return score;
		
	def getMIPSignal5(self, sequence):
		score = self.getAliphaticHelix(sequence[1:100]);
		return score;
		
	def getMIPSignal6(self, sequence):
		score = self.getAliphaticHelix(sequence[1:200]);
		return score;
		
	def getMIPSignal7(self, sequence):
		score = self.getAliphaticHelix(sequence[-100:]);
		return score;

	def getMIPSignal8(self, sequence):
		sequence = sequence[1:200];
		# Eisenberg hydrophobic moment -> modified to basic moment
		# Eisenberg paper Hydrophobic Moments and Protein structure Hydrophobicity Values:
		charge = {"A" : 6.107, "C" : 5.02, "D" : 2.98, "E" : 3.08, "F" : 5.91,
				"G" : 6.064, "H" : 7.64, "I" : 6.038, "K" : 9.47, "L" : 6.038,
				"M" : 5.74, "N" : 5.41, "P" : 6.3, "Q" : 5.65, "R" : 10.76,
				"S" : 5.68, "T" : 5.64, "V" : 6.002, "W" : 5.88, "Y" : 5.63};
		max_pattern = 0;
		for pos in range(0,len(sequence)-18):
			sum_cos = 0;
			sum_sin = 0;
			aa_index = 1;
			for aa in sequence[pos:pos+18]:
				if aa in charge.keys():
					sum_cos += (charge[aa]-5)*cos(1.745329252*aa_index);
					sum_sin += (charge[aa]-5)*sin(1.745329252*aa_index);
				aa_index += 1;
			c = sqrt( sum_cos**2 + sum_sin**2 );	
			if c>max_pattern:
				max_pattern=c;
		return max_pattern;
		
	def getMIPSignal9(self, sequence):
		# this and others from "Cleavage-site motifs in mitochondrial targeting peptides"
		score = len(re.findall("R.[FLI]",sequence[1:200]));
		return score;
	
	def getMIPSignal10(self, sequence):
		score = len(re.findall("R.Y[SA]",sequence[1:200]));
		return score;
		
	def getMIPSignal11(self, sequence):
		score = len(re.findall("R..S",sequence[1:200]));
		return score;
		
	def getMIPSignal12(self, sequence):
		# discriminant function Nakai'92  A knowledge base for predicting protein localization sites
		coeff = { 'R': 0.116, 'D' : -0.238 , 'P' : -0.253, 'E' : -0.233, 'G' : -0.25,
			'Q' : -0.55, 'H' : -0.239, 'K' : -0.113, 'N' : -0.134, 'Y' : -0.157};
		score = 0;
		for aa in sequence[1:20]:
			if aa in coeff.keys():
				score += coeff[aa];
				
		return score;
		
	def getMIPSignal13(self, sequence):
		# discriminant function Nakai'92  A knowledge base for predicting protein localization sites + Sebs Rule: only few S
		coeff = { 'R': 0.116, 'D' : -0.238 , 'P' : -0.253, 'E' : -0.233, 'G' : -0.25,
			'Q' : -0.55, 'H' : -0.239, 'K' : -0.113, 'N' : -0.134, 'Y' : -0.157, 'S' : -0.1};
		score = 0;
		for aa in sequence[1:20]:
			if aa in coeff.keys():
				score += coeff[aa];
				
		return score;
	
	def getAliphaticHelix(self, sequence):
		# Eisenberg hydrophobic moment
		# Eisenberg paper Hydrophobic Moments and Protein structure Hydrophobicity Values:
		hydr = {"A" : 0.25, "C" : 0.04, "D" : -0.72, "E" : -0.62, "F" : 0.61,
				"G" : 0.16, "H" : -0.4, "I" : 0.73, "K" : -1.1, "L" : 0.53,
				"M" : 0.26, "N" : -0.64, "P" : -0.07, "Q" : -0.69, "R" : -1.8,
				"S" : -0.26, "T" : -0.18, "V" : 0.54, "W" : 0.37, "Y" : 0.02};
		max_pattern = 0;
		for pos in range(0,len(sequence)-18):
			sum_cos = 0;
			sum_sin = 0;
			aa_index = 1;
			for aa in sequence[pos:pos+18]:
				if aa in hydr.keys():
					sum_cos += hydr[aa]*cos(1.745329252*aa_index);
					sum_sin += hydr[aa]*sin(1.745329252*aa_index);
				aa_index += 1;
			c = sqrt( sum_cos**2 + sum_sin**2 );	
			if c>max_pattern:
				max_pattern=c;
		return max_pattern;


	def getAliphaticHelixNormed(self, sequence):
		# Eisenberg hydrophobic moment
		# other hydrophobic moments
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880};
		max_pattern = 0;
		for pos in range(0,len(sequence)-18):
			sum_cos = 0;
			sum_sin = 0;
			aa_index = 1;
			overall = 0.0;
			for aa in sequence[pos:pos+18]:
				if aa in hydr.keys():
					sum_cos += hydr[aa]*cos(2.792526803*aa_index);
					sum_sin += hydr[aa]*sin(2.792526803*aa_index);
					overall += hydr[aa];
				aa_index += 1;
			c = sqrt( sum_cos**2 + sum_sin**2 );
			if overall != 0:
				c = c / overall;
			else:
				c = 0;
			if c>max_pattern:
				max_pattern=c;
		return max_pattern;

	def getAliphaticSheet(self, sequence):
		# Eisenberg hydrophobic moment
		# Eisenberg paper Hydrophobic Moments and Protein structure Hydrophobicity Values:
		hydr = {"A" : 0.25, "C" : 0.04, "D" : -0.72, "E" : -0.62, "F" : 0.61,
				"G" : 0.16, "H" : -0.4, "I" : 0.73, "K" : -1.1, "L" : 0.53,
				"M" : 0.26, "N" : -0.64, "P" : -0.07, "Q" : -0.69, "R" : -1.8,
				"S" : -0.26, "T" : -0.18, "V" : 0.54, "W" : 0.37, "Y" : 0.02};
		max_pattern = 0;
		for pos in range(0,len(sequence)-18):
			sum_cos = 0;
			sum_sin = 0;
			aa_index = 1;
			for aa in sequence[pos:pos+18]:
				if aa in hydr.keys():
					sum_cos += hydr[aa]*cos(2.792526803*aa_index);
					sum_sin += hydr[aa]*sin(2.792526803*aa_index);
				aa_index += 1;
			c = sqrt( sum_cos**2 + sum_sin**2 );	
			if c>max_pattern:
				max_pattern=c;
		return max_pattern;


	def getAliphaticSheetNormed(self, sequence):
		# Eisenberg hydrophobic moment
		# other hydrophobic moments
		hydr = {"A" : 0.616, "C" : 0.680, "D" : 0.028, "E" : 0.043, "F" : 1.0,
				"G" : 0.501, "H" : 0.165, "I" : 0.943, "K" : 0.283, "L" : 0.943,
				"M" : 0.738, "N" : 0.236, "P" : 0.711, "Q" : 0.251, "R" : 0.0,
				"S" : 0.359, "T" : 0.450, "V" : 0.825, "W" : 0.878, "Y" : 0.880};
		max_pattern = 0;
		for pos in range(0,len(sequence)-18):
			sum_cos = 0;
			sum_sin = 0;
			aa_index = 1;
			overall = 0.0;
			for aa in sequence[pos:pos+18]:
				if aa in hydr.keys():
					sum_cos += hydr[aa]*cos(1.745329252*aa_index);
					sum_sin += hydr[aa]*sin(1.745329252*aa_index);
					overall += hydr[aa];
				aa_index += 1;
			c = sqrt( sum_cos**2 + sum_sin**2 );
			if overall != 0:
				c = c / overall;
			else:
				c = 0;
			if c>max_pattern:
				max_pattern=c;
		return max_pattern;

	

	def getChlorSignal(self, sequence):
		# few acids, often hydroxylated,[Ema02] many S and T [Nak00]
		# at N-terminus often M or A [Ema02]
		# motif at cleavage site VRA | AAV
		pos1 = ['P','S','T'];
		pos2 = ['M','A'];
		score = 0;
		score += 3*len(re.findall("VRA",sequence));
		score += 3*len(re.findall("AAV",sequence));
		for aa in sequence:
			if aa in pos1:
				score += 1;
		for aa in sequence[1:100]:
			if aa in pos2:
				score += 2;
		return score;
		
	def getChlorSignal2(self, sequence):
		score = 0;
		score += len(re.findall("^MA",sequence));
		return score;
		
	def getVesicularSignal(self, sequence):
		score = len(re.findall("D.E",sequence[1:50]));
		score += len(re.findall("D.E",sequence[-50:]));
		return score;

	def getCleavageSite(self, sequence):
		"""
		vanHeinje 1986 "A new method for cleavage site prediction"
		weight_matrix = { "A": 
		[16,13,14,15,20,18,18,17,25,15,47,6,80,18,6,14.5],
		"C":[3,6,9,7,9,14,6,8,5,6,19,3,9,8,3,4.5],
		"D":[0,0,0,0,0,0,0,0,5,3,0,5,0,10,11,8.9],
		"E":[0,0,0,1,0,0,0,0,3,7,0,7,0,13,14,10.0],
		"F":[13,9,11,11,6,7,18,13,4,5,0,13,0,6,4,5.6],
		"G":[4,4,3,6,3,13,3,2,19,34,5,7,39,10,7,12.1],
		"H":[0,0,0,0,0,1,1,0,5,0,0,6,0,4,2,3.4],
		"I":[15,15,8,6,11,5,4,8,5,1,10,5,0,8,7,7.4],
		"K":[0,0,0,1,0,0,1,0,0,4,0,2,0,11,9,11.3],
		"L":[71,68,72,79,78,45,64,49,10,23,8,20,1,8,4,12.1],
		"M":[0,3,7,4,1,6,2,2,0,0,0,1,0,1,2,2.7],
		"N":[0,1,0,1,1,0,0,0,3,3,0,10,0,4,7,7.1],
		"P":[2,0,2,0,0,4,1,8,20,14,0,1,3,0,22,7.4],
		"Q":[0,0,0,1,0,6,1,0,10,8,0,18,3,19,10,6.3],
		"R":[2,0,0,0,0,1,0,0,7,4,0,15,0,12,9,7.6],
		"S":[9,3,8,6,13,10,15,16,26,11,23,17,20,15,10,11.4],
		"T":[2,10,5,4,5,13,7,7,12,6,17,8,6,3,10,9.7],
		"V":[20,25,15,18,13,15,11,27,0,12,32,3,0,8,17,11.1],
		"W":[4,3,3,1,1,2,6,3,1,3,0,9,0,2,0,1.8],
		"Y":[0,1,4,0,0,1,3,1,1,2,0,5,0,1,7,5.6]};

		for aa in weight_matrix.keys():
			s = "\""+aa+"\":[";
			for elem in weight_matrix[aa][:-1]:
				s += str( log(float(elem+1) / float(weight_matrix[aa][-1])) ) + ", ";
			s += "],";
			print(s);
		print("####");
		"""

		w = {"A":[0.15906469463, -0.0350913198113, 0.0339015516757, 0.0984400728133, 0.370373788297, 0.27029032974, 
		0.27029032974, 0.21622310847, 0.583947888595, 0.0984400728133, 1.19705236148, -0.728238500371, 
		1.72030050525, 0, 0.27029032974, -0.728238500371 ],
		"C":[-0.117783035656, 0.441832752279, 0.798507696218, 0.575364144904, 0.798507696218, 1.20397280433, 
		0.441832752279, 0.69314718056, 0.287682072452, 0.441832752279, 1.49165487678, -0.117783035656, 
		0.798507696218, 0, 0.69314718056, -0.117783035656 ],
		"E":[-2.30258509299, -2.30258509299, -2.30258509299, -1.60943791243, -2.30258509299, -2.30258509299, 
		-2.30258509299, -2.30258509299, -0.916290731874, -0.223143551314, -2.30258509299, -0.223143551314, 
		-2.30258509299, 0, 0.336472236621, 0.405465108108 ],
		"D":[-2.18605127674, -2.18605127674, -2.18605127674, -2.18605127674, -2.18605127674, -2.18605127674, 
		-2.18605127674, -2.18605127674, -0.39429180751, -0.799756915618, -2.18605127674, -0.39429180751, 
		-2.18605127674, 0, 0.21184399606, 0.29885537305 ],
		"G":[-0.883767540169, -0.883767540169, -1.10691109148, -0.547295303547, -1.10691109148, 0.145851877013, 
		-1.10691109148, -1.39459316393, 0.502526820951, 1.06214260889, -0.701445983375, -0.413763910923, 
		1.19567400151, 0, -0.0953101798043, -0.413763910923 ],
		"F":[0.916290731874, 0.579818495253, 0.762140052047, 0.762140052047, 0.223143551314, 0.356674943939, 
		1.22167238143, 0.916290731874, -0.113328685307, 0.068992871487, -1.72276659774, 0.916290731874, 
		-1.72276659774, 0, 0.223143551314, -0.113328685307 ],
		"I":[0.77110872203, 0.77110872203, 0.195744577126, -0.0555698511548, 0.483426649578, -0.209720530982, 
		-0.392042087776, 0.195744577126, -0.209720530982, -1.30833281965, 0.396415272588, -0.209720530982, 			-2.00148000021, 0, 0.195744577126, 0.0779615414697 ],
		"H":[-1.22377543162, -1.22377543162, -1.22377543162, -1.22377543162, -1.22377543162, -0.530628251062, 
		-0.530628251062, -1.22377543162, 0.567984037606, -1.22377543162, -1.22377543162, 0.722134717433, 
		-1.22377543162, 0, 0.385662480812, -0.125163142954 ],
		"K":[-2.42480272572, -2.42480272572, -2.42480272572, -1.73165554516, -2.42480272572, -2.42480272572, 
		-1.73165554516, -2.42480272572, -2.42480272572, -0.815364813284, -2.42480272572, -1.32619043705, 
		-2.42480272572, 0, 0.0601039240697, -0.122217632724 ],
		"M":[-0.99325177301, 0.39304258811, 1.08618976867, 0.616186139424, -0.30010459245, 0.952658376045, 
		0.105360515658, 0.105360515658, -0.99325177301, -0.99325177301, -0.99325177301, -0.30010459245, 
		-0.99325177301, 0, -0.30010459245, 0.105360515658 ],
		"L":[1.78346066641, 1.74090105199, 1.79725398855, 1.88882118207, 1.87624239986, 1.33543594389, 
		1.68118181729, 1.41881755283, -0.0953101798043, 0.684848377745, -0.295980875266, 0.551316985121, 
		-1.80005827204, 0, -0.295980875266, -0.883767540169],
		"N":[-1.96009478405, -1.26694760349, -1.96009478405, -1.26694760349, -1.26694760349, -1.96009478405, 
		-1.96009478405, -1.96009478405, -0.573800422927, -0.573800422927, -1.96009478405, 0.437800488751, 
		-1.96009478405, 0, -0.350656871613, 0.119346757633 ],
		"Q":[-1.8405496334, -1.8405496334, -1.8405496334, -1.14740245284, -1.8405496334, 0.105360515658, 
		-1.14740245284, -1.8405496334, 0.557345639401, 0.356674943939, -1.8405496334, 1.10388934577, 
		-0.454255272278, 0, 1.15518264016, 0.557345639401 ],
		"P":[-0.902867711542, -2.00148000021, -0.902867711542, -2.00148000021, -2.00148000021, -0.392042087776, 
		-1.30833281965, 0.195744577126, 1.04304243751, 0.706570200892, -2.00148000021, -1.30833281965, 
		-0.61518563909, 0, -2.00148000021, 1.13401421572 ],
		"S":[-0.131028262406, -1.04731899428, -0.236388778064, -0.487703206345, 0.205443974215, -0.0357180826021, 
		0.338975366839, 0.399599988656, 0.862223510604, 0.0512932943876, 0.744440474947, 0.456758402496, 
		0.610909082323, 0, 0.338975366839, -0.0357180826021 ],
		"R":[-0.929535958624, -2.02814824729, -2.02814824729, -2.02814824729, -2.02814824729, -1.33500106673, 
		-2.02814824729, -2.02814824729, 0.0512932943876, -0.418710334858, -2.02814824729, 0.744440474947, 
		-2.02814824729, 0, 0.536801110169, 0.274436845702 ],
		"T":[-1.17351359684, 0.125769387289, -0.480366416281, -0.662687973075, -0.480366416281, 0.366931444106, 
		-0.19268434383, -0.19268434383, 0.292823471952, -0.326215736454, 0.618245872387, -0.0749013081731, 
		-0.326215736454, 0, -0.885831524389, 0.125769387289 ],
		"W":[1.02165124753, 0.798507696218, 0.798507696218, 0.105360515658, 0.105360515658, 0.510825623766, 
		1.35812348415, 0.798507696218, 0.105360515658, 0.798507696218, -0.587786664902, 1.71479842809, 
		-0.587786664902, 0, 0.510825623766, -0.587786664902 ],
		"V":[0.637577329405, 0.851151429703, 0.365643613921, 0.537493870848, 0.232112221297, 0.365643613921, 
		0.0779615414697, 0.925259401857, -2.40694510832, 0.158004249143, 1.08956245315, -1.0206507472, 
		-2.40694510832, 0, -0.209720530982, 0.483426649578 ],
		"Y":[-1.72276659774, -1.02961941718, -0.113328685307, -1.72276659774, -1.72276659774, -1.02961941718, 
		-0.336472236621, -1.02961941718, -1.02961941718, -0.624154309073, -1.72276659774, 0.068992871487, 			-1.72276659774, 0, -1.02961941718, 0.356674943939 ]};

		score_max = 0;
		pos_max = 0;
		end_pos = len(sequence)-2;
		if len(sequence)>202:
			end_pos = 200;
		for aa_pos in range(14,end_pos):
			score_pos = 0;
			index_pos = 0;
			for k in range(aa_pos-14,aa_pos-1):
				if sequence[k] in w.keys():
					score_pos += w[sequence[k]][index_pos];
				index_pos += 1;
			index_pos +=1;
			for k in range(aa_pos+1,aa_pos+2):
				if sequence[k] in w.keys():
					score_pos += w[sequence[k]][index_pos];
				index_pos += 1;
			if score_pos > score_max:
				score_max = score_pos;
				pos_max = aa_pos;
			#test = str(aa_pos)+" ";
			#for k in range(1,ceil(score_max)):
			#	test += "|";
			#print(test);
		#print(pos_max);
		#print(score_max);
		#print("--");
		return pos_max;

	def getSecretedSignal1(self, sequence):
		score = 0;
		g = self.getGlycoSiteSignal(sequence);
		if g == 0:
			score += 3;
		s = len(re.findall("S",sequence[1:20]));
		if s == 0:
			score += 2;
		g = len(re.findall("W.W",sequence));
		if g == 0:
			score += 2;
		p = self.getPeroxiSignal8(sequence);
		if p <= 25:
			score +=4;
		elif p<=50:
			score +=1;
		p = len(re.findall("P",sequence));
		if p<16:
			score += 2;
		return score;
