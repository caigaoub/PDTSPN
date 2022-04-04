import json
import ast
import glob2
import re
import pandas as pd
import numpy as np
def get_statistics_mode11():

	path = 'C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/python/'

	''' combine all result files into one '''
	implement_combination = True
	if implement_combination:
		filespool = glob2.glob(path + "ResultsMode_1-1/" + "ret_cvxp_*")
		with open(path + "combined_Mode11.py", "w") as outfile:
		    for f in filespool:
		        with open(f, "r") as infile:
		            outfile.write( infile.read() + "\n")

	rile = open(path + 'combined_Mode11.py', "r")
	line_ = rile.readline()
	RET = []
	while line_ != '':
		str_ = re.split('{|}|:|=|,|\'| |\n', line_)
		str_ = list(filter(None, str_))
		print(str_)
		RET.append(str_)
		line_ = rile.readline()

	df = pd.DataFrame(RET)
	df.to_csv('stats_11.csv', sep=',')

def get_statistics_modeMM():
	path = 'C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/python/'
	''' combine all result files into one '''
	implement_combination = True
	if implement_combination:
		filespool = glob2.glob(path + "ResultsMode_M-M/" + "ret_cvxp_*")
		with open(path + "combined.py", "w") as outfile:
		    for f in filespool:
		        with open(f, "r") as infile:
		            outfile.write( infile.read() + "\n")

	rile = open(path + 'combined_ModeMM.py', "r")
	line_ = rile.readline()
	RET = []
	while line_ != '':
		str_ = re.split('{|}|:|=|,|\'| |\n', line_)
		str_ = list(filter(None, str_))
		print(str_)
		RET.append(str_)
		line_ = rile.readline()

	df = pd.DataFrame(RET)
	df.to_csv('stats_MM.csv', sep=',')

def get_statistics_separators():
	path = 'C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/python/'
	implement_combination = True
	if implement_combination:
		filespool = glob2.glob(path + "ResultsSeparatorsEvaluation/" + "ret_cvxp_*")
		with open(path + "combined_separators.py", "w") as outfile:
		    for f in filespool:
		        with open(f, "r") as infile:
		            outfile.write(infile.read() + "\n")

	rile = open(path + 'combined_separators.py', "r")
	line_ = rile.readline()
	RET = []
	while line_ != '':
		# print(line_.replace('\n',''))
		str_ = re.split('{|}|:|=|,|\'| |\[|\]|\{|\}|\(|\)|\n', line_)
		str_ = list(filter(None, str_))
		# print(str_)
		RET.append(str_)
		line_ = rile.readline()
	TIME = {'cvxp_6_1_6': 0.06796940000000001, 'cvxp_6_1_9': 0.0993252, 'cvxp_6_1_12': 0.13763120000000006, 'cvxp_6_1_15': 0.16048839999999998, 'cvxp_6_2_6': 0.08585369999999992, 'cvxp_6_2_9': 0.12542830000000005, 'cvxp_6_2_12': 0.18243089999999995, 'cvxp_6_2_15': 0.22482970000000013, 'cvxp_6_3_6': 0.06251479999999998, 'cvxp_6_3_9': 0.10257439999999995, 'cvxp_6_3_12': 0.13755289999999998, 'cvxp_6_3_15': 0.19219410000000003, 'cvxp_6_4_6': 0.08057619999999988, 'cvxp_6_4_9': 0.12203169999999997, 'cvxp_6_4_12': 0.16732710000000006, 'cvxp_6_4_15': 0.21061209999999964, 'cvxp_6_5_6': 0.07894159999999983, 'cvxp_6_5_9': 0.1232527000000001, 'cvxp_6_5_12': 0.1707506000000003, 'cvxp_6_5_15': 0.21475890000000009, 'cvxp_6_11_6': 0.07576360000000015, 'cvxp_6_11_9': 0.12863360000000013, 'cvxp_6_11_12': 0.17337130000000034, 'cvxp_6_11_15': 0.21685109999999996, 'cvxp_6_12_6': 0.11215909999999996, 'cvxp_6_12_9': 0.19367709999999994, 'cvxp_6_12_12': 0.2376541999999997, 'cvxp_6_12_15': 0.28845840000000056, 'cvxp_6_13_6': 0.10591279999999959, 'cvxp_6_13_9': 0.1705918999999998, 'cvxp_6_13_12': 0.22969989999999996, 'cvxp_6_13_15': 0.2985459000000006, 'cvxp_6_14_6': 0.05824050000000014, 'cvxp_6_14_9': 0.09739729999999991, 'cvxp_6_14_12': 0.13065090000000001, 'cvxp_6_14_15': 0.16432380000000002, 'cvxp_6_15_6': 0.06965500000000002, 'cvxp_6_15_9': 0.1188760000000002, 'cvxp_6_15_12': 0.1390311000000004, 'cvxp_6_15_15': 0.18555299999999963, 'cvxp_6_21_6': 0.08596709999999952, 'cvxp_6_21_9': 0.14286270000000023, 'cvxp_6_21_12': 0.20431689999999936, 'cvxp_6_21_15': 0.26463780000000003, 'cvxp_6_22_6': 0.06301609999999958, 'cvxp_6_22_9': 0.10426629999999992, 'cvxp_6_22_12': 0.14244009999999996, 'cvxp_6_22_15': 0.18945860000000003, 'cvxp_6_23_6': 0.08686859999999985, 'cvxp_6_23_9': 0.13239770000000028, 'cvxp_6_23_12': 0.19324779999999908, 'cvxp_6_23_15': 0.24289369999999977, 'cvxp_6_24_6': 0.08505689999999966, 'cvxp_6_24_9': 0.13133019999999895, 'cvxp_6_24_12': 0.18112490000000037, 'cvxp_6_24_15': 0.23454389999999847, 'cvxp_6_25_6': 0.054114600000000124, 'cvxp_6_25_9': 0.08730139999999942, 'cvxp_6_25_12': 0.13167439999999964, 'cvxp_6_25_15': 0.17793439999999983, 'cvxp_8_1_6': 0.10913069999999969, 'cvxp_8_1_9': 0.17588189999999848, 'cvxp_8_1_12': 0.23764270000000032, 'cvxp_8_1_15': 0.30784069999999986, 'cvxp_8_2_6': 0.0857600999999999, 'cvxp_8_2_9': 0.14126940000000054, 'cvxp_8_2_12': 0.19203790000000076, 'cvxp_8_2_15': 0.2358009999999986, 'cvxp_8_3_6': 0.10584839999999929, 'cvxp_8_3_9': 0.16420859999999848, 'cvxp_8_3_12': 0.22925859999999965, 'cvxp_8_3_15': 0.3025196000000001, 'cvxp_8_4_6': 0.1271678000000005, 'cvxp_8_4_9': 0.20128389999999996, 'cvxp_8_4_12': 0.27511559999999946, 'cvxp_8_4_15': 0.3408718000000004, 'cvxp_8_5_6': 0.07327940000000055, 'cvxp_8_5_9': 0.13192480000000018, 'cvxp_8_5_12': 0.17923689999999937, 'cvxp_8_5_15': 0.22481370000000034, 'cvxp_8_11_6': 0.08222540000000045, 'cvxp_8_11_9': 0.14193820000000024, 'cvxp_8_11_12': 0.1844184000000002, 'cvxp_8_11_15': 0.2474607000000013, 'cvxp_8_12_6': 0.11487669999999994, 'cvxp_8_12_9': 0.1733722000000011, 'cvxp_8_12_12': 0.2627967000000009, 'cvxp_8_12_15': 0.30573799999999984, 'cvxp_8_13_6': 0.12654989999999877, 'cvxp_8_13_9': 0.19771189999999983, 'cvxp_8_13_12': 0.2831140000000012, 'cvxp_8_13_15': 0.3431139999999999, 'cvxp_8_14_6': 0.10157699999999892, 'cvxp_8_14_9': 0.15708380000000233, 'cvxp_8_14_12': 0.2271745999999979, 'cvxp_8_14_15': 0.2768000999999991, 'cvxp_8_15_6': 0.12295519999999982, 'cvxp_8_15_9': 0.19986469999999912, 'cvxp_8_15_12': 0.2651237000000002, 'cvxp_8_15_15': 0.34439820000000054, 'cvxp_8_21_6': 0.10747229999999774, 'cvxp_8_21_9': 0.16529060000000229, 'cvxp_8_21_12': 0.22953180000000017, 'cvxp_8_21_15': 0.29324370000000144, 'cvxp_8_22_6': 0.10867939999999976, 'cvxp_8_22_9': 0.17210589999999826, 'cvxp_8_22_12': 0.24346889999999988, 'cvxp_8_22_15': 0.29554170000000113, 'cvxp_8_23_6': 0.09469019999999873, 'cvxp_8_23_9': 0.1680141000000006, 'cvxp_8_23_12': 0.2321163999999989, 'cvxp_8_23_15': 0.3031063999999972, 'cvxp_8_24_6': 0.11908049999999903, 'cvxp_8_24_9': 0.18873909999999938, 'cvxp_8_24_12': 0.2491139000000011, 'cvxp_8_24_15': 0.31074999999999875, 'cvxp_8_25_6': 0.09977159999999685, 'cvxp_8_25_9': 0.15846379999999982, 'cvxp_8_25_12': 0.22603379999999973, 'cvxp_8_25_15': 0.2820323999999985, 'cvxp_10_1_6': 0.09917079999999956, 'cvxp_10_1_9': 0.1695500999999986, 'cvxp_10_1_12': 0.23239940000000203, 'cvxp_10_1_15': 0.30824989999999985, 'cvxp_10_2_6': 0.142725200000001, 'cvxp_10_2_9': 0.21387490000000042, 'cvxp_10_2_12': 0.300608399999998, 'cvxp_10_2_15': 0.38185889999999745, 'cvxp_10_3_6': 0.1428244999999997, 'cvxp_10_3_9': 0.23009469999999865, 'cvxp_10_3_12': 0.31268030000000024, 'cvxp_10_3_15': 0.3989244999999997, 'cvxp_10_4_6': 0.08766319999999794, 'cvxp_10_4_9': 0.12722289999999958, 'cvxp_10_4_12': 0.1796717000000001, 'cvxp_10_4_15': 0.23536590000000146, 'cvxp_10_5_6': 0.12213570000000118, 'cvxp_10_5_9': 0.20520819999999773, 'cvxp_10_5_12': 0.2722371999999993, 'cvxp_10_5_15': 0.36566809999999705, 'cvxp_10_11_6': 0.14490590000000125, 'cvxp_10_11_9': 0.2142382999999981, 'cvxp_10_11_12': 0.2965106999999989, 'cvxp_10_11_15': 0.37106119999999976, 'cvxp_10_12_6': 0.12067130000000148, 'cvxp_10_12_9': 0.2012122000000005, 'cvxp_10_12_12': 0.273325100000001, 'cvxp_10_12_15': 0.3472705000000005, 'cvxp_10_13_6': 0.10475890000000021, 'cvxp_10_13_9': 0.142414500000001, 'cvxp_10_13_12': 0.20256749999999712, 'cvxp_10_13_15': 0.2530731000000017, 'cvxp_10_14_6': 0.14020100000000113, 'cvxp_10_14_9': 0.22733129999999946, 'cvxp_10_14_12': 0.27658879999999897, 'cvxp_10_14_15': 0.37242789999999815, 'cvxp_10_15_6': 0.10653879999999916, 'cvxp_10_15_9': 0.18330740000000034, 'cvxp_10_15_12': 0.2418143999999991, 'cvxp_10_15_15': 0.3135323999999997, 'cvxp_10_21_6': 0.15995969999999815, 'cvxp_10_21_9': 0.253223000000002, 'cvxp_10_21_12': 0.3402902999999995, 'cvxp_10_21_15': 0.40108809999999906, 'cvxp_10_22_6': 0.13921609999999873, 'cvxp_10_22_9': 0.2139351000000005, 'cvxp_10_22_12': 0.301054900000004, 'cvxp_10_22_15': 0.38799960000000056, 'cvxp_10_23_6': 0.1621787999999995, 'cvxp_10_23_9': 0.22305709999999834, 'cvxp_10_23_12': 0.32308439999999905, 'cvxp_10_23_15': 0.4024023999999997, 'cvxp_10_24_6': 0.13159009999999682, 'cvxp_10_24_9': 0.19886460000000028, 'cvxp_10_24_12': 0.2733922000000035, 'cvxp_10_24_15': 0.34307460000000134, 'cvxp_10_25_6': 0.10996409999999912, 'cvxp_10_25_9': 0.18589749999999583, 'cvxp_10_25_12': 0.24776640000000327, 'cvxp_10_25_15': 0.3174457000000004, 'cvxp_12_1_6': 0.14457389999999748, 'cvxp_12_1_9': 0.19378180000000356, 'cvxp_12_1_12': 0.27106030000000203, 'cvxp_12_1_15': 0.3512211999999977, 'cvxp_12_2_6': 0.19390419999999864, 'cvxp_12_2_9': 0.3026389999999992, 'cvxp_12_2_12': 0.4036104999999992, 'cvxp_12_2_15': 0.5285646000000028, 'cvxp_12_3_6': 0.15250359999999574, 'cvxp_12_3_9': 0.23387000000000313, 'cvxp_12_3_12': 0.3264299999999949, 'cvxp_12_3_15': 0.40544009999999986, 'cvxp_12_4_6': 0.1556434999999965, 'cvxp_12_4_9': 0.2530378000000013, 'cvxp_12_4_12': 0.35458369999999917, 'cvxp_12_4_15': 0.4411198000000027, 'cvxp_12_5_6': 0.1511698999999993, 'cvxp_12_5_9': 0.22391609999999673, 'cvxp_12_5_12': 0.3117984000000007, 'cvxp_12_5_15': 0.39209100000000063, 'cvxp_12_11_6': 0.1406244000000001, 'cvxp_12_11_9': 0.2142281000000068, 'cvxp_12_11_12': 0.2957338000000007, 'cvxp_12_11_15': 0.36545449999999846, 'cvxp_12_12_6': 0.13431149999999548, 'cvxp_12_12_9': 0.21221690000000137, 'cvxp_12_12_12': 0.29268839999999585, 'cvxp_12_12_15': 0.3537770000000009, 'cvxp_12_13_6': 0.17670640000000049, 'cvxp_12_13_9': 0.28582380000000285, 'cvxp_12_13_12': 0.37933060000000296, 'cvxp_12_13_15': 0.4808540000000008, 'cvxp_12_14_6': 0.19295879999999954, 'cvxp_12_14_9': 0.2955455000000029, 'cvxp_12_14_12': 0.3981986000000006, 'cvxp_12_14_15': 0.5104145999999972, 'cvxp_12_15_6': 0.17266899999999907, 'cvxp_12_15_9': 0.28402089999999447, 'cvxp_12_15_12': 0.37781900000000235, 'cvxp_12_15_15': 0.48799410000000165, 'cvxp_12_21_6': 0.16357659999999896, 'cvxp_12_21_9': 0.25231589999999926, 'cvxp_12_21_12': 0.3554302000000007, 'cvxp_12_21_15': 0.4471418999999983, 'cvxp_12_22_6': 0.14276290000000103, 'cvxp_12_22_9': 0.21652639999999934, 'cvxp_12_22_12': 0.29601209999999867, 'cvxp_12_22_15': 0.3816096999999985, 'cvxp_12_23_6': 0.13940589999999986, 'cvxp_12_23_9': 0.22523379999999804, 'cvxp_12_23_12': 0.31260239999999584, 'cvxp_12_23_15': 0.38859029999999706, 'cvxp_12_24_6': 0.11356150000000298, 'cvxp_12_24_9': 0.19340429999999742, 'cvxp_12_24_12': 0.26175500000000085, 'cvxp_12_24_15': 0.32885989999999765, 'cvxp_12_25_6': 0.10665799999999592, 'cvxp_12_25_9': 0.16269600000000395, 'cvxp_12_25_12': 0.2216364000000013, 'cvxp_12_25_15': 0.292595900000002, 'cvxp_14_1_6': 0.20343469999999542, 'cvxp_14_1_9': 0.30521110000000107, 'cvxp_14_1_12': 0.427149600000007, 'cvxp_14_1_15': 0.5391557999999961, 'cvxp_14_2_6': 0.17467830000000362, 'cvxp_14_2_9': 0.2615735999999984, 'cvxp_14_2_12': 0.36918960000000567, 'cvxp_14_2_15': 0.45569090000000045, 'cvxp_14_3_6': 0.21001479999999617, 'cvxp_14_3_9': 0.32944260000000014, 'cvxp_14_3_12': 0.44655669999999503, 'cvxp_14_3_15': 0.5653381000000053, 'cvxp_14_4_6': 0.16592769999999746, 'cvxp_14_4_9': 0.26708160000000447, 'cvxp_14_4_12': 0.3546674999999979, 'cvxp_14_4_15': 0.46410029999999836, 'cvxp_14_5_6': 0.14038520000000432, 'cvxp_14_5_9': 0.23697570000000212, 'cvxp_14_5_12': 0.32236540000000247, 'cvxp_14_5_15': 0.42588779999999815, 'cvxp_14_11_6': 0.17182289999999512, 'cvxp_14_11_9': 0.2628977999999975, 'cvxp_14_11_12': 0.3790164000000047, 'cvxp_14_11_15': 0.4759555999999989, 'cvxp_14_12_6': 0.17762410000000273, 'cvxp_14_12_9': 0.26486650000000367, 'cvxp_14_12_12': 0.3772616000000042, 'cvxp_14_12_15': 0.4506332999999998, 'cvxp_14_13_6': 0.1975145999999981, 'cvxp_14_13_9': 0.32538499999999715, 'cvxp_14_13_12': 0.4357545000000016, 'cvxp_14_13_15': 0.5683524000000091, 'cvxp_14_14_6': 0.20393339999999682, 'cvxp_14_14_9': 0.3257940000000019, 'cvxp_14_14_12': 0.43930860000000393, 'cvxp_14_14_15': 0.5523799000000054, 'cvxp_14_15_6': 0.19113180000000796, 'cvxp_14_15_9': 0.3109125000000006, 'cvxp_14_15_12': 0.43507470000000126, 'cvxp_14_15_15': 0.552634299999994, 'cvxp_14_21_6': 0.2129620999999986, 'cvxp_14_21_9': 0.3446892000000048, 'cvxp_14_21_12': 0.4569694000000055, 'cvxp_14_21_15': 0.5984260999999975, 'cvxp_14_22_6': 0.17300430000000233, 'cvxp_14_22_9': 0.27002040000000704, 'cvxp_14_22_12': 0.3573181999999946, 'cvxp_14_22_15': 0.4566073000000017, 'cvxp_14_23_6': 0.14332780000000866, 'cvxp_14_23_9': 0.23828429999998946, 'cvxp_14_23_12': 0.3321536999999921, 'cvxp_14_23_15': 0.44548360000000287, 'cvxp_14_24_6': 0.2203970000000055, 'cvxp_14_24_9': 0.35305570000001296, 'cvxp_14_24_12': 0.4736200999999909, 'cvxp_14_24_15': 0.5969478000000095, 'cvxp_14_25_6': 0.18890249999999753, 'cvxp_14_25_9': 0.31070569999999975, 'cvxp_14_25_12': 0.4152562999999958, 'cvxp_14_25_15': 0.5525683999999984, 'cvxp_16_1_6': 0.14209119999999587, 'cvxp_16_1_9': 0.25475380000000314, 'cvxp_16_1_12': 0.3368217999999956, 'cvxp_16_1_15': 0.4452609999999879, 'cvxp_16_2_6': 0.18985649999999055, 'cvxp_16_2_9': 0.3150856000000033, 'cvxp_16_2_12': 0.44357580000000496, 'cvxp_16_2_15': 0.5682450000000046, 'cvxp_16_3_6': 0.1862384000000077, 'cvxp_16_3_9': 0.32372999999999763, 'cvxp_16_3_12': 0.4187404999999984, 'cvxp_16_3_15': 0.567776699999996, 'cvxp_16_4_6': 0.21226089999998976, 'cvxp_16_4_9': 0.32395389999999225, 'cvxp_16_4_12': 0.4611289000000056, 'cvxp_16_4_15': 0.5735648999999938, 'cvxp_16_5_6': 0.2294085999999993, 'cvxp_16_5_9': 0.34230390000000455, 'cvxp_16_5_12': 0.48904020000000514, 'cvxp_16_5_15': 0.6126418000000058, 'cvxp_16_11_6': 0.21556280000000072, 'cvxp_16_11_9': 0.3315419999999989, 'cvxp_16_11_12': 0.46460280000000864, 'cvxp_16_11_15': 0.586514099999988, 'cvxp_16_12_6': 0.14088380000001166, 'cvxp_16_12_9': 0.2372569999999996, 'cvxp_16_12_12': 0.33325169999999105, 'cvxp_16_12_15': 0.42002750000000333, 'cvxp_16_13_6': 0.21068590000000142, 'cvxp_16_13_9': 0.3114731000000006, 'cvxp_16_13_12': 0.4422855000000112, 'cvxp_16_13_15': 0.5568351000000007, 'cvxp_16_14_6': 0.19226130000001262, 'cvxp_16_14_9': 0.3197361000000001, 'cvxp_16_14_12': 0.4156418999999971, 'cvxp_16_14_15': 0.5264985999999965, 'cvxp_16_15_6': 0.20114119999999502, 'cvxp_16_15_9': 0.34135700000000213, 'cvxp_16_15_12': 0.4528001000000046, 'cvxp_16_15_15': 0.5983943999999894, 'cvxp_16_21_6': 0.21330310000000452, 'cvxp_16_21_9': 0.3216977999999955, 'cvxp_16_21_12': 0.44684719999999345, 'cvxp_16_21_15': 0.5665528999999907, 'cvxp_16_22_6': 0.19435540000000628, 'cvxp_16_22_9': 0.3244232999999923, 'cvxp_16_22_12': 0.46002630000000977, 'cvxp_16_22_15': 0.567067300000005, 'cvxp_16_23_6': 0.1726846999999907, 'cvxp_16_23_9': 0.2991522999999887, 'cvxp_16_23_12': 0.3981621000000075, 'cvxp_16_23_15': 0.5049776999999978, 'cvxp_16_24_6': 0.136308800000009, 'cvxp_16_24_9': 0.2327615000000094, 'cvxp_16_24_12': 0.31500400000000184, 'cvxp_16_24_15': 0.39528149999999584, 'cvxp_16_25_6': 0.2251579999999933, 'cvxp_16_25_9': 0.37418430000001024, 'cvxp_16_25_12': 0.49784490000000403, 'cvxp_16_25_15': 0.6500283000000024, 'cvxp_18_1_6': 0.20762190000000658, 'cvxp_18_1_9': 0.33525990000001116, 'cvxp_18_1_12': 0.4472535000000022, 'cvxp_18_1_15': 0.5811117000000081, 'cvxp_18_2_6': 0.24929139999998995, 'cvxp_18_2_9': 0.4004354000000063, 'cvxp_18_2_12': 0.5493316999999962, 'cvxp_18_2_15': 0.7059601999999927, 'cvxp_18_3_6': 0.22303440000000307, 'cvxp_18_3_9': 0.36346370000001116, 'cvxp_18_3_12': 0.487219500000009, 'cvxp_18_3_15': 0.6155439000000058, 'cvxp_18_4_6': 0.16931150000000628, 'cvxp_18_4_9': 0.2742211999999995, 'cvxp_18_4_12': 0.3731801999999931, 'cvxp_18_4_15': 0.4770522000000028, 'cvxp_18_5_6': 0.18495270000001085, 'cvxp_18_5_9': 0.3064593000000002, 'cvxp_18_5_12': 0.42337589999999636, 'cvxp_18_5_15': 0.5281523999999962, 'cvxp_18_11_6': 0.19070529999999053, 'cvxp_18_11_9': 0.35491519999999355, 'cvxp_18_11_12': 0.48311879999999974, 'cvxp_18_11_15': 0.6276531000000034, 'cvxp_18_12_6': 0.2531452999999999, 'cvxp_18_12_9': 0.40815079999998716, 'cvxp_18_12_12': 0.5365608000000037, 'cvxp_18_12_15': 0.6860175999999996, 'cvxp_18_13_6': 0.1967458000000022, 'cvxp_18_13_9': 0.30438530000000696, 'cvxp_18_13_12': 0.43358950000001073, 'cvxp_18_13_15': 0.5531950999999964, 'cvxp_18_14_6': 0.20518409999999676, 'cvxp_18_14_9': 0.3219807999999915, 'cvxp_18_14_12': 0.4606078999999994, 'cvxp_18_14_15': 0.5770127999999914, 'cvxp_18_15_6': 0.21197749999998905, 'cvxp_18_15_9': 0.3462945999999931, 'cvxp_18_15_12': 0.44231850000001316, 'cvxp_18_15_15': 0.5712434000000002, 'cvxp_18_21_6': 0.18072930000001008, 'cvxp_18_21_9': 0.26111279999999226, 'cvxp_18_21_12': 0.37011419999998907, 'cvxp_18_21_15': 0.5026732000000038, 'cvxp_18_22_6': 0.19605969999999218, 'cvxp_18_22_9': 0.3086460000000102, 'cvxp_18_22_12': 0.4326880999999929, 'cvxp_18_22_15': 0.5611023000000017, 'cvxp_18_23_6': 0.27265549999999905, 'cvxp_18_23_9': 0.4327744999999936, 'cvxp_18_23_12': 0.6105833999999959, 'cvxp_18_23_15': 0.7801412000000028, 'cvxp_18_24_6': 0.21749499999999955, 'cvxp_18_24_9': 0.3528646999999978, 'cvxp_18_24_12': 0.4937520999999947, 'cvxp_18_24_15': 0.5928591000000125, 'cvxp_18_25_6': 0.23806940000000054, 'cvxp_18_25_9': 0.37166960000000415, 'cvxp_18_25_12': 0.49513539999999523, 'cvxp_18_25_15': 0.6327136000000024, 'cvxp_20_1_6': 0.28336990000001094, 'cvxp_20_1_9': 0.452736299999998, 'cvxp_20_1_12': 0.6229386999999917, 'cvxp_20_1_15': 0.7875062000000099, 'cvxp_20_2_6': 0.22961059999998668, 'cvxp_20_2_9': 0.34416810000000453, 'cvxp_20_2_12': 0.4954853000000128, 'cvxp_20_2_15': 0.6276027999999911, 'cvxp_20_3_6': 0.2561908999999929, 'cvxp_20_3_9': 0.38729179999999985, 'cvxp_20_3_12': 0.5090600000000052, 'cvxp_20_3_15': 0.6328815000000105, 'cvxp_20_4_6': 0.2233748000000162, 'cvxp_20_4_9': 0.3539150999999947, 'cvxp_20_4_12': 0.5018684000000064, 'cvxp_20_4_15': 0.6491088999999874, 'cvxp_20_5_6': 0.14047209999998245, 'cvxp_20_5_9': 0.2330748999999912, 'cvxp_20_5_12': 0.2961585000000184, 'cvxp_20_5_15': 0.38242270000000644, 'cvxp_20_11_6': 0.23044369999999503, 'cvxp_20_11_9': 0.3613843999999915, 'cvxp_20_11_12': 0.5000662999999861, 'cvxp_20_11_15': 0.6270660999999791, 'cvxp_20_12_6': 0.16734700000000657, 'cvxp_20_12_9': 0.2852435000000071, 'cvxp_20_12_12': 0.36421450000000277, 'cvxp_20_12_15': 0.5035172000000046, 'cvxp_20_13_6': 0.2545063999999968, 'cvxp_20_13_9': 0.41632549999999924, 'cvxp_20_13_12': 0.5603586000000007, 'cvxp_20_13_15': 0.7128103999999951, 'cvxp_20_14_6': 0.2011109000000033, 'cvxp_20_14_9': 0.3380035999999791, 'cvxp_20_14_12': 0.4589413999999863, 'cvxp_20_14_15': 0.6105176000000085, 'cvxp_20_15_6': 0.28701019999999744, 'cvxp_20_15_9': 0.45555459999999925, 'cvxp_20_15_12': 0.627375700000016, 'cvxp_20_15_15': 0.7919190000000071, 'cvxp_20_21_6': 0.2250613999999871, 'cvxp_20_21_9': 0.34283110000001216, 'cvxp_20_21_12': 0.5039911000000075, 'cvxp_20_21_15': 0.6098235000000045, 'cvxp_20_22_6': 0.22818269999999075, 'cvxp_20_22_9': 0.3788553000000263, 'cvxp_20_22_12': 0.5152132999999992, 'cvxp_20_22_15': 0.6347621999999831, 'cvxp_20_23_6': 0.22563529999999332, 'cvxp_20_23_9': 0.3764340999999831, 'cvxp_20_23_12': 0.5102567999999792, 'cvxp_20_23_15': 0.6525622000000055, 'cvxp_20_24_6': 0.22932160000002, 'cvxp_20_24_9': 0.36783760000000143, 'cvxp_20_24_12': 0.4986176000000171, 'cvxp_20_24_15': 0.6842673999999818, 'cvxp_20_25_6': 0.19952039999998306, 'cvxp_20_25_9': 0.35574989999997797, 'cvxp_20_25_12': 0.4879839999999831, 'cvxp_20_25_15': 0.6450759000000232, 'cvxp_22_1_6': 0.3298141000000214, 'cvxp_22_1_9': 0.5227117999999962, 'cvxp_22_1_12': 0.7221351000000027, 'cvxp_22_1_15': 0.9102034999999944, 'cvxp_22_2_6': 0.26939709999999195, 'cvxp_22_2_9': 0.4594654999999932, 'cvxp_22_2_12': 0.6122761000000025, 'cvxp_22_2_15': 0.8239868999999942, 'cvxp_22_3_6': 0.20460249999999292, 'cvxp_22_3_9': 0.3436776000000066, 'cvxp_22_3_12': 0.4816471999999976, 'cvxp_22_3_15': 0.603566700000016, 'cvxp_22_4_6': 0.28601280000000884, 'cvxp_22_4_9': 0.44836469999998485, 'cvxp_22_4_12': 0.6403907000000117, 'cvxp_22_4_15': 0.7805610999999999, 'cvxp_22_5_6': 0.21741950000000543, 'cvxp_22_5_9': 0.3365374999999915, 'cvxp_22_5_12': 0.4751009999999951, 'cvxp_22_5_15': 0.5794466000000114, 'cvxp_22_11_6': 0.2611114999999984, 'cvxp_22_11_9': 0.43464480000000094, 'cvxp_22_11_12': 0.5692196999999908, 'cvxp_22_11_15': 0.7546453000000213, 'cvxp_22_12_6': 0.2689517999999964, 'cvxp_22_12_9': 0.4310471000000007, 'cvxp_22_12_12': 0.5867163000000062, 'cvxp_22_12_15': 0.7376252999999906, 'cvxp_22_13_6': 0.3269859999999767, 'cvxp_22_13_9': 0.5235887999999989, 'cvxp_22_13_12': 0.7087787000000105, 'cvxp_22_13_15': 0.9057423999999799, 'cvxp_22_14_6': 0.2853168000000039, 'cvxp_22_14_9': 0.40855700000000184, 'cvxp_22_14_12': 0.5802374000000157, 'cvxp_22_14_15': 0.7096566000000166, 'cvxp_22_15_6': 0.2510481999999854, 'cvxp_22_15_9': 0.42054469999999355, 'cvxp_22_15_12': 0.5565781999999899, 'cvxp_22_15_15': 0.7198573999999951, 'cvxp_22_21_6': 0.2830381999999929, 'cvxp_22_21_9': 0.4548778999999854, 'cvxp_22_21_12': 0.6381834000000026, 'cvxp_22_21_15': 0.8133301999999958, 'cvxp_22_22_6': 0.27285630000000083, 'cvxp_22_22_9': 0.4478249000000005, 'cvxp_22_22_12': 0.5944618999999989, 'cvxp_22_22_15': 0.7719111000000112, 'cvxp_22_23_6': 0.24204860000000394, 'cvxp_22_23_9': 0.4251307999999767, 'cvxp_22_23_12': 0.5849475999999925, 'cvxp_22_23_15': 0.7522120000000143, 'cvxp_22_24_6': 0.29993930000000546, 'cvxp_22_24_9': 0.4899096000000043, 'cvxp_22_24_12': 0.6482844999999884, 'cvxp_22_24_15': 0.7968238999999926, 'cvxp_22_25_6': 0.26528899999999567, 'cvxp_22_25_9': 0.39047740000000886, 'cvxp_22_25_12': 0.5832144000000028, 'cvxp_22_25_15': 0.730111700000009, 'cvxp_24_1_6': 0.26228779999999574, 'cvxp_24_1_9': 0.43004159999998137, 'cvxp_24_1_12': 0.5649328000000082, 'cvxp_24_1_15': 0.7658709000000101, 'cvxp_24_2_6': 0.3184733999999878, 'cvxp_24_2_9': 0.49369340000001216, 'cvxp_24_2_12': 0.6966060999999968, 'cvxp_24_2_15': 0.8595353999999986, 'cvxp_24_3_6': 0.3137777999999969, 'cvxp_24_3_9': 0.5130314999999825, 'cvxp_24_3_12': 0.7317242000000022, 'cvxp_24_3_15': 0.938545199999993, 'cvxp_24_4_6': 0.3521296999999777, 'cvxp_24_4_9': 0.547101199999986, 'cvxp_24_4_12': 0.7431278999999904, 'cvxp_24_4_15': 0.9536917000000074, 'cvxp_24_5_6': 0.24223440000000096, 'cvxp_24_5_9': 0.3658959999999922, 'cvxp_24_5_12': 0.49843590000000404, 'cvxp_24_5_15': 0.6498284999999839, 'cvxp_24_11_6': 0.257619300000016, 'cvxp_24_11_9': 0.40962970000001064, 'cvxp_24_11_12': 0.5576441999999986, 'cvxp_24_11_15': 0.7023696999999913, 'cvxp_24_12_6': 0.3202588999999989, 'cvxp_24_12_9': 0.530182300000007, 'cvxp_24_12_12': 0.7014099000000158, 'cvxp_24_12_15': 0.881979699999988, 'cvxp_24_13_6': 0.34226150000000644, 'cvxp_24_13_9': 0.5330677000000037, 'cvxp_24_13_12': 0.7619674000000032, 'cvxp_24_13_15': 0.9565833999999995, 'cvxp_24_14_6': 0.3201215999999931, 'cvxp_24_14_9': 0.5004347999999936, 'cvxp_24_14_12': 0.6996279999999899, 'cvxp_24_14_15': 0.8925992000000065, 'cvxp_24_15_6': 0.3179606000000206, 'cvxp_24_15_9': 0.5061680000000024, 'cvxp_24_15_12': 0.7104605000000106, 'cvxp_24_15_15': 0.8897731000000135, 'cvxp_24_21_6': 0.24933329999998932, 'cvxp_24_21_9': 0.4092314000000101, 'cvxp_24_21_12': 0.5909948999999983, 'cvxp_24_21_15': 0.7070671999999831, 'cvxp_24_22_6': 0.25243950000000837, 'cvxp_24_22_9': 0.41057470000001217, 'cvxp_24_22_12': 0.5503779000000009, 'cvxp_24_22_15': 0.6917126000000167, 'cvxp_24_23_6': 0.2467551999999955, 'cvxp_24_23_9': 0.4032903999999746, 'cvxp_24_23_12': 0.5468122999999991, 'cvxp_24_23_15': 0.7182433000000117, 'cvxp_24_24_6': 0.26206879999998023, 'cvxp_24_24_9': 0.41221810000001824, 'cvxp_24_24_12': 0.5810053999999809, 'cvxp_24_24_15': 0.7258867999999836, 'cvxp_24_25_6': 0.3604411999999968, 'cvxp_24_25_9': 0.5872032000000047, 'cvxp_24_25_12': 0.8117747999999949, 'cvxp_24_25_15': 1.0018266000000153, 'cvxp_26_1_6': 0.3474723999999867, 'cvxp_26_1_9': 0.5589129999999898, 'cvxp_26_1_12': 0.731170800000001, 'cvxp_26_1_15': 0.9149987000000124, 'cvxp_26_2_6': 0.2676232000000027, 'cvxp_26_2_9': 0.4372754000000043, 'cvxp_26_2_12': 0.5875181999999768, 'cvxp_26_2_15': 0.7934904999999901, 'cvxp_26_3_6': 0.31557639999999765, 'cvxp_26_3_9': 0.5082222000000058, 'cvxp_26_3_12': 0.6646360000000016, 'cvxp_26_3_15': 0.8437718000000132, 'cvxp_26_4_6': 0.3278639000000112, 'cvxp_26_4_9': 0.5249557999999865, 'cvxp_26_4_12': 0.7292897999999752, 'cvxp_26_4_15': 0.8686229999999853, 'cvxp_26_5_6': 0.39229950000000713, 'cvxp_26_5_9': 0.5976344999999981, 'cvxp_26_5_12': 0.8333937999999819, 'cvxp_26_5_15': 1.0406068000000062, 'cvxp_26_11_6': 0.2709983999999963, 'cvxp_26_11_9': 0.4371241999999995, 'cvxp_26_11_12': 0.6210851999999818, 'cvxp_26_11_15': 0.7893706000000122, 'cvxp_26_12_6': 0.2404214999999965, 'cvxp_26_12_9': 0.3845597999999768, 'cvxp_26_12_12': 0.564783900000009, 'cvxp_26_12_15': 0.7138758000000109, 'cvxp_26_13_6': 0.3052853000000084, 'cvxp_26_13_9': 0.4447747000000106, 'cvxp_26_13_12': 0.6715688000000171, 'cvxp_26_13_15': 0.7831562000000076, 'cvxp_26_14_6': 0.2576035999999817, 'cvxp_26_14_9': 0.4577547000000095, 'cvxp_26_14_12': 0.6273985999999923, 'cvxp_26_14_15': 0.8036075999999923, 'cvxp_26_15_6': 0.32545190000001867, 'cvxp_26_15_9': 0.5361288000000002, 'cvxp_26_15_12': 0.730468099999996, 'cvxp_26_15_15': 0.9253480000000138, 'cvxp_26_21_6': 0.31815890000001446, 'cvxp_26_21_9': 0.5132897000000014, 'cvxp_26_21_12': 0.7073099000000127, 'cvxp_26_21_15': 0.8638214000000062, 'cvxp_26_22_6': 0.3861769999999751, 'cvxp_26_22_9': 0.6129234999999937, 'cvxp_26_22_12': 0.8451690999999926, 'cvxp_26_22_15': 1.020910299999997, 'cvxp_26_23_6': 0.27661179999998353, 'cvxp_26_23_9': 0.44167849999999476, 'cvxp_26_23_12': 0.5712737000000061, 'cvxp_26_23_15': 0.747193899999985, 'cvxp_26_24_6': 0.2686783000000048, 'cvxp_26_24_9': 0.47645740000001524, 'cvxp_26_24_12': 0.63576070000002, 'cvxp_26_24_15': 0.7769346999999982, 'cvxp_26_25_6': 0.25959349999999404, 'cvxp_26_25_9': 0.42763610000000085, 'cvxp_26_25_12': 0.5868662000000029, 'cvxp_26_25_15': 0.7426750999999854, 'cvxp_28_1_6': 0.3671327000000002, 'cvxp_28_1_9': 0.5648051999999999, 'cvxp_28_1_12': 0.7808045000000003, 'cvxp_28_1_15': 1.0089396000000002, 'cvxp_28_2_6': 0.3758400000000002, 'cvxp_28_2_9': 0.6142614000000002, 'cvxp_28_2_12': 0.8420207, 'cvxp_28_2_15': 1.0591285, 'cvxp_28_3_6': 0.32806219999999975, 'cvxp_28_3_9': 0.5000837000000011, 'cvxp_28_3_12': 0.7024057000000017, 'cvxp_28_3_15': 0.8837887999999996, 'cvxp_28_4_6': 0.29125070000000086, 'cvxp_28_4_9': 0.4777553000000001, 'cvxp_28_4_12': 0.6330477999999999, 'cvxp_28_4_15': 0.821916100000001, 'cvxp_28_5_6': 0.25004329999999975, 'cvxp_28_5_9': 0.4438441000000015, 'cvxp_28_5_12': 0.5932241999999999, 'cvxp_28_5_15': 0.7349683000000002, 'cvxp_28_11_6': 0.3359053000000003, 'cvxp_28_11_9': 0.5068840999999988, 'cvxp_28_11_12': 0.7113095000000005, 'cvxp_28_11_15': 0.903092400000002, 'cvxp_28_12_6': 0.38677960000000056, 'cvxp_28_12_9': 0.6185025999999993, 'cvxp_28_12_12': 0.8325735000000023, 'cvxp_28_12_15': 1.0871762999999994, 'cvxp_28_13_6': 0.3537575000000004, 'cvxp_28_13_9': 0.5856705999999967, 'cvxp_28_13_12': 0.7958914999999998, 'cvxp_28_13_15': 0.9863622999999997, 'cvxp_28_14_6': 0.23669860000000043, 'cvxp_28_14_9': 0.4060522000000013, 'cvxp_28_14_12': 0.5247757000000028, 'cvxp_28_14_15': 0.6933924000000005, 'cvxp_28_15_6': 0.29309820000000286, 'cvxp_28_15_9': 0.48387650000000093, 'cvxp_28_15_12': 0.6422154000000013, 'cvxp_28_15_15': 0.8456987999999974, 'cvxp_28_21_6': 0.3496210000000026, 'cvxp_28_21_9': 0.5396942999999972, 'cvxp_28_21_12': 0.7558570000000024, 'cvxp_28_21_15': 0.9728233999999993, 'cvxp_28_22_6': 0.32666650000000175, 'cvxp_28_22_9': 0.5198384999999988, 'cvxp_28_22_12': 0.6997319999999974, 'cvxp_28_22_15': 0.9013106000000022, 'cvxp_28_23_6': 0.2832040000000049, 'cvxp_28_23_9': 0.494128400000001, 'cvxp_28_23_12': 0.6520758000000058, 'cvxp_28_23_15': 0.8378298000000015, 'cvxp_28_24_6': 0.32816940000000017, 'cvxp_28_24_9': 0.5247332, 'cvxp_28_24_12': 0.7152189000000035, 'cvxp_28_24_15': 0.9220611000000005, 'cvxp_28_25_6': 0.23900739999999843, 'cvxp_28_25_9': 0.4165464999999955, 'cvxp_28_25_12': 0.5501908000000029, 'cvxp_28_25_15': 0.6859202999999994, 'cvxp_30_1_6': 0.38655109999999837, 'cvxp_30_1_9': 0.6513360000000006, 'cvxp_30_1_12': 0.8564316000000005, 'cvxp_30_1_15': 1.1143701999999962, 'cvxp_30_2_6': 0.31955159999999694, 'cvxp_30_2_9': 0.5148260000000064, 'cvxp_30_2_12': 0.707427200000005, 'cvxp_30_2_15': 0.9097918000000007, 'cvxp_30_3_6': 0.3744337000000044, 'cvxp_30_3_9': 0.5967374000000021, 'cvxp_30_3_12': 0.7841239999999985, 'cvxp_30_3_15': 0.9826465000000013, 'cvxp_30_4_6': 0.34114909999999554, 'cvxp_30_4_9': 0.5951872999999992, 'cvxp_30_4_12': 0.8048335000000009, 'cvxp_30_4_15': 1.0574974999999966, 'cvxp_30_5_6': 0.385088599999996, 'cvxp_30_5_9': 0.6149442000000036, 'cvxp_30_5_12': 0.8572083999999975, 'cvxp_30_5_15': 1.0985749999999967, 'cvxp_30_11_6': 0.3695738999999989, 'cvxp_30_11_9': 0.5748251000000053, 'cvxp_30_11_12': 0.8109503000000018, 'cvxp_30_11_15': 1.0002588000000046, 'cvxp_30_12_6': 0.3893872000000016, 'cvxp_30_12_9': 0.6110350000000011, 'cvxp_30_12_12': 0.8460971000000015, 'cvxp_30_12_15': 1.0488445000000013, 'cvxp_30_13_6': 0.3381036000000037, 'cvxp_30_13_9': 0.5121451999999991, 'cvxp_30_13_12': 0.699727799999998, 'cvxp_30_13_15': 0.9058762000000016, 'cvxp_30_14_6': 0.42660550000000086, 'cvxp_30_14_9': 0.6499140000000025, 'cvxp_30_14_12': 0.9230649, 'cvxp_30_14_15': 1.1428075999999976, 'cvxp_30_15_6': 0.39017509999999334, 'cvxp_30_15_9': 0.6283399999999943, 'cvxp_30_15_12': 0.879097299999998, 'cvxp_30_15_15': 1.0901097000000135, 'cvxp_30_21_6': 0.3470628999999974, 'cvxp_30_21_9': 0.4947786999999977, 'cvxp_30_21_12': 0.7452897000000007, 'cvxp_30_21_15': 0.9204871999999966, 'cvxp_30_22_6': 0.39367539999999224, 'cvxp_30_22_9': 0.6226434000000012, 'cvxp_30_22_12': 0.8720157999999998, 'cvxp_30_22_15': 1.057648499999999, 'cvxp_30_23_6': 0.3901735999999971, 'cvxp_30_23_9': 0.6425093000000004, 'cvxp_30_23_12': 0.8586162999999942, 'cvxp_30_23_15': 1.1320070999999956, 'cvxp_30_24_6': 0.36108159999999145, 'cvxp_30_24_9': 0.5713279, 'cvxp_30_24_12': 0.8117275000000035, 'cvxp_30_24_15': 1.0530251999999933, 'cvxp_30_25_6': 0.3517772000000008, 'cvxp_30_25_9': 0.5288498999999973, 'cvxp_30_25_12': 0.7569031000000024, 'cvxp_30_25_15': 0.9213778000000019}
	newRet = []
	RetForStats = []
	for el in RET:
		nb_cvxp = re.split('_', el[0])[2]
		
		instance_idx = re.split('_', el[0])[3]
		singelrow = [el[0], 'cvxps_size', el[1], 
					'cutoff-w-sep-10', el[3], 'time', TIME['cvxp_' + nb_cvxp + '_' + instance_idx + '_6'],
					'cutoff-w-sep-16', el[6], 'time', TIME['cvxp_' + nb_cvxp + '_' + instance_idx + '_9'],
					'cutoff-w-sep-22', el[9], 'time', TIME['cvxp_' + nb_cvxp + '_' + instance_idx + '_12'],
					'cutoff-w-sep-28', el[12],'time', TIME['cvxp_' + nb_cvxp + '_' + instance_idx + '_15']]
		singelrow2 = [el[0], float(el[1]), 
					float(el[3])/float(el[1]), float(TIME['cvxp_' + nb_cvxp + '_' + instance_idx + '_6']),
					float(el[6])/float(el[1]), float(TIME['cvxp_' + nb_cvxp + '_' + instance_idx + '_9']),
					float(el[9])/float(el[1]), float(TIME['cvxp_' + nb_cvxp + '_' + instance_idx + '_12']),
					float(el[12])/float(el[1]),float(TIME['cvxp_' + nb_cvxp + '_' + instance_idx + '_15'])]
		
		newRet.append(singelrow)
		RetForStats.append(singelrow2)
		# print(singelrow)
	df = pd.DataFrame(newRet)
	df.to_csv('stats_separators.csv', sep=',')


	# newRetStats = []
	# k = 1
	# for ridx in range(0, 165, 5):
	# 	nb_cvxp = re.split('_', RetForStats[ridx][0])[2]
		
		
	# 	if k%3 == 1:
	# 		avg = (np.array(RetForStats[ridx][1::]) + np.array(RetForStats[ridx+1][1::])+ np.array(RetForStats[ridx+2][1::])
	# 	 + np.array(RetForStats[ridx+3][1::]) + np.array(RetForStats[ridx+4][1::]))/5.0
	# 		avglist = list(avg)
	# 		avglist.insert(0,nb_cvxp)
	# 		avglist.insert(1,'easy')
	# 	if k%3 == 2:
	# 		avg = (np.array(RetForStats[ridx][1::]) + np.array(RetForStats[ridx+1][1::])+ np.array(RetForStats[ridx+2][1::])
	# 	 + np.array(RetForStats[ridx+3][1::]) + np.array(RetForStats[ridx+4][1::]))/5.0
	# 		avglist = list(avg)
	# 		avglist.insert(0,nb_cvxp)
	# 		avglist.insert(1,'medium')
	# 	if k%3 == 0:
	# 		avg = (np.array(RetForStats[ridx][1::]) + np.array(RetForStats[ridx+1][1::])+ np.array(RetForStats[ridx+2][1::])
	# 	 + np.array(RetForStats[ridx+3][1::]) + np.array(RetForStats[ridx+4][1::]))/5.0
	# 		avglist = list(avg)
	# 		avglist.insert(0,nb_cvxp)
	# 		avglist.insert(1,'hard')
	# 	k += 1
	# 	newRetStats.append(avglist)
	# df = pd.DataFrame(newRetStats)
	# df.to_csv('stats_separators_agg.csv', sep=',')
	# print(newRetStats)
if __name__ == "__main__":
	# mode = 'separator-result-integration'

	# if mode == 'separator-result-integration':
	# 	print("separator-result-integration running ...")
	# 	get_statistics_separators()

	# elif mode == '1-1-result-integration':
	# 	get_statistics_mode11()

	# else:
	# 	print("wrong mode input")
	# get_statistics_mode11()
	get_statistics_mode11()