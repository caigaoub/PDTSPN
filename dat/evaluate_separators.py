

from InstanceReader import CvxPolygon
import timeit

def evaluate_single_set(cvp, nb_cvxp):
	path = "C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/dat/Cai2/"
	
	ret = []
	for idx in [1,2,3,4,5,10,11,12,13,14,15,20,21,22,23,24,25]:
		instance = path + "cvxp_" + str(nb_cvxp) + "_" + str(idx)
		starttime = timeit.default_timer()
		cvp.read_cvxp_instance(instance)
		sep_set= cvp.generate_separators(25)
		endtime = timeit.default_timer()
		# print(endtime - starttime)
		truepolygonarea, estimatedpolygonarea, cutoffarea = cvp.evaluate_separators(sep_set)
		res = [nb_cvxp, idx, endtime - starttime, truepolygonarea, estimatedpolygonarea, cutoffarea]
		print(res)
		ret.append(res)
		# cvp.plot_hull()
	print(ret)


if __name__ == "__main__":
	cvp = CvxPolygon('cut off test')
	evaluate_single_set(cvp, 10)
