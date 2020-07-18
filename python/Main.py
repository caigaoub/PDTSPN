import gurobipy as gp
from gurobipy import GRB
from itertools import combinations
from sys import argv
import re
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pch
import networkx as nx
#--------------------------------
import GBCuts
import CircularNeighborhoods as CNgb
import ConvexPolygonNeighborhoods as CvxPolyNgb
def read_instance_Behdani(instancefile):
	file_ = open(instancefile, 'r')
	X = []
	Y = []
	R = []
	line_ = file_.readline()
	list_ = re.split(" |\t|\n", line_)
	depot = [float(list_[0]), float(list_[1])]
	while True:
		line_ = file_.readline()
		list_ = re.split(" |\t|\n", line_)
		if list_[0] == '':
			break
		else:
			list_ = list(filter(None, list_))
			# print(list_)
			X.append(float(list_[0]))
			Y.append(float(list_[1]))
			R.append(1.5)

	file_.close()
	return depot, X, Y, R

	
def plot_Behdani_instance(depot, X, Y, R):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	# plt.axis('off')
	XX = [x for x in X]
	XX.append(depot[0])
	YY = [y for y in Y]
	YY.append(depot[1])
	ax.set_xlim([min(XX)-max(R), max(XX)+max(R)])
	ax.set_ylim([min(YY)-max(R), max(YY)+max(R)])
	ax.set_aspect('equal', adjustable='box')
	rect = pch.Rectangle((depot[0], depot[1]), 0.2, 0.2, edgecolor='k',fill=True) 
	ax.add_patch(rect) 
	for i in range(len(X)):
		circle = plt.Circle((X[i], Y[i]), radius=R[i], edgecolor='r',fill=False, alpha = 0.5)
		ax.add_artist(circle)
		ax.annotate(str(i+1), xy=(X[i], Y[i]), xytext=(X[i], Y[i]),size=11)
	plt.xticks(np.arange(min(XX)-max(R),max(XX)+max(R),0.5))
	plt.yticks(np.arange(min(YY)-max(R),max(YY)+max(R),0.5))
	# plt.grid(alpha=.5)
	plt.show()


def  solve_MasterProb(depot, Ox, Oy, Or, separators):
	nb_tars = len(Ox)
	try:
		# Create a new model
		model = gp.Model("GBD_CETSP")
		model.setParam(GRB.Param.OutputFlag, 0)
		model.setParam(GRB.Param.TimeLimit, 1500.0)

		model._depot1 = depot
		model._depot2 = depot  

		model._Ox = Ox
		model._Oy = Oy
		model._Or = Or
		model._separators = separators
		# zbar = np.full((nb_tars+1, nb_tars+1), np.inf)
		zbar = np.zeros((nb_tars+2, nb_tars+2))
		for i in range(nb_tars):
			dist = max([math.sqrt((Ox[i]-depot[0])**2 + (Oy[i]-depot[1])**2) -Or[i],0])
			zbar[0,i+1] = dist
			zbar[i+1,0] = zbar[0,i+1]
			zbar[i+1,nb_tars+1] = zbar[0,i+1]
			zbar[nb_tars+1,i+1] = zbar[i+1,nb_tars+1]

		for i in range(nb_tars):
			for j in range(nb_tars):
				if i != j:
					dist = max([math.sqrt((Ox[i]-Ox[j])**2 + (Oy[i]-Oy[j])**2) - Or[i] - Or[j],0])
					zbar[i+1,j+1] = dist
					zbar[j+1,i+1] = zbar[i+1,j+1]

		# print(zbar)
		model._zbar = zbar
		# binary variables
		M = 1000.0
		varE = model.addVars(nb_tars+2, nb_tars+2, vtype=GRB.BINARY,name="E")
		theta = model.addVar(lb=0.0,vtype=GRB.CONTINUOUS,name="theta")
		obj = 0
		for i in range(nb_tars+2):
			for j in range(nb_tars+2):
				obj += zbar[i,j] * varE[i,j] 
		obj += 1.0 * theta
		model.setObjective(obj, GRB.MINIMIZE)
		# model.modelSense = GRB.MINIMIZE
		
		# constraints
		for i in range(0, nb_tars+2):
			constr1 = 0
			constr2 = 0
			for j in range(0, nb_tars+2):
				constr1 += varE[i,j]
				constr2 += varE[j,i]
			model.addConstr(constr1 == 1, 'r_'+str(i)) 
			model.addConstr(constr2 == 1, 'c_'+str(i)) 
			model.update()
		# remove TSP sequence of inverse order
		exprA = 0
		exprB = 0
		for i in range(1,nb_tars+1):
			exprA += i * varE[0,i]
			exprB += i * varE[i,nb_tars+1]
		model.addConstr(exprA <= exprB, 'inv') 

		varE[nb_tars+1,0].lb = 1.0
		varE[0,nb_tars+1].ub = 0.0

		for i in range(0, nb_tars+2):
			varE[i,i].ub = 0.0

		# model.write('master.lp')

		model._varE = varE
		model._theta = theta
		model._size = nb_tars + 2
		model._nb_SECs = 0
		model._nb_GBCs = 0
		model.Params.lazyConstraints = 1

		# ADD_USERCUT(model)
		model.optimize(callback_GBD)
		''' ------------- model output  ----------------------'''
		if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
			curObj = np.zeros((model._size, model._size))
			for i in range(0, model._size):
				for j in range(0, model._size):
					curObj[i,j] = varE[i,j].x
					# print(varE[i,j].x, end=' ')
				# print('\n',end='')

			cirset, curBestSeq =  find_all_subtours(curObj, model._size)
			# print(curBestSeq)
			optTour, objLen = GBCuts.solve_SOCP_Disc(depot, Ox, Oy, Or, curBestSeq)
			CNgb.plot_Behdani_instance(depot, Ox, Oy, Or, optTour)
			# print('Theta: %g' % theta.x)
			# print('Gurobi Time: %g' % model.Runtime)
			# print('OBJ: %g' % model.objVal)
			print(model._nb_SECs, '& ', model._nb_GBCs, '& ', len(separators), '& ', '{:.2f}'.format(model.Runtime), '& ', '{:.2f}'.format(model.MIPGap), '& ', '{:.4f}'.format(objLen))

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))
	except AttributeError:
		print('Encountered an attribute error')



def  solve_MasterProb_CvxPolyNgbs(depot, LPC, HULLs):
	nb_tars = len(LPC)
	try:
		# Create a new model
		model = gp.Model("GBD_CETSP_CvxPolyNgbs")
		model.setParam(GRB.Param.OutputFlag, 1)
		model.setParam(GRB.Param.TimeLimit, 10.0)
		model._depot1 = depot
		model._depot2 = depot  
		model._LPC = LPC
		zbar = CvxPolyNgb.calc_zbar(depot, HULLs)
		model._zbar = zbar
		# binary variables
		varE = model.addVars(nb_tars+2, nb_tars+2, vtype=GRB.BINARY,name="E")
		theta = model.addVar(lb=0.0,vtype=GRB.CONTINUOUS,name="theta")
		obj = 0
		for i in range(nb_tars+2):
			for j in range(nb_tars+2):
				obj += zbar[i,j] * varE[i,j] 
		obj += 1.0 * theta
		model.setObjective(obj, GRB.MINIMIZE)
		# model.modelSense = GRB.MINIMIZE
		
		# constraints
		for i in range(0, nb_tars+2):
			constr1 = 0
			constr2 = 0
			for j in range(0, nb_tars+2):
				constr1 += varE[i,j]
				constr2 += varE[j,i]
			model.addConstr(constr1 == 1, 'r_'+str(i)) 
			model.addConstr(constr2 == 1, 'c_'+str(i)) 
			model.update()
		# remove TSP sequence of inverse order
		exprA = 0
		exprB = 0
		for i in range(1,nb_tars+1):
			exprA += i * varE[0,i]
			exprB += i * varE[i,nb_tars+1]
		model.addConstr(exprA <= exprB, 'inv') 

		varE[nb_tars+1,0].lb = 1.0
		varE[0,nb_tars+1].ub = 0.0

		for i in range(0, nb_tars+2):
			varE[i,i].ub = 0.0

		# model.write('master.lp')

		model._varE = varE
		model._theta = theta
		model._size = nb_tars + 2
		model._nb_SECs = 0
		model._nb_GBCs = 0
		model.Params.lazyConstraints = 1

		# ADD_USERCUT(model)
		model.optimize(callback_GBD)
		''' ------------- model output  ----------------------'''
		if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
			curObj = np.zeros((model._size, model._size))
			for i in range(0, model._size):
				for j in range(0, model._size):
					curObj[i,j] = varE[i,j].x
					# print(varE[i,j].x, end=' ')
				# print('\n',end='')

			cirset, curBestSeq =  find_all_subtours(curObj, model._size)
			# print(curBestSeq)
			optTour, objLen = GBCuts.solve_SOCP_CvxPolyNgbs(model, curBestSeq)
			# CvxPolyNgb.plot_CvxPoly_Instance_wOptTour(model._depot1, HULLs, optTour)
			# print('Theta: %g' % theta.x)
			# print('Gurobi Time: %g' % model.Runtime)
			# print('OBJ: %g' % model.objVal)
			print(model._nb_SECs, '& ', model._nb_GBCs, '& ', '{:.2f}'.format(model.Runtime), '& ', '{:.2f}'.format(model.MIPGap), '& ', '{:.4f}'.format(objLen))

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))
	except AttributeError:
		print('Encountered an attribute error')

def callback_GBD(model, where):
	if where == GRB.Callback.MIPSOL:
		vals = model.cbGetSolution(model._varE)
		if False:
			for i in range(model._size):
				for j in range(model._size):
					print(vals[i,j],end=' ')
				print(end='\n')

		# find_Subtour(vals, model._size)
		Cirset,fseq = find_all_subtours(vals, model._size)
		# print(Cirset,fseq)
		if (len(Cirset) != 0):
			for circle in Cirset:
				constr = 0
				for i in range(0,len(circle)-1):
					constr += model._varE[circle[i],circle[i+1]]
				constr += model._varE[circle[-1],circle[0]]

				model.cbLazy(constr <= len(circle)-1)
				model._nb_SECs  += 1
		else:
			# mu0, mu1, obj = GBCuts.generate_GenOptimalityCut(model, fseq)
			mu0, mu1, obj = GBCuts.generate_GBC_CvxPolyNgbs(model, fseq)

			# print(objbar,obj, objbar + obj)

			constr = 0
			for el in mu1:
				constr += el[2] * (model._varE[el[0],el[1]] - 1.0)
			# for el in mu0:
			# 	constr += el[2] * (model._varE[el[0],el[1]] - 0.0)
			model.cbLazy(obj + constr <= model._theta)
			model._nb_GBCs += 1

def find_all_subtours(curSol, matsize):
	visited = [False]*matsize
	cur_node = 0
	visited[cur_node] = True
	Cirset = []
	fseq = []
	circle = [cur_node]
	while True:
		zeronext = True
		i = 0
		while i < matsize:
			if curSol[cur_node, i] > 0.5:
				zeronext = False
			if curSol[cur_node, i] > 0.5 and visited[i] == False:
				circle.append(i)
				visited[i] = True
				cur_node = i
				break
			i += 1			
		if i == matsize: # could be a cycle or a path. 
			# print(circle)
			if zeronext != True and len(circle) > 1 and len(circle) < matsize:
				Cirset.append(circle)
			if len(circle) == matsize:
				fseq = circle
			j = 0
			while j < matsize: # start from another unvisited node
				if visited[j] == False:
					cur_node = j
					visited[j] = True
					circle = [cur_node]
					break
				j += 1
			if j == matsize: # no unvisited node left; terminate
				break
	
	return Cirset,fseq
 		
def find_all_subcomponents(curSol, matsize):
	# visited = [False]*matsize
	# cur_node = 0
	# visited[cur_node] = True
	G = nx.Graph()
	for i in range(matsize):
		for j in range(matsize):
			if curSol[i,j] > 0:
				G.add_edge(i, j, weight=curSol[i,j])
				G.add_edge(j, i, weight=curSol[i,j])


	# list(nx.connected_components(G))
	all_subcomps = []
	for compo in list(nx.connected_components(G)):
		emptylst = []
		for e in compo:
			emptylst.append(e)
		all_subcomps.append(emptylst)
	return all_subcomps

if __name__ == "__main__":
	instance = argv[1]
	''' discs '''
	# depot, Ox, Oy, Or = read_instance_Behdani(instance)
	# # plot_Behdani_instance(depot, Ox, Oy, Or)
	# separators = CNgb.find_separators(depot, Ox, Oy, Or)
	# solve_MasterProb(depot, Ox, Oy, Or, separators)

	''' Convex polygon Neighborhoods '''
	depot, HULLs = CvxPolyNgb.read_cvxp_instance(instance)
	LPC = CvxPolyNgb.convert_HULL_LPConsts(HULLs)
	solve_MasterProb_CvxPolyNgbs(depot, LPC, HULLs)

	# for i in range(15,21):
	# 	instance = '/home/latte/Dropbox/Box_Research/Github/CETSP/dat/Cai/cvxp_18_' + str(i)
	# 	depot, HULLs = CvxPolyNgb.read_cvxp_instance(instance)
	# 	LPC = CvxPolyNgb.convert_HULL_LPConsts(HULLs)
	# 	solve_MasterProb_CvxPolyNgbs(depot, LPC, HULLs)
	# print('------------------------------------------------------')

	# for i in range(14,21):
	# 	instance = '/home/latte/Dropbox/Box_Research/Github/CETSP/dat/Cai/cvxp_20_' + str(i)
	# 	depot, HULLs = CvxPolyNgb.read_cvxp_instance(instance)
	# 	LPC = CvxPolyNgb.convert_HULL_LPConsts(HULLs)
	# 	solve_MasterProb_CvxPolyNgbs(depot, LPC, HULLs)












