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
import time
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
			R.append(1)

	file_.close()
	return depot, X, Y, R

def read_Mennell_instance(instancefile):
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
			R.append(float(list_[3])* 1)

	file_.close()
	return depot, X, Y, R


def  solve_MasterProb(depot, Ox, Oy, Or, zbar, gbc_set):
	nb_tars = len(Ox)
	try:
		# Create a new model
		model = gp.Model("GBD_CETSP")
		model.setParam(GRB.Param.OutputFlag, 0)
		model.setParam(GRB.Param.TimeLimit, 15000.0)
		model._nb_GBCs = 0
		model._size = nb_tars + 2
		model._nb_SECs = 0
		model._gbc_set = gbc_set
		model._gbc_set_added = False
		# binary variables
		varE = model.addVars(nb_tars+2, nb_tars+2, vtype=GRB.BINARY, name="E")
		theta = model.addVar(lb=0, vtype=GRB.CONTINUOUS,name="theta")
		obj = 0
		for i in range(nb_tars+2):
			for j in range(nb_tars+2):
				obj += zbar[i,j] * varE[i,j] 
		obj += 1.0 * theta
		model.setObjective(obj, GRB.MINIMIZE)
		
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

		if len(model._gbc_set) != 0:
				# model._gbc_set_added = True
				for gbc in gbc_set:
					obj = gbc[0]
					mu1 = gbc[1]
					constr = 0
					for el in mu1:
						constr += el[2] * (varE[el[0],el[1]] - 1.0)
					model.addConstr(obj + constr <= theta)
					# model._nb_GBCs += 1
		
		model._varE = varE
		model._theta = theta

		model.Params.lazyConstraints = 1
		model.optimize(callback_subtour)
		if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
			curObj = np.zeros((model._size, model._size))
			tmp = 0
			for i in range(0, model._size):
				for j in range(0, model._size):
					curObj[i,j] = varE[i,j].x
					tmp += curObj[i,j] * zbar[i][j]
					# print(varE[i,j].x, end=' ')
				# print('\n',end='')
			objval = model.objVal
			# print('theta', theta.x, 'tmp', tmp, 'model obj', objval)
			cirset, curBestSeq =  find_all_subtours(curObj, model._size)
			optTour, curtour_Len = GBCuts.solve_SOCP_Disc(depot, Ox, Oy, Or, curBestSeq)
			print(curBestSeq,curtour_Len )

			return curBestSeq, tmp, objval
	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))
	except AttributeError:
		print('Encountered an attribute error')


def callback_subtour(model, where):
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
			# for circle in Cirset:
			smalllens = float('inf')
			smallindex = -1

			for i in range(len(Cirset)):
				if len(Cirset[i]) < smalllens:
					smalllens = len(Cirset[i])
					smallindex = i
			circle = Cirset[smallindex]	
			constr = 0
			for i in range(0,len(circle)-1):
				constr += model._varE[circle[i],circle[i+1]]
			constr += model._varE[circle[-1],circle[0]]
			model.cbLazy(constr <= len(circle)-1)
			model._nb_SECs  += 1
		# else:
		# 	if not model._gbc_set_added and  len(model._gbc_set) != 0:
		# 		model._gbc_set_added = True
		# 		for gbc in model._gbc_set:
		# 			obj = gbc[0]
		# 			mu1 = gbc[1]
		# 			constr = 0
		# 			for el in mu1:
		# 				constr += el[2] * (model._varE[el[0],el[1]] - 1.0)
		# 			model.cbLazy(obj + constr <= model._theta)
		# 			model._nb_GBCs += 1


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


def generate_GenOptimalityCut(fseq, size, depot, Ox, Oy, Or, zbar):
	# print(fseq)
	depot1 = depot
	depot2 = depot
	nb_tars =  size - 2
	vismat = np.zeros(( size, size))
	for i in range(len(fseq)-1):
		vismat[fseq[i],fseq[i+1]] = 1
	# print(vismat)
	try:
		SP_m = gp.Model("SOCP")
		SP_m.setParam(GRB.Param.OutputFlag, 0)
		SP_m.setParam(GRB.Param.TimeLimit, 1000.0)

		varZ = SP_m.addVars( size,  size, lb=0,ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Z")

		varX = SP_m.addVars( size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="X")
		varY = SP_m.addVars( size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Y")

		varS = SP_m.addVars( size,  size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="S")
		varT = SP_m.addVars( size,  size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="T")

		varK = SP_m.addVars( size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="K")
		varL = SP_m.addVars( size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="L")

		SP_m.update()

		# objective function
		obj = 0
		for i in range( size):
			for j in range( size):
				obj += vismat[i,j] * (varZ[i,j] -  zbar[i,j])
		SP_m.setObjective(obj, GRB.MINIMIZE)

		# constraints - depot1
		for i in range( size):
			for j in range(i+1,  size):
				SP_m.addConstr(varS[j,i]==varS[i,j])
				SP_m.addConstr(varT[j,i]==varT[i,j])
				SP_m.addConstr(varZ[j,i]==varZ[i,j])
			SP_m.update()

		for i in range(0, nb_tars):
			SP_m.addConstr(varS[0,i+1] ==  depot1[0]-varX[i], 'cS_'+str(0)+str(i+1))
			SP_m.addConstr(varS[i+1,0]==varS[0,i+1])
			SP_m.addConstr(varT[0,i+1] ==  depot1[1]-varY[i], 'cT_'+str(0)+str(i+1))
			SP_m.addConstr(varT[i+1,0]==varT[0,i+1])

		for i in range(0, nb_tars):
			SP_m.addConstr(varS[i+1,nb_tars+1] == varX[i] -  depot2[0], 'cS_'+str(i+1)+str(nb_tars+1))
			SP_m.addConstr(varS[nb_tars+1,i+1]==varS[i+1,nb_tars+1])
			SP_m.addConstr(varT[i+1,nb_tars+1] == varY[i] -  depot2[1], 'cT_'+str(i+1)+str(nb_tars+1))
			SP_m.addConstr(varT[nb_tars+1,i+1]==varT[i+1,nb_tars+1])

		for i in range(1, nb_tars+1):
			for j in range(i+1, nb_tars+1):
				SP_m.addConstr(varS[i,j] == varX[i-1]-varX[j-1], 'cS_'+str(i)+str(j))
				SP_m.addConstr(varS[j,i] == varS[i,j])
				SP_m.addConstr(varT[i,j] == varY[i-1]-varY[j-1], 'cT_'+str(i)+str(j))
				SP_m.addConstr(varT[j,i] == varT[i,j])
			SP_m.update()
		
		for i in range( size):
			for j in range( size):
				if i!=j:
					SP_m.addConstr(varS[i,j]*varS[i,j] + varT[i,j]*varT[i,j] <= varZ[i,j]*varZ[i,j], 'cZ_'+str(i)+str(j)) #constraint S
			SP_m.update()

		for i in range(0, nb_tars):
			SP_m.addConstr(varX[i]- Ox[i] == varK[i], 'cK_'+str(i)) 
			SP_m.addConstr(varY[i]- Oy[i] == varL[i], 'cL_'+str(i)) 
		SP_m.update()

		for i in range(0, nb_tars):
			SP_m.addConstr(varK[i]*varK[i] + varL[i]*varL[i] <=  Or[i]**2, 'cR_'+str(i)) 
		SP_m.update()


		SP_m.optimize()

		''' ------------- model output  ----------------------'''
		if SP_m.status == GRB.OPTIMAL or SP_m.status == GRB.TIME_LIMIT:
			# print('here2')
			mu1 = []
			for i in range( size):
				for j in range( size):
					if i != j and vismat[i,j] > 0:
						# mu1.append([i,j,varZ[i,j].x])
						mu1.append([i,j,varZ[i,j].x -  zbar[i,j]])
						# mu1.append([j,i,varZ[i,j].x -  zbar[i,j]])					


			mu0 = []
			for i in range(nb_tars):
				for j in range(i+2, nb_tars+1):
					if i == 0:
						val =  math.sqrt(( depot1[0] - varX[fseq[j]-1].x)**2 + ( depot1[1]-varY[fseq[j]-1].x)**2)
						val = val -  zbar[0,fseq[j]]
						mu0.append([0, fseq[j],val])
						# mu0.append([fseq[j], 0,val])

					if i != 0:
						val =  math.sqrt((varX[fseq[i]-1].x - varX[fseq[j]-1].x)**2 + (varY[fseq[i]-1].x-varY[fseq[j]-1].x)**2)	
						val = val -  zbar[fseq[i],fseq[j]]
						mu0.append([fseq[i], fseq[j], val])	
						# mu0.append([fseq[j], fseq[i], val])

			for i in range(1, nb_tars):
				val =  math.sqrt(( depot2[0] - varX[fseq[i]-1].x)**2 + ( depot2[1]-varY[fseq[i]-1].x)**2)	
				val = val -  zbar[fseq[i],nb_tars+1]
				mu0.append([fseq[i], nb_tars+1, val])
				# mu0.append([nb_tars+1, fseq[i], val])

			return mu0, mu1, SP_m.objVal

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))
	except AttributeError:
		print('Encountered an attribute error')

def main_loop(instance):
	depot, Ox, Oy, Or = read_instance_Behdani(instance)
	# depot, Ox, Oy, Or = read_Mennell_instance(instance)

	# separators = CNgb.find_separators(depot, Ox, Oy, Or)
	''' set coefs for master problem '''
	nb_tars = len(Ox)
	size = nb_tars + 2
	zbar = np.zeros((nb_tars+2, nb_tars+2))
	for i in range(nb_tars):
		dist = max([math.sqrt((Ox[i]-depot[0])**2 + (Oy[i]-depot[1])**2) -Or[i],0])
		# dist = math.sqrt((Ox[i]-depot[0])**2 + (Oy[i]-depot[1])**2) 
		zbar[0,i+1] = dist
		zbar[i+1,0] = zbar[0,i+1]
		zbar[i+1,nb_tars+1] = zbar[0,i+1]
		zbar[nb_tars+1,i+1] = zbar[i+1,nb_tars+1]
	for i in range(nb_tars):
		for j in range(nb_tars):
			if i != j:
				dist = max([math.sqrt((Ox[i]-Ox[j])**2 + (Oy[i]-Oy[j])**2) - Or[i] - Or[j],0])
				# dist = math.sqrt((Ox[i]-Ox[j])**2 + (Oy[i]-Oy[j])**2)
				zbar[i+1,j+1] = dist
				zbar[j+1,i+1] = zbar[i+1,j+1]
	
	UB = float('inf')
	LB = 0.000001
	tolerance = 0.001
	gbc_set = []
	itr = 0
	start_time = time.time()

	# local_search_seqs(nb_tars, depot, Ox, Oy, Or)
	# exit(1)
	while UB - LB > tolerance:
		# print(itr, ':', UB, LB, (UB-LB)/LB * 100.0)

		fseq, tmp, modelobjvalue = solve_MasterProb(depot, Ox, Oy, Or, zbar,  gbc_set)
		mu0, mu1, obj = generate_GenOptimalityCut(fseq, size, depot, Ox, Oy, Or, zbar)

		UB = min(UB, tmp + obj)
		# print(tmp + obj)
		# print('after solving master', time.time() - start_time, 'secs')

		LB = max(LB, modelobjvalue)
		# print( '--------------', UB, LB)

		# print('after solving sub', time.time() - start_time, 'secs')

		# print('sub obj = ', obj)
		gbc_set.append([obj, mu1])

		# print('--------------------------------------------------------------------------------------')
		itr += 1
	print('Total time', time.time() - start_time, 'secs')
	
# def local_search_seqs(nb_tars, depot, Ox, Oy, Or):
# 	sorted_distmat = []
# 	tmp = {}
# 	for i in range(nb_tars):
# 		dist = math.sqrt((Ox[i]-depot[0])**2 + (Oy[i]-depot[1])**2) -Or[i]
# 		tmp[i+1] = dist
# 	sorted_tmp = {k: v for k, v in sorted(tmp.items(), key=lambda item: item[1])}
# 	sorted_distmat.append(sorted_tmp)
# 	for i in range(nb_tars):
# 		tmp = {}
# 		for j in range(nb_tars):
# 			if i != j:
# 				dist =  math.sqrt((Ox[i]-Ox[j])**2 + (Oy[i]-Oy[j])**2) - Or[i] - Or[j]
# 				tmp[j+1] = dist
# 		sorted_tmp = {k: v for k, v in sorted(tmp.items(), key=lambda item: item[1])}
# 		sorted_distmat.append(sorted_tmp)
# 	# print(sorted_distmat)
	# for idx in range(1, nb_tars+1):


# def BFS_simplified(sorted_distmat, nb_tars):
	
# 	one_seq = []
# 	breadth_limit = 3
# 	while True:
# 		for idx in range(1, nb_tars+1):
# 			while nb_next < breadth_limit:
# 				one_seq.append()








if __name__ == "__main__":
	instance = argv[1]
	main_loop(instance)
