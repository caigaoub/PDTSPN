import gurobipy as gp
from gurobipy import GRB
from itertools import combinations

from sys import argv
import re
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pch
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, PathPatch
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
import mpl_toolkits.mplot3d.art3d as art3d



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
			R.append(0.5)

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
	plt.grid(alpha=.5)
	plt.show()


def  solve_MasterProb(depot, Ox, Oy, Or ):
	nb_tars = len(Ox)
	try:
		# Create a new model
		model = gp.Model("GBD_CETSP")
		model.setParam(GRB.Param.OutputFlag, 1)
		model.setParam(GRB.Param.TimeLimit, 1000.0)

		model._nb_SECs = 0
		model._depot1 = depot
		model._depot2 = depot  

		model._Ox = Ox
		model._Oy = Oy
		model._Or = Or
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
		model.Params.lazyConstraints = 1
		# ADD_USERCUT(model)
		model.optimize(cb_GenBendersCuts)
		''' ------------- model output  ----------------------'''
		if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
			for i in range(0, model._size):
				for j in range(0, model._size):
					print(int(varE[i,j].x), end=' ')
				print('\n',end='')
			print('Theta: %g' % theta.x)
			print('Gurobi Time: %g' % model.Runtime)
			print('OBJ: %g' % model.objVal)

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))
	except AttributeError:
		print('Encountered an attribute error')

def cb_GenBendersCuts(model, where):
	if where == GRB.Callback.MIPSOL:
		vals = model.cbGetSolution(model._varE)
		if False:
			for i in range(model._size):
				for j in range(model._size):
					print(vals[i,j],end=' ')
				print('\n')

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
			mu0, mu1, obj = generate_GenOptimalityCut(model, fseq)
			# print(objbar,obj, objbar + obj)

			constr = 0
			for el in mu1:
				constr += el[2] * (model._varE[el[0],el[1]] - 1.0)
			# for el in mu0:
			# 	constr += el[2] * (model._varE[el[0],el[1]] - 0.0)
			model.cbLazy(obj + constr <= model._theta)

# def ADD_USERCUT(model):
# 	for i in range(model._size):
# 		for j in range(model._size):
# 			model._varE[i,j].vtype = GRB.CONTINUOUS
# 	model.optimize(cb_fractionalSolCut)
# 	for i in range(model._size):
# 		for j in range(model._size):
# 			model._varE[i,j].vtype = GRB.BINARY

# def cb_fractionalSolCut(model, where):
# 	if where == GRB.Callback.MIPNODE:
# 		print("************************")

# 		fracSol = model.cbGetSolution(model._varE)
# 		if False:
# 			for i in range(model._size):
# 				for j in range(model._size):
# 					print(fracSol[i,j],end=' ')
# 				print('\n')

# 		# find_Subtour(vals, model._size)
# 		Cirset,fseq = find_all_subcomponents(fracSol, model._size)
# 		# print(Cirset,fseq)
# 		print("************************")

# 		if (len(Cirset) != 0):
# 			for circle in Cirset:
# 				constr = 0
# 				for i in range(0,len(circle)-1):
# 					constr += model._varE[circle[i],circle[i+1]]
# 				constr += model._varE[circle[-1],circle[0]]
# 				model.cbCut(constr <= len(circle)-1)
# 				model._nb_SECs  += 1
# 				print("************************")

# 		else:
# 			mu0, mu1, obj = generate_fractionalSolCut(model, fracSol)
# 			# print(objbar,obj, objbar + obj)
# 			constr = 0
# 			for el in mu1:
# 				constr += el[2] * (model._varE[el[0],el[1]] - 1.0)
# 			# for el in mu0:
# 			# 	constr += el[2] * (model._varE[el[0],el[1]] - 0.0)
# 			model.cbCut(obj + constr <= model._theta)
# 			print("*&&&&&&&&&&&&&&&&&&&&&&&&")

# def generate_fractionalSolCut(model, fracSol):
# 	nb_tars = model._size - 2
# 	try:
# 		SP_m = gp.Model("SOCP")
# 		SP_m.setParam(GRB.Param.OutputFlag, 0)
# 		SP_m.setParam(GRB.Param.TimeLimit, 1000.0)

# 		varZ = SP_m.addVars(model._size, model._size, lb=0.001,ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Z")

# 		varX = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="X")
# 		varY = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Y")

# 		varS = SP_m.addVars(model._size, model._size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="S")
# 		varT = SP_m.addVars(model._size, model._size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="T")

# 		varK = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="K")
# 		varL = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="L")

# 		SP_m.update()

# 		# objective function
# 		obj = 0
# 		for i in range(model._size):
# 			for j in range(model._size):
# 				obj += fracSol[i,j] * (varZ[i,j] - model._zbar[i,j])
# 		SP_m.setObjective(obj, GRB.MINIMIZE)

# 		# constraints - depot1
# 		for i in range(model._size):
# 			for j in range(i+1, model._size):
# 				SP_m.addConstr(varS[j,i]==varS[i,j])
# 				SP_m.addConstr(varT[j,i]==varT[i,j])
# 				SP_m.addConstr(varZ[j,i]==varZ[i,j])
# 			SP_m.update()

# 		for i in range(0, nb_tars):
# 			SP_m.addConstr(varS[0,i+1] == model._depot1[0]-varX[i], 'cS_'+str(0)+str(i+1))
# 			SP_m.addConstr(varS[i+1,0]==varS[0,i+1])
# 			SP_m.addConstr(varT[0,i+1] == model._depot1[1]-varY[i], 'cT_'+str(0)+str(i+1))
# 			SP_m.addConstr(varT[i+1,0]==varT[0,i+1])

# 		for i in range(0, nb_tars):
# 			SP_m.addConstr(varS[i+1,nb_tars+1] == varX[i] - model._depot2[0], 'cS_'+str(i+1)+str(nb_tars+1))
# 			SP_m.addConstr(varS[nb_tars+1,i+1]==varS[i+1,nb_tars+1])
# 			SP_m.addConstr(varT[i+1,nb_tars+1] == varY[i] - model._depot2[1], 'cT_'+str(i+1)+str(nb_tars+1))
# 			SP_m.addConstr(varT[nb_tars+1,i+1]==varT[i+1,nb_tars+1])

# 		for i in range(1, nb_tars+1):
# 			for j in range(i+1, nb_tars+1):
# 				SP_m.addConstr(varS[i,j] == varX[i-1]-varX[j-1], 'cS_'+str(i)+str(j))
# 				SP_m.addConstr(varS[j,i] == varS[i,j])
# 				SP_m.addConstr(varT[i,j] == varY[i-1]-varY[j-1], 'cT_'+str(i)+str(j))
# 				SP_m.addConstr(varT[j,i] == varT[i,j])
# 			SP_m.update()
		
# 		for i in range(model._size):
# 			for j in range(model._size):
# 				if i!=j:
# 					SP_m.addConstr(varS[i,j]*varS[i,j] + varT[i,j]*varT[i,j] <= varZ[i,j]*varZ[i,j], 'cZ_'+str(i)+str(j)) #constraint S
# 			SP_m.update()

# 		for i in range(0, nb_tars):
# 			SP_m.addConstr(varX[i]-model._Ox[i] == varK[i], 'cK_'+str(i)) 
# 			SP_m.addConstr(varY[i]-model._Oy[i] == varL[i], 'cL_'+str(i)) 
# 		SP_m.update()

# 		for i in range(0, nb_tars):
# 			SP_m.addConstr(varK[i]*varK[i] + varL[i]*varL[i] <= model._Or[i]**2, 'cR_'+str(i)) 
# 		SP_m.update()

# 		# SP_m.write('socp.lp')
# 		SP_m.optimize()

# 		''' ------------- model output  ----------------------'''
# 		if SP_m.status == GRB.OPTIMAL or SP_m.status == GRB.TIME_LIMIT:
		
# 			mu1 = []
# 			for i in range(model._size):
# 				for j in range(model._size):
# 					if i!=j and fracSol[i,j] > 0:
# 						mu1.append([i,j,varZ[i,j].x])					
# 			mu0 = []
# 			# for i in range(len(fseq)-2):
# 			# 	for j in range(i+2, len(fseq)-1):
# 			# 		if i == 0:
# 			# 			val =  math.sqrt((model._depot1[0] - varX[fseq[j]-1].x)**2 + (model._depot1[1]-varY[fseq[j]-1].x)**2)
# 			# 			val = val - model._zbar[0,fseq[j]]
# 			# 			mu0.append([0, fseq[j],val])
# 			# 		if i != 0:
# 			# 			val =  math.sqrt((varX[fseq[i]-1].x - varX[fseq[j]-1].x)**2 + (varY[fseq[i]-1].x-varY[fseq[j]-1].x)**2)	
# 			# 			val = val - model._zbar[fseq[i],fseq[j]]
# 			# 			mu0.append([fseq[i], fseq[j], val])	


# 			# for i in range(1, len(fseq)-2):
# 			# 	val =  math.sqrt((model._depot2[0] - varX[fseq[i]-1].x)**2 + (model._depot2[1]-varY[fseq[i]-1].x)**2)	
# 			# 	val = val - model._zbar[fseq[i],model._size-1]
# 			# 	mu0.append([fseq[i], model._size-1,val])

# 			# print(mu1)
# 			# print(mu0)
# 			# for i in range(model._size):
# 			# 	for j in range(model._size):
# 			# 		if i!=j:
# 			# 			print(varZ[i,j].x, end=' ')
# 			# 	print('')
# 			# print("(", model._depot1[0], model._depot1[1],")")
# 			# for i in range(nb_tars):
# 			# 	print("(", varX[fseq[i+1]-1].x,varY[fseq[i+1]-1].x,")")
# 			# print('\n')
# 			# print(' ===>>> Gurobi Time: %g' % SP_m.Runtime)
# 			# print(' ===>>> SOCP OBJ: %g' % SP_m.objVal)

# 			return mu0, mu1, SP_m.objVal

# 	except gp.GurobiError as e:
# 		print('Error code ' + str(e.errno) + ': ' + str(e))
# 	except AttributeError:
# 		print('Encountered an attribute error')




def generate_GenOptimalityCut(model, fseq):
	nb_tars = model._size - 2
	vismat = np.zeros((model._size,model._size))
	for i in range(len(fseq)-1):
		vismat[fseq[i],fseq[i+1]] = 1
	# print(vismat)
	try:
		SP_m = gp.Model("SOCP")
		SP_m.setParam(GRB.Param.OutputFlag, 0)
		SP_m.setParam(GRB.Param.TimeLimit, 1000.0)

		varZ = SP_m.addVars(model._size, model._size, lb=0.001,ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Z")

		varX = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="X")
		varY = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Y")

		varS = SP_m.addVars(model._size, model._size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="S")
		varT = SP_m.addVars(model._size, model._size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="T")

		varK = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="K")
		varL = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="L")

		SP_m.update()

		# objective function
		obj = 0
		for i in range(model._size):
			for j in range(model._size):
				obj += vismat[i,j] * (varZ[i,j] - model._zbar[i,j])
		SP_m.setObjective(obj, GRB.MINIMIZE)

		# constraints - depot1
		for i in range(model._size):
			for j in range(i+1, model._size):
				SP_m.addConstr(varS[j,i]==varS[i,j])
				SP_m.addConstr(varT[j,i]==varT[i,j])
				SP_m.addConstr(varZ[j,i]==varZ[i,j])
			SP_m.update()

		for i in range(0, nb_tars):
			SP_m.addConstr(varS[0,i+1] == model._depot1[0]-varX[i], 'cS_'+str(0)+str(i+1))
			SP_m.addConstr(varS[i+1,0]==varS[0,i+1])
			SP_m.addConstr(varT[0,i+1] == model._depot1[1]-varY[i], 'cT_'+str(0)+str(i+1))
			SP_m.addConstr(varT[i+1,0]==varT[0,i+1])

		for i in range(0, nb_tars):
			SP_m.addConstr(varS[i+1,nb_tars+1] == varX[i] - model._depot2[0], 'cS_'+str(i+1)+str(nb_tars+1))
			SP_m.addConstr(varS[nb_tars+1,i+1]==varS[i+1,nb_tars+1])
			SP_m.addConstr(varT[i+1,nb_tars+1] == varY[i] - model._depot2[1], 'cT_'+str(i+1)+str(nb_tars+1))
			SP_m.addConstr(varT[nb_tars+1,i+1]==varT[i+1,nb_tars+1])

		for i in range(1, nb_tars+1):
			for j in range(i+1, nb_tars+1):
				SP_m.addConstr(varS[i,j] == varX[i-1]-varX[j-1], 'cS_'+str(i)+str(j))
				SP_m.addConstr(varS[j,i] == varS[i,j])
				SP_m.addConstr(varT[i,j] == varY[i-1]-varY[j-1], 'cT_'+str(i)+str(j))
				SP_m.addConstr(varT[j,i] == varT[i,j])
			SP_m.update()
		
		for i in range(model._size):
			for j in range(model._size):
				if i!=j:
					SP_m.addConstr(varS[i,j]*varS[i,j] + varT[i,j]*varT[i,j] <= varZ[i,j]*varZ[i,j], 'cZ_'+str(i)+str(j)) #constraint S
			SP_m.update()

		for i in range(0, nb_tars):
			SP_m.addConstr(varX[i]-model._Ox[i] == varK[i], 'cK_'+str(i)) 
			SP_m.addConstr(varY[i]-model._Oy[i] == varL[i], 'cL_'+str(i)) 
		SP_m.update()

		for i in range(0, nb_tars):
			SP_m.addConstr(varK[i]*varK[i] + varL[i]*varL[i] <= model._Or[i]**2, 'cR_'+str(i)) 
		SP_m.update()

		# SP_m.write('socp.lp')
		SP_m.optimize()

		''' ------------- model output  ----------------------'''
		if SP_m.status == GRB.OPTIMAL or SP_m.status == GRB.TIME_LIMIT:
		
			mu1 = []
			for i in range(model._size):
				for j in range(model._size):
					if i!=j and vismat[i,j] > 0:
						mu1.append([i,j,varZ[i,j].x])					
			mu0 = []
			for i in range(len(fseq)-2):
				for j in range(i+2, len(fseq)-1):
					if i == 0:
						val =  math.sqrt((model._depot1[0] - varX[fseq[j]-1].x)**2 + (model._depot1[1]-varY[fseq[j]-1].x)**2)
						val = val - model._zbar[0,fseq[j]]
						mu0.append([0, fseq[j],val])
					if i != 0:
						val =  math.sqrt((varX[fseq[i]-1].x - varX[fseq[j]-1].x)**2 + (varY[fseq[i]-1].x-varY[fseq[j]-1].x)**2)	
						val = val - model._zbar[fseq[i],fseq[j]]
						mu0.append([fseq[i], fseq[j], val])	


			for i in range(1, len(fseq)-2):
				val =  math.sqrt((model._depot2[0] - varX[fseq[i]-1].x)**2 + (model._depot2[1]-varY[fseq[i]-1].x)**2)	
				val = val - model._zbar[fseq[i],model._size-1]
				mu0.append([fseq[i], model._size-1,val])
				# mu0.append([ model._size-1,fseq[i],val])

			# print(mu1)
			# print(mu0)
			# for i in range(model._size):
			# 	for j in range(model._size):
			# 		if i!=j:
			# 			print(varZ[i,j].x, end=' ')
			# 	print('')
			# print("(", model._depot1[0], model._depot1[1],")")
			# for i in range(nb_tars):
			# 	print("(", varX[fseq[i+1]-1].x,varY[fseq[i+1]-1].x,")")
			# print('\n')
			# print(' ===>>> Gurobi Time: %g' % SP_m.Runtime)
			# print(' ===>>> SOCP OBJ: %g' % SP_m.objVal)

			return mu0, mu1, SP_m.objVal

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))
	except AttributeError:
		print('Encountered an attribute error')


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
			if curSol[cur_node, i] > 0.0001:
				zeronext = False
			if curSol[cur_node, i] > 0.0001 and visited[i] == False:
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

if __name__ == "__main__":
	instance = argv[1]
	depot, Ox, Oy, Or = read_instance_Behdani(instance)
	# plot_Behdani_instance(depot, Ox, Oy, Or)

	solve_MasterProb(depot, Ox, Oy, Or)























# def generate_GenOptimalityCut(model, fseq):
# 	nb_tars = model._size - 1
# 	# re-order the circles by following the given sequence 
# 	OX = []
# 	OY = []
# 	OR = []
# 	# print(fseq)
# 	# print(model._Ox)
# 	# print(fseq[1:])
# 	for i in fseq[1:]:
# 		OX.append(model._Ox[i-1])
# 		OY.append(model._Oy[i-1])
# 		OR.append(model._Or[i-1])

# 	# print(model._Ox, model._Oy)
# 	# print(fseq)
# 	# print(OX, OY,OR)
# 	try:
# 		SP_m = gp.Model("SOCP")
# 		SP_m.setParam(GRB.Param.OutputFlag, 0)
# 		SP_m.setParam(GRB.Param.TimeLimit, 1000.0)
# 		LBZ = []
# 		LBZ.append(model._zbar[fseq[nb_tars],0])
# 		for i in range(nb_tars):
# 			LBZ.append(model._zbar[fseq[i],fseq[i+1]])
# 		# print(LBZ)
# 		# distance between two turning points
# 		Z = SP_m.addVars(nb_tars+1, lb=LBZ, vtype=GRB.CONTINUOUS,name="spZ")

# 		X = SP_m.addVars(nb_tars,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="spX")
# 		Y = SP_m.addVars(nb_tars,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="spY")

# 		W = SP_m.addVars(nb_tars+1,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="spW")
# 		U = SP_m.addVars(nb_tars+1,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="spU")

# 		S = SP_m.addVars(nb_tars,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="spS")
# 		T = SP_m.addVars(nb_tars,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="spT")


# 		obj = 0
# 		obj += 1.0 * (Z[0]-model._zbar[fseq[nb_tars]][0])
# 		for i in range(0,nb_tars):
# 			obj += 1.0 * (Z[i]-model._zbar[fseq[i]][fseq[i+1]])
# 		SP_m.setObjective(obj, GRB.MINIMIZE)
# 		# SP_m.modelSense = GRB.MINIMIZE
	

# 		SP_m.addConstr(X[nb_tars-1] - model._depot[0] == W[0], 'cw'+str(0)) #constraint w
# 		SP_m.addConstr(model._depot[0] - X[0] == W[1], 'cw'+str(1)) 		
# 		for i in range(1, nb_tars):		
# 			SP_m.addConstr(X[i-1]-X[i] == W[i+1], 'cw'+str(i+1)) 
# 		SP_m.update()

# 		SP_m.addConstr(Y[nb_tars-1] - model._depot[1] == U[0], 'cu'+str(0)) 
# 		SP_m.addConstr(model._depot[1] - Y[0] == U[1], 'cu'+str(1)) 		
# 		for i in range(1, nb_tars):		
# 			SP_m.addConstr(Y[i-1]-Y[i] == U[i+1], 'cu'+str(i+1)) 
# 		SP_m.update()

# 		for i in range(0, nb_tars):		
# 			SP_m.addConstr(OX[i]-X[i] == S[i], 'Dx_'+str(i)) #Delta x: distance to the center
# 		SP_m.update()

# 		for i in range(0, nb_tars):		
# 			SP_m.addConstr(OY[i]-Y[i] == T[i], 'Dy_'+str(i)) 
# 		SP_m.update()


# 		for i in range(0, nb_tars + 1):		
# 			SP_m.addConstr(W[i]*W[i] + U[i]*U[i] <= Z[i]* Z[i], 'Dist_'+str(i)) 
# 		SP_m.update()		

# 		for i in range(0, nb_tars):		
# 			SP_m.addConstr(S[i]*S[i] + T[i]*T[i] <= OR[i]**2, 'InCir_'+str(i)) 
# 		SP_m.update()		
	

# 		# model.write('socp.lp')
# 		SP_m.optimize()

# 		''' ------------- model output  ----------------------'''
# 		if SP_m.status == GRB.OPTIMAL or SP_m.status == GRB.TIME_LIMIT:
# 			mu = []
# 			mu.append([fseq[nb_tars],0,-(Z[0].x -model._zbar[fseq[nb_tars],0])])
# 			for i in range(nb_tars):
# 				mu.append([fseq[i],fseq[i+1],-(Z[i+1].x-model._zbar[fseq[i],fseq[i+1]])])




# 			# for i in range(nb_tars):
# 			# 	print("(", X[i].x,Y[i].x,")")
# 			# print('\n')
# 			# print(' ===>>> Gurobi Time: %g' % SP_m.Runtime)
# 			# print(' ===>>> SOCP OBJ: %g' % SP_m.objVal)

# 			return mu, SP_m.objVal



# 	except gp.GurobiError as e:
# 		print('Error code ' + str(e.errno) + ': ' + str(e))

# 	except AttributeError:
# 		print('Encountered an attribute error')
