import gurobipy as gp
from gurobipy import GRB
from itertools import combinations
from sys import argv
import re
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pch

def generate_fractionalSolCut(model, curSol):
	nb_tars = model._size - 2
	# print(curSol)
	try:
		SP_m = gp.Model("SOCP")
		SP_m.setParam(GRB.Param.OutputFlag, 0)
		SP_m.setParam(GRB.Param.TimeLimit, 1000.0)

		varZ = SP_m.addVars(model._size, model._size, lb=0.00001,ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Z")

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
				obj += curSol[i,j] * (varZ[i,j] - model._zbar[i,j])
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
					if i!=j and curSol[i,j] > 0:
						# mu1.append([i,j,varZ[i,j].x])
						mu1.append([i,j,varZ[i,j].x - model._zbar[i,j]])					

			mu0 = []
			return mu0, mu1, SP_m.objVal

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))
	except AttributeError:
		print('Encountered an attribute error')


def generate_GenOptimalityCut(model, fseq):
	# print(fseq)

	nb_tars = model._size - 2
	vismat = np.zeros((model._size,model._size))
	for i in range(len(fseq)-1):
		vismat[fseq[i],fseq[i+1]] = 1
	# print(vismat)
	try:
		SP_m = gp.Model("SOCP")
		SP_m.setParam(GRB.Param.OutputFlag, 0)
		SP_m.setParam(GRB.Param.TimeLimit, 1000.0)

		varZ = SP_m.addVars(model._size, model._size, lb=0.000000,ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Z")

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

		# separator 
		# itr = 1
		# for sep in model._separators:
		# 	for i in range(0, nb_tars):
		# 		SP_m.addConstr( sep[0] * varX[i] + sep[1] * varY[i] + sep[2] <= 0, 'sep_'+str(itr)+str(i)) 
		# 	itr += 1
		# SP_m.update()
		# SP_m.write('socp2.lp')
		SP_m.optimize()
		# print('here1')

		''' ------------- model output  ----------------------'''
		if SP_m.status == GRB.OPTIMAL or SP_m.status == GRB.TIME_LIMIT:
			# print('here2')
			mu1 = []
			for i in range(model._size):
				for j in range(model._size):
					if i != j and vismat[i,j] > 0:
						# mu1.append([i,j,varZ[i,j].x])
						mu1.append([i,j,varZ[i,j].x - model._zbar[i,j]])
						# mu1.append([j,i,varZ[i,j].x - model._zbar[i,j]])					


			mu0 = []
			for i in range(nb_tars):
				for j in range(i+2, nb_tars+1):
					if i == 0:
						val =  math.sqrt((model._depot1[0] - varX[fseq[j]-1].x)**2 + (model._depot1[1]-varY[fseq[j]-1].x)**2)
						val = val - model._zbar[0,fseq[j]]
						mu0.append([0, fseq[j],val])
						# mu0.append([fseq[j], 0,val])

					if i != 0:
						val =  math.sqrt((varX[fseq[i]-1].x - varX[fseq[j]-1].x)**2 + (varY[fseq[i]-1].x-varY[fseq[j]-1].x)**2)	
						val = val - model._zbar[fseq[i],fseq[j]]
						mu0.append([fseq[i], fseq[j], val])	
						# mu0.append([fseq[j], fseq[i], val])

			for i in range(1, nb_tars):
				val =  math.sqrt((model._depot2[0] - varX[fseq[i]-1].x)**2 + (model._depot2[1]-varY[fseq[i]-1].x)**2)	
				val = val - model._zbar[fseq[i],nb_tars+1]
				mu0.append([fseq[i], nb_tars+1, val])
				# mu0.append([nb_tars+1, fseq[i], val])

			# print(mu1)
			# print(fseq)
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



def generate_GBC_CvxPolyNgbs(model, fseq):
	nb_tars = model._size - 2
	vismat = np.zeros((model._size,model._size))
	for i in range(len(fseq)-1):
		vismat[fseq[i],fseq[i+1]] = 1
	# print(vismat)
	try:
		SP_m = gp.Model("SOCP_CvxPolyNgbs")
		SP_m.setParam(GRB.Param.OutputFlag, 0)
		SP_m.setParam(GRB.Param.TimeLimit, 1000.0)
		varZ = SP_m.addVars(model._size, model._size, lb=0.00001,ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Z")
		varX = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="X")
		varY = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Y")
		varS = SP_m.addVars(model._size, model._size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="S")
		varT = SP_m.addVars(model._size, model._size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="T")
		SP_m.update()

		# objective function
		obj = 0
		for i in range(model._size):
			for j in range(model._size):
				obj += vismat[i,j] * (varZ[i,j] - model._zbar[i,j])
				# obj += vismat[i,j] * varZ[i,j]

		SP_m.setObjective(obj, GRB.MINIMIZE)
		SP_m.update()

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

		# add convex polygon constraints 
		idx = 0
		for lpc in model._LPC:
			for c in lpc:
				SP_m.addConstr(c[0] * varX[idx] + c[1] * varY[idx] + c[2] <= 0)
			idx += 1
		SP_m.update()

		# add boundary-projection-closed separators if exists
		if model._BPC_sep_flag:
			for sep in model._BPC_separators:
				for k in range(nb_tars):
					SP_m.addConstr(sep[0] * varX[k] + sep[1] * varY[k] + sep[2] <= 0)
		SP_m.update()

		SP_m.optimize()

		''' ------------- model output  ----------------------'''
		if SP_m.status == GRB.OPTIMAL or SP_m.status == GRB.TIME_LIMIT:
		
			mu1 = []
			for i in range(model._size):
				for j in range(model._size):
					if i!=j and vismat[i,j] > 0:
						# mu1.append([i,j,varZ[i,j].x])
						mu1.append([i,j,varZ[i,j].x - model._zbar[i,j]])
						# mu1.append([j,i,varZ[i,j].x - model._zbar[i,j]])					


			mu0 = []
			# for i in range(nb_tars):
			# 	for j in range(i+2, nb_tars+1):
			# 		if i == 0:
			# 			val =  math.sqrt((model._depot1[0] - varX[fseq[j]-1].x)**2 + (model._depot1[1]-varY[fseq[j]-1].x)**2)
			# 			val = val - model._zbar[0,fseq[j]]
			# 			mu0.append([0, fseq[j],val])
			# 			# mu0.append([fseq[j], 0,val])

			# 		if i != 0:
			# 			val =  math.sqrt((varX[fseq[i]-1].x - varX[fseq[j]-1].x)**2 + (varY[fseq[i]-1].x-varY[fseq[j]-1].x)**2)	
			# 			val = val - model._zbar[fseq[i],fseq[j]]
			# 			mu0.append([fseq[i], fseq[j], val])	
			# 			# mu0.append([fseq[j], fseq[i], val])

			# for i in range(1, nb_tars):
			# 	val =  math.sqrt((model._depot2[0] - varX[fseq[i]-1].x)**2 + (model._depot2[1]-varY[fseq[i]-1].x)**2)	
			# 	val = val - model._zbar[fseq[i],nb_tars+1]
			# 	mu0.append([fseq[i], nb_tars+1, val])
				# mu0.append([nb_tars+1, fseq[i], val])

			# print(mu1)
			# print(fseq)
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

''' generate Generalized Benders cuts when neigbhorhood regions are concave polygon'''
def generate_GBC_MixedPolyNgbs(model, fseq):
	nb_tars = model._size - 2
	vismat = np.zeros((model._size,model._size))
	for i in range(len(fseq)-1):
		vismat[fseq[i],fseq[i+1]] = 1
	# print(vismat)
	max_num_ploygonedges = 0
	for polygon in model._POLYGONs:
		vertices = list(polygon.exterior.coords)
		if len(vertices)-1 > max_num_ploygonedges:
			max_num_ploygonedges = len(vertices)-1

	try:
		SP_m = gp.Model("SubProblem_MixedPolyNgbs")
		SP_m.setParam(GRB.Param.OutputFlag, 0)
		SP_m.setParam(GRB.Param.TimeLimit, 1000.0)
		varZ = SP_m.addVars(model._size, model._size, lb=0.00001,ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Z")
		varX = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="X")
		varY = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Y")
		varS = SP_m.addVars(model._size, model._size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="S")
		varT = SP_m.addVars(model._size, model._size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="T")
		varLambda = SP_m.addVars(model._size-2, max_num_ploygonedges, lb=0, ub=1.0, vtype=GRB.CONTINUOUS,name="Lambda")
		varPhi 	  = SP_m.addVars(model._size-2, max_num_ploygonedges, vtype=GRB.BINARY, name="Phi")

		SP_m.update()

		# objective function
		obj = 0
		for i in range(model._size):
			for j in range(model._size):
				obj += vismat[i,j] * (varZ[i,j] - model._zbar[i,j])
				# obj += vismat[i,j] * varZ[i,j]

		SP_m.setObjective(obj, GRB.MINIMIZE)
		SP_m.update()

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

		for pindex, polygon in enumerate(model._POLYGONs):
			if model._neighborhood_types[pindex] == 'concave':
				# add concave polygon constraints 
				constrX, constrY, constrSumPhi = 0, 0, 0
				vertices = list(polygon.exterior.coords)
				for index in range(len(vertices)-1):
					constrX += varLambda[pindex,index]*(vertices[index][0] - vertices[index+1][0]) + varPhi[pindex,index]*vertices[index+1][0]
					constrY += varLambda[pindex,index]*(vertices[index][1] - vertices[index+1][1]) + varPhi[pindex,index]*vertices[index+1][1]
					constrSumPhi += varPhi[pindex,index]
				SP_m.addConstr(constrX==varX[pindex])
				SP_m.addConstr(constrY==varY[pindex])
				SP_m.addConstr(constrSumPhi==1, "sum_Phi_constr")
				SP_m.update()
			else:
				# add convex polygon constraints
				lpc = model._LPC[pindex]
				for c in lpc:
					SP_m.addConstr(c[0] * varX[pindex] + c[1] * varY[pindex] + c[2] <= 0)
				SP_m.update()


		for pindex, polygon in enumerate(model._POLYGONs):
			if model._neighborhood_types[pindex] == 'concave':
				vertices = list(polygon.exterior.coords)
				for index in range(len(vertices)-1):
					SP_m.addConstr(varLambda[pindex,index]<=varPhi[pindex,index])				
				SP_m.update()

		# add boundary-projection-closed separators if exists
		if True:
			for sep in model._BPC_separators:
				for k in range(nb_tars):
					SP_m.addConstr(sep[0] * varX[k] + sep[1] * varY[k] + sep[2] <= 0)
		SP_m.update()

		SP_m.optimize()

		''' ------------- model output  ----------------------'''
		if SP_m.status == GRB.OPTIMAL or SP_m.status == GRB.TIME_LIMIT:
		
			mu1 = []
			for i in range(model._size):
				for j in range(model._size):
					if i!=j and vismat[i,j] > 0:
						# mu1.append([i,j,varZ[i,j].x])
						mu1.append([i,j,varZ[i,j].x - model._zbar[i,j]])
						# mu1.append([j,i,varZ[i,j].x - model._zbar[i,j]])					


			mu0 = []
			return mu0, mu1, SP_m.objVal

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))
	except AttributeError:
		print('Subproblem encountered an attribute error')




def solve_SOCP_Disc(_depot, _Ox, _Oy, _Or, fseq):
	nb_tars = len(fseq) -2
	# re-order the circles by following the given sequence 
	OX = [_depot[0]]
	OY = [_depot[1]]
	OR = [0]
	for i in range(nb_tars):
		OX.append(_Ox[fseq[i+1]-1])
		OY.append(_Oy[fseq[i+1]-1])
		OR.append(_Or[fseq[i+1]-1])
	OX.append(_depot[0])
	OY.append(_depot[1])
	OR.append(0)
	# print(model._Ox, model._Oy)
	# print(fseq)
	# print(OX, OY,OR)
	try:
		SP_m = gp.Model("SOCP_Disc")
		SP_m.setParam(GRB.Param.OutputFlag, 0)
		SP_m.setParam(GRB.Param.TimeLimit, 1000.0)
		
		# distance between two turning points
		Z = SP_m.addVars(nb_tars+1, lb=0, vtype=GRB.CONTINUOUS,name="spZ")

		X = SP_m.addVars(nb_tars+2,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="spX")
		Y = SP_m.addVars(nb_tars+2,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="spY")

		W = SP_m.addVars(nb_tars+1,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="spW")
		U = SP_m.addVars(nb_tars+1,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="spU")

		S = SP_m.addVars(nb_tars+2,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="spS")
		T = SP_m.addVars(nb_tars+2,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="spT")


		obj = 0
		for i in range(0,nb_tars+1):
			obj += 1.0 * Z[i]
		SP_m.setObjective(obj, GRB.MINIMIZE)
		# SP_m.modelSense = GRB.MINIMIZE
			
		X[0].lb = OX[0]
		X[0].ub = OX[0]
		Y[0].lb = OY[0]
		Y[0].ub = OY[0]
		X[nb_tars+1].lb = OX[nb_tars+1]
		X[nb_tars+1].ub = OX[nb_tars+1]
		Y[nb_tars+1].lb = OY[nb_tars+1]
		Y[nb_tars+1].ub = OY[nb_tars+1]
		S[0].lb = 0
		S[0].ub = 0
		T[0].lb = 0
		T[0].ub = 0
		S[nb_tars+1].lb = 0
		S[nb_tars+1].ub = 0
		T[nb_tars+1].lb = 0
		T[nb_tars+1].ub = 0

		for i in range(1, nb_tars+2):
			SP_m.addConstr(X[i]-X[i-1] == W[i-1], 'cW'+str(i-1)) 
		SP_m.update()
		for i in range(1, nb_tars+2):
			SP_m.addConstr(Y[i]-Y[i-1] == U[i-1], 'cU'+str(i-1)) 
		SP_m.update()

		for i in range(0, nb_tars + 1):		
			SP_m.addConstr(W[i]*W[i] + U[i]*U[i] <= Z[i]* Z[i], 'Dist_'+str(i)) 
		SP_m.update()		
		
		for i in range(1, nb_tars+1):		
			SP_m.addConstr(OX[i]-X[i] == S[i], 'Dx_'+str(i)) #Delta x: distance to the center
		SP_m.update()

		for i in range(1, nb_tars+1):		
			SP_m.addConstr(OY[i]-Y[i] == T[i], 'Dy_'+str(i)) #Delta x: distance to the center
		SP_m.update()

		for i in range(1, nb_tars+1):		
			SP_m.addConstr(S[i]*S[i] + T[i]*T[i] <= OR[i]**2, 'InCir_'+str(i)) 
		SP_m.update()		
	

		# model.write('socp.lp')
		SP_m.optimize()

		''' ------------- model output  ----------------------'''
		if SP_m.status == GRB.OPTIMAL or SP_m.status == GRB.TIME_LIMIT:
			# print(' ===>>> Gurobi Time: %g' % SP_m.Runtime)
			# print(' ===>>> SOCP OBJ: %g' % SP_m.objVal)
			optTour = []
			for i in range(0, nb_tars+2):
				optTour.append([X[i].x, Y[i].x])
			optTour = np.array(optTour)
			# print(optTour)
			return optTour, SP_m.objVal



	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	except AttributeError:
		print('Encountered an attribute error')


def solve_SOCP_CvxPolyNgbs(model, fseq):
	nb_tars = model._size - 2
	vismat = np.zeros((model._size,model._size))
	for i in range(len(fseq)-1):
		vismat[fseq[i],fseq[i+1]] = 1
	# print(vismat)
	try:
		SP_m = gp.Model("SOCP_CvxPolyNgbs")
		SP_m.setParam(GRB.Param.OutputFlag, 0)
		SP_m.setParam(GRB.Param.TimeLimit, 1000.0)

		varZ = SP_m.addVars(model._size, model._size, lb=0.00001,ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Z")
		varX = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="X")
		varY = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Y")
		varS = SP_m.addVars(model._size, model._size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="S")
		varT = SP_m.addVars(model._size, model._size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="T")
		SP_m.update()

		# objective function
		obj = 0
		for i in range(model._size):
			for j in range(model._size):
				# obj += vismat[i,j] * (varZ[i,j] - model._zbar[i,j])
				obj += vismat[i,j] * varZ[i,j]

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

		# add convex polygon constraints 
		idx = 0
		# print(model._LPC)
		for lpc in False:
			for c in lpc:
				SP_m.addConstr(c[0]*varX[idx]+c[1]*varY[idx] + c[2] <= 0)
			idx += 1
		SP_m.update()
		SP_m.optimize()

		''' ------------- model output  ----------------------'''
		if SP_m.status == GRB.OPTIMAL or SP_m.status == GRB.TIME_LIMIT:
			optTour = [model._depot1]
			for i in range(1, nb_tars+1):
				optTour.append([varX[fseq[i]-1].x, varY[fseq[i]-1].x])
			optTour.append(model._depot2)
			# print(optTour)
			optTour = np.array(optTour)
			return optTour, SP_m.objVal

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))
	except AttributeError:
		print('Encountered an attribute error')


def solve_ConcavePolyNgbs(model, fseq):
	nb_tars = model._size - 2
	vismat = np.zeros((model._size,model._size))
	for i in range(len(fseq)-1):
		vismat[fseq[i],fseq[i+1]] = 1
	# print(vismat)
	max_num_ploygonedges = 0
	for polygon in model._POLYGONs:
		vertices = list(polygon.exterior.coords)
		if len(vertices)-1 > max_num_ploygonedges:
			max_num_ploygonedges = len(vertices)-1
	try:
		SP_m = gp.Model("SOCP_CvxPolyNgbs")
		SP_m.setParam(GRB.Param.OutputFlag, 0)
		SP_m.setParam(GRB.Param.TimeLimit, 1000.0)

		varZ = SP_m.addVars(model._size, model._size, lb=0.00001,ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Z")
		varX = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="X")
		varY = SP_m.addVars(model._size-2, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="Y")
		varS = SP_m.addVars(model._size, model._size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="S")
		varT = SP_m.addVars(model._size, model._size, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,name="T")
		varLambda = SP_m.addVars(model._size-2, max_num_ploygonedges, lb=0, ub=1.0, vtype=GRB.CONTINUOUS,name="Lambda")
		varPhi 	  = SP_m.addVars(model._size-2, max_num_ploygonedges, vtype=GRB.BINARY, name="Phi")
		SP_m.update()

		# objective function
		obj = 0
		for i in range(model._size):
			for j in range(model._size):
				# obj += vismat[i,j] * (varZ[i,j] - model._zbar[i,j])
				obj += vismat[i,j] * varZ[i,j]

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

		# add concave polygon constraints 
		for pindex, polygon in enumerate(model._POLYGONs):
			constrX, constrY, constrSumPhi = 0, 0, 0
			vertices = list(polygon.exterior.coords)
			for index in range(len(vertices)-1):
				constrX += varLambda[pindex,index]*(vertices[index][0] - vertices[index+1][0]) + varPhi[pindex,index]*vertices[index+1][0]
				constrY += varLambda[pindex,index]*(vertices[index][1] - vertices[index+1][1]) + varPhi[pindex,index]*vertices[index+1][1]
				constrSumPhi += varPhi[pindex,index]
			SP_m.addConstr(constrX==varX[pindex], "x_in_concavepoly_constr")
			SP_m.addConstr(constrY==varY[pindex], "y_in_concavepoly_constr")
			SP_m.addConstr(constrSumPhi==1, "sum_Phi_constr")
			SP_m.update()

		for pindex, polygon in enumerate(model._POLYGONs):
			vertices = list(polygon.exterior.coords)
			for index in range(len(vertices)-1):
				SP_m.addConstr(varLambda[pindex,index]<=varPhi[pindex,index], "lambda<Phi_constr")				
			SP_m.update()
		SP_m.update()
		SP_m.optimize()

		''' ------------- model output  ----------------------'''
		if SP_m.status == GRB.OPTIMAL or SP_m.status == GRB.TIME_LIMIT:
			optTour = [model._depot1]
			for i in range(1, nb_tars+1):
				optTour.append([varX[fseq[i]-1].x, varY[fseq[i]-1].x])
			optTour.append(model._depot2)
			# print(optTour)
			optTour = np.array(optTour)
			return optTour, SP_m.objVal

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))
	except AttributeError:
		print('Solve subproblem directly: error')