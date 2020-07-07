import gurobipy as gp
from gurobipy import GRB
from itertools import combinations

from sys import argv
import re
import math
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, PathPatch
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
import mpl_toolkits.mplot3d.art3d as art3d


def read_toy(instancefile):
	file_ = open(instancefile, 'r')
	X = []
	Y = []
	R = []
	line_ = file_.readline()
	list_ = re.split(" |\t|\n", line_)
	depot = [float(list_[0]), float(list_[1])]
	# line_ = file_.readline()
	# list_ = re.split(" |\t|\n", line_)
	# depot2 = [float(list_[0]), float(list_[1])]
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
			R.append(float(list_[2]))
	file_.close()
	return depot, X, Y, R


def plot_toy(X, Y, R, path):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	# plt.axis('off')
	ax.set_xlim([min(X)-max(R), max(X)+max(R)])
	ax.set_ylim([min(Y)-max(R), max(Y)+max(R)])
	ax.set_aspect('equal', adjustable='box')
	# circle = plt.Circle((depot1[0], depot1[1]), radius=0.2, edgecolor='k', alpha = 0.5)
	# ax.add_artist(circle)	
	# circle = plt.Circle((depot2[0], depot2[1]), radius=0.2, edgecolor='k', alpha = 0.5)
	# ax.add_artist(circle)	

	
	for i in range(len(X)):
		circle = plt.Circle((X[i], Y[i]), radius=R[i], edgecolor='r',fill=False, alpha = 0.5)
		ax.add_artist(circle)
	plt.plot(path[0],path[1])

	plt.show()


''' ----------------------------------------- '''
''' solve SOCP model '''
''' ----------------------------------------- '''
def solve_model(depot, oX, oY, oR):
	nb_tars = len(oX)
	try:
		# Create a new model
		model = gp.Model("SOCP")
		model.setParam(GRB.Param.OutputFlag, 1)
		model.setParam(GRB.Param.TimeLimit, 1000.0)

		# distance between two turning points
		varZ = model.addVars(nb_tars+1, lb=[0]*(nb_tars+1), vtype=GRB.CONTINUOUS,name="Z")

		varX = model.addVars(nb_tars,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="X")
		varY = model.addVars(nb_tars,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="Y")

		varW = model.addVars(nb_tars+1,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="W")
		varU = model.addVars(nb_tars+1,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="U")

		varS = model.addVars(nb_tars,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="S")
		varT = model.addVars(nb_tars,lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="T")


		obj = 0
		for i in range(0,nb_tars+1):
			obj += 1.0 * varZ[i]
		model.setObjective(obj, GRB.MINIMIZE)
		model.modelSense = GRB.MINIMIZE
	

		model.addConstr(varX[nb_tars-1] - depot[0] == varW[0], 'cw'+str(0)) 
		model.addConstr(depot[0] - varX[0] == varW[1], 'cw'+str(1)) 		
		for i in range(1, nb_tars):		
			model.addConstr(varX[i-1]-varX[i] == varW[i+1], 'cw'+str(i+1)) 
		model.update()

		model.addConstr(varY[nb_tars-1] - depot[1] == varU[0], 'cu'+str(0)) 
		model.addConstr(depot[1] - varY[0] == varU[1], 'cu'+str(1)) 		
		for i in range(1, nb_tars):		
			model.addConstr(varY[i-1]-varY[i] == varU[i+1], 'cu'+str(i+1)) 
		model.update()

		for i in range(0, nb_tars):		
			model.addConstr(oX[i]-varX[i] == varS[i], 'Dx_'+str(i)) 
		model.update()

		for i in range(0, nb_tars):		
			model.addConstr(oY[i]-varY[i] == varT[i], 'Dy_'+str(i)) 
		model.update()


		for i in range(0, nb_tars + 1):		
			model.addConstr(varW[i]*varW[i] + varU[i]*varU[i] <= varZ[i]* varZ[i], 'Dist_'+str(i)) 
		model.update()		

		for i in range(0, nb_tars):		
			model.addConstr(varS[i]*varS[i] + varT[i]*varT[i] <= oR[i]**2, 'InCir_'+str(i)) 
		model.update()		

		# separators 
		# for i in range(0, nb_tars-1):
		# 	model.addConstr(varX[i] <= varY[i], 'sepa_'+str(i)) 
		# model.update()		

		# model.write('socp.lp')
		model.optimize()

		''' ------------- model output  ----------------------'''
		if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
			pathX = []
			pathY = []
			pathX.append(depot[0])
			pathY.append(depot[1])
			for i in range(nb_tars):
				pathX.append(varX[i].x)
				pathY.append(varY[i].x)
				# print('(', varX[i].x, ',', varY[i].x, ')', end=' ')
			# print('\n')
			pathX.append(depot[0])
			pathY.append(depot[1])
			# for i in range(nb_tars+1):
			# 	print(varZ[i].x, end=' ')
			print('\n')
			print('Gurobi Time: %g' % model.Runtime)
			print('OBJ: %g' % obj.getValue())
			path = [pathX, pathY]
			return path


	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	except AttributeError:
		print('Encountered an attribute error')

''' ----------------------------------------- '''
''' solve SOCP model '''
''' ----------------------------------------- '''
def solve_model2(depot, oX, oY, oR):
	nb_tars = len(oX)
	try:
		# Create a new model
		model = gp.Model("SOCP")
		model.setParam(GRB.Param.OutputFlag, 1)
		model.setParam(GRB.Param.TimeLimit, 1000.0)

		# distance between two turning points
		z0 = model.addVar(lb=0.0, vtype=GRB.CONTINUOUS,name="z0")
		z1 = model.addVar(lb=0.0, vtype=GRB.CONTINUOUS,name="z1")
		z2 = model.addVar(lb=0.0, vtype=GRB.CONTINUOUS,name="z2")

		x1 = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="x1")
		y1 = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="y1")
		x2 = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="x2")
		y2 = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="y2")

		dx0 = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="dx0")
		dx1 = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="dx1")
		dx2 = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="dx2")

		dy0 = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="dy0")
		dy1 = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="dy1")
		dy2 = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="dy2")
		
		dr1x = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="dr1x")
		dr1y = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="dr1y")

		dr2x = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="dr2x")
		dr2y = model.addVar(lb=-gp.GRB.INFINITY,vtype=GRB.CONTINUOUS,name="dr2y")

		obj = z0 + z1 + z2
		model.setObjective(obj, GRB.MINIMIZE)
	

		model.addConstr(dx0 == x1 - depot[0], '1_dx') 
		model.addConstr(dy0 == y1 - depot[1], '1_dy') 

		model.addConstr(dx1 == x2 - x1, '2_dx') 
		model.addConstr(dy1 == y2 - y1, '2_dy') 

		model.addConstr(dx2 == depot[0] - x2, '3_dx') 
		model.addConstr(dy2 == depot[1] - y2, '3_dy') 

		model.addConstr(dx0*dx0+ dy0*dy0 <= z0*z0, 'dz0')
		model.addConstr(dx1*dx1+ dy1*dy1 <= z1*z1, 'dz1')
		model.addConstr(dx2*dx2+ dy2*dy2 <= z2*z2, 'dz2')

        #---
		model.addConstr(dr1x == x1 - oX[0], 'rx1') 
		model.addConstr(dr1y == y1 - oY[0], 'ry1') 

		model.addConstr(dr2x == x2 - oX[1], 'rx2') 
		model.addConstr(dr2y == y2 - oY[1], 'ry2') 
		model.addConstr(dr1x*dr1x + dr1y * dr1y <= 1, 'dr1')
		model.addConstr(dr2x*dr2x + dr2y * dr2y <= 1, 'dr2')


		model.write('socp.lp')
		model.optimize()

		''' ------------- model output  ----------------------'''
		if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
			print(x1.x,y1.x)
			print(x2.x,y2.x)
			
			print('Gurobi Time: %g' % model.Runtime)
			print('OBJ: %g' % obj.getValue())


	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	except AttributeError:
		print('Encountered an attribute error')


if __name__ == "__main__":
	instance = argv[1]
	depot, oX, oY, oR = read_toy(instance)
	path = solve_model(depot, oX, oY, oR)
	plot_toy(oX, oY, oR, path)
