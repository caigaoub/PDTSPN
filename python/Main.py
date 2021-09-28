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


import InstanceReader
import CutPool


class MINLP:

	def __init__(self, inst, SD, mode='1-M-1'):
		self._instance = inst
		self._SD = SD
		self._mode = mode
		self._depot = inst._depot
		self._HULLs = inst._HULLs
		self._SEPs = inst._SEPs
		self._nb_tars = inst._nb_cvxps

	def _solve(self):
		try:
			'''Create a new model'''
			model = gp.Model("GBD_PDTSPN")

			'''Set up GRBmodel parameters '''
			self._set_parameter(model)

			'''Initialize GRBmodel input '''
			model._nb_SECs = 0
			model._nb_GBCs = 0
			model._size = self._nb_tars + 2  
			model._depot1, model._depot2 = self._depot, self._depot
			model._zbar = self._calc_zbar()
			model._LPC = self._convert_HULL_LPConsts()

			'''Add binary variables to the model'''
			self._add_vars(model)

			'''Create objective function '''
			self._set_objective(model)
			
			'''Add constraints to the model'''
			self._set_constraints(model)			

			'''Optimize the model with Callback function ''' 
			model.optimize(MINLP.callback_GBD)

			'''Terminate the model and output the results '''
			self._terminate(model)
			
		except gp.GurobiError as e:
			print('Error code ' + str(e.errno) + ': ' + str(e))
		except AttributeError:
			print('Encountered an attribute error')

	def _set_parameter(self, model):
		model.setParam(GRB.Param.OutputFlag, 1)
		# model.setParam(GRB.Param.TimeLimit, 10.0)
		model.Params.lazyConstraints = 1

	def _add_vars(self, model):
		model._varE = model.addVars(model._size, model._size, vtype=GRB.BINARY, name="E")
		model._theta = model.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="theta")
		model.update()

	def _set_objective(self, model):
		obj = 0
		for i in range(model._size):
			for j in range(model._size):
				obj += model._zbar[i,j] * model._varE[i,j] 
		obj += 1.0 * model._theta
		model.setObjective(obj, GRB.MINIMIZE)
		model.update()

	def _set_constraints(self, model):
		for i in range(0, model._size):
			constr1 = 0
			constr2 = 0
			for j in range(0, model._size):
				constr1 += model._varE[i,j]
				constr2 += model._varE[j,i]
			model.addConstr(constr1 == 1, 'r_'+str(i)) 
			model.addConstr(constr2 == 1, 'c_'+str(i)) 
		model.update()

		model._varE[self._nb_tars+1,0].lb = 1.0
		model._varE[0,self._nb_tars+1].ub = 0.0

		for i in range(0, model._size):
			model._varE[i,i].ub = 0.0
		model.update()
	
	def _terminate(self, model):
		if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
			curObj = np.zeros((model._size, model._size))
			for i in range(0, model._size):
				for j in range(0, model._size):
					curObj[i,j] = model._varE[i,j].x
					# print(varE[i,j].x, end=' ')
				# print('\n',end='')

			cirset, curBestSeq = MINLP.find_all_subtours(curObj, model._size)
			# print(curBestSeq)

			optTour, objLen = CutPool.solve_SOCP_CvxPolyNgbs(model, curBestSeq)
			self._plot_CvxPoly_Instance_wOptTour(optTour)
			# print(model._nb_SECs, '& ', 
			# 	  model._nb_GBCs, '& ', 
			# 	  '{:.2f}'.format(model.Runtime), '& ', 
			# 	  '{:.2f}'.format(model.MIPGap), '& ', 
			# 	  '{:.4f}'.format(objLen))

	@staticmethod
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

	def _calc_zbar(self):
		# self._nb_tars = len(self._HULLs)
		if self._nb_tars <= 0:
			print("the number of polygons is wrong!! ")
			exit()
			
		zbar = np.zeros((self._nb_tars+2,self._nb_tars+2))
		for i in range(0, self._nb_tars):
			hull = self._HULLs[i]
			a0 = np.array([self._depot[0], self._depot[1], 0])
			a1 = np.array([self._depot[0]+0.0000001, self._depot[1]+0.00000001, 0])
			dist = []
			for simplex in hull.simplices:
				b0 = np.array([hull.points[simplex[0]][0], hull.points[simplex[0]][1], 0])
				b1 = np.array([hull.points[simplex[1]][0], hull.points[simplex[1]][1], 0])
				pa, pb, md = self._closestDistanceBetweenLines(a0, a1, b0, b1, clampAll=True)
				dist.append(md)
			# print(dist)
			zbar[0, i+1] = min(dist)
			zbar[i+1, 0] = zbar[0, i+1]
			zbar[i+1, self._nb_tars+1] = zbar[0, i+1]
			zbar[self._nb_tars+1, i+1] = zbar[0, i+1]

		for i in range(0, self._nb_tars):
			for j in range(0, self._nb_tars):
				if i != j:
					hulli = self._HULLs[i]
					hullj = self._HULLs[j]
					dist = []
					for ksplx in hulli.simplices:
						a0 = np.array([hulli.points[ksplx[0]][0], hulli.points[ksplx[0]][1], 0])
						a1 = np.array([hulli.points[ksplx[1]][0], hulli.points[ksplx[1]][1], 0])
						for lsplx in hullj.simplices:
							b0 = np.array([hullj.points[lsplx[0]][0], hullj.points[lsplx[0]][1], 0])
							b1 = np.array([hullj.points[lsplx[1]][0], hullj.points[lsplx[1]][1], 0])
							pa, pb, md = self._closestDistanceBetweenLines(a0, a1, b0, b1, clampAll=True)
							dist.append(md)
					zbar[i+1, j+1] = min(dist)
					zbar[j+1, i+1] = zbar[i+1, j+1]
		return zbar

	def _closestDistanceBetweenLines(self, a0,a1,b0,b1,clampAll=False,clampA0=False,clampA1=False,clampB0=False,clampB1=False):
		''' Given two lines defined by numpy.array pairs (a0,a1,b0,b1)
			Return the closest points on each segment and their distance
		'''
		# If clampAll=True, set all clamps to True
		if clampAll:
			clampA0=True
			clampA1=True
			clampB0=True
			clampB1=True
		# Calculate denomitator
		A = a1 - a0
		B = b1 - b0
		magA = np.linalg.norm(A)
		magB = np.linalg.norm(B)
		_A = A / magA
		_B = B / magB
		cross = np.cross(_A, _B);
		denom = np.linalg.norm(cross)**2
		# If lines are parallel (denom=0) test if lines overlap.
		# If they don't overlap then there is a closest point solution.
		# If they do overlap, there are infinite closest positions, but there is a closest distance
		if not denom:
			d0 = np.dot(_A,(b0-a0))

			# Overlap only possible with clamping
			if clampA0 or clampA1 or clampB0 or clampB1:
				d1 = np.dot(_A,(b1-a0))

				# Is segment B before A?
				if d0 <= 0 >= d1:
					if clampA0 and clampB1:
						if np.absolute(d0) < np.absolute(d1):
							return a0,b0,np.linalg.norm(a0-b0)
						return a0,b1,np.linalg.norm(a0-b1)
				# Is segment B after A?
				elif d0 >= magA <= d1:
					if clampA1 and clampB0:
						if np.absolute(d0) < np.absolute(d1):
							return a1,b0,np.linalg.norm(a1-b0)
						return a1,b1,np.linalg.norm(a1-b1)
			# Segments overlap, return distance between parallel segments
			return None,None,np.linalg.norm(((d0*_A)+a0)-b0)
		# Lines criss-cross: Calculate the projected closest points
		t = (b0 - a0);
		detA = np.linalg.det([t, _B, cross])
		detB = np.linalg.det([t, _A, cross])
		t0 = detA/denom;
		t1 = detB/denom;
		pA = a0 + (_A * t0) # Projected closest point on segment A
		pB = b0 + (_B * t1) # Projected closest point on segment B
		# Clamp projections
		if clampA0 or clampA1 or clampB0 or clampB1:
			if clampA0 and t0 < 0:
				pA = a0
			elif clampA1 and t0 > magA:
				pA = a1
			if clampB0 and t1 < 0:
				pB = b0
			elif clampB1 and t1 > magB:
				pB = b1
			# Clamp projection A
			if (clampA0 and t0 < 0) or (clampA1 and t0 > magA):
				dot = np.dot(_B,(pA-b0))
				if clampB0 and dot < 0:
					dot = 0
				elif clampB1 and dot > magB:
					dot = magB
				pB = b0 + (_B * dot)
			# Clamp projection B
			if (clampB0 and t1 < 0) or (clampB1 and t1 > magB):
				dot = np.dot(_A,(pB-a0))
				if clampA0 and dot < 0:
					dot = 0
				elif clampA1 and dot > magA:
					dot = magA
				pA = a0 + (_A * dot)
		return pA,pB,np.linalg.norm(pA-pB)

	def _plot_CvxPoly_Instance_wOptTour(self, optTour):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_aspect('equal', adjustable='box')
		itr = 1
		cx = 0
		cy = 0
		# ax.set_xlim([0, 17])
		# ax.set_ylim([0, 14])
		plt.xlabel('x', fontsize=18)
		plt.ylabel('y', fontsize=18)
		ax.tick_params(axis='both', labelsize=16)
		ax.set_aspect('equal', adjustable='box')

		plt.plot(self._depot[0], self._depot[1], 'rs', ms=10)
		for hull in self._HULLs:
			for simplex in hull.simplices:
				plt.plot(hull.points[simplex, 0], hull.points[simplex, 1], 'k-')
				# centroid
				cx = np.mean(hull.points[hull.vertices,0])
				cy = np.mean(hull.points[hull.vertices,1])
					#Plot centroid
				plt.plot(cx, cy,'ko',ms=2)
			ax.annotate(str(itr), xy=(cx, cy),textcoords="offset points", xytext=(cx, cy),size=14)
			itr += 1
		plt.plot(optTour[:,0],optTour[:,1], linestyle='-', color = 'blue', markersize=3, lw=2)
		plt.show()
	
	@staticmethod		
	def callback_GBD(model, where):
		if where == GRB.Callback.MIPSOL:
			vals = model.cbGetSolution(model._varE)
			if False:
				for i in range(model._size):
					for j in range(model._size):
						print(vals[i,j],end=' ')
					print(end='\n')

			# find_Subtour(vals, model._size)
			Cirset,fseq = MINLP.find_all_subtours(vals, model._size)
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
			else:
				# '''supply_demand balance check '''
				# SD_balanced = True
				# for i in range(1, model._size):
				# 	sum_sd = 0
				# 	for j in range(0, i):
				# 		sum_sd += model._SD[fseq[j]]
				# 	if sum_sd < 0:
				# 		# print(sum_sd)
				# 		SD_balanced = False
				# 		break
				# if not SD_balanced:
				# # if False:
				# 	constr = 0
				# 	# print(fseq)
				# 	for i in range(0, model._size-1):
				# 		constr += model._varE[fseq[i],fseq[i+1]]
				# 		# print(fseq[i],fseq[i+1], end=', ')
				# 	# print('')
				# 	model.cbLazy(constr <= model._size-2)
				

				# 	# print('....addded ')
				# else:
				# optTour, objLen = CutPool.solve_SOCP_Disc(self._depot, Ox, Oy, Or, fseq)
				# # print(fseq, objLen)
				# mu0, mu1, obj = CutPool.generate_GenOptimalityCut(model, fseq)
				# print(mu0,mu1)
				mu0, mu1, obj = CutPool.generate_GBC_CvxPolyNgbs(model, fseq)
				constr = 0
				for el in mu1:
					constr += el[2] * (model._varE[el[0],el[1]] - 1.0)
				# for el in mu0:
				# 	constr += el[2] * (model._varE[el[0],el[1]] - 0.0)
				model.cbLazy(obj + constr <= model._theta)
				model._nb_GBCs += 1
	
	def _convert_HULL_LPConsts(self):
		LPC = []
		for hull in self._HULLs:
			lpc = []
			# centroid of current convex hull
			cx = np.mean(hull.points[hull.vertices,0])
			cy = np.mean(hull.points[hull.vertices,1])
			for simplex in hull.simplices:
				# print(simplex)
				AB = [(hull.points[simplex[0],0],hull.points[simplex[0],1]),(hull.points[simplex[1],0],hull.points[simplex[1],1])]
				x_coords, y_coords = zip(*AB)
				A = np.vstack([x_coords,np.ones(len(x_coords))]).T
				m, c = np.linalg.lstsq(A, y_coords,rcond=None)[0]
				# print("Line Solution is y = {m}x + {c}".format(m=m,c=c))
				
				if m * cx + c < cy:
					lpc.append([m, -1, c])
				else:
					lpc.append([-m, +1, -c])
			LPC.append(lpc)
		# print(LPC[0])
		return LPC

if __name__ == "__main__":

	filepath = 'C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/dat/Cai2/'
	instancename = 'cvxp_10_22'
	inst = InstanceReader.CvxPolygon(instancename)
	inst.read_cvxp_instance(filepath + instancename)
	SD = []
	mdl = MINLP(inst, SD, mode='M-M')
	mdl._solve()