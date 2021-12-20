''' python package '''
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
import timeit

''' self-defined python class/funcions ''' 
import InstanceReader
import CutPool

'''
	@Class PDTSPNModel is defined for solving 1-1/M-M mode Pickup and Delivery Traveling Salesman Problem 
	with Neighborhoods (PDTSPN). 
	@Author Cai Gao
	@Date 9/1/2020
 '''

class PDTSPNModel:
	''' Constructor '''
	def __init__(self, inst, SD, mode='M-M'):
		self._instance 		= inst 				# instance name
		self._SD 			= SD                # supply-demand list
		self._mode 			= mode              # problem mode: CETSP, 1-1 or M-M
		self._depot 		= inst._depot       # depot position 
		self._HULLs 		= inst._HULLs       # convex neighborhoods  
		self._SEPs 			= inst._SEPs        # boundary-projection-closed separators
		self._nb_tars 		= inst._nb_cvxps    # number of neighborhoods/targets
		self._BPC_sep_flag 	= False 			# 
		

	''' add boundary-projection-closed separators to the model '''
	def _add_BPC_separators(self, sep_set):
		if len(sep_set) != 0:
			self._BPC_sep_flag		= True		# flag of adding separators 
			self._BPC_separators 	= sep_set

	''' Mixed-integer nonlinear programming model, solved by Generalized Benders Decomposition` '''
	def _solve(self):
		try:
			'''Create a new model'''
			model = gp.Model("PDTSPNModel_GBD")

			'''Set up GRBmodel parameters '''
			self._set_parameter(model)

			'''Initialize GRBmodel output'''
			model._nb_SECs, model._nb_SDCs, model._nb_GBCs = 0, 0, 0
			model._time_SECs, model._time_SDCs, model._time_GBCs = 0, 0, 0

			model._depot1 		= self._depot
			model._depot2 		= self._depot
			model._size 		= self._nb_tars + 2  
			model._SD 			= self._SD
			model._zbar 		= self._calc_zbar()
			model._LPC 			= self._convert_HULL_LPConsts()
			model._mode 		= self._mode
			model._BPC_sep_flag = self._BPC_sep_flag

			if self._BPC_sep_flag:
				model._BPC_sep_flag = self._BPC_sep_flag
				model._BPC_separators = self._BPC_separators

			'''Add binary variables to the model'''
			self._add_vars(model)

			'''Create objective function '''
			self._set_objective(model)
			
			'''Add constraints to the model'''
			self._set_constraints(model)
			
			if model._mode == '1-1': self._set_precedence_constraints(model) 			

			'''Optimize the model with Callback function ''' 
			model.optimize(PDTSPNModel.callback_GBD)

			'''Terminate the model and output the results '''
			self._terminate(model)
			'''Plot the optimal tour '''

			# self._plot_wOptTour_M_M(model._opt_seq, model._opt_tour) if model._mode == 'M-M' else self._plot_wOptTour_1_1(model._opt_seq, model._opt_tour)

		except gp.GurobiError as e:
			print('Error code ' + str(e.errno) + ': ' + str(e))
		except AttributeError:
			print('Encountered an attribute error')

	''' Set Gurobi model parameters '''
	def _set_parameter(self, model):
		model.setParam(GRB.Param.OutputFlag, 1)
		model.setParam(GRB.Param.TimeLimit, 3600.0)
		model.Params.lazyConstraints = 1

	''' add model variables '''
	def _add_vars(self, model):
		model._varE = model.addVars(model._size, model._size, vtype=GRB.BINARY, name="E")
		model._theta = model.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="theta")
		model.update()

	''' set model objective function '''
	def _set_objective(self, model):
		obj = 0
		for i in range(model._size):
			for j in range(model._size):
				obj += model._zbar[i,j] * model._varE[i,j] 
		obj += 1.0 * model._theta
		model.setObjective(obj, GRB.MINIMIZE)
		model.update()

	''' add TSP-related constraints  '''
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
	
	''' add preceduence constraints to the model if model is in 1-1 mode'''
	def _set_precedence_constraints(self, model):
		nb_pairs = len(model._SD) - 2					
		model._varF = model.addVars(nb_pairs, model._size, model._size, lb=0.0, ub=1.0, vtype= GRB.CONTINUOUS, name="F")
		for k in range(nb_pairs):
			for i in range(1, model._size-1):
				constrIN = 0
				constrOUT = 0
				for j in range(1, model._size-1):
					constrOUT += model._varF[k, i, j]
					constrIN += model._varF[k, j, i]
				if i == model._SD[k+1][0]:
					model.addConstr(constrOUT == 1, 's_F'+str(k)+'_'+str(i))
				if i == model._SD[k+1][1]:
					model.addConstr(constrIN  == 1, 't_F'+str(k)+'_'+str(i))
				if i != model._SD[k+1][0] and i != model._SD[k+1][1]:
					model.addConstr(constrIN - constrOUT == 0, '!st_F'+str(k)+'_'+str(i))
			model.update()

		for k in range(nb_pairs):
			for i in range(model._size):
				for j in range(model._size):
					model.addConstr(model._varF[k,i,j] <= model._varE[i,j], 'F_UB'+ str(k)+'_'+str(i)+'_'+str(j))
			model.update()


	''' terminate the model if the model reaches optimality or time limit '''
	def _terminate(self, model):
		if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
			curObj = np.zeros((model._size, model._size))
			for i in range(0, model._size):
				for j in range(0, model._size):
					curObj[i,j] = model._varE[i,j].x
					# print(varE[i,j].x, end=' ')
				# print('\n',end='')

			cirset, curBestSeq = PDTSPNModel.find_all_subtours(curObj, model._size)
			# print(curBestSeq)

			optTour, objLen = CutPool.solve_SOCP_CvxPolyNgbs(model, curBestSeq)
			model._opt_seq = curBestSeq
			model._opt_tour = optTour
			model._opt_objval = objLen
			self._write_results_to_file(model, self._instance._name)
			# self._print_results(model)
			
	''' use a dictionary to hold all results '''
	def _print_results(self, model):
		ret = {'nb_SECs': model._nb_SECs, 
			   'nb_SDCs': model._nb_SDCs, 
			   'nb_GBCs': model._nb_GBCs,
			   'Time_SECs': model._time_SECs,
			   'Time_SDCs': model._time_SDCs,
			   'Time_GBCs': model._time_GBCs,
			   'Model Time': model.Runtime,
			   'MIPGap': model.MIPGap,
			   'ObjVal': model._opt_objval}
		print(ret)

	def _write_results_to_file(self, model, filename):
		out = {'nb_SECs': model._nb_SECs, 
			   'nb_SDCs': model._nb_SDCs, 
			   'nb_GBCs': model._nb_GBCs,
			   'Time_SECs': model._time_SECs,
			   'Time_SDCs': model._time_SDCs,
			   'Time_GBCs': model._time_GBCs,
			   'Model Time': model.Runtime,
			   'MIPGap': model.MIPGap,
			   'ObjVal': model._opt_objval}
		with open('./ResultsMode_M-M/ret_' + filename,'w') as data:
			data.write('ret_' +filename + ' = ' + str(out))

	''' Generalized benders cuts added for different modes - CETSP, 1-1 and M-M '''
	@staticmethod		
	def callback_GBD(model, where):
		if where == GRB.Callback.MIPSOL:
			vals = model.cbGetSolution(model._varE)
			if False:
				for i in range(model._size):
					for j in range(model._size):
						print(vals[i,j],end=' ')
					print(end='\n')
			
			
			# Cirset is the list of subtours if exists. Otherwise, it is empty list 
			# fseq is a list of full feasible hamiltonian path if exists, otherwise, it is empty list
			
			Cirset,fseq = PDTSPNModel.find_all_subtours(vals, model._size)			
			
			if model._mode == 'CETSP':
				if len(Cirset) != 0:
					PDTSPNModel.generate_smallest_SECs(Cirset, model)
				if len(fseq) != 0:
					mu0, mu1, obj = CutPool.generate_GBC_CvxPolyNgbs(model, fseq)
					constr = 0
					for el in mu1:
						constr += el[2] * (model._varE[el[0],el[1]] - 1.0)
					model.cbLazy(obj + constr <= model._theta)
					model._nb_GBCs += 1

			if model._mode == 'M-M':
				feasible_tour = PDTSPNModel.generate_SEC_SDC_M_M(Cirset, fseq, model)
				if feasible_tour:
					start_time = timeit.default_timer()
					mu0, mu1, obj = CutPool.generate_GBC_CvxPolyNgbs(model, fseq)
					constr = 0
					for el in mu1:
						constr += el[2] * (model._varE[el[0],el[1]] - 1.0)
					model.cbLazy(obj + constr <= model._theta)
					end_time = timeit.default_timer()
					model._time_GBCs += end_time - start_time
					model._nb_GBCs += 1

			if model._mode == '1-1':
				if len(Cirset) != 0:
					start_time = timeit.default_timer()
					PDTSPNModel.generate_smallest_SECs(Cirset, model)
					end_time = timeit.default_timer()
					model._time_SECs += end_time - start_time
					model._nb_SECs += 1
				if len(fseq) != 0:
					start_time = timeit.default_timer()
					mu0, mu1, obj = CutPool.generate_GBC_CvxPolyNgbs(model, fseq)
					constr = 0
					for el in mu1:
						constr += el[2] * (model._varE[el[0],el[1]] - 1.0)
					model.cbLazy(obj + constr <= model._theta)
					end_time = timeit.default_timer()
					model._time_GBCs += end_time - start_time
					model._nb_GBCs += 1


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

	@staticmethod
	def generate_SEC_SDC_M_M(Cirset, fseq, model):
		feasible_tour = True
		''' case I: subtour exists '''
		if len(Cirset) != 0: 
			start_time = timeit.default_timer()
			subtour0 = Cirset[0]
			feasible_tour = False
			sum_sd, constr = 0, 0
			s = 0
			while s < len(subtour0)-1:
				sum_sd += model._SD[subtour0[s+1]]
				constr += model._varE[subtour0[s],subtour0[s+1]]
				if sum_sd < 0:
					model.cbLazy(constr <= s)
					model._nb_SDCs += 1
					break
				s += 1
			end_time = timeit.default_timer()
			model._time_SDCs += end_time - start_time
			# if s == len(subtour0)-1:
			# 	constr += model._varE[subtour0[-1], subtour0[0]]
			# 	model.cbLazy(constr <= len(subtour0)-1)
			# 	model._nb_SECs += 1	

			start_time = timeit.default_timer()
			PDTSPNModel.generate_smallest_SECs(Cirset, model)
			end_time = timeit.default_timer()
			model._time_SECs += end_time - start_time
			model._nb_SECs += 1
		''' case II: a full-length Hamiltonian path '''
		if len(fseq) != 0:
			sum_sd, constr = 0, 0
			s = 0
			start_time = timeit.default_timer()

			while s < model._size-1:
				sum_sd += model._SD[fseq[s+1]]
				constr += model._varE[fseq[s],fseq[s+1]]
				if sum_sd < 0:
					feasible_tour = False
					model.cbLazy(constr <= s)
					model._nb_SDCs += 1
					break
				s += 1
			end_time = timeit.default_timer()
			model._time_SDCs += end_time - start_time
		return feasible_tour

	@staticmethod
	def generate_SEC_SDC_1_1(Cirset, fseq, model):
		feasible_tour = True
		''' case I: subtour exists '''
		if len(Cirset) != 0: 
			start_time = timeit.default_timer()
			subtour0 = Cirset[0]
			feasible_tour = False
			sum_sd, constr = 0, 0
			s = 0
			while s < len(subtour0)-1:
				sum_sd += model._SD[subtour0[s+1]]
				constr += model._varE[subtour0[s],subtour0[s+1]]
				if sum_sd < 0:
					model.cbLazy(constr <= s)
					model._nb_SDCs += 1
					break
				s += 1
			end_time = timeit.default_timer()
			model._time_SDCs += end_time - start_time

			start_time = timeit.default_timer()
			PDTSPNModel.generate_smallest_SECs(Cirset, model)
			end_time = timeit.default_timer()
			model._time_SECs += end_time - start_time

		''' case II: a full-length Hamiltonian path '''
		if len(fseq) != 0:
			sum_sd, constr = 0, 0
			s = 0
			start_time = timeit.default_timer()

			while s < model._size-1:
				sum_sd += model._SD[fseq[s+1]]
				constr += model._varE[fseq[s],fseq[s+1]]
				if sum_sd < 0:
					feasible_tour = False
					model.cbLazy(constr <= s)
					model._nb_SDCs += 1
					break
				s += 1
			end_time = timeit.default_timer()
			model._time_SDCs += end_time - start_time
		return feasible_tour

	@staticmethod
	def generate_smallest_SECs(Cirset, model):
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

	''' calculate the distance between neighborhoods, which is used as objective coefficient '''
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

	def _plot_wOptTour_M_M(self, curBestSeq, optTour):
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
		accLoad = []
		sum_sd = 0
		for el in curBestSeq:
			sum_sd += self._SD[el]
			accLoad.append(sum_sd)
		plt.title('Optimal tour = ' + str(curBestSeq))
		for hull in self._HULLs:
			for simplex in hull.simplices:
				plt.plot(hull.points[simplex, 0], hull.points[simplex, 1], 'k-')
				# centroid
				cx = np.mean(hull.points[hull.vertices,0])
				cy = np.mean(hull.points[hull.vertices,1])
					#Plot centroid
				plt.plot(cx, cy,'ko',ms=2)
			# ax.annotate(str(itr), xy=(cx, cy),textcoords="offset points", xytext=(cx, cy),size=14)
			bbox = dict(boxstyle ="round", fc ="0.8")
			cl = accLoad[curBestSeq.index(itr)]
			ax.annotate('id='+str(itr)+'\nsd=' + str(self._SD[itr])+'\ncl='+str(cl), xy=(cx+0.1, cy+0.1),textcoords="offset points", bbox = bbox, xytext=(cx, cy),size=10)

			itr += 1
		plt.plot(optTour[:,0],optTour[:,1], linestyle='-', color = 'blue', markersize=3, lw=2)
		plt.show()

	def _plot_wOptTour_1_1(self, curBestSeq, optTour):
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
		plt.title('Optimal tour = ' + str(curBestSeq) + '\nSD pair=' + str(self._SD))
		for hull in self._HULLs:
			for simplex in hull.simplices:
				plt.plot(hull.points[simplex, 0], hull.points[simplex, 1], 'k-')
				# centroid
				cx = np.mean(hull.points[hull.vertices,0])
				cy = np.mean(hull.points[hull.vertices,1])
					#Plot centroid
				plt.plot(cx, cy,'ko',ms=2)
			# ax.annotate(str(itr), xy=(cx, cy),textcoords="offset points", xytext=(cx, cy),size=14)
			bbox = dict(boxstyle ="round", fc ="0.8")
			ax.annotate('id='+str(itr), xy=(cx+0.1, cy+0.1),textcoords="offset points", bbox = bbox, xytext=(cx, cy),size=10)

			itr += 1
		plt.plot(optTour[:,0],optTour[:,1], linestyle='-', color = 'blue', markersize=3, lw=2)
		plt.show()
	
	''' Convert the convex hull as linear constraints '''
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
	mode = 'run-all-instances-of-Mode-M-M'

	if mode  == 'test':
		problem_type = 'M-M'
		''' Step 1: read instances '''
		filepath = 'C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/dat/Cai2/'
		instancename = 'cvxp_20_4'
		inst = InstanceReader.CvxPolygon(instancename)
		inst.read_cvxp_instance(filepath + instancename)

		SD_cvxp10 = [0,3, 2, 3, -1, -3, 5, -4, 4, -5, -3, 0]
		SD_cvxp20 = [0, 5, 2, 3, -1, -3, 5, -10, 4, -8, 3, 5, 2, 3, -1, -3, 5, -10, 4, -8, 3, 0]
		SD_cvxp24 = [0, 5, 2, 3, 1, -5, 5, -10, 4, -8, 3, 5, 2, 3, +1, -5, 5, -10, 4, -16, 3, 6, 4, -4, 2, 0]
		SD_cvxp10_pair = [0, (1,4), (9,2), (3, 7), (5, 8), (6,10), 11]
		SD_cvxp20_pair = [0, (1,4), (2,9), (3, 7), (5, 8), (6,10), (11,14), (12,19), (13, 17), (15, 18), (16,20), 21]
		SD_cvxp30_pair = [0, (1,4), (2,9), (3, 7), (5, 8), (6,10), (11,14), (12,19), (13, 17), (15, 18), (16,20), (21,24), (22,29), (23, 27), (25, 28), (26,30), 31]
		SD_cvxp26_pair = [0, (1,4), (2,9), (3, 7), (5, 8), (6,10), (11,14), (12,19), (13, 17), (15, 18), (16,20), (21,24), (22,25), (23, 26), 27]
		SD_cvxp24_pair = [0, (1,4), (2,9), (3, 7), (5, 8), (6,10), (11,14), (12,19), (13, 17), (15, 18), (16,20), (21,24), (22,23), 25]

		''' Step 2: create model '''
		mdl = PDTSPNModel(inst, SD_cvxp20, problem_type)

		''' Step 3: generate boundary-projection-closed separators '''
		sep_set = inst.generate_separators(7)
		mdl._add_BPC_separators(sep_set)
		
		''' Step 4: solve the model '''
		mdl._solve()

	if mode == 'run-single-instance-Mode1-1':
		problem_type = '1-1'

		''' Step 1: read instances '''
		nb_cvxp = 6
		instance_idx = 4
		
		instancename = 'cvxp_'+str(nb_cvxp)+'_' + str(instance_idx)
		inst = InstanceReader.CvxPolygon(instancename)
		filepath = 'C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/dat/Cai2/'
		inst.read_cvxp_instance(filepath + instancename)
		import SupplyDemandPairs

		SD_cvxp10 = SupplyDemandPairs.AllSDPairs['cp_'+str(nb_cvxp)+'_idx_'+str(instance_idx)]
		''' Step 2: create model '''
		mdl = PDTSPNModel(inst, SD_cvxp10, problem_type)

		''' Step 3: generate boundary-projection-closed separators '''
		sep_set = inst.generate_separators(7)
		mdl._add_BPC_separators(sep_set)
		
		''' Step 4: solve the model '''
		mdl._solve()

	if mode == 'run-all-instances-of-Mode-1-1':
		problem_type = '1-1'
		for nb_cvxp in range(20, 28, 2):
			for instance_idx in [1,2,3,4,5,11,12,13,14,15,21,22,23,24,25]:
				print('running instance ', nb_cvxp, ' - ', instance_idx)
				''' Step 1: read instances '''				
				instancename = 'cvxp_'+str(nb_cvxp)+'_' + str(instance_idx)
				inst = InstanceReader.CvxPolygon(instancename)
				filepath = 'C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/dat/Cai2/'
				inst.read_cvxp_instance(filepath + instancename)
				import SupplyDemandPairs

				SD_cvxp10 = SupplyDemandPairs.AllSDPairs['cp_'+str(nb_cvxp)+'_idx_'+str(instance_idx)]
				''' Step 2: create model '''
				mdl = PDTSPNModel(inst, SD_cvxp10, problem_type)

				''' Step 3: generate boundary-projection-closed separators '''
				sep_set = inst.generate_separators(7)
				mdl._add_BPC_separators(sep_set)
				
				''' Step 4: solve the model '''
				mdl._solve()

	if mode == 'run-all-instances-of-Mode-M-M':
		problem_type = 'M-M'
		for nb_cvxp in range(24, 26, 2):
			for instance_idx in [24,25]:
				print('running instance ', nb_cvxp, ' - ', instance_idx)
				''' Step 1: read instances '''				
				instancename = 'cvxp_'+str(nb_cvxp)+'_' + str(instance_idx)
				inst = InstanceReader.CvxPolygon(instancename)
				filepath = 'C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/dat/Cai2/'
				inst.read_cvxp_instance(filepath + instancename)
				import SupplyDemand_M_M

				SD_cvxp10 = SupplyDemand_M_M.AllSDPairs['cp_'+str(nb_cvxp)+'_idx_'+str(instance_idx)]
				''' Step 2: create model '''
				mdl = PDTSPNModel(inst, SD_cvxp10, problem_type)

				''' Step 3: generate boundary-projection-closed separators '''
				sep_set = inst.generate_separators(7)
				mdl._add_BPC_separators(sep_set)
				
				''' Step 4: solve the model '''
				mdl._solve()