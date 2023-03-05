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
import alphashape
from descartes import PolygonPatch
from shapely.geometry import Polygon, mapping
''' self-defined python class/funcions ''' 
import PolygonNeighborhoodsInstanceReader
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
		self._POLYGONs 		= inst._POLYGONs       # convex neighborhoods
		self._PLOYGON_POINTS = inst._POINTS 
		self._nb_tars 		= inst._nb_ploygons    # number of neighborhoods/targets
		self._BPC_sep_flag 	= False 
		self._BPC_simle_separators = inst._SIMPLE_SEPERATORS
		self._BPC_separators 	= []
		self._BPC_separators_plot = []
		self._neighborhood_types = []
		self._LPC = []


	''' add boundary-projection-closed separators to the model '''
	def _add_BPC_separators(self, sep_set, simple_separators, sep_plot):
		self._BPC_sep_flag		= True		# flag of adding separators 
		self._BPC_separators 	= sep_set
		self._BPC_simle_separators = simple_separators
		self._BPC_separators_plot = sep_plot

	def _update_convexity(self):
		self._neighborhood_types = []
		for ployindex, points in enumerate(self._PLOYGON_POINTS):
			if self._is_convex(points):
				self._neighborhood_types.append('convex')
			else:
				self._neighborhood_types.append('concave')

	''' Convert the convex ploygon as linear constraints '''
	def _convert_LPConsts_if_convex(self):
		self._LPC = []
		for pindex, polygon in enumerate(self._POLYGONs):
			lpc = []
			if self._neighborhood_types[pindex] == 'convex':
				vertices = self._PLOYGON_POINTS[pindex]
				cx, cy = [sum(x)/len(vertices) for x in zip(*vertices)]
				for vindex in range(len(self._PLOYGON_POINTS[pindex]) - 1):
					curr_v = (vertices[vindex][0], vertices[vindex][1])
					next_v = (vertices[vindex+1][0], vertices[vindex+1][1])
					AB = [curr_v, next_v]			
					x_coords, y_coords = zip(*AB)
					A = np.vstack([x_coords,np.ones(len(x_coords))]).T
					m, c = np.linalg.lstsq(A, y_coords,rcond=None)[0]
					# print("Line Solution is y = {m}x + {c}".format(m=m,c=c))
					if m * cx + c < cy:
						lpc.append([m, -1, c])
					else:
						lpc.append([-m, +1, -c])
			self._LPC.append(lpc)


	def _is_convex(self, points):
		N = len(points)
		prev = 0
		curr = 0
		for i in range(N):
			temp = [points[i], points[(i + 1) % N], points[(i + 2) % N]]
			curr = self._crossproduct(temp)
			if (curr != 0):
				if (curr * prev < 0):
					return False
				else:
					prev = curr
		return True

	def _crossproduct(self, A):
		X1 = (A[1][0] - A[0][0])
		Y1 = (A[1][1] - A[0][1])
		X2 = (A[2][0] - A[0][0])
		Y2 = (A[2][1] - A[0][1])
		return (X1 * Y2 - Y1 * X2)
	
	''' Mixed-integer nonlinear programming model, solved by Generalized Benders Decomposition` '''
	def _solve(self, use_seperators):
		try:
			'''Create a new model'''
			model = gp.Model("PDTSPNModel_GBD")

			'''Set up GRBmodel parameters '''
			self._set_parameter(model)

			'''Initialize GRBmodel output'''
			model._nb_SECs, model._nb_SDCs, model._nb_GBCs = 0, 0, 0
			model._time_SECs, model._time_SDCs, model._time_GBCs = 0, 0, 0

			if use_seperators:
				self._cut_neighborhoods_by_seperators()
			self._update_convexity()
			self._convert_LPConsts_if_convex()


			model._depot1 		= self._depot
			model._depot2 		= self._depot
			model._POLYGONs 	= self._POLYGONs
			model._size 		= self._nb_tars + 2  
			model._SD 			= self._SD
			model._zbar 		= self._calc_zbar()
			model._mode 		= self._mode
			model._BPC_sep_flag = self._BPC_sep_flag
			model._neighborhood_types = self._neighborhood_types
			model._LPC = self._LPC
			model._use_seperators = use_seperators

			print(model._neighborhood_types)
			if self._BPC_sep_flag:
				model._BPC_sep_flag = self._BPC_sep_flag
				model._BPC_separators = self._BPC_separators

			'''Add binary variables to the model'''
			self._add_vars(model)

			'''Create objective function '''
			self._set_objective(model)
			
			'''Add constraints to the model'''
			self._set_constraints(model)		

			self._set_precedence_constraints(model) 			

			'''Optimize the model with Callback function ''' 
			model.optimize(PDTSPNModel.callback_GBD)

			'''Terminate the model and output the results '''
			self._terminate(model)
			'''Plot the optimal tour '''

			# self._plot_wOptTour_M_M(model._opt_seq, model._opt_tour) if model._mode == 'M-M' else self._plot_wOptTour_1_1(model._opt_seq, model._opt_tour)

		except gp.GurobiError as e:
			print('Error code ' + str(e.errno) + ': ' + str(e))
		except AttributeError:
			print('Master problem encountered an attribute error')

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
			print("TSP sequence:", curBestSeq)

			optTour, objLen = CutPool.solve_ConcavePolyNgbs(model, curBestSeq)
			# print(optTour)
			model._opt_seq = curBestSeq
			model._opt_tour = optTour
			model._opt_objval = objLen
			self._write_results_to_file(model, self._instance._name)
			self._plot_wOptTour_generic(curBestSeq, optTour)
			self._print_results(model)
			
	''' use a dictionary to hold all results '''
	def _print_results(self, model):
		ret = {'use_seperators': model._use_seperators,
			   'nb_SECs': model._nb_SECs, 
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
		out = {'use_seperators': model._use_seperators,
			   'nb_SECs': model._nb_SECs, 
			   'nb_SDCs': model._nb_SDCs, 
			   'nb_GBCs': model._nb_GBCs,
			   'Time_SECs': model._time_SECs,
			   'Time_SDCs': model._time_SDCs,
			   'Time_GBCs': model._time_GBCs,
			   'Model Time': model.Runtime,
			   'MIPGap': model.MIPGap,
			   'ObjVal': model._opt_objval}
		with open('./MixedPolygons/ret','a') as data:
			data.write('ret_' +filename + ' = ' + str(out) + '\n')

	''' 
	Generalized benders cuts added for different modes - CETSP, 1-1 and M-M '''
	@staticmethod		
	def callback_GBD(model, where):
		if where == GRB.Callback.MIPSOL:
			vals = model.cbGetSolution(model._varE)
			if False:
				for i in range(model._size):
					for j in range(model._size):
						print(vals[i,j],end=' ')
					print(end='\n')
			
			Cirset,fseq = PDTSPNModel.find_all_subtours(vals, model._size)			
					
			if model._mode == 'MixedPolygons':
				if len(Cirset) != 0:
					start_time = timeit.default_timer()
					PDTSPNModel.generate_smallest_SECs(Cirset, model)
					end_time = timeit.default_timer()
					model._time_SECs += end_time - start_time
					model._nb_SECs += 1
				if len(fseq) != 0:
					start_time = timeit.default_timer()
					mu0, mu1, obj = CutPool.generate_GBC_MixedPolyNgbs(model, fseq)
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
		# self._nb_tars = len(self._POLYGONs)
		if self._nb_tars <= 0:
			print("the number of polygons is wrong!! ")
			exit()
			
		zbar = np.zeros((self._nb_tars+2,self._nb_tars+2))
		for i in range(0, self._nb_tars):
			polygon = self._POLYGONs[i]
			a0 = np.array([self._depot[0], self._depot[1], 0])
			a1 = np.array([self._depot[0]+0.0000001, self._depot[1]+0.00000001, 0])
			dist = []
			polygon_simplices = list(polygon.exterior.coords)
			for index in range(len(polygon_simplices) - 1):
				b0 = np.array([polygon_simplices[index][0], polygon_simplices[index][1], 0])
				b1 = np.array([polygon_simplices[index + 1][0], polygon_simplices[index + 1][1], 0])
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
					polygoni = self._POLYGONs[i]
					polygonj = self._POLYGONs[j]
					dist = []
					polygoni_simplices = list(polygoni.exterior.coords)
					polygonj_simplices = list(polygonj.exterior.coords)
					for indexi in range(len(polygoni_simplices) - 1):
						a0 = np.array([polygoni_simplices[indexi][0], polygoni_simplices[indexi][1], 0])
						a1 = np.array([polygoni_simplices[indexi + 1][0], polygoni_simplices[indexi + 1][1], 0])
						for indexj in range(len(polygonj_simplices) - 1):
							b0 = np.array([polygonj_simplices[indexj][0], polygonj_simplices[indexj][1], 0])
							b1 = np.array([polygonj_simplices[indexj + 1][0], polygonj_simplices[indexj + 1][1], 0])
							pa, pb, md = self._closestDistanceBetweenLines(a0, a1, b0, b1, clampAll=True)
							dist.append(md)
					zbar[i+1, j+1] = min(dist)
					zbar[j+1, i+1] = zbar[i+1, j+1]
		# print(zbar)
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

	def _plot_wOptTour_generic(self, curBestSeq, optTour):
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
		plt.title('Optimal tour = ' + str(curBestSeq))

		minx, maxx, miny, maxy = self._depot[0], self._depot[0], self._depot[1], self._depot[1]
		for polygon in self._POLYGONs:
			polygon_vertices = list(polygon.exterior.coords)
			for point in polygon_vertices:
				maxx = maxx if maxx >= np.max(point[0]) else np.max(point[0])
				minx = minx if minx <= np.min(point[0]) else np.min(point[0])
				maxy = maxy if maxy >= np.max(point[1]) else np.max(point[1])
				miny = miny if miny <= np.min(point[1]) else np.min(point[1])
		ax.set_xlim([minx-1, maxx+1])
		ax.set_ylim([miny-1, maxy+1])
		plt.plot(self._depot[0], self._depot[1],'ko',ms=5)
		itr, cx, cy = 1, 0, 0
		for polygon in self._POLYGONs:
			vertices = list(polygon.exterior.coords)
			cx = [p[0] for p in vertices]
			cy = [p[1] for p in vertices]
			plt.plot(cx, cy,'k-')
			itr += 1
		plt.plot(optTour[:,0],optTour[:,1], linestyle='-', color = 'blue', markersize=3, lw=2)
		for sep in self._BPC_separators_plot:
			plt.plot((sep[0][0],sep[0][1]), (sep[1][0],sep[1][1]), 'k-.', linewidth=0.8)
		plt.show()

	def _cut_neighborhoods_by_seperators(self):
		updated_ploygons = []
		updated_points = []
		# print(self._BPC_simle_separators)
		for polyindex, points in enumerate(self._PLOYGON_POINTS):
			points_to_remove = []
			cloned_points = points
			# print(cloned_points)
			simple_separators = self._BPC_simle_separators[polyindex]
			if len(simple_separators) != 0:
				# print(polyindex, simple_separators)
				while True:
					itr = 0
					for point in cloned_points:
						if self._is_outside_halfspaces(point, simple_separators):
							cloned_points.remove(point)
							break
						itr += 1

					if itr == len(cloned_points):
						break
			# print(cloned_points)
			updated_points.append(cloned_points)
			updated_ploygons.append(Polygon(cloned_points))

		self._POLYGONs = updated_ploygons
		self._PLOYGON_POINTS = updated_points


	def _is_outside_halfspaces(self, p, halfspaces):
		for halfspace in halfspaces:
			if p[0] * halfspace[0] + p[1] * halfspace[1] + halfspace[2] > 0.00001:
				return True
		return False

if __name__ == "__main__":
	# mode = 'run-all-instances-of-Mode-1-1'
	mode = 'test'

	if mode  == 'test':
		problem_type = 'MixedPolygons'
		nb_plgn_ngbrs = argv[1]
		cluster_coef = argv[2]
		index = argv[3]
		''' Step 1: read instances '''
		filepath = 'C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/dat/Polygons/'
		instancename = 'cavp_'+nb_plgn_ngbrs+'_' + cluster_coef +'_'+index
		inst = PolygonNeighborhoodsInstanceReader.PolygonNeighborhoods(instancename)
		inst.read_instance(filepath + instancename)

		import SupplyDemandPairs
		convert_2_SD_index = (int(cluster_coef)-1)*5 + int(index)
		SD = SupplyDemandPairs.AllSDPairs['cp_'+str(nb_plgn_ngbrs)+'_idx_'
		+str(convert_2_SD_index)]
		''' Step 2: create model '''
		mdl = PDTSPNModel(inst, SD, problem_type)

		''' Step 3: generate boundary-projection-closed separators '''
		# inst.generate_separators(7)
		inst._generate_simple_separators()
		use_seperators = False
		mdl._add_BPC_separators(inst._SEPERATORS, inst._SIMPLE_SEPERATORS, inst._SEPERATORS_PLOT)
		
		''' Step 4: solve the model '''
		mdl._solve(use_seperators)

	if mode  == 'run-all-instances-of-Mode-1-1':
		problem_type = 'MixedPolygons'
		for nb_plgn_ngbrs in [16]:
			for cluster_coef in [1,2,3,4,5]:
				for index in [1,2,3,4,5]:
					try:
						''' Step 1: read instances '''
						filepath = 'C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/dat/Polygons/'
						instancename = 'cavp_'+str(nb_plgn_ngbrs)+'_' + str(cluster_coef) +'_'+str(index)
						print('******************************\n')
						print(instancename)
						print('******************************\n')

						inst = PolygonNeighborhoodsInstanceReader.PolygonNeighborhoods(instancename)
						inst.read_instance(filepath + instancename)

						import SupplyDemandPairs
						convert_2_SD_index = (int(cluster_coef)-1)*5 + int(index)
						SD = SupplyDemandPairs.AllSDPairs['cp_'+str(nb_plgn_ngbrs)+'_idx_'
						+str(convert_2_SD_index)]
						''' Step 2: create model '''
						mdl = PDTSPNModel(inst, SD, problem_type)

						''' Step 3: generate boundary-projection-closed separators '''
						# inst.generate_separators(7)
						inst._generate_simple_separators()
						use_seperators = True
						mdl._add_BPC_separators(inst._SEPERATORS, inst._SIMPLE_SEPERATORS, inst._SEPERATORS_PLOT)
						
						''' Step 4: solve the model '''
						mdl._solve(use_seperators)
					except: 
						pass