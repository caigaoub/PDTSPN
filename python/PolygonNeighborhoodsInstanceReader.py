import gurobipy as gp
from gurobipy import GRB
from itertools import combinations
from sys import argv
import re
import math
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.patches as pch
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
from matplotlib.path import Path
from shapely.geometry import LineString
from shapely.geometry import Point
from shapely.ops import nearest_points
import random
import timeit
import alphashape
from descartes import PolygonPatch
from shapely.geometry import Polygon, mapping
''' -------------------------------------------------
	Things done: create a class concave polygon
	Date: 01/14/2023
	-------------------------------------------------
'''

class PolygonNeighborhoods:

	def __init__(self, instance):
		self._name = instance
		self._POINTS = []
		self._POLYGONs = []
		self._SEPERATORS = []
		self._SIMPLE_SEPERATORS = {}
		self._SEPERATORS_PLOT = []
		self._nb_ploygons = 0
		
	def read_instance(self, instance):
		nb_cvxp = int(re.split('/', instance)[-1].split('_')[1])
		rile = open(instance, "r")
		line_ = rile.readline()
		str_ = re.split('\t|\n', line_)
		self._depot = [float(re.split(':', str_[0])[0]),float(re.split(':', str_[0])[1])]
		# print(depot)
		line_ = rile.readline()
		self._POINTS = []
		self._POLYGONs = []
		self._nb_ploygons = 0
		while line_ != '':
			str_ = re.split('\t|\n', line_)
			pts = []
			for el in str_[1:]:
				if el != '':
					pts.append([float(re.split(':', el)[0]),float(re.split(':', el)[1])])
			if len(pts) == 0:
				break 
			polygon = Polygon(pts)
			polygon_vertices = list(polygon.exterior.coords)
			# print(polygon_vertices)
			self._POINTS.append(polygon_vertices)
			self._POLYGONs.append(polygon)

			self._nb_ploygons += 1
			line_ = rile.readline()

	def plot_polygon(self):
		fig, ax = plt.subplots()
		ax.set_aspect('equal', adjustable='box')
		# print(seps_plot)
		for sep in self._SEPERATORS_PLOT:
			plt.plot((sep[0][0],sep[0][1]), (sep[1][0],sep[1][1]), 'k-.', linewidth=0.8)

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
		plt.show()


	def generate_separators(self, nb_intvals = 20):
		# ---> Step 1: find the center of the smallest horizontal rectangle that contains all polygon vertices
		minx, maxx, miny, maxy = self._depot[0], self._depot[0], self._depot[1], self._depot[1]
		pset = [] # store all vertices of the convex polygons
		pset.append(self._depot)
		for polygon_vertices in self._POINTS:
			for point in polygon_vertices:
				maxx = maxx if maxx >= np.max(point[0]) else np.max(point[0])
				minx = minx if minx <= np.min(point[0]) else np.min(point[0])
				maxy = maxy if maxy >= np.max(point[1]) else np.max(point[1])
				miny = miny if miny <= np.min(point[1]) else np.min(point[1])
			for v in polygon_vertices:
				pset.append(v)
		center = [(maxx + minx)/2.0, (maxy+miny)/2.0]

		# ---> Step 2:		
		sep_set, seps_plot = [], []
		for i in range(1, nb_intvals):
			theta = math.pi/nb_intvals* i + 0.00001
			x1 = center[0]  + 1 * math.cos(theta)
			y1 = center[1]  + 1 * math.sin(theta)
			x2 = center[0]  + 1 * math.cos(theta + math.pi)
			y2 = center[1]  + 1 * math.sin(theta + math.pi)
			a = (y2- y1)/(x2 - x1)
			b = -1 
			c = -a * x2 + y2
			halfspace = [a, b, c]
			farp1, fardist1, farp2, fardist2 = self.farthest_points_twosides(pset, halfspace)

			if len(farp1) != 0:
				bestsep, endpoints = self.binary_search_separator(theta, farp1, center, 0, 1, sign=True)
				seps_plot.append(endpoints)
				sep_set.append(bestsep)
			if len(farp2) != 0:
				bestsep, endpoints = self.binary_search_separator(theta, farp2, center, 0, 1, sign=False)
				seps_plot.append(endpoints)
				sep_set.append(bestsep)
		# self.plot_polygon(seps_plot=seps_plot)
		self._SEPERATORS_PLOT = seps_plot
		self._SEPERATORS = sep_set
		return sep_set

	def _generate_simple_separators(self):
		# ---> Step 1: find the center of the smallest horizontal rectangle that contains all polygon vertices
		minx, maxx, miny, maxy = self._depot[0], self._depot[0], self._depot[1], self._depot[1]
		pset = [] # store all vertices of the convex polygons
		pset.append(self._depot)
		for polygon_vertices in self._POINTS:
			for point in polygon_vertices:
				maxx = maxx if maxx >= np.max(point[0]) else np.max(point[0])
				minx = minx if minx <= np.min(point[0]) else np.min(point[0])
				maxy = maxy if maxy >= np.max(point[1]) else np.max(point[1])
				miny = miny if miny <= np.min(point[1]) else np.min(point[1])
			for v in polygon_vertices:
				pset.append(v)
		center = [(maxx + minx)/2.0, (maxy+miny)/2.0]

		# ---> Step 2:		
		seps_plot = []
		linesegscale = 20
		self._SIMPLE_SEPERATORS = {}
		for ployindex, points in enumerate(self._POINTS):
			sep_set = []
			for idx1, p1 in enumerate(points):
				for idx2, p2 in enumerate(points):
					if idx2 > idx1 + 1:
						if p2[0] - p1[0] == 0:
							continue
						a = (p2[1]- p1[1])/(p2[0] - p1[0])
						b = -1 
						c = -a * p2[0] + p2[1]
						halfspace = [a, b, c]
						if not self.is_in_halfspace(center, halfspace) or not self.is_in_halfspace(self._depot, halfspace) :
							halfspace = [-a, -b, -c]
						point_on_line = [(p1[0]+p2[0])/2.0, (p1[1]+p2[1])/2.0]
						theta = math.atan(a)
						if self.is_projclosed(halfspace):
							x0 = point_on_line[0] + linesegscale * math.cos(theta)
							x1 = point_on_line[0] + linesegscale * math.cos(theta + math.pi)
							y0 = point_on_line[1] + linesegscale * math.sin(theta)
							y1 = point_on_line[1] + linesegscale * math.sin(theta + math.pi)
							seps_plot.append([[x0,x1], [y0,y1]])
							sep_set.append(halfspace)
			
			self._SIMPLE_SEPERATORS[ployindex] = sep_set

		self._SEPERATORS_PLOT = seps_plot
		self.plot_polygon()

		return sep_set	


	def binary_search_separator(self, theta, farp, center, lval, rval, sign=True):
		pivot = [0 , 0]
		pivot[0] = farp[0] 
		pivot[1] = farp[1]
		linesegscale = 20 
		x3 = pivot[0] + linesegscale * math.cos(theta)
		y3 = pivot[1] + linesegscale * math.sin(theta)
		x4 = pivot[0] + linesegscale * math.cos(theta + math.pi)
		y4 = pivot[1] + linesegscale * math.sin(theta + math.pi)
		a = (y4- y3)/(x4 - x3)
		b = -1 
		c = -a * x4 + y4
		halfspace = [a, b, c] if sign else [-a, -b, -c]
		bestsep = halfspace
		bestX = [x3, x4]
		bestY = [y3, y4]
		while rval-lval > 0.000001:
			alp = (lval + rval)/2.0
			pivot[0] = farp[0] * (1 - alp) + alp*center[0]
			pivot[1] = farp[1] * (1 - alp) + alp*center[1]
			x3 = pivot[0] + linesegscale * math.cos(theta)
			y3 = pivot[1] + linesegscale * math.sin(theta)
			x4 = pivot[0] + linesegscale * math.cos(theta + math.pi)
			y4 = pivot[1] + linesegscale * math.sin(theta + math.pi)
			a = (y4- y3)/(x4 - x3)
			b = -1 
			c = -a * x4 + y4
			halfspace = [a, b, c] if sign else [-a, -b, -c]
			flag = self.is_projclosed(halfspace)
			if flag:
				lval = alp
				bestsep = halfspace
				bestX = [x3, x4]
				bestY = [y3, y4]
			else:
				rval = alp	
		return bestsep, [bestX, bestY]


	def farthest_points_twosides(self, pset, halfspace):
		fardist1, fardist2 = 0, 0
		farp1, farp2 = [], []
		for p in pset:
			if not self.is_in_halfspace(p, halfspace):
				tmp = self.shortest_distance(p, halfspace)
				if fardist1 < tmp:
					fardist1 = tmp
					farp1 = p
			if self.is_in_halfspace(p, halfspace):
				tmp = self.shortest_distance(p, halfspace)
				if fardist2 < tmp:
					fardist2 = tmp
					farp2 = p
		return farp1, fardist1,farp2, fardist2

	def is_projclosed(self, halfspace):
		proj_closed = True
		if not self.is_in_halfspace(self._depot, halfspace):
			return False
		for polygon in self._POLYGONs:
			if not self.is_projclosed_single(polygon, halfspace):
				proj_closed = False
				break
		return proj_closed

	def is_projclosed_single(self, polygon, halfspace):	
		line = self.halfspace_2_line(halfspace)
		# print(line)
		projection_closed = True
		vertices = list(polygon.exterior.coords)
		for v in vertices:
			p = Point(v)
			if not self.is_in_halfspace([p.x, p.y], halfspace):
				nearp = self.nearestpoint_on_line(np.array(line.coords[0]), np.array(line.coords[1]), np.array([p.x, p.y]))
				# print(p.x, p.y, nearp)
				if not self.is_in_polygon(Point(nearp), polygon):
					projection_closed = False
					# print('point ', p.x, p.y, 'cannot be projected closed')
					break
		return projection_closed

	def nearestpoint_on_line(self, a, b, p):
	    ap = p - a
	    ab = b - a
	    result = a + np.dot(ap, ab) / np.dot(ab, ab) * ab
	    return result

	def halfspace_2_line(self, halfspace):
		x1 = 1
		y1 =  -(halfspace[0] * x1 + halfspace[2])/halfspace[1]
		x2 = 2
		y2 =  -(halfspace[0] * x2 + halfspace[2])/halfspace[1]
		line = LineString([(x1, y1), (x2, y2)])
		return line

	def is_in_polygon(self, point, polygon, tolerance=1e-12):
		return polygon.contains(point)

	def is_in_halfspace(self, p, halfspace):
		return True if p[0] * halfspace[0] + p[1] * halfspace[1] + halfspace[2] <= 0 else False

	def shortest_distance(self, p,halfspace):  
	    d = abs((halfspace[0] * p[0] + halfspace[1] * p[1] + halfspace[2])) / (math.sqrt(halfspace[0] **2 + halfspace[1] ** 2)) 
	    return d

		
	def is_in_separators(self, point, sep_set):
		flag = True
		for hs in sep_set:
			# print(point, hs)
			if not self.is_in_halfspace(point, hs):
				flag = False
				break
		return flag

	def create_instances_set(self):
		path = "C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/dat/Polygons/"
		NBCVXP = [14]
		for nb in NBCVXP:
			for clustercoef in range(1, 6):
				for i in range(1, 6):
					cvxpsize = 1
					self.generate_convex_polygons(nb, cvxpsize, clustercoef, path+'cavp_'+str(nb)+'_'
					 + str(clustercoef) +'_' +str(i))

	def generate_convex_polygons(self, nb_cvxp, cvxpsize, clustercoef, filename):
		x_range = 16
		y_range = 10
		max_radius  = 5
		Hulls = []
		POINTS = []
		# posx = random.uniform(0, 1)*x_range
		# posy = random.uniform(0, 1)*y_range
		angle = random.uniform(0, math.pi * 2)
		radius = random.uniform(max_radius* (0.08 * clustercoef), max_radius)
		posx = 8 + math.cos(angle)* radius
		posy = 5 + math.sin(angle)* radius

		depot = [posx, posy]
		for i in range(nb_cvxp):
			# posx = random.uniform(0, 1)*x_range
			# posy = random.uniform(0, 1)* y_range
			angle = random.uniform(0, math.pi * 2)
			radius = random.uniform(max_radius* (0.08 * clustercoef), max_radius)
			posx = 8 + math.cos(angle)* radius
			posy = 5 + math.sin(angle)* radius

			points = np.random.rand(10,2) * cvxpsize
			points[:,0] += posx
			points[:,1] += posy
			hull = ConvexHull(points)
			hull_vertices = points[hull.vertices, :]
			hull_vertices = hull_vertices.tolist()
			if i % 2 == 0:
				xs = [el[0] for el in hull_vertices]
				ys = [el[1] for el in hull_vertices]
				avgx = np.mean(xs)
				avgy = np.mean(ys)
				if len(xs) < 3:
					continue
				insertposes = random.sample(range(1, len(xs)-1), min([2,len(xs)-1]))
				sorted(insertposes, reverse=True)
				for pos in insertposes:
					currp = hull_vertices[pos]
					nextp = hull_vertices[pos + 1]
					midx  = (currp[0] + nextp[0]) / 2.0
					midy  = (currp[1] + nextp[1]) / 2.0
					newmidx = [midx * 0.85 + 0.15 * avgx]
					newmidy = [midy * 0.85 + 0.15 * avgy]
					mid = [newmidx[0], newmidy[0]]
					hull_vertices.insert(pos+1, mid)
			POINTS.append(hull_vertices)
		self.write_cvxp_instance(filename, depot, POINTS)

	def write_cvxp_instance(self, instance, depot, POINTS):
		file = open(instance,'w')
		file.write(str('{0:.3f}'.format(depot[0]))+':'+str('{0:.3f}'.format(depot[1]))+'\n')
		index = 1
		for pts in POINTS:
			file.write('cp'+str(index)+':\t')
			for p in pts:
				file.write(str('{0:.3f}'.format(p[0])) + ':' + str('{0:.3f}'.format(p[1])) + '\t')
			file.write('\n')
			index += 1
		file.close()


if __name__ == "__main__":

	'''choose which instance by two parameters: nb of convex polygons and instance index  '''
	nb_plgn_ngbrs = argv[1]
	cluster_coef = argv[2]
	index = argv[3]
	# ''' create instance object '''
	instance = "C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/dat/Polygons/cavp_" + str(nb_plgn_ngbrs) + "_" + str(cluster_coef) + "_" + str(index)
	polygons = PolygonNeighborhoods('polyon_'+ 'nb_polygon= '+ nb_plgn_ngbrs +' index= '+index)
	polygons.create_instances_set()
	exit(1)

	polygons.read_instance(instance)
	# # sep_set= polygons.generate_separaors(10)
	polygons._generate_simple_separators()
	polygons.plot_polygon()


	# generate_supply_demand_1_1()
	# generate_supply_demand_M_M()


	'''ddd '''
	# TIME = {}
	# for nb_cvxp in range(28, 32, 2):
	# 	for instance_idx in [1,2,3,4,5,11,12,13,14,15,21,22,23,24,25]:
	# 		filename = 'cvxp_' + str(nb_cvxp) + '_' + str(instance_idx)
	# 		instance = "C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/dat/Cai2/" + filename
	# 		result = []
	# 		for nb_seps in [6, 9, 12, 15]:
	# 			print(nb_cvxp, instance_idx, nb_seps)
	# 			cvp = CvxPolygon('convex_polyon_'+ 'nb_cvx_polygon= '+ str(nb_cvxp) +' index= '+str(instance_idx))
	# 			cvp.read_cvxp_instance(instance)
	# 			start_time = timeit.default_timer()

	# 			sep_set= cvp.generate_separators(nb_seps)
	# 			stop_time = timeit.default_timer()
	# 			elapse_time = stop_time - start_time
	# 			TIME[filename + '_' + str(nb_seps)] = elapse_time
	# 			result.append(cvp.evaluate_separators(sep_set))
	# 		with open('./ResultsSeparatorsEvaluation/ret_' + filename,'w') as data:
	# 			data.write('ret_' +filename + ' = ' + str(result))

	# print(TIME)











			