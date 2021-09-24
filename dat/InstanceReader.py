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

''' -------------------------------------------------
	Things done: create a class CvxPolygon
	Date: 9/21/2021
	-------------------------------------------------
'''

class CvxPolygon:

	def __init__(self, instance):
		self._name = instance
		self.name = ''
		self._POINTS = []
		self._HULLs = []
		self._SEPs = []
		self._totalCvxPolygonarea = 0
	def read_cvxp_instance(self, instance):
		nb_cvxp = int(re.split('/', instance)[-1].split('_')[1])
		rile = open(instance, "r")
		line_ = rile.readline()
		str_ = re.split('\t|\n', line_)
		self._depot = [float(re.split(':', str_[0])[0]),float(re.split(':', str_[0])[1])]
		# print(depot)
		line_ = rile.readline()
		self._POINTS = []
		self._HULLs = []
		self._totalCvxPolygonarea = 0
		while line_ != '':
			str_ = re.split('\t|\n', line_)
			pts = []
			for el in str_[1:]:
				if el != '':
					pts.append([float(re.split(':', el)[0]),float(re.split(':', el)[1])])
			
			self._POINTS.append(np.array(pts))
			hull = ConvexHull(np.array(pts))
			self._totalCvxPolygonarea += hull.volume
			self._HULLs.append(hull)
			line_ = rile.readline()

	def plot_hull(self, seps_plot=[]):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_aspect('equal', adjustable='box')
		for sep in seps_plot:
			plt.plot(sep[0], sep[1],'r', alpha=0.5)

		minx, maxx, miny, maxy = self._depot[0], self._depot[0], self._depot[1], self._depot[1]
		pset = [] # store all vertices of the convex polygons
		for hull in self._HULLs:
			for simplex in hull.simplices:
				maxx = maxx if maxx >= np.max(hull.points[hull.vertices,0]) else np.max(hull.points[hull.vertices,0])
				minx = minx if minx <= np.min(hull.points[hull.vertices,0]) else np.min(hull.points[hull.vertices,0])
				maxy = maxy if maxy >= np.max(hull.points[hull.vertices,1]) else np.max(hull.points[hull.vertices,1])
				miny = miny if miny <= np.min(hull.points[hull.vertices,1]) else np.min(hull.points[hull.vertices,1])
		ax.set_xlim([minx-1, maxx+1])
		ax.set_ylim([miny-1, maxy+1])
		plt.plot(self._depot[0], self._depot[1],'rs',ms=10)
		itr, cx, cy = 1, 0, 0
		for hull in self._HULLs:
			#Plot convex hull
			for simplex in hull.simplices:
				plt.plot(hull.points[simplex, 0], hull.points[simplex, 1], 'k-')
				# centroid
				cx = np.mean(hull.points[hull.vertices,0])
				cy = np.mean(hull.points[hull.vertices,1])
					#Plot centroid
				plt.plot(cx, cy,'ko',ms=2)
			ax.annotate(str(itr), xy=(cx, cy),textcoords="offset points", xytext=(cx, cy),size=14)
			itr += 1
		plt.show()

	

	def generate_separators(self, nb_intvals = 20 ):
		# ---> Step 1: find the center of the smallest horizontal rectangle that contains all polygon vertices
		minx, maxx, miny, maxy = self._depot[0], self._depot[0], self._depot[1], self._depot[1]
		pset = [] # store all vertices of the convex polygons
		pset.append(self._depot)
		for hull in self._HULLs:
			for simplex in hull.simplices:
				maxx = maxx if maxx >= np.max(hull.points[hull.vertices,0]) else np.max(hull.points[hull.vertices,0])
				minx = minx if minx <= np.min(hull.points[hull.vertices,0]) else np.min(hull.points[hull.vertices,0])
				maxy = maxy if maxy >= np.max(hull.points[hull.vertices,1]) else np.max(hull.points[hull.vertices,1])
				miny = miny if miny <= np.min(hull.points[hull.vertices,1]) else np.min(hull.points[hull.vertices,1])
			for v in hull.points[hull.vertices]:
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
		self.plot_hull(seps_plot=seps_plot)
		self._SEPs = sep_set
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
			flag = self.is_projclosed(self._HULLs, self._depot, halfspace)
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

	def is_projclosed(self, HULL, depot, halfspace):
		proj_closed = True
		if not self.is_in_halfspace(depot, halfspace):
			return False
		for hull in HULL:
			if not self.is_projclosed_single(hull, halfspace):
				proj_closed = False
				break
		return proj_closed

	def is_projclosed_single(self, hull, halfspace):	
		line = self.halfspace_2_line(halfspace)
		# print(line)
		projection_closed = True
		for v in hull.vertices:
			p = Point(hull.points[v])
			if not self.is_in_halfspace([p.x, p.y], halfspace):
				nearp = self.nearestpoint_on_line(np.array(line.coords[0]), np.array(line.coords[1]), np.array([p.x, p.y]))
				# print(p.x, p.y, nearp)
				if not self.is_in_hull(nearp, hull):
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

	def is_in_hull(self, point, hull, tolerance=1e-12):
		return all((np.dot(eq[:-1], point) + eq[-1] <= tolerance) for eq in hull.equations)

	def is_in_halfspace(self, p, halfspace):
		return True if p[0] * halfspace[0] + p[1] * halfspace[1] + halfspace[2] <= 0 else False

	def shortest_distance(self, p,halfspace):  
	    d = abs((halfspace[0] * p[0] + halfspace[1] * p[1] + halfspace[2])) / (math.sqrt(halfspace[0] **2 + halfspace[1] ** 2)) 
	    return d

	def evaluate_separators(self, sep_set):
		gridsize = 1000
		minx, maxx, miny, maxy = self._depot[0], self._depot[0], self._depot[1], self._depot[1]
		for hull in self._HULLs:
			for simplex in hull.simplices:
				maxx = maxx if maxx >= np.max(hull.points[hull.vertices,0]) else np.max(hull.points[hull.vertices,0])
				minx = minx if minx <= np.min(hull.points[hull.vertices,0]) else np.min(hull.points[hull.vertices,0])
				maxy = maxy if maxy >= np.max(hull.points[hull.vertices,1]) else np.max(hull.points[hull.vertices,1])
				miny = miny if miny <= np.min(hull.points[hull.vertices,1]) else np.min(hull.points[hull.vertices,1])
		gridarea = (maxx -minx)*(maxy-miny)
		unitgridarea = gridarea/gridsize**2
		estimated_area = 0
		cutoff_area = 0
		for i in range(0,gridsize):
			px = minx + (maxx - minx)/gridsize * i
			for j in range(0, gridsize):
				py = miny + (maxy - miny)/gridsize * j
				for hull in self._HULLs:
					if self.is_in_hull([px,py], hull):
						estimated_area += unitgridarea
						if not self.is_in_separators([px, py], sep_set):
							cutoff_area += unitgridarea

		# print("grid area=", gridarea, "True polygon area = ", self._totalCvxPolygonarea, "estimated polygon area = ",estimated_area, "cutff = ", cutoff_area)
		return self._totalCvxPolygonarea, estimated_area, cutoff_area

		
	def is_in_separators(self, point, sep_set):
		flag = True
		for hs in sep_set:
			# print(point, hs)
			if not self.is_in_halfspace(point, hs):
				flag = False
				break
		return flag

	def generate_convex_polygons(self, nb_cvxp, cvxpsize, filename):
		x_range = 16
		y_range = 10
		Hulls = []
		POINTS = []
		posx = random.uniform(0, 1)*x_range
		posy = random.uniform(0, 1)*y_range
		depot = [posx, posy]
		for i in range(nb_cvxp):
			posx = random.uniform(0, 1)*x_range
			posy = random.uniform(0, 1)* y_range
			points = np.random.rand(10,2) * cvxpsize
			points[:,0] += posx
			points[:,1] += posy
			# print(points)
			POINTS.append(points)
			hull = ConvexHull(points)
			Hulls.append(hull)
		pset = POINTS[0]
		self.write_cvxp_instance(filename, depot, POINTS)
		# print(pset)
		for pts in POINTS[1:]:
			# print(pts)
			pset = np.concatenate((pset,pts), axis=0)
		# print(pset)
		big_hull = ConvexHull(pset)
		# plot_hull(Hulls, big_hull)
		# self.plot_hull_temp(depot, Hulls)

	def plot_hull_temp(self, depot, Hulls):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_aspect('equal', adjustable='box')
		minx, maxx, miny, maxy = depot[0], depot[0], depot[1], depot[1]
		pset = [] # store all vertices of the convex polygons
		for hull in Hulls:
			for simplex in hull.simplices:
				maxx = maxx if maxx >= np.max(hull.points[hull.vertices,0]) else np.max(hull.points[hull.vertices,0])
				minx = minx if minx <= np.min(hull.points[hull.vertices,0]) else np.min(hull.points[hull.vertices,0])
				maxy = maxy if maxy >= np.max(hull.points[hull.vertices,1]) else np.max(hull.points[hull.vertices,1])
				miny = miny if miny <= np.min(hull.points[hull.vertices,1]) else np.min(hull.points[hull.vertices,1])
		ax.set_xlim([minx-1, maxx+1])
		ax.set_ylim([miny-1, maxy+1])
		plt.plot(depot[0], depot[1],'rs',ms=10)
		itr, cx, cy = 1, 0, 0
		for hull in Hulls:
			#Plot convex hull
			for simplex in hull.simplices:
				plt.plot(hull.points[simplex, 0], hull.points[simplex, 1], 'k-')
				# centroid
				cx = np.mean(hull.points[hull.vertices,0])
				cy = np.mean(hull.points[hull.vertices,1])
					#Plot centroid
				plt.plot(cx, cy,'ko',ms=2)
			ax.annotate(str(itr), xy=(cx, cy),textcoords="offset points", xytext=(cx, cy),size=14)
			itr += 1
		plt.show()

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

def create_instances_set():
	path = "C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/dat/Cai2/"
	NBCVXP = [6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
	for nb in NBCVXP:
		for i in range(1, 11):
			cvxpsize = 1
			cvp = CvxPolygon('convex_polyon')
			cvp.generate_convex_polygons(nb, cvxpsize, path+'cvxp_'+str(nb)+'_'+str(i))
		for i in range(11, 21):
			cvxpsize = 2
			cvp = CvxPolygon('convex_polyon')

			cvp.generate_convex_polygons(nb, cvxpsize, path+'cvxp_'+str(nb)+'_'+str(i))
		for i in range(21, 31):
			cvxpsize = 3
			cvp = CvxPolygon('convex_polyon')
			cvp.generate_convex_polygons(nb, cvxpsize, path+'cvxp_'+str(nb)+'_'+str(i))




if __name__ == "__main__":

	'''choose which instance by two parameters: nb of convex polygons and instance index  '''
	nb_cvxp = argv[1]
	index = argv[2]
	''' create instance object '''
	instance = "C:/Users/caiga/Dropbox/Box_Research/Projects/CETSP/CETSP_Code/CETSP/dat/Cai2/cvxp_" + nb_cvxp + "_" + index
	cvp = CvxPolygon('convex_polyon_'+ 'nb_cvx_polygon= '+ nb_cvxp +' index= '+index)
	cvp.read_cvxp_instance(instance)
	sep_set= cvp.generate_separators(10)
	cvp.evaluate_separators(sep_set)
	cvp.plot_hull()


	# create_instances_set()