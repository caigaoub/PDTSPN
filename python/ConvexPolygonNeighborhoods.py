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

def generate_convex_polygons(nb_cvxp, filename):
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
		points = np.random.rand(10,2)*4.0
		points[:,0] += posx
		points[:,1] += posy
		# print(points)
		POINTS.append(points)
		hull = ConvexHull(points)
		Hulls.append(hull)
	pset = POINTS[0]
	write_cvxp_instance(filename, depot, POINTS)
	# print(pset)
	for pts in POINTS[1:]:
		# print(pts)
		pset = np.concatenate((pset,pts), axis=0)
	# print(pset)
	big_hull = ConvexHull(pset)
	# plot_hull(Hulls, big_hull)
	# plot_hull(depot, Hulls)

def read_cvxp_instance(instance):
	nb_cvxp = int(re.split('/', instance)[-1].split('_')[1])
	rile = open(instance, "r")
	line_ = rile.readline()
	# print(re.split('\t|\n', line_))
	str_ = re.split('\t|\n', line_)
	depot = [float(re.split(':', str_[0])[0]),float(re.split(':', str_[0])[1])]
	# print(depot)
	line_ = rile.readline()
	POINTS = []
	HULLs = []
	while line_ != '':
		str_ = re.split('\t|\n', line_)
		pts = []
		for el in str_[1:]:
			if el != '':
				pts.append([float(re.split(':', el)[0]),float(re.split(':', el)[1])])
		POINTS.append(np.array(pts))
		hull = ConvexHull(np.array(pts))
		HULLs.append(hull)
		line_ = rile.readline()
	# print(HULL)

	return depot, HULLs

def convert_HULL_LPConsts(HULLs):
	LPC = []
	for hull in HULLs:
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


def write_cvxp_instance(instance, depot, POINTS):
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

# def plot_hull(Hulls, big_hull):
def plot_hull(depot, Hulls, seps_plot):

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_aspect('equal', adjustable='box')
	# ax.set_xlim([0, 4])
	# ax.set_ylim([0, 4])
	for sep in seps_plot:
		plt.plot(sep[0], sep[1])


	# for simplex in big_hull.simplices:
	# 	plt.plot(big_hull.points[simplex, 0], big_hull.points[simplex, 1], 'r--')
	# 	cx = np.mean(big_hull.points[big_hull.vertices,0])
	# 	cy = np.mean(big_hull.points[big_hull.vertices,1])
	# 	plt.plot(cx, cy,'rs',ms=10)
	plt.plot(depot[0], depot[1],'rs',ms=10)

	itr = 1
	cx = 0
	cy = 0
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

def plot_CvxPoly_Instance_wOptTour(depot, Hulls, optTour):
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

	plt.plot(depot[0], depot[1],'rs',ms=10)
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
	plt.plot(optTour[:,0],optTour[:,1], linestyle='-', color = 'blue', markersize=3, lw=2)
	plt.show()


def calc_zbar(depot, HULLs):
	nb_cvxp = len(HULLs)
	zbar = np.zeros((nb_cvxp+2,nb_cvxp+2))
	for i in range(0, nb_cvxp):
		hull = HULLs[i]

		a0 = np.array([depot[0], depot[1], 0])
		a1 = np.array([depot[0]+0.0001, depot[1]+0.0001, 0])
		dist = []
		for simplex in hull.simplices:
			b0 = np.array([hull.points[simplex[0]][0], hull.points[simplex[0]][1], 0])
			b1 = np.array([hull.points[simplex[1]][0], hull.points[simplex[1]][1], 0])
			pa, pb, md = closestDistanceBetweenLines(a0, a1, b0, b1, clampAll=True)
			dist.append(md)
		# print(dist)
		zbar[0, i+1] = min(dist)
		zbar[i+1, 0] = zbar[0, i+1]
		zbar[i+1, nb_cvxp+1] = zbar[0, i+1]
		zbar[nb_cvxp+1, i+1] = zbar[0, i+1]

	for i in range(0, nb_cvxp):
		for j in range(0, nb_cvxp):
			if i != j:
				hulli = HULLs[i]
				hullj = HULLs[j]
				dist = []
				for ksplx in hulli.simplices:
					a0 = np.array([hulli.points[ksplx[0]][0], hulli.points[ksplx[0]][1], 0])
					a1 = np.array([hulli.points[ksplx[1]][0], hulli.points[ksplx[1]][1], 0])
					for lsplx in hullj.simplices:
						b0 = np.array([hullj.points[lsplx[0]][0], hullj.points[lsplx[0]][1], 0])
						b1 = np.array([hullj.points[lsplx[1]][0], hullj.points[lsplx[1]][1], 0])
						pa, pb, md = closestDistanceBetweenLines(a0, a1, b0, b1, clampAll=True)
						dist.append(md)
				zbar[i+1, j+1] = min(dist)
				zbar[j+1, i+1] = zbar[i+1, j+1]
	# print(zbar)	
	return zbar


def closestDistanceBetweenLines(a0,a1,b0,b1,clampAll=False,clampA0=False,clampA1=False,clampB0=False,clampB1=False):

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



def generate_separators(HULL, depot):
	minx = depot[0]
	maxx  = depot[0]
	miny = depot[1]
	maxy = depot[1]
	pset = []
	for hull in HULL:
		for simplex in hull.simplices:
			maxx = maxx if maxx >= np.max(hull.points[hull.vertices,0]) else np.max(hull.points[hull.vertices,0])
			minx = minx if minx <= np.max(hull.points[hull.vertices,0]) else np.min(hull.points[hull.vertices,0])
			maxy = maxy if maxy >= np.max(hull.points[hull.vertices,1]) else np.max(hull.points[hull.vertices,1])
			miny = miny if miny <= np.max(hull.points[hull.vertices,1]) else np.min(hull.points[hull.vertices,1])
		for v in hull.points[hull.vertices]:
			# print(v)
			pset.append(v)
	# print(minx, maxx, miny, maxy)
	center = [(maxx + minx)/2.0, (maxy+miny)/2.0]

	# c_range = math.sqrt((maxx-minx)**2/4.0 + (maxy-miny)**2/4.0)

	nb_intvals = 30
	sep_set = []
	seps_plot = []
	for i in range(1, nb_intvals):
		theta = math.pi/nb_intvals* i + 0.00001

		x1 = center[0]  + 1 * math.cos(theta)
		y1 = center[1]  + 1 * math.sin(theta)
		x2 = center[0]  + 1 * math.cos(theta + math.pi)
		y2 = center[1]  + 1 * math.sin(theta + math.pi)
		# line = LineString([(x1, y1), (x2, y2)])
		a = (y2- y1)/(x2 - x1)
		b = -1 
		c = -a * x2 + y2

		halfspace = [a, b, c]
		nearp1, neardist1, nearp2, neardist2 = nearest_point(pset, halfspace)

		pivot1 = [0 , 0]
		pivot2 = [0 , 0]
		bestX = []
		bestY = []
		bestsep = []
		for alp in np.arange(0, 0.2,0.01):
			pivot1[0] = nearp1[0] * (1 - alp) + alp*center[0]
			pivot1[1] = nearp1[1] * (1 - alp) + alp*center[1]
			x3 = pivot1[0]  + 8 * math.cos(theta)
			y3 = pivot1[1]  +8 * math.sin(theta)
			x4 = pivot1[0]  + 8 * math.cos(theta + math.pi)
			y4 = pivot1[1]  + 8* math.sin(theta + math.pi)
			a = (y4- y3)/(x4 - x3)
			b = -1 
			c = -a * x4 + y4
			halfspace = [a, b, c]
			flag = is_projclosed(HULL, depot, halfspace)
			
			print('projection closed: ', flag)
			if flag:
				bestsep = [a, b, c]
				bestX = [x3, x4]
				bestY = [y3, y4]
			else:
				break
			# plot_hull(depot, HULL, [x3, x4], [y3, y4])
		seps_plot.append([bestX, bestY])
		sep_set.append(bestsep)
		for alp in np.arange(0, 0.1,0.01):
			pivot2[0] = nearp2[0] * (1 - alp) + alp*center[0]
			pivot2[1] = nearp2[1] * (1 - alp) + alp*center[1]
			x3 = pivot2[0]  + 4 * math.cos(theta)
			y3 = pivot2[1]  +4 * math.sin(theta)
			x4 = pivot2[0]  + 4 * math.cos(theta + math.pi)
			y4 = pivot2[1]  + 4 * math.sin(theta + math.pi)
			a = (y4- y3)/(x4 - x3)
			b = -1 
			c = -a * x4 + y4
			halfspace = [-a, -b, -c]
			flag = is_projclosed(HULL, depot, halfspace)
			
			print('projection closed: ', flag)
			if flag:
				bestsep = [-a, -b, -c]
				bestX = [x3, x4]
				bestY = [y3, y4]
			else:
				break
			# plot_hull(depot, HULL, [x3, x4], [y3, y4])
		sep_set.append(bestsep)
		seps_plot.append([bestX, bestY])
	plot_hull(depot, HULL, seps_plot)

	# print(center)


def nearest_point(pset, halfspace):
	neardist1 = 0
	neardist2 = 0
	nearp1 = []
	nearp2 = []

	for p in pset:
		if not is_in_halfspace(p, halfspace):
			tmp = shortest_distance(p, halfspace)
			if neardist1 < tmp:
				neardist1 = tmp
				nearp1 = p
		if is_in_halfspace(p, halfspace):
			tmp = shortest_distance(p, halfspace)
			if neardist2 < tmp:
				neardist2 = tmp
				nearp2 = p

	return nearp1, neardist1,nearp2, neardist2

def shortest_distance(p,halfspace):  
    d = abs((halfspace[0] * p[0] + halfspace[1] * p[1] + halfspace[2])) / (math.sqrt(halfspace[0] **2 + halfspace[1] ** 2)) 
    return d

def is_projclosed(HULL, depot, halfspace):
	proj_closed = True
	if not is_in_halfspace(depot, halfspace):
		return False
	# i = 1
	for hull in HULL:
		if not is_projclosed_single(hull, halfspace):
			proj_closed = False
			# print('hull',i, ' is NOT projection closed')
			break
		# else:
			# print('hull',i, 'is projection closed')
		# i += 1
	return proj_closed

def is_projclosed_single(hull, halfspace):	
	line = halfspace_2_line(halfspace)
	# print(line)
	projection_closed = True
	for v in hull.vertices:
		p = Point(hull.points[v])
		if not is_in_halfspace([p.x, p.y], halfspace):
			nearp = nearestpoint_on_line(np.array(line.coords[0]), np.array(line.coords[1]), np.array([p.x, p.y]))
			# print(p.x, p.y, nearp)
			if not is_in_hull(nearp, hull):
				projection_closed = False
				# print('point ', p.x, p.y, 'cannot be projected closed')
				break
	return projection_closed

def nearestpoint_on_line(a, b, p):
    ap = p - a
    ab = b - a
    result = a + np.dot(ap, ab) / np.dot(ab, ab) * ab
    return result

def halfspace_2_line(halfspace):
	x1 = 1
	y1 =  -(halfspace[0] * x1 + halfspace[2])/halfspace[1]
	x2 = 2
	y2 =  -(halfspace[0] * x2 + halfspace[2])/halfspace[1]
	line = LineString([(x1, y1), (x2, y2)])
	return line

def is_in_hull(point, hull, tolerance=1e-12):
	return all((np.dot(eq[:-1], point) + eq[-1] <= tolerance) for eq in hull.equations)

def is_in_halfspace(p, halfspace):
	return True if p[0] * halfspace[0] + p[1] * halfspace[1] + halfspace[2] <= 0 else False



if __name__ == "__main__":

	'''Generate instances '''
	nb_cvxp = argv[1]
	index = argv[2]
	# instance = "/home/latte/Dropbox/Box_Rese arch/Github/CETSP/dat/Cai/cvxp_" + nb_cvxp + "_" + str(index)
	# generate_convex_polygons(int(nb_cvxp), instance)
	
	''' plot insrtance'''
	instance = "/home/latte/Dropbox/Box_Research/Github/CETSP/dat/Cai/cvxp_" + nb_cvxp + "_" + index
	depot, HULL = read_cvxp_instance(instance)
	convert_HULL_LPConsts(HULL)
	# print(HULL[0].vertices)
	# halfspace = [-1, -1, 6.7]
	# is_projclosed_single(HULL[3], halfspace)
	# is_projclosed(HULL, halfspace)
	generate_separators(HULL,depot)
	# plot_hull(depot, HULL)
	# calc_zbar(depot, HULL)