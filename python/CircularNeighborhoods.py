import gurobipy as gp
from gurobipy import GRB
from itertools import combinations
from sys import argv
import re
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pch
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay



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

def plot_Behdani_instance(depot, X, Y, R, optTour):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	# plt.axis('off')
	XX = [x for x in X]
	XX.append(depot[0])
	YY = [y for y in Y]
	YY.append(depot[1])
	ax.set_xlim([min(XX)-max(R)-1, max(XX)+max(R)+1])
	ax.set_ylim([min(YY)-max(R), max(YY)+max(R)+1])
	# ax.set_xlim([-2, 15])
	# ax.set_ylim([-2, 12])
	plt.xlabel('x', fontsize=18)
	plt.ylabel('y', fontsize=18)
	ax.tick_params(axis='both', labelsize=16)
	ax.set_aspect('equal', adjustable='box')
	# plot depot
	sidelens = 0.2
	rect = pch.Rectangle((depot[0]-sidelens/2.0, depot[1]-sidelens/2.0), sidelens, sidelens, edgecolor='r',fill=True,facecolor='r') 
	ax.add_patch(rect) 
	# plot circles
	for i in range(len(X)):
		circle = plt.Circle((X[i], Y[i]), radius=R[i], edgecolor='k',fill=False, alpha = 1)
		ax.add_artist(circle)
		ax.annotate(str(i+1), xy=(X[i], Y[i]),textcoords="offset points", xytext=(X[i], Y[i]),size=14)
	# convex hull
	ls_points = [depot]
	for x, y in zip(X, Y):
		ls_points.append([x,y])
	points = np.array(ls_points)
	# print(points)
	hull = ConvexHull(points)
	# print(hull.simplices)
	# plot circle centers
	plt.plot(points[:,0], points[:,1], 'ko', markersize=3)
	# plot convex hull
	for simplex in hull.simplices:
		plt.plot(points[simplex, 0], points[simplex, 1], color='red', linestyle='--', lw=2, alpha=1)
	# plot the optimal tour 
	plt.plot(optTour[:,0],optTour[:,1], linestyle='-', color = 'blue', markersize=3, lw=2)

	# plt.xticks(np.arange(min(XX)-max(R),max(XX)+max(R),0.5))
	# plt.yticks(np.arange(min(YY)-max(R),max(YY)+max(R),0.5))
	plt.grid(alpha=.1)
	# plt.axis('square')
	plt.savefig('/home/latte/Dropbox/Box_Research/Github/CETSP/python/cetsp_Binst_8_10.pdf',dpi=300)
	plt.show()

def find_separators(depot, X, Y, R):
	ls_points = [depot]
	for x, y in zip(X, Y):
		ls_points.append([x,y])
	points = np.array(ls_points)
	# print(points)
	hull = ConvexHull(points)

	# hull.simplices -> facet
	# hull.vertices
	innerpoint  = []
	
	for e in range(len(X)+1):
		if e not in hull.vertices:
			innerpoint = points[e,:]
			break
	# print(innerpoint)

	if len(innerpoint) == 0: #all points are vertices of the convex hull
		xum = 0
		yum = 0
		for e in hull.vertices:
			xum += points[e,0]
			yum += points[e,1]
		innerpoint = [xum/2.0, yum/2.0]	
	seperators = []
	for simplex in hull.simplices:
		# print(simplex)
		AB = [(points[simplex[0],0],points[simplex[0],1]),(points[simplex[1],0],points[simplex[1],1])]
		x_coords, y_coords = zip(*AB)
		A = np.vstack([x_coords,np.ones(len(x_coords))]).T
		m, c = np.linalg.lstsq(A, y_coords,rcond=None)[0]
		# print("Line Solution is y = {m}x + {c}".format(m=m,c=c))
		
		if m * innerpoint[0] + c < innerpoint[1]:
			seperators.append([m, -1, c])
		else:
			seperators.append([-m, +1, -c])

	# plot_Behdani_instance(depot, X, Y, R)

	return seperators

if __name__ == "__main__":
	instance = argv[1]
	depot, Ox, Oy, Or = read_instance_Behdani(instance)
	plot_Behdani_instance(depot, Ox, Oy, Or)
	# construct_convexhull(depot, Ox, Oy, Or)
