import gurobipy as gp
from gurobipy import GRB
from itertools import combinations

from sys import argv
import re
import math
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, PathPatch
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
import mpl_toolkits.mplot3d.art3d as art3d

def plot_instance(X, Y, R, Tour):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	# plt.axis('off')
	ax.set_xlim([min(X)-max(R), max(X)+max(R)])
	ax.set_ylim([min(Y)-max(R), max(Y)+max(R)])
	ax.set_aspect('equal', adjustable='box')
	if True:
		for i in range(len(X)):
			if True:
				circle = plt.Circle((X[i], Y[i]), radius=R[i], edgecolor='r',fill=False, alpha = 0.5)
				ax.add_artist(circle)
		plt.plot(Tour[0], Tour[1])
	plt.show()

def test_III():
	X = [0, 2, 1]
	Y = [0, 0, 2]
	R = [1.2, 1.2, 1.2]

	pset1 = []
	pset2 = []
	pset3 = []
	nbp = 100
	angle = 2.0*math.pi/float(nbp)
	for i in range(3):
		for j in range(nbp):
			x = R[i] * math.cos(j*angle) + X[i]
			y = R[i] * math.sin(j*angle) + Y[i]
			if i == 0:
				pset1.append([x,y])
			if i == 1:
				pset2.append([x,y])
			if i == 2:
				pset3.append([x,y])
	# print(pset1)
	# print(pset2)
	# print(pset3)

	SD = 10000
	bestp1 = []
	bestp2 = []
	bestp3 = []
	for p1 in pset1:
		for p2 in pset2:
			for p3 in pset3:
				d12 = eucl_dist(p1, p2)
				d23 = eucl_dist(p2, p3)
				d31 = eucl_dist(p3, p1)
				if d12 + d23 + d31 < SD:
					SD = d12 + d23 + d31
					bestp1 = p1
					bestp2 = p2
					bestp3 = p3

	print(bestp1)
	print(bestp2)
	print(bestp3)
	tourx = [bestp1[0],bestp2[0],bestp3[0],bestp1[0]]
	toury = [bestp1[1],bestp2[1],bestp3[1],bestp1[1]]
	Tour = [tourx, toury]
	plot_instance(X, Y, R, Tour)


def test_IV():
	X = [0, 2, 0.5, 3]
	Y = [0, 0, 1, 3]
	R = [0.7, 0.7, 1.2, 1]

	pset1 = []
	pset2 = []
	pset3 = []
	pset4 = []

	nbp = 20
	angle = 2.0*math.pi/float(nbp)
	for i in range(4):
		for j in range(nbp):
			x = R[i] * math.cos(j*angle) + X[i]
			y = R[i] * math.sin(j*angle) + Y[i]
			if i == 0:
				pset1.append([x,y])
			if i == 1:
				pset2.append([x,y])
			if i == 2:
				pset3.append([x,y])
			if i == 3:
				pset4.append([x,y])
	# print(pset1)
	# print(pset2)
	# print(pset3)

	SD = 10000
	bestp1 = []
	bestp2 = []
	bestp3 = []
	bestp4 = []

	for p1 in pset1:
		for p2 in pset2:
			for p3 in pset3:
				for p4 in pset4:
					d12 = eucl_dist(p1, p2)
					d23 = eucl_dist(p2, p3)
					d34 = eucl_dist(p3, p4)
					d41 = eucl_dist(p4, p1)
					if d12 + d23 + d34 + d41 < SD:
						SD = d12 + d23 + d34 + d41
						bestp1 = p1
						bestp2 = p2
						bestp3 = p3
						bestp4 = p4

	# print(bestp1)
	# print(bestp2)
	# print(bestp3)
	# print(bestp4)

	tourx = [bestp1[0],bestp2[0],bestp3[0],bestp4[0],bestp1[0]]
	toury = [bestp1[1],bestp2[1],bestp3[1],bestp4[1],bestp1[1]]
	Tour = [tourx, toury]
	plot_instance(X, Y, R, Tour)


def eucl_dist(p1, p2):
	return math.sqrt((p1[0]-p2[0])**2 +(p1[1]-p2[1])**2)





if __name__ == "__main__":
	test_III()
	# test_IV()