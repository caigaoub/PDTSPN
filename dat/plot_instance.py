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

def read_instance(instancefile):
	file_ = open(instancefile, 'r')
	X = []
	Y = []
	R = []
	while True:
		line_ = file_.readline()
		list_ = re.split(" |\t|\n", line_)

		if list_[0] == '':
			break
		else:
			list_ = list(filter(None, list_))
			print(list_)
			X.append(float(list_[0]))
			Y.append(float(list_[1]))
			R.append(float(list_[3]))

	file_.close()
	return X, Y, R

def plot_instance(X, Y, R):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	# plt.axis('off')
	ax.set_xlim([min(X)-max(R), max(X)+max(R)])
	ax.set_ylim([min(Y)-max(R), max(Y)+max(R)])
	ax.set_aspect('equal', adjustable='box')
	if True:
		for i in range(len(X)):
			circle = plt.Circle((X[i], Y[i]), radius=R[i], edgecolor='r',fill=False, alpha = 0.5)
			ax.add_artist(circle)
	plt.show()


if __name__ == "__main__":
	instance = argv[1]
	X, Y, R = read_instance(instance)
	plot_instance(X, Y, R)
