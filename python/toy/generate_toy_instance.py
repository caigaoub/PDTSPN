import gurobipy as gp
from gurobipy import GRB
from itertools import combinations

from sys import argv
import re
import math
import numpy as np
import random as rd
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, PathPatch
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
import mpl_toolkits.mplot3d.art3d as art3d




def write_instance(filename):
	nb_tars = 50000
	file = open(filename, "w")
	file.write(str(0) + '\t' + str(0) + '\n')
	for i in range(1,nb_tars):
		file.write(str(i*3) + '\t' + str(i*3) + '\t' + str(2) + '\n')
	file.write(str(nb_tars*3) + '\t' + str(0) + '\t' + str(2) + '\n')
	file.close()

def write_rdminstance(filename):
	nb_tars = 1000
	file = open(filename, "w")
	file.write(str(0) + '\t' + str(0) + '\n')
	for i in range(1,nb_tars):
		file.write(str(rd.randint(0,1000)) + '\t' + str(rd.randint(0,1000)) + '\t' + str(50) + '\n')
		# file.write(str(i*3) + '\t' + str(i*3) + '\t' + str(2) + '\n')
	file.write(str(nb_tars*3) + '\t' + str(0) + '\t' + str(2) + '\n')
	file.close()


if __name__ == "__main__":
	# filename = '/home/cai/Dropbox/Box_Research/Github/CETSP/python/toy/toy1.cetsp'
	# write_instance(filename)

	filename = '/home/latte/Dropbox/Box_Research/Github/CETSP/python/toy/toy1rdm.cetsp'
	write_rdminstance(filename)
