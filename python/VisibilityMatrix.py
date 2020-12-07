import Main
from sys import argv
import math
import numpy as np
import shapely.geometry as spg


def get_blockmatrix(depot, Ox, Oy, Or):
	nb_targets = len(Ox)
	blkmat = np.zeros((nb_targets+2, nb_targets+2))

	for i in range(nb_targets):
		circle = [Ox[i], Oy[i], Or[i]]
		if is_point_circle_blocked(depot, circle, Ox, Oy, Or, i):
			blkmat[0, i+1] = 1
			blkmat[i+1, 0] = 1
			blkmat[nb_targets+1, i+1] = 1
			blkmat[i+1, nb_targets+1] = 1

	for i in range(nb_targets):
		for j in range(i+1, nb_targets):
			cira = [Ox[i], Oy[i], Or[i]]
			cirb = [Ox[j], Oy[j], Or[j]]
			if  not is_two_circles_intersected(cira, cirb):
				if is_two_circles_blocked(cira, cirb, Ox, Oy, Or, i, j):
					blkmat[i+1, j+1] = 1
					blkmat[j+1, i+1] = 1
			# exit()
	for i in range(nb_targets+2):
		for j in range(nb_targets+2):
			print(int(blkmat[i,j]), end=' ')
		print(end='\n')
	return blkmat

def is_two_circles_blocked(cira, cirb, Ox, Oy, Or, idxa, idxb):
	if cira[2] == cirb[2]:
		maxtheta = math.pi/2.0
	else:
		dist = math.sqrt((cira[0]-cirb[0])**2 + (cira[1]-cirb[1])**2)
		ratio = cira[2]/cirb[2] if cira < cirb else cirb[2]/cira[2]
		x = dist * ratio /(1-ratio)
		maxtheta = math.acos(cira[2]/x) if cira < cirb else math.acos(cirb[2]/x)

	#
	itp_a = get_lineseg_circle_intersection( [cira[0],cira[1]],  [cirb[0],cirb[1]], cira)
	itp_a = itp_a.coords[0]
	obse_points = 5
	max_theta_a = maxtheta if cira > cirb else math.pi - maxtheta
	angleset_a = np.linspace(-max_theta_a, max_theta_a, obse_points).tolist()
	potentialvispoints_a = []
	for ang in angleset_a: 
		qx, qy = rotate([cira[0],cira[1]], itp_a,  ang)
		potentialvispoints_a.append([qx,qy])
	# print(potentialvispoints_a)

	#
	itp_on_cirb = get_lineseg_circle_intersection( [cira[0],cira[1]],  [cirb[0],cirb[1]], cirb)
	itp_on_cirb = itp_on_cirb.coords[0]
	maxtheta_cirb = math.pi - maxtheta if cira > cirb else maxtheta
	angleset_cirb = np.linspace(-maxtheta_cirb, maxtheta_cirb, obse_points).tolist()
	potentialvispoints_b = []
	for ang in angleset_cirb: 
		qx, qy = rotate([cirb[0],cirb[1]], itp_on_cirb,  ang)
		potentialvispoints_b.append([qx,qy])

	# print(potentialvispoints_b)
	seen = np.ones((obse_points, obse_points))
	for idxpa in range(obse_points):
		for idxpb in range(obse_points):
			for i in range(0, len(Ox)):
				if not (i == idxa or i == idxb):
					if is_lineseg_circle_intersected(potentialvispoints_a[idxpa], potentialvispoints_b[idxpb], [Ox[i], Oy[i], Or[i]]):
						# print(potentialvispoints_a[idxpa], potentialvispoints_b[idxpb], i+1)
						seen[idxpa, idxpb] = 0
						break
	# print(seen)
	for i in range(obse_points):
		for j in range(obse_points):	
			# print(seen[i,j], potentialvispoints_a[i], potentialvispoints_b[j])
			if seen[i,j] == 1:
				return False
	return True

def is_point_circle_blocked(p, circle, Ox, Oy, Or, idx):
	dist = math.sqrt((p[0]-circle[0])**2 + (p[1]-circle[1])**2)

	maxtheta = math.acos(circle[2]/dist)
	origin = [circle[0],circle[1]]
	itp = get_lineseg_circle_intersection(p, origin, circle)
	itp = itp.coords[0]
	obse_points = 5
	angleset = np.linspace(-maxtheta, maxtheta, obse_points).tolist()
	potentialvispoints = []
	for ang in angleset: 
		qx, qy = rotate([circle[0],circle[1]], itp,  ang)
		potentialvispoints.append([qx,qy])
	# print(potentialvispoints)	
	seen = np.ones(obse_points)
	for idxp in range(len(potentialvispoints)):
		for i in range(0, len(Ox)):
			if i != idx:
				if is_lineseg_circle_intersected(potentialvispoints[idxp], p, [Ox[i], Oy[i], Or[i]]):
					seen[idxp] = 0
					break
	# print(seen)
	for e in seen:
		if e == 1:
			return False
	return True

def rotate(origin, point, angle):
	ox, oy = origin
	px, py = point
	qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
	qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
	return qx, qy


def is_lineseg_circle_intersected(pa, pb, cir):
	p = spg.Point(cir[0],cir[1])
	c = p.buffer(cir[2]).boundary
	ls = spg.LineString([(pa[0],pa[1]), (pb[0],pb[1])])
	if c.intersects(ls):
		return True
	else:
		return False

def get_lineseg_circle_intersection(pa, pb, cir):
	p = spg.Point(cir[0],cir[1])
	c = p.buffer(cir[2]).boundary
	ls = spg.LineString([(pa[0],pa[1]), (pb[0],pb[1])])
	itp = c.intersection(ls)
	return itp

def is_two_circles_intersected(circle1, circle2):
	dist = 	math.sqrt((circle1[0]-circle2[0])**2 + (circle1[1]-circle2[1])**2)
	if dist - circle1[2] -circle2[2] <= 0:
		return True
	else:
		return False

def remove_useless_circles(depot, X, Y, R):
	X.insert(0, depot[0])
	Y.insert(0, depot[1])
	R.insert(0, 0.00001)
	nb_cirs = len(X)
	
	is_del = [False] * nb_cirs	
	for i in range(0, nb_cirs):
		for j in range(0, nb_cirs):
			if i != j:
				d = math.sqrt((X[i]-X[j])**2 + (Y[i]-Y[j])**2)
				if R[i] < R[j]:
					if d + R[i] <= R[j]:
						is_del[j] = True
				else:
					if d + R[j] <= R[i]:
						is_del[i] = True

	nb_removes = 0
	for e in is_del:
		if e == True:
			nb_removes += 1
	print("nb of removes ", nb_removes)
	return is_del

def get_bubble_1():
	block_mat = np.zeros((38, 38))
	block_mat[1,3] = 1
	block_mat[1,4] = 1
	block_mat[1,5] = 1
	block_mat[1,6] = 1
	block_mat[1,7] = 1
	block_mat[1,8] = 1
	block_mat[1,9] = 1
	block_mat[1,10] = 1

	block_mat[2,4] = 1
	block_mat[2,5] = 1
	block_mat[2,6] = 1
	block_mat[2,7] = 1
	block_mat[2,8] = 1
	block_mat[2,9] = 1
	block_mat[2,10] = 1
	
	block_mat[3,1] = 1
	block_mat[3,5] = 1
	block_mat[3,6] = 1
	block_mat[3,7] = 1
	block_mat[3,8] = 1
	block_mat[3,9] = 1
	block_mat[3,10] = 1
	block_mat[4,1] = 1
	block_mat[4,2] = 1
	block_mat[4,6] = 1
	block_mat[4,7] = 1
	block_mat[4,8] = 1
	block_mat[4,9] = 1
	block_mat[4,10] = 1
	block_mat[5,1] = 1
	block_mat[5,2] = 1
	block_mat[5,3] = 1
	block_mat[5,7] = 1
	block_mat[5,8] = 1
	block_mat[5,9] = 1
	block_mat[5,10] = 1
	block_mat[6,1] = 1
	block_mat[6,2] = 1
	block_mat[6,3] = 1
	block_mat[6,4] = 1
	block_mat[6,8] = 1
	block_mat[6,9] = 1
	block_mat[6,10] = 1
	block_mat[7,1] = 1
	block_mat[7,2] = 1
	block_mat[7,3] = 1
	block_mat[7,4] = 1
	block_mat[7,5] = 1
	block_mat[7,9] = 1
	block_mat[7,10] = 1
	block_mat[8,1] = 1
	block_mat[8,2] = 1
	block_mat[8,3] = 1
	block_mat[8,4] = 1
	block_mat[8,5] = 1
	block_mat[8,6] = 1
	block_mat[8,10] = 1
	block_mat[9,1] = 1
	block_mat[9,2] = 1
	block_mat[9,3] = 1
	block_mat[9,4] = 1
	block_mat[9,5] = 1
	block_mat[9,6] = 1
	block_mat[9,7] = 1
	block_mat[10,1] = 1
	block_mat[10,2] = 1
	block_mat[10,3] = 1
	block_mat[10,4] = 1
	block_mat[10,5] = 1
	block_mat[10,6] = 1
	block_mat[10,7] = 1
	block_mat[10,8] = 1

	block_mat[10,22] = 1
	block_mat[10,23] = 1
	block_mat[10,24] = 1
	block_mat[10,25] = 1
	block_mat[10,26] = 1
	block_mat[10,27] = 1
	block_mat[10,28] = 1
	block_mat[10,20] = 1

	block_mat[21,23] = 1
	block_mat[21,24] = 1
	block_mat[21,25] = 1
	block_mat[21,26] = 1
	block_mat[21,27] = 1
	block_mat[21,28] = 1
	block_mat[21,20] = 1

	block_mat[22,10] = 1
	block_mat[22,24] = 1
	block_mat[22,25] = 1
	block_mat[22,26] = 1
	block_mat[22,27] = 1
	block_mat[22,28] = 1
	block_mat[22,20] = 1

	block_mat[23,10] = 1
	block_mat[23,21] = 1
	block_mat[23,25] = 1
	block_mat[23,26] = 1
	block_mat[23,27] = 1
	block_mat[23,28] = 1
	block_mat[23,20] = 1

	block_mat[24,10] = 1
	block_mat[24,21] = 1
	block_mat[24,22] = 1
	block_mat[24,26] = 1
	block_mat[24,27] = 1
	block_mat[24,28] = 1
	block_mat[24,20] = 1

	block_mat[25,10] = 1
	block_mat[25,21] = 1
	block_mat[25,22] = 1
	block_mat[25,23] = 1
	block_mat[25,27] = 1
	block_mat[25,28] = 1
	block_mat[25,20] = 1

	block_mat[26,10] = 1
	block_mat[26,21] = 1
	block_mat[26,22] = 1
	block_mat[26,23] = 1
	block_mat[26,24] = 1
	block_mat[26,28] = 1
	block_mat[26,20] = 1

	block_mat[27,10] = 1
	block_mat[27,21] = 1
	block_mat[27,22] = 1
	block_mat[27,23] = 1
	block_mat[27,24] = 1
	block_mat[27,25] = 1
	block_mat[27,20] = 1

	block_mat[28,10] = 1
	block_mat[28,21] = 1
	block_mat[28,22] = 1
	block_mat[28,23] = 1
	block_mat[28,24] = 1
	block_mat[28,25] = 1
	block_mat[28,26] = 1

	block_mat[20,10] = 1
	block_mat[20,21] = 1
	block_mat[20,22] = 1
	block_mat[20,23] = 1
	block_mat[20,24] = 1
	block_mat[20,25] = 1
	block_mat[20,26] = 1
	block_mat[20,27] = 1

	block_mat[20,18] = 1
	block_mat[20,17] = 1
	block_mat[20,16] = 1
	block_mat[20,15] = 1
	block_mat[20,14] = 1
	block_mat[20,13] = 1
	block_mat[20,12] = 1
	block_mat[20,11] = 1

	block_mat[19,17] = 1
	block_mat[19,16] = 1
	block_mat[19,15] = 1
	block_mat[19,14] = 1
	block_mat[19,13] = 1
	block_mat[19,12] = 1
	block_mat[19,11] = 1

	block_mat[18,20] = 1
	block_mat[18,16] = 1
	block_mat[18,15] = 1
	block_mat[18,14] = 1
	block_mat[18,13] = 1
	block_mat[18,12] = 1
	block_mat[18,11] = 1

	block_mat[17,20] = 1
	block_mat[17,19] = 1
	block_mat[17,15] = 1
	block_mat[17,14] = 1
	block_mat[17,13] = 1
	block_mat[17,12] = 1
	block_mat[17,11] = 1

	block_mat[16,20] = 1
	block_mat[16,19] = 1
	block_mat[16,18] = 1
	block_mat[16,14] = 1
	block_mat[16,13] = 1
	block_mat[16,12] = 1
	block_mat[16,11] = 1

	block_mat[15,20] = 1
	block_mat[15,19] = 1
	block_mat[15,18] = 1
	block_mat[15,17] = 1
	block_mat[15,13] = 1
	block_mat[15,12] = 1
	block_mat[15,11] = 1

	block_mat[14,20] = 1
	block_mat[14,19] = 1
	block_mat[14,18] = 1
	block_mat[14,17] = 1
	block_mat[14,16] = 1
	block_mat[14,12] = 1
	block_mat[14,11] = 1

	block_mat[13,20] = 1
	block_mat[13,19] = 1
	block_mat[13,18] = 1
	block_mat[13,17] = 1
	block_mat[13,16] = 1
	block_mat[13,15] = 1
	block_mat[13,11] = 1

	block_mat[12,20] = 1
	block_mat[12,19] = 1
	block_mat[12,18] = 1
	block_mat[12,17] = 1
	block_mat[12,16] = 1
	block_mat[12,15] = 1
	block_mat[12,14] = 1

	block_mat[11,20] = 1
	block_mat[11,19] = 1
	block_mat[11,18] = 1
	block_mat[11,17] = 1
	block_mat[11,16] = 1
	block_mat[11,15] = 1
	block_mat[11,14] = 1
	block_mat[11,13] = 1


	block_mat[1,30] = 1
	block_mat[1,31] = 1
	block_mat[1,32] = 1
	block_mat[1,33] = 1
	block_mat[1,34] = 1
	block_mat[1,35] = 1
	block_mat[1,36] = 1
	block_mat[1,11] = 1

	block_mat[29,31] = 1
	block_mat[29,32] = 1
	block_mat[29,33] = 1
	block_mat[29,34] = 1
	block_mat[29,35] = 1
	block_mat[29,36] = 1
	block_mat[29,11] = 1


	block_mat[30,1] = 1
	block_mat[30,32] = 1
	block_mat[30,33] = 1
	block_mat[30,34] = 1
	block_mat[30,35] = 1
	block_mat[30,36] = 1
	block_mat[30,11] = 1

	block_mat[31,1] = 1
	block_mat[31,29] = 1
	block_mat[31,33] = 1
	block_mat[31,34] = 1
	block_mat[31,35] = 1
	block_mat[31,36] = 1
	block_mat[31,11] = 1

	block_mat[32,1] = 1
	block_mat[32,29] = 1
	block_mat[32,30] = 1
	block_mat[32,34] = 1
	block_mat[32,35] = 1
	block_mat[32,36] = 1
	block_mat[32,11] = 1

	block_mat[33,1] = 1
	block_mat[33,29] = 1
	block_mat[33,30] = 1
	block_mat[33,31] = 1
	block_mat[33,35] = 1
	block_mat[33,36] = 1
	block_mat[33,11] = 1

	block_mat[34,1] = 1
	block_mat[34,29] = 1
	block_mat[34,30] = 1
	block_mat[34,31] = 1
	block_mat[34,32] = 1
	block_mat[34,36] = 1
	block_mat[34,11] = 1

	block_mat[35,1] = 1
	block_mat[35,29] = 1
	block_mat[35,30] = 1
	block_mat[35,31] = 1
	block_mat[35,32] = 1
	block_mat[35,33] = 1
	block_mat[35,11] = 1

	block_mat[36,1] = 1
	block_mat[36,29] = 1
	block_mat[36,30] = 1
	block_mat[36,31] = 1
	block_mat[36,32] = 1
	block_mat[36,33] = 1
	block_mat[36,34] = 1

	block_mat[11,1] = 1
	block_mat[11,29] = 1
	block_mat[11,30] = 1
	block_mat[11,31] = 1
	block_mat[11,32] = 1
	block_mat[11,33] = 1
	block_mat[11,34] = 1
	block_mat[11,35] = 1

	# for i in range(38):
	# 	for j in range(38):
	# 		if block_mat[i,j] == 1:
	# 			print(i,j, end='\t')
	return block_mat


def main_blkmat(instance):
	# depot, Ox, Oy, Or = Main.read_Mennell_instance(instance)
	depot, Ox, Oy, Or = Main.read_instance_Behdani(instance)
	
	is_del = remove_useless_circles(depot, Ox, Oy, Or)
	new_Ox = []
	new_Oy = []
	new_Or = []
	for i, val in enumerate(is_del):
		print(i,val, end=', ')
	for i, val in enumerate(is_del[1:]):
		if  not val:
			new_Ox.append(Ox[i+1])
			new_Oy.append(Oy[i+1])
			new_Or.append(Or[i+1])
	Main.plot_instance(depot, new_Ox, new_Oy, new_Or)

	blkmat = get_blockmatrix(depot, new_Ox, new_Oy, new_Or)

	return blkmat, depot, new_Ox, new_Oy, new_Or
	# Main.plot_instance(depot, Ox, Oy, Or)
	# Main.plot_instance(depot, new_Ox, new_Oy, new_Or)


if __name__ == "__main__":
	instance = argv[1]
	''' discs '''
	# depot, Ox, Oy, Or = Main.read_instance_Behdani(instance)
	depot, Ox, Oy, Or = Main.read_Mennell_instance(instance)
	is_del = remove_useless_circles(depot, Ox, Oy, Or)
	new_Ox = []
	new_Oy = []
	new_Or = []
	# Main.plot_instance(depot, Ox, Oy, Or)
	for i, val in enumerate(is_del):
		print(i,val, end=', ')
	for i, val in enumerate(is_del[1:]):
		if  not val:
			new_Ox.append(Ox[i+1])
			new_Oy.append(Oy[i+1])
			new_Or.append(Or[i+1])

	get_blockmatrix(depot, new_Ox, new_Oy, new_Or)

	# Main.plot_instance(depot, Ox, Oy, Or)
	Main.plot_instance(depot, new_Ox, new_Oy, new_Or)
	