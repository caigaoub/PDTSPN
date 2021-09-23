


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


