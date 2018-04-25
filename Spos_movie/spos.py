import numpy as np
import pandas as pd
from numpy import linalg as la
from sklearn.cluster import KMeans
from time import time
import sys

time_start = time()

start = int(sys.argv[1])
end = int(sys.argv[3])
step = int(sys.argv[2])
# n_cluster = input('n_clusters:')
n_cluster = int(sys.argv[4])
# R = input('R:')
R = int(sys.argv[5])
print(n_cluster, R)

center_df = pd.DataFrame(columns = ['d' + str(i) + str(j) for i in range(n_cluster) for j in range(3)])
angle_df = pd.DataFrame(columns = ['a' + str(i) for i in range(n_cluster * (n_cluster - 1) / 2)])

with open('param.in', 'rb') as f:
	for line in f.readlines():
		if line[:2] == 'Nx':
			nx = int(line[3:])	
		elif line[:2] == 'Ny':
			ny = int(line[3:])	
		elif line[:2] == 'Nz':
			nz = int(line[3:])	
		else:
			break
print('nx: {}\nny: {}\nnz: {}'.format(nx, ny, nz))
 
with open('grid.bin', 'rb') as f:
	grid = np.fromfile(f, dtype=np.int32)
print('Grid size: {}'.format(grid.shape[0]))

for file_idx in range(start, end, step):
	filename = 'Q' + str(file_idx).zfill(5) + '.bin'
	with open('Q00000.bin', 'rb') as f:
		qtensor = np.fromfile(f, dtype=float)
	q_idx = 0
	g_idx = 0
	defects = np.array([])
	for k in range(nz):
		for j in range(ny):
			for i in range(nx):
				if grid[g_idx] == 0 or grid[g_idx] == 1:
					t = np.zeros((3,3))
					t[0, 0] = qtensor[q_idx]
					t[0, 1] = qtensor[q_idx+1]
					t[1, 0] = qtensor[q_idx+1]
					t[0, 2] = qtensor[q_idx+2]
					t[2, 0] = qtensor[q_idx+2]
					t[1, 1] = qtensor[q_idx+3]
					t[1, 2] = qtensor[q_idx+4]
					t[2, 1] = qtensor[q_idx+4]
					t[2, 2] = qtensor[q_idx+5]
					w, v = la.eig(t)
					w = np.sort(w)
					s = 0.5*3*w[2]
					if s < 0.4:
						defects = np.append(defects, [i, j, k])
					q_idx += 6
				g_idx += 1
			
	n_defect = defects.shape[0] / 3
	defects = defects.reshape((n_defect, 3))

	# defects=np.array([[280,280,280], [0,0,0], [298, 278, 270], [3,3,3], [0.7, 0.9, 1.1]])

	rc = np.array([R + 2, R + 2, R + 2])
	model = KMeans(n_clusters = n_cluster)
	model.fit(defects)

	centers= model.cluster_centers_
	for i in range(n_cluster):
		r = centers[i] - rc
		centers[i] = r / np.sqrt(r.dot(r)) * R
	center_df.loc[file_idx] = centers.ravel()

	angles = []
	for i in range(n_cluster):
		for j in range(i + 1, n_cluster):
			dis = centers[i] - centers[j]
			dis = np.sqrt(dis.dot(dis))
			angle = 2 * np.arcsin(dis * 0.5 / R) / np.pi * 180
			angles.append(angle)
	angle_df.loc[file_idx] = angles

center_df.to_csv('centers.csv')
angle_df.to_csv('angles.csv')

time_end = time()
print('Time: {} min'.format((time_end-time_start)/60))
