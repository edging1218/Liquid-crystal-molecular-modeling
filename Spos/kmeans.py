from sklearn.cluster import KMeans
import numpy as np


filename = 'Spos.out'
with open(filename) as f:
	content = f.readlines()
content = np.array([x.strip().split('\t') for x in content])


n_cluster = input('n_clusters:')
R = input('R:')
rc = np.array([R + 2, R + 2, R + 2])
model = KMeans(n_clusters = n_cluster)
model.fit(content)
centers= model.cluster_centers_

for i in range(n_cluster):
	r = centers[i] - rc
	centers[i] = r / np.sqrt(r.dot(r)) * R
for i in range(n_cluster):
	for j in range(i + 1, n_cluster):
		dis = centers[i] - centers[j]
		dis = np.sqrt(dis.dot(dis))
#		print dis
		angle = 2 * np.arcsin(dis * 0.5 / R) / np.pi * 180
		print i, j, angle




