import subprocess
import numpy as np
import pandas as pd
from numpy import linalg as la
from sklearn.cluster import KMeans
from time import time
import sys

def cal_dis(a, b):
    try:
        dis = a - b
        return np.sqrt(dis.dot(dis))
    except:
        print('Wrong input for distance calculation.')
        

time_start = time()

# cmd = '/scratch/midway2/yezhou/LC/Spos/iso'
cmd = '~/midway2/LC/Spos/iso'
filename = 'Spos.out'
start = int(sys.argv[1])
end = int(sys.argv[3])
step = int(sys.argv[2])
n_cluster = int(sys.argv[4])
R = int(sys.argv[5])
print(n_cluster, R)
n_angle = n_cluster * (n_cluster - 1) / 2
center_df = pd.DataFrame(columns = ['d' + str(i) + str(j) for i in range(n_cluster) for j in range(3)])
angle_df = pd.DataFrame(columns = ['a' + str(i) for i in range(n_angle)])

for file_idx in range(start, end, step):
	qtensor_name = 'Q' + str(file_idx).zfill(5) + '.bin'
        subprocess.call(['mv', qtensor_name, 'Qtensor.bin'])
        subprocess.call(cmd, shell=True)
        subprocess.call(['mv', 'Qtensor.bin', qtensor_name])

        with open(filename) as f:
                content = f.readlines()
        content = np.array([x.strip().split('\t') for x in content])

        rc = np.array([R + 2, R + 2, R + 2])
        if file_idx > start:
            model = KMeans(n_clusters = n_cluster, init=pre_centers)
        else:
            model = KMeans(n_clusters = n_cluster)
        model.fit(content)
        centers= model.cluster_centers_
        pre_centers=centers
        
	for i in range(n_cluster):
		r = centers[i] - rc
		centers[i] = r / np.sqrt(r.dot(r)) * R

	center_df.loc[file_idx] = centers.ravel()
	angles = []
	for i in range(n_cluster):
		for j in range(i + 1, n_cluster):
                        dis = cal_dis(centers[i], centers[j])
			angle = 2 * np.arcsin(dis * 0.5 / R) / np.pi * 180
			angles.append(angle)
	angle_df.loc[file_idx] = angles
        pre_centers = centers
        
center_df.to_csv('centers.csv')
angle_df.to_csv('angles.csv')

time_end = time()
print('Time: {} min'.format((time_end-time_start)/60))
