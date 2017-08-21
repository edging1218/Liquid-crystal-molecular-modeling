import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
markers = ['<', 'o', '*', '>', 'd', 's', '^', 'p','h']
colors = ['b', 'm', 'c', 'r', 'g', 'y', 'k']
filenames = ['1', '2', '3', '4', '5', '6']
labels = ['2.4', '2.6', '2.8', '3.0', '3.2', '3.4']
label_size = 24
linewidth = 0
angle = np.linspace(0,90,19)
for i, filename in enumerate(filenames):
	lines = open('c45d' + filename + '.out').readlines()
	data = np.array([x.strip().split() for x in lines])
	data = data.astype(np.float)
	data[:, 2] *= 10.438
	plt.plot(angle, data[:,2], mew = 0, label = 'd/R = ' + labels[i], linestyle = 'dashed', lw = linewidth, marker = markers[i], color = colors[i])
	tck = interpolate.splrep(angle, data[:,2], s=10)
	en = interpolate.splev(angle, tck, der=0)
	plt.plot(angle, en, mew = 0, linestyle = 'dashed', lw = linewidth + 1, color = colors[i])
plt.legend(bbox_to_anchor=(0.6, 1), loc=2, borderaxespad = 0.)
plt.xlabel('Theta (degrees)', fontsize = label_size)
plt.ylabel('Free energy (KT)', fontsize = label_size)
plt.tick_params(labelsize = 20)
plt.grid()
#plt.ylim([-10, 250])
plt.savefig('c45dis.png', bbox_inches='tight')
plt.show()
