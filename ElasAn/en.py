import numpy as np
import matplotlib.pyplot as plt
filename = 'elas.out'
with open(filename) as f:
	content = f.readlines()
content = np.array([x.strip().split('\t') for x in content])
plt.plot(content[:, -3])
plt.plot(content[:, -2])
plt.plot(content[:, -1])
plt.show()
