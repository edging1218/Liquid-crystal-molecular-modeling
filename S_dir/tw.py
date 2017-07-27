import matplotlib.pyplot as plt
import numpy as np
LW = 1.5
MS = 4
FS = 25


filename = 'TwAn.out'
xxs = []
pos = []
S = []
tw = []

lines = open(filename).readlines()
for l in lines:
        if len(l.split()):
                xxs.append([float(x) for x in l.split()])
for x in xxs:
	if x[1] != 0:
		pos.append(x[0])
		S.append(x[1])
		tw.append(abs(x[2]))

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(pos, S, color='m', linestyle='dashed', marker='^', markersize=MS, lw = LW, mew= 0, label='Order parameter')
ax.legend(loc='upper left', fontsize=17)
plt.xticks(fontsize=FS)
plt.yticks(fontsize=FS)
plt.ylim([0.766, 0.772])

ax2 = ax.twinx()
ax2.plot(pos, tw, color='b', linestyle='dashed', marker='d', markersize=MS, lw = LW, mew= 0, label='Tw')
ax2.legend(loc='upper right', fontsize=17)
plt.ylim([-0.1, 1.3])
ax.set_xlabel('Anchoring Strength $W$ ($J/m^2$)', fontsize=25)
ax.set_ylabel('$F$ - $F_{min}$ ($10^3$ $k_BT$)', fontsize=25)
ax2.set_ylabel('$F_{surf}$ ($10^3$ $k_BT$)', fontsize=25)
#plt.gca().set_xscale('log')
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

plt.savefig('Tw', bbox_inches='tight');
plt.show();
