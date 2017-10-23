import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
name = ["polx", "poly", "polz"]; 

for filename in name:
    polx = [];
    lines = open(filename + '.out').readlines()
    for l in lines:
            if len(l.split()):
                    polx.append([float(x) for x in l.split()])
    polx = np.array(polx).T
   # im = plt.imshow(polx, cmap=plt.get_cmap('Blues'), interpolation='bicubic', vmin=0, vmax=1.5);
    im = plt.imshow(polx, cmap=plt.get_cmap('bone'), interpolation='bicubic', vmax = 1.5);
    #plt.colorbar();
    im.axes.get_xaxis().set_visible(False);
    im.axes.get_yaxis().set_visible(False);
    plt.savefig(filename + ".png");
    plt.clf();
