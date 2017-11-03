import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
start = input('File number start:')
increment = input('File number increment:')
end = input('File number end:')
for num in range(start, end, increment):
    sx = 'polx%.5d' % num
    sy = 'poly%.5d' % num
    sz = 'polz%.5d' % num
    name = [sx, sy, sz]; 

    polx = []

    for filename in name:
        polx = [];
        lines = open(filename + '.out').readlines()
        for l in lines:
                if len(l.split()):
                        polx.append([float(x) for x in l.split()])
    	polx = np.array(polx).T
        im = plt.imshow(polx, cmap=plt.get_cmap('bone'), interpolation='bicubic', vmin=0, vmax=1.5);
       # im = plt.imshow(polx, cmap=plt.get_cmap('bone'), interpolation='bicubic', vmax=1.5);
        #plt.colorbar();
        im.axes.get_xaxis().set_visible(False);
        im.axes.get_yaxis().set_visible(False);
        plt.savefig(filename + ".png");
        plt.clf();
