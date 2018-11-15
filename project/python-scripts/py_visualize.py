# -*- coding: utf-8 -*-

from ims_core import ims, peak_list
import sys
import matplotlib.patches as pat
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import itertools

## main script
if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print("usage: python3 py_visualize.py input-measurement.csv [input-peak-list.csv [input-peak-list.csv]]")
        exit()
    m = ims(sys.argv[1])
    
    peak_lists = []
    for i in range(2, len(sys.argv)):
        read_peak_list = peak_list(sys.argv[i])
        pl = []
        pl.append([ipl.t for ipl in read_peak_list.ims_peak_list])
        pl.append([ipl.r for ipl in read_peak_list.ims_peak_list])
        peak_lists.append(pl)
    
    cdict = {'red':   ((0.0, 1.0, 1.0),
                (0.05, 0.0, 0.0),
                (0.1, 1.0, 1.0),
                (0.2, 1.0, 1.0),
                (0.6, 1.0, 1.0),
                (1.0, 1.0, 1.0)),

        'green': ((0.0, 1.0, 1.0),
                (0.05,0.0, 0.0),
                (0.1, 0.0, 0.0),
                (0.2, 0.0, 0.0),
                (0.6, 1.0, 1.0),
                (1.0, 1.0, 1.0)),

        'blue':  ((0.0, 1.0, 1.0),
                (0.05, 1.0, 1.0),
                (0.1, 1.0, 1.0),
                (0.2, 0.0, 0.0),
                (0.6, 0.0, 0.0),
                (1.0, 0.0, 0.0))}

    cmap = mplc.LinearSegmentedColormap('ims_colormap',cdict,256)
    axis = plt.gca().set_color_cycle(['black', 'DarkGray', 'green', 'brown'])
    marker = ['.', '*', 'x', '+']
    plt.imshow(m.points, interpolation="nearest", origin="lower",vmin=0, vmax=100, cmap=cmap, extent = m.extent, aspect="auto")
    for i, pl in enumerate(peak_lists):
        plt.plot(pl[0], pl[1], linestyle = '', marker = marker[i % 4], markersize=10.0, label = str(i + 1) + ". peak list")
        
    plt.subplots_adjust(left=0.02, bottom=0.02, right=1, top=1, wspace=None, hspace=None)
    plt.legend()
    plt.show(block=True)