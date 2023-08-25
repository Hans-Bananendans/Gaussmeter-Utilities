# Imports ====================================================================
import numpy as np
from time import time

from tools.local_emf import local_emf
from tools.plotting import VectorPlot
from tools.GaussmeterLib import (
    EulerRotation,
    read_data,
    data_normal_distribution,
    data_savgol_filter,
    data_rotate,
    data_add_vector,
)

time0 = time()

# Modelling

t = np.linspace(15, 59, 12)
t2 = np.linspace(0, 60, 601)  # For plotting dotted line only

def bxt(t_array, bmax=1, period=60):
    bx = np.zeros(len(t_array))
    for i in range(len(t_array)):
        bx[i] = bmax * np.sin(2*np.pi / period * t_array[i])
    return bx


def bzt(t_array, bmax=1, period=60):
    bz = np.zeros(len(t_array))
    for i in range(len(t_array)):
        bz[i] = bmax * np.cos(2*np.pi / period * t_array[i])
    return bz


bx = bxt(t, bmax=-241, period=60)
bz = bzt(t, bmax=241, period=60)

sx = bxt(t2, bmax=-241, period=60)
sz = bzt(t2, bmax=241, period=60)
sa = np.zeros(len(t2))
for i in range(len(t2)):
    sa[i] = i/len(t2)*1

print(max(sa), min(sa))

# Vector plot
plot_title = """
    Title
    """
vectorplot = VectorPlot()
vectorplot.autoplot(tripod=True, coils=True, walls=False,
                    table=False)
# vectorplot.plot_vector(local_emf, linewidth=3, scaling=0.02, alr=0.2, color="orange")

# B Arrows
for i in range(len(t)):
    vectorplot.plot_vector([bx[i], 0, bz[i]], color="black",
                           scaling=0.004, linewidth=1+1.5*i/len(t), alr=0.15,
                           alpha=i/len(t)*1)  # alpha=i/max(t)*1

# Circle shape
scale = 0.004
vectorplot.ax.scatter(scale*sx, 0, scale*sz, alpha=sa, s=0.5, c='#555')

# Red sinusoid
linear = np.linspace(250, 550, len(t2))
vectorplot.ax.scatter(scale*sx, 0, -scale*linear, alpha=0.8, s=0.5, c='#F00')

# Blue sinusoid
vectorplot.ax.scatter(-scale*linear, 0, scale*sz, alpha=0.8, s=0.5, c='#00F')

# Labels
vectorplot.ax.text(scale*(bx[-1]+10), 0, scale*(bz[-1]+10),
           "B(t)", color="k", size=16, weight="bold")

vectorplot.ax.text(0, 0, scale*-430,
           "I_X (t)", color="#F00", weight="bold")

vectorplot.ax.text(scale*-540, 0, 0,
           "I_Z (t)", color="#00F", weight="bold")

# Title
vectorplot.ax.text2D(0.05, 0.9, plot_title,
                     fontsize=14,
                     multialignment="center",
                     transform=vectorplot.ax.transAxes)
vectorplot.ax.view_init(elev=0, azim=-90)
vectorplot.show()
# vectorplot.fig.savefig("./test.png", dpi=150)

print("Elapsed time:", time()-time0)
