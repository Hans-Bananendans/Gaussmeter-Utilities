"""
This is a script that hardcodes the results of the Weekday4 data and plots the
vectorplot. This can be used to more easily visualize the net disturbance that
was measured in the lab, since it is a 3D view rather than an image.
"""

# Imports ====================================================================
import numpy as np
from tools.plotting import VectorPlot

# Local EMF vector in the cage frame (C)
local_emf = np.array([8.1609506, 17.36678791, -45.45745])
local_emf_abs = round(np.linalg.norm(local_emf), 1)
print("Local EMF:".ljust(20, " "), np.round(local_emf, 3),
      "uT   (abs: {} uT)".format(local_emf_abs))

# Mean of the raw measured data vector, but rotated to the cage frame (C)
mdv = np.array([ -8.19666212, 9.02900061, -45.70576045])
mdv_abs = round(np.linalg.norm(mdv), 1)
print("Raw measured:".ljust(20, " "), np.round(mdv, 3),
      "uT   (abs: {} uT)".format(mdv_abs))

# Mean of the data vector, normalized to the local EMF, rotated to C.
# This is referred to as the (net) disturbance vector
mnv = np.array([-16.35761272, -8.3377873, -0.248310450])
mnv_abs = round(np.linalg.norm(mnv), 3)
print("Disturbance vector:".ljust(20, " "), np.round(mnv, 3),
      "uT   (abs: {} uT)".format(mnv_abs))

# Vector plot

plot_title = """
    Magnetic disturbance vector at the test site during a regular day
    Measured on 13-07-2023, part of a 24+ hour measurement"""

vp = VectorPlot()
vp.ax.view_init(elev=20, azim=-65.5)

# Vector scaling factor: 1 uT corresponds to 2 cm on the 3D plot
vscf = 0.02
vp.autoplot(tripod=True, coils=True, walls=True,
                    table=True)

# Local EMF vector
vp.plot_vector(local_emf, linewidth=3, scaling=vscf, alr=0.2, color="orange")

# Raw data vector
vp.plot_vector(mdv, scaling=vscf, linewidth=3, alr=0.2, color="purple")

# Disturbance vector
vp.plot_vector(mnv, color="black", scaling=vscf, linewidth=3, alr=0.3)


# Transparent line from mdv to emf vector
vp.plot_vector((local_emf-mdv)*vscf, origin=mdv*vscf, alpha=0.5, linewidth=1,
               scaling=1, alr=0.001, color="#555")
# Transparent line from mdv to disturbance vector
vp.plot_vector((mnv-mdv)*vscf, origin=mdv*vscf, alpha=0.5, linewidth=1,
               scaling=1, alr=0.001, color="#555")

# Title
vp.ax.text2D(0.15, 0.9, plot_title,
             fontsize=14,
             multialignment="center",
             transform=vp.ax.transAxes)

# Vector labels
vp.ax.text(local_emf[0]*vscf/2, local_emf[1]*vscf/2, local_emf[2]*vscf/2,
           "Local EMF\n{} uT".format(local_emf_abs),
           color="orange", weight="bold")
vp.ax.text(mdv[0]*vscf*1.2, mdv[1]*vscf*1.2, mdv[2]*vscf*1.2,
           "Raw data\n{} uT".format(mdv_abs),
           color="purple", weight="bold")
vp.ax.text(mnv[0]*vscf*1.1, mnv[1]*vscf*1.1, mnv[2]*vscf*1.1,
           "Disturbance\n{} uT".format(mnv_abs),
           color="black", ha="right", weight="bold")

# Wall labels
vp.ax.text(-1.5, 0, 0, "Hallway wall", color="#555", ha="center")
vp.ax.text(0, 1.5, 0, "Back wall", color="#555", ha="center")

vp.show()

vp.fig.savefig("./figures/disturbance_visualization.png", dpi=150)
