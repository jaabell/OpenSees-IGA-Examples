import numpy as np


uCuts = np.array([0.011, 0.251, 0.696, 1.141, 1.586])  # Places to define sections


from aspa_geometry import *
aspa_geometry(uCuts)

# Flag to decide whether to use bending strips or not (for debugging)
use_bendingStrip = True

import aspa_loaded
from aspa_loaded import *

# Angle of each ply for each zone
θ_0 = deg2rad*np.array([45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0])
θ_1 = deg2rad*np.array([45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 45.0, -45.0, 0.0, 90.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0])
θ_2 = deg2rad*np.array([45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 45.0, -45, -45.0, 45.0, -45.0, 45.0])

θs = [θ_0, θ_1, θ_2]


aspa_loaded(use_bendingStrip, θs)