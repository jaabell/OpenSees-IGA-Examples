#  refinement testing

import numpy as np
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, knotvector, utilities
from surfVisualize import *

# These are given in v,u
controlPts = np.array([
    [0, 0, 1],
   	[0.5,1,1],
   	[1,0,1]
])

# Create a B-Spline curve instance
curve = NURBS.Curve()

# Set up the curve
curve.degree = 2
# curve.ctrlpts = # Setting control points for surface
curve.set_ctrlpts(controlPts.tolist(), 3)

# Auto-generate knot vector
curve.knotvector = utilities.generate_knot_vector(curve.degree, len(curve.ctrlpts))

# Refine surface
operations.refine_knotvector(curve, [1])

# Set evaluation delta
curve.delta = 0.01

# Evaluate curve
curve.evaluate()

# Plot the control point polygon and the evaluated curve
vis_comp = VisMPL.VisCurve2D()
curve.vis = vis_comp
curve.render()

# Good to have something here to put a breakpoint
pass