# Function to visualize nurbs surfaces using geomdl
# Felipe Elgueta
# Universidad de los Andes, Chile

from geomdl.visualization import VisMPL, VisVTK
from geomdl import NURBS
from multiprocessing import Process
from matplotlib import cm


def surfVisualize(surf, visualizer="VTK", sample_size=25):

    if surf.type=="container":
        for surface in surf:
            surface.sample_size=sample_size
            # surface.delta=delta
            surface.evaluate()
    else:
        surf.sample_size=sample_size
        # surf.sample_size_v=100
        # surf.delta=delta
        surf.evaluate()

    # Import colormaps
    if visualizer == 'MPL':
        vis_config = VisMPL.VisConfig(
            ctrlpts=False, axes=True, legend=False)
        vis_comp = VisMPL.VisSurface(vis_config)

    elif visualizer == "VTK":
        vis_config = VisVTK.VisConfig(
            ctrlpts=False, axes=True, legend=False)
        vis_comp = VisVTK.VisSurface(vis_config)

    surf.vis = vis_comp
    surf.render(colormap=cm.coolwarm)
