# Function to visualize nurbs surfaces using geomdl
# Felipe Elgueta
# Universidad de los Andes, Chile

from geomdl.visualization import VisMPL
from geomdl.visualization import VisVTK
from geomdl import NURBS
from multiprocessing import Process


def surfVisualize(surf,i=0, sample_size=25,hold=False,save=False):
         # Set sample size and evaluate surface
    surf.sample_size = 25
    surf.evaluate()

    # Import colormaps

    # if is_interactive():
    #     ioff()

    # Plot the surface

    # ion()
    from matplotlib import cm

    vis_config = VisMPL.VisConfig(
        ctrlpts=True, axes=True, legend=True)
    vis_comp = VisMPL.VisSurface(vis_config)

    # vis_config = VisVTK.VisConfig(
    #     ctrlpts=True, axes=True, legend=True)
    # vis_comp = VisVTK.VisSurface(vis_config)
    
    surf.vis = vis_comp

    # show(block=True)
    if hold:
        p = Process(target=surfVisualize, args=(surf,))
        p.start()
    else:
        if save:
            string="mode{}.pdf".format(i)
            print("Saving: ", string)
            surf.render(filename=string,colormap=cm.coolwarm,plot=False,label=string.strip(".pdf"))
            print("Finished Saving")

        else:
            surf.render(colormap=cm.coolwarm)


    # show(block=False)


    # ioff()

    # pause(1)
    # Good to have something here to put a breakpoint
    pass

# def surfVisualize_thread(surf):
#     from multiprocessing import Process
#     p = Process(target=surfVisualize, args=(surf,))
#     p.start()

