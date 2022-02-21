import glob
from PIL import Image

# filepaths
fp_in = "figuresLinear8HZ/8HZ*.png"
fp_out = "aspa_8Hz.gif"

import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

# alist=[
#     "something1",
#     "something12",
#     "something17",
#     "something2",
#     "something25",
#     "something29"]

# alist.sort(key=natural_keys)
# print(alist)
listaImagenes=sorted(glob.glob(fp_in))
# https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif
print("sorted(glob.glob(fp_in)): ", listaImagenes)
listaImagenes.sort(key=natural_keys)
print("listaImagenes: ", listaImagenes)
img, *imgs = [Image.open(f) for f in listaImagenes]
img.save(fp=fp_out, format='GIF', append_images=imgs,
         save_all=True, duration=2e1, loop=0)