import numpy as np

K_ops=np.load("K_ops.npy")
K_python=np.load("K_python.npy")

diff=K_ops-K_python

import matplotlib.pyplot as plt

fig=plt.figure()

plt.subplot(311)
plt.imshow(K_ops)
plt.title("K_ops")
plt.colorbar()

plt.subplot(312)
plt.imshow(K_python)
plt.title("K_python")
plt.colorbar()

plt.subplot(313)
plt.imshow(diff)
plt.title("diff")
plt.colorbar()

plt.show()


