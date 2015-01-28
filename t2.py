import numpy as np
import matplotlib.pyplot as plt


xbins = 11
ybins = 11


grid = np.zeros(shape=(xbins,ybins))

grid[5,3] += 1

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
circ_1 = plt.Circle((.5, 0), radius = .05, color = 'b')
circ_2 = plt.Circle((-.5, 0), radius = .05, color = 'g')
circ_mid = plt.Circle((0, 0), radius = .05, color = 'r')

x = np.linspace(-3, 3, xbins)
y = np.linspace(-3, 3, ybins)
X, Y = np.meshgrid(x, y)
plt.contour(x, y, grid)
#plt.scatter(xs, ys)
plt.xlabel('[R_sep]')
plt.ylabel('[R_sep]')
ax.add_artist(circ_1)
ax.add_artist(circ_2)
ax.add_artist(circ_mid)
plt.show()
