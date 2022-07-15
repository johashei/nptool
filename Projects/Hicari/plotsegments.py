import numpy as np
import matplotlib.pyplot as plt

txt = np.loadtxt('segmentpositions.dat')
x = txt[:,3]
y = txt[:,4]
z = txt[:,5]

#fig, ax = plt.subplots(2,2)
#ax[0,0].scatter(x,y, marker='o')
#ax[0,0].set_xlabel('X (mm)')
#ax[0,0].set_ylabel('Y (mm)')
#ax[0,1].scatter(z,x, color='red')
#ax[0,1].set_xlabel('Z(mm)')
#ax[0,1].set_ylabel('X(mm)')
#fig.tight_layout()

fig2 = plt.figure()
ax2 = fig2.add_subplot(projection='3d')
ax2.scatter(x, z, y, marker='o')
ax2.set_xlabel('X(mm)')
ax2.set_ylabel('Z(mm)')
ax2.set_zlabel('Y(mm)')


plt.show()
