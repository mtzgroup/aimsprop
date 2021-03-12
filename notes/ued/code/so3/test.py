from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import lebedev  # E.g., the one in aimsprop/iam/lebedev.py

nlebedev = 302
nomega = 12

leb = lebedev.Lebedev.build(nlebedev)

phis = leb.phi
zetas = leb.theta
omegas = np.linspace(0.0, 2.0 * np.pi, nomega, endpoint=False)

xyz = np.array([0.0, 0.0, 1.0])
uvw = np.array([1.1, 0.0, 1.0])

xyzs = []
uvws = []
ws = []
for phi, zeta in zip(phis, zetas):
    for omega in omegas:
        R1 = np.array(
            [
                [np.cos(phi), -np.sin(phi), 0.0],
                [np.sin(phi), np.cos(phi), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        R2 = np.array(
            [
                [np.cos(zeta), 0.0, -np.sin(zeta)],
                [0.0, 1.0, 0.0],
                [np.sin(zeta), 0.0, np.cos(zeta)],
            ]
        )
        R3 = np.array(
            [
                [np.cos(omega), -np.sin(omega), 0.0],
                [np.sin(omega), np.cos(omega), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        R = np.dot(R1, np.dot(R2, R3))
        xyz1 = np.dot(R, xyz)
        uvw1 = np.dot(R, uvw)
        xyzs.append(xyz1)
        uvws.append(uvw1)
        ws.append(np.sin(zeta))

xyzs = np.array(xyzs)
uvws = np.array(uvws)
ws = np.array(ws)

uvws -= xyzs

x = xyzs[:, 0]
y = xyzs[:, 1]
z = xyzs[:, 2]
u = uvws[:, 0]
v = uvws[:, 1]
w = uvws[:, 2]

fig = plt.figure()
ax = fig.gca(projection="3d")

h = ax.quiver(x, y, z, u, v, w, length=0.1)
print((type(h)))

plt.axis("off")
plt.show()
