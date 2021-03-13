import numpy as np

p = np.random.rand() * 2.0 * np.pi
z = np.random.rand() * np.pi
o = np.random.rand() * 2.0 * np.pi

R1 = np.array(
    [
        [np.cos(p), -np.sin(p), 0.0],
        [np.sin(p), np.cos(p), 0.0],
        [0.0, 0.0, 1.0],
    ]
)
R2 = np.array(
    [
        [np.cos(z), 0.0, -np.sin(z)],
        [0.0, 1.0, 0.0],
        [np.sin(z), 0.0, np.cos(z)],
    ]
)
R3 = np.array(
    [
        [np.cos(o), -np.sin(o), 0.0],
        [np.sin(o), np.cos(o), 0.0],
        [0.0, 0.0, 1.0],
    ]
)

R = np.dot(R1, np.dot(R2, R3))

S = np.array(
    [
        [
            np.cos(p) * np.cos(z) * np.cos(o) - np.sin(p) * np.sin(o),
            -np.cos(p) * np.cos(z) * np.sin(o) - np.sin(p) * np.cos(o),
            -np.cos(p) * np.sin(z),
        ],
        [
            np.sin(p) * np.cos(z) * np.cos(o) + np.cos(p) * np.sin(o),
            -np.sin(p) * np.cos(z) * np.sin(o) + np.cos(p) * np.cos(o),
            -np.sin(p) * np.sin(z),
        ],
        [np.sin(z) * np.cos(o), -np.sin(z) * np.sin(o), np.cos(z)],
    ]
)

print((R - S))
print((np.max(np.abs(R - S))))
