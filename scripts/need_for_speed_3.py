from aimsprop import need_for_speed3 as nfs
import numpy as np
import time

nfs.first_try()

def take_this_geom_p(x):
    a = x.shape[0]
    b = x.shape[1]
    for i in range(a):
        for j in range(b):
            x[i,j] = multiply_plus_ten_p(x[i,j])
    return x

def multiply_plus_ten_p(x):
    return x*2+10

a = np.random.rand(3000000,3)
t0 = time.time()
pyth = take_this_geom_p(a)
t1 = time.time()
num = a*2+10
t2 = time.time()
b = nfs.take_this_geom(a)
t3 = time.time()

print(np.all(b==num) and np.all(pyth==num), t0-t1, t2-t1, t3-t2)
