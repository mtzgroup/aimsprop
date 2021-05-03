cimport cython

def first_try():
    print("hello to your empty world")

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cpdef double [:,:] take_this_geom(double [:,:] x) nogil:
    cdef int a,b,i,j
    a = x.shape[0]
    b = x.shape[1]
    for i in range(a):
        for j in range(b):
           x[i,j] = multiply_plus_ten(x[i,j])
    return x

cdef double multiply_plus_ten(double a) nogil:
    return a*2+10
