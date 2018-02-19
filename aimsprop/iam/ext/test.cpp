#include <boost/python.hpp>
#include <lightspeed/tensor.hpp>

using namespace lightspeed;
using namespace boost::python;

/**
 * Compute the UED 
 *  
 *
 * 
 **/
lightspeed::shared_ptr<Tensor> compute_diffraction(
    const lightspeed::shared_ptr<Tensor>& ss,
    const lightspeed::shared_ptr<Tensor>& xyzs,
    const lightspeed::shared_ptr<Tensor>& Zs,
    const lightspeed::shared_ptr<Tensor>& as,
    const lightspeed::shared_ptr<Tensor>& bs,
    const lightspeed::shared_ptr<Tensor>& Rs,
    const lightspeed::shared_ptr<Tensor>& ws,
    bool ued,
    bool mod,
    bool cos2
    )
{
    // Validity errors
    ss->ndim_error(2);
    ss->shape_error({ss->shape()[0], 3});
    xyzs->ndim_error(2);
    xyzs->shape_error({xyzs->shape()[0], 3});
    Zs->shape_error({xyzs->shape()[0]});
    Rs->ndim_error(3);
    Rs->shape_error({Rs->shape()[0], 3, 3});
    ws->shape_error({Rs->shape()[0]});

    // Target
    lightspeed::shared_ptr<Tensor> Is(new Tensor({ss->shape()[0]}));
    double* Ip = Is->data().data();

    // Pointers
    const double* sp = ss->data().data();
    const double* xyzp = xyzs->data().data();
    const double* Zp = Zs->data().data();
    const double* ap = as->data().data();
    const double* bp = bs->data().data();
    const double* Rp = Rs->data().data();
    const double* wp = ws->data().data();

    lightspeed::shared_ptr<Tensor> xyzs2(new Tensor(xyzs->shape()));
    double* xyz2p = xyzs2->data().data();
    
    for (size_t Rind = 0; Rind < Rs->shape()[0]; Rind++) {
        // Rotation quadrature
        const double* R2p = Rp + Rind * 9;
        double w = wp[Rind];
        // Rotate coordinates to frame
        for (size_t A = 0; A < xyzs->shape()[0]; A++) {
            xyz2p[3*A + 0] = R2p[0] * xyzp[3*A + 0] + R2p[1] * xyzp[3*A + 1] + R2p[2] * xyzp[3*A + 2];
            xyz2p[3*A + 1] = R2p[3] * xyzp[3*A + 0] + R2p[4] * xyzp[3*A + 1] + R2p[5] * xyzp[3*A + 2];
            xyz2p[3*A + 2] = R2p[6] * xyzp[3*A + 0] + R2p[7] * xyzp[3*A + 1] + R2p[8] * xyzp[3*A + 2];
        }
        // Anisotropy weight
        if (cos2) {
            w *= pow(R2p[8], 2);
        }
        // Diffraction Intensity computation
        for (size_t P = 0; P < ss->shape()[0]; P++) {
            double sx = sp[3*P + 0];
            double sy = sp[3*P + 1];
            double sz = sp[3*P + 2];
            double s2 = pow(sx, 2) + pow(sy, 2) + pow(sz, 2);
            double D = 0.0;
            double NR = 0.0;
            double NI = 0.0;
            for (size_t A = 0; A < xyzs->shape()[0]; A++) {
                const double* a2p = ap + A * as->shape()[1];
                const double* b2p = bp + A * bs->shape()[1];
                double fA = 0.0;
                for (size_t J = 0; J < as->shape()[1]; J++) {
                    fA += a2p[J] * exp(-b2p[J] * s2 / pow(4.0 * M_PI, 2));
                }
                if (ued) {
                    fA = 1.0 / s2 * (Zp[A] - fA);
                }
                D += pow(fA, 2);
                double x = xyz2p[3*A + 0];
                double y = xyz2p[3*A + 1];
                double z = xyz2p[3*A + 2];
                double theta = sx * x + sy * y + sz * z;
                NR += fA * cos(theta);
                NI += fA * sin(theta);
            }
            double F = NR * NR + NI * NI;
            if (mod) {
                F = (F - D) / D;
            }
            Ip[P] += w * F;
        }
    } 
    
    return Is;
}

BOOST_PYTHON_MODULE(pyplugin)
{
    def("compute_diffraction", compute_diffraction);
}

