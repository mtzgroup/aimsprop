#include <boost/python.hpp>
#include <lightspeed/tensor.hpp>

using namespace lightspeed;
using namespace boost::python;

/**
 * Compute the UED or X-Ray elastic scattering intensity for a molecule,
 * averaging over a molecular orientation quadrature grid, and accounting for
 * anisotropy and/or modifications to the scattering intensity definition. Uses
 * the independent atom model.
 *
 * First, the scattering angles theta are computed according to the elastic
 * scattering formula,
 *
 * s = 4 pi / lambda * sin(theta / 2) => theta = 2 arcsin(s * lambda / (4 * pi))
 *
 * This will throw if s * lambda / (4 * pi) is outside [+1, -1]
 *
 * Next, the scattering vectors are computed according to,
 *
 * \vec s = {sx, sy, sz} = s { cos(theta/2)sin(eta), sin(theta/2), cos(theta/2)cos(eta) }
 *  
 * Next the (isotropic) atomic scattering intensity is computed according to,
 *
 * D(s) = \sum_{A} | f_A (s) |^2
 *
 * Next, the scattering intensity is computed for each point in the molecular
 * orientation quadrature { R, w_R }, and summed up, including the possible
 * anisotropy weighting factor a(R).
 *
 * I(s, eta) = \sum_{R} w_{R} a(R) I_R (s, eta)
 *
 * Here the raw scattering intensity for a given molecular orientation is computed as,
 *
 * I_R (s, eta) = | \sum_{A} exp(i \vec s * (R^\dagger r_{A})) |^2
 *
 * e.g., the diffraction pattern of the rotated molecular geometry. 
 *
 * If the "anisotropy" option is set to true, the factor a(R) will be computed as
 * the square of the projection of the rotated z axis on the original z axis
 * (cos^2 zeta perpendicular one-photon pump-probe anisotropy)
 *
 * If the "mod" option is set to true, the modified scattering intensity will
 * replace the raw scattering intensity in the result:
 *
 * M_R (s, eta) = [I_R (s, eta) - D(s)] / D(s)
 *
 * @param lambda (double) - probe wavelength (e.g., deBroglie wavelength) [Angstrom].
 * @param ss (Tensor of shape (ns,)) - s values to collocate to [Angstrom^-1]
 * @param etas (Tensor of shape (neta,)) - eta values to collocate to [Radian]
 * @param xyzs (Tensor of shape (natom, 3)) - x, y, z centers of atoms [Anstrom]
 * @param fs (Tensor of shape (natom, ns)) - atomic form factors f_A (s)
 *  computed for the above atoms and s values. Specialization to UED or X-Ray
 *  can be made outside of the present code by computing these factors
 *  appropriately. E.g., f_A^UED (s) = 1 / s^2 (Z_A - f_A^XRAY (s)).
 * @param Rs (Tensor of shape (nrot, 3, 3)) - rotation matrices in the
 *  molecular orientation quadrature grid.
 * @param ws (Tensor of shape (nrot,)) - weights of rotation matrices in the
 *  molecular orientation quadrature grid.
 * @param anisotropy (bool) - a(R) will be computed as R_zz'^2 if true
 *  (perpendicular one-photon anisotropy), else a(R) will be set to 1 if false
 *  (parallel arrangement or isotropic pump-probe experiment).
 * @param mod (bool) - compute M(s,eta) if true, else compute I(s,eta) if false.
 *
 * @return I (Tensor of shape (ns, neta)) - molecular scattering intensity I(s, eta) 
 *  or M(s, eta)
 **/
lightspeed::shared_ptr<Tensor> compute_diffraction(
    double lambda,
    const lightspeed::shared_ptr<Tensor>& ss,
    const lightspeed::shared_ptr<Tensor>& etas,
    const lightspeed::shared_ptr<Tensor>& xyzs,
    const lightspeed::shared_ptr<Tensor>& fs,
    const lightspeed::shared_ptr<Tensor>& Rs,
    const lightspeed::shared_ptr<Tensor>& ws,
    bool anisotropy,
    bool mod
    )
{
    // Validity checks
    ss->ndim_error(1);
    etas->ndim_error(1);
    xyzs->ndim_error(2);
    xyzs->shape_error({xyzs->shape()[0], 3});
    fs->shape_error({xyzs->shape()[0], ss->shape()[0]});
    Rs->ndim_error(3);
    Rs->shape_error({Rs->shape()[0], 3, 3});
    ws->shape_error({Rs->shape()[0]});

    // Sizes
    size_t ns = ss->shape()[0];
    size_t neta = etas->shape()[0];
    size_t nA = xyzs->shape()[0];
    size_t nR = Rs->shape()[0];

    // Pointers
    const double* sp = ss->data().data();
    const double* etap = etas->data().data();
    const double* xyzp = xyzs->data().data();
    const double* fp = fs->data().data();
    const double* Rp = Rs->data().data();
    const double* wp = ws->data().data();

    // Scattering vectors
    lightspeed::shared_ptr<Tensor> sxyz(new Tensor({ns, neta, 3}));
    double* sxyzp = sxyz->data().data();
    for (size_t sind = 0; sind < ns; sind++) {
        double arg = lambda * sp[sind] / (4.0 * M_PI);
        if (arg > 1.0 or arg < -1.0) throw std::runtime_error("Invalid s:" + std::to_string(sp[sind]));
        double theta = 2.0 * asin(arg);
        for (size_t eind = 0; eind < neta; eind++) {
            sxyzp[3 * (sind * neta + eind) + 0] = sp[sind] * cos(theta / 2.0) * sin(etap[eind]);
            sxyzp[3 * (sind * neta + eind) + 1] = sp[sind] * sin(theta / 2.0);
            sxyzp[3 * (sind * neta + eind) + 2] = sp[sind] * cos(theta / 2.0) * cos(etap[eind]);
        }
    }

    // Atomic scattering intensity
    lightspeed::shared_ptr<Tensor> Ds(new Tensor({ns}));
    double* Dp = Ds->data().data();
    for (size_t sind = 0; sind < ns; sind++) {
        for (size_t A = 0; A < nA; A++) {
            Dp[sind] += pow(fp[A * ns + sind], 2);
        }
    }

    // Target
    lightspeed::shared_ptr<Tensor> Is(new Tensor({ns, neta}));
    double* Ip = Is->data().data();

    // #pragma omp parallel for schedule(static, 16)
    #pragma omp parallel for schedule(static)
    for (size_t P = 0; P < ns * neta; P++) {
        size_t sind = P / neta;
        size_t eind = P % neta;
        double sx = sxyzp[3*P + 0];
        double sy = sxyzp[3*P + 1];
        double sz = sxyzp[3*P + 2];
        double D = Dp[sind];
        double I = 0.0;
        for (size_t Rind = 0; Rind < Rs->shape()[0]; Rind++) {
            // Rotation quadrature
            const double* R2p = Rp + Rind * 9;
            double w = wp[Rind];
            // Anisotropy weight
            if (anisotropy) {
                w *= pow(R2p[8], 2);
            }
            // Diffraction cross section
            double NR = 0.0;
            double NI = 0.0;
            for (size_t A = 0; A < xyzs->shape()[0]; A++) {
                // Rotate coordinates to frame
                double x = R2p[0] * xyzp[3*A + 0] + R2p[1] * xyzp[3*A + 1] + R2p[2] * xyzp[3*A + 2];
                double y = R2p[3] * xyzp[3*A + 0] + R2p[4] * xyzp[3*A + 1] + R2p[5] * xyzp[3*A + 2];
                double z = R2p[6] * xyzp[3*A + 0] + R2p[7] * xyzp[3*A + 1] + R2p[8] * xyzp[3*A + 2];
                double theta = sx * x + sy * y + sz * z;
                double fA = fp[A * ns + sind];
                NR += fA * cos(theta);
                NI += fA * sin(theta);
            }
            double F = NR * NR + NI * NI;
            if (mod) {
                F = (F - D) / D;
            }
            I += w * F;
        }
        Ip[P] = I;
    }

    return Is;
}

BOOST_PYTHON_MODULE(pyplugin)
{
    def("compute_diffraction", compute_diffraction);
}

