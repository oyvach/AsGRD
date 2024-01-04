#ifndef GW_HEADER
#define GW_HEADER
#endif

#ifndef PI
#define PI 3.1415
#endif

#ifndef Cplx
#define Cplx Imag
#endif
using namespace std;
using namespace LATfield2;

//////////////////////////
// evolveFTtensor
//////////////////////////
// Description:
//   projection of the Fourier image of a tensor field on the transverse
//   trace-free tensor component
//
// Arguments:
//   SijFT      reference to the Fourier image of the input tensor field
//   hijFT      reference to the Fourier image of the current hij
//   hijprimeFT reference to the Fourier image of hij' that will be updated
//   hubble     conformal Hubble rate
//   dtau       next time step
//   dtau_old   previous time step
//
// Returns:
//
//////////////////////////
void (*evolveFTtensor)(Field<Cplx> &SijFT, Field<Cplx> &hijFT, Field<Cplx> &hijprimeFT, const double hubble, const double dtau, const double dtau_old);

void evolveFTtensorAdiabatic(Field<Cplx> &SijFT, Field<Cplx> &hijFT, Field<Cplx> &hijprimeFT, const double hubble, const double dtau, const double dtau_old)
{
    const int linesize = hijFT.lattice().size(1);
    int i;
    Real *gridk2;
    Cplx *kshift;
    rKSite k(hijFT.lattice());
    Real k2, k4;
    double dtau_mean = 0.5 * (dtau_old + dtau);

    gridk2 = (Real *)malloc(linesize * sizeof(Real));
    kshift = (Cplx *)malloc(linesize * sizeof(Cplx));

    for (i = 0; i < linesize; i++)
    {
        gridk2[i] = 2. * (Real)linesize * sin(M_PI * (Real)i / (Real)linesize);
        kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real)i / (Real)linesize), -sin(M_PI * (Real)i / (Real)linesize));
        gridk2[i] *= gridk2[i];
    }

    k.first();
    if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
    {
        for (i = 0; i < hijprimeFT.components(); i++)
        {
            hijprimeFT(k, i) = Cplx(0., 0.);
            hijFT(k, i) = Cplx(0., 0.);
        }

        k.next();
    }

    for (; k.test(); k.next())
    {
        k2 = gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)];
        k4 = k2 * k2 * linesize;

        if (k2 * dtau_mean * dtau_mean < 1.)
        {
            for (i = 0; i < hijprimeFT.components(); i++)
                hijFT(k, i) += hijprimeFT(k, i) * dtau_old;

            hijprimeFT(k, 0, 0) = ((1. - hubble * dtau_mean) * hijprimeFT(k, 0, 0) + (((gridk2[k.coord(0)] - k2) * ((gridk2[k.coord(0)] - k2) * SijFT(k, 0, 0) + 2. * kshift[k.coord(0)] * (kshift[k.coord(1)] * SijFT(k, 0, 1) + kshift[k.coord(2)] * SijFT(k, 0, 2))) + ((gridk2[k.coord(0)] + k2) * (gridk2[k.coord(1)] + k2) - 2. * k2 * k2) * SijFT(k, 1, 1) + ((gridk2[k.coord(0)] + k2) * (gridk2[k.coord(2)] + k2) - 2. * k2 * k2) * SijFT(k, 2, 2) + 2. * (gridk2[k.coord(0)] + k2) * kshift[k.coord(1)] * kshift[k.coord(2)] * SijFT(k, 1, 2)) / k4 - k2 * hijFT(k, 0, 0)) * dtau_mean) / (1. + hubble * dtau_mean);

            hijprimeFT(k, 0, 1) = ((1. - hubble * dtau_mean) * hijprimeFT(k, 0, 1) + ((2. * (gridk2[k.coord(0)] - k2) * (gridk2[k.coord(1)] - k2) * SijFT(k, 0, 1) + (gridk2[k.coord(2)] + k2) * kshift[k.coord(0)].conj() * kshift[k.coord(1)].conj() * SijFT(k, 2, 2) + (gridk2[k.coord(0)] - k2) * kshift[k.coord(1)].conj() * (kshift[k.coord(0)].conj() * SijFT(k, 0, 0) + 2. * kshift[k.coord(2)] * SijFT(k, 0, 2)) + (gridk2[k.coord(1)] - k2) * kshift[k.coord(0)].conj() * (kshift[k.coord(1)].conj() * SijFT(k, 1, 1) + 2. * kshift[k.coord(2)] * SijFT(k, 1, 2))) / k4 - k2 * hijFT(k, 0, 1)) * dtau_mean) / (1. + hubble * dtau_mean);

            hijprimeFT(k, 0, 2) = ((1. - hubble * dtau_mean) * hijprimeFT(k, 0, 2) + ((2. * (gridk2[k.coord(0)] - k2) * (gridk2[k.coord(2)] - k2) * SijFT(k, 0, 2) + (gridk2[k.coord(1)] + k2) * kshift[k.coord(0)].conj() * kshift[k.coord(2)].conj() * SijFT(k, 1, 1) + (gridk2[k.coord(0)] - k2) * kshift[k.coord(2)].conj() * (kshift[k.coord(0)].conj() * SijFT(k, 0, 0) + 2. * kshift[k.coord(1)] * SijFT(k, 0, 1)) + (gridk2[k.coord(2)] - k2) * kshift[k.coord(0)].conj() * (kshift[k.coord(2)].conj() * SijFT(k, 2, 2) + 2. * kshift[k.coord(1)] * SijFT(k, 1, 2))) / k4 - k2 * hijFT(k, 0, 2)) * dtau_mean) / (1. + hubble * dtau_mean);

            hijprimeFT(k, 1, 1) = ((1. - hubble * dtau_mean) * hijprimeFT(k, 1, 1) + (((gridk2[k.coord(1)] - k2) * ((gridk2[k.coord(1)] - k2) * SijFT(k, 1, 1) + 2. * kshift[k.coord(1)] * (kshift[k.coord(0)] * SijFT(k, 0, 1) + kshift[k.coord(2)] * SijFT(k, 1, 2))) + ((gridk2[k.coord(1)] + k2) * (gridk2[k.coord(0)] + k2) - 2. * k2 * k2) * SijFT(k, 0, 0) + ((gridk2[k.coord(1)] + k2) * (gridk2[k.coord(2)] + k2) - 2. * k2 * k2) * SijFT(k, 2, 2) + 2. * (gridk2[k.coord(1)] + k2) * kshift[k.coord(0)] * kshift[k.coord(2)] * SijFT(k, 0, 2)) / k4 - k2 * hijFT(k, 1, 1)) * dtau_mean) / (1. + hubble * dtau_mean);

            hijprimeFT(k, 1, 2) = ((1. - hubble * dtau_mean) * hijprimeFT(k, 1, 2) + ((2. * (gridk2[k.coord(1)] - k2) * (gridk2[k.coord(2)] - k2) * SijFT(k, 1, 2) + (gridk2[k.coord(0)] + k2) * kshift[k.coord(1)].conj() * kshift[k.coord(2)].conj() * SijFT(k, 0, 0) + (gridk2[k.coord(1)] - k2) * kshift[k.coord(2)].conj() * (kshift[k.coord(1)].conj() * SijFT(k, 1, 1) + 2. * kshift[k.coord(0)] * SijFT(k, 0, 1)) + (gridk2[k.coord(2)] - k2) * kshift[k.coord(1)].conj() * (kshift[k.coord(2)].conj() * SijFT(k, 2, 2) + 2. * kshift[k.coord(0)] * SijFT(k, 0, 2))) / k4 - k2 * hijFT(k, 1, 2)) * dtau_mean) / (1. + hubble * dtau_mean);

            hijprimeFT(k, 2, 2) = ((1. - hubble * dtau_mean) * hijprimeFT(k, 2, 2) + (((gridk2[k.coord(2)] - k2) * ((gridk2[k.coord(2)] - k2) * SijFT(k, 2, 2) + 2. * kshift[k.coord(2)] * (kshift[k.coord(0)] * SijFT(k, 0, 2) + kshift[k.coord(1)] * SijFT(k, 1, 2))) + ((gridk2[k.coord(2)] + k2) * (gridk2[k.coord(0)] + k2) - 2. * k2 * k2) * SijFT(k, 0, 0) + ((gridk2[k.coord(2)] + k2) * (gridk2[k.coord(1)] + k2) - 2. * k2 * k2) * SijFT(k, 1, 1) + 2. * (gridk2[k.coord(2)] + k2) * kshift[k.coord(0)] * kshift[k.coord(1)] * SijFT(k, 0, 1)) / k4 - k2 * hijFT(k, 2, 2)) * dtau_mean) / (1. + hubble * dtau_mean);
        }
        else
        {

            // store temporarily
            hijprimeFT(k, 0, 0) = (((gridk2[k.coord(0)] - k2) * ((gridk2[k.coord(0)] - k2) * SijFT(k, 0, 0) + 2. * kshift[k.coord(0)] * (kshift[k.coord(1)] * SijFT(k, 0, 1) + kshift[k.coord(2)] * SijFT(k, 0, 2))) + ((gridk2[k.coord(0)] + k2) * (gridk2[k.coord(1)] + k2) - 2. * k2 * k2) * SijFT(k, 1, 1) + ((gridk2[k.coord(0)] + k2) * (gridk2[k.coord(2)] + k2) - 2. * k2 * k2) * SijFT(k, 2, 2) + 2. * (gridk2[k.coord(0)] + k2) * kshift[k.coord(1)] * kshift[k.coord(2)] * SijFT(k, 1, 2)) * dtau_old / k4 + 2. * hubble * hijFT(k, 0, 0)) / (k2 * dtau_old + 2. * hubble);

            hijprimeFT(k, 0, 1) = ((2. * (gridk2[k.coord(0)] - k2) * (gridk2[k.coord(1)] - k2) * SijFT(k, 0, 1) + (gridk2[k.coord(2)] + k2) * kshift[k.coord(0)].conj() * kshift[k.coord(1)].conj() * SijFT(k, 2, 2) + (gridk2[k.coord(0)] - k2) * kshift[k.coord(1)].conj() * (kshift[k.coord(0)].conj() * SijFT(k, 0, 0) + 2. * kshift[k.coord(2)] * SijFT(k, 0, 2)) + (gridk2[k.coord(1)] - k2) * kshift[k.coord(0)].conj() * (kshift[k.coord(1)].conj() * SijFT(k, 1, 1) + 2. * kshift[k.coord(2)] * SijFT(k, 1, 2))) * dtau_old / k4 + 2. * hubble * hijFT(k, 0, 1)) / (k2 * dtau_old + 2. * hubble);

            hijprimeFT(k, 0, 2) = ((2. * (gridk2[k.coord(0)] - k2) * (gridk2[k.coord(2)] - k2) * SijFT(k, 0, 2) + (gridk2[k.coord(1)] + k2) * kshift[k.coord(0)].conj() * kshift[k.coord(2)].conj() * SijFT(k, 1, 1) + (gridk2[k.coord(0)] - k2) * kshift[k.coord(2)].conj() * (kshift[k.coord(0)].conj() * SijFT(k, 0, 0) + 2. * kshift[k.coord(1)] * SijFT(k, 0, 1)) + (gridk2[k.coord(2)] - k2) * kshift[k.coord(0)].conj() * (kshift[k.coord(2)].conj() * SijFT(k, 2, 2) + 2. * kshift[k.coord(1)] * SijFT(k, 1, 2))) * dtau_old / k4 + 2. * hubble * hijFT(k, 0, 2)) / (k2 * dtau_old + 2. * hubble);

            hijprimeFT(k, 1, 1) = (((gridk2[k.coord(1)] - k2) * ((gridk2[k.coord(1)] - k2) * SijFT(k, 1, 1) + 2. * kshift[k.coord(1)] * (kshift[k.coord(0)] * SijFT(k, 0, 1) + kshift[k.coord(2)] * SijFT(k, 1, 2))) + ((gridk2[k.coord(1)] + k2) * (gridk2[k.coord(0)] + k2) - 2. * k2 * k2) * SijFT(k, 0, 0) + ((gridk2[k.coord(1)] + k2) * (gridk2[k.coord(2)] + k2) - 2. * k2 * k2) * SijFT(k, 2, 2) + 2. * (gridk2[k.coord(1)] + k2) * kshift[k.coord(0)] * kshift[k.coord(2)] * SijFT(k, 0, 2)) * dtau_old / k4 + 2. * hubble * hijFT(k, 1, 1)) / (k2 * dtau_old + 2. * hubble);

            hijprimeFT(k, 1, 2) = ((2. * (gridk2[k.coord(1)] - k2) * (gridk2[k.coord(2)] - k2) * SijFT(k, 1, 2) + (gridk2[k.coord(0)] + k2) * kshift[k.coord(1)].conj() * kshift[k.coord(2)].conj() * SijFT(k, 0, 0) + (gridk2[k.coord(1)] - k2) * kshift[k.coord(2)].conj() * (kshift[k.coord(1)].conj() * SijFT(k, 1, 1) + 2. * kshift[k.coord(0)] * SijFT(k, 0, 1)) + (gridk2[k.coord(2)] - k2) * kshift[k.coord(1)].conj() * (kshift[k.coord(2)].conj() * SijFT(k, 2, 2) + 2. * kshift[k.coord(0)] * SijFT(k, 0, 2))) * dtau_old / k4 + 2. * hubble * hijFT(k, 1, 2)) / (k2 * dtau_old + 2. * hubble);

            hijprimeFT(k, 2, 2) = (((gridk2[k.coord(2)] - k2) * ((gridk2[k.coord(2)] - k2) * SijFT(k, 2, 2) + 2. * kshift[k.coord(2)] * (kshift[k.coord(0)] * SijFT(k, 0, 2) + kshift[k.coord(1)] * SijFT(k, 1, 2))) + ((gridk2[k.coord(2)] + k2) * (gridk2[k.coord(0)] + k2) - 2. * k2 * k2) * SijFT(k, 0, 0) + ((gridk2[k.coord(2)] + k2) * (gridk2[k.coord(1)] + k2) - 2. * k2 * k2) * SijFT(k, 1, 1) + 2. * (gridk2[k.coord(2)] + k2) * kshift[k.coord(0)] * kshift[k.coord(1)] * SijFT(k, 0, 1)) * dtau_old / k4 + 2. * hubble * hijFT(k, 2, 2)) / (k2 * dtau_old + 2. * hubble);

            for (i = 0; i < hijprimeFT.components(); i++)
            {
                hijprimeFT(k, i) = (hijprimeFT(k, i) - hijFT(k, i)) / dtau_old; // adiabatic velocity
                hijFT(k, i) += hijprimeFT(k, i) * dtau_old;                     // from above
            }
        }
    }

    free(gridk2);
    free(kshift);
}


void evolveFTtensorZero(Field<Cplx> &SijFT, Field<Cplx> &hijFT, Field<Cplx> &hijprimeFT, const double hubble, const double dtau, const double dtau_old)
{
    const int linesize = hijFT.lattice().size(1);
    int i;
    Real *gridk2;
    Cplx *kshift;
    rKSite k(hijFT.lattice());
    Real k2, k4;
    double dtau_mean = 0.5 * (dtau_old + dtau);

    gridk2 = (Real *)malloc(linesize * sizeof(Real));
    kshift = (Cplx *)malloc(linesize * sizeof(Cplx));

    for (i = 0; i < linesize; i++)
    {
        gridk2[i] = 2. * (Real)linesize * sin(M_PI * (Real)i / (Real)linesize);
        kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real)i / (Real)linesize), -sin(M_PI * (Real)i / (Real)linesize));
        gridk2[i] *= gridk2[i];
    }

    k.first();
    if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
    {
        for (i = 0; i < hijprimeFT.components(); i++)
        {
            hijprimeFT(k, i) = Cplx(0., 0.);
            hijFT(k, i) = Cplx(0., 0.);
        }

        k.next();
    }

    for (; k.test(); k.next())
    {
        k2 = gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)];
        k4 = k2 * k2 * linesize;

        if (k2 * dtau_mean * dtau_mean < 1.)
        {
            for (i = 0; i < hijprimeFT.components(); i++)
                hijFT(k, i) += hijprimeFT(k, i) * dtau_old;

            hijprimeFT(k, 0, 0) = ((1. - hubble * dtau_mean) * hijprimeFT(k, 0, 0) + (((gridk2[k.coord(0)] - k2) * ((gridk2[k.coord(0)] - k2) * SijFT(k, 0, 0) + 2. * kshift[k.coord(0)] * (kshift[k.coord(1)] * SijFT(k, 0, 1) + kshift[k.coord(2)] * SijFT(k, 0, 2))) + ((gridk2[k.coord(0)] + k2) * (gridk2[k.coord(1)] + k2) - 2. * k2 * k2) * SijFT(k, 1, 1) + ((gridk2[k.coord(0)] + k2) * (gridk2[k.coord(2)] + k2) - 2. * k2 * k2) * SijFT(k, 2, 2) + 2. * (gridk2[k.coord(0)] + k2) * kshift[k.coord(1)] * kshift[k.coord(2)] * SijFT(k, 1, 2)) / k4 - k2 * hijFT(k, 0, 0)) * dtau_mean) / (1. + hubble * dtau_mean);

            hijprimeFT(k, 0, 1) = ((1. - hubble * dtau_mean) * hijprimeFT(k, 0, 1) + ((2. * (gridk2[k.coord(0)] - k2) * (gridk2[k.coord(1)] - k2) * SijFT(k, 0, 1) + (gridk2[k.coord(2)] + k2) * kshift[k.coord(0)].conj() * kshift[k.coord(1)].conj() * SijFT(k, 2, 2) + (gridk2[k.coord(0)] - k2) * kshift[k.coord(1)].conj() * (kshift[k.coord(0)].conj() * SijFT(k, 0, 0) + 2. * kshift[k.coord(2)] * SijFT(k, 0, 2)) + (gridk2[k.coord(1)] - k2) * kshift[k.coord(0)].conj() * (kshift[k.coord(1)].conj() * SijFT(k, 1, 1) + 2. * kshift[k.coord(2)] * SijFT(k, 1, 2))) / k4 - k2 * hijFT(k, 0, 1)) * dtau_mean) / (1. + hubble * dtau_mean);

            hijprimeFT(k, 0, 2) = ((1. - hubble * dtau_mean) * hijprimeFT(k, 0, 2) + ((2. * (gridk2[k.coord(0)] - k2) * (gridk2[k.coord(2)] - k2) * SijFT(k, 0, 2) + (gridk2[k.coord(1)] + k2) * kshift[k.coord(0)].conj() * kshift[k.coord(2)].conj() * SijFT(k, 1, 1) + (gridk2[k.coord(0)] - k2) * kshift[k.coord(2)].conj() * (kshift[k.coord(0)].conj() * SijFT(k, 0, 0) + 2. * kshift[k.coord(1)] * SijFT(k, 0, 1)) + (gridk2[k.coord(2)] - k2) * kshift[k.coord(0)].conj() * (kshift[k.coord(2)].conj() * SijFT(k, 2, 2) + 2. * kshift[k.coord(1)] * SijFT(k, 1, 2))) / k4 - k2 * hijFT(k, 0, 2)) * dtau_mean) / (1. + hubble * dtau_mean);

            hijprimeFT(k, 1, 1) = ((1. - hubble * dtau_mean) * hijprimeFT(k, 1, 1) + (((gridk2[k.coord(1)] - k2) * ((gridk2[k.coord(1)] - k2) * SijFT(k, 1, 1) + 2. * kshift[k.coord(1)] * (kshift[k.coord(0)] * SijFT(k, 0, 1) + kshift[k.coord(2)] * SijFT(k, 1, 2))) + ((gridk2[k.coord(1)] + k2) * (gridk2[k.coord(0)] + k2) - 2. * k2 * k2) * SijFT(k, 0, 0) + ((gridk2[k.coord(1)] + k2) * (gridk2[k.coord(2)] + k2) - 2. * k2 * k2) * SijFT(k, 2, 2) + 2. * (gridk2[k.coord(1)] + k2) * kshift[k.coord(0)] * kshift[k.coord(2)] * SijFT(k, 0, 2)) / k4 - k2 * hijFT(k, 1, 1)) * dtau_mean) / (1. + hubble * dtau_mean);

            hijprimeFT(k, 1, 2) = ((1. - hubble * dtau_mean) * hijprimeFT(k, 1, 2) + ((2. * (gridk2[k.coord(1)] - k2) * (gridk2[k.coord(2)] - k2) * SijFT(k, 1, 2) + (gridk2[k.coord(0)] + k2) * kshift[k.coord(1)].conj() * kshift[k.coord(2)].conj() * SijFT(k, 0, 0) + (gridk2[k.coord(1)] - k2) * kshift[k.coord(2)].conj() * (kshift[k.coord(1)].conj() * SijFT(k, 1, 1) + 2. * kshift[k.coord(0)] * SijFT(k, 0, 1)) + (gridk2[k.coord(2)] - k2) * kshift[k.coord(1)].conj() * (kshift[k.coord(2)].conj() * SijFT(k, 2, 2) + 2. * kshift[k.coord(0)] * SijFT(k, 0, 2))) / k4 - k2 * hijFT(k, 1, 2)) * dtau_mean) / (1. + hubble * dtau_mean);

            hijprimeFT(k, 2, 2) = ((1. - hubble * dtau_mean) * hijprimeFT(k, 2, 2) + (((gridk2[k.coord(2)] - k2) * ((gridk2[k.coord(2)] - k2) * SijFT(k, 2, 2) + 2. * kshift[k.coord(2)] * (kshift[k.coord(0)] * SijFT(k, 0, 2) + kshift[k.coord(1)] * SijFT(k, 1, 2))) + ((gridk2[k.coord(2)] + k2) * (gridk2[k.coord(0)] + k2) - 2. * k2 * k2) * SijFT(k, 0, 0) + ((gridk2[k.coord(2)] + k2) * (gridk2[k.coord(1)] + k2) - 2. * k2 * k2) * SijFT(k, 1, 1) + 2. * (gridk2[k.coord(2)] + k2) * kshift[k.coord(0)] * kshift[k.coord(1)] * SijFT(k, 0, 1)) / k4 - k2 * hijFT(k, 2, 2)) * dtau_mean) / (1. + hubble * dtau_mean);
        }
        else
        {

            // store temporarily
            hijFT(k, 0, 0) = (((gridk2[k.coord(0)] - k2) * ((gridk2[k.coord(0)] - k2) * SijFT(k, 0, 0) + 2. * kshift[k.coord(0)] * (kshift[k.coord(1)] * SijFT(k, 0, 1) + kshift[k.coord(2)] * SijFT(k, 0, 2))) + ((gridk2[k.coord(0)] + k2) * (gridk2[k.coord(1)] + k2) - 2. * k2 * k2) * SijFT(k, 1, 1) + ((gridk2[k.coord(0)] + k2) * (gridk2[k.coord(2)] + k2) - 2. * k2 * k2) * SijFT(k, 2, 2) + 2. * (gridk2[k.coord(0)] + k2) * kshift[k.coord(1)] * kshift[k.coord(2)] * SijFT(k, 1, 2)) * dtau_old / k4 + 2. * hubble * hijFT(k, 0, 0)) / (k2 * dtau_old + 2. * hubble);

            hijFT(k, 0, 1) = ((2. * (gridk2[k.coord(0)] - k2) * (gridk2[k.coord(1)] - k2) * SijFT(k, 0, 1) + (gridk2[k.coord(2)] + k2) * kshift[k.coord(0)].conj() * kshift[k.coord(1)].conj() * SijFT(k, 2, 2) + (gridk2[k.coord(0)] - k2) * kshift[k.coord(1)].conj() * (kshift[k.coord(0)].conj() * SijFT(k, 0, 0) + 2. * kshift[k.coord(2)] * SijFT(k, 0, 2)) + (gridk2[k.coord(1)] - k2) * kshift[k.coord(0)].conj() * (kshift[k.coord(1)].conj() * SijFT(k, 1, 1) + 2. * kshift[k.coord(2)] * SijFT(k, 1, 2))) * dtau_old / k4 + 2. * hubble * hijFT(k, 0, 1)) / (k2 * dtau_old + 2. * hubble);

            hijFT(k, 0, 2) = ((2. * (gridk2[k.coord(0)] - k2) * (gridk2[k.coord(2)] - k2) * SijFT(k, 0, 2) + (gridk2[k.coord(1)] + k2) * kshift[k.coord(0)].conj() * kshift[k.coord(2)].conj() * SijFT(k, 1, 1) + (gridk2[k.coord(0)] - k2) * kshift[k.coord(2)].conj() * (kshift[k.coord(0)].conj() * SijFT(k, 0, 0) + 2. * kshift[k.coord(1)] * SijFT(k, 0, 1)) + (gridk2[k.coord(2)] - k2) * kshift[k.coord(0)].conj() * (kshift[k.coord(2)].conj() * SijFT(k, 2, 2) + 2. * kshift[k.coord(1)] * SijFT(k, 1, 2))) * dtau_old / k4 + 2. * hubble * hijFT(k, 0, 2)) / (k2 * dtau_old + 2. * hubble);

            hijFT(k, 1, 1) = (((gridk2[k.coord(1)] - k2) * ((gridk2[k.coord(1)] - k2) * SijFT(k, 1, 1) + 2. * kshift[k.coord(1)] * (kshift[k.coord(0)] * SijFT(k, 0, 1) + kshift[k.coord(2)] * SijFT(k, 1, 2))) + ((gridk2[k.coord(1)] + k2) * (gridk2[k.coord(0)] + k2) - 2. * k2 * k2) * SijFT(k, 0, 0) + ((gridk2[k.coord(1)] + k2) * (gridk2[k.coord(2)] + k2) - 2. * k2 * k2) * SijFT(k, 2, 2) + 2. * (gridk2[k.coord(1)] + k2) * kshift[k.coord(0)] * kshift[k.coord(2)] * SijFT(k, 0, 2)) * dtau_old / k4 + 2. * hubble * hijFT(k, 1, 1)) / (k2 * dtau_old + 2. * hubble);

            hijFT(k, 1, 2) = ((2. * (gridk2[k.coord(1)] - k2) * (gridk2[k.coord(2)] - k2) * SijFT(k, 1, 2) + (gridk2[k.coord(0)] + k2) * kshift[k.coord(1)].conj() * kshift[k.coord(2)].conj() * SijFT(k, 0, 0) + (gridk2[k.coord(1)] - k2) * kshift[k.coord(2)].conj() * (kshift[k.coord(1)].conj() * SijFT(k, 1, 1) + 2. * kshift[k.coord(0)] * SijFT(k, 0, 1)) + (gridk2[k.coord(2)] - k2) * kshift[k.coord(1)].conj() * (kshift[k.coord(2)].conj() * SijFT(k, 2, 2) + 2. * kshift[k.coord(0)] * SijFT(k, 0, 2))) * dtau_old / k4 + 2. * hubble * hijFT(k, 1, 2)) / (k2 * dtau_old + 2. * hubble);

            hijFT(k, 2, 2) = (((gridk2[k.coord(2)] - k2) * ((gridk2[k.coord(2)] - k2) * SijFT(k, 2, 2) + 2. * kshift[k.coord(2)] * (kshift[k.coord(0)] * SijFT(k, 0, 2) + kshift[k.coord(1)] * SijFT(k, 1, 2))) + ((gridk2[k.coord(2)] + k2) * (gridk2[k.coord(0)] + k2) - 2. * k2 * k2) * SijFT(k, 0, 0) + ((gridk2[k.coord(2)] + k2) * (gridk2[k.coord(1)] + k2) - 2. * k2 * k2) * SijFT(k, 1, 1) + 2. * (gridk2[k.coord(2)] + k2) * kshift[k.coord(0)] * kshift[k.coord(1)] * SijFT(k, 0, 1)) * dtau_old / k4 + 2. * hubble * hijFT(k, 2, 2)) / (k2 * dtau_old + 2. * hubble);

            // for (i = 0; i < hijprimeFT.components(); i++)
            // {
            //     hijprimeFT(k, i) =         kept zero          
            // }
        }
    }

    free(gridk2);
    free(kshift);
}
