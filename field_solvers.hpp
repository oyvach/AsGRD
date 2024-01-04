#ifndef CONVTEST_HEADER
#define CONVTEST_HEADER
#endif

#include "prng_engine.hpp" // for use with ICs
#ifndef PI
#define PI 3.1415
#endif

#ifndef Cplx
#define Cplx Imag
#endif
using namespace std;
using namespace LATfield2;



template <class FieldType>
void update_aq_ptr(
    double dtau, double boxsize, double dx, double a,
    double fourpiG, const cosmology cosmo, const metadata sim,
    Field<FieldType> &nT_cdm, Field<FieldType> *achi,
    Field<FieldType> *aq)
{
    double Laplacian_achi, dqdtau, aeta, a2;
    double achif, aetanorm, prefac, dAcnf;

    dAcnf = cosmo.dAcnf;
    Site x(aq->lattice());
    a2 = a * a;

    aetanorm = cosmo.aetanorm / (a2 * a);
    prefac = cosmo.prefac * a2 * a2;
    for (x.first(); x.test(); x.next())
    {
        achif = (*achi)(x);
        Laplacian_achi = ((*achi)(x + 0) + (*achi)(x - 0) + (*achi)(x + 1) + (*achi)(x - 1) + (*achi)(x + 2) + (*achi)(x - 2) - 6. * achif) / (dx * dx);

        aeta = nT_cdm(x) * aetanorm / (1 + dAcnf * achif * achif);

        dqdtau = a2 * Laplacian_achi - prefac * ((pow(achif, 3.0) + cosmo.kappa2_as * achif * achif) + (aeta - 1.0) * achif); // extra terms
        (*aq)(x) = (*aq)(x) + dqdtau * dtau;
    }
}

template <class FieldType>
void update_achi_ptr(double dtau, double boxsize, double dx, double a, Field<FieldType> *achi,
                     Field<FieldType> *aq)
{
    Site x(achi->lattice());
    for (x.first(); x.test(); x.next())
    {
        (*achi)(x) = (*achi)(x) + (*aq)(x)*dtau / pow(a, 2.0);
    }
}

void doEvolveLF(double dtau, int N_as, double a_As, double fourpiG,
                metadata sim, double dx, double avgsource,
                cosmology cosmo, Field<Real> *achi, Field<Real> *aq,
                Field<Real> &nT_cdm)
{
    double dtauN = dtau / N_as;
    double dtauN_2 = dtauN / 2.;
    update_aq_ptr(
        dtauN_2, sim.boxsize, dx, a_As, fourpiG, cosmo, sim,
        nT_cdm, achi, aq);
    rungekutta4bg(a_As, fourpiG, dtauN_2, cosmo, avgsource);
    for (int i = 0; i < N_as; i++)
    {
        update_achi_ptr(dtauN, sim.boxsize, dx, a_As, achi, aq);
        achi->updateHalo();
        rungekutta4bg(a_As, fourpiG, dtauN_2, cosmo, avgsource);
        update_aq_ptr(
            dtauN, sim.boxsize, dx, a_As, fourpiG, cosmo, sim,
            nT_cdm, achi, aq);
        rungekutta4bg(a_As, fourpiG, dtauN_2, cosmo, avgsource);
    }
    rungekutta4bg(a_As, fourpiG, -dtauN_2, cosmo, avgsource);
    update_aq_ptr(
        -dtauN_2, sim.boxsize, dx, a_As, fourpiG, cosmo, sim,
        nT_cdm, achi, aq);
    aq->updateHalo();
}

void doEvolveEC(double dtau, int N_as, double a_As, double fourpiG,
                metadata sim, double dx, double avgsource,
                cosmology cosmo, Field<Real> *achi, Field<Real> *aq,
                Field<Real> &nT_cdm)
{
    double dtauN = dtau / N_as;
    // update velocity, then background, then field
    for (int i = 0; i < N_as; i++)
    {
        update_aq_ptr(
            dtauN, sim.boxsize, dx, a_As, fourpiG, cosmo, sim,
            nT_cdm, achi, aq);
        rungekutta4bg(a_As, fourpiG, dtauN, cosmo, avgsource);
        update_achi_ptr(dtauN, sim.boxsize,
                        dx, a_As, achi, aq);
        achi->updateHalo();
    }
    aq->updateHalo();
}

void doEvolveRK4(double dtau, int N_as, double a, double fourpiG,
                 metadata sim, double dx, double avgsource,
                 cosmology cosmo, Field<Real> *achi, Field<Real> *aq,
                 Field<Real> &nT_cdm, Field<Real> &achi2,
                 Field<Real> &aq2, Field<Real> &achi3, Field<Real> &aq3, Field<Real> &achiLap)
{
    double Laplacian_achi, dqdtau, aeta;
    double achif, aetanorm, prefac, dAcnf;
    dAcnf = 0.5 * pow(cosmo.mu_as / cosmo.M_as, 2) / cosmo.lambda_as;
    Site x(aq->lattice());
    double k1, k2, k3, k4;
    double a2 = a * a, k1a, k2a, k3a, k4a, ak2, ak3, ak4;
    double a2k2, a2k3, a2k4;
    double dtauN = dtau / N_as;
    aetanorm = cosmo.aetanorm;
    prefac = cosmo.prefac;
    double maxLap = 0, tmp;
    // store final solution in achi
    // initial field in achi3
    // store temporary solutions in achi2 / achiLap (need two for keeping neighbour points saved)
    for (x.first(); x.test(); x.next())
    {
        achi3(x) = (*achi)(x);
        aq3(x) = (*aq)(x);
    }
    achi3.updateHalo();

    for (int i = 0; i < N_as; i++)
    {
        k1a = a * Hconf(a, fourpiG, cosmo, avgsource);
        ak2 = a + dtauN * k1a / 2.;
        a2k2 = ak2 * ak2;
        k2a = ak2 * Hconf(ak2, fourpiG, cosmo, avgsource);
        ak3 = a + dtauN * k2a / 2.;
        a2k3 = ak3 * ak3;
        k3a = ak3 * Hconf(ak3, fourpiG, cosmo, avgsource);
        ak4 = a + dtauN * k3a;
        a2k4 = ak4 * ak4;
        k4a = ak4 * Hconf(ak4, fourpiG, cosmo, avgsource);

        // k1 + evolve to k2
        for (x.first(); x.test(); x.next())
        {
            achif = achi3(x);
            Laplacian_achi = (achi3(x + 0) + achi3(x - 0) + achi3(x + 1) + achi3(x - 1) + achi3(x + 2) + achi3(x - 2) - 6. * achi3(x)) / (dx * dx);
            aeta = nT_cdm(x) * aetanorm / a2 / a / (1 + dAcnf * achif * achif);
            k1 = a2 * Laplacian_achi - prefac * a2 * a2 * ((pow(achif, 3.0) + cosmo.kappa2_as * achif * achif) + (aeta - 1.0) * achif);
            (*achi)(x) = achi3(x) + (1. / 6.) * aq3(x) * dtauN / a2;
            (*aq)(x) = aq3(x) + (1. / 6.) * k1 * dtauN;
            achi2(x) = achi3(x) + dtauN * aq3(x) / 2 / a2;
            aq2(x) = aq3(x) + dtauN * k1 / 2;
            achiLap(x) = achi2(x);
        }
        achiLap.updateHalo();
        // k2 + evolve to k3
        for (x.first(); x.test(); x.next())
        {
            achif = achiLap(x);
            Laplacian_achi = (achiLap(x + 0) + achiLap(x - 0) + achiLap(x + 1) + achiLap(x - 1) + achiLap(x + 2) + achiLap(x - 2) - 6. * achiLap(x)) / (dx * dx);
            aeta = nT_cdm(x) * aetanorm / a2k2 / ak2 / (1 + dAcnf * achif * achif);
            k2 = a2k2 * Laplacian_achi - prefac * a2k2 * a2k2 * ((pow(achif, 3.0) + cosmo.kappa2_as * achif * achif) + (aeta - 1.0) * achif);
            (*achi)(x) = (*achi)(x) + (1. / 3.) * aq2(x) * dtauN / a2k2;
            (*aq)(x) = (*aq)(x) + (1. / 3.) * k2 * dtauN;
            achi2(x) = achi3(x) + dtauN * aq2(x) / 2 / a2k2;
            aq2(x) = aq3(x) + dtauN * k2 / 2;
        }
        achi2.updateHalo();
        // k3 + evolve to k4
        for (x.first(); x.test(); x.next())
        {
            achif = achi2(x);
            Laplacian_achi = (achi2(x + 0) + achi2(x - 0) + achi2(x + 1) + achi2(x - 1) + achi2(x + 2) + achi2(x - 2) - 6. * achif) / (dx * dx);
            aeta = nT_cdm(x) * aetanorm / a2k3 / ak3 / (1 + dAcnf * achif * achif);
            k3 = a2k3 * Laplacian_achi - prefac * a2k3 * a2k3 * ((pow(achif, 3.0) + cosmo.kappa2_as * achif * achif) + (aeta - 1.0) * achif);
            (*achi)(x) = (*achi)(x) + (1. / 3.) * aq2(x) * dtauN / a2k3;
            (*aq)(x) = (*aq)(x) + (1. / 3.) * k3 * dtauN;
            achiLap(x) = achi3(x) + dtauN * aq2(x) / a2k3;
            aq2(x) = aq3(x) + dtauN * k3;
        }
        achiLap.updateHalo();
        // k4 + final evolution
        for (x.first(); x.test(); x.next())
        {
            achif = achiLap(x);
            Laplacian_achi = (achiLap(x + 0) + achiLap(x - 0) + achiLap(x + 1) + achiLap(x - 1) + achiLap(x + 2) + achiLap(x - 2) - 6. * achif) / (dx * dx);
            aeta = nT_cdm(x) * aetanorm / a2k4 / ak4 / (1 + dAcnf * achif * achif);
            k4 = a2k4 * Laplacian_achi - prefac * a2k4 * a2k4 * ((pow(achif, 3.0) + cosmo.kappa2_as * achif * achif) + (aeta - 1.0) * achif);
            (*achi)(x) = (*achi)(x) + (1. / 6.) * aq2(x) * dtauN / a2k4;
            (*aq)(x) = (*aq)(x) + (1. / 6.) * k4 * dtauN;
            achi3(x) = (*achi)(x);
            aq3(x) = (*aq)(x);
        }
        achi->updateHalo();
        achi3.updateHalo();
        a += dtauN * (k1a + 2. * k2a + 2. * k3a + k4a) / 6.;
        a2 = a * a;
    }
    aq->updateHalo();
}
