#ifndef INTERFACE_F_HEADER
#define INTERFACE_F_HEADER
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

void doSetBackground(const metadata sim)
{
    if (strcmp(sim.bckopt, "sym") == 0 || strcmp(sim.bckopt, "symdyn") == 0)
        Hconf = &Hconf_sym;
    else if (strcmp(sim.bckopt, "rho") == 0)
        Hconf = &Hconf_rho;
    else if (strcmp(sim.bckopt, "mink") == 0)
        Hconf = &Hconf_mink;
    else // i.e. lcdm
        Hconf = &Hconf_old;
}

void doUpdateBackgroundEnergy(cosmology &cosmo, const metadata sim, double a,
                              double achiB, double aqB)
{
    // NB: only made for symmetron (deltaBeta=0)
    if (sim.As_use_dynbck == 0)
    {
        double phib;
        if (a <= cosmo.astar)
            cosmo.Omega_sym = 0.;
        else
        {
            double a3 = pow(a, 3);
            phib = (cosmo.mu_as / sqrt(cosmo.lambda_as)) * sqrt(1 - pow(cosmo.astar, 3) / a3);
            cosmo.Omega_sym = 0.5 * cosmo.Omega_m / a3 * pow(phib / cosmo.M_as, 2) + (-0.5 * pow(cosmo.mu_as * phib, 2) + 0.25 * cosmo.lambda_as * pow(phib, 4)) * 2998 * 2998 / 3;
        }
    }
    else
    {
        double a3 = a * a * a;
        double phifac = cosmo.mu_as / sqrt(cosmo.lambda_as);
        cosmo.Omega_sym =
            +0.5 * cosmo.Omega_m / a3 * pow(achiB * phifac / cosmo.M_as, 2) + (0.5 * pow(aqB * phifac / (a3 * sim.boxsize), 2) - 0.5 * pow(cosmo.mu_as * achiB * phifac, 2) + 0.25 * cosmo.lambda_as * pow(achiB * phifac, 4)) * 2998 * 2998 / 3;
    }
}

double sign(double val)
{
    if (val < 0)
        return -1;
    else
        return 1;
}

void doAchiQSA(
    double dtau, double a_As, double fourpiG,
    metadata sim, double dx, double avgsource,
    cosmology cosmo, Field<Real> &achi, Field<Real> &aq,
    Field<Real> &nT_cdm, double achiB, int cycle, int ICopt, int bckgRho)
{
    Site x(achi.lattice());
    // Quasi-static approximation using Gauss-Seidel iterative solution of Poisson
    double convcheck = 0.9; // pass if GStol = 1
    double old, QSAsource, QSAsource2;
    int count = 0;
    double prefac = pow(a_As * cosmo.mu_as * sim.boxsize * dx, 2.0);
    double aetanorm = 2. * fourpiG / (a_As * pow(sim.boxsize * a_As * cosmo.M_as * cosmo.mu_as, 2.0));
    double dAcnf = 0.5 * pow(cosmo.mu_as, 2) / cosmo.lambda_as;
    double ffac = 1. / sim.numpts;
    double half = 1. / 2., cx;
    double T = cosmo.Omega_m / pow(a_As, 3); // use if bckgRho
    double GStol = sim.GStol;

    if (a_As > cosmo.astar)
    {
        if (sim.QSA_guess_dyn == 0)
            old = sqrt(1.0 - pow(cosmo.astar / a_As, 3.0));
        else
            old = achiB;
    }
    else
        old = 1e-4; // done in ISIS
    if (ICopt == 0) // homogeneous ICs
    {
        for (x.first(); x.test(); x.next())
        {
            achi(x) = old;
        }
    }
    else if (ICopt == 1) // try to have domain wall across box
    {
        for (x.first(); x.test(); x.next())
        {
            cx = x.coord(2) * ffac;
            if (cx < half)
                achi(x) = -old;
            else
            {
                achi(x) = old;
            }
        }
    }
    achi.updateHalo();
    if (parallel.rank() == 0 and cycle % CYCLE_INFO_INTERVAL == 0)
        cout << "Starting Gauss-Seidel sweeps" << endl;
    while (convcheck > GStol)
    {
        convcheck = 0.;
        count++;
        for (x.first(); x.test(); x.next())
        {
            old = achi(x);
            if (bckgRho < 1)
                T = nT_cdm(x);
            QSAsource = pow(old, 3.0) + (aetanorm * T - 1.0) * old + cosmo.kappa2_as * old * old;
            QSAsource2 = 3. * old * old + (aetanorm * T / (1 + dAcnf * old * old) - 1.0) //
                         + 2. * cosmo.kappa2_as * old;

            achi(x) = old - (achi(x + 0) + achi(x + 1) + achi(x + 2) + achi(x - 0) + achi(x - 1) + achi(x - 2) //
                             - 6. * old - QSAsource * prefac) /
                                (-6. - prefac * QSAsource2);
            
            if (fabs(achi(x) - old) > convcheck)
                convcheck = fabs(achi(x) - old);
        }
        parallel.max(convcheck);
        achi.updateHalo();
        if (parallel.rank() == 0 and cycle % CYCLE_INFO_INTERVAL == 0 and count % 1000 == 0)
            cout << "finished GS sweep " << count << endl
                 << "   error = " << convcheck << ", tol = " << sim.GStol << endl;
    }
    if (parallel.rank() == 0 and cycle % CYCLE_INFO_INTERVAL == 0)
        cout << "finished GS sweeps\n";
    // find velocity
    if (sim.QSAq > 0)
    {
        for (x.first(); x.test(); x.next())
        {
            aq(x) = a_As * a_As * (achi(x) - aq(x)) / dtau;
        }
    }
    if (sim.QSAmodulateAmp != 1)
    {
        for (x.first(); x.test(); x.next())
        {
            achi(x) *= sim.QSAmodulateAmp;
        }
        achi.updateHalo();
    }
}

void doAchiICs(metadata sim, icsettings ic, cosmology cosmo,
               Field<Real> &achi, Field<Real> &aq, Field<Cplx> &scalarFT,
               PlanFFT<Cplx> &plan_achi,
               double fourpiG, double avgsource, double a,
               double &achiB, double &aqB, int numpts3d,
               double dtau, double dx, int cycle,
               Field<Real> &nT_cdm)
{
    Site x(achi.lattice());
    if (strcmp(sim.achifile, "no") != 0)
    {
        COUT << "Reading achi/aq from file " << sim.achifile << endl;
        string filename = sim.achifile;
        achi.loadHDF5(filename);
        if (strcmp(sim.aqfile, "no") != 0)
        {
            filename = sim.aqfile;
            aq.loadHDF5(filename);
        }
        double mx = 0, mn = 0;
        for (x.first(); x.test(); x.next())
        {
            if (achi(x) > mx)
                mx = achi(x);
            if (achi(x) < mn)
                mn = achi(x);
        }
        parallel.max<double>(mx);
        parallel.min<double>(mn);

        COUT << "achi max/min = " << mx << " / " << mn << endl;
    }
    else if (strcmp(sim.As_ICopt, "homogeneous") == 0)
    {
        COUT << "Initialising achi as random small perturbations." << endl;

        for (x.first(); x.test(); x.next())
        {
            achi(x) = sim.ICamp;
            aq(x) = sim.ICamp2;
        }
    }
    // some of below has read achi first
    else if (strcmp(sim.As_ICopt, "random") == 0)
    {
        COUT << "Initialising achi as random small perturbations." << endl;
        // Set x-grid small real numbers
        sitmo::prng_engine prng;
        prng.seed(ic.seed * parallel.rank());
        for (x.first(); x.test(); x.next())
        {
            achi(x) = (2.0 * ((double)prng()) / (double)sitmo::prng_engine::max() - 1.0) * sim.ICamp;
            aq(x) = (2.0 * ((double)prng()) / (double)sitmo::prng_engine::max() - 1.0) * sim.ICamp2 * sim.boxsize; // multiply because dimensionfull
        }
    }
    else if (strcmp(sim.As_ICopt, "relaxed") == 0)
    { // Homogeneous field relaxed around ptcl input
        // solve field
        COUT << "Starting field with relaxed ICs\n";
        doAchiQSA(
            dtau, a, fourpiG,
            sim, dx, avgsource,
            cosmo, achi, aq,
            nT_cdm, achiB, cycle, 0, 0);
    }
    else if (strcmp(sim.As_ICopt, "relaxedWall") == 0)
    { // have initial domain wall across x-axis
        // solve field
        COUT << "Starting field with relaxed ICs and a wall\n";
        doAchiQSA(
            dtau, a, fourpiG,
            sim, dx, avgsource,
            cosmo, achi, aq,
            nT_cdm, achiB, cycle, 1, 0);
    }
    else if (strcmp(sim.As_ICopt, "relaxedSnap") == 0)
    { // relax previously loaded snapshot
        COUT << "Starting field with relaxing input from file\n";
        // solve field
        doAchiQSA(
            dtau, a, fourpiG,
            sim, dx, avgsource,
            cosmo, achi, aq,
            nT_cdm, achiB, cycle, 2, 0);
    }
    else // proper ICs from scale-invariant power spectrum
    {    // should be correspondent at different resolutions
        COUT << "Initialising achi from scale-invariant power spectrum" << endl;
        // Set x-grid small real numbers
        gsl_spline *pkspline = NULL; // transfer function?
        int num = 10000;             // increase dim>ngrid
        double k[num], pk[num];
        k[0] = 0.;
        pk[0] = sim.ICamp;
        for (int i = 1; i < num; i++)
        {
            k[i] = (double)((i + 1) * sim.numpts) * 2. * M_PI / (double)(num + 1);
            pk[i] = sim.ICamp / sqrt(2.) / M_PI; // scale-invariant
        }
        pkspline = gsl_spline_alloc(gsl_interp_cspline, num);
        gsl_spline_init(pkspline, k, pk, num);
        generateRealization(scalarFT, 1,
                            pkspline, (unsigned int)ic.seed, ic.flags & ICFLAG_KSPHERE);

        plan_achi.execute(FFT_BACKWARD); // leave velocity zero
    }

    achi.updateHalo();
    aq.updateHalo();
    if (sim.As_compute_dynbck > 0)
    {
        for (x.first(); x.test(); x.next())
        {
            achiB += achi(x);
            aqB += aq(x);
        }
        parallel.sum<double>(achiB);
        parallel.sum<double>(aqB);
        achiB /= numpts3d;
        aqB /= numpts3d;
    }
}

void doAchiFieldUpdate(
    double dtau, int &N_as, double a_As, double fourpiG,
    metadata sim, double dx, double avgsource,
    cosmology cosmo, Field<Real> &achi, Field<Real> &aq,
    Field<Real> &scalar, Field<Real> &scalar2,
    Field<Real> &scalar3, Field<Real> &scalar4,
    Field<Real> &scalar5,
    Field<Real> &nT_cdm,
    int cycle, double achiB)
{
    if (strcmp(sim.As_solver, "leapfrog") == 0)
        doEvolveLF(dtau, N_as, a_As, fourpiG,
                   sim, dx, avgsource,
                   cosmo, &achi, &aq,
                   nT_cdm);
    else if (strcmp(sim.As_solver, "euler") == 0)
        doEvolveEC(dtau, N_as, a_As, fourpiG,
                   sim, dx, avgsource, cosmo, &achi, &aq, nT_cdm);
    else if (strcmp(sim.As_solver, "QSA") == 0)
        doAchiQSA(
            dtau, a_As, fourpiG, sim, dx,
            avgsource, cosmo, achi, aq, nT_cdm, achiB, cycle, 0, 0);
    else if (strcmp(sim.As_solver, "RK4") == 0)
        doEvolveRK4(
            dtau, N_as, a_As, fourpiG, sim, dx, avgsource, cosmo,
            &achi, &aq, nT_cdm, scalar, scalar2, scalar3, scalar4, scalar5);
}

void doFindBackgroundQuantities(
    Field<Real> &source, Field<Real> &T0i, Field<Real> &Sij,
    double dx, double boxsize, double a, const cosmology cosmo,
    Field<Real> &achi, Field<Real> &aq,
    const metadata sim, Field<Real> &nT_cdm,
    double &avgsource, double &eospar, double &eosparT, double &Vav,
    double &T00J, double &T00E, double &T00hom, double &T00As_hom,
    double &achiav, double &achi2av, double &aqav, double &aq2av,
    double &fracp, double &fracm, int cycle)
{
    Site xField(achi.lattice());
    double psif, phif, aphif, daphidt;
    double gradientachi_squared, Dx_achi_Dx_achi, Dy_achi_Dy_achi, Dz_achi_Dz_achi;
    double Dx_achi, Dy_achi, Dz_achi;
    double termtt, termii, termti, pot, lag;
    double phi0 = cosmo.mu_as / sqrt(cosmo.lambda_as);
    double phi02 = 0.5 * phi0 / dx / a;
    double boxsize2 = boxsize * boxsize;
    double rhocnorm = pow(2998., 2) / (3. * boxsize2);
    double a2 = a * a;
    double aetanorm = cosmo.aetanorm / (a * a2);
    double a3 = a * a2 * rhocnorm;
    double a4 = a2 * a2 * rhocnorm;
    double pressavgT = 0, pressavg = 0;
    double a6 = a2 * a4;
    double Acnf, Acnf4, src, temp;
    T00hom = 0.;
    T00As_hom = 0.;
    T00J = cosmo.Omega_m;
    T00E = 0;
    Vav = 0;
    achiav = 0;
    achi2av = 0;
    aqav = 0;
    aq2av = 0;
    fracp = 0;
    fracm = 0;
    for (xField.first(); xField.test(); xField.next())
    {
        Acnf = 1. + cosmo.dAcnf * pow(achi(xField), 2.);
        Acnf4 = pow(Acnf, 4);
        T00E += source(xField) * Acnf4;
        //------------------------
        //(D_i achi)^2
        //------------------------
        Dx_achi = phi02 * (achi(xField + 0) - achi(xField - 0));
        Dy_achi = phi02 * (achi(xField + 1) - achi(xField - 1));
        Dz_achi = phi02 * (achi(xField + 2) - achi(xField - 2));
        Dx_achi_Dx_achi = Dx_achi * Dx_achi;
        Dy_achi_Dy_achi = Dy_achi * Dy_achi;
        Dz_achi_Dz_achi = Dz_achi * Dz_achi;
        gradientachi_squared = Dx_achi_Dx_achi + Dy_achi_Dy_achi + Dz_achi_Dz_achi;
        //------------------------
        // Useful quantities
        aphif = phi0 * achi(xField);
        daphidt = phi0 * aq(xField) / a2 / a;
        //------------------------
        // STRESS TENSOR COMPONENTS
        //------------------------
        pot = boxsize2 * ((-cosmo.mu_as * cosmo.mu_as / 2.0 + (cosmo.kappa_as / 3.0 + cosmo.lambda_as / 4.0 * aphif) * aphif) * aphif * aphif);
        lag = -0.5 * (-daphidt * daphidt + gradientachi_squared) - pot;
        //
        src = source(xField) * Acnf4 - (-0.5 * daphidt * daphidt - pot - 0.5 * gradientachi_squared) * a3;
        // do misc
        T00hom += src;
        T00As_hom -= (-0.5 * daphidt * daphidt - pot - 0.5 * gradientachi_squared) * a3;
        pressavgT += (+Sij(xField, 0, 0) + Sij(xField, 1, 1) + Sij(xField, 2, 2)) / 6.;
        pressavg += (+Dx_achi_Dx_achi + Dy_achi_Dy_achi + Dz_achi_Dz_achi + 3. * lag) * a3 / 3.;
        Vav += pot * rhocnorm;
        achiav += achi(xField);
        achi2av += achi(xField) * achi(xField);
        aqav += aq(xField);
        aq2av += aq(xField) * aq(xField);
        temp = max(max(fabs(Dx_achi), fabs(Dy_achi)), fabs(Dz_achi));

        if (achi(xField) > 0.01)
        {
            fracp += 1;
        }
        if (achi(xField) < -0.01)
            fracm += 1;
        if (strcmp(sim.As_solver, "QSA") == 0 && sim.QSAq > 0 && a > sim.aMG)
            aq(xField) = achi(xField); // prepare
    }
    double numpts3d = pow((double)sim.numpts, 3.);

    parallel.sum<double>(T00hom);
    parallel.sum<double>(T00E);
    parallel.sum<double>(T00As_hom);
    parallel.sum<double>(pressavg);
    parallel.sum<double>(pressavgT);
    parallel.sum<double>(Vav);
    parallel.sum<double>(achiav);
    parallel.sum<double>(achi2av);
    parallel.sum<double>(aqav);
    parallel.sum<double>(aq2av);
    parallel.sum<double>(fracp);
    parallel.sum<double>(fracm);

    T00As_hom /= numpts3d;
    T00hom /= numpts3d;
    T00E /= numpts3d;
    Vav /= numpts3d;
    achiav /= numpts3d;
    achi2av /= numpts3d;
    aqav /= numpts3d;
    aq2av /= numpts3d;
    pressavg /= numpts3d;
    pressavgT /= numpts3d;
    fracp /= numpts3d;
    fracm /= numpts3d;
    eospar = pressavg / T00As_hom;
    eosparT = pressavgT / T00hom;

    avgsource = (T00hom + cosmo.Omega_rad / a + cosmo.Omega_Lambda * pow(a, 3.0));
    if (cycle % CYCLE_INFO_INTERVAL == 0)
    {
        COUT << " cycle " << cycle << ", background information: z = " << (1. / a) - 1. << ", average T00 = " << T00hom << ", background model = " << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) << endl;
    }
}

template <class FieldType>
void doConstructAchiSEtensor(
    Field<FieldType> &source, Field<FieldType> &T0i, Field<FieldType> &Sij,
    double dx, double boxsize, double a, const cosmology cosmo,
    Field<FieldType> &achi, Field<FieldType> &aq,
    const metadata sim, Field<Real> &nT_cdm, Field<Cplx> &scalarFT,
    PlanFFT<Cplx> &plan_Tcdm_As,
    double &avgsource, double &eospar, double &eosparT, double &Vav,
    double &T00J, double &T00E, double &T00hom, double &T00As_hom,
    double &achiav, double &achi2av, double &aqav, double &aq2av,
    double &fracp, double &fracm, int cycle, double fourpiG)
{
    Site xField(achi.lattice());
    double psif, phif, aphif, daphidt;
    double gradientachi_squared, Dx_achi_Dx_achi, Dy_achi_Dy_achi, Dz_achi_Dz_achi;
    double Dx_achi, Dy_achi, Dz_achi;
    double termtt, termii, termti, pot, lag;
    double phi0 = cosmo.mu_as / sqrt(cosmo.lambda_as);
    double phi02 = 0.5 * phi0 / dx / a;
    double boxsize2 = boxsize * boxsize;
    double rhocnorm = pow(2998., 2) / (3. * boxsize2);
    double a2 = a * a;
    double aetanorm = cosmo.aetanorm / (a * a2);
    double a3 = a * a2 * rhocnorm;
    double a4 = a2 * a2 * rhocnorm;
    double pressavgT = 0, pressavg = 0;
    double a6 = a2 * a4;
    double Acnf, Acnf4, temp;
    T00hom = 0.;
    T00As_hom = 0.;
    T00J = cosmo.Omega_m;
    T00E = 0;
    Vav = 0;
    achiav = 0;
    achi2av = 0;
    aqav = 0;
    aq2av = 0;
    fracp = 0;
    fracm = 0;
    for (xField.first(); xField.test(); xField.next())
    {
        Acnf = 1. + cosmo.dAcnf * pow(achi(xField), 2.);
        Acnf4 = pow(Acnf, 4);
        T00E += source(xField) * Acnf4;
        //------------------------
        //(D_i achi)^2
        //------------------------
        Dx_achi = phi02 * (achi(xField + 0) - achi(xField - 0));
        Dy_achi = phi02 * (achi(xField + 1) - achi(xField - 1));
        Dz_achi = phi02 * (achi(xField + 2) - achi(xField - 2));
        Dx_achi_Dx_achi = Dx_achi * Dx_achi;
        Dy_achi_Dy_achi = Dy_achi * Dy_achi;
        Dz_achi_Dz_achi = Dz_achi * Dz_achi;
        gradientachi_squared = Dx_achi_Dx_achi + Dy_achi_Dy_achi + Dz_achi_Dz_achi;
        //------------------------
        // Useful quantities
        aphif = phi0 * achi(xField);
        daphidt = phi0 * aq(xField) / a2 / a;
        //------------------------
        // STRESS TENSOR COMPONENTS
        //------------------------
        pot = boxsize2 * ((-cosmo.mu_as * cosmo.mu_as / 2.0 + (cosmo.kappa_as / 3.0 + cosmo.lambda_as / 4.0 * aphif) * aphif) * aphif * aphif);
        lag = -0.5 * (-daphidt * daphidt + gradientachi_squared) - pot;
        //
        source(xField) = source(xField) * Acnf4 - (-0.5 * daphidt * daphidt - pot - 0.5 * gradientachi_squared) * a3;
        if (sim.vector_flag == VECTOR_ELLIPTIC)
        {
            T0i(xField, 0) = Acnf4 * T0i(xField, 0) - Dx_achi * daphidt * a4;
            T0i(xField, 1) = Acnf4 * T0i(xField, 1) - Dy_achi * daphidt * a4;
            T0i(xField, 2) = Acnf4 * T0i(xField, 2) - Dz_achi * daphidt * a4;
        }
        Sij(xField, 0, 0) = Acnf4 * Sij(xField, 0, 0) + (Dx_achi_Dx_achi + lag) * a3;
        Sij(xField, 1, 1) = Acnf4 * Sij(xField, 1, 1) + (Dy_achi_Dy_achi + lag) * a3;
        Sij(xField, 2, 2) = Acnf4 * Sij(xField, 2, 2) + (Dz_achi_Dz_achi + lag) * a3;
        Sij(xField, 0, 1) = Acnf4 * Sij(xField, 0, 1) + Dy_achi * Dx_achi * a3;
        Sij(xField, 0, 2) = Acnf4 * Sij(xField, 0, 2) + Dz_achi * Dx_achi * a3;
        Sij(xField, 1, 2) = Acnf4 * Sij(xField, 1, 2) + Dz_achi * Dy_achi * a3;
        nT_cdm(xField) = source(xField) - (double)(sim.As_full_trace) * (Sij(xField, 0, 0) + Sij(xField, 1, 1) + Sij(xField, 2, 2));

        // do misc
        T00hom += source(xField);
        T00As_hom -= (-0.5 * daphidt * daphidt - pot - 0.5 * gradientachi_squared) * a3;
        pressavgT += (+Sij(xField, 0, 0) + Sij(xField, 1, 1) + Sij(xField, 2, 2)) / 3.;
        pressavg += (+Dx_achi_Dx_achi + Dy_achi_Dy_achi + Dz_achi_Dz_achi + 3. * lag) * a3 / 3.;
        Vav += pot * rhocnorm;
        achiav += achi(xField);
        achi2av += achi(xField) * achi(xField);
        aqav += aq(xField);
        aq2av += aq(xField) * aq(xField);
        if (achi(xField) > 0.01)
        {
            fracp += 1;
        }
        if (achi(xField) < -0.01)
            fracm += 1;
        if (strcmp(sim.As_solver, "QSA") == 0 && sim.QSAq > 0 && a > sim.aMG)
            aq(xField) = achi(xField); // prepare
    }
    source.updateHalo();
    if (sim.vector_flag == VECTOR_ELLIPTIC)
    {
        T0i.updateHalo();
    }
    Sij.updateHalo();
    double numpts3d = pow((double)sim.numpts, 3.);
    parallel.sum<double>(T00hom);
    parallel.sum<double>(T00E);
    parallel.sum<double>(T00As_hom);
    parallel.sum<double>(pressavg);
    parallel.sum<double>(pressavgT);
    parallel.sum<double>(Vav);
    parallel.sum<double>(achiav);
    parallel.sum<double>(achi2av);
    parallel.sum<double>(aqav);
    parallel.sum<double>(aq2av);
    parallel.sum<double>(fracp);
    parallel.sum<double>(fracm);

    T00As_hom /= numpts3d;
    T00hom /= numpts3d;
    T00E /= numpts3d;
    Vav /= numpts3d;
    achiav /= numpts3d;
    achi2av /= numpts3d;
    aqav /= numpts3d;
    aq2av /= numpts3d;
    pressavg /= numpts3d;
    pressavgT /= numpts3d;
    fracp /= numpts3d;
    fracm /= numpts3d;
    eospar = pressavg / T00As_hom;
    eosparT = pressavgT / T00hom;

    avgsource = (T00hom + cosmo.Omega_rad / a + cosmo.Omega_Lambda * pow(a, 3.0));

    if (cycle % CYCLE_INFO_INTERVAL == 0)
    {
        COUT << " cycle " << cycle << ", background information: z = " << (1. / a) - 1. << ", average T00 = " << T00hom << ", background model = " << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) << endl;
    }
}

void doUpdateAchiBackground(double &achiB, double &aqB, double a, const double fourpiG,
                            const double dtau, const cosmology cosmo, const metadata sim, double avgsource)
{
    int i;
    double tmpa = a;
    for (i = 0; i < sim.nAs_numsteps_dynbck; i++)
        rungekutta4achi(achiB, aqB, tmpa, fourpiG, dtau / sim.nAs_numsteps_dynbck, cosmo, sim, avgsource);
}

