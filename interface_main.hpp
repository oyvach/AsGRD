#ifndef INTERFACE_HEADER
#define INTERFACE_HEADER
#endif

#ifndef PI
#define PI 3.1415
#endif

#ifndef Cplx
#define Cplx Imag
#endif
using namespace std;
using namespace LATfield2;

void doReadCommandLine(char **argv, int argc, char *&settingsfile, char *&precisionfile, int &n, int &m, int &io_size, int &io_group_size)
{
    int i;
    for (i = 1; i < argc; i++)
    {
        if (argv[i][0] != '-')
            continue;
        switch (argv[i][1])
        {
        case 's':
            settingsfile = argv[++i]; // settings file name
            break;
        case 'n':
            n = atoi(argv[++i]); // size of the dim 1 of the processor grid
            break;
        case 'm':
            m = atoi(argv[++i]); // size of the dim 2 of the processor grid
            break;
        case 'p':
#ifndef HAVE_CLASS
            cout << "HAVE_CLASS needs to be set at compilation to use CLASS precision files" << endl;
            exit(-100);
#endif
            precisionfile = argv[++i];
            break;
        case 'i':
#ifndef EXTERNAL_IO
            cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server" << endl;
            exit(-1000);
#endif
            io_size = atoi(argv[++i]);
            break;
        case 'g':
#ifndef EXTERNAL_IO
            cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server" << endl;
            exit(-1000);
#endif
            io_group_size = atoi(argv[++i]);
        }
    }
}

void doTitleAndSettings(int n, int m, char *settingsfile,
                        int &numparam, double &start_time,
                        char *filename, cosmology &cosmo,
                        metadata &sim, icsettings &ic, parameter *params)
{
    COUT << COLORTEXT_WHITE << endl;

    COUT
        << "         /\\                  .-..-.          .-.  \n   "
        << "     _  / |           .;;.`-'  (_) )-.      (_) )-.    \n"
        << "    (  /  |  .    .  ;; (_;      .:   \\       .:   \\   \n"
        << "     `/.__|_.'  .'; ;;          .::.   )     .:'    \\  \n"
        << " .:' /    |   .' .';;    `;;' .-:. `:-'    .-:.      ) \n"
        << "(__.'     `-''     `;.___.'  (_/     `:._.(_/  `----'  \n"
        << "                  .-.    .  _   *     _   .\n"
        << "      *          /   \\     ((       _/ \\       *    .\n"
        << "    _    .   .--'\\/\\_ \\     `      /    \\  *    ___\n"
        << "*  / \\_    _/ ^      \\/\'__        /\\/\\  /\\  __/   \\ *\n"
        << "  /    \\  /    .'   _/  /  \\  *' /    \\/  \\/ .`'\\_/\\   .\n"
        << " /\\/\\  /\\/ :' __  ^/  ^/    `--./.'  ^  `-.\\ _    _:\\ _\n"
        << COLORTEXT_RESET << endl;
    COUT << "version 1.0         running on " << n * m << " cores." << endl;

    if (settingsfile == NULL)
    {
        COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no settings file specified!" << endl;
        parallel.abortForce();
    }

    COUT << " initializing..." << endl;

    start_time = MPI_Wtime();

    numparam = loadParameterFile(settingsfile, params);

    int usedparams = parseMetadata(params, numparam, sim, cosmo, ic);

    COUT << " parsing of settings file completed. " << numparam << " parameters found, " << usedparams << " were used." << endl;

    snprintf(filename, SIZFILENAME, "%s%s_settings_used.ini", sim.output_path, sim.basename_generic);
    saveParameterFile(filename, params, numparam);

    free(params);
}

void doICs(cosmology &cosmo, double &a,
           metadata &sim, icsettings &ic,
           Field<Real> &achi, Field<Real> &aq,
           double avgsource, double fourpiG,
           double tau, double dtau, double dtau_old,
           string h5filename,
           Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_cdm,
           Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_b,
           Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_ncdm,
           double *maxvel, Field<Real> &phi, Field<Real> &chi, Field<Real> &Bi,
           Field<Real> &hij, Field<Real> &source, Field<Cplx> &scalarFT,
           Field<Cplx> &BiFT, Field<Cplx> &SijFT, PlanFFT<Cplx> &plan_phi,
           PlanFFT<Cplx> &plan_chi, PlanFFT<Cplx> &plan_Bi, PlanFFT<Cplx> &plan_Sij,
           PlanFFT<Cplx> &plan_source,
           parameter *params, int &numparam, int &numspecies, int cycle, int snapcount,
           int pkcount, int restartcount, set<long> *IDbacklog,
           PlanFFT<Cplx> &plan_achi, PlanFFT<Cplx> &plan_aq, Field<Real> &nT_cdm
#ifdef GW
           ,
           Field<Cplx> &hijFT, Field<Cplx> &hijprimeFT, PlanFFT<Cplx> &plan_hij, PlanFFT<Cplx> &plan_hijprime
#endif
)
{
    int i;
    rKSite kFT(scalarFT.lattice());

#ifdef GW // initialise properly at cycle=2
    if (sim.GWhighkAdiabatic > 0)
        evolveFTtensor = &evolveFTtensorAdiabatic;
    else
        evolveFTtensor = &evolveFTtensorZero;
    if (sim.GWreadFields == 0)
    {
        for (kFT.first(); kFT.test(); kFT.next())
        {
            for (i = 0; i < hijprimeFT.components(); i++)
                hijprimeFT(kFT, i) = Cplx(0., 0.);
        }
    }
    else // initialise from reading fields
    {
        if (strcmp(sim.hijfile, "no") == 0 || strcmp(sim.hijprimefile, "no") == 0)
        {
            COUT << COLORTEXT_RED << "Have not provided files to read hij fields from" << COLORTEXT_RESET << endl;
            return;
        }
        else
        {
            COUT << "Reading hij/hijprime from files " << sim.hijfile
                 << " and " << sim.hijprimefile << endl;
            string filename = sim.hijfile;
            hij.loadHDF5(filename);
            plan_hij.execute(FFT_FORWARD);
            filename = sim.hijprimefile;
            hij.loadHDF5(filename);
            plan_hijprime.execute(FFT_FORWARD);
            // normalise (unormalised FFTW convention)
            double numpts3d = (double)pow(sim.numpts, 3);
            for (kFT.first(); kFT.test(); kFT.next())
            {
                for (i = 0; i < hijprimeFT.components(); i++)
                {
                    hijprimeFT(kFT, i) /= numpts3d;
                    hijFT(kFT, i) /= numpts3d;
                }
            }
        }
    }
#endif

    if (ic.generator == ICGEN_BASIC)
        generateIC_basic(sim, ic, cosmo, avgsource, fourpiG, pcls_cdm, pcls_b, pcls_ncdm, maxvel, &phi,
                         &chi, &Bi, &source, &hij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi,
                         &plan_source, &plan_Sij, params, numparam); // generates ICs on the fly
    else if (ic.generator == ICGEN_READ_FROM_DISK)
        readIC(sim, ic, cosmo, avgsource, fourpiG, a, tau, dtau, dtau_old, pcls_cdm, pcls_b, pcls_ncdm, maxvel,
               &phi, &chi, &Bi, &source, &hij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source,
               &plan_Sij, cycle, snapcount, pkcount, restartcount, IDbacklog);
#ifdef ICGEN_PREVOLUTION
    else if (ic.generator == ICGEN_PREVOLUTION)
        generateIC_prevolution(sim, ic, cosmo, avgsource, fourpiG, a, tau, dtau, dtau_old, pcls_cdm, pcls_b, pcls_ncdm,
                               maxvel, &phi, &chi, &Bi, &source, &hij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi,
                               &plan_Bi, &plan_source, &plan_Sij, params, numparam);
#endif
#ifdef ICGEN_FALCONIC
    else if (ic.generator == ICGEN_FALCONIC)
        maxvel[0] = generateIC_FalconIC(sim, ic, cosmo, fourpiG, dtau, pcls_cdm, pcls_ncdm, maxvel + 1, &phi, &source,
                                        &chi, &Bi, &source, &hij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_source,
                                        &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif
    else
    {
        COUT << " error: IC generator not implemented!" << endl;
        parallel.abortForce();
    }

    if (sim.baryon_flag > 1)
    {
        COUT << " error: baryon_flag > 1 after IC generation, something went wrong in IC generator!" << endl;
        parallel.abortForce();
    }

    numspecies = 1 + sim.baryon_flag + cosmo.num_ncdm;
    parallel.max<double>(maxvel, numspecies);

    if (sim.gr_flag > 0)
    {
        for (i = 0; i < numspecies; i++)
            maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
    }
}

void doConstructSEtensor(
    const cosmology cosmo, double a, const metadata sim,
    Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_cdm,
    Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_b,
    Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_ncdm,
    Field<Real> &phi, Field<Real> &phiJ, Field<Real> &chi,
    Field<Real> &achi, Field<Real> &T0i,
    Field<Real> &source, Field<Real> &Sij)
{
    int i;
    double tmp;
    double fac = 0.5 * pow(cosmo.mu_as / (sqrt(cosmo.lambda_as) * cosmo.M_as), 2) * (double)sim.As_use_phiJ_SE;
    double fac2;
    Site x(phi.lattice());

    // construct stress-energy tensor
    projection_init(&source);
    if (sim.gr_flag > 0)
    {
        if (a > sim.aMG)
            fac2 = 1.;
        else
            fac2 = 0.;
        // convert to Jordan frame potential
        for (x.first(); x.test(); x.next())
        {
            phiJ(x) = phi(x) - fac2 * fac * pow(achi(x), 2);
        }
        phiJ.updateHalo();
        projection_T00_project(pcls_cdm, &source, a, &phiJ);
        if (sim.baryon_flag)
            projection_T00_project(pcls_b, &source, a, &phiJ);
        for (i = 0; i < cosmo.num_ncdm; i++)
        {
            if (a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1 + sim.baryon_flag + i] > 0)
                projection_T00_project(pcls_ncdm + i, &source, a, &phiJ);
            else if (sim.radiation_flag == 0 || (a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1 + sim.baryon_flag + i] == 0))
            {
                tmp = bg_ncdm(a, cosmo, i);
                for (x.first(); x.test(); x.next())
                    source(x) += tmp;
            }
        }
    }
    else
    {
        scalarProjectionCIC_project(pcls_cdm, &source);
        if (sim.baryon_flag)
            scalarProjectionCIC_project(pcls_b, &source);
        for (i = 0; i < cosmo.num_ncdm; i++)
        {
            if (a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1 + sim.baryon_flag + i] > 0)
                scalarProjectionCIC_project(pcls_ncdm + i, &source);
        }
    }
    projection_T00_comm(&source);

    if (sim.vector_flag == VECTOR_ELLIPTIC)
    {
        projection_init(&T0i);
        projection_T0i_project(pcls_cdm, &T0i, &phiJ);
        if (sim.baryon_flag)
            projection_T0i_project(pcls_b, &T0i, &phiJ);
        for (i = 0; i < cosmo.num_ncdm; i++)
        {
            if (a >= 1. / (sim.z_switch_Bncdm[i] + 1.) && sim.numpcl[1 + sim.baryon_flag + i] > 0)
                projection_T0i_project(pcls_ncdm + i, &T0i, &phiJ);
        }
        projection_T0i_comm(&T0i);
    }

    projection_init(&Sij);
    projection_Tij_project(pcls_cdm, &Sij, a, &phiJ);
    if (sim.baryon_flag)
        projection_Tij_project(pcls_b, &Sij, a, &phiJ);
    if (a >= 1. / (sim.z_switch_linearchi + 1.))
    {
        for (i = 0; i < cosmo.num_ncdm; i++)
        {
            if (sim.numpcl[1 + sim.baryon_flag + i] > 0)
                projection_Tij_project(pcls_ncdm + i, &Sij, a, &phiJ);
        }
    }
    projection_Tij_comm(&Sij);
}

void doBackgroundOutput(
    double phik0, const cosmology cosmo, double a, const metadata sim,
    char *filename, FILE *outfile, int cycle, int snapcount,
    double tau, double fourpiG, double avgsource, double T00hom, int N_as,
    double T00As_hom, double eospar, double eosparT, double Vav, double T00J, double T00E,
    double achiav, double achi2av, double aqav, double aq2av,
    double achiB, double aqB, double tau2,
    double fracp, double fracm, double maxvelCDM,
    int pkcount
#ifdef GW
    ,
    double avhijprimenorm
#endif
)
{
    int snapoutb = 0, pkoutb = 0;
    string tmp;

    // record some background data
    if (parallel.isRoot())
    {
        snprintf(filename, SIZFILENAME, "%s%s_background.dat", sim.output_path, sim.basename_generic);
        // modification -- fix annoying issue:
        if (cycle == 0)
            outfile = fopen(filename, "w"); // make output files new
        else
            outfile = fopen(filename, "a");
        if (outfile == NULL)
        {
            cout << " error opening file for background output!" << endl;
        }
        else // modification -- add columns
        {
            tmp = "# background statistics\n# cycle   tau/boxsize   a   conformal H/H0   ";
            tmp += "phi(k=0)   T00(k=0)   n_as   eta   T00As(k=0)   eospar   eosparT   ";
            tmp += "Vav   T00J   T00E   achiav   achi2av   aqav   aq2av   achib   aqb   tau2   avgsource   ";
            tmp += "fracp   fracm   ";
            tmp += "maxvelCDM   snapoutbool   pkoutbool";
#ifdef GW
            tmp += +"   avhijprimenorm";
#endif
            tmp += "\n";
            if (cycle == 0)
            {
                fprintf(outfile, "%s", tmp.c_str());
            }
            if (snapcount < sim.num_snapshot && 1. / a <= sim.z_snapshot[snapcount] + 1.0000001)
                snapoutb = 1;
            if (pkcount < sim.num_pk && 1. / a <= sim.z_pk[pkcount] + 1.0000001)
                pkoutb = 1;
            tmp = " %6d   %e   %e   %e   %e   %e   %d   %e   %e   %e   %e   %e   %e   %e   %e   %e   %e   %e   ";
            tmp += "%e   %e   %e   %e   %e   %e   %e   %d   %d";
#ifdef GW
            tmp += "   %e";
#endif
            tmp += "\n";
            fprintf(outfile, tmp.c_str(),
                    cycle, tau, a, Hconf(a, fourpiG, cosmo, avgsource) / sqrt(2. * fourpiG / 3.),
                    phik0, T00hom, N_as,
                    T00hom / (a * a * a) / pow(cosmo.M_as * cosmo.mu_as, 2.0), T00As_hom, eospar,
                    eosparT, Vav, T00J, T00E, achiav, achi2av, aqav, aq2av, achiB, aqB, tau2, avgsource,
                    fracp, fracm, maxvelCDM, snapoutb, pkoutb
#ifdef GW
                    ,
                    avhijprimenorm
#endif
            );
            fclose(outfile);
        }
    }
}

void doSolveEinsteinEquation(
    const cosmology cosmo, double a, const metadata sim,
    Field<Real> &phi, Field<Real> &chi, Field<Real> &Bi,
    Field<Real> &Sij, Field<Real> &source, Field<Real> &source2, double T00hom, double avgsource,
    Field<Cplx> &scalarFT, Field<Cplx> &BiFT, Field<Cplx> &SijFT,
    PlanFFT<Cplx> &plan_source, PlanFFT<Cplx> &plan_source2, PlanFFT<Cplx> &plan_phi, PlanFFT<Cplx> &plan_chi,
    PlanFFT<Cplx> &plan_Bi, PlanFFT<Cplx> &plan_Sij,
    double dx, double fourpiG, double dtau_old, int cycle, double &phik0, int &done_hij
#ifdef GW
    ,
    Field<Cplx> &hijFT, Field<Cplx> &hijprimeFT, double dtau,
    PlanFFT<Cplx> &plan_hij, PlanFFT<Cplx> &plan_hijprime, double &avhijprimenorm
#endif
)
{
    if (cycle > 0)
    {
        if (sim.gr_flag > 0)
        {

            if (dtau_old > 0.)
            {
                // use source2 to not have to redo projection at pk output (since we are keeping scalar2 anyway)
                prepareFTsource<Real>(phi, chi, source, T00hom, // cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) + cosmo.Omega_sym,
                                      source2, 3. * Hconf(a, fourpiG, cosmo, avgsource) * dx * dx / dtau_old,
                                      fourpiG * dx * dx / a,
                                      3. * Hconf(a, fourpiG, cosmo, avgsource) * Hconf(a, fourpiG, cosmo, avgsource) * dx * dx); // prepare nonlinear source for phi update

                plan_source2.execute(FFT_FORWARD);                                                                               // go to k-space
                solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx), 3. * Hconf(a, fourpiG, cosmo, avgsource) / dtau_old); // phi update (k-space)
                plan_phi.execute(FFT_BACKWARD);
            }
        }
        else
        {
            plan_source.execute(FFT_FORWARD);                        // Newton: directly go to k-space
            solveModifiedPoissonFT(scalarFT, scalarFT, fourpiG / a); // Newton: phi update (k-space)
            plan_phi.execute(FFT_BACKWARD);
        }
        phi.updateHalo();
    }
    rKSite kFT(scalarFT.lattice());
    if (kFT.setCoord(0, 0, 0))
    { // prepare for background output
        phik0 = scalarFT(kFT).real();
    }
    prepareFTsource<Real>(phi, Sij, Sij, 2. * fourpiG * dx * dx / a); // prepare nonlinear source for additional equations
    plan_Sij.execute(FFT_FORWARD);
    if (cycle > 0)
    {
        projectFTscalar(SijFT, scalarFT); // construct chi by scalar projection (k-space)
        plan_chi.execute(FFT_BACKWARD);

        chi.updateHalo();

        if (sim.gr_flag > 0)
        {
            if (sim.vector_flag == VECTOR_ELLIPTIC)
            {
                plan_Bi.execute(FFT_FORWARD);
                projectFTvector(BiFT, BiFT, fourpiG * dx * dx); // solve B using elliptic constraint (k-space)
            }
            else
                evolveFTvector(SijFT, BiFT, a * a * dtau_old); // evolve B using vector projection (k-space)

            plan_Bi.execute(FFT_BACKWARD);
            Bi.updateHalo();
        }
#ifdef GW
        avhijprimenorm = 0.;
        if (cycle == sim.GWinitialiseCycle - 1 && sim.GWreadFields == 0 && sim.GWstartAdiabatic > 0)
        {
            projectFTtensor(SijFT, hijFT);
            rKSite kFT(hijprimeFT.lattice());
            // prepare for initialisation from adiabatic condition
            for (kFT.first(); kFT.test(); kFT.next())
            {
                for (int i = 0; i < hijprimeFT.components(); i++)
                    hijprimeFT(kFT, i) = hijFT(kFT, i);
            }
        }
        else if (cycle < sim.GWinitialiseCycle && sim.GWreadFields == 0 && sim.GWstartAdiabatic == 0)
        {
            projectFTtensor(SijFT, hijFT);
            rKSite kFT(hijprimeFT.lattice());
            // prepare for initialisation from adiabatic condition
            for (kFT.first(); kFT.test(); kFT.next())
            {
                for (int i = 0; i < hijprimeFT.components(); i++)
                    hijprimeFT(kFT, i) = Cplx(0., 0.);
            }
        }
        else if (cycle == sim.GWinitialiseCycle && sim.GWreadFields == 0 && sim.GWstartAdiabatic > 0)
        {
            projectFTtensor(SijFT, hijFT);
            rKSite kFT(hijprimeFT.lattice());
            for (kFT.first(); kFT.test(); kFT.next())
            {
                for (int i = 0; i < hijprimeFT.components(); i++)
                    hijprimeFT(kFT, i) = (hijFT(kFT, i) - hijprimeFT(kFT, i)) / dtau_old;
            }
        }
        else if (cycle > 0 && sim.GWreadFields > 0)
        {
            // set special dtau_old
            evolveFTtensor(SijFT, hijFT, hijprimeFT, Hconf(a, fourpiG, cosmo, avgsource), dtau, dtau);
        }
        else if (cycle > 0)
        {
            evolveFTtensor(SijFT, hijFT, hijprimeFT, Hconf(a, fourpiG, cosmo, avgsource), dtau, dtau_old);
        }
        plan_hijprime.execute(FFT_BACKWARD);
        Site x(Sij.lattice());
        for (x.first(); x.test(); x.next())
        {
            avhijprimenorm += +Sij(x, 0, 0) * Sij(x, 0, 0) + Sij(x, 1, 1) * Sij(x, 1, 1)       /* */
                              + Sij(x, 2, 2) * Sij(x, 2, 2) + 2. * Sij(x, 0, 1) * Sij(x, 0, 1) /* */
                              + 2. * Sij(x, 0, 2) * Sij(x, 0, 2) + 2. * Sij(x, 1, 2) * Sij(x, 1, 2);
        }
        parallel.sum<double>(avhijprimenorm);
        avhijprimenorm /= pow((double)sim.numpts, 3.);
        // remember to plan_hij.execute(FFT_BACKWARD before use in real space)
#else
        // don't do anything since hij is only used for output
        done_hij = 0;
#endif
    }
}

void doWriteSpectra(double dtau, int cycle, double tau,
                    metadata &sim, cosmology &cosmo, double avgsource, const double fourpiG, const double a, int &pkcount,
                    Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_cdm,
                    Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_b,
                    Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_ncdm,
                    Field<Real> *phi, Field<Real> *achi, Field<Real> *aq, Field<Real> *chi,
                    Field<Real> *Bi, Field<Real> *hij, Field<Real> *source,
                    Field<Real> *scalar,
                    Field<Cplx> *scalarFT, Field<Cplx> *BiFT, Field<Cplx> *SijFT,
                    PlanFFT<Cplx> *plan_phi, PlanFFT<Cplx> *plan_achi, PlanFFT<Cplx> *plan_aq,
                    PlanFFT<Cplx> *plan_chi, PlanFFT<Cplx> *plan_Bi,
                    PlanFFT<Cplx> *plan_source, PlanFFT<Cplx> *plan_Sij,
                    PlanFFT<Cplx> *plan_scalar, int &done_hij
#ifdef GW
                    ,
                    Field<Cplx> *hijFT, Field<Cplx> *hijprimeFT,
                    PlanFFT<Cplx> *plan_hij, PlanFFT<Cplx> *plan_hijprime
#endif
)
{
    COUT << COLORTEXT_CYAN << " writing power spectra" << COLORTEXT_RESET << " at z = " << ((1. / a) - 1.) << " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

    writeSpectra(sim, cosmo, avgsource, fourpiG, a, pkcount,
                 pcls_cdm, pcls_b, pcls_ncdm, phi, achi, aq, chi, Bi, hij, source, scalar,
                 scalarFT, BiFT, SijFT, plan_phi, plan_achi, plan_aq, plan_chi, plan_Bi,
                 plan_source, plan_Sij, plan_scalar, done_hij
#ifdef GW
                 ,
                 hijFT, hijprimeFT,
                 plan_hij, plan_hijprime
#endif
    );

    pkcount++;
}
// debug
/*To use:
Field<Real> *fields[2];
fields[0] = update_cdm_fields[0];
fields[1] = update_cdm_fields[1];
printMaxMinFields(fields, 2);
*/
void printMaxMinFields(Field<Real> **fields, int length)
{
    Site x(fields[0]->lattice());
    double vmin, vmax, vval;
    for (int i = 0; i < length; i++)
    {
        vmin = 1e10;
        vmax = 0;
        for (x.first(); x.test(); x.next())
        {
            vval = abs((*fields[i])(x));
            if (vval > vmax)
                vmax = vval;
            if (vval < vmin)
                vmin = vval;
        }
        parallel.max<double>(vmax);
        parallel.min<double>(vmin);
        COUT << "Field " << i << " has min/max " << vmin << " / " << vmax << endl;
    }
}
void doUpdateParticlesAndBackground(
    const metadata sim, const cosmology cosmo, double avgsource, double fourpiG, double &a,
    Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> &pcls_cdm,
    Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> &pcls_b,
    Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_ncdm,
    double dtau, double dtau_old, double dx, double *maxvel, double *f_params,
    icsettings ic, Field<Real> **update_cdm_fields, Field<Real> **update_b_fields,
    Field<Real> **update_ncdm_fields,
    int nfields, int numspecies, int *numsteps_ncdm)
{
    int i, j;
    double tmp;
    // compute number of step subdivisions for ncdm particle updates
    for (i = 0; i < cosmo.num_ncdm; i++)
    {
        if (dtau * maxvel[i + 1 + sim.baryon_flag] > dx * sim.movelimit)
            numsteps_ncdm[i] = (int)ceil(dtau * maxvel[i + 1 + sim.baryon_flag] / dx / sim.movelimit);
        else
            numsteps_ncdm[i] = 1;
    }
    for (i = 0; i < cosmo.num_ncdm; i++) // non-cold DM particle update
    {
        if (sim.numpcl[1 + sim.baryon_flag + i] == 0)
            continue;

        tmp = a;
        for (j = 0; j < numsteps_ncdm[i]; j++)
        {
            f_params[0] = tmp;
            f_params[1] = tmp * tmp * (double)sim.numpts;
            if (sim.As_5th_force * sim.fifth > 0) // modification
                maxvel[i + 1 + sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q_5th, (dtau + dtau_old) / 2. / (double)numsteps_ncdm[i], update_ncdm_fields, nfields, f_params);
            else if (sim.As_5th_force_Newton * sim.fifth > 0)
                maxvel[i + 1 + sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q_Newton_5th, (dtau + dtau_old) / 2. / (double)numsteps_ncdm[i], update_ncdm_fields, nfields, f_params);
            else if (sim.gr_flag > 0)
                maxvel[i + 1 + sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q, (dtau + dtau_old) / 2. / (double)numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
            else
                maxvel[i + 1 + sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q_Newton, (dtau + dtau_old) / 2. / (double)numsteps_ncdm[i], update_ncdm_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);

            rungekutta4bg(tmp, fourpiG, 0.5 * dtau / (double)numsteps_ncdm[i], cosmo, avgsource);
            f_params[0] = tmp;
            f_params[1] = tmp * tmp * (double)sim.numpts;
            if (sim.As_5th_force * sim.fifth > 0) // modification
                pcls_ncdm[i].moveParticles(update_pos_5th, dtau / (double)numsteps_ncdm[i], update_ncdm_fields, nfields, f_params);
            else if (sim.As_5th_force_Newton * sim.fifth > 0)
                pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / (double)numsteps_ncdm[i], NULL, 0, f_params);
            else if (sim.gr_flag > 0)
                pcls_ncdm[i].moveParticles(update_pos, dtau / (double)numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
            else
                pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / (double)numsteps_ncdm[i], NULL, 0, f_params);

            rungekutta4bg(tmp, fourpiG, 0.5 * dtau / (double)numsteps_ncdm[i], cosmo, avgsource);
        }
    }

    // cdm and baryon particle update
    f_params[0] = a;
    f_params[1] = a * a * (double)sim.numpts;
    f_params[3] = Hconf(a, fourpiG, cosmo, avgsource);
    if (sim.As_5th_force * sim.fifth > 0) // modification
    {
        maxvel[0] = pcls_cdm.updateVel(update_q_5th, (dtau + dtau_old) / 2., update_cdm_fields, nfields, f_params);
        if (sim.baryon_flag)
            maxvel[1] = pcls_b.updateVel(update_q_5th, (dtau + dtau_old) / 2., update_b_fields, nfields, f_params);
    }
    else if (sim.As_5th_force_Newton * sim.fifth > 0)
    {
        maxvel[0] = pcls_cdm.updateVel(update_q_Newton_5th, (dtau + dtau_old) / 2., update_cdm_fields, nfields, f_params);
        if (sim.baryon_flag)
        {
            maxvel[1] = pcls_b.updateVel(update_q_Newton_5th, (dtau + dtau_old) / 2., update_b_fields, nfields, f_params);
        }
    }
    else if (sim.gr_flag > 0)
    {

        maxvel[0] = pcls_cdm.updateVel(update_q, (dtau + dtau_old) / 2., update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
        if (sim.baryon_flag)
            maxvel[1] = pcls_b.updateVel(update_q, (dtau + dtau_old) / 2., update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
    }
    else
    {
        maxvel[0] = pcls_cdm.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_cdm_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
        if (sim.baryon_flag)
            maxvel[1] = pcls_b.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_b_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
    }

    rungekutta4bg(a, fourpiG, 0.5 * dtau, cosmo, avgsource); // evolve background by half a time step

    f_params[0] = a;
    f_params[1] = a * a * (double)sim.numpts;
    f_params[3] = Hconf(a, fourpiG, cosmo, avgsource);
    if (sim.As_5th_force * sim.fifth > 0) // modification
    {
        pcls_cdm.moveParticles(update_pos_5th, dtau, update_cdm_fields, nfields, f_params);
        if (sim.baryon_flag)
            pcls_b.moveParticles(update_pos_5th, dtau, update_b_fields, nfields, f_params);
    }
    else if (sim.As_5th_force_Newton * sim.fifth > 0)
    {
        pcls_cdm.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
        if (sim.baryon_flag)
            pcls_b.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
    }
    else if (sim.gr_flag > 0)
    {
        pcls_cdm.moveParticles(update_pos, dtau, update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
        if (sim.baryon_flag)
            pcls_b.moveParticles(update_pos, dtau, update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
    }
    else
    {
        pcls_cdm.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
        if (sim.baryon_flag)
            pcls_b.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
    }

    rungekutta4bg(a, fourpiG, 0.5 * dtau, cosmo, avgsource); // evolve background by half a time step
    parallel.max<double>(maxvel, numspecies);

    if (sim.gr_flag > 0)
    {
        for (i = 0; i < numspecies; i++)
            maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
    }
    // done particle update
}

bool startedMG = false;
void doSetNewTimestep(double &dtau, double dx, double a, int &N_as, metadata &sim, const cosmology cosmo,
                      double avgsource, double fourpiG, int pkcount, int snapcount, int snapcount2, bool &writeSnap, bool &writeSnap2, bool &writeSpec,
                      int cycle, double maxvelcdm)
{
    double zmax;

    double Cf;
    if (!startedMG && a > sim.aMG && sim.fifth > 0)
    {
        COUT << "Adjusting Courant factor for MG\n";
        Cf = sim.Cf_MG;
        N_as = max(int(Cf / sim.As_Courant + 0.5), 1);
        startedMG = true;
    }
    else if (!startedMG)
    {
        Cf = sim.Cf;
    }
    else if (a > sim.aMG)
    {
        Cf = sim.Cf_MG;
    }
    double tol = 1e-1 * Cf * dx;
    // consider if now is output
    writeSpec = (pkcount < sim.num_pk && 1. / a - 1 - sim.z_pk[pkcount] <= tol);
    writeSnap = (snapcount < sim.num_snapshot && 1. / a - 1 - sim.z_snapshot[snapcount] <= tol);
    writeSnap2 = (snapcount2 < sim.num_snapshot2 && 1. / a - 1 - sim.z_snapshot2[snapcount2] <= tol);
    // standard choice of dtau
    dtau = Cf * dx;

    if (!writeSnap && (sim.writefirstsnap > 0))
    {
        sim.num_snapshot++;
        writeSnap = true;
        sim.writefirstsnap = 0;
        // fix sim.z_snapcount
        double tmp[sim.num_snapshot];
        tmp[snapcount] = sim.z_snapshot[snapcount];
        sim.z_snapshot[snapcount] = 1. / a - 1.;
        for (int i = snapcount + 1; i < sim.num_snapshot; i++)
        {
            tmp[i] = sim.z_snapshot[i];
            sim.z_snapshot[i] = tmp[i - 1];
        }
    }
    if (!writeSpec && (sim.writefirstspec > 0))
    {
        sim.num_pk++;
        writeSpec = true;
        sim.writefirstspec = 0;
        // fix sim.z_snapcount
        double tmp[sim.num_pk];
        tmp[pkcount] = sim.z_pk[pkcount];
        sim.z_pk[pkcount] = 1. / a - 1.;
        for (int i = pkcount + 1; i < sim.num_pk; i++)
        {
            tmp[i] = sim.z_pk[i];
            sim.z_pk[i] = tmp[i - 1];
        }
    }

    if (writeSpec)
        pkcount++;
    if (writeSnap)
        snapcount++;
    if (writeSnap2)
        snapcount2++;

    double anext = a;
    rungekutta4bg(anext, fourpiG, dtau, cosmo, avgsource);
    // consider if the iteration after next dtau is output and adjust dtau
    bool specnext = (pkcount < sim.num_pk && 1. / anext - 1 - sim.z_pk[pkcount] <= tol);
    bool snapnext = (snapcount < sim.num_snapshot && 1. / anext - 1 - sim.z_snapshot[snapcount] <= tol);
    bool snapnext2 = (snapcount2 < sim.num_snapshot2 && 1. / anext - 1 - sim.z_snapshot2[snapcount2] <= tol);

    if (snapnext || snapnext2 || specnext || anext > 1)
    {
        if (anext < 1)
            zmax = max(sim.z_pk[pkcount], max(sim.z_snapshot[snapcount], sim.z_snapshot2[snapcount2]));
        else
            zmax = 0;
        double anext = 1. / (1 + zmax);
        rungekutta4dtau(dtau, a, fourpiG, anext - a, cosmo, avgsource);
        COUT << COLORTEXT_YELLOW << "z = " << 1 / a - 1. << COLORTEXT_RESET << ": Adding one iteration for redshift precision outputs.\n";
    }
}

void doCoutCycleInfo(int cycle, double dtau, double dx, double a,
                     double fourpiG, const cosmology cosmo, const metadata sim,
                     int *numsteps_ncdm, double *maxvel, double avgsource)
{
    int i;
    if (cycle % CYCLE_INFO_INTERVAL == 0)
    {
        COUT << " cycle " << cycle << ", time integration information: max |v| = " << maxvel[0] << " (cdm Courant factor = " << maxvel[0] * dtau / dx;
        if (sim.baryon_flag)
        {
            COUT << "), baryon max |v| = " << maxvel[1] << " (Courant factor = " << maxvel[1] * dtau / dx;
        }

        COUT << "), time step / Hubble time = " << Hconf(a, fourpiG, cosmo, avgsource) * dtau;

        for (i = 0; i < cosmo.num_ncdm; i++)
        {
            if (i == 0)
            {
                COUT << endl
                     << " time step subdivision for ncdm species: ";
            }
            COUT << numsteps_ncdm[i] << " (max |v| = " << maxvel[i + 1 + sim.baryon_flag] << ")";
            if (i < cosmo.num_ncdm - 1)
            {
                COUT << ", ";
            }
        }

        COUT << endl;
    }
}
