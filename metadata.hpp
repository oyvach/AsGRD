//////////////////////////
// metadata.hpp
//////////////////////////
//
// Constants and metadata structures
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: April 2019
//
//////////////////////////

#ifndef METADATA_HEADER
#define METADATA_HEADER

#define GEVOLUTION_VERSION 1.2

#ifndef MAX_OUTPUTS
#define MAX_OUTPUTS 10000 // for animations
#endif

#ifndef PARAM_MAX_LENGTH
#define PARAM_MAX_LENGTH 512
#endif

#ifndef PARAM_MAX_LINESIZE
#define PARAM_MAX_LINESIZE 2048
#endif

#ifndef MAX_INTERSECTS
#define MAX_INTERSECTS 12
#endif

#ifndef LIGHTCONE_IDCHECK_ZONE
#define LIGHTCONE_IDCHECK_ZONE 0.05
#endif

#define LIGHTCONE_PHI_OFFSET 0
#define LIGHTCONE_CHI_OFFSET 1
#define LIGHTCONE_B_OFFSET 2
#define LIGHTCONE_HIJ_OFFSET 5
// #define LIGHTCONE_MAX_FIELDS 10
// Modification
#define LIGHTCONE_ACHI_OFFSET 10
#define LIGHTCONE_AQ_OFFSET 11
#define LIGHTCONE_SGW_OFFSET 12
#define LIGHTCONE_MU_TOT_AC_OFFSET 13
#define LIGHTCONE_MU_PAR_OFFSET 14
#define LIGHTCONE_MU_TOT_OFFSET 15
#define LIGHTCONE_T00_OFFSET 16
#define LIGHTCONE_MAX_FIELDS 17
//-------------------

#ifndef MAX_PCL_SPECIES
#define MAX_PCL_SPECIES 3
#endif

#ifndef CYCLE_INFO_INTERVAL
#define CYCLE_INFO_INTERVAL 10
#endif
//----------Modification------- In below
#define MASK_PHI 1
#define MASK_CHI 2
#define MASK_POT 4
#define MASK_B 8
#define MASK_T00 16
#define MASK_TIJ 32
#define MASK_HIJ 64
#define MASK_GADGET 128
#define MASK_PCLS 1024
#define MASK_DELTA 256
#define MASK_ACHI 512
#define MASK_AQ 2048
#define MASK_MULTI 4096
#define MASK_T00_AS 8192
#define MASK_EOS_AS 16384
#define MASK_SOURCE 32768
#define MASK_HIJPRIME 65536
#define MASK_HIJPRIMENORM 131072

#define ICFLAG_CORRECT_DISPLACEMENT 1
#define ICFLAG_KSPHERE 2

// Identifiers for IC generator modules
#define ICGEN_BASIC 0
#define ICGEN_READ_FROM_DISK 1
#ifdef ICGEN_PREVOLUTION
#undef ICGEN_PREVOLUTION
#define ICGEN_PREVOLUTION 2
#endif
#ifdef ICGEN_SONG
#undef ICGEN_SONG
#define ICGEN_SONG 3
#endif
#ifdef ICGEN_FALCONIC
#undef ICGEN_FALCONIC
#define ICGEN_FALCONIC 4
#endif
#define ICGEN_CUSTOM 5

#define VECTOR_PARABOLIC 0
#define VECTOR_ELLIPTIC 1

// Physical constants
#define C_PLANCK_LAW 4.48147e-7		// omega_g / (T_cmb [K])^4
#define C_BOLTZMANN_CST 8.61733e-5	// Boltzmann constant [eV/K]
#define C_SPEED_OF_LIGHT 2997.92458 // speed of light [100 km/s]
#define C_RHO_CRIT 2.77459457e11	// critical density [M_sun h^2 / Mpc^3]
#define C_FD_NORM 1.80308535		// Integral[q*q/(exp(q)+1), 0, infinity]

// default physical parameters (used in parser.hpp)
#define P_HUBBLE 0.67556		// default value for h
#define P_T_NCDM 0.71611		// default value for T_ncdm
#define P_NCDM_MASS_OMEGA 93.14 // m_ncdm / omega_ncdm [eV]
#define P_N_UR 3.046			// default value for N_ur (= N_eff)
#define P_SPECTRAL_AMP 2.215e-9 // default value for A_s
#define P_SPECTRAL_INDEX 0.9619 // default value for n_s
#define P_PIVOT_SCALE 0.05		// default pivot scale [Mpc^-1]

#ifndef GADGET_LENGTH_CONVERSION
#define GADGET_LENGTH_CONVERSION 0.001 // Gadget length unit in Mpc / h
#endif
#ifndef GADGET_MASS_CONVERSION
#define GADGET_MASS_CONVERSION 1.0e10 // Gadget mass unit in M_sun / h
#endif
#ifndef GADGET_VELOCITY_CONVERSION
#define GADGET_VELOCITY_CONVERSION 3.335640952e-6 // Gadget velocity unit / speed of light
#endif
#ifndef GADGET_ID_BYTES
#define GADGET_ID_BYTES 8
#endif

#ifdef EXTERNAL_IO
#ifndef NUMBER_OF_IO_FILES
#define NUMBER_OF_IO_FILES 4
#endif
#endif

/* For use in snprintf.. */
#define SIZ 1000
#define SIZBUFFER 64
#define SIZFILENAME 200
#define SIZPARSTRING 512

// color escape sequences for terminal highlighting (enable with -DCOLORTERMINAL)
#ifdef COLORTERMINAL
#define COLORTEXT_WHITE "\033[37;1m"
#define COLORTEXT_CYAN "\033[36;1m"
#define COLORTEXT_GREEN "\033[32;1m"
#define COLORTEXT_RED "\033[31;1m"
#define COLORTEXT_YELLOW "\033[33;1m"
#define COLORTEXT_RESET "\033[0m"
#else
#define COLORTEXT_WHITE ""
#define COLORTEXT_CYAN ""
#define COLORTEXT_GREEN ""
#define COLORTEXT_RED ""
#define COLORTEXT_YELLOW ""
#define COLORTEXT_RESET ""
#endif

// header structure for GADGET-2 files [V. Springel, N. Yoshida, and S.D. White, New Astron. 6 (2001) 79
// and V. Springel, Mon. Not. R. Astron. Soc. 364 (2005) 1105]

#ifndef GADGET2_HEADER
#define GADGET2_HEADER
struct gadget2_header
{
	uint32_t npart[6];
	double mass[6];
	double time;
	double redshift;
	int32_t flag_sfr;
	int32_t flag_feedback;
	uint32_t npartTotal[6];
	int32_t flag_cooling;
	int32_t num_files;
	double BoxSize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	int32_t flag_age;
	int32_t flag_metals;
	uint32_t npartTotalHW[6];
	char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4]; /* fills to 256 Bytes */
};
#endif

#ifdef HAVE_HEALPIX
#include "chealpix.h"
#ifndef PIXBUFFER
#define PIXBUFFER 1048576
#endif

struct healpix_header
{
	uint32_t Nside;
	uint32_t Npix;
	uint32_t precision;
	uint32_t Ngrid;
	double direction[3];
	double distance;
	double boxsize;
	uint32_t Nside_ring;
	char fill[256 - 5 * 4 - 5 * 8]; /* fills to 256 Bytes */
};
#endif

struct lightcone_geometry
{
	double vertex[3];
	double z;
	double direction[3];
	double opening;
	double distance[2];
	double tauobs;
};

struct metadata
{
	int numpts;
	int downgrade_factor;
	int downgrade_factor2;
	long numpcl[MAX_PCL_SPECIES];
	int tracer_factor[MAX_PCL_SPECIES];
	int baryon_flag;
	int gr_flag;
	int vector_flag;
	int radiation_flag;
	int fluid_flag;
	int out_pk;
	int out_snapshot;
	int out_snapshot2;
	int out_animation;
	int out_lightcone[MAX_OUTPUTS];
	int num_pk;
	int numbins;
	int num_snapshot;
	int num_snapshot2;
	int num_lightcone;
	int num_restart;
	int Nside[MAX_OUTPUTS][2];
	double Cf;
	double movelimit;
	double steplimit;
	double boxsize;
	double dx;
	double wallclocklimit;
	double pixelfactor[MAX_OUTPUTS];
	double shellfactor[MAX_OUTPUTS];
	double covering[MAX_OUTPUTS];
	double z_in;
	double z_snapshot[MAX_OUTPUTS];
	double z_snapshot2[MAX_OUTPUTS];
	double z_pk[MAX_OUTPUTS];
	double z_restart[MAX_OUTPUTS];
	double z_switch_deltarad;
	double z_switch_linearchi;
	double z_switch_deltancdm[MAX_PCL_SPECIES - 2];
	double z_switch_Bncdm[MAX_PCL_SPECIES - 2];
	lightcone_geometry lightcone[MAX_OUTPUTS];
	char basename_lightcone[PARAM_MAX_LENGTH];
	char basename_snapshot[PARAM_MAX_LENGTH];
	char basename_pk[PARAM_MAX_LENGTH];
	char basename_generic[PARAM_MAX_LENGTH];
	char basename_animation[PARAM_MAX_LENGTH];
	char output_path[PARAM_MAX_LENGTH];
	char restart_path[PARAM_MAX_LENGTH];
	char basename_restart[PARAM_MAX_LENGTH];
	/*added in asevolution*/

	// Whether or not to evolve MG
	int fifth;
	// Which background to use: lcdm, mink, rho, etc
	char bckopt[PARAM_MAX_LENGTH];
	// .h5 snapshot to initialise scalar field from
	char achifile[PARAM_MAX_LENGTH];
	// .h5 snapshot initialise scalar velocity from
	char aqfile[PARAM_MAX_LENGTH];
	// Which field solver should be used: RK4, leapfrog, euler
	char As_solver[PARAM_MAX_LENGTH];
	// Which initialisation option for scalar: scaleinv, random, relaxedWall, gaussian, etc
	char As_ICopt[PARAM_MAX_LENGTH];
	// How many numbers of steps to use in evolving scalar for each external CDM evolution loop
	int nAs_numsteps;
	// Courant factor to use for scalar (instead of above parameter)
	double As_Courant;
	// Whether or not to add SE tensor of scalar to the one used in Einstein equations and background
	int As_source_gravity;
	// Whether to compute scalar SE tensor even if it is not used for evolution
	int As_compute_SE;
	// Whether to compute the scalar dynamically at background level
	int As_compute_dynbck;
	// Whether to use dynamical background computation for scalar
	int As_use_dynbck;
	// How many steps should be used for the dynamical background computation
	int nAs_numsteps_dynbck;
	// Whether to use fully relativistic geodesic w/ 5th force (w/ Bi)
	int As_5th_force;
	// Whether to use only Psi geodesic w/ 5th force
	int As_5th_force_Newton;
	// Whether to use the full T_cdm trace in the scalar equations of motion
	int As_full_trace;
	// Whether to use the Jordan frame potential when solving the Einstein equations
	int As_use_phiJ_SE;
	// Whether to update the particles
	int updatePtcls;
	// Whether to update the scalar
	int updateScalar;
	// Whether to solve Einstein equations
	int solveEinstein;
	// Whether to find the adiabatic velocity of scalar when solving quasi-statically
	int QSAq;
	// Whether the QSA initial guess should be from dynamical background calculation
	int QSA_guess_dyn;
	// Debug: Whether to modulate QSA result by some factor
	double QSAmodulateAmp;
	// Constraint solver tolerance
	double GStol;
	// Second constraint solve tolerance for initialisation
	double GStol2;
	// Scale factor at time of starting to evolve scalar
	double aMG;
	// New CDM Courant factor to use when starting to evolve scalar
	double Cf_MG;
	// Amplitude of scalar ICs, whether random or scaleinv
	double ICamp;
	// Amplitude of scalar ICs for velocity
	double ICamp2;
	// IC choice for custom particle setup initialisation:
	// centreWall
	char custom_pcl_IC[PARAM_MAX_LENGTH];
	// IC type choice for custom..: field or particles
	char custom_pcl_IC_type[PARAM_MAX_LENGTH];
	// Debug: whether to stop program after initialisation
	int breakafterzini;
	// Whether to stop program after redshift z
	double breakafterz;
	// Whether to stop program after domain walls have vanished
	int breakWhenNoDW;
	// Whether or not to output snapshots for cycle=0
	int writefirstsnap;
	int writefirstspec;

	// animation
	// whether to have animation
	int animation;
	// what z to start animation
	double zanimationstart;
	// what z to end animation
	double zanimationend;
	// output animation every nth cycle
	int animationeveryn;
	// what downgrade factor to use for animation outputs
	int animation_downgrade_factor;
	// what thickness to use for animation outputs
	int animation_thickness;
	// which x coordinate (node nr) to use for the 2D snapshot output
	int animation_xcoord;

	// GW
	// whether to have adiabatic difference between current and previous solutions for beyond Courant condition modes
	int GWhighkAdiabatic;
	// whether to have adiabatic difference for ICs
	int GWstartAdiabatic;
	// which cycle to initialise GWs at
	int GWinitialiseCycle;
	// whether to read initial dynamic fields (hij and hijprime) ! NB: read in real space than transformed
	int GWreadFields;
	char hijfile[PARAM_MAX_LENGTH];
	char hijprimefile[PARAM_MAX_LENGTH];
};

struct icsettings
{
	int numtile[MAX_PCL_SPECIES];
	int seed;
	int flags;
	int generator;
	int restart_cycle;
	char pclfile[MAX_PCL_SPECIES][PARAM_MAX_LENGTH];
	char pkfile[PARAM_MAX_LENGTH];
	char tkfile[PARAM_MAX_LENGTH];
	char metricfile[3][PARAM_MAX_LENGTH];
	double restart_tau;
	double restart_dtau;
	double restart_version;
	double z_ic;
	double z_relax;
	double Cf;
	double A_s;
	double n_s;
	double k_pivot;
};

//----------Modification------- adding parameters below
struct cosmology
{
	double Omega_cdm;
	double Omega_b;
	double Omega_m;
	// Set by flatness criterion at background level
	double Omega_Lambda;
	// Constant term energy density for changing Omega_Lambda effectively
	double Omega_asymmetron;
	// Energy density for symmetron background
	double Omega_sym; // evaluated in the code
	// Constant term energy density for symmetron background
	double Omega_sym0;
	// must be in order that they are found

	// Compton wavelength as a fraction of Hubble scale
	double xistar;
	// scale factor at symmetry breaking
	double astar;
	// fifth force strength relative to Newtonian
	double betastar;
	// Second fifth force strength in case of asymmetron
	double betastar2;
	// Compton wavelength in Mpc/h
	double Lc_as;
	// Lagrangian mass parameter
	double mu_as;
	// redshift of symmetry breaking
	double zssb_as;
	// Conformal coupling scale
	double M_as;
	// Lagrangian coupling strength parameter
	double lambda_as;
	// Lagrangian cubic term parameter
	double kappa_as;
	// Dimensionless parameter
	double kappa2_as;
	// Sign of cubic term in potential
	double kappa_sign;
	// Prefactor of potential in equations of motion
	double prefac;
	// Factor in baryonic coupling term in equations of motion
	double aetanorm;
	// term in conformal coupling factor
	double dAcnf;
	// for particle horizon integrand
	double avgsource;
	//--
	double Omega_fld;
	double w0_fld;
	double wa_fld;
	double cs2_fld;
	double Omega_g;
	double Omega_ur;
	double Omega_rad;
	double Omega_ncdm[MAX_PCL_SPECIES - 2];
	double h;
	double m_ncdm[MAX_PCL_SPECIES - 2];
	double T_ncdm[MAX_PCL_SPECIES - 2];
	double deg_ncdm[MAX_PCL_SPECIES - 2];
	int num_ncdm;
};

#endif
