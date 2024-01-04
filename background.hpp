//////////////////////////
// background.hpp
//////////////////////////
//
// code components related to background evolution
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: September 2018
//
//////////////////////////

#ifndef BACKGROUND_HEADER
#define BACKGROUND_HEADER

#include <gsl/gsl_integration.h>

double FermiDiracIntegrand(double q, void *w)
{
	return q * q * sqrt(q * q + *(double *)w) / (exp(q) + 1.0l);
}

//////////////////////////
// FermiDiracIntegral
//////////////////////////
// Description:
//   computes the integral of the relativistic Fermi-Dirac distribution
//
// Arguments:
//   w          parameter in the F-D distribution, "(m a / kB T)^2"
//
// Returns: value for the integral
//
//////////////////////////

double FermiDiracIntegral(double &w)
{
	double result;
	gsl_function f;
	double err;
	size_t n;

	f.function = &FermiDiracIntegrand;
	f.params = &w;

	gsl_integration_qng(&f, 0.0l, 24.0l, 5.0e-7, 1.0e-7, &result, &err, &n);

	return result;
}

//////////////////////////
// bg_ncdm (1)
//////////////////////////
// Description:
//   computes the background model for one ncdm species by integrating the relativistic
//   Fermi-Dirac distribution
//
// Arguments:
//   a          scale factor at which to compute the background model
//   cosmo      structure containing the cosmological parameters
//   p          index of the ncdm species
//
// Returns: value for the background model
//
//////////////////////////

double bg_ncdm(const double a, const cosmology cosmo, const int p)
{
	if (p < 0 || p >= cosmo.num_ncdm)
		return 0;
	else
	{
		double w = a * cosmo.m_ncdm[p] / (pow(cosmo.Omega_g * cosmo.h * cosmo.h / C_PLANCK_LAW, 0.25) * cosmo.T_ncdm[p] * C_BOLTZMANN_CST);
		w *= w;

		return FermiDiracIntegral(w) * cosmo.Omega_ncdm[p] * pow(cosmo.Omega_g * cosmo.h * cosmo.h / C_PLANCK_LAW, 0.25) * cosmo.T_ncdm[p] * C_BOLTZMANN_CST / cosmo.m_ncdm[p] / C_FD_NORM / a;
	}
}

//////////////////////////
// bg_ncdm (2)
//////////////////////////
// Description:
//   computes the background model for all ncdm species by integrating the relativistic
//   Fermi-Dirac distribution
//
// Arguments:
//   a          scale factor at which to compute the background model
//   cosmo      structure containing the cosmological parameters
//
// Note:
//   For optimization, the last value of a is stored in a static variable such that
//   multiple calls at the same value of a will not result in multiple integrations
//   being carried out. This assumes that the cosmological model should not change!
//
// Returns: value for the background model
//
//////////////////////////

double bg_ncdm(const double a, const cosmology cosmo)
{
	double w;
	static double result = -1.0;
	static double a_prev = -1.0;

	if (a != a_prev)
	{
		result = 0.0;
		a_prev = a;

		for (int p = 0; p < cosmo.num_ncdm; p++)
			result += bg_ncdm(a, cosmo, p);
	}

	return result;
}

//////////////////////////
// Hconf
//////////////////////////
// Description:
//   computes the conformal Hubble rate at given scale factor
//
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//
// Returns: conformal Hubble rate
//
//////////////////////////
//---------------Modification------------------------
// NB, since the intergrator evaluates the derivative for a, we really should find Omega_sym inside of the functions
// This is a shortcoming of this procedure

double (*Hconf)(const double a, const double fourpiG, const cosmology cosmo, double avgsource);

double Hconf_old(const double a, const double fourpiG, const cosmology cosmo, double avgsource)
{ // use lcdm in ICs anyway
	return sqrt((2. * fourpiG / 3.) * (((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) / a) + (cosmo.Omega_Lambda * a * a) + (cosmo.Omega_rad / a / a) + (cosmo.Omega_fld * exp(3. * cosmo.wa_fld * (a - 1.)) / pow(a, 1. + 3. * (cosmo.w0_fld + cosmo.wa_fld)))));
}
double Hconf_sym(const double a, const double fourpiG, const cosmology cosmo, double avgsource)
{
	return sqrt((2. * fourpiG / 3.) * (((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) / a) + ((cosmo.Omega_sym + cosmo.Omega_Lambda) * a * a) + (cosmo.Omega_rad / a / a)));
}
double Hconf_rho(const double a, const double fourpiG, const cosmology cosmo, double avgsource)
{
	return sqrt((2. * fourpiG / 3.) * avgsource / a);
}
double Hconf_mink(const double a, const double fourpiG, const cosmology cosmo, double avgsource)
{
	return 0;
}

double daqdtau(double achiB, double a, const cosmology cosmo, const metadata sim)
{
	return -pow(a * a * cosmo.mu_as * sim.boxsize, 2) * (achiB * (achiB * (achiB + cosmo.kappa_as / (cosmo.mu_as * sqrt(cosmo.lambda_as))) + (cosmo.Omega_m * 3. / (a * pow(2998 * a * cosmo.mu_as * cosmo.M_as, 2)) - 1.)));
}

void rungekutta4dtau(double &dtau, double a, const double fourpiG, const double da, const cosmology cosmo, double avgsource)
{
	double k1a, k2a, k4a;

	k1a = 1. / (a * Hconf(a, fourpiG, cosmo, avgsource));
	k2a = 1. / ((a + da / 2.) * Hconf(a + da / 2., fourpiG, cosmo, avgsource));
	k4a = 1. / ((a + da) * Hconf(a + da, fourpiG, cosmo, avgsource));

	dtau = da * (k1a + 4. * k2a + k4a) / 6.;
}

void rungekutta4achi(double &achiB, double &aqB, double &a, const double fourpiG,
					 const double dtau, const cosmology cosmo, const metadata sim, double avgsource)
{
	double k1a, k2a, k3a, k4a;
	double k1b, k2b, k3b, k4b;

	k1a = a * Hconf(a, fourpiG, cosmo, avgsource);
	k2a = (a + k1a * dtau / 2.) * Hconf(a + k1a * dtau / 2., fourpiG, cosmo, avgsource);
	k3a = (a + k2a * dtau / 2.) * Hconf(a + k2a * dtau / 2., fourpiG, cosmo, avgsource);
	k4a = (a + k3a * dtau) * Hconf(a + k3a * dtau, fourpiG, cosmo, avgsource);

	a += dtau * (k1a + 2. * k2a + 2. * k3a + k4a) / 6.;

	// Euler-Cromer variant
	k1b = 1 * daqdtau(achiB, a, cosmo, sim);
	k2b = (1 + k1a * dtau / 2.) * daqdtau(achiB, a + k1a * dtau / 2., cosmo, sim);
	k3b = (1 + k2a * dtau / 2.) * daqdtau(achiB, a + k2a * dtau / 2., cosmo, sim);
	k4b = (1 + k3a * dtau) * daqdtau(achiB, a + k3a * dtau, cosmo, sim);
	aqB += dtau * (k1b + 2. * k2b + 2. * k3b + k4b) / 6.;

	achiB += dtau * aqB / (a * a);
}

// These are only used in ic_basic.hpp, therefore it is fair to treat asymmetron as either zero or some negligible cosmological constant
double Omega_m(const double a, const cosmology cosmo) { return cosmo.Omega_m / (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) + cosmo.Omega_asymmetron * a * a * a + cosmo.Omega_Lambda * a * a * a + cosmo.Omega_rad / a + cosmo.Omega_fld * exp(3. * cosmo.wa_fld * (a - 1.)) / pow(a, 3. * (cosmo.w0_fld + cosmo.wa_fld))); }
double Omega_rad(const double a, const cosmology cosmo) { return (cosmo.Omega_rad + (bg_ncdm(a, cosmo) + cosmo.Omega_cdm + cosmo.Omega_b - cosmo.Omega_m) * a) / ((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * a + cosmo.Omega_asymmetron * a * a * a * a + cosmo.Omega_Lambda * a * a * a * a + cosmo.Omega_rad + cosmo.Omega_fld * exp(3. * cosmo.wa_fld * (a - 1.)) / pow(a, 3. * (cosmo.w0_fld + cosmo.wa_fld) - 1.)); }
double Omega_Lambda(const double a, const cosmology cosmo) { return cosmo.Omega_Lambda / ((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) / a / a / a + cosmo.Omega_Lambda + cosmo.Omega_asymmetron + cosmo.Omega_rad / a / a / a / a + cosmo.Omega_fld * exp(3. * cosmo.wa_fld * (a - 1.)) / pow(a, 3. + 3. * (cosmo.w0_fld + cosmo.wa_fld))); }

//////////////////////////
// rungekutta4bg
//////////////////////////
// Description:
//   integrates the Friedmann equation for the background model using a fourth-order
//   Runge-Kutta method
//
// Arguments:
//   a          scale factor (will be advanced by dtau)
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//   dtau       time step by which the scale factor should be advanced
//
// Returns:
//
//////////////////////////
// ------------modification--- in parameters and Hconf call

void rungekutta4bg(double &a, const double fourpiG, const double dtau, const cosmology cosmo, double avgsource)
{
	double k1a, k2a, k3a, k4a;

	k1a = a * Hconf(a, fourpiG, cosmo, avgsource);
	k2a = (a + k1a * dtau / 2.) * Hconf(a + k1a * dtau / 2., fourpiG, cosmo, avgsource);
	k3a = (a + k2a * dtau / 2.) * Hconf(a + k2a * dtau / 2., fourpiG, cosmo, avgsource);
	k4a = (a + k3a * dtau) * Hconf(a + k3a * dtau, fourpiG, cosmo, avgsource);

	a += dtau * (k1a + 2. * k2a + 2. * k3a + k4a) / 6.;
}

void rungekutta4bg_old(double &a, const double fourpiG, const double dtau, const cosmology cosmo, double avgsource)
{
	double k1a, k2a, k3a, k4a;

	k1a = a * Hconf_old(a, fourpiG, cosmo, avgsource);
	k2a = (a + k1a * dtau / 2.) * Hconf_old(a + k1a * dtau / 2., fourpiG, cosmo, avgsource);
	k3a = (a + k2a * dtau / 2.) * Hconf_old(a + k2a * dtau / 2., fourpiG, cosmo, avgsource);
	k4a = (a + k3a * dtau) * Hconf_old(a + k3a * dtau, fourpiG, cosmo, avgsource);

	a += dtau * (k1a + 2. * k2a + 2. * k3a + k4a) / 6.;
}

// To fix the below, integrate particleHorizon inside of the code
double particleHorizonIntegrand(double sqrta, void *cosmo)
{
	return 2. / (sqrta * Hconf_old(sqrta * sqrta, 1., (*(cosmology *)cosmo), (*(cosmology *)cosmo).avgsource)); // NB: must be Hconf_old since should be function of a // only used early redshift..
}
double particleHorizonIntegrand_old(double sqrta, void *cosmo)
{
	return 2. / (sqrta * Hconf_old(sqrta * sqrta, 1., (*(cosmology *)cosmo), (*(cosmology *)cosmo).avgsource));
}

//////////////////////////
// particleHorizon
//////////////////////////
// Description:
//   computes the particle horizon (tau) at given scale factor
//
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//
// Returns: particle horizon (tau)
//
//////////////////////////

double particleHorizon(const double a, const double fourpiG, cosmology &cosmo, double &avgsource)
{
	double result;
	gsl_function f;
	double err;
	size_t n;
	cosmo.avgsource = avgsource;
	f.function = &particleHorizonIntegrand;
	f.params = &cosmo;

	gsl_integration_qng(&f, sqrt(a) * 1.0e-7, sqrt(a), 5.0e-7, 1.0e-7, &result, &err, &n);

	return result / sqrt(fourpiG);
}
double particleHorizon_old(const double a, const double fourpiG, cosmology &cosmo, double &avgsource)
{
	double result;
	gsl_function f;
	double err;
	size_t n;
	f.function = &particleHorizonIntegrand_old;
	f.params = &cosmo;

	gsl_integration_qng(&f, sqrt(a) * 1.0e-7, sqrt(a), 5.0e-7, 1.0e-7, &result, &err, &n);

	return result / sqrt(fourpiG);
}
double particleHorizon_old2(const double a, const double fourpiG, cosmology &cosmo)
{
	double result;
	gsl_function f;
	double err;
	size_t n;
	f.function = &particleHorizonIntegrand_old;
	f.params = &cosmo;

	gsl_integration_qng(&f, sqrt(a) * 1.0e-7, sqrt(a), 5.0e-7, 1.0e-7, &result, &err, &n);

	return result / sqrt(fourpiG);
}

#endif
