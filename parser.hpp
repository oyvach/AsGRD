//////////////////////////
// parser.hpp
//////////////////////////
//
// Parser for settings file
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: April 2019
//
//////////////////////////
// then modifed.. etc.. see main file.. applies to all files

#ifndef PARSER_HEADER
#define PARSER_HEADER

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "metadata.hpp"

//--modification--
#ifndef PI
#define PI 3.141593
#endif
// for conversion from eV to code units
#define hbarSI 1.0545718e-34
#define eVSI 1.60217662e-19
#define cSI 2.99792458e8
#define GSI 6.67430e-11

using namespace std;

struct parameter
{
	char name[PARAM_MAX_LENGTH];
	char value[PARAM_MAX_LENGTH];
	bool used;
};

int sort_descending(const void *z1, const void *z2)
{
	return (*(double *)z1 > *(double *)z2) ? -1 : ((*(double *)z1 < *(double *)z2) ? 1 : 0);
}

//////////////////////////
// readline
//////////////////////////
// Description:
//   reads a line of characters and checks if it declares a parameter; if yes,
//   i.e. the line has the format "<parameter name> = <parameter value>" (with
//   an optional comment added, preceded by a hash-symbol, '#'), the parameter
//   name and value are copied to the corresponding arrays, and 'true' is returned.
//   If the format is not recognized (or the line is commented using the hash-symbol)
//   'false' is returned instead.
//
// Arguments:
//   line       string containing the line to be read
//   pname      will contain the name of the declared parameter (if found)
//   pvalue     will contain the value of the declared parameter (if found)
//
// Returns:
//   'true' if a parameter is declared in the line, 'false' otherwise.
//
//////////////////////////

bool readline(char *line, char *pname, char *pvalue)
{
	char *pequal;
	char *phash;
	char *l;
	char *r;

	pequal = strchr(line, '=');

	if (pequal == NULL || pequal == line)
		return false;

	phash = strchr(line, '#');

	if (phash != NULL && phash < pequal)
		return false;

	l = line;
	while (*l == ' ' || *l == '\t')
		l++;

	r = pequal - 1;
	while ((*r == ' ' || *r == '\t') && r > line)
		r--;

	if (r < l)
		return false;

	if (r - l + 1 >= PARAM_MAX_LENGTH)
		return false;

	strncpy(pname, l, r - l + 1);
	pname[r - l + 1] = '\0';

	l = pequal + 1;
	while (*l == ' ' || *l == '\t')
		l++;

	if (phash == NULL)
		r = line + strlen(line) - 1;
	else
		r = phash - 1;

	while (*r == ' ' || *r == '\t' || *r == '\n' || *r == '\r')
		r--;

	if (r < l)
		return false;

	if (r - l + 1 >= PARAM_MAX_LENGTH)
		return false;

	strncpy(pvalue, l, r - l + 1);
	pvalue[r - l + 1] = '\0';

	return true;
}

//////////////////////////
// loadParameterFile
//////////////////////////
// Description:
//   loads a parameter file and creates an array of parameters declared therein
//
// Arguments:
//   filename   string containing the path to the parameter file
//   params     will contain the array of parameters (memory will be allocated)
//
// Returns:
//   number of parameters defined in the parameter file (= length of parameter array)
//
//////////////////////////

int loadParameterFile(const char *filename, parameter *&params)
{
	int numparam = 0;
	int i = 0;

#ifdef LATFIELD2_HPP
	if (parallel.grid_rank()[0] == 0) // read file
	{
#endif
		FILE *paramfile;
		char line[PARAM_MAX_LINESIZE];
		char pname[PARAM_MAX_LENGTH];
		char pvalue[PARAM_MAX_LENGTH];

		paramfile = fopen(filename, "r");

		if (paramfile == NULL)
		{
#ifdef LATFIELD2_HPP
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! Unable to open parameter file " << filename << "." << endl;
			parallel.abortForce();
#else
		cerr << " error in loadParameterFile! Unable to open parameter file " << filename << "." << endl;
		return -1;
#endif
		}

		while (!feof(paramfile) && !ferror(paramfile))
		{
			if (fgets(line, PARAM_MAX_LINESIZE, paramfile) == NULL)
				break;

			if (readline(line, pname, pvalue) == true)
				numparam++;
		}

		if (numparam == 0)
		{
			fclose(paramfile);
#ifdef LATFIELD2_HPP
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! No valid data found in file " << filename << "." << endl;
			parallel.abortForce();
#else
		cerr << " error in loadParameterFile! No valid data found in file " << filename << "." << endl;
		return -1;
#endif
		}

		params = (parameter *)malloc(sizeof(parameter) * numparam);

		if (params == NULL)
		{
			fclose(paramfile);
#ifdef LATFIELD2_HPP
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! Memory error." << endl;
			parallel.abortForce();
#else
		cerr << " error in loadParameterFile! Memory error." << endl;
		return -1;
#endif
		}

		rewind(paramfile);

		while (!feof(paramfile) && !ferror(paramfile) && i < numparam)
		{
			if (fgets(line, PARAM_MAX_LINESIZE, paramfile) == NULL)
				break;

			if (readline(line, params[i].name, params[i].value) == true)
			{
				params[i].used = false;
				i++;
			}
		}

		fclose(paramfile);

		if (i < numparam)
		{
			free(params);
#ifdef LATFIELD2_HPP
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! File may have changed or file pointer corrupted." << endl;
			parallel.abortForce();
#else
		cerr << " error in loadParameterFile! File may have changed or file pointer corrupted." << endl;
		return -1;
#endif
		}

#ifdef LATFIELD2_HPP
		parallel.broadcast_dim0<int>(numparam, 0);
	}
	else
	{
		parallel.broadcast_dim0<int>(numparam, 0);

		params = (parameter *)malloc(sizeof(parameter) * numparam);

		if (params == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! Memory error." << endl;
			parallel.abortForce();
		}
	}

	parallel.broadcast_dim0<parameter>(params, numparam, 0);
#endif

	return numparam;
}

//////////////////////////
// saveParameterFile
//////////////////////////
// Description:
//   saves a parameter file
//
// Arguments:
//   filename   string containing the path to the parameter file
//   params     array of parameters
//   numparam   length of parameter array
//   used_only  if 'true', only the used parameters will be written (default)
//
// Returns:
//
//////////////////////////

void saveParameterFile(const char *filename, parameter *params, const int numparam, bool used_only = true)
{
#ifdef LATFIELD2_HPP
	if (parallel.isRoot())
#endif
	{
		FILE *paramfile;

		paramfile = fopen(filename, "w");

		if (paramfile == NULL)
		{
			cout << " error in saveParameterFile! Unable to open file " << filename << "." << endl;
		}
		else
		{
			for (int i = 0; i < numparam; i++)
			{
				if (!used_only || params[i].used)
					fprintf(paramfile, "%s = %s\n", params[i].name, params[i].value);
			}

			fclose(paramfile);
		}
	}
}

//////////////////////////
// parseParameter (int)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value as integer
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to integer which will contain the parsed parameter value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter *&params, const int numparam, const char *pname, int &pvalue)
{
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			if (sscanf(params[i].value, "%d", &pvalue) == 1)
			{
				params[i].used = true;
				return true;
			}
		}
	}

	return false;
}

//////////////////////////
// parseParameter (long)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value as integer
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to integer which will contain the parsed parameter value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter *&params, const int numparam, const char *pname, long &pvalue)
{
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			if (sscanf(params[i].value, "%ld", &pvalue) == 1)
			{
				params[i].used = true;
				return true;
			}
		}
	}

	return false;
}

//////////////////////////
// parseParameter (double)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value as double
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to double which will contain the parsed parameter value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter *&params, const int numparam, const char *pname, double &pvalue)
{
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			if (sscanf(params[i].value, "%lf", &pvalue) == 1)
			{
				params[i].used = true;
				return true;
			}
		}
	}

	return false;
}

//////////////////////////
// parseParameter (char *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and retrieves its value as string
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     character string which will contain a copy of the parameter value (if found)
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter *&params, const int numparam, const char *pname, char *pvalue)
{
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			strcpy(pvalue, params[i].value);
			params[i].used = true;
			return true;
		}
	}

	return false;
}

//////////////////////////
// parseParameter (double *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated double values
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of double values which will contain the list of parsed parameters (if found)
//   nmax       maximum size of array; will be set to the actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter *&params, const int numparam, const char *pname, double *pvalue, int &nmax)
{
	char *start;
	char *comma;
	char item[PARAM_MAX_LENGTH];
	int n = 0;

	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			start = params[i].value;
			if (nmax > 1)
			{
				while ((comma = strchr(start, ',')) != NULL)
				{
					strncpy(item, start, comma - start);
					item[comma - start] = '\0';
					if (sscanf(item, " %lf ", pvalue + n) != 1)
					{
						nmax = n;
						return false;
					}
					start = comma + 1;
					if (++n > nmax - 2)
						break;
				}
			}
			if (sscanf(start, " %lf ", pvalue + n) != 1)
			{
				nmax = n;
				return false;
			}
			nmax = ++n;
			params[i].used = true;
			return true;
		}
	}

	nmax = 0;
	return false;
}

//////////////////////////
// parseParameter (int *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated integer values
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of integer values which will contain the list of parsed parameters (if found)
//   nmax       maximum size of array; will be set to the actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter *&params, const int numparam, const char *pname, int *pvalue, int &nmax)
{
	char *start;
	char *comma;
	char item[PARAM_MAX_LENGTH];
	int n = 0;

	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			start = params[i].value;
			if (nmax > 1)
			{
				while ((comma = strchr(start, ',')) != NULL)
				{
					strncpy(item, start, comma - start);
					item[comma - start] = '\0';
					if (sscanf(item, " %d ", pvalue + n) != 1)
					{
						nmax = n;
						return false;
					}
					start = comma + 1;
					if (++n > nmax - 2)
						break;
				}
			}
			if (sscanf(start, " %d ", pvalue + n) != 1)
			{
				nmax = n;
				return false;
			}
			nmax = ++n;
			params[i].used = true;
			return true;
		}
	}

	nmax = 0;
	return false;
}

//////////////////////////
// parseParameter (char **)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated strings
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of character strings which will contain the list of parsed parameters (if found)
//   nmax       maximum size of array; will be set to the actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter *&params, const int numparam, const char *pname, char **pvalue, int &nmax)
{
	char *start;
	char *comma;
	char *l;
	char *r;
	char item[PARAM_MAX_LENGTH];
	int n = 0;

	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			start = params[i].value;
			if (nmax > 1)
			{
				while ((comma = strchr(start, ',')) != NULL)
				{
					l = start;
					while (*l == ' ' || *l == '\t')
						l++;
					r = comma - 1;
					while ((*r == ' ' || *r == '\t') && r > start)
						r--;

					if (r < l)
					{
						nmax = n;
						return false;
					}

					strncpy(pvalue[n], l, r - l + 1);
					pvalue[n][r - l + 1] = '\0';

					start = comma + 1;
					if (++n > nmax - 2)
						break;
				}
			}
			l = start;
			while (*l == ' ' || *l == '\t')
				l++;
			r = l;
			while (*r != ' ' && *r != '\t' && *r != '\0')
				r++;
			r--;

			if (r < l)
			{
				nmax = n;
				return false;
			}

			strncpy(pvalue[n], l, r - l + 1);
			pvalue[n][r - l + 1] = '\0';

			nmax = ++n;
			params[i].used = true;
			return true;
		}
	}

	nmax = 0;
	return false;
}

//////////////////////////
// parseFieldSpecifiers
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated field specifiers
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     integer which will contain the binary-encoded list of parsed specifiers (if found)
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseFieldSpecifiers(parameter *&params, const int numparam, const char *pname, int &pvalue)
{
	char *start;
	char *comma;
	int pos;
	char item[PARAM_MAX_LENGTH];

	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			pvalue = 0;
			start = params[i].value;
			while ((comma = strchr(start, ',')) != NULL)
			{
				strncpy(item, start, comma - start);
				for (pos = comma - start; pos > 0; pos--)
				{
					if (item[pos - 1] != ' ' && item[pos - 1] != '\t')
						break;
				}
				item[pos] = '\0';

				if (strcmp(item, "Phi") == 0 || strcmp(item, "phi") == 0)
					pvalue |= MASK_PHI;
				/*--------------Modification------------------*/
				else if (strcmp(item, "achi") == 0 || strcmp(item, "Achi") == 0)
					pvalue |= MASK_ACHI;
				else if (strcmp(item, "aq") == 0 || strcmp(item, "Aq") == 0)
					pvalue |= MASK_AQ;
				else if (strcmp(item, "T00_As") == 0 || strcmp(item, "T00_as") == 0)
					pvalue |= MASK_T00_AS;
				else if (strcmp(item, "eos_As") == 0 || strcmp(item, "eos_as") == 0)
					pvalue |= MASK_EOS_AS;
				else if (strcmp(item, "source") == 0 || strcmp(item, "rhocdm") == 0)
					pvalue |= MASK_SOURCE;
				else if (strcmp(item, "hij_prime") == 0 || strcmp(item, "HIJ_PRIME") == 0)
					pvalue |= MASK_HIJPRIME;
				else if (strcmp(item, "hij_prime_norm") == 0 || strcmp(item, "HIJ_PRIME_NORM") == 0)
					pvalue |= MASK_HIJPRIMENORM;
				/*--------------------------------------------*/
				else if (strcmp(item, "Chi") == 0 || strcmp(item, "chi") == 0)
					pvalue |= MASK_CHI;
				else if (strcmp(item, "Pot") == 0 || strcmp(item, "pot") == 0 || strcmp(item, "Psi_N") == 0 || strcmp(item, "psi_N") == 0 || strcmp(item, "PsiN") == 0 || strcmp(item, "psiN") == 0)
					pvalue |= MASK_POT;
				else if (strcmp(item, "B") == 0 || strcmp(item, "Bi") == 0)
					pvalue |= MASK_B;
				else if (strcmp(item, "T00") == 0 || strcmp(item, "rho") == 0)
					pvalue |= MASK_T00;
				else if (strcmp(item, "Tij") == 0)
					pvalue |= MASK_TIJ;
				else if (strcmp(item, "hij") == 0 || strcmp(item, "GW") == 0)
					pvalue |= MASK_HIJ;
				else if (strcmp(item, "Gadget") == 0 || strcmp(item, "Gadget2") == 0 || strcmp(item, "gadget") == 0 || strcmp(item, "gadget2") == 0)
					pvalue |= MASK_GADGET;
				else if (strcmp(item, "multi-Gadget") == 0 || strcmp(item, "multi-Gadget2") == 0 || strcmp(item, "multi-gadget") == 0 || strcmp(item, "multi-gadget2") == 0)
					pvalue |= MASK_GADGET | MASK_MULTI;
				else if (strcmp(item, "Particles") == 0 || strcmp(item, "particles") == 0 || strcmp(item, "pcls") == 0 || strcmp(item, "part") == 0)
					pvalue |= MASK_PCLS;
				else if (strcmp(item, "delta") == 0 || strcmp(item, "Ds") == 0 || strcmp(item, "D_s") == 0)
					pvalue |= MASK_DELTA;
				start = comma + 1;
				while (*start == ' ' || *start == '\t')
					start++;
			}

			if (strcmp(start, "Phi") == 0 || strcmp(start, "phi") == 0)
				pvalue |= MASK_PHI;
			/*---------------Modification-----------------*/
			else if (strcmp(start, "Achi") == 0 || strcmp(start, "achi") == 0)
				pvalue |= MASK_ACHI;
			else if (strcmp(start, "Aq") == 0 || strcmp(start, "aq") == 0)
				pvalue |= MASK_AQ;
			else if (strcmp(start, "T00_As") == 0 || strcmp(start, "T00_as") == 0)
				pvalue |= MASK_T00_AS;
			else if (strcmp(start, "eos_As") == 0 || strcmp(start, "eos_as") == 0)
				pvalue |= MASK_EOS_AS;
			else if (strcmp(start, "source") == 0 || strcmp(start, "rhocdm") == 0)
				pvalue |= MASK_SOURCE;
			else if (strcmp(start, "hij_prime") == 0 || strcmp(start, "HIJ_PRIME") == 0)
				pvalue |= MASK_HIJPRIME;
			else if (strcmp(start, "hij_prime_norm") == 0 || strcmp(start, "HIJ_PRIME_NORM") == 0)
				pvalue |= MASK_HIJPRIMENORM;
			/*--------------------------------------------*/
			else if (strcmp(start, "Chi") == 0 || strcmp(start, "chi") == 0)
				pvalue |= MASK_CHI;
			else if (strcmp(start, "Pot") == 0 || strcmp(start, "pot") == 0 || strcmp(start, "Psi_N") == 0 || strcmp(start, "psi_N") == 0 || strcmp(start, "PsiN") == 0 || strcmp(start, "psiN") == 0)
				pvalue |= MASK_POT;
			else if (strcmp(start, "B") == 0 || strcmp(start, "Bi") == 0)
				pvalue |= MASK_B;
			else if (strcmp(start, "T00") == 0 || strcmp(start, "rho") == 0)
				pvalue |= MASK_T00;
			else if (strcmp(start, "Tij") == 0)
				pvalue |= MASK_TIJ;
			else if (strcmp(start, "hij") == 0 || strcmp(start, "GW") == 0)
				pvalue |= MASK_HIJ;
			else if (strcmp(start, "Gadget") == 0 || strcmp(start, "Gadget2") == 0 || strcmp(start, "gadget") == 0 || strcmp(start, "gadget2") == 0)
				pvalue |= MASK_GADGET;
			else if (strcmp(start, "multi-Gadget") == 0 || strcmp(start, "multi-Gadget2") == 0 || strcmp(start, "multi-gadget") == 0 || strcmp(start, "multi-gadget2") == 0)
				pvalue |= MASK_GADGET | MASK_MULTI;
			else if (strcmp(start, "Particles") == 0 || strcmp(start, "particles") == 0 || strcmp(start, "pcls") == 0 || strcmp(start, "part") == 0)
				pvalue |= MASK_PCLS;
			else if (strcmp(start, "delta") == 0 || strcmp(start, "Ds") == 0 || strcmp(start, "D_s") == 0)
				pvalue |= MASK_DELTA;

			params[i].used = true;
			return true;
		}
	}

	return false;
}

//////////////////////////
// parseMetadata
//////////////////////////
// Description:
//   parses all metadata from the parameter array
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   sim        reference to metadata stucture (holds simulation parameters)
//   cosmo      reference to cosmology structure (holds cosmological parameters)
//   ic         reference to icsettings structure (holds settings for IC generation)
//
// Returns:
//   number of parameters parsed
//
//////////////////////////

#ifndef LATFIELD2_HPP
#define COUT cout
#endif

int parseMetadata(parameter *&params, const int numparam, metadata &sim, cosmology &cosmo, icsettings &ic)
{
	char par_string[PARAM_MAX_LENGTH];
	char *ptr;
	char *pptr[MAX_PCL_SPECIES];
	int usedparams = 0;
	int i;
	double tmp;

	// parse settings for IC generator

	ic.pkfile[0] = '\0';
	ic.tkfile[0] = '\0';
	ic.metricfile[0][0] = '\0';
	ic.metricfile[1][0] = '\0';
	ic.metricfile[2][0] = '\0';
	ic.seed = 0;
	ic.flags = 0;
	ic.z_ic = -2.;
	ic.z_relax = -2.;
	ic.Cf = -1.0;
	ic.A_s = P_SPECTRAL_AMP;
	ic.n_s = P_SPECTRAL_INDEX;
	ic.k_pivot = P_PIVOT_SCALE;
	ic.restart_cycle = -1;
	ic.restart_tau = 0.;
	ic.restart_dtau = 0.;
	ic.restart_version = -1.;

	parseParameter(params, numparam, "seed", ic.seed);

	if (parseParameter(params, numparam, "IC generator", par_string))
	{
		if (par_string[0] == 'B' || par_string[0] == 'b')
			ic.generator = ICGEN_BASIC;
		else if (par_string[0] == 'R' || par_string[0] == 'r')
			ic.generator = ICGEN_READ_FROM_DISK;
		else if (par_string[0] == 'C' || par_string[0] == 'c')
			ic.generator = ICGEN_CUSTOM;
#ifdef ICGEN_PREVOLUTION
		else if (par_string[0] == 'P' || par_string[0] == 'p')
			ic.generator = ICGEN_PREVOLUTION;
#endif
#ifdef ICGEN_SONG
		else if (par_string[0] == 'S' || par_string[0] == 's')
			ic.generator = ICGEN_SONG;
#endif
#ifdef ICGEN_FALCONIC
		else if (par_string[0] == 'F' || par_string[0] == 'f')
			ic.generator = ICGEN_FALCONIC;
#endif
		else
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": IC generator not recognized!" << endl;
#ifdef LATFIELD2_HPP
			parallel.abortForce();
#endif
		}
	}
	else
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": IC generator not specified, selecting default (basic)" << endl;
		ic.generator = ICGEN_BASIC;
	}

	for (i = 0; i < MAX_PCL_SPECIES; i++)
		pptr[i] = ic.pclfile[i];
	if (ic.generator == ICGEN_READ_FROM_DISK)
	{
		if (!parseParameter(params, numparam, "particle file", pptr, i))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no particle file specified!" << endl;
#ifdef LATFIELD2_HPP
			parallel.abortForce();
#endif
		}
	}
	else if (!parseParameter(params, numparam, "template file", pptr, i))
	{
		COUT << COLORTEXT_RED << " warning: " << COLORTEXT_RESET << ": no template file specified!" << endl;
#ifdef LATFIELD2_HPP
#endif
	}

	for (; i < MAX_PCL_SPECIES; i++)
	{
		if (ic.generator == ICGEN_READ_FROM_DISK)
			strcpy(ic.pclfile[i], "/dev/null");
		else
			strcpy(ic.pclfile[i], ic.pclfile[i - 1]);
	}

#ifdef ICGEN_FALCONIC
	if ((!parseParameter(params, numparam, "mPk file", ic.pkfile) && !parseParameter(params, numparam, "Tk file", ic.tkfile) && ic.generator != ICGEN_READ_FROM_DISK && ic.generator != ICGEN_FALCONIC)
#else
	if ((!parseParameter(params, numparam, "mPk file", ic.pkfile) && !parseParameter(params, numparam, "Tk file", ic.tkfile) && ic.generator != ICGEN_READ_FROM_DISK)
#endif
#ifdef ICGEN_PREVOLUTION
		|| ic.generator == ICGEN_PREVOLUTION)
#else
	)
#endif
	{
#ifdef HAVE_CLASS
		COUT << " initial transfer functions will be computed by calling CLASS" << endl;
#else
		COUT << COLORTEXT_RED << " warning: " << COLORTEXT_RESET << ": no power spectrum file nor transfer function file specified!" << endl;
#ifdef LATFIELD2_HPP
#endif
#endif
	}

	if (parseParameter(params, numparam, "correct displacement", par_string))
	{
		if (par_string[0] == 'Y' || par_string[0] == 'y')
			ic.flags |= ICFLAG_CORRECT_DISPLACEMENT;
		else if (par_string[0] != 'N' && par_string[0] != 'n')
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": setting chosen for deconvolve displacement option not recognized, using default (no)" << endl;
	}

	if (parseParameter(params, numparam, "k-domain", par_string))
	{
		if (par_string[0] == 'S' || par_string[0] == 's')
			ic.flags |= ICFLAG_KSPHERE;
		else if (par_string[0] != 'C' && par_string[0] != 'c')
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": setting chosen for k-domain option not recognized, using default (cube)" << endl;
	}

	for (i = 0; i < MAX_PCL_SPECIES; i++)
		ic.numtile[i] = 0;

	if (!parseParameter(params, numparam, "tiling factor", ic.numtile, i) && ic.generator != ICGEN_READ_FROM_DISK)
	{
		if (parseParameter(params, numparam, "Ngrid", sim.numpts))
		{
			ic.numtile[0] = (int)(sim.numpts / 4 + 0.5);
			COUT << COLORTEXT_YELLOW << "Setting tiling factor from Ngrid\n"
				 << COLORTEXT_RESET;
		}
		else
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": tiling factor not specified, using default value for all species (1)" << endl;
			ic.numtile[0] = 1;
			i = 1;
		}
	}

	for (; i < MAX_PCL_SPECIES; i++)
		ic.numtile[i] = ic.numtile[i - 1];

	if (ic.numtile[0] <= 0 && ic.generator != ICGEN_READ_FROM_DISK)
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": tiling number for cdm particle template not set properly; using default value (1)" << endl;
		ic.numtile[0] = 1;
	}

	for (i = 1; i < MAX_PCL_SPECIES; i++)
	{
		if (ic.numtile[i] < 0 && ic.generator != ICGEN_READ_FROM_DISK)
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": tiling number for particle template not set properly; using default value (1)" << endl;
			ic.numtile[i] = 1;
		}
		else if (ic.generator == ICGEN_READ_FROM_DISK && strcmp(ic.pclfile[i], "/dev/null") != 0)
			ic.numtile[i] = 1;
	}

	if (ic.pkfile[0] != '\0')
	{
		sim.baryon_flag = 0;
		if (parseParameter(params, numparam, "baryon treatment", par_string))
		{
			if (par_string[0] != 'i' && par_string[0] != 'I')
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": using mPk file forces baryon treatment = ignore" << endl;
			}
		}
	}
	else if (parseParameter(params, numparam, "baryon treatment", par_string))
	{
		if (par_string[0] == 'i' || par_string[0] == 'I')
		{
			COUT << " baryon treatment set to: " << COLORTEXT_CYAN << "ignore" << COLORTEXT_RESET << endl;
			sim.baryon_flag = 0;
		}
		else if (par_string[0] == 's' || par_string[0] == 'S')
		{
			COUT << " baryon treatment set to: " << COLORTEXT_CYAN << "sample" << COLORTEXT_RESET << endl;
			sim.baryon_flag = 1;
		}
		else if (par_string[0] == 'b' || par_string[0] == 'B')
		{
			COUT << " baryon treatment set to: " << COLORTEXT_CYAN << "blend" << COLORTEXT_RESET << endl;
			sim.baryon_flag = 2;
		}
		else if (par_string[0] == 'h' || par_string[0] == 'H')
		{
			COUT << " baryon treatment set to: " << COLORTEXT_CYAN << "hybrid" << COLORTEXT_RESET << endl;
			sim.baryon_flag = 3;
		}
		else
		{
			COUT << COLORTEXT_RED << " warning: " << COLORTEXT_RESET << ": baryon treatment not supported!" << endl;
#ifdef LATFIELD2_HPP
#endif
		}
	}
	else if (ic.generator == ICGEN_READ_FROM_DISK)
	{
		sim.baryon_flag = 0;
	}
	else
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": baryon treatment not specified, using default (blend)" << endl;
		sim.baryon_flag = 2;
	}

	if (sim.baryon_flag == 1 && ic.numtile[1] <= 0 && ic.generator != ICGEN_READ_FROM_DISK)
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": tiling number for baryon particle template not set properly; using default value (1)" << endl;
		ic.numtile[1] = 1;
	}

	if (parseParameter(params, numparam, "radiation treatment", par_string))
	{
		if (par_string[0] == 'i' || par_string[0] == 'I' || par_string[0] == 'b' || par_string[0] == 'B')
		{
			sim.radiation_flag = 0;
			COUT << " radiation treatment set to: " << COLORTEXT_CYAN << "background" << COLORTEXT_RESET << endl;
		}
#ifdef HAVE_CLASS
		else if (par_string[0] == 'c' || par_string[0] == 'C')
		{
			sim.radiation_flag = 1;
			COUT << " radiation treatment set to: " << COLORTEXT_CYAN << "CLASS" << COLORTEXT_RESET << endl;
			if ((ic.pkfile[0] != '\0' || ic.tkfile[0] != '\0')
#ifdef ICGEN_PREVOLUTION
				&& ic.generator != ICGEN_PREVOLUTION
#endif
			)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": using radiation treatment = CLASS and providing initial power spectra / transfer functions independently" << endl;
				COUT << "              is dangerous! In order to ensure consistency, it is recommended to call CLASS directly." << endl;
			}
		}
		else
		{
			COUT << COLORTEXT_RED << " warning: " << COLORTEXT_RESET << ": radiation treatment not supported!" << endl;
#ifdef LATFIELD2_HPP
#endif
		}
#else
		else
		{
			sim.radiation_flag = 0;
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": CLASS is not available, setting radiation treatment = background" << endl;
		}
#endif
	}
	else
		sim.radiation_flag = 0;

	if (parseParameter(params, numparam, "fluid treatment", par_string))
	{
		if (par_string[0] == 'b' || par_string[0] == 'B')
		{
			sim.fluid_flag = 0;
			COUT << " fluid treatment set to: " << COLORTEXT_CYAN << "background" << COLORTEXT_RESET << endl;
		}
#ifdef HAVE_CLASS
		else if (par_string[0] == 'c' || par_string[0] == 'C')
		{
			sim.fluid_flag = 1;
			COUT << " fluid treatment set to: " << COLORTEXT_CYAN << "CLASS" << COLORTEXT_RESET << endl;
			if ((ic.pkfile[0] != '\0' || ic.tkfile[0] != '\0')
#ifdef ICGEN_PREVOLUTION
				&& ic.generator != ICGEN_PREVOLUTION
#endif
			)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": using fluid treatment = CLASS and providing initial power spectra / transfer functions independently" << endl;
				COUT << "              is dangerous! In order to ensure consistency, it is recommended to call CLASS directly." << endl;
			}
		}
		else
		{
			COUT << COLORTEXT_RED << " warning: " << COLORTEXT_RESET << ": fluid treatment not supported!" << endl;
#ifdef LATFIELD2_HPP
#endif
		}
#else
		else
		{
			sim.fluid_flag = 0;
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": CLASS is not available, setting fluid treatment = background" << endl;
		}
#endif
	}
	else
		sim.fluid_flag = 0;

	parseParameter(params, numparam, "relaxation redshift", ic.z_relax);

	if (ic.generator == ICGEN_READ_FROM_DISK)
	{
		parseParameter(params, numparam, "restart redshift", ic.z_ic);
		if (!parseParameter(params, numparam, "cycle", ic.restart_cycle))
		{
			if (sim.radiation_flag)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": using radiation treatment = CLASS and IC generator = read from disk does not guarantee" << endl;
				COUT << "              that the realization of the radiation perturbations is consistent!" << endl;
			}
		}
		parseParameter(params, numparam, "tau", ic.restart_tau);
		parseParameter(params, numparam, "dtau", ic.restart_dtau);
		for (i = 0; i < 3; i++)
			pptr[i] = ic.metricfile[i];
		parseParameter(params, numparam, "metric file", pptr, i);
		if (parseParameter(params, numparam, "gevolution version", ic.restart_version))
		{
			if (ic.restart_version - GEVOLUTION_VERSION > 0.0001)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": version number of settings file (" << ic.restart_version << ") is higher than version of executable (" << GEVOLUTION_VERSION << ")!" << endl;
			}
		}
	}
#ifdef ICGEN_PREVOLUTION
	else if (ic.generator == ICGEN_PREVOLUTION)
	{
		if (!parseParameter(params, numparam, "prevolution redshift", ic.z_ic))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no starting redshift specified for IC generator = prevolution" << endl;
#ifdef LATFIELD2_HPP
			parallel.abortForce();
#endif
		}

		parseParameter(params, numparam, "prevolution Courant factor", ic.Cf);
	}
#endif

	if (!parseParameter(params, numparam, "A_s", ic.A_s) && (
#ifdef ICGEN_FALCONIC
																ic.generator == ICGEN_FALCONIC ||
#endif
#ifdef ICGEN_PREVOLUTION
																ic.generator == ICGEN_PREVOLUTION ||
#endif
																sim.radiation_flag > 0 || (ic.pkfile[0] == '\0' && ic.generator != ICGEN_READ_FROM_DISK)))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": power spectrum normalization not specified, using default value (2.215e-9)" << endl;
	}

	if (!parseParameter(params, numparam, "n_s", ic.n_s) && (
#ifdef ICGEN_FALCONIC
																ic.generator == ICGEN_FALCONIC ||
#endif
#ifdef ICGEN_PREVOLUTION
																ic.generator == ICGEN_PREVOLUTION ||
#endif
																sim.radiation_flag > 0 || (ic.pkfile[0] == '\0' && ic.generator != ICGEN_READ_FROM_DISK)))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": scalar spectral index not specified, using default value (0.9619)" << endl;
	}

	if (!parseParameter(params, numparam, "k_pivot", ic.k_pivot) && (
#ifdef ICGEN_FALCONIC
																		ic.generator == ICGEN_FALCONIC ||
#endif
#ifdef ICGEN_PREVOLUTION
																		ic.generator == ICGEN_PREVOLUTION ||
#endif
																		sim.radiation_flag > 0 || (ic.pkfile[0] == '\0' && ic.generator != ICGEN_READ_FROM_DISK)))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": pivot scale not specified, using default value (0.05 / Mpc)" << endl;
	}

	// parse metadata

	sim.numpts = 0;
	sim.downgrade_factor = 1;
	sim.downgrade_factor2 = 1;
	sim.animation_downgrade_factor = 1;
	for (i = 0; i < MAX_PCL_SPECIES; i++)
		sim.numpcl[i] = 0;
	sim.vector_flag = VECTOR_PARABOLIC;
	sim.gr_flag = 0;
	sim.out_pk = 0;
	sim.out_snapshot = 0;
	sim.out_snapshot2 = 0;
	sim.out_animation = 0;
	sim.out_lightcone[0] = 0;
	sim.num_pk = MAX_OUTPUTS;
	sim.numbins = 0;
	sim.num_snapshot = MAX_OUTPUTS;
	sim.num_snapshot2 = MAX_OUTPUTS;
	sim.num_lightcone = 0;
	sim.num_restart = MAX_OUTPUTS;
	for (i = 0; i < MAX_PCL_SPECIES; i++)
		sim.tracer_factor[i] = 1;
	sim.Cf = 1.;
	sim.steplimit = 1.;
	sim.boxsize = -1.;
	sim.wallclocklimit = -1.;
	sim.z_in = 0.;

	if (parseParameter(params, numparam, "vector method", par_string))
	{
		if (par_string[0] == 'p' || par_string[0] == 'P')
		{
			COUT << " vector method set to: " << COLORTEXT_CYAN << "parabolic" << COLORTEXT_RESET << endl;
			sim.vector_flag = VECTOR_PARABOLIC;
		}
		else if (par_string[0] == 'e' || par_string[0] == 'E')
		{
			COUT << " vector method set to: " << COLORTEXT_CYAN << "elliptic" << COLORTEXT_RESET << endl;
			sim.vector_flag = VECTOR_ELLIPTIC;
		}
		else
		{
			COUT << COLORTEXT_RED << " warning: " << COLORTEXT_RESET << ": vector method not supported!" << endl;
#ifdef LATFIELD2_HPP
#endif
		}
	}

	if (!parseParameter(params, numparam, "generic file base", sim.basename_generic))
		sim.basename_generic[0] = '\0';

	if (!parseParameter(params, numparam, "snapshot file base", sim.basename_snapshot))
		strcpy(sim.basename_snapshot, "snapshot");

	if (!parseParameter(params, numparam, "animation file base", sim.basename_animation))
		strcpy(sim.basename_animation, "animation");

	if (!parseParameter(params, numparam, "Pk file base", sim.basename_pk))
		strcpy(sim.basename_pk, "pk");

	if (!parseParameter(params, numparam, "lightcone file base", sim.basename_lightcone))
		strcpy(sim.basename_lightcone, "lightcone");

	if (!parseParameter(params, numparam, "output path", sim.output_path))
		sim.output_path[0] = '\0';

	if (!parseParameter(params, numparam, "hibernation path", sim.restart_path))
		strcpy(sim.restart_path, sim.output_path);

	if (!parseParameter(params, numparam, "hibernation file base", sim.basename_restart))
		strcpy(sim.basename_restart, "restart");

	parseParameter(params, numparam, "Ngrid", sim.numpts);
	if (sim.numpts < 2 || !isfinite(sim.numpts))
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": number of grid points not set properly!" << endl;
#ifdef LATFIELD2_HPP
		parallel.abortForce();
#endif
	}
	if (!parseParameter(params, numparam, "boxsize", sim.boxsize))
	{
		if (parseParameter(params, numparam, "dx", sim.dx))
		{
			COUT << COLORTEXT_YELLOW << "Setting boxsize from dx\n"
				 << COLORTEXT_RESET;
			sim.boxsize = sim.numpts / sim.dx;
		}
	}
	if (!parseParameter(params, numparam, "dx", sim.dx))
	{
		sim.dx = sim.boxsize / sim.numpts;
	}
	if (sim.boxsize <= 0. || !isfinite(sim.boxsize))
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": simulation box size not set properly!" << endl;
#ifdef LATFIELD2_HPP
		parallel.abortForce();
#endif
	}

	if (parseParameter(params, numparam, "downgrade factor", sim.downgrade_factor))
	{
		if (sim.downgrade_factor < 1 || sim.downgrade_factor >= sim.numpts || !isfinite(sim.downgrade_factor))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": downgrade factor makes no sense!" << endl;
#ifdef LATFIELD2_HPP
			parallel.abortForce();
		}
		if (sim.downgrade_factor > 1 && (sim.numpts % parallel.grid_size()[0] || sim.numpts % parallel.grid_size()[1] || (sim.numpts / parallel.grid_size()[0]) % sim.downgrade_factor || (sim.numpts / parallel.grid_size()[1]) % sim.downgrade_factor))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": downgrade factor does not appear to be compatible with process layout at given Ngrid!" << endl;
			parallel.abortForce();
#endif
		}
	}
	if (parseParameter(params, numparam, "downgrade factor 2", sim.downgrade_factor2))
	{
		if (sim.downgrade_factor2 < 1 || sim.downgrade_factor2 >= sim.numpts || !isfinite(sim.downgrade_factor2))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": downgrade factor2 makes no sense!" << endl;
#ifdef LATFIELD2_HPP
			parallel.abortForce();
		}
		if (sim.downgrade_factor2 > 1 && (sim.numpts % parallel.grid_size()[0] || sim.numpts % parallel.grid_size()[1] || (sim.numpts / parallel.grid_size()[0]) % sim.downgrade_factor2 || (sim.numpts / parallel.grid_size()[1]) % sim.downgrade_factor2))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": downgrade factor2 does not appear to be compatible with process layout at given Ngrid!" << endl;
			parallel.abortForce();
#endif
		}
	}
	if (parseParameter(params, numparam, "animation downgrade factor", sim.animation_downgrade_factor))
	{
		if (sim.animation_downgrade_factor < 1 || sim.animation_downgrade_factor >= sim.numpts || !isfinite(sim.animation_downgrade_factor))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": animation downgrade factor makes no sense!" << endl;
#ifdef LATFIELD2_HPP
			parallel.abortForce();
		}
		if (sim.animation_downgrade_factor > 1 && (sim.numpts % parallel.grid_size()[0] || sim.numpts % parallel.grid_size()[1] || (sim.numpts / parallel.grid_size()[0]) % sim.animation_downgrade_factor || (sim.numpts / parallel.grid_size()[1]) % sim.animation_downgrade_factor))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": animation downgrade factor does not appear to be compatible with process layout at given Ngrid!" << endl;
			parallel.abortForce();
#endif
		}
	}

	if (!parseParameter(params, numparam, "Courant factor", sim.Cf))
	{
		double temp;
		if (parseParameter(params, numparam, "dtaucdm", temp))
		{
			sim.Cf = temp * sim.numpts / sim.boxsize;
		}
		else
		{
			COUT << COLORTEXT_RED << "Did not specify time resolution!\n"
				 << COLORTEXT_RESET;
		}
	}

	if (ic.Cf < 0.)
		ic.Cf = sim.Cf;

	parseParameter(params, numparam, "time step limit", sim.steplimit);

	if (!parseParameter(params, numparam, "move limit", sim.movelimit))
		sim.movelimit = (double)sim.numpts;

	if (!parseParameter(params, numparam, "initial redshift", sim.z_in))
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": initial redshift not specified!" << endl;
#ifdef LATFIELD2_HPP
		parallel.abortForce();
#endif
	}

	if (ic.z_relax < -1.)
		ic.z_relax = sim.z_in;
#ifdef ICGEN_PREVOLUTION
	else if (ic.generator == ICGEN_PREVOLUTION && ic.z_relax < sim.z_in)
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": relaxation redshift cannot be below initial redshift for IC generator = prevolution; reset to initial redshift!" << endl;
		ic.z_relax = sim.z_in;
	}
#endif

	// if (ic.z_ic < sim.z_in && ic.generator != ICGEN_READ_FROM_DISK)
	ic.z_ic = sim.z_in;

	parseParameter(params, numparam, "snapshot redshifts", sim.z_snapshot, sim.num_snapshot);
	if (sim.num_snapshot > 0)
	{
		int tmpnum = sim.num_snapshot;
		for (int i = 0; i < sim.num_snapshot; i++)
		{
			if (sim.z_snapshot[i] > ic.z_ic)
			{
				sim.num_snapshot--;
				tmpnum++;
				COUT << "zic = " << ic.z_ic << ": ";
				COUT << "Removing snapshot output at z = " << sim.z_snapshot[i] << endl;
				sim.z_snapshot[i] = -1.;
			}
		}
		qsort((void *)sim.z_snapshot, (size_t)tmpnum, sizeof(double), sort_descending);
	}
	parseParameter(params, numparam, "snapshot redshifts 2", sim.z_snapshot2, sim.num_snapshot2);
	if (sim.num_snapshot2 > 0)
	{
		int tmpnum = sim.num_snapshot2;
		for (int i = 0; i < sim.num_snapshot2; i++)
		{
			if (sim.z_snapshot2[i] > ic.z_ic)
			{
				sim.num_snapshot2--;
				tmpnum++;
				COUT << "zic = " << ic.z_ic << ": ";
				COUT << "Removing snapshot output at z = " << sim.z_snapshot2[i] << endl;
				sim.z_snapshot2[i] = -1.;
			}
		}
		qsort((void *)sim.z_snapshot2, (size_t)tmpnum, sizeof(double), sort_descending);
	}
	parseParameter(params, numparam, "Pk redshifts", sim.z_pk, sim.num_pk);
	if (sim.num_pk > 0)
	{
		int tmpnum = sim.num_pk;
		for (int i = 0; i < sim.num_pk; i++)
		{
			if (sim.z_pk[i] > ic.z_ic)
			{
				sim.num_pk--;
				tmpnum++;
				COUT << "Removing spectra output at z = " << sim.z_pk[i] << endl;
				sim.z_pk[i] = -1.;
			}
		}
		qsort((void *)sim.z_pk, (size_t)tmpnum, sizeof(double), sort_descending);
	}

	parseParameter(params, numparam, "hibernation redshifts", sim.z_restart, sim.num_restart);
	if (sim.num_restart > 0)
		qsort((void *)sim.z_restart, (size_t)sim.num_restart, sizeof(double), sort_descending);

	parseParameter(params, numparam, "hibernation wallclock limit", sim.wallclocklimit);

	parseFieldSpecifiers(params, numparam, "lightcone outputs", sim.out_lightcone[0]);
	parseFieldSpecifiers(params, numparam, "snapshot outputs", sim.out_snapshot);
	parseFieldSpecifiers(params, numparam, "animation outputs", sim.out_animation);
	parseFieldSpecifiers(params, numparam, "snapshot outputs 2", sim.out_snapshot2);
	parseFieldSpecifiers(params, numparam, "Pk outputs", sim.out_pk);

	if (!parseParameter(params, numparam, "lightcone pixel factor", sim.pixelfactor[0]))
		sim.pixelfactor[0] = 0.5;

	if (!parseParameter(params, numparam, "lightcone shell factor", sim.shellfactor[0]))
		sim.shellfactor[0] = 1.;

	if (!parseParameter(params, numparam, "lightcone covering", sim.covering[0]))
		sim.covering[0] = 2. + 4. / sim.Cf / sim.shellfactor[0];

	i = 2;
	if (parseParameter(params, numparam, "lightcone Nside", sim.Nside[0], i))
	{
		if (i < 1)
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": parsing of lightcone Nside parameter failed, assuming minimum value Nside=2" << endl;
			sim.Nside[0][0] = 2;
		}
		if (i < 2)
			sim.Nside[0][1] = sim.Nside[0][0];
		if (sim.Nside[0][0] < 2)
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": lightcone Nside parameter out of bounds, assuming minimum value Nside=2" << endl;
			sim.Nside[0][0] = 2;
		}
		else
		{
			for (i = 2; i < sim.Nside[0][0]; i *= 2)
				;
			if (i != sim.Nside[0][0])
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": lightcone Nside parameter was found to be no power of two, assuming nearest value Nside=" << i << endl;
				sim.Nside[0][0] = i;
			}
		}
		if (sim.Nside[0][1] < sim.Nside[0][0])
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": lightcone Nside parameter out of bounds, assuming minimum value Nside=" << sim.Nside[0][0] << endl;
			sim.Nside[0][1] = sim.Nside[0][0];
		}
		else
		{
			for (i = 2; i < sim.Nside[0][1]; i *= 2)
				;
			if (i != sim.Nside[0][1])
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": lightcone Nside parameter was found to be no power of two, assuming nearest value Nside=" << i << endl;
				sim.Nside[0][1] = i;
			}
		}
	}
	else
	{
		sim.Nside[0][0] = 2;
		for (sim.Nside[0][1] = 2; sim.Nside[0][1] < sim.numpts; sim.Nside[0][1] *= 2)
			;
	}

	for (i = 1; i < MAX_OUTPUTS; i++)
	{
		sim.out_lightcone[i] = sim.out_lightcone[0];
		sim.Nside[i][0] = sim.Nside[0][0];
		sim.Nside[i][1] = sim.Nside[0][1];
		sim.pixelfactor[i] = sim.pixelfactor[0];
		sim.shellfactor[i] = sim.shellfactor[0];
		sim.covering[i] = sim.covering[0];
	}

	i = 3;
	if (parseParameter(params, numparam, "lightcone vertex", sim.lightcone[0].vertex, i))
	{
		if (i != 3 || sim.lightcone[0].vertex[0] < 0. || sim.lightcone[0].vertex[0] >= sim.boxsize || sim.lightcone[0].vertex[1] < 0. || sim.lightcone[0].vertex[1] >= sim.boxsize || sim.lightcone[0].vertex[2] < 0. || sim.lightcone[0].vertex[2] >= sim.boxsize)
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": parsing of lightcone geometry failed, no lightcone will be written!" << endl;
		}
		else
		{
			sim.num_lightcone = 1;
			for (i = 0; i < 3; i++)
				sim.lightcone[0].vertex[i] /= sim.boxsize;

			if (parseParameter(params, numparam, "lightcone opening half-angle", sim.lightcone[0].opening))
			{
				if (sim.lightcone[0].opening > 180. || sim.lightcone[0].opening <= 0.)
				{
					COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": lightcone opening half-angle out of bounds, assuming full sky" << endl;
					sim.lightcone[0].opening = -1.;
				}
				else
					sim.lightcone[0].opening = cos(sim.lightcone[0].opening * M_PI / 180.);
			}
			else
				sim.lightcone[0].opening = -1.;

			i = 2;
			if (parseParameter(params, numparam, "lightcone distance", sim.lightcone[0].distance, i))
			{
				sim.lightcone[0].distance[0] /= sim.boxsize;
				if (i > 1)
				{
					sim.lightcone[0].distance[1] /= sim.boxsize;
					qsort((void *)sim.lightcone[0].distance, 2, sizeof(double), sort_descending);
				}
				else
					sim.lightcone[0].distance[1] = 0.;
			}
			else
			{
				sim.lightcone[0].distance[0] = 0.5;
				sim.lightcone[0].distance[1] = 0.;
			}

			if (!parseParameter(params, numparam, "lightcone redshift", sim.lightcone[0].z))
				sim.lightcone[0].z = 0.;

			i = 3;
			if (parseParameter(params, numparam, "lightcone direction", sim.lightcone[0].direction, i))
			{
				if (i == 2)
				{
					tmp = sin(sim.lightcone[0].direction[0] * M_PI / 180.);
					sim.lightcone[0].direction[2] = cos(sim.lightcone[0].direction[0] * M_PI / 180.);
					sim.lightcone[0].direction[0] = tmp * cos(sim.lightcone[0].direction[1] * M_PI / 180.);
					sim.lightcone[0].direction[1] = tmp * sin(sim.lightcone[0].direction[1] * M_PI / 180.);
				}
				else if (i == 3)
				{
					tmp = sqrt(sim.lightcone[0].direction[0] * sim.lightcone[0].direction[0] + sim.lightcone[0].direction[1] * sim.lightcone[0].direction[1] + sim.lightcone[0].direction[2] * sim.lightcone[0].direction[2]);
					sim.lightcone[0].direction[0] /= tmp;
					sim.lightcone[0].direction[1] /= tmp;
					sim.lightcone[0].direction[2] /= tmp;
				}
				else
				{
					COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": parsing of lightcone geometry failed, no lightcone will be written!" << endl;
					sim.num_lightcone = 0;
				}
			}
			else
			{
				sim.lightcone[0].direction[0] = 0.;
				sim.lightcone[0].direction[1] = 0.;
				sim.lightcone[0].direction[2] = 1.;
			}
		}
	}
	else
	{
		for (sim.num_lightcone = 0; sim.num_lightcone < MAX_OUTPUTS; sim.num_lightcone++)
		{
			snprintf(par_string, SIZPARSTRING, "lightcone %d vertex", sim.num_lightcone);
			i = 3;
			if (parseParameter(params, numparam, par_string, sim.lightcone[sim.num_lightcone].vertex, i))
			{
				if (i != 3 || sim.lightcone[sim.num_lightcone].vertex[0] < 0. || sim.lightcone[sim.num_lightcone].vertex[0] >= sim.boxsize || sim.lightcone[sim.num_lightcone].vertex[1] < 0. || sim.lightcone[sim.num_lightcone].vertex[1] >= sim.boxsize || sim.lightcone[sim.num_lightcone].vertex[2] < 0. || sim.lightcone[sim.num_lightcone].vertex[2] >= sim.boxsize)
				{
					COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": parsing of lightcone geometry failed, not all lightcones will be written!" << endl;
					break;
				}
				else
				{
					for (i = 0; i < 3; i++)
						sim.lightcone[sim.num_lightcone].vertex[i] /= sim.boxsize;

					snprintf(par_string, SIZPARSTRING, "lightcone %d outputs", sim.num_lightcone);
					parseFieldSpecifiers(params, numparam, par_string, sim.out_lightcone[sim.num_lightcone]);

					snprintf(par_string, SIZPARSTRING, "lightcone %d opening half-angle", sim.num_lightcone);
					if (parseParameter(params, numparam, par_string, sim.lightcone[sim.num_lightcone].opening))
					{
						if (sim.lightcone[sim.num_lightcone].opening > 180. || sim.lightcone[sim.num_lightcone].opening <= 0.)
						{
							COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": lightcone opening half-angle out of bounds, assuming full sky" << endl;
							sim.lightcone[sim.num_lightcone].opening = -1.;
						}
						else
							sim.lightcone[sim.num_lightcone].opening = cos(sim.lightcone[sim.num_lightcone].opening * M_PI / 180.);
					}
					else
						sim.lightcone[sim.num_lightcone].opening = -1.;

					snprintf(par_string, SIZPARSTRING, "lightcone %d distance", sim.num_lightcone);
					i = 2;
					if (parseParameter(params, numparam, par_string, sim.lightcone[sim.num_lightcone].distance, i))
					{
						sim.lightcone[sim.num_lightcone].distance[0] /= sim.boxsize;
						if (i > 1)
						{
							sim.lightcone[sim.num_lightcone].distance[1] /= sim.boxsize;
							qsort((void *)sim.lightcone[sim.num_lightcone].distance, 2, sizeof(double), sort_descending);
						}
						else
							sim.lightcone[sim.num_lightcone].distance[1] = 0.;
					}
					else
					{
						sim.lightcone[sim.num_lightcone].distance[0] = 0.5;
						sim.lightcone[sim.num_lightcone].distance[1] = 0.;
					}

					snprintf(par_string, SIZPARSTRING, "lightcone %d redshift", sim.num_lightcone);
					if (!parseParameter(params, numparam, par_string, sim.lightcone[sim.num_lightcone].z))
						sim.lightcone[sim.num_lightcone].z = 0.;

					snprintf(par_string, SIZPARSTRING, "lightcone %d direction", sim.num_lightcone);
					i = 3;
					if (parseParameter(params, numparam, par_string, sim.lightcone[sim.num_lightcone].direction, i))
					{
						if (i == 2)
						{
							tmp = sin(sim.lightcone[sim.num_lightcone].direction[0] * M_PI / 180.);
							sim.lightcone[sim.num_lightcone].direction[2] = cos(sim.lightcone[sim.num_lightcone].direction[0] * M_PI / 180.);
							sim.lightcone[sim.num_lightcone].direction[0] = tmp * cos(sim.lightcone[sim.num_lightcone].direction[1] * M_PI / 180.);
							sim.lightcone[sim.num_lightcone].direction[1] = tmp * sin(sim.lightcone[sim.num_lightcone].direction[1] * M_PI / 180.);
						}
						else if (i == 3)
						{
							tmp = sqrt(sim.lightcone[sim.num_lightcone].direction[0] * sim.lightcone[sim.num_lightcone].direction[0] + sim.lightcone[sim.num_lightcone].direction[1] * sim.lightcone[sim.num_lightcone].direction[1] + sim.lightcone[sim.num_lightcone].direction[2] * sim.lightcone[sim.num_lightcone].direction[2]);
							sim.lightcone[sim.num_lightcone].direction[0] /= tmp;
							sim.lightcone[sim.num_lightcone].direction[1] /= tmp;
							sim.lightcone[sim.num_lightcone].direction[2] /= tmp;
						}
						else
						{
							COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": parsing of lightcone geometry failed, not all lightcones will be written!" << endl;
							break;
						}
					}
					else
					{
						sim.lightcone[sim.num_lightcone].direction[0] = 0.;
						sim.lightcone[sim.num_lightcone].direction[1] = 0.;
						sim.lightcone[sim.num_lightcone].direction[2] = 1.;
					}
					snprintf(par_string, SIZPARSTRING, "lightcone %d pixel factor", sim.num_lightcone);
					parseParameter(params, numparam, par_string, sim.pixelfactor[sim.num_lightcone]);
					snprintf(par_string, SIZPARSTRING, "lightcone %d shell factor", sim.num_lightcone);
					parseParameter(params, numparam, par_string, sim.shellfactor[sim.num_lightcone]);
					snprintf(par_string, SIZPARSTRING, "lightcone %d covering", sim.num_lightcone);
					parseParameter(params, numparam, par_string, sim.covering[sim.num_lightcone]);
					snprintf(par_string, SIZPARSTRING, "lightcone %d Nside", sim.num_lightcone);

					i = 2;
					if (parseParameter(params, numparam, par_string, sim.Nside[sim.num_lightcone], i))
					{
						if (i < 1)
						{
							COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": parsing of lightcone Nside parameter failed, assuming minimum value Nside=2" << endl;
							sim.Nside[sim.num_lightcone][0] = 2;
						}
						if (i < 2)
							sim.Nside[sim.num_lightcone][1] = sim.Nside[sim.num_lightcone][0];
						if (sim.Nside[sim.num_lightcone][0] < 2)
						{
							COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": lightcone Nside parameter out of bounds, assuming minimum value Nside=2" << endl;
							sim.Nside[sim.num_lightcone][0] = 2;
						}
						else
						{
							for (i = 2; i < sim.Nside[sim.num_lightcone][0]; i *= 2)
								;
							if (i != sim.Nside[sim.num_lightcone][0])
							{
								COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": lightcone Nside parameter was found to be no power of two, assuming nearest value Nside=" << i << endl;
								sim.Nside[sim.num_lightcone][0] = i;
							}
						}
						if (sim.Nside[sim.num_lightcone][1] < sim.Nside[sim.num_lightcone][0])
						{
							COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": lightcone Nside parameter out of bounds, assuming minimum value Nside=" << sim.Nside[sim.num_lightcone][0] << endl;
							sim.Nside[sim.num_lightcone][1] = sim.Nside[sim.num_lightcone][0];
						}
						else
						{
							for (i = 2; i < sim.Nside[sim.num_lightcone][1]; i *= 2)
								;
							if (i != sim.Nside[sim.num_lightcone][1])
							{
								COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": lightcone Nside parameter was found to be no power of two, assuming nearest value Nside=" << i << endl;
								sim.Nside[sim.num_lightcone][1] = i;
							}
						}
					}
				}
			}
			else
				break;
		}
	}

	i = MAX_PCL_SPECIES;
	parseParameter(params, numparam, "tracer factor", sim.tracer_factor, i);
	for (; i > 0; i--)
	{
		if (sim.tracer_factor[i - 1] < 0)
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": tracer factor not set properly; using default value (1)" << endl;
			sim.tracer_factor[i - 1] = 1;
		}
	}

	if ((sim.num_snapshot <= 0 || sim.out_snapshot == 0) && (sim.num_snapshot2 <= 0 || sim.out_snapshot2 == 0) && (sim.num_pk <= 0 || sim.out_pk == 0) && (sim.out_animation == 0 || sim.out_animation <= 0))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": no output specified!" << endl;
	}

	if (!parseParameter(params, numparam, "Pk bins", sim.numbins))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of Pk bins not set properly; using default value (64)" << endl;
		sim.numbins = 64;
	}

	if (parseParameter(params, numparam, "gravity theory", par_string))
	{
		if (par_string[0] == 'N' || par_string[0] == 'n')
		{
			COUT << " gravity theory set to: " << COLORTEXT_CYAN << "Newtonian" << COLORTEXT_RESET << endl;
			sim.gr_flag = 0;
			if (ic.pkfile[0] == '\0' && ic.tkfile[0] != '\0'
#ifdef ICGEN_PREVOLUTION
				&& ic.generator != ICGEN_PREVOLUTION
#endif
			)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": gauge transformation to N-body gauge can only be performed for the positions; the transformation for" << endl;
				COUT << "              the velocities requires time derivatives of transfer functions. Call CLASS directly to avoid this issue." << endl;
			}
		}
		else if (par_string[0] == 'G' || par_string[0] == 'g')
		{
			COUT << " gravity theory set to: " << COLORTEXT_CYAN << "General Relativity" << COLORTEXT_RESET << endl;
			sim.gr_flag = 1;
		}
		else
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": gravity theory unknown, using default (General Relativity)" << endl;
			sim.gr_flag = 1;
		}
	}
	else
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": gravity theory not selected, using default (General Relativity)" << endl;
		sim.gr_flag = 1;
	}

	// parse cosmological parameters

	if (!parseParameter(params, numparam, "h", cosmo.h))
	{
		cosmo.h = P_HUBBLE;
	}

	cosmo.num_ncdm = MAX_PCL_SPECIES - 2;
	if (!parseParameter(params, numparam, "m_ncdm", cosmo.m_ncdm, cosmo.num_ncdm))
	{
		for (i = 0; i < MAX_PCL_SPECIES - 2; i++)
			cosmo.m_ncdm[i] = 0.;
		cosmo.num_ncdm = 0;
	}

	if (parseParameter(params, numparam, "N_ncdm", i))
	{
		if (i < 0 || !isfinite(i))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": number of ncdm species not set properly!" << endl;
#ifdef LATFIELD2_HPP
			parallel.abortForce();
#endif
		}
		if (i > cosmo.num_ncdm)
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": N_ncdm = " << i << " is larger than the number of mass parameters specified (" << cosmo.num_ncdm << ")!" << endl;
#ifdef LATFIELD2_HPP
			parallel.abortForce();
#endif
		}
		cosmo.num_ncdm = i;
	}
	else if (cosmo.num_ncdm > 0)
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": N_ncdm not specified, inferring from number of mass parameters in m_ncdm (" << cosmo.num_ncdm << ")!" << endl;
	}

	for (i = 0; i < MAX_PCL_SPECIES - 2; i++)
	{
		cosmo.T_ncdm[i] = P_T_NCDM;
		cosmo.deg_ncdm[i] = 1.0;
	}
	parseParameter(params, numparam, "T_ncdm", cosmo.T_ncdm, i);
	i = MAX_PCL_SPECIES - 2;
	parseParameter(params, numparam, "deg_ncdm", cosmo.deg_ncdm, i);

	for (i = 0; i < cosmo.num_ncdm; i++)
	{
		cosmo.Omega_ncdm[i] = cosmo.m_ncdm[i] * cosmo.deg_ncdm[i] / P_NCDM_MASS_OMEGA / cosmo.h / cosmo.h;
	}

	if (parseParameter(params, numparam, "T_cmb", cosmo.Omega_g))
	{
		cosmo.Omega_g = cosmo.Omega_g * cosmo.Omega_g / cosmo.h;
		cosmo.Omega_g = cosmo.Omega_g * cosmo.Omega_g * C_PLANCK_LAW; // Planck's law
	}
	else if (parseParameter(params, numparam, "omega_g", cosmo.Omega_g))
	{
		cosmo.Omega_g /= cosmo.h * cosmo.h;
	}
	else if (!parseParameter(params, numparam, "Omega_g", cosmo.Omega_g))
	{
		cosmo.Omega_g = 0.;
	}

	if (parseParameter(params, numparam, "N_ur", cosmo.Omega_ur))
	{
		cosmo.Omega_ur *= (7. / 8.) * pow(4. / 11., 4. / 3.) * cosmo.Omega_g;
	}
	else if (parseParameter(params, numparam, "N_eff", cosmo.Omega_ur))
	{
		cosmo.Omega_ur *= (7. / 8.) * pow(4. / 11., 4. / 3.) * cosmo.Omega_g;
	}
	else if (parseParameter(params, numparam, "omega_ur", cosmo.Omega_ur))
	{
		cosmo.Omega_ur /= cosmo.h * cosmo.h;
	}
	else if (!parseParameter(params, numparam, "Omega_ur", cosmo.Omega_ur))
	{
		cosmo.Omega_ur = P_N_UR * (7. / 8.) * pow(4. / 11., 4. / 3.) * cosmo.Omega_g;
	}

	cosmo.Omega_rad = cosmo.Omega_g + cosmo.Omega_ur;

	if (parseParameter(params, numparam, "omega_fld", cosmo.Omega_fld))
		cosmo.Omega_fld /= cosmo.h * cosmo.h;
	else if (!parseParameter(params, numparam, "Omega_fld", cosmo.Omega_fld))
		cosmo.Omega_fld = 0.;
	if (!parseParameter(params, numparam, "w0_fld", cosmo.w0_fld))
		cosmo.w0_fld = -1.;
	if (!parseParameter(params, numparam, "wa_fld", cosmo.wa_fld))
		cosmo.wa_fld = 0.;
	if (!parseParameter(params, numparam, "cs2_fld", cosmo.cs2_fld))
		cosmo.cs2_fld = 1.;

	if (cosmo.Omega_fld > 0 && cosmo.w0_fld == -1.)
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": w0_fld = -1 is singular, setting Omega_fld = 0." << endl;
		cosmo.Omega_fld = 0.;
	}

	if (parseParameter(params, numparam, "omega_b", cosmo.Omega_b))
	{
		cosmo.Omega_b /= cosmo.h * cosmo.h;
	}
	else if (!parseParameter(params, numparam, "Omega_b", cosmo.Omega_b))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": Omega_b not found in settings file, setting to default (0)." << endl;
		cosmo.Omega_b = 0.;
	}

	if (parseParameter(params, numparam, "omega_cdm", cosmo.Omega_cdm))
	{
		cosmo.Omega_cdm /= cosmo.h * cosmo.h;
	}
	else if (!parseParameter(params, numparam, "Omega_cdm", cosmo.Omega_cdm))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": Omega_cdm not found in settings file, setting to default (1)." << endl;
		cosmo.Omega_cdm = 1.;
	}

	cosmo.Omega_m = cosmo.Omega_cdm + cosmo.Omega_b;
	for (i = 0; i < cosmo.num_ncdm; i++)
		cosmo.Omega_m += cosmo.Omega_ncdm[i];

	//-----------Modification---------------
	//
	double H0fac = 1 / (cSI / 1e5);
	double fourpiG2 = 1.5 * H0fac * H0fac;
	COUT << "Parser defined fourpiG2 " << fourpiG2 << endl;

	if (!parseParameter(params, numparam, "Omega_asymmetron", cosmo.Omega_asymmetron))
	{
		cosmo.Omega_asymmetron = 0.0; // just add to cosmological constant
	}
	if (!parseParameter(params, numparam, "mu_as", cosmo.mu_as))
	{ // if not defined
		if (parseParameter(params, numparam, "Lc_as", cosmo.Lc_as))
		{
			// cosmo.mu_as  = 1. / (sqrt(2)*cosmo.Lc_as);
			COUT << "Defined mu from Compton physical Compton in Mpc/h" << endl;
		}
		else if (parseParameter(params, numparam, "xistar", cosmo.xistar))
		{
			COUT << "Defined mu from xistar ";
		}
		else
		{
			cosmo.xistar = 3.33e-4;
			COUT << "Defined mu from default (Lc=1 Mpc/h)" << endl;
		}
		cosmo.mu_as = H0fac / (sqrt(2.) * cosmo.xistar);
	}
	else
	{
		COUT << "Defined mu directly\n";
		cosmo.xistar = H0fac / (sqrt(2.) * cosmo.mu_as);
	}
	COUT << "Compton is " << 1. / cosmo.mu_as / sqrt(2)
		 << "Mpc/h" << endl;
	double deltastar;
	if (!parseParameter(params, numparam, "M_as", cosmo.M_as))
	{
		if (parseParameter(params, numparam, "zssb_as", cosmo.zssb_as))
		{
			COUT << "Defined M from zssb FIX" << endl;
		}
		else if (parseParameter(params, numparam, "astar", cosmo.astar))
		{
			cosmo.M_as = cosmo.xistar * sqrt(6. * cosmo.Omega_m / pow(cosmo.astar, 3.0));
			COUT << "Defined M from astar" << endl;
		}
		else if (parseParameter(params, numparam, "deltastar",deltastar))
		{
			double rhocrit = 3./2998./2998.;
			double rhostar = (deltastar+1)*rhocrit;
			cosmo.M_as = sqrt(rhostar)/cosmo.mu_as;
			cosmo.astar = pow(pow(cosmo.xistar/cosmo.M_as,2)*6.*cosmo.Omega_m,0.33);
			COUT << "Defined M from deltastar. astar is " << cosmo.astar << endl;
		}
		else
		{
			// cosmo.M_as = sqrt( cosmo.Omega_m * pow(1.+3.,3.0) ) / cosmo.mu_as;
			COUT << "Defined M from default (a=0.5)" << endl;
			cosmo.astar = 0.5;
		}
		cosmo.M_as = cosmo.xistar * sqrt(6 * cosmo.Omega_m / pow(cosmo.astar, 3.));
	}
	else
	{
		COUT << "Defined M directly\n";
		cosmo.astar = pow(pow(cosmo.xistar / cosmo.M_as, 2.) * 6 * cosmo.Omega_m, 1. / 3.);
	}
	//------------asymmetric part-------- check kappa/betas
	if (!parseParameter(params, numparam, "kappa_as", cosmo.kappa_as))
	{
		cosmo.kappa_as = 0.;
	}
	else
		COUT << "Defined kappa directly\n";
	if (!parseParameter(params, numparam, "kappa_sign", cosmo.kappa_sign))
	{
		cosmo.kappa_sign = -1; // this makes the positive minimum deepest / largest value
	}
	else
	{
		cosmo.kappa_sign = -cosmo.kappa_sign / fabs(cosmo.kappa_sign);
	}
	if (!parseParameter(params, numparam, "betastar", cosmo.betastar))
	{
		cosmo.betastar = 1.;
		COUT << "Defined betastar from default (beta=1)"
			 << "only relevant if not setting manually" << endl;
	}
	if (!parseParameter(params, numparam, "betastar2", cosmo.betastar2))
	{
		cosmo.betastar2 = cosmo.betastar;
		COUT << "Defined betastar2 from default (beta+=beta-)"
			 << "\tonly relevant if not setting manually" << endl;
	}

	// for now, use cosmo.betastar and cosmo.betastar2
	if (cosmo.kappa_as == 0)
	{
		double betav = (cosmo.betastar + cosmo.betastar2) / 2.;
		double dbet = fabs(cosmo.betastar - cosmo.betastar2);
		cosmo.lambda_as = cosmo.mu_as * cosmo.mu_as / (pow(cosmo.M_as, 4) * (pow(betav, 2) - 0.25 * pow(dbet, 2)));
		cosmo.kappa_as = cosmo.lambda_as * dbet * pow(cosmo.M_as, 2) * cosmo.kappa_sign;
	}
	else
	{
		if (parseParameter(params, numparam, "lambda_as", cosmo.lambda_as))
		{
			COUT << "Defined lambda_as directly\n";
		}
		else
			COUT << "Did not define lambda! Fix\n";
	}

	cosmo.kappa2_as = cosmo.kappa_as / (cosmo.mu_as * sqrt(cosmo.lambda_as));
	cosmo.aetanorm = 2. * fourpiG2 / (pow(cosmo.M_as * cosmo.mu_as, 2.0));
	cosmo.prefac = pow(cosmo.mu_as * sim.boxsize, 2.0);
	cosmo.dAcnf = 0.5 * pow(cosmo.mu_as / cosmo.M_as, 2) / cosmo.lambda_as;

	//-----------------------------------------------
	if (!parseParameter(params, numparam, "fifth", sim.fifth))
	{
		sim.fifth = 1; 
	}
	if (!parseParameter(params, numparam, "updateScalar", sim.updateScalar))
	{
		sim.updateScalar = 1;
	}
	if (!parseParameter(params, numparam, "updatePtcls", sim.updatePtcls))
	{
		sim.updatePtcls = 1;
	}
	if (!parseParameter(params, numparam, "solveEinstein", sim.solveEinstein))
	{
		sim.solveEinstein = 1;
	}
	if (!parseParameter(params, numparam, "As_solver", sim.As_solver))
	{
		// old convention parameters
		if (parseParameter(params, numparam, "leapfrog", sim.As_solver))
		{
			if (strcmp(sim.As_solver, "1"))
			{
				strcpy(sim.As_solver, "leapfrog");
			}
			else if (parseParameter(params, numparam, "euler", sim.As_solver))
			{
				if (strcmp(sim.As_solver, "1"))
				{
					strcpy(sim.As_solver, "euler");
				}
				else if (parseParameter(params, numparam, "QSA", sim.As_solver))
				{
					if (strcmp(sim.As_solver, "1"))
					{
						strcpy(sim.As_solver, "QSA");
					}
					else
						strcpy(sim.As_solver, "leapfrog");
				}
			}
		}
	}
	if (!parseParameter(params, numparam, "bckopt", sim.bckopt))
	{
		strcpy(sim.bckopt, "rho");
	}
	if (!parseParameter(params, numparam, "As_compute_dynbck", sim.As_compute_dynbck))
	{
		sim.As_compute_dynbck = 0;
	}
	if (!parseParameter(params, numparam, "As_use_dynbck", sim.As_use_dynbck))
	{
		if (strcmp(sim.bckopt, "symdyn") == 0)
		{
			sim.As_use_dynbck = 1;
			sim.As_compute_dynbck = 1;
		}
		else
			sim.As_use_dynbck = 0;
		COUT << "Using dynamic background " << sim.As_use_dynbck << endl;
	}
	if (!parseParameter(params, numparam, "nAs_numsteps_dynbck", sim.nAs_numsteps_dynbck))
	{
		sim.nAs_numsteps_dynbck = 1;
	}
	if (!parseParameter(params, numparam, "achifile", sim.achifile))
	{
		strcpy(sim.achifile, "no");
	}
	if (!parseParameter(params, numparam, "aqfile", sim.aqfile))
	{
		strcpy(sim.aqfile, "no");
	}
	if (!parseParameter(params, numparam, "As_ICopt", sim.As_ICopt))
	{
		strcpy(sim.As_ICopt, "no");
	}
	if (!parseParameter(params, numparam, "nAs_numsteps", sim.nAs_numsteps))
	{
		double temp;
		if (parseParameter(params, numparam, "As_Courant", sim.As_Courant))
		{
			sim.nAs_numsteps = (int)(sim.Cf / sim.As_Courant + 0.5);
		}
		else if (parseParameter(params, numparam, "As_dtau", temp))
		{
			sim.As_Courant = (temp * sim.numpts / sim.boxsize);
			sim.nAs_numsteps = (int)(sim.Cf / sim.As_Courant + 0.5);
		}
		else
		{
			sim.As_Courant = sim.Cf;
			sim.nAs_numsteps = 1; // relative to external loop number
		}
	}
	else
	{
		sim.As_Courant = sim.Cf / sim.nAs_numsteps;
	}
	if (!parseParameter(params, numparam, "Cf_MG", sim.Cf_MG))
	{
		sim.Cf_MG = sim.Cf;
	}
	COUT << "Courant factor is " << sim.Cf << endl
		 << "During MG Courant is " << sim.Cf_MG << endl
		 << "While for scalar field " << sim.As_Courant << endl;
	if (!parseParameter(params, numparam, "As_source_gravity", sim.As_source_gravity))
	{
		sim.As_source_gravity = 1; // 0 is false
	}
	if (!parseParameter(params, numparam, "As_compute_SE", sim.As_compute_SE))
	{						   // whether to compute SE tensor even if not needed
		sim.As_compute_SE = 1; // 0 is false
	}
	if (!parseParameter(params, numparam, "As_use_phiJ_SE", sim.As_use_phiJ_SE))
	{
		sim.As_use_phiJ_SE = 1; // 0 is false
	}
	if (!parseParameter(params, numparam, "animation", sim.animation))
	{
		sim.animation = 0;
	}
	if (!parseParameter(params, numparam, "zanimationstart", sim.zanimationstart))
	{
		sim.zanimationstart = 1.5;
	}
	if (!parseParameter(params, numparam, "zanimationend", sim.zanimationend))
	{
		sim.zanimationend = 1.5;
	}
	if (!parseParameter(params, numparam, "animationeveryn", sim.animationeveryn))
	{
		sim.animationeveryn = 1;
	}
	if (!parseParameter(params, numparam, "animation_downgrade_factor", sim.animation_downgrade_factor))
	{
		sim.animation_downgrade_factor = 1;
	}
	if (!parseParameter(params, numparam, "downgrade factor 2", sim.downgrade_factor2))
	{
		sim.downgrade_factor2 = 1;
	}
	if (!parseParameter(params, numparam, "animation_thickness", sim.animation_thickness))
	{
		sim.animation_thickness = 1;
	}
	if (!parseParameter(params, numparam, "animation_xcoord", sim.animation_xcoord))
	{
		sim.animation_xcoord = 1;
	}
	if (!parseParameter(params, numparam, "As_5th_force", sim.As_5th_force))
	{
		sim.As_5th_force = 1;
	}
	if (!parseParameter(params, numparam, "As_5th_force_Newton", sim.As_5th_force_Newton))
	{
		sim.As_5th_force_Newton = 0;
	}
	if (!parseParameter(params, numparam, "QSAq", sim.QSAq))
	{ // whether to have q as the velocity between loops a^2 Delta achi / dtau
		sim.QSAq = 0;
	}
	if (!parseParameter(params, numparam, "QSA_guess_dyn", sim.QSA_guess_dyn))
	{ // whether to have dynamical background starting guess for iterations
		sim.QSA_guess_dyn = 0;
	}
	if (!parseParameter(params, numparam, "QSAmodulateAmp", sim.QSAmodulateAmp))
	{ // whether to have dynamical background starting guess for iterations
		sim.QSAmodulateAmp = 1;
	}
	else if (sim.As_compute_dynbck == 0)
	{
		sim.As_compute_dynbck = 1;
	}
	if (!parseParameter(params, numparam, "GStol", sim.GStol))
	{
		sim.GStol = 1;
	}
	if (!parseParameter(params, numparam, "GStol2", sim.GStol2))
	{
		sim.GStol2 = 1;
	}
	if (!parseParameter(params, numparam, "aMG", sim.aMG))
	{
		if (strcmp(sim.As_solver, "QSA") == 0)
			sim.aMG = cosmo.astar;
		else
			sim.aMG = 1. / (sim.z_in + 1);
	}
	if (!parseParameter(params, numparam, "ICamp", sim.ICamp))
	{
		sim.ICamp = 1e-20;
	}
	if (!parseParameter(params, numparam, "ICamp2", sim.ICamp2))
	{
		sim.ICamp2 = 0;
	}
	if (!parseParameter(params, numparam, "As_full_trace", sim.As_full_trace))
	{
		sim.As_full_trace = 1; // have full trace in eom
	}
	if (!parseParameter(params, numparam, "breakafterzini", sim.breakafterzini))
	{
		sim.breakafterzini = 0;
	}
	if (!parseParameter(params, numparam, "breakWhenNoDW", sim.breakWhenNoDW))
	{
		sim.breakWhenNoDW = 0;
	}
	if (!parseParameter(params, numparam, "GWhighkAdiabatic", sim.GWhighkAdiabatic))
	{
		sim.GWhighkAdiabatic = 0;
	}
	if (!parseParameter(params, numparam, "GWstartAdiabatic", sim.GWstartAdiabatic))
	{
		sim.GWstartAdiabatic = 0;
	}
	if (!parseParameter(params, numparam, "GWinitialiseCycle", sim.GWinitialiseCycle))
	{
		sim.GWinitialiseCycle = 3;
	}
	else
	{
		if (sim.GWinitialiseCycle < 3)
		{
			COUT << COLORTEXT_RED << "NB, initialising GWs at cycle<3 contaminates the GWs. Setting to 3." << COLORTEXT_RESET << endl;
			sim.GWinitialiseCycle = 3;
		}
	}
	if (!parseParameter(params, numparam, "GWreadFields", sim.GWreadFields))
	{
		sim.GWreadFields = 0;
	}
	if (sim.GWreadFields > 0)
	{
		if (!parseParameter(params, numparam, "hijfile", sim.hijfile))
		{
			strcpy(sim.hijfile, "no");
		}
		if (!parseParameter(params, numparam, "hijprimefile", sim.hijprimefile))
		{
			strcpy(sim.hijprimefile, "no");
		}
		if (strcmp(sim.hijfile, "no") == 0 || strcmp(sim.hijprimefile, "no") == 0)
		{
			COUT << COLORTEXT_RED << "Did not provide initialisation files for GWs" << COLORTEXT_RESET << endl;
#ifdef LATFIELD2_HPP
			parallel.abortForce();
#endif
		}
	}
	if (!parseParameter(params, numparam, "writefirstsnap", sim.writefirstsnap))
	{
		sim.writefirstsnap = 0; // output ICs before evolving
	}
	if (!parseParameter(params, numparam, "writefirstspec", sim.writefirstspec))
	{
		sim.writefirstspec = 0; // output ICs before evolving
	}
	if (!parseParameter(params, numparam, "breakafterz", sim.breakafterz))
	{
		sim.breakafterz = -1;
	}
	if (strcmp(sim.bckopt, "sym") == 0 || strcmp(sim.bckopt, "symdyn") == 0)
	{
		double phib = (cosmo.mu_as / sqrt(cosmo.lambda_as)) * sqrt(1 - pow(cosmo.astar, 3));
		cosmo.Omega_sym0 = (0.5 * cosmo.Omega_m * pow(phib / cosmo.M_as, 2) + (-0.5 * pow(cosmo.mu_as * phib, 2) + 0.25 * cosmo.lambda_as * pow(phib, 4)) * 2998. * 2998 / 3.) * (double)sim.As_source_gravity;
		cosmo.Omega_asymmetron += cosmo.Omega_sym0;
	}
	else
		cosmo.Omega_sym0 = 0.;
	if (sim.fifth == 0)
	{
		// still allow to set to adjust Omega_Lambda if desired
		// cosmo.Omega_asymmetron = 0;
		cosmo.Omega_sym0 = 0;
	}
	cosmo.Omega_Lambda = 1. - cosmo.Omega_m - cosmo.Omega_rad - cosmo.Omega_fld - cosmo.Omega_asymmetron;

	// lightcone added
	COUT << "Number of lightcones: " << sim.num_lightcone << endl;
	for (i = 0; i < sim.num_lightcone; i++)
	{
		snprintf(par_string, SIZPARSTRING, "lightcone %d conformal time", i);
		if (!parseParameter(params, numparam, par_string, sim.lightcone[i].tauobs))
		{
			double avgsource = 0.; // not used in _old
			double fourpiG = fourpiG2 * sim.boxsize * sim.boxsize;

			sim.lightcone[i].tauobs = particleHorizon_old2(1. / (1. + sim.lightcone[i].z), fourpiG, cosmo);
			COUT << "Using LCDM tauobs " << sim.lightcone[i].tauobs << " for lightcone " << i << endl;
			COUT << "Is corresponding to lightcone observer at redshift " << sim.lightcone[i].z << endl;
		}
		else
			COUT << "Using parameter tauobs " << sim.lightcone[i].tauobs << " for lightcone " << i << endl;
	}

	//====================================================================

	if (cosmo.Omega_m <= 0. || cosmo.Omega_m > 1.)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": total matter density out of range!" << endl;
#ifdef LATFIELD2_HPP
		parallel.abortForce();
#endif
	}
	else if (cosmo.Omega_rad < 0. || cosmo.Omega_rad > 1. - cosmo.Omega_m)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": total radiation energy density out of range!" << endl;
#ifdef LATFIELD2_HPP
		parallel.abortForce();
#endif
	}
	else
	{
		//-------------Modification-----------------------
		COUT << "Asymmetron parameters: \n"
			 << "Asymmetron source gravity = " << sim.As_source_gravity << ", 5th force turned on = " << sim.As_5th_force << ", Number of asymmetron updates = " << sim.nAs_numsteps
			 << ", phenomenological parameters (xistar,astar,betastar,betastar2) = ("
			 << cosmo.xistar << ", " << cosmo.astar << ", " << cosmo.betastar << ", " << cosmo.betastar2 << ")" << endl;
		COUT << " cosmological parameters are: Omega_m0 = " << cosmo.Omega_m << ", Omega_rad0 = " << cosmo.Omega_rad
			 << ", Omega_lambda = " << cosmo.Omega_Lambda << ", h = " << cosmo.h << ", Omega_asymmetron0 = " << cosmo.Omega_asymmetron << ", M_as = " << cosmo.M_as << ", mu_as = " << cosmo.mu_as
			 << ", lambda_as = " << cosmo.lambda_as << ", kappa_as = " << cosmo.kappa_as << endl;
		//-------------------------------------------------
	}

	if (!parseParameter(params, numparam, "switch delta_rad", sim.z_switch_deltarad))
		sim.z_switch_deltarad = 0.;

	i = MAX_PCL_SPECIES - 2;
	if (!parseParameter(params, numparam, "switch delta_ncdm", sim.z_switch_deltancdm, i))
	{
		for (i = 0; i < MAX_PCL_SPECIES - 2; i++)
			sim.z_switch_deltancdm[i] = (ic.numtile[(sim.baryon_flag == 1) ? 2 + i : 1 + i] > 0) ? sim.z_in : 0.;
	}
	else
	{
		for (; i < MAX_PCL_SPECIES - 2; i++)
			sim.z_switch_deltancdm[i] = sim.z_switch_deltancdm[i - 1];
	}

	if (!parseParameter(params, numparam, "switch linear chi", sim.z_switch_linearchi))
	{
		if (sim.gr_flag > 0)
			sim.z_switch_linearchi = 0.;
		else
			sim.z_switch_linearchi = 0.011;

		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.z_switch_linearchi < sim.z_switch_deltancdm[i])
				sim.z_switch_linearchi = sim.z_switch_deltancdm[i];
		}
	}
	else if (sim.gr_flag == 0 && sim.z_switch_linearchi <= 0.01)
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": with garavity theory = Newton the switch linear chi redshift must be larger than 0.01." << endl;
		COUT << "              setting switch linear chi = 0.011" << endl;
		sim.z_switch_linearchi = 0.011;
	}

	i = MAX_PCL_SPECIES - 2;
	if (!parseParameter(params, numparam, "switch B ncdm", sim.z_switch_Bncdm, i))
	{
		for (i = 0; i < MAX_PCL_SPECIES - 2; i++)
			sim.z_switch_Bncdm[i] = sim.z_switch_deltancdm[i];
	}
	for (; i < MAX_PCL_SPECIES - 2; i++)
		sim.z_switch_Bncdm[i] = sim.z_switch_Bncdm[i - 1];

	for (i = 0; i < numparam; i++)
	{
		if (params[i].used)
			usedparams++;
	}

	return usedparams;
}

#endif
