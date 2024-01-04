//////////////////////////
// output.hpp
//////////////////////////
//
// Output of snapshots, light cones and spectra
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: March 2020
//
//////////////////////////

#ifndef OUTPUT_HEADER
#define OUTPUT_HEADER

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

using namespace std;

//////////////////////////
// writeSnapshots
//////////////////////////
// Description:
//   output of snapshots
//
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   a              scale factor
//   snapcount      snapshot index
//   h5filename     base name for HDF5 output file
//   pcls_cdm       pointer to particle handler for CDM
//   pcls_b         pointer to particle handler for baryons
//   pcls_ncdm      array of particle handlers for
//                  non-cold DM (may be set to NULL)
//   phi            pointer to allocated field
//   achi           pointer to allocated field
//   aq            pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   source         pointer to allocated field
//   hij            pointer to allocated field
//   T00_As		  	pointer to allocated field
//   T0i_As       	pointer to allocated field
//   Tij_As       	pointer to allocated field
//   scalarFT       pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_phi       pointer to FFT planner
//   plan_chi       pointer to FFT planner
//   plan_Bi        pointer to FFT planner
//   plan_source    pointer to FFT planner
//   plan_hij       pointer to FFT planner
//   Bi_check       pointer to allocated field
//   BiFT_check     pointer to allocated field
//   plan_Bi_check  pointer to FFT planner
//   vi             pointer to allocated field
//
// Returns:
//
//////////////////////////
// ----------Modification-------------- Parameters below and description above
void writeSnapshots(int isSecondary, metadata &sim, cosmology &cosmo, double avgsource, const double fourpiG,
					const double a, const double dtau_old, const int snapcount, string h5filename,
					Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_cdm,
					Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_b,
					Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_ncdm,
					Field<Real> *phi, Field<Real> *achi, Field<Real> *aq, Field<Real> *chi, Field<Real> *Bi,
					Field<Real> *hij, Field<Real> *source, Field<Real> *scalar,
					PlanFFT<Cplx> *plan_source, PlanFFT<Cplx> *plan_scalar,
					PlanFFT<Cplx> *plan_Sij, Field<Cplx> *SijFT, int &done_hij
#ifdef GW
					,
					Field<Cplx> *hijFT, Field<Cplx> *hijprimeFT, PlanFFT<Cplx> *plan_hij, PlanFFT<Cplx> *plan_hijprime
#endif
)
{
	char filename[2 * PARAM_MAX_LENGTH + 24];
	char buffer[64];
	int i;
	gadget2_header hdr;
	Site x(phi->lattice());
	Real divB, curlB, divh, traceh, normh;
	double dtau_pos = 0.;
	string secondMod = "Second";
	int downgrade_factor = sim.downgrade_factor;
	int out_snapshot = sim.out_snapshot;
	double z_snapshot[MAX_OUTPUTS];

	if (isSecondary > 0)
	{
		h5filename += secondMod;
		downgrade_factor = sim.downgrade_factor2;
		out_snapshot = sim.out_snapshot2;
		for (i = 0; i < sim.num_snapshot2; i++)
		{
			z_snapshot[i] = sim.z_snapshot2[i];
		}
	}
	else
	{
		for (i = 0; i < sim.num_snapshot; i++)
			z_snapshot[i] = sim.z_snapshot[i];
	}
	if (snapcount < 1e3)
		snprintf(filename, SIZFILENAME, "%03d", snapcount);
	else
		snprintf(filename, SIZFILENAME, "%04d", snapcount);
#ifdef EXTERNAL_IO
	while (ioserver.openOstream() == OSTREAM_FAIL)
		;
	if (out_snapshot & MASK_PCLS)
	{
		pcls_cdm->saveHDF5_server_open(h5filename + filename + "_cdm");
		if (sim.baryon_flag)
			pcls_b->saveHDF5_server_open(h5filename + filename + "_b");
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.numpcl[1 + sim.baryon_flag + i] == 0)
				continue;
			snprintf(buffer, SIZBUFFER, "_ncdm%d", i);
			pcls_ncdm[i].saveHDF5_server_open(h5filename + filename + buffer);
		}
	}

	if (out_snapshot & MASK_T00)
		source->saveHDF5_server_open(h5filename + filename + "_T00");

	if (out_snapshot & MASK_B)
		Bi->saveHDF5_server_open(h5filename + filename + "_B");

	if (out_snapshot & MASK_PHI)
		phi->saveHDF5_server_open(h5filename + filename + "_phi");
	//-------------Modification----------------------------
	if (out_snapshot & MASK_ACHI && sim.aMG <= a)
		achi->saveHDF5_server_open(h5filename + filename + "_achi");
	if (out_snapshot & MASK_AQ && sim.aMG <= a)
		aq->saveHDF5_server_open(h5filename + filename + "_aq");
	//------------------------------------------------------
	if (out_snapshot & MASK_CHI)
		chi->saveHDF5_server_open(h5filename + filename + "_chi");

	if (out_snapshot & MASK_HIJ)
	{
#ifdef GW
		plan_hij.execute(FFT_BACKWARD);
		hij->updateHalo();
#endif
		hij->saveHDF5_server_open(h5filename + filename + "_hij");
	}
#ifdef GW
	if (out_snapshot & MASK_HIJPRIME)
	{
		plan_hijprime->execute(FFT_BACKWARD);
		hij->updateHalo();
		hij->saveHDF5_server_open(h5filename + filename + "_hij_prime");
	}
	if (out_snapshot & MASK_HIJPRIMENORM)
	{
		if (!(out_snapshot & MASK_HIJPRIME))
			plan_hijprime->execute(FFT_BACKWARD);
		for (x.first(); x.test(); x.next())
		{
			(*scalar)(x) = +(*hij)(x, 0, 0) * (*hij)(x, 0, 0) + (*hij)(x, 1, 1) * (*hij)(x, 1, 1)		/* */
						   + (*hij)(x, 2, 2) * (*hij)(x, 2, 2) + 2. * (*hij)(x, 0, 1) * (*hij)(x, 0, 1) /* */
						   + 2. * (*hij)(x, 0, 2) * (*hij)(x, 0, 2) + 2. * (*hij)(x, 1, 2) * (*hij)(x, 1, 2);
		}
		scalar->updateHalo();
		scalar->saveHDF5_server_open(h5filename + filename + "_hij_prime_norm");
	}
#endif
#endif

	if (out_snapshot & MASK_T00)
	{
#ifdef EXTERNAL_IO
		source->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (downgrade_factor > 1)
			source->saveHDF5_coarseGrain3D(h5filename + filename + "_T00.h5", downgrade_factor);
		else
			source->saveHDF5(h5filename + filename + "_T00.h5");
#endif
	}

	if (out_snapshot & MASK_B)
	{
		computeVectorDiagnostics(*Bi, divB, curlB);
		double toPhysicalFactor = 1. / (a * a * sim.numpts);
		COUT << " B diagnostics: max |divB| = " << divB * toPhysicalFactor << ", max |curlB| = " << curlB * toPhysicalFactor << endl;

#ifdef EXTERNAL_IO
		Bi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (downgrade_factor > 1)
			Bi->saveHDF5_coarseGrain3D(h5filename + filename + "_B.h5", downgrade_factor);
		else
			Bi->saveHDF5(h5filename + filename + "_B.h5");
#endif
	}
	if (out_snapshot & MASK_PHI)
	{
#ifdef EXTERNAL_IO
		phi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (downgrade_factor > 1)
			phi->saveHDF5_coarseGrain3D(h5filename + filename + "_phi.h5", downgrade_factor);
		else
			phi->saveHDF5(h5filename + filename + "_phi.h5");
#endif
		//-------------------Modification----------------------------
		if (out_snapshot & MASK_ACHI && sim.aMG <= a)
		{
#ifdef EXTERNAL_IO
			achi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
			if (downgrade_factor > 1)
				achi->saveHDF5_coarseGrain3D(h5filename + filename + "_achi.h5", downgrade_factor);
			else
				achi->saveHDF5(h5filename + filename + "_achi.h5");
#endif
		}

		if (out_snapshot & MASK_AQ && sim.aMG <= a)
		{
#ifdef EXTERNAL_IO
			aq->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
			if (downgrade_factor > 1)
				aq->saveHDF5_coarseGrain3D(h5filename + filename + "_aq.h5", downgrade_factor);
			else
				aq->saveHDF5(h5filename + filename + "_aq.h5");
#endif
		}

		//------------------------------------------------------------
		if (out_snapshot & MASK_CHI)
		{
#ifdef EXTERNAL_IO
			chi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
			if (downgrade_factor > 1)
				chi->saveHDF5_coarseGrain3D(h5filename + filename + "_chi.h5", downgrade_factor);
			else
				chi->saveHDF5(h5filename + filename + "_chi.h5");
#endif
		}
#ifdef GW
		if (out_snapshot & MASK_HIJ)
		{
			plan_hij->execute(FFT_BACKWARD);
			hij->updateHalo();

			computeTensorDiagnostics(*hij, divh, traceh, normh);
			COUT << " GW diagnostics: max |divh| = " << divh << ", max |traceh| = " << traceh << ", max |h| = " << normh << endl;

#ifdef EXTERNAL_IO
			hij->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
			if (downgrade_factor > 1)
				hij->saveHDF5_coarseGrain3D(h5filename + filename + "_hij.h5", downgrade_factor);
			else
				hij->saveHDF5(h5filename + filename + "_hij.h5");
#endif
		}
		if (out_snapshot & MASK_HIJPRIME)
		{
			plan_hijprime->execute(FFT_BACKWARD);
			hij->updateHalo();
#ifdef EXTERNAL_IO
			hij->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
			if (downgrade_factor > 1)
				hij->saveHDF5_coarseGrain3D(h5filename + filename + "_hij_prime.h5", downgrade_factor);
			else
				hij->saveHDF5(h5filename + filename + "_hij_prime.h5");
#endif
		}
		if (out_snapshot & MASK_HIJPRIMENORM)
		{
			if (!(out_snapshot & MASK_HIJPRIME))
				plan_hijprime->execute(FFT_BACKWARD);
			for (x.first(); x.test(); x.next())
			{
				(*scalar)(x) = +(*hij)(x, 0, 0) * (*hij)(x, 0, 0) + (*hij)(x, 1, 1) * (*hij)(x, 1, 1)		/* */
							   + (*hij)(x, 2, 2) * (*hij)(x, 2, 2) + 2. * (*hij)(x, 0, 1) * (*hij)(x, 0, 1) /* */
							   + 2. * (*hij)(x, 0, 2) * (*hij)(x, 0, 2) + 2. * (*hij)(x, 1, 2) * (*hij)(x, 1, 2);
			}
			scalar->updateHalo();
#ifdef EXTERNAL_IO
			scalar->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
			if (downgrade_factor > 1)
				scalar->saveHDF5_coarseGrain3D(h5filename + filename + "_hij_prime_norm.h5", downgrade_factor);
			else
				scalar->saveHDF5(h5filename + filename + "_hij_prime_norm.h5");
#endif
		}
#else
		if (out_snapshot & MASK_HIJ)
		{
			if (done_hij == 0)
			{
				projectFTtensor(*SijFT, *SijFT); // adiabatic
				plan_Sij->execute(FFT_BACKWARD);
				hij->updateHalo();
				done_hij = 1;
			}

			computeTensorDiagnostics(*hij, divh, traceh, normh);
			COUT << " GW diagnostics: max |divh| = " << divh << ", max |traceh| = " << traceh << ", max |h| = " << normh << endl;

#ifdef EXTERNAL_IO
			hij->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
			if (downgrade_factor > 1)
				hij->saveHDF5_coarseGrain3D(h5filename + filename + "_hij.h5", downgrade_factor);
			else
				hij->saveHDF5(h5filename + filename + "_hij.h5");
#endif
		}

#endif

		if (out_snapshot & MASK_GADGET)
		{
			if (out_snapshot & MASK_MULTI)
				hdr.num_files = parallel.grid_size()[1];
			else
				hdr.num_files = 1;
			hdr.Omega0 = cosmo.Omega_m;
			hdr.OmegaLambda = cosmo.Omega_Lambda;
			hdr.HubbleParam = cosmo.h;
			hdr.BoxSize = sim.boxsize / GADGET_LENGTH_CONVERSION;
			hdr.flag_sfr = 0;
			hdr.flag_cooling = 0;
			hdr.flag_feedback = 0;
			hdr.flag_age = 0;
			hdr.flag_metals = 0;
			for (i = 0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4; i++)
				hdr.fill[i] = 0;
			for (i = 0; i < 6; i++)
			{
				hdr.npart[i] = 0;
				hdr.npartTotal[i] = 0;
				hdr.npartTotalHW[i] = 0;
				hdr.mass[i] = 0.;
			}

#ifdef EXACT_OUTPUT_REDSHIFTS
			hdr.time = 1. / (z_snapshot[snapcount] + 1.);
			hdr.redshift = z_snapshot[snapcount];
			dtau_pos = (hdr.time - a) / a / Hconf(a, fourpiG, cosmo, avgsource);
#else
			hdr.time = a;
			hdr.redshift = (1. / a) - 1.;
#endif
			if (sim.tracer_factor[0] > 0)
			{
				hdr.npart[1] = (uint32_t)(((sim.numpcl[0] % sim.tracer_factor[0]) ? (1 + (sim.numpcl[0] / sim.tracer_factor[0])) : (sim.numpcl[0] / sim.tracer_factor[0])) % (1ll << 32));
				hdr.npartTotal[1] = hdr.npart[1];
				hdr.npartTotalHW[1] = (uint32_t)(((sim.numpcl[0] % sim.tracer_factor[0]) ? (1 + (sim.numpcl[0] / sim.tracer_factor[0])) : (sim.numpcl[0] / sim.tracer_factor[0])) / (1ll << 32));
				if (sim.baryon_flag)
				{
					hdr.mass[1] = (double)sim.tracer_factor[0] * C_RHO_CRIT * cosmo.Omega_cdm * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;
				}
				else
					hdr.mass[1] = (double)sim.tracer_factor[0] * C_RHO_CRIT * (cosmo.Omega_cdm + cosmo.Omega_b) * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;

				if (hdr.npartTotalHW[1] != 0 && hdr.num_files == 1)
				{
					COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of particles (" << (sim.numpcl[0] / sim.tracer_factor[0]) << ") in Gadget2 file exceeds limit (4294967295). Try using multi-Gadget2 output format." << endl;
				}
				else
				{
					pcls_cdm->saveGadget2(h5filename + filename + "_cdm", hdr, sim.tracer_factor[0], dtau_pos, dtau_pos + 0.5 * dtau_old, phi);
				}
			}
			if (sim.baryon_flag && sim.tracer_factor[1] > 0)
			{
				hdr.npart[1] = (uint32_t)(((sim.numpcl[1] % sim.tracer_factor[1]) ? (1 + (sim.numpcl[1] / sim.tracer_factor[1])) : (sim.numpcl[1] / sim.tracer_factor[1])) % (1ll << 32));
				hdr.npartTotal[1] = hdr.npart[1];
				hdr.npartTotalHW[1] = (uint32_t)(((sim.numpcl[1] % sim.tracer_factor[1]) ? (1 + (sim.numpcl[1] / sim.tracer_factor[1])) : (sim.numpcl[1] / sim.tracer_factor[1])) / (1ll << 32));
				hdr.mass[1] = (double)sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
				if (hdr.npartTotalHW[1] != 0 && hdr.num_files == 1)
				{
					COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of particles (" << (sim.numpcl[1] / sim.tracer_factor[1]) << ") in Gadget2 file exceeds limit (4294967295). Try using multi-Gadget2 output format." << endl;
				}
				else
					pcls_b->saveGadget2(h5filename + filename + "_b", hdr, sim.tracer_factor[1], dtau_pos, dtau_pos + 0.5 * dtau_old, phi);
			}
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (sim.numpcl[1 + sim.baryon_flag + i] == 0 || sim.tracer_factor[i + 1 + sim.baryon_flag] == 0)
					continue;
				snprintf(buffer, SIZBUFFER, "_ncdm%d", i);
				hdr.npart[1] = (uint32_t)(((sim.numpcl[i + 1 + sim.baryon_flag] % sim.tracer_factor[i + 1 + sim.baryon_flag]) ? (1 + (sim.numpcl[i + 1 + sim.baryon_flag] / sim.tracer_factor[i + 1 + sim.baryon_flag])) : (sim.numpcl[i + 1 + sim.baryon_flag] / sim.tracer_factor[i + 1 + sim.baryon_flag])) % (1ll << 32));
				hdr.npartTotal[1] = hdr.npart[1];
				hdr.npartTotalHW[1] = (uint32_t)(((sim.numpcl[i + 1 + sim.baryon_flag] % sim.tracer_factor[i + 1 + sim.baryon_flag]) ? (1 + (sim.numpcl[i + 1 + sim.baryon_flag] / sim.tracer_factor[i + 1 + sim.baryon_flag])) : (sim.numpcl[i + 1 + sim.baryon_flag] / sim.tracer_factor[i + 1 + sim.baryon_flag])) / (1ll << 32));
				hdr.mass[1] = (double)sim.tracer_factor[i + 1 + sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[i] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[i + 1 + sim.baryon_flag] / GADGET_MASS_CONVERSION;
				if (hdr.npartTotalHW[1] != 0 && hdr.num_files == 1)
				{
					COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of particles (" << (sim.numpcl[i + 1 + sim.baryon_flag] / sim.tracer_factor[i + 1 + sim.baryon_flag]) << ") in Gadget2 file exceeds limit (4294967295). Try using multi-Gadget2 output format." << endl;
				}
				else
					pcls_ncdm[i].saveGadget2(h5filename + filename + buffer, hdr, sim.tracer_factor[i + 1 + sim.baryon_flag], dtau_pos, dtau_pos + 0.5 * dtau_old, phi);
			}
		}

		if (out_snapshot & MASK_PCLS)
		{
#ifdef EXTERNAL_IO
			pcls_cdm->saveHDF5_server_write();
			if (sim.baryon_flag)
				pcls_b->saveHDF5_server_write();
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (sim.numpcl[1 + sim.baryon_flag + i] == 0)
					continue;
				pcls_ncdm[i].saveHDF5_server_write();
			}
#else
			pcls_cdm->saveHDF5(h5filename + filename + "_cdm", 1);
			if (sim.baryon_flag)
				pcls_b->saveHDF5(h5filename + filename + "_b", 1);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (sim.numpcl[1 + sim.baryon_flag + i] == 0)
					continue;
				snprintf(buffer, SIZBUFFER, "_ncdm%d", i);
				pcls_ncdm[i].saveHDF5(h5filename + filename + buffer, 1);
			}
#endif
		}
#ifdef EXTERNAL_IO
		ioserver.closeOstream();
#endif
	}
}
void writeAnimation(metadata &sim, cosmology &cosmo, double avgsource, const double fourpiG,
					const double a, const double dtau_old, const int animcount, string h5filename,
					Field<Real> *phi, Field<Real> *achi, Field<Real> *aq, Field<Real> *chi, Field<Real> *Bi,
					Field<Real> *hij, Field<Real> *source, Field<Real> *scalar,
					PlanFFT<Cplx> *plan_source, PlanFFT<Cplx> *plan_scalar,
					PlanFFT<Cplx> *plan_Sij, Field<Cplx> *SijFT, int &done_hij
#ifdef GW
					,
					Field<Cplx> *hijFT, Field<Cplx> *hijprimeFT, PlanFFT<Cplx> *plan_hij, PlanFFT<Cplx> *plan_hijprime
#endif
)
{
	char filename[2 * PARAM_MAX_LENGTH + 24];
	char buffer[64];
	int i;
	Site x(phi->lattice());
	double dtau_pos = 0.;
	if (animcount < 1e3)
		snprintf(filename, SIZFILENAME, "%03d", animcount);
	else
		snprintf(filename, SIZFILENAME, "%04d", animcount);

	if (sim.out_animation & MASK_T00)
	{
		if (sim.animation_downgrade_factor > 1)
			source->saveSliceHDF5_coarseGrain2D(h5filename + filename + "_T00.h5", sim.animation_xcoord,
												sim.animation_thickness, sim.animation_downgrade_factor);
		else
			source->saveSliceHDF5(h5filename + filename + "_T00.h5", sim.animation_xcoord,
								  sim.animation_thickness);
	}

	if (sim.out_animation & MASK_B)
	{

		if (sim.animation_downgrade_factor > 1)
			Bi->saveSliceHDF5_coarseGrain2D(h5filename + filename + "_B.h5", sim.animation_xcoord,
											sim.animation_thickness, sim.animation_downgrade_factor);
		else
			Bi->saveSliceHDF5(h5filename + filename + "_B.h5", sim.animation_xcoord,
							  sim.animation_thickness);
	}
	if (sim.out_animation & MASK_PHI)
	{

		if (sim.animation_downgrade_factor > 1)
			phi->saveSliceHDF5_coarseGrain2D(h5filename + filename + "_phi.h5", sim.animation_xcoord,
											 sim.animation_thickness, sim.animation_downgrade_factor);
		else
			phi->saveSliceHDF5(h5filename + filename + "_phi.h5", sim.animation_xcoord,
							   sim.animation_thickness);
	}
	//-------------------Modification----------------------------
	if (sim.out_animation & MASK_ACHI && sim.aMG <= a)
	{
		if (sim.animation_downgrade_factor > 1)
			achi->saveSliceHDF5_coarseGrain2D(h5filename + filename + "_achi.h5", sim.animation_xcoord,
											  sim.animation_thickness, sim.animation_downgrade_factor);
		else
			achi->saveSliceHDF5(h5filename + filename + "_achi.h5", sim.animation_xcoord,
								sim.animation_thickness);
	}

	if (sim.out_animation & MASK_AQ && sim.aMG <= a)
	{
		if (sim.animation_downgrade_factor > 1)
			aq->saveSliceHDF5_coarseGrain2D(h5filename + filename + "_aq.h5", sim.animation_xcoord,
											sim.animation_thickness, sim.animation_downgrade_factor);
		else
			aq->saveSliceHDF5(h5filename + filename + "_aq.h5", sim.animation_xcoord,
							  sim.animation_thickness);
	}

	//------------------------------------------------------------
	if (sim.out_animation & MASK_CHI)
	{
		if (sim.animation_downgrade_factor > 1)
			chi->saveSliceHDF5_coarseGrain2D(h5filename + filename + "_chi.h5", sim.animation_xcoord,
											 sim.animation_thickness, sim.animation_downgrade_factor);
		else
			chi->saveSliceHDF5(h5filename + filename + "_chi.h5", sim.animation_xcoord,
							   sim.animation_thickness);
	}

#ifdef GW
	if (sim.out_animation & MASK_HIJ)
	{
		plan_hij->execute(FFT_BACKWARD);
		hij->updateHalo();
		if (sim.animation_downgrade_factor > 1)
			hij->saveSliceHDF5_coarseGrain2D(h5filename + filename + "_hij.h5", sim.animation_xcoord,
											 sim.animation_thickness, sim.animation_downgrade_factor);
		else
			hij->saveSliceHDF5(h5filename + filename + "_hij.h5", sim.animation_xcoord,
							   sim.animation_thickness);
	}
	if (sim.out_animation & MASK_HIJPRIME)
	{
		plan_hijprime->execute(FFT_BACKWARD);
		hij->updateHalo();
		if (sim.animation_downgrade_factor > 1)
			hij->saveSliceHDF5_coarseGrain2D(h5filename + filename + "_hij_prime.h5", sim.animation_xcoord,
											 sim.animation_thickness, sim.animation_downgrade_factor);
		else
			hij->saveSliceHDF5(h5filename + filename + "_hij_prime.h5", sim.animation_xcoord,
							   sim.animation_thickness);
	}
	if (sim.out_animation & MASK_HIJPRIMENORM)
	{
		if (!(sim.out_animation & MASK_HIJPRIME))
			plan_hijprime->execute(FFT_BACKWARD);
		for (x.first(); x.test(); x.next())
		{
			(*scalar)(x) = +(*hij)(x, 0, 0) * (*hij)(x, 0, 0) + (*hij)(x, 1, 1) * (*hij)(x, 1, 1)		/* */
						   + (*hij)(x, 2, 2) * (*hij)(x, 2, 2) + 2. * (*hij)(x, 0, 1) * (*hij)(x, 0, 1) /* */
						   + 2. * (*hij)(x, 0, 2) * (*hij)(x, 0, 2) + 2. * (*hij)(x, 1, 2) * (*hij)(x, 1, 2);
		}
		scalar->updateHalo();

		if (sim.animation_downgrade_factor > 1)
			scalar->saveSliceHDF5_coarseGrain2D(h5filename + filename + "_hij_prime_norm.h5", sim.animation_xcoord,
												sim.animation_thickness, sim.animation_downgrade_factor);
		else
			scalar->saveSliceHDF5(h5filename + filename + "_hij_prime_norm.h5", sim.animation_xcoord,
								  sim.animation_thickness);
	}
#else
	if (sim.out_animation & MASK_HIJ)
	{
		if (done_hij == 0)
		{
			projectFTtensor(*SijFT, *SijFT); // adiabatic
			plan_Sij->execute(FFT_BACKWARD);
			hij->updateHalo();
			done_hij = 1;
		}

		if (sim.animation_downgrade_factor > 1)
			hij->saveSliceHDF5_coarseGrain2D(h5filename + filename + "_hij.h5", sim.animation_xcoord,
											 sim.animation_thickness, sim.animation_downgrade_factor);
		else
			hij->saveSliceHDF5(h5filename + filename + "_hij.h5", sim.animation_xcoord,
							   sim.animation_thickness);
	}
#endif
}

double ninjhijprime(double *vec, Field<Real> *hij_prime, Site &xsim)
{
	return 0.5 * (vec[0] * vec[0] * (*hij_prime)(xsim, 0, 0) + vec[1] * vec[1] * (*hij_prime)(xsim, 1, 1) + vec[2] * vec[2] * (*hij_prime)(xsim, 2, 2) /**/
				  + 2 * (vec[0] * vec[1] * (*hij_prime)(xsim, 0, 1) + vec[0] * vec[2] * (*hij_prime)(xsim, 0, 2) + vec[1] * vec[2] * (*hij_prime)(xsim, 1, 2)));
}

//////////////////////////
// writeLightcones
//////////////////////////
// Description:
//   output of light cones
//
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   a              scale factor
//   tau            conformal time
//   dtau           conformal time step
//   dtau_old       conformal time step of previous cycle
//   maxvel         maximum cdm velocity
//   cycle          current simulation cycle
//   h5filename     base name for HDF5 output file
//   pcls_cdm       pointer to particle handler for CDM
//   pcls_b         pointer to particle handler for baryons
//   pcls_ncdm      array of particle handlers for
//                  non-cold DM (may be set to NULL)
//   phi            pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   hij            pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_Bi        pointer to FFT planner
//   plan_hij       pointer to FFT planner
//   IDbacklog      IDs of particles written in previous cycle
//
// Returns:
//
//////////////////////////

void writeLightcones(metadata &sim, cosmology &cosmo, double avgsource, const double fourpiG,
					 const double a, const double tau, const double dtau, const double dtau_old,
					 const double maxvel, const int cycle, string h5filename,
					 Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_cdm,
					 Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_b,
					 Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_ncdm,
					 Field<Real> *phi, Field<Real> *chi, Field<Real> *Bi, Field<Real> *hij_prime,
					 Field<Real> *achi, Field<Real> *aq,
					 set<long> *IDbacklog,
					 PlanFFT<Cplx> *plan_Sij, Field<Cplx> *SijFT, int &done_hij
#ifdef GW
					 ,
					 PlanFFT<Cplx> *plan_hijprime, Field<Cplx> *hijprimeFT
#endif
)
{
	// NB: edited hij part to output GW observable instead (Delta P/Delta tau /P), for pulsar period P
	// removed ability to output hij field on lightcone because not needed. Easily modified to include again.
	int i, j, n, p;
	double d;
	double vertex[MAX_INTERSECTS][3];
	double domain[6];
	double pos[3];
	double s[2];
	char filename[2 * PARAM_MAX_LENGTH + 24];
	char buffer[268];
	FILE *outfile;
	gadget2_header hdr;
	set<long> IDprelog[MAX_PCL_SPECIES];
	long *IDcombuf;
	long *IDcombuf2;
	Site xsim;

#ifdef HAVE_HEALPIX

	if (sim.out_lightcone[i] & MASK_HIJ)
	{
#ifdef GW
		plan_hijprime->execute(FFT_BACKWARD);
		hij_prime->updateHalo();
#endif
		// #ifdef GW
		// 		plan_hij->execute(FFT_BACKWARD);
		// 		hij->updateHalo();
		// #else
		// 		if (done_hij == 0)
		// 		{
		// 			projectFTtensor(*SijFT, *SijFT); // adiabatic
		// 			plan_Sij->execute(FFT_BACKWARD);
		// 			hij->updateHalo();
		// 			done_hij = 1;
		// 		}
		// #endif
	}

	xsim.initialize(phi->lattice());
	int64_t pix, pix2, q;
	vector<int> pixbatch_id;
	vector<int> sender_proc;
	vector<int> pixbatch_size[3];
	vector<int> pixbatch_delim[3];
	int pixbatch_type;
	int commdir[2];
	Real *pixbuf[LIGHTCONE_MAX_FIELDS][9];
	Real *commbuf;
	int pixbuf_size[9];
	int pixbuf_reserve[9];
	int64_t bytes, bytes2, offset2;
	vector<MPI_Offset> offset;
	char **outbuf = new char *[LIGHTCONE_MAX_FIELDS];
	healpix_header maphdr;
	double R[3][3];
	double w[3];
	double temp;
	int base_pos[3];
	int shell, shell_inner, shell_outer, shell_write;
	uint32_t blocksize;
	MPI_File mapfile;
	MPI_Status status;
	MPI_Datatype patch;
	int io_group_size;
	long pixcount = 0;

	for (j = 0; j < 9 * LIGHTCONE_MAX_FIELDS; j++)
		pixbuf[j / 9][j % 9] = NULL;

	for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
		outbuf[j] = NULL;
#endif
	domain[0] = -0.5;
	domain[1] = phi->lattice().coordSkip()[1] - 0.5;
	domain[2] = phi->lattice().coordSkip()[0] - 0.5;
	for (i = 0; i < 3; i++)
		domain[i + 3] = domain[i] + phi->lattice().sizeLocal(i) + 1.;

	for (i = 0; i < 6; i++)
		domain[i] /= (double)sim.numpts;

	for (i = 0; i < sim.num_lightcone; i++)
	{
		// Modification
		//  this is LCDM estimate if not inputted parameter in settings file
		d = sim.lightcone[i].tauobs;

		if (parallel.isRoot())
		{
			if (sim.num_lightcone > 1)
				snprintf(filename, SIZFILENAME, "%s%s%d_info.dat", sim.output_path, sim.basename_lightcone, i);
			else
				snprintf(filename, SIZFILENAME, "%s%s_info.dat", sim.output_path, sim.basename_lightcone);

			if (cycle == 0)
				outfile = fopen(filename, "w");
			else
				outfile = fopen(filename, "a");
			if (outfile == NULL)
			{
				cout << " error opening file for lightcone info!" << endl;
			}
			else if (cycle == 0)
			{
				if (sim.num_lightcone > 1) // now printing lcdm distance
					fprintf(outfile, "# information file for lightcone %d\n# geometric parameters:\n# vertex = (%f, %f, %f) Mpc/h\n# redshift = %f\n# distance = (%f - %f) Mpc/h\n# lcdm distance = %f\n# opening half-angle = %f degrees\n# direction = (%f, %f, %f)\n# cycle   tau/boxsize    a              pcl_inner        pcl_outer        metric_inner     metric_outer\n", i, sim.lightcone[i].vertex[0] * sim.boxsize, sim.lightcone[i].vertex[1] * sim.boxsize, sim.lightcone[i].vertex[2] * sim.boxsize, sim.lightcone[i].z, sim.lightcone[i].distance[0] * sim.boxsize, sim.lightcone[i].distance[1] * sim.boxsize, d, (sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180., sim.lightcone[i].direction[0], sim.lightcone[i].direction[1], sim.lightcone[i].direction[2]);
				else
					fprintf(outfile, "# information file for lightcone\n# geometric parameters:\n# vertex = (%f, %f, %f) Mpc/h\n# redshift = %f\n# distance = (%f - %f) Mpc/h\n# lcdm distance = %f\n# opening half-angle = %f degrees\n# direction = (%f, %f, %f)\n# cycle   tau/boxsize    a              pcl_inner        pcl_outer        metric_inner     metric_outer\n", sim.lightcone[i].vertex[0] * sim.boxsize, sim.lightcone[i].vertex[1] * sim.boxsize, sim.lightcone[i].vertex[2] * sim.boxsize, sim.lightcone[i].z, sim.lightcone[i].distance[0] * sim.boxsize, sim.lightcone[i].distance[1] * sim.boxsize, d, (sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180., sim.lightcone[i].direction[0], sim.lightcone[i].direction[1], sim.lightcone[i].direction[2]);
			}
		}
		// using correct tau, lcdm d amounting to error on z_observer to be quantified later
		s[0] = d - tau - 0.5 * sim.covering[i] * dtau;
		s[1] = d - tau + 0.5 * sim.covering[i] * dtau_old;
#ifdef HAVE_HEALPIX
		shell_inner = (s[0] > sim.lightcone[i].distance[1]) ? ceil(s[0] * sim.numpts * sim.shellfactor[i]) : ceil(sim.lightcone[i].distance[1] * sim.numpts * sim.shellfactor[i]);
		shell_outer = (sim.lightcone[i].distance[0] > s[1]) ? floor(s[1] * sim.numpts * sim.shellfactor[i]) : (ceil(sim.lightcone[i].distance[0] * sim.numpts * sim.shellfactor[i]) - 1);
		if (shell_outer < shell_inner && s[1] > 0) // removed && s[1]>0 criterion
		{
			shell_outer = shell_inner;
		}
		maphdr.precision = sizeof(Real);
		maphdr.Ngrid = sim.numpts;
		maphdr.direction[0] = sim.lightcone[i].direction[0];
		maphdr.direction[1] = sim.lightcone[i].direction[1];
		maphdr.direction[2] = sim.lightcone[i].direction[2];
		maphdr.boxsize = sim.boxsize;
		memset((void *)maphdr.fill, 0, 256 - 5 * 4 - 5 * 8);

		// xsim.initialize(phi->lattice()); // do earlier

		if (sim.lightcone[i].direction[0] == 0 && sim.lightcone[i].direction[1] == 0)
		{
			R[0][0] = sim.lightcone[i].direction[2];
			R[0][1] = 0;
			R[0][2] = 0;
			R[1][0] = 0;
			R[1][1] = 1;
			R[1][2] = 0;
			R[2][0] = 0;
			R[2][1] = 0;
			R[2][2] = sim.lightcone[i].direction[2];
		}
		else
		{
			temp = atan2(sim.lightcone[i].direction[1], sim.lightcone[i].direction[0]);
			R[0][0] = cos(temp) * sim.lightcone[i].direction[2];
			R[0][1] = -sin(temp);
			R[0][2] = sim.lightcone[i].direction[0];
			R[1][0] = sin(temp) * sim.lightcone[i].direction[2];
			R[1][1] = cos(temp);
			R[1][2] = sim.lightcone[i].direction[1];
			R[2][0] = -sqrt(1. - sim.lightcone[i].direction[2] * sim.lightcone[i].direction[2]);
			R[2][1] = 0;
			R[2][2] = sim.lightcone[i].direction[2];
		}
#endif
		if (sim.lightcone[i].distance[0] > s[0] && sim.lightcone[i].distance[1] <= s[1] && s[1] > 0.)
		{
			if (parallel.isRoot() && outfile != NULL)
			{
				fprintf(outfile, "%6d   %e   %e   %2.12f   %2.12f   %2.12f   %2.12f\n", cycle, tau, a, d - tau - 0.5 * dtau, d - tau + 0.5 * dtau_old, s[0], s[1]);
				fclose(outfile);

				if (sim.num_lightcone > 1)
					snprintf(filename, SIZFILENAME, "%s%s%d_info.bin", sim.output_path, sim.basename_lightcone, i);
				else
					snprintf(filename, SIZFILENAME, "%s%s_info.bin", sim.output_path, sim.basename_lightcone);
				// fix for overwrite on old info
				outfile = fopen(filename, "a");
				if (outfile == NULL)
				{
					cout << " error opening file for lightcone info!" << endl;
				}
				else
				{
					((double *)buffer)[0] = tau;
					((double *)buffer)[1] = a;
					((double *)buffer)[2] = d - tau - 0.5 * dtau;
					((double *)buffer)[3] = d - tau + 0.5 * dtau_old;
					// since d is unknown, fix post-process
					// tau and taulcdm are in background file

					fwrite((const void *)&cycle, sizeof(int), 1, outfile);
					fwrite((const void *)buffer, sizeof(double), 4, outfile);
					fwrite((const void *)s, sizeof(double), 2, outfile);

					fclose(outfile);
				}
			}

#ifdef HAVE_HEALPIX

			bytes = 0;
			bytes2 = 0;

			for (j = 0; j < 9; j++)
				pixbuf_reserve[j] = PIXBUFFER;

			if (sim.out_lightcone[i] & MASK_PHI)
			{
				for (j = 0; j < 9; j++)
					pixbuf[LIGHTCONE_PHI_OFFSET][j] = (Real *)malloc(sizeof(Real) * PIXBUFFER);
			}

			if (sim.out_lightcone[i] & MASK_CHI)
			{
				for (j = 0; j < 9; j++)
					pixbuf[LIGHTCONE_CHI_OFFSET][j] = (Real *)malloc(sizeof(Real) * PIXBUFFER);
			}

			if (sim.out_lightcone[i] & MASK_B)
			{
				for (j = 0; j < 9; j++)
				{
					pixbuf[LIGHTCONE_B_OFFSET][j] = (Real *)malloc(sizeof(Real) * PIXBUFFER);
					pixbuf[LIGHTCONE_B_OFFSET + 1][j] = (Real *)malloc(sizeof(Real) * PIXBUFFER);
					pixbuf[LIGHTCONE_B_OFFSET + 2][j] = (Real *)malloc(sizeof(Real) * PIXBUFFER);
				}
			}

#ifdef GW
			// Stochastic GWs
			if (sim.out_lightcone[i] & MASK_HIJ)
			{
				for (j = 0; j < 9; j++)
					pixbuf[LIGHTCONE_SGW_OFFSET][j] = (Real *)malloc(sizeof(Real) * PIXBUFFER);
			}
#endif

			// Modification
			if (sim.out_lightcone[i] & MASK_ACHI && sim.aMG <= a)
			{
				for (j = 0; j < 9; j++)
					pixbuf[LIGHTCONE_ACHI_OFFSET][j] = (Real *)malloc(sizeof(Real) * PIXBUFFER);
			}
			if (sim.out_lightcone[i] & MASK_AQ && sim.aMG <= a)
			{
				for (j = 0; j < 9; j++)
					pixbuf[LIGHTCONE_AQ_OFFSET][j] = (Real *)malloc(sizeof(Real) * PIXBUFFER);
			}

			//------------------

			if ((shell_outer + 1 - shell_inner) > parallel.size())
			{
				shell_write = ((shell_outer + 1 - shell_inner) * parallel.rank() + parallel.size() - 1) / parallel.size();
				io_group_size = 0;
			}
			else
			{
				shell_write = ((shell_outer + 1 - shell_inner) * parallel.rank()) / parallel.size();
				io_group_size = (((shell_write + 1) * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner)) - ((shell_write * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner));
			}

			for (shell = shell_inner; shell <= shell_outer; shell++)
			{
				maphdr.distance = (double)shell / (double)sim.numpts / sim.shellfactor[i];

				for (maphdr.Nside = sim.Nside[i][0]; maphdr.Nside < sim.Nside[i][1]; maphdr.Nside *= 2)
				{
					if (12. * maphdr.Nside * maphdr.Nside > sim.pixelfactor[i] * 4. * M_PI * maphdr.distance * maphdr.distance * sim.numpts * sim.numpts)
						break;
				}

				for (maphdr.Nside_ring = 2; 2.137937882409166 * sim.numpts * maphdr.distance / maphdr.Nside_ring > phi->lattice().sizeLocal(1) && 2.137937882409166 * sim.numpts * maphdr.distance / maphdr.Nside_ring > phi->lattice().sizeLocal(2) && maphdr.Nside_ring < maphdr.Nside; maphdr.Nside_ring *= 2)
					;

				if (sim.lightcone[i].opening > 2. / 3.)
				{
					p = 1 + (int)floor(maphdr.Nside * sqrt(3. - 3. * sim.lightcone[i].opening));
					maphdr.Npix = 2 * p * (p + 1);
				}
				else if (sim.lightcone[i].opening > -2. / 3.)
				{
					p = 1 + (int)floor(maphdr.Nside * (2. - 1.5 * sim.lightcone[i].opening));
					maphdr.Npix = 2 * maphdr.Nside * (maphdr.Nside + 1) + (p - maphdr.Nside) * 4 * maphdr.Nside;
				}
				else if (sim.lightcone[i].opening > -1.)
				{
					p = (int)floor(maphdr.Nside * sqrt(3. + 3. * sim.lightcone[i].opening));
					maphdr.Npix = 12 * maphdr.Nside * maphdr.Nside - 2 * p * (p + 1);
					p = 4 * maphdr.Nside - 1 - p;
				}
				else
				{
					maphdr.Npix = 12 * maphdr.Nside * maphdr.Nside;
					p = 4 * maphdr.Nside - 1;
				}

				pixbatch_size[0].push_back(maphdr.Nside / maphdr.Nside_ring);

				pixbatch_delim[1].push_back(p / pixbatch_size[0].back());
				pixbatch_delim[0].push_back((pixbatch_delim[1].back() > 0) ? pixbatch_delim[1].back() - 1 : 0);
				pixbatch_delim[2].push_back(pixbatch_delim[1].back() + 1);
				pixbatch_size[1].push_back((pixbatch_size[0].back() * (pixbatch_size[0].back() + 1) + (2 * pixbatch_size[0].back() - 1 - p % pixbatch_size[0].back()) * (p % pixbatch_size[0].back())) / 2);
				pixbatch_size[2].push_back(((p % pixbatch_size[0].back() + 1) * (p % pixbatch_size[0].back())) / 2);
				pixbatch_size[0].back() *= pixbatch_size[0].back();
				for (p = 0; p < 3; p++)
				{
					if (pixbatch_delim[p].back() <= maphdr.Nside_ring)
						pixbatch_delim[p].back() = 2 * pixbatch_delim[p].back() * (pixbatch_delim[p].back() + 1);
					else if (pixbatch_delim[p].back() <= 3 * maphdr.Nside_ring)
						pixbatch_delim[p].back() = 2 * maphdr.Nside_ring * (maphdr.Nside_ring + 1) + (pixbatch_delim[p].back() - maphdr.Nside_ring) * 4 * maphdr.Nside_ring;
					else if (pixbatch_delim[p].back() < 4 * maphdr.Nside_ring)
						pixbatch_delim[p].back() = 12 * maphdr.Nside_ring * maphdr.Nside_ring - 2 * (4 * maphdr.Nside_ring - 1 - pixbatch_delim[p].back()) * (4 * maphdr.Nside_ring - pixbatch_delim[p].back());
					else
						pixbatch_delim[p].back() = 12 * maphdr.Nside_ring * maphdr.Nside_ring;
				}

				if (pixbatch_size[1].back() == pixbatch_size[0].back())
					pixbatch_delim[0].back() = pixbatch_delim[1].back();

				for (j = 0; j < 9; j++)
					pixbuf_size[j] = 0;

				for (p = 0; p < pixbatch_delim[2].back(); p++)
				{
					pix2vec_ring64(maphdr.Nside_ring, p, w);

					base_pos[1] = (int)floor((maphdr.distance * (R[1][0] * w[0] + R[1][1] * w[1] + R[1][2] * w[2]) + sim.lightcone[i].vertex[1]) * sim.numpts) % sim.numpts;
					if (base_pos[1] < 0)
						base_pos[1] += sim.numpts;

					commdir[1] = phi->lattice().getRankDim1(base_pos[1]);
					j = commdir[1] * parallel.grid_size()[0];
					commdir[1] -= parallel.grid_rank()[1];

					if (commdir[1] < -1)
						commdir[1] += parallel.grid_size()[1];
					else if (commdir[1] > 1)
						commdir[1] -= parallel.grid_size()[1];

					base_pos[2] = (int)floor((maphdr.distance * (R[2][0] * w[0] + R[2][1] * w[1] + R[2][2] * w[2]) + sim.lightcone[i].vertex[2]) * sim.numpts) % sim.numpts;
					if (base_pos[2] < 0)
						base_pos[2] += sim.numpts;

					commdir[0] = phi->lattice().getRankDim0(base_pos[2]);
					j += commdir[0];
					commdir[0] -= parallel.grid_rank()[0];

					if (commdir[0] < -1)
						commdir[0] += parallel.grid_size()[0];
					else if (commdir[0] > 1)
						commdir[0] -= parallel.grid_size()[0];

					if ((io_group_size == 0 && parallel.rank() == ((shell - shell_inner) * parallel.size()) / (shell_outer + 1 - shell_inner)) || (io_group_size > 0 && shell - shell_inner == shell_write && ((pixbatch_delim[2].back() >= io_group_size && p / (pixbatch_delim[2].back() / io_group_size) < io_group_size && p / (pixbatch_delim[2].back() / io_group_size) == parallel.rank() - (shell_write * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner)) || (parallel.rank() - (shell_write * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner) == io_group_size - 1 && (pixbatch_delim[2].back() < io_group_size || p / (pixbatch_delim[2].back() / io_group_size) >= io_group_size)))))
					{
						sender_proc.push_back(j);
					}

					if (commdir[0] * commdir[0] > 1 || commdir[1] * commdir[1] > 1)
						continue;

					ring2nest64(maphdr.Nside_ring, p, &pix);
					pix *= pixbatch_size[0].back();

					if (p < pixbatch_delim[0].back())
						pixbatch_type = 0;
					else if (p < pixbatch_delim[1].back())
						pixbatch_type = 1;
					else
						pixbatch_type = 2;

					j = 3 * commdir[0] + commdir[1] + 4;

					if (pixbuf_size[j] + pixbatch_size[pixbatch_type].back() > pixbuf_reserve[j])
					{
						do
						{
							pixbuf_reserve[j] += PIXBUFFER;
						} while (pixbuf_size[j] + pixbatch_size[pixbatch_type].back() > pixbuf_reserve[j]);

						for (q = 0; q < LIGHTCONE_MAX_FIELDS; q++)
						{
							if (pixbuf[q][j] != NULL)
							{
								pixbuf[q][j] = (Real *)realloc((void *)pixbuf[q][j], sizeof(Real) * pixbuf_reserve[j]);
								if (pixbuf[q][j] == NULL)
								{
									cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " unable to allocate memory for pixelisation!" << endl;
									parallel.abortForce();
								}
							}
						}
					}

					for (q = 0; q < pixbatch_size[pixbatch_type].back(); pix++)
					{
						pixcount++;

						if (pixbatch_type)
						{
							nest2ring64(maphdr.Nside, pix, &pix2);
							if (pix2 >= maphdr.Npix)
								continue;
						}

						pix2vec_nest64(maphdr.Nside, pix, w);
						/* x-component of nvec to pixel from vertex R[0][0] * w[0] + R[0][1] * w[1] + R[0][2] * w[2] etc */
						double nvec[3] = {R[0][0] * w[0] + R[0][1] * w[1] + R[0][2] * w[2],
										  R[1][0] * w[0] + R[1][1] * w[1] + R[1][2] * w[2],
										  R[2][0] * w[0] + R[2][1] * w[1] + R[2][2] * w[2]};
						pos[0] = (maphdr.distance * (nvec[0]) + sim.lightcone[i].vertex[0]) * sim.numpts;
						pos[1] = (maphdr.distance * (nvec[1]) + sim.lightcone[i].vertex[1]) * sim.numpts;
						pos[2] = (maphdr.distance * (nvec[2]) + sim.lightcone[i].vertex[2]) * sim.numpts;

						if (pos[0] >= 0)
						{
							w[0] = modf(pos[0], &temp);
							base_pos[0] = (int)temp % sim.numpts;
						}
						else
						{
							w[0] = 1. + modf(pos[0], &temp);
							base_pos[0] = sim.numpts - 1 - (((int)-temp) % sim.numpts);
						}
						if (pos[1] >= 0)
						{
							w[1] = modf(pos[1], &temp);
							base_pos[1] = (int)temp % sim.numpts;
						}
						else
						{
							w[1] = 1. + modf(pos[1], &temp);
							base_pos[1] = sim.numpts - 1 - (((int)-temp) % sim.numpts);
						}
						if (pos[2] >= 0)
						{
							w[2] = modf(pos[2], &temp);
							base_pos[2] = (int)temp % sim.numpts;
						}
						else
						{
							w[2] = 1. + modf(pos[2], &temp);
							base_pos[2] = sim.numpts - 1 - (((int)-temp) % sim.numpts);
						}

						if (xsim.setCoord(base_pos))
						{
							if (sim.out_lightcone[i] & MASK_PHI)
							{
								*(pixbuf[LIGHTCONE_PHI_OFFSET][j] + pixbuf_size[j] + q) = (1. - w[0]) * (1. - w[1]) * ((1. - w[2]) * (*phi)(xsim) + w[2] * (*phi)(xsim + 2));
								*(pixbuf[LIGHTCONE_PHI_OFFSET][j] + pixbuf_size[j] + q) += w[0] * (1. - w[1]) * ((1. - w[2]) * (*phi)(xsim + 0) + w[2] * (*phi)(xsim + 0 + 2));
								*(pixbuf[LIGHTCONE_PHI_OFFSET][j] + pixbuf_size[j] + q) += w[0] * w[1] * ((1. - w[2]) * (*phi)(xsim + 0 + 1) + w[2] * (*phi)(xsim + 0 + 1 + 2));
								*(pixbuf[LIGHTCONE_PHI_OFFSET][j] + pixbuf_size[j] + q) += (1. - w[0]) * w[1] * ((1. - w[2]) * (*phi)(xsim + 1) + w[2] * (*phi)(xsim + 1 + 2));
							}
							if (sim.out_lightcone[i] & MASK_CHI)
							{
								*(pixbuf[LIGHTCONE_CHI_OFFSET][j] + pixbuf_size[j] + q) = (1. - w[0]) * (1. - w[1]) * ((1. - w[2]) * (*chi)(xsim) + w[2] * (*chi)(xsim + 2));
								*(pixbuf[LIGHTCONE_CHI_OFFSET][j] + pixbuf_size[j] + q) += w[0] * (1. - w[1]) * ((1. - w[2]) * (*chi)(xsim + 0) + w[2] * (*chi)(xsim + 0 + 2));
								*(pixbuf[LIGHTCONE_CHI_OFFSET][j] + pixbuf_size[j] + q) += w[0] * w[1] * ((1. - w[2]) * (*chi)(xsim + 0 + 1) + w[2] * (*chi)(xsim + 0 + 1 + 2));
								*(pixbuf[LIGHTCONE_CHI_OFFSET][j] + pixbuf_size[j] + q) += (1. - w[0]) * w[1] * ((1. - w[2]) * (*chi)(xsim + 1) + w[2] * (*chi)(xsim + 1 + 2));
							}
							if (sim.out_lightcone[i] & MASK_B)
							{
#ifdef LIGHTCONE_INTERPOLATE
								*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) = (1. - w[2]) * (1. - w[1]) * ((1. - w[0]) * (*Bi)(xsim - 0, 0) + w[0] * (*Bi)(xsim + 0, 0) + (*Bi)(xsim, 0));
								*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) += w[2] * (1. - w[1]) * ((1. - w[0]) * (*Bi)(xsim - 0 + 2, 0) + w[0] * (*Bi)(xsim + 0 + 2, 0) + (*Bi)(xsim + 2, 0));
								*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) += w[2] * w[1] * ((1. - w[0]) * (*Bi)(xsim - 0 + 1 + 2, 0) + w[0] * (*Bi)(xsim + 0 + 1 + 2, 0) + (*Bi)(xsim + 1 + 2, 0));
								*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) += (1. - w[2]) * w[1] * ((1. - w[0]) * (*Bi)(xsim - 0 + 1, 0) + w[0] * (*Bi)(xsim + 0 + 1, 0) + (*Bi)(xsim + 1, 0));

								*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) = (1. - w[2]) * (1. - w[0]) * ((1. - w[1]) * (*Bi)(xsim - 1, 1) + w[1] * (*Bi)(xsim + 1, 1) + (*Bi)(xsim, 1));
								*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) += w[2] * (1. - w[0]) * ((1. - w[1]) * (*Bi)(xsim - 1 + 2, 1) + w[1] * (*Bi)(xsim + 1 + 2, 1) + (*Bi)(xsim + 2, 1));
								*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) += w[2] * w[0] * ((1. - w[1]) * (*Bi)(xsim + 0 - 1 + 2, 1) + w[1] * (*Bi)(xsim + 0 + 1 + 2, 1) + (*Bi)(xsim + 0 + 2, 1));
								*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) += (1. - w[2]) * w[0] * ((1. - w[1]) * (*Bi)(xsim + 0 - 1, 1) + w[1] * (*Bi)(xsim + 0 + 1, 1) + (*Bi)(xsim + 0, 1));

								*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) = (1. - w[0]) * (1. - w[1]) * ((1. - w[2]) * (*Bi)(xsim - 2, 2) + w[2] * (*Bi)(xsim + 2, 2) + (*Bi)(xsim, 2));
								*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) += w[0] * (1. - w[1]) * ((1. - w[2]) * (*Bi)(xsim + 0 - 2, 2) + w[2] * (*Bi)(xsim + 0 + 2, 2) + (*Bi)(xsim + 0, 2));
								*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) += w[0] * w[1] * ((1. - w[2]) * (*Bi)(xsim + 0 + 1 - 2, 2) + w[2] * (*Bi)(xsim + 0 + 1 + 2, 2) + (*Bi)(xsim + 0 + 1, 2));
								*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) += (1. - w[0]) * w[1] * ((1. - w[2]) * (*Bi)(xsim + 1 - 2, 2) + w[2] * (*Bi)(xsim + 1 + 2, 2) + (*Bi)(xsim + 1, 2));

								*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) /= 2. * a * a * sim.numpts;
								*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) /= 2. * a * a * sim.numpts;
								*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) /= 2. * a * a * sim.numpts;
#else
								if (w[0] > 0.5)
								{
									*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) = (1.5 - w[0]) * (1. - w[1]) * ((1. - w[2]) * (*Bi)(xsim, 0) + w[2] * (*Bi)(xsim + 2, 0));
									*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) += (1.5 - w[0]) * w[1] * ((1. - w[2]) * (*Bi)(xsim + 1, 0) + w[2] * (*Bi)(xsim + 1 + 2, 0));
									*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) += (w[0] - 0.5) * w[1] * ((1. - w[2]) * (*Bi)(xsim + 0 + 1, 0) + w[2] * (*Bi)(xsim + 0 + 1 + 2, 0));
									*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) += (w[0] - 0.5) * (1. - w[1]) * ((1. - w[2]) * (*Bi)(xsim + 0, 0) + w[2] * (*Bi)(xsim + 0 + 2, 0));
								}
								else
								{
									*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) = (0.5 + w[0]) * (1. - w[1]) * ((1. - w[2]) * (*Bi)(xsim, 0) + w[2] * (*Bi)(xsim + 2, 0));
									*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) += (0.5 + w[0]) * w[1] * ((1. - w[2]) * (*Bi)(xsim + 1, 0) + w[2] * (*Bi)(xsim + 1 + 2, 0));
									*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) += (0.5 - w[0]) * w[1] * ((1. - w[2]) * (*Bi)(xsim - 0 + 1, 0) + w[2] * (*Bi)(xsim - 0 + 1 + 2, 0));
									*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) += (0.5 - w[0]) * (1. - w[1]) * ((1. - w[2]) * (*Bi)(xsim - 0, 0) + w[2] * (*Bi)(xsim - 0 + 2, 0));
								}
								if (w[1] > 0.5)
								{
									*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) = (1. - w[0]) * (1.5 - w[1]) * ((1. - w[2]) * (*Bi)(xsim, 1) + w[2] * (*Bi)(xsim + 2, 1));
									*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) += (1. - w[0]) * (w[1] - 0.5) * ((1. - w[2]) * (*Bi)(xsim + 1, 1) + w[2] * (*Bi)(xsim + 1 + 2, 1));
									*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) += w[0] * (w[1] - 0.5) * ((1. - w[2]) * (*Bi)(xsim + 0 + 1, 1) + w[2] * (*Bi)(xsim + 0 + 1 + 2, 1));
									*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) += w[0] * (1.5 - w[1]) * ((1. - w[2]) * (*Bi)(xsim + 0, 1) + w[2] * (*Bi)(xsim + 0 + 2, 1));
								}
								else
								{
									*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) = (1. - w[0]) * (0.5 + w[1]) * ((1. - w[2]) * (*Bi)(xsim, 1) + w[2] * (*Bi)(xsim + 2, 1));
									*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) += (1. - w[0]) * (0.5 - w[1]) * ((1. - w[2]) * (*Bi)(xsim - 1, 1) + w[2] * (*Bi)(xsim - 1 + 2, 1));
									*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) += w[0] * (0.5 - w[1]) * ((1. - w[2]) * (*Bi)(xsim + 0 - 1, 1) + w[2] * (*Bi)(xsim + 0 - 1 + 2, 1));
									*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) += w[0] * (0.5 + w[1]) * ((1. - w[2]) * (*Bi)(xsim + 0, 1) + w[2] * (*Bi)(xsim + 0 + 2, 1));
								}
								if (w[2] > 0.5)
								{
									*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) = (1. - w[0]) * (1. - w[1]) * ((1.5 - w[2]) * (*Bi)(xsim, 2) + (w[2] - 0.5) * (*Bi)(xsim + 2, 2));
									*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) += (1. - w[0]) * w[1] * ((1.5 - w[2]) * (*Bi)(xsim + 1, 2) + (w[2] - 0.5) * (*Bi)(xsim + 1 + 2, 2));
									*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) += w[0] * w[1] * ((1.5 - w[2]) * (*Bi)(xsim + 0 + 1, 2) + (w[2] - 0.5) * (*Bi)(xsim + 0 + 1 + 2, 2));
									*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) += w[0] * (1. - w[1]) * ((1.5 - w[2]) * (*Bi)(xsim + 0, 2) + (w[2] - 0.5) * (*Bi)(xsim + 0 + 2, 2));
								}
								else
								{
									*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) = (1. - w[0]) * (1. - w[1]) * ((0.5 + w[2]) * (*Bi)(xsim, 2) + (0.5 - w[2]) * (*Bi)(xsim - 2, 2));
									*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) += (1. - w[0]) * w[1] * ((0.5 + w[2]) * (*Bi)(xsim + 1, 2) + (0.5 - w[2]) * (*Bi)(xsim + 1 - 2, 2));
									*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) += w[0] * w[1] * ((0.5 + w[2]) * (*Bi)(xsim + 0 + 1, 2) + (0.5 - w[2]) * (*Bi)(xsim + 0 + 1 - 2, 2));
									*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) += w[0] * (1. - w[1]) * ((0.5 + w[2]) * (*Bi)(xsim + 0, 2) + (0.5 - w[2]) * (*Bi)(xsim + 0 - 2, 2));
								}
								*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) /= a * a * sim.numpts;
								*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) /= a * a * sim.numpts;
								*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) /= a * a * sim.numpts;
#endif
							}
#ifdef GW
							if (sim.out_lightcone[i] & MASK_HIJ)
							{
								
								*(pixbuf[LIGHTCONE_SGW_OFFSET][j] + pixbuf_size[j] + q) =
									0.5 * (1. - w[0]) * (1. - w[1]) * ((1. - w[2]) *											/*xsim*/
																		   (nvec[0] * nvec[0] * (*hij_prime)(xsim, 0, 0) +		/**/
																			nvec[1] * nvec[1] * (*hij_prime)(xsim, 1, 1) +		/**/
																			nvec[2] * nvec[2] * (*hij_prime)(xsim, 2, 2) +		/**/
																			2 * (nvec[0] * nvec[1] * (*hij_prime)(xsim, 0, 1) + /**/
																				 nvec[0] * nvec[2] * (*hij_prime)(xsim, 0, 2) + /**/
																				 nvec[1] * nvec[2] * (*hij_prime)(xsim, 1, 2))) +
																	   w[2] *														/*xsim+2*/
																		   (nvec[0] * nvec[0] * (*hij_prime)(xsim + 2, 0, 0) +		/**/
																			nvec[1] * nvec[1] * (*hij_prime)(xsim + 2, 1, 1) +		/**/
																			nvec[2] * nvec[2] * (*hij_prime)(xsim + 2, 2, 2) +		/**/
																			2 * (nvec[0] * nvec[1] * (*hij_prime)(xsim + 2, 0, 1) + /**/
																				 nvec[0] * nvec[2] * (*hij_prime)(xsim + 2, 0, 2) + /**/
																				 nvec[1] * nvec[2] * (*hij_prime)(xsim + 2, 1, 2))));
								*(pixbuf[LIGHTCONE_SGW_OFFSET][j] + pixbuf_size[j] + q) +=
									0.5 * w[0] * (1. - w[1]) * ((1. - w[2]) *												 /*xsim+0*/
																	(nvec[0] * nvec[0] * (*hij_prime)(xsim + 0, 0, 0) +		 /**/
																	 nvec[1] * nvec[1] * (*hij_prime)(xsim + 0, 1, 1) +		 /**/
																	 nvec[2] * nvec[2] * (*hij_prime)(xsim + 0, 2, 2) +		 /**/
																	 2 * (nvec[0] * nvec[1] * (*hij_prime)(xsim + 0, 0, 1) + /**/
																		  nvec[0] * nvec[2] * (*hij_prime)(xsim + 0, 0, 2) + /**/
																		  nvec[1] * nvec[2] * (*hij_prime)(xsim + 0, 1, 2))) +
																w[2] *															 /*xsim+0+2*/
																	(nvec[0] * nvec[0] * (*hij_prime)(xsim + 0 + 2, 0, 0) +		 /**/
																	 nvec[1] * nvec[1] * (*hij_prime)(xsim + 0 + 2, 1, 1) +		 /**/
																	 nvec[2] * nvec[2] * (*hij_prime)(xsim + 0 + 2, 2, 2) +		 /**/
																	 2 * (nvec[0] * nvec[1] * (*hij_prime)(xsim + 0 + 2, 0, 1) + /**/
																		  nvec[0] * nvec[2] * (*hij_prime)(xsim + 0 + 2, 0, 2) + /**/
																		  nvec[1] * nvec[2] * (*hij_prime)(xsim + 0 + 2, 1, 2))));
								*(pixbuf[LIGHTCONE_SGW_OFFSET][j] + pixbuf_size[j] + q) +=
									0.5 * w[0] * w[1] * ((1. - w[2]) *													  /*xsim+0+1*/
															 (nvec[0] * nvec[0] * (*hij_prime)(xsim + 0 + 1, 0, 0) +	  /**/
															  nvec[1] * nvec[1] * (*hij_prime)(xsim + 0 + 1, 1, 1) +	  /**/
															  nvec[2] * nvec[2] * (*hij_prime)(xsim + 0 + 1, 2, 2) +	  /**/
															  2 * (nvec[0] * nvec[1] * (*hij_prime)(xsim + 0 + 1, 0, 1) + /**/
																   nvec[0] * nvec[2] * (*hij_prime)(xsim + 0 + 1, 0, 2) + /**/
																   nvec[1] * nvec[2] * (*hij_prime)(xsim + 0 + 1, 1, 2))) +
														 w[2] *																  /*xsim+0+1+2*/
															 (nvec[0] * nvec[0] * (*hij_prime)(xsim + 0 + 1 + 2, 0, 0) +	  /**/
															  nvec[1] * nvec[1] * (*hij_prime)(xsim + 0 + 1 + 2, 1, 1) +	  /**/
															  nvec[2] * nvec[2] * (*hij_prime)(xsim + 0 + 1 + 2, 2, 2) +	  /**/
															  2 * (nvec[0] * nvec[1] * (*hij_prime)(xsim + 0 + 1 + 2, 0, 1) + /**/
																   nvec[0] * nvec[2] * (*hij_prime)(xsim + 0 + 1 + 2, 0, 2) + /**/
																   nvec[1] * nvec[2] * (*hij_prime)(xsim + 0 + 1 + 2, 1, 2))));
								*(pixbuf[LIGHTCONE_SGW_OFFSET][j] + pixbuf_size[j] + q) +=
									0.5 * (1. - w[0]) * w[1] * ((1. - w[2]) *												 /*xsim+1*/
																	(nvec[0] * nvec[0] * (*hij_prime)(xsim + 1, 0, 0) +		 /**/
																	 nvec[1] * nvec[1] * (*hij_prime)(xsim + 1, 1, 1) +		 /**/
																	 nvec[2] * nvec[2] * (*hij_prime)(xsim + 1, 2, 2) +		 /**/
																	 2 * (nvec[0] * nvec[1] * (*hij_prime)(xsim + 1, 0, 1) + /**/
																		  nvec[0] * nvec[2] * (*hij_prime)(xsim + 1, 0, 2) + /**/
																		  nvec[1] * nvec[2] * (*hij_prime)(xsim + 1, 1, 2))) +
																w[2] *															 /*xsim+1+2*/
																	(nvec[0] * nvec[0] * (*hij_prime)(xsim + 1 + 2, 0, 0) +		 /**/
																	 nvec[1] * nvec[1] * (*hij_prime)(xsim + 1 + 2, 1, 1) +		 /**/
																	 nvec[2] * nvec[2] * (*hij_prime)(xsim + 1 + 2, 2, 2) +		 /**/
																	 2 * (nvec[0] * nvec[1] * (*hij_prime)(xsim + 1 + 2, 0, 1) + /**/
																		  nvec[0] * nvec[2] * (*hij_prime)(xsim + 1 + 2, 0, 2) + /**/
																		  nvec[1] * nvec[2] * (*hij_prime)(xsim + 1 + 2, 1, 2))));
							}
#endif
							
							// Modification
							if (sim.out_lightcone[i] & MASK_ACHI && sim.aMG <= a)
							{
								*(pixbuf[LIGHTCONE_ACHI_OFFSET][j] + pixbuf_size[j] + q) = (1. - w[0]) * (1. - w[1]) * ((1. - w[2]) * (*achi)(xsim) + w[2] * (*achi)(xsim + 2));
								*(pixbuf[LIGHTCONE_ACHI_OFFSET][j] + pixbuf_size[j] + q) += w[0] * (1. - w[1]) * ((1. - w[2]) * (*achi)(xsim + 0) + w[2] * (*achi)(xsim + 0 + 2));
								*(pixbuf[LIGHTCONE_ACHI_OFFSET][j] + pixbuf_size[j] + q) += w[0] * w[1] * ((1. - w[2]) * (*achi)(xsim + 0 + 1) + w[2] * (*achi)(xsim + 0 + 1 + 2));
								*(pixbuf[LIGHTCONE_ACHI_OFFSET][j] + pixbuf_size[j] + q) += (1. - w[0]) * w[1] * ((1. - w[2]) * (*achi)(xsim + 1) + w[2] * (*achi)(xsim + 1 + 2));
							}
							if (sim.out_lightcone[i] & MASK_AQ && sim.aMG <= a)
							{
								*(pixbuf[LIGHTCONE_AQ_OFFSET][j] + pixbuf_size[j] + q) = (1. - w[0]) * (1. - w[1]) * ((1. - w[2]) * (*aq)(xsim) + w[2] * (*aq)(xsim + 2));
								*(pixbuf[LIGHTCONE_AQ_OFFSET][j] + pixbuf_size[j] + q) += w[0] * (1. - w[1]) * ((1. - w[2]) * (*aq)(xsim + 0) + w[2] * (*aq)(xsim + 0 + 2));
								*(pixbuf[LIGHTCONE_AQ_OFFSET][j] + pixbuf_size[j] + q) += w[0] * w[1] * ((1. - w[2]) * (*aq)(xsim + 0 + 1) + w[2] * (*aq)(xsim + 0 + 1 + 2));
								*(pixbuf[LIGHTCONE_AQ_OFFSET][j] + pixbuf_size[j] + q) += (1. - w[0]) * w[1] * ((1. - w[2]) * (*aq)(xsim + 1) + w[2] * (*aq)(xsim + 1 + 2));
							}
						}
						else
						{
							if (sim.out_lightcone[i] & MASK_PHI)
								*(pixbuf[LIGHTCONE_PHI_OFFSET][j] + pixbuf_size[j] + q) = 0;
							if (sim.out_lightcone[i] & MASK_CHI)
								*(pixbuf[LIGHTCONE_CHI_OFFSET][j] + pixbuf_size[j] + q) = 0;
							if (sim.out_lightcone[i] & MASK_B)
							{
								*(pixbuf[LIGHTCONE_B_OFFSET][j] + pixbuf_size[j] + q) = 0;
								*(pixbuf[LIGHTCONE_B_OFFSET + 1][j] + pixbuf_size[j] + q) = 0;
								*(pixbuf[LIGHTCONE_B_OFFSET + 2][j] + pixbuf_size[j] + q) = 0;
							}
							
#ifdef GW
							if (sim.out_lightcone[i] & MASK_HIJ)
								*(pixbuf[LIGHTCONE_SGW_OFFSET][j] + pixbuf_size[j] + q) = 0;
#endif
							// Modification
							if (sim.out_lightcone[i] & MASK_ACHI && sim.aMG <= a)
								*(pixbuf[LIGHTCONE_ACHI_OFFSET][j] + pixbuf_size[j] + q) = 0;
							if (sim.out_lightcone[i] & MASK_AQ && sim.aMG <= a)
								*(pixbuf[LIGHTCONE_AQ_OFFSET][j] + pixbuf_size[j] + q) = 0;
						}

						q++;
					} // q-loop

					pixbuf_size[j] += pixbatch_size[pixbatch_type].back();

					if (j == 4)
					{
						pixbatch_id.push_back(p);
					}
				} // p-loop
				p = 0;
				for (j = 0; j < 3; j++)
				{
					if (pixbuf_size[3 * j] + pixbuf_size[3 * j + 1] + pixbuf_size[3 * j + 2] > p)
						p = pixbuf_size[3 * j] + pixbuf_size[3 * j + 1] + pixbuf_size[3 * j + 2];
				}

				if (p > 0)
				{
					commbuf = (Real *)malloc(sizeof(Real) * p);

					for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
					{
						if (pixbuf[j][0] != NULL)
						{
							if (parallel.grid_rank()[0] % 2 == 0)
							{
								if (pixbuf_size[0] + pixbuf_size[1] + pixbuf_size[2] > 0)
								{
									if (pixbuf_size[0] > 0)
										memcpy((void *)commbuf, (void *)pixbuf[j][0], pixbuf_size[0] * sizeof(Real));
									if (pixbuf_size[1] > 0)
										memcpy((void *)(commbuf + pixbuf_size[0]), (void *)pixbuf[j][1], pixbuf_size[1] * sizeof(Real));
									if (pixbuf_size[2] > 0)
										memcpy((void *)(commbuf + pixbuf_size[0] + pixbuf_size[1]), (void *)pixbuf[j][2], pixbuf_size[2] * sizeof(Real));
									parallel.send_dim0<Real>(commbuf, pixbuf_size[0] + pixbuf_size[1] + pixbuf_size[2], (parallel.grid_size()[0] + parallel.grid_rank()[0] - 1) % parallel.grid_size()[0]);
								}
								if (pixbuf_size[3] + pixbuf_size[4] + pixbuf_size[5] > 0)
								{
									parallel.receive_dim0<Real>(commbuf, pixbuf_size[3] + pixbuf_size[4] + pixbuf_size[5], (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
									for (q = 0; q < pixbuf_size[3]; q++)
										*(pixbuf[j][3] + q) += commbuf[q];
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4] + q) += commbuf[pixbuf_size[3] + q];
									for (q = 0; q < pixbuf_size[5]; q++)
										*(pixbuf[j][5] + q) += commbuf[pixbuf_size[3] + pixbuf_size[4] + q];
								}
							}
							else
							{
								if (pixbuf_size[3] + pixbuf_size[4] + pixbuf_size[5] > 0 && parallel.grid_size()[0] > 2)
								{
									parallel.receive_dim0<Real>(commbuf, pixbuf_size[3] + pixbuf_size[4] + pixbuf_size[5], (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
									for (q = 0; q < pixbuf_size[3]; q++)
										*(pixbuf[j][3] + q) += commbuf[q];
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4] + q) += commbuf[pixbuf_size[3] + q];
									for (q = 0; q < pixbuf_size[5]; q++)
										*(pixbuf[j][5] + q) += commbuf[pixbuf_size[3] + pixbuf_size[4] + q];
								}
								if (pixbuf_size[0] + pixbuf_size[1] + pixbuf_size[2] > 0)
								{
									if (pixbuf_size[0] > 0)
										memcpy((void *)commbuf, (void *)pixbuf[j][0], pixbuf_size[0] * sizeof(Real));
									if (pixbuf_size[1] > 0)
										memcpy((void *)(commbuf + pixbuf_size[0]), (void *)pixbuf[j][1], pixbuf_size[1] * sizeof(Real));
									if (pixbuf_size[2] > 0)
										memcpy((void *)(commbuf + pixbuf_size[0] + pixbuf_size[1]), (void *)pixbuf[j][2], pixbuf_size[2] * sizeof(Real));
									parallel.send_dim0<Real>(commbuf, pixbuf_size[0] + pixbuf_size[1] + pixbuf_size[2], (parallel.grid_size()[0] + parallel.grid_rank()[0] - 1) % parallel.grid_size()[0]);
								}
							}

							if (parallel.grid_rank()[0] % 2 == 0)
							{
								if (pixbuf_size[6] + pixbuf_size[7] + pixbuf_size[8] > 0)
								{
									if (pixbuf_size[6] > 0)
										memcpy((void *)commbuf, (void *)pixbuf[j][6], pixbuf_size[6] * sizeof(Real));
									if (pixbuf_size[7] > 0)
										memcpy((void *)(commbuf + pixbuf_size[6]), (void *)pixbuf[j][7], pixbuf_size[7] * sizeof(Real));
									if (pixbuf_size[8] > 0)
										memcpy((void *)(commbuf + pixbuf_size[6] + pixbuf_size[7]), (void *)pixbuf[j][8], pixbuf_size[8] * sizeof(Real));
									parallel.send_dim0<Real>(commbuf, pixbuf_size[6] + pixbuf_size[7] + pixbuf_size[8], (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
								}
								if (pixbuf_size[3] + pixbuf_size[4] + pixbuf_size[5] > 0 && parallel.grid_size()[0] > 2)
								{
									parallel.receive_dim0<Real>(commbuf, pixbuf_size[3] + pixbuf_size[4] + pixbuf_size[5], (parallel.grid_size()[0] + parallel.grid_rank()[0] - 1) % parallel.grid_size()[0]);
									for (q = 0; q < pixbuf_size[3]; q++)
										*(pixbuf[j][3] + q) += commbuf[q];
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4] + q) += commbuf[pixbuf_size[3] + q];
									for (q = 0; q < pixbuf_size[5]; q++)
										*(pixbuf[j][5] + q) += commbuf[pixbuf_size[3] + pixbuf_size[4] + q];
								}
							}
							else
							{
								if (pixbuf_size[3] + pixbuf_size[4] + pixbuf_size[5] > 0)
								{
									parallel.receive_dim0<Real>(commbuf, pixbuf_size[3] + pixbuf_size[4] + pixbuf_size[5], (parallel.grid_size()[0] + parallel.grid_rank()[0] - 1) % parallel.grid_size()[0]);
									for (q = 0; q < pixbuf_size[3]; q++)
										*(pixbuf[j][3] + q) += commbuf[q];
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4] + q) += commbuf[pixbuf_size[3] + q];
									for (q = 0; q < pixbuf_size[5]; q++)
										*(pixbuf[j][5] + q) += commbuf[pixbuf_size[3] + pixbuf_size[4] + q];
								}
								if (pixbuf_size[6] + pixbuf_size[7] + pixbuf_size[8] > 0)
								{
									if (pixbuf_size[6] > 0)
										memcpy((void *)commbuf, (void *)pixbuf[j][6], pixbuf_size[6] * sizeof(Real));
									if (pixbuf_size[7] > 0)
										memcpy((void *)(commbuf + pixbuf_size[6]), (void *)pixbuf[j][7], pixbuf_size[7] * sizeof(Real));
									if (pixbuf_size[8] > 0)
										memcpy((void *)(commbuf + pixbuf_size[6] + pixbuf_size[7]), (void *)pixbuf[j][8], pixbuf_size[8] * sizeof(Real));
									parallel.send_dim0<Real>(commbuf, pixbuf_size[6] + pixbuf_size[7] + pixbuf_size[8], (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
								}
							}

							if (parallel.grid_rank()[1] % 2 == 0)
							{
								if (pixbuf_size[3] > 0)
									parallel.send_dim1<Real>(pixbuf[j][3], pixbuf_size[3], (parallel.grid_size()[1] + parallel.grid_rank()[1] - 1) % parallel.grid_size()[1]);
								if (pixbuf_size[4] > 0)
								{
									parallel.receive_dim1<Real>(commbuf, pixbuf_size[4], (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4] + q) += commbuf[q];
								}
							}
							else
							{
								if (pixbuf_size[4] > 0 && parallel.grid_size()[1] > 2)
								{
									parallel.receive_dim1<Real>(commbuf, pixbuf_size[4], (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4] + q) += commbuf[q];
								}
								if (pixbuf_size[3] > 0)
									parallel.send_dim1<Real>(pixbuf[j][3], pixbuf_size[3], (parallel.grid_size()[1] + parallel.grid_rank()[1] - 1) % parallel.grid_size()[1]);
							}

							if (parallel.grid_rank()[1] % 2 == 0)
							{
								if (pixbuf_size[5] > 0)
									parallel.send_dim1<Real>(pixbuf[j][5], pixbuf_size[5], (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
								if (pixbuf_size[4] > 0 && parallel.grid_size()[1] > 2)
								{
									parallel.receive_dim1<Real>(commbuf, pixbuf_size[4], (parallel.grid_size()[1] + parallel.grid_rank()[1] - 1) % parallel.grid_size()[1]);
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4] + q) += commbuf[q];
								}
							}
							else
							{
								if (pixbuf_size[4] > 0)
								{
									parallel.receive_dim1<Real>(commbuf, pixbuf_size[4], (parallel.grid_size()[1] + parallel.grid_rank()[1] - 1) % parallel.grid_size()[1]);
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4] + q) += commbuf[q];
								}
								if (pixbuf_size[5] > 0)
									parallel.send_dim1<Real>(pixbuf[j][5], pixbuf_size[5], (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
							}
						}
					}

					free(commbuf);
				}

				if (io_group_size == 0 && parallel.rank() == ((shell - shell_inner) * parallel.size() / (shell_outer + 1 - shell_inner)))
				{
					for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
					{
						if (pixbuf[j][0] != NULL)
						{
							if (bytes2 == 0)
							{
								outbuf[j] = (char *)malloc(maphdr.Npix * maphdr.precision + 272);

								if (outbuf[j] == NULL)
								{
									cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " unable to allocate " << maphdr.Npix * maphdr.precision + 272 << " bytes of memory for pixelisation!" << endl;
									parallel.abortForce();
								}
							}
							else
							{
								outbuf[j] = (char *)realloc((void *)outbuf[j], bytes2 + maphdr.Npix * maphdr.precision + 272);

								if (outbuf[j] == NULL)
								{
									cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " unable to reallocate " << bytes2 + maphdr.Npix * maphdr.precision + 272 << " bytes of memory (" << maphdr.Npix * maphdr.precision + 272 << " additional bytes) for pixelisation!" << endl;
									parallel.abortForce();
								}
							}

							blocksize = 256;
							memcpy((void *)(outbuf[j] + bytes2), (void *)&blocksize, 4);
							memcpy((void *)(outbuf[j] + bytes2 + 4), (void *)&maphdr, 256);
							memcpy((void *)(outbuf[j] + bytes2 + 260), (void *)&blocksize, 4);
							blocksize = maphdr.precision * maphdr.Npix;
							memcpy((void *)(outbuf[j] + bytes2 + 264), (void *)&blocksize, 4);
							memcpy((void *)(outbuf[j] + bytes2 + 268 + blocksize), (void *)&blocksize, 4);
						}
					}
					offset2 = bytes2 + 268;
					bytes2 += maphdr.Npix * maphdr.precision + 272;
					p = 0;
					q = pixbatch_delim[2].back();
				}
				else if (io_group_size > 0 && shell - shell_inner == shell_write)
				{
					q = pixbatch_delim[2].back() / io_group_size;
					p = parallel.rank() - (shell_write * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner);

					for (j = 0; p * q >= pixbatch_delim[j].back(); j++)
						;

					if ((p + 1) * q >= pixbatch_delim[j].back() && j < 2)
					{
						bytes2 = (pixbatch_delim[j].back() - p * q) * pixbatch_size[j].back() * maphdr.precision;
						if ((p + 1) * q >= pixbatch_delim[j + 1].back() && j < 1)
						{
							bytes2 += (pixbatch_delim[j + 1].back() - pixbatch_delim[j].back()) * pixbatch_size[j + 1].back() * maphdr.precision;
							bytes2 += ((p + 1) * q - pixbatch_delim[j + 1].back()) * pixbatch_size[j + 2].back() * maphdr.precision;
						}
						else
							bytes2 += ((p + 1) * q - pixbatch_delim[j].back()) * pixbatch_size[j + 1].back() * maphdr.precision;
					}
					else
						bytes2 = q * pixbatch_size[j].back() * maphdr.precision;

					if (p == 0)
					{
						bytes2 += 268;
						offset2 = 268;
					}
					else
						offset2 = 0;

					if (p == io_group_size - 1)
					{
						bytes2 += 4;
						q = pixbatch_delim[2].back() % io_group_size;
						if (pixbatch_delim[2].back() - pixbatch_delim[1].back() < q)
						{
							bytes2 += (pixbatch_delim[2].back() - pixbatch_delim[1].back()) * pixbatch_size[2].back() * maphdr.precision;
							if (pixbatch_delim[2].back() - pixbatch_delim[0].back() < q)
							{
								bytes2 += (pixbatch_delim[1].back() - pixbatch_delim[0].back()) * pixbatch_size[1].back() * maphdr.precision;
								bytes2 += (q - pixbatch_delim[2].back() + pixbatch_delim[0].back()) * pixbatch_size[0].back() * maphdr.precision;
							}
							else
								bytes2 += (q - pixbatch_delim[2].back() + pixbatch_delim[1].back()) * pixbatch_size[1].back() * maphdr.precision;
						}
						else
							bytes2 += q * pixbatch_size[2].back() * maphdr.precision;

						q += pixbatch_delim[2].back() / io_group_size;
					}

					for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
					{
						if (pixbuf[j][0] != NULL)
						{
							if (bytes2 > 0)
							{
								outbuf[j] = (char *)malloc(bytes2);

								if (outbuf[j] == NULL)
								{
									cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " unable to allocate " << bytes2 << " bytes of memory for pixelisation!" << endl;
									parallel.abortForce();
								}
							}

							if (p == 0)
							{
								blocksize = 256;
								memcpy((void *)outbuf[j], (void *)&blocksize, 4);
								memcpy((void *)(outbuf[j] + 4), (void *)&maphdr, 256);
								memcpy((void *)(outbuf[j] + 260), (void *)&blocksize, 4);
								blocksize = maphdr.precision * maphdr.Npix;
								memcpy((void *)(outbuf[j] + 264), (void *)&blocksize, 4);
							}

							if (p == io_group_size - 1)
							{
								blocksize = maphdr.precision * maphdr.Npix;
								memcpy((void *)(outbuf[j] + bytes2 - 4), (void *)&blocksize, 4);
							}
						}
					}

					p *= pixbatch_delim[2].back() / io_group_size;
				}

				pix = 0;
				pix2 = 0;

				if ((io_group_size == 0 && parallel.rank() == ((shell - shell_inner) * parallel.size()) / (shell_outer + 1 - shell_inner)) || (io_group_size > 0 && shell - shell_inner == shell_write))
				{
					if (q != sender_proc.size())
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " pixel batch count mismatch! expecting " << q << " but sender list contains " << sender_proc.size() << " entries!" << endl;
						exit(-99);
					}

					for (int64_t p2 = p; p2 < p + q; p2 += n)
					{
						while (pix < pixbatch_id.size() && pixbatch_id[pix] < p2)
						{
							for (pixbatch_type = 0; pixbatch_delim[pixbatch_type].back() <= pixbatch_id[pix]; pixbatch_type++)
								;
							if (io_group_size > 0 && pixbatch_delim[2].back() >= io_group_size && pixbatch_id[pix] / (pixbatch_delim[2].back() / io_group_size) < io_group_size)
							{
								for (n = 1; pix + n < pixbatch_id.size() && pixbatch_id[pix + n] == pixbatch_id[pix + n - 1] + 1 && pixbatch_id[pix + n] < pixbatch_delim[pixbatch_type].back() && (pixbatch_id[pix + n] / (pixbatch_delim[2].back() / io_group_size) == pixbatch_id[pix] / (pixbatch_delim[2].back() / io_group_size) || pixbatch_id[pix] / (pixbatch_delim[2].back() / io_group_size) == io_group_size - 1); n++)
									;
								for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
								{
									if (pixbuf[j][4] != NULL && pixbatch_size[pixbatch_type].back() > 0)
										parallel.send<Real>(pixbuf[j][4] + pix2, n * pixbatch_size[pixbatch_type].back(), (pixbatch_id[pix] / (pixbatch_delim[2].back() / io_group_size)) + ((shell - shell_inner) * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner));
								}
							}
							else
							{
								for (n = 1; pix + n < pixbatch_id.size() && pixbatch_id[pix + n] == pixbatch_id[pix + n - 1] + 1 && pixbatch_id[pix + n] < pixbatch_delim[pixbatch_type].back(); n++)
									;
								for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
								{
									if (pixbuf[j][4] != NULL && pixbatch_size[pixbatch_type].back() > 0)
										parallel.send<Real>(pixbuf[j][4] + pix2, n * pixbatch_size[pixbatch_type].back(), (io_group_size ? io_group_size - 1 : 0) + ((shell - shell_inner) * parallel.size() + (io_group_size ? shell_outer - shell_inner : 0)) / (shell_outer + 1 - shell_inner));
								}
							}
							pix += n;
							pix2 += n * pixbatch_size[pixbatch_type].back();
						}

						for (pixbatch_type = 0; pixbatch_delim[pixbatch_type].back() <= p2; pixbatch_type++)
							;

						for (n = 1; p2 + n < p + q && sender_proc[p2 + n - p] == sender_proc[p2 - p] && p2 + n < pixbatch_delim[pixbatch_type].back(); n++)
							;

						if (sender_proc[p2 - p] == parallel.rank())
						{
							if (pix + n - 1 >= pixbatch_id.size())
							{
								cerr << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " pixel batch index mismatch! expecting " << p2 << " but ID list contains not enough elements!" << endl;
								exit(-99);
							}
							else if (pixbatch_id[pix] != p2)
							{
								cerr << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " pixel batch index mismatch! expecting " << p2 << " but ID list says " << pixbatch_id[pix] << "!" << endl;
								exit(-99);
							}
							for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
							{
								if (pixbuf[j][4] != NULL && pixbatch_size[pixbatch_type].back() > 0)
									memcpy((void *)(outbuf[j] + offset2), (void *)(pixbuf[j][4] + pix2), n * pixbatch_size[pixbatch_type].back() * maphdr.precision);
							}
							pix += n;
							pix2 += n * pixbatch_size[pixbatch_type].back();
						}
						else
						{
							for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
							{
								if (outbuf[j] != NULL && pixbatch_size[pixbatch_type].back() > 0)
									parallel.receive<Real>((Real *)(outbuf[j] + offset2), n * pixbatch_size[pixbatch_type].back(), sender_proc[p2 - p]);
							}
						}

						offset2 += n * pixbatch_size[pixbatch_type].back() * maphdr.precision;
					}

					if (io_group_size > 0)
					{
						if (p > 0)
						{
							if (p >= pixbatch_delim[0].back())
							{
								offset2 = 268 + pixbatch_delim[0].back() * pixbatch_size[0].back() * maphdr.precision;
								if (p >= pixbatch_delim[1].back())
								{
									offset2 += (pixbatch_delim[1].back() - pixbatch_delim[0].back()) * pixbatch_size[1].back() * maphdr.precision;
									offset2 += (p - pixbatch_delim[1].back()) * pixbatch_size[2].back() * maphdr.precision;
								}
								else
									offset2 += (p - pixbatch_delim[0].back()) * pixbatch_size[1].back() * maphdr.precision;
							}
							else
								offset2 = 268 + p * pixbatch_size[0].back() * maphdr.precision;
						}
						else if (parallel.rank() == (shell_write * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner))
							offset2 = 0;
						else
							offset2 = 268;
					}
				}

				p = ((((shell + 1 - shell_inner) * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner)) - (((shell - shell_inner) * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner)));

				while (pix < pixbatch_id.size())
				{
					for (pixbatch_type = 0; pixbatch_delim[pixbatch_type].back() <= pixbatch_id[pix]; pixbatch_type++)
						;

					if (p > 0 && pixbatch_delim[2].back() >= p && pixbatch_id[pix] / (pixbatch_delim[2].back() / p) < p)
					{
						for (n = 1; pix + n < pixbatch_id.size() && pixbatch_id[pix + n] == pixbatch_id[pix + n - 1] + 1 && pixbatch_id[pix + n] < pixbatch_delim[pixbatch_type].back() && (pixbatch_id[pix + n] / (pixbatch_delim[2].back() / p) == pixbatch_id[pix] / (pixbatch_delim[2].back() / p) || pixbatch_id[pix] / (pixbatch_delim[2].back() / p) == p - 1); n++)
							;
						for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
						{
							if (pixbuf[j][4] != NULL && pixbatch_size[pixbatch_type].back() > 0)
								parallel.send<Real>(pixbuf[j][4] + pix2, n * pixbatch_size[pixbatch_type].back(), (pixbatch_id[pix] / (pixbatch_delim[2].back() / p)) + ((shell - shell_inner) * parallel.size() + (io_group_size ? shell_outer - shell_inner : 0)) / (shell_outer + 1 - shell_inner));
						}
					}
					else
					{
						for (n = 1; pix + n < pixbatch_id.size() && pixbatch_id[pix + n] == pixbatch_id[pix + n - 1] + 1 && pixbatch_id[pix + n] < pixbatch_delim[pixbatch_type].back(); n++)
							;
						for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
						{
							if (pixbuf[j][4] != NULL && pixbatch_size[pixbatch_type].back() > 0)
								parallel.send<Real>(pixbuf[j][4] + pix2, n * pixbatch_size[pixbatch_type].back(), (p ? p - 1 : 0) + ((shell - shell_inner) * parallel.size() + (p ? shell_outer - shell_inner : 0)) / (shell_outer + 1 - shell_inner));
						}
					}

					pix += n;
					pix2 += n * pixbatch_size[pixbatch_type].back();
				}

				offset.push_back(bytes);
				bytes += maphdr.Npix * maphdr.precision + 272;

				pixbatch_id.clear();
				sender_proc.clear();
			} // shell-loop

			if (io_group_size == 0)
				offset2 = 0;

			for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
			{
				if (pixbuf[j][0] == NULL || shell_outer < shell_inner)
					continue;

				if (sim.num_lightcone > 1)
				{
					if (j == LIGHTCONE_PHI_OFFSET)
						snprintf(filename, SIZFILENAME, "%s%s%d_%04d_phi.map", sim.output_path, sim.basename_lightcone, i, cycle);
					else if (j == LIGHTCONE_CHI_OFFSET)
						snprintf(filename, SIZFILENAME, "%s%s%d_%04d_chi.map", sim.output_path, sim.basename_lightcone, i, cycle);
					else if (j >= LIGHTCONE_B_OFFSET && j < LIGHTCONE_B_OFFSET + 3)
						snprintf(filename, SIZFILENAME, "%s%s%d_%04d_B%d.map", sim.output_path, sim.basename_lightcone, i, cycle, j + 1 - LIGHTCONE_B_OFFSET);
					if (j == LIGHTCONE_SGW_OFFSET)
						snprintf(filename, SIZFILENAME, "%s%s%d_%04d_sgw.map", sim.output_path, sim.basename_lightcone, i, cycle);
					// Modification
					else if (j == LIGHTCONE_ACHI_OFFSET)
						snprintf(filename, SIZFILENAME, "%s%s%d_%04d_achi.map", sim.output_path, sim.basename_lightcone, i, cycle);
					else if (j == LIGHTCONE_AQ_OFFSET)
						snprintf(filename, SIZFILENAME, "%s%s%d_%04d_aq.map", sim.output_path, sim.basename_lightcone, i, cycle);
					else if (j == LIGHTCONE_T00_OFFSET)
						snprintf(filename, SIZFILENAME, "%s%s%d_%04d_T00.map", sim.output_path, sim.basename_lightcone, i, cycle);
				}
				else
				{
					if (j == LIGHTCONE_PHI_OFFSET)
						snprintf(filename, SIZFILENAME, "%s%s_%04d_phi.map", sim.output_path, sim.basename_lightcone, cycle);
					else if (j == LIGHTCONE_CHI_OFFSET)
						snprintf(filename, SIZFILENAME, "%s%s_%04d_chi.map", sim.output_path, sim.basename_lightcone, cycle);
					else if (j >= LIGHTCONE_B_OFFSET && j < LIGHTCONE_B_OFFSET + 3)
						snprintf(filename, SIZFILENAME, "%s%s_%04d_B%d.map", sim.output_path, sim.basename_lightcone, cycle, j + 1 - LIGHTCONE_B_OFFSET);
					if (j == LIGHTCONE_SGW_OFFSET)
						snprintf(filename, SIZFILENAME, "%s%s_%04d_sgw.map", sim.output_path, sim.basename_lightcone, cycle);
					// Modification
					else if (j == LIGHTCONE_ACHI_OFFSET)
						snprintf(filename, SIZFILENAME, "%s%s_%04d_achi.map", sim.output_path, sim.basename_lightcone, cycle);
					else if (j == LIGHTCONE_AQ_OFFSET)
						snprintf(filename, SIZFILENAME, "%s%s_%04d_aq.map", sim.output_path, sim.basename_lightcone, cycle);
				}

				MPI_File_open(parallel.lat_world_comm(), filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mapfile);
				MPI_File_set_size(mapfile, (MPI_Offset)bytes);
				MPI_File_write_at_all(mapfile, (MPI_Offset)offset[shell_write] + offset2, (void *)outbuf[j], bytes2, MPI_BYTE, &status);
				MPI_File_close(&mapfile);
			}

			for (j = 0; j < 3; j++)
			{
				pixbatch_size[j].clear();
				pixbatch_delim[j].clear();
			}

			offset.clear();

			for (j = 0; j < 9 * LIGHTCONE_MAX_FIELDS; j++)
			{
				if (pixbuf[j / 9][j % 9] != NULL)
				{
					free(pixbuf[j / 9][j % 9]);
					pixbuf[j / 9][j % 9] = NULL;
				}
			}

			for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
			{
				if (outbuf[j] != NULL)
				{
					free(outbuf[j]);
					outbuf[j] = NULL;
				}
			}
#endif // HAVE_HEALPIX
		}
		else if (parallel.isRoot() && outfile != NULL)
			fclose(outfile);

		if (sim.out_lightcone[i] & MASK_GADGET && sim.lightcone[i].distance[0] > d - tau + 0.5 * dtau_old && sim.lightcone[i].distance[1] <= d - tau + 0.5 * dtau_old && d - tau + 0.5 * dtau_old > 0.)
		{
			n = findIntersectingLightcones(sim.lightcone[i], d - tau + (0.5 + LIGHTCONE_IDCHECK_ZONE) * dtau_old, d - tau - 0.5 * dtau, domain, vertex);

			hdr.num_files = 1;
			hdr.Omega0 = cosmo.Omega_m;
			hdr.OmegaLambda = cosmo.Omega_Lambda;
			hdr.HubbleParam = cosmo.h;
			hdr.BoxSize = sim.boxsize / GADGET_LENGTH_CONVERSION;
			hdr.flag_sfr = 0;
			hdr.flag_cooling = 0;
			hdr.flag_feedback = 0;
			hdr.flag_age = 0;
			hdr.flag_metals = 0;
			for (p = 0; p < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4; p++)
				hdr.fill[p] = 0;
			for (p = 0; p < 6; p++)
			{
				hdr.npart[p] = 0;
				hdr.npartTotal[p] = 0;
				hdr.npartTotalHW[p] = 0;
				hdr.mass[p] = 0.;
			}

			hdr.time = a;
			hdr.redshift = (1. / a) - 1.;

			if (sim.baryon_flag)
				hdr.mass[1] = (double)sim.tracer_factor[0] * C_RHO_CRIT * cosmo.Omega_cdm * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;
			else
				hdr.mass[1] = (double)sim.tracer_factor[0] * C_RHO_CRIT * (cosmo.Omega_cdm + cosmo.Omega_b) * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;

			if (sim.num_lightcone > 1)
				snprintf(filename, SIZFILENAME, "%d_%04d", i, cycle);
			else
				snprintf(filename, SIZFILENAME, "_%04d", cycle);

			if (sim.tracer_factor[0] > 0)
				pcls_cdm->saveGadget2(h5filename + filename + "_cdm", hdr, sim.lightcone[i], d - tau, dtau, dtau_old, a * Hconf(a, fourpiG, cosmo, avgsource), vertex, n, IDbacklog[0], IDprelog[0], phi, sim.tracer_factor[0]);

			if (sim.baryon_flag && sim.tracer_factor[1] > 0)
			{
				hdr.mass[1] = (double)sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
				pcls_b->saveGadget2(h5filename + filename + "_b", hdr, sim.lightcone[i], d - tau, dtau, dtau_old, a * Hconf(a, fourpiG, cosmo, avgsource), vertex, n, IDbacklog[1], IDprelog[1], phi, sim.tracer_factor[1]);
			}

			for (p = 0; p < cosmo.num_ncdm; p++)
			{
				if (sim.numpcl[1 + sim.baryon_flag + p] == 0 || sim.tracer_factor[p + 1 + sim.baryon_flag] == 0)
					continue;
				snprintf(buffer, SIZBUFFER, "_ncdm%d", p);
				hdr.mass[1] = (double)sim.tracer_factor[p + 1 + sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[p] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[p + 1 + sim.baryon_flag] / GADGET_MASS_CONVERSION;
				pcls_ncdm[p].saveGadget2(h5filename + filename + buffer, hdr, sim.lightcone[i], d - tau, dtau, dtau_old, a * Hconf(a, fourpiG, cosmo, avgsource), vertex, n, IDbacklog[p + 1 + sim.baryon_flag], IDprelog[p + 1 + sim.baryon_flag], phi, sim.tracer_factor[p + 1 + sim.baryon_flag]);
			}
		}
	}

#ifdef HAVE_HEALPIX
	delete[] outbuf;
#endif

	for (p = 0; p <= cosmo.num_ncdm + sim.baryon_flag; p++)
	{
		IDbacklog[p] = IDprelog[p];
		IDprelog[p].clear();

		n = IDbacklog[p].size();
		// dim 0 send/rec
		if (parallel.grid_rank()[0] % 2 == 0)
		{
			parallel.send_dim0<int>(n, (parallel.grid_size()[0] + parallel.grid_rank()[0] - 1) % parallel.grid_size()[0]);
			parallel.receive_dim0<int>(i, (parallel.grid_size()[0] + parallel.grid_rank()[0] - 1) % parallel.grid_size()[0]);
			parallel.send_dim0<int>(n, (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
			parallel.receive_dim0<int>(j, (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
		}
		else
		{
			parallel.receive_dim0<int>(i, (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
			parallel.send_dim0<int>(n, (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
			parallel.receive_dim0<int>(j, (parallel.grid_size()[0] + parallel.grid_rank()[0] - 1) % parallel.grid_size()[0]);
			parallel.send_dim0<int>(n, (parallel.grid_size()[0] + parallel.grid_rank()[0] - 1) % parallel.grid_size()[0]);
		}

		if (n + i + j > 0)
		{
			IDcombuf = (long *)malloc((n + i + j) * sizeof(long));

			n = 0;
			for (std::set<long>::iterator it = IDbacklog[p].begin(); it != IDbacklog[p].end(); it++)
				IDcombuf[n++] = *it;

			if (parallel.grid_rank()[0] % 2 == 0)
			{
				if (n > 0)
					parallel.send_dim0<long>(IDcombuf, n, (parallel.grid_size()[0] + parallel.grid_rank()[0] - 1) % parallel.grid_size()[0]);
				if (i > 0)
					parallel.receive_dim0<long>(IDcombuf + n, i, (parallel.grid_size()[0] + parallel.grid_rank()[0] - 1) % parallel.grid_size()[0]);
				if (n > 0)
					parallel.send_dim0<long>(IDcombuf, n, (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
				if (j > 0)
					parallel.receive_dim0<long>(IDcombuf + n + i, j, (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
			}
			else
			{
				if (i > 0)
					parallel.receive_dim0<long>(IDcombuf + n, i, (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
				if (n > 0)
					parallel.send_dim0<long>(IDcombuf, n, (parallel.grid_rank()[0] + 1) % parallel.grid_size()[0]);
				if (j > 0)
					parallel.receive_dim0<long>(IDcombuf + n + i, j, (parallel.grid_size()[0] + parallel.grid_rank()[0] - 1) % parallel.grid_size()[0]);
				if (n > 0)
					parallel.send_dim0<long>(IDcombuf, n, (parallel.grid_size()[0] + parallel.grid_rank()[0] - 1) % parallel.grid_size()[0]);
			}

			n += i + j;

			for (i = IDbacklog[p].size(); i < n; i++)
				IDbacklog[p].insert(IDcombuf[i]);
		}

		// dim 1 send/rec
		if (parallel.grid_rank()[1] % 2 == 0)
		{
			parallel.send_dim1<int>(n, (parallel.grid_size()[1] + parallel.grid_rank()[1] - 1) % parallel.grid_size()[1]);
			parallel.receive_dim1<int>(i, (parallel.grid_size()[1] + parallel.grid_rank()[1] - 1) % parallel.grid_size()[1]);
			parallel.send_dim1<int>(n, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
			parallel.receive_dim1<int>(j, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);

			if (n > 0)
				parallel.send_dim1<long>(IDcombuf, n, (parallel.grid_size()[1] + parallel.grid_rank()[1] - 1) % parallel.grid_size()[1]);

			if (i > 0)
			{
				IDcombuf2 = (long *)malloc(i * sizeof(long));
				parallel.receive_dim1<long>(IDcombuf2, i, (parallel.grid_size()[1] + parallel.grid_rank()[1] - 1) % parallel.grid_size()[1]);
				while (i > 0)
					IDbacklog[p].insert(IDcombuf2[--i]);
				free(IDcombuf2);
			}

			if (n > 0)
			{
				parallel.send_dim1<long>(IDcombuf, n, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
				free(IDcombuf);
			}

			if (j > 0)
			{
				IDcombuf2 = (long *)malloc(j * sizeof(long));
				parallel.receive_dim1<long>(IDcombuf2, j, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
				while (j > 0)
					IDbacklog[p].insert(IDcombuf2[--j]);
				free(IDcombuf2);
			}
		}
		else
		{
			parallel.receive_dim1<int>(i, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
			parallel.send_dim1<int>(n, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
			parallel.receive_dim1<int>(j, (parallel.grid_size()[1] + parallel.grid_rank()[1] - 1) % parallel.grid_size()[1]);
			parallel.send_dim1<int>(n, (parallel.grid_size()[1] + parallel.grid_rank()[1] - 1) % parallel.grid_size()[1]);

			if (i > 0)
			{
				IDcombuf2 = (long *)malloc(i * sizeof(long));
				parallel.receive_dim1<long>(IDcombuf2, i, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);
				while (i > 0)
					IDbacklog[p].insert(IDcombuf2[--i]);
				free(IDcombuf2);
			}

			if (n > 0)
				parallel.send_dim1<long>(IDcombuf, n, (parallel.grid_rank()[1] + 1) % parallel.grid_size()[1]);

			if (j > 0)
			{
				IDcombuf2 = (long *)malloc(j * sizeof(long));
				parallel.receive_dim1<long>(IDcombuf2, j, (parallel.grid_size()[1] + parallel.grid_rank()[1] - 1) % parallel.grid_size()[1]);
				while (j > 0)
					IDbacklog[p].insert(IDcombuf2[--j]);
				free(IDcombuf2);
			}

			if (n > 0)
			{
				parallel.send_dim1<long>(IDcombuf, n, (parallel.grid_size()[1] + parallel.grid_rank()[1] - 1) % parallel.grid_size()[1]);
				free(IDcombuf);
			}
		}
	}
}

//////////////////////////
// writeSpectra
//////////////////////////
// Description:
//   output of spectra
//
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   a              scale factor
//   pkcount        spectrum output index
//   pcls_cdm       pointer to particle handler for CDM
//   pcls_b         pointer to particle handler for baryons
//   pcls_ncdm      array of particle handlers for
//                  non-cold DM (may be set to NULL)
//   phi            pointer to allocated field
//	 achi			pointer to allocated field
//	 aq				pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   T00_As         pointer to allocated field
//   T0i_As         pointer to allocated field
//   Tij_As         pointer to allocated field
//   source         pointer to allocated field
//   hij            pointer to allocated field
//   scalarFT       pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_phi       pointer to FFT planner
//   plan_achi      pointer to FFT planner
//   plan_aq        pointer to FFT planner
//   plan_chi       pointer to FFT planner
//   plan_Bi        pointer to FFT planner
//   plan_T00_As    pointer to FFT planner
//   plan_T0i_As    pointer to FFT planner
//   plan_Tij_As    pointer to FFT planner
//   plan_source    pointer to FFT planner
//   plan_hij       pointer to FFT planner
//   Bi_check       pointer to allocated field (or NULL)
//   BiFT_check     pointer to allocated field (or NULL)	e
//   plan_Bi_check  pointer to FFT planner (or NULL)
//
// Returns:
//
//////////////////////////

void writeSpectra(metadata &sim, cosmology &cosmo, double avgsource, const double fourpiG, const double a, const int pkcount,
#ifdef HAVE_CLASS
				  background &class_background, perturbs &class_perturbs, icsettings &ic,
#endif
				  Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_cdm,
				  Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_b,
				  Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> *pcls_ncdm,
				  Field<Real> *phi, Field<Real> *achi, Field<Real> *aq, Field<Real> *chi,
				  Field<Real> *Bi, Field<Real> *hij, Field<Real> *source, Field<Real> *scalar,
				  Field<Cplx> *scalarFT,
				  Field<Cplx> *BiFT, Field<Cplx> *SijFT, PlanFFT<Cplx> *plan_phi,
				  PlanFFT<Cplx> *plan_achi, PlanFFT<Cplx> *plan_aq,
				  PlanFFT<Cplx> *plan_chi, PlanFFT<Cplx> *plan_Bi, PlanFFT<Cplx> *plan_source,
				  PlanFFT<Cplx> *plan_Sij, PlanFFT<Cplx> *plan_scalar,
				  int &done_hij
#ifdef GW
				  ,
				  Field<Cplx> *hijFT, Field<Cplx> *hijprimeFT,
				  PlanFFT<Cplx> *plan_hij, PlanFFT<Cplx> *plan_hijprime
#endif
)
{
	char filename[2 * PARAM_MAX_LENGTH + 24];
	char buffer[64];
	int i, j;
	Site x(phi->lattice());
	rKSite kFT(scalarFT->lattice());
	long numpts3d = (long)sim.numpts * (long)sim.numpts * (long)sim.numpts;
	Cplx tempk;
	double Omega_ncdm;

	Real *kbin;
	Real *power;
	Real *kscatter;
	Real *pscatter;
	int *occupation;

	kbin = (Real *)malloc(sim.numbins * sizeof(Real));
	power = (Real *)malloc(sim.numbins * sizeof(Real));
	kscatter = (Real *)malloc(sim.numbins * sizeof(Real));
	pscatter = (Real *)malloc(sim.numbins * sizeof(Real));
	occupation = (int *)malloc(sim.numbins * sizeof(int));

	if (sim.out_pk & MASK_PHI)
	{
		plan_phi->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		snprintf(filename, SIZFILENAME, "%s%s%03d_phi.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real)numpts3d * (Real)numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of phi", a, sim.z_pk[pkcount]);
	}

	//---------------------------Modification----------------------------
	if (sim.out_pk & MASK_ACHI && sim.aMG <= a)
	{
		plan_achi->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		snprintf(filename, SIZFILENAME, "%s%s%03d_achi.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real)numpts3d * (Real)numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of achi", a, sim.z_pk[pkcount]);
	}

	if (sim.out_pk & MASK_AQ && sim.aMG <= a)
	{
		plan_aq->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		snprintf(filename, SIZFILENAME, "%s%s%03d_aq.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real)numpts3d * (Real)numpts3d * 2. * M_PI * M_PI * Hconf(a, fourpiG, cosmo, avgsource) * Hconf(a, fourpiG, cosmo, avgsource), filename, "power spectrum of aq", a, sim.z_pk[pkcount]);
	}
	if (sim.out_pk & MASK_SOURCE)
	{
		plan_source->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
		snprintf(filename, SIZFILENAME, "%s%s%03d_source.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real)numpts3d * (Real)numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of source", a, sim.z_pk[pkcount]);
	}
	//------------------------------------------------------------------------------

	if (sim.out_pk & MASK_CHI)
	{
		plan_chi->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		snprintf(filename, SIZFILENAME, "%s%s%03d_chi.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real)numpts3d * (Real)numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of chi", a, sim.z_pk[pkcount]);
	}

	if (sim.out_pk & MASK_HIJ)
	{
#ifdef GW
		//	SijFT should not have been used since Tij Fourier
		projectFTtensor(*SijFT, *SijFT); // adiabatic
#else
		if (done_hij == 0)
		{
			projectFTtensor(*SijFT, *SijFT); // adiabatic
		}
#endif
		extractPowerSpectrum(*SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		snprintf(filename, SIZFILENAME, "%s%s%03d_hij.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, 2. * M_PI * M_PI, filename, "power spectrum of hij", a, sim.z_pk[pkcount]);
#ifdef GW
		// dynamical dofs:
		extractPowerSpectrum(*hijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		snprintf(filename, SIZFILENAME, "%s%s%03d_hij_dyn.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, 2. * M_PI * M_PI, filename, "power spectrum of hij_dyn", a, sim.z_pk[pkcount]);

		extractPowerSpectrum(*hijprimeFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		snprintf(filename, SIZFILENAME, "%s%s%03d_hij_prime.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, 2. * M_PI * M_PI * Hconf(a, fourpiG, cosmo, avgsource) * Hconf(a, fourpiG, cosmo, avgsource), filename, "power spectrum of Hconf * hij_prime'", a, sim.z_pk[pkcount]);
#endif
	}

	if ((sim.out_pk & MASK_T00 || sim.out_pk & MASK_DELTA) && sim.gr_flag > 0)
	{

		plan_source->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);

		if (sim.out_pk & MASK_T00)
		{
			snprintf(filename, SIZFILENAME, "%s%s%03d_T00.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real)numpts3d * (Real)numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00", a, sim.z_pk[pkcount]);
		}

		if (sim.out_pk & MASK_DELTA)
		{
			snprintf(filename, SIZFILENAME, "%s%s%03d_delta.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real)numpts3d * (Real)numpts3d * 2. * M_PI * M_PI * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)), filename, "power spectrum of delta", a, sim.z_pk[pkcount]);
		}

		if (cosmo.num_ncdm > 0 || sim.baryon_flag || sim.radiation_flag > 0 || sim.fluid_flag > 0)
		{
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
			if (sim.out_pk & MASK_T00)
			{
				snprintf(filename, SIZFILENAME, "%s%s%03d_T00cdm.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real)numpts3d * (Real)numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00 for cdm", a, sim.z_pk[pkcount]);
			}
			if (sim.out_pk & MASK_DELTA)
			{
				snprintf(filename, SIZFILENAME, "%s%s%03d_deltacdm.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real)numpts3d * (Real)numpts3d * 2. * M_PI * M_PI * (sim.baryon_flag ? (cosmo.Omega_cdm * cosmo.Omega_cdm) : ((cosmo.Omega_cdm + cosmo.Omega_b) * (cosmo.Omega_cdm + cosmo.Omega_b))), filename, "power spectrum of delta for cdm", a, sim.z_pk[pkcount]);
			}
			
		}
	}

	if (sim.out_pk & MASK_B)
	{
		extractPowerSpectrum(*BiFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		snprintf(filename, SIZFILENAME, "%s%s%03d_B.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, a * a * a * a * sim.numpts * sim.numpts * 2. * M_PI * M_PI, filename, "power spectrum of B", a, sim.z_pk[pkcount]);
	}

	free(kbin);
	free(power);
	free(kscatter);
	free(pscatter);
	free(occupation);
}

#endif
