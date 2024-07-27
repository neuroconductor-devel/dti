
#include "Fibertracking.h"
#include "Converter.h"
#include <iostream>
#include <string>
#include <R.h>
#include <Rinternals.h>

extern "C"{
	SEXP interface_tracking(SEXP data_dir_coords, SEXP data_FA_values, 
                                SEXP data_mask,
			        SEXP dim_x, SEXP dim_y, SEXP dim_z,
			        SEXP roi_x_s, SEXP roi_x_e,
			        SEXP roi_y_s, SEXP roi_y_e,
			        SEXP roi_z_s, SEXP roi_z_e,
			        SEXP voxelext_x, SEXP voxelext_y, SEXP voxelext_z,
			        SEXP min_anisotropy,
			        SEXP max_angle
			       )
	{
		SEXP retVec;
		int ii, length, nProtected = 0;
		
		Converter converter(REAL(data_dir_coords), REAL(data_FA_values), INTEGER(data_mask), *INTEGER(dim_x), *INTEGER(dim_y), *INTEGER(dim_z));
		
		int marked[] = {*INTEGER(roi_x_s), *INTEGER(roi_x_e), *INTEGER(roi_y_s), *INTEGER(roi_y_e), *INTEGER(roi_z_s), *INTEGER(roi_z_e)};
		
		Fibertracking *tester = new Fibertracking(converter.getVoxels(), *INTEGER(dim_x), *INTEGER(dim_y), *INTEGER(dim_z), *REAL(voxelext_x), *REAL(voxelext_y), *REAL(voxelext_z), *REAL(min_anisotropy), *REAL(max_angle));
		tester->findMarkedFibers(marked);
		
		length = tester->getLength();
		
		double *vals = tester->convertToDouble();
		
//		Rprintf("length: %d\n", length);
	
		delete tester;
			
		PROTECT(retVec = Rf_allocVector(REALSXP, length));
		++nProtected;
		
		for (ii = 0; ii < length; ++ii)
		{
			REAL(retVec)[ii] = vals[ii];
		}
		
		UNPROTECT(nProtected);
		
		delete vals;
		
		return retVec;
	}
} // extern C end

extern "C"{

	SEXP interface_tracking_mixtensor(
					SEXP data_dir_coords, SEXP data_order, 
					SEXP data_FA_values, SEXP data_mask,
					SEXP data_mix, 
					SEXP maxorder, 
			        SEXP dim_x, SEXP dim_y, SEXP dim_z,
			        SEXP roi_x_s, SEXP roi_x_e,
			        SEXP roi_y_s, SEXP roi_y_e,
			        SEXP roi_z_s, SEXP roi_z_e,
			        SEXP voxelext_x, SEXP voxelext_y, SEXP voxelext_z,
			        SEXP min_anisotropy,
			        SEXP max_angle
			       )
	{		
		SEXP retVec;

//		Rprintf("retVec allociert\n");

		int ii, length, nProtected = 0;
		
		// data_dir_num in Konstruktor
		Converter *converter = new Converter(REAL(data_dir_coords), REAL(data_FA_values), INTEGER(data_mask), INTEGER(data_order), *INTEGER(maxorder), *INTEGER(dim_x), *INTEGER(dim_y), *INTEGER(dim_z));
		
//		Rprintf("R-Daten konvertiert\n");

		int marked[] = {*INTEGER(roi_x_s), *INTEGER(roi_x_e), *INTEGER(roi_y_s), *INTEGER(roi_y_e), *INTEGER(roi_z_s), *INTEGER(roi_z_e)};
		
		Fibertracking *tracks = new Fibertracking(converter->getVoxels(), *INTEGER(dim_x), *INTEGER(dim_y), *INTEGER(dim_z), *REAL(voxelext_x), *REAL(voxelext_y), *REAL(voxelext_z), *REAL(min_anisotropy), *REAL(max_angle));

//		Rprintf("Fibertracking-Objekt angelegt\n");

		tracks->findMarkedFibers(marked);
		
//		Rprintf("Fibertracking durchgefuehrt\n");

		length = tracks->getLength();
		
		double *vals = tracks->convertToDouble();
		
//		Rprintf("verkettete Liste konvertiert mit laenge %d\n", length);

//		Rprintf("length: %d\n", length);
	
		delete tracks;
			
//		Rprintf("delete tracks ausgefuehrt\n");

		PROTECT(retVec = Rf_allocVector(REALSXP, length));
		++nProtected;

//		Rprintf("Protect ausgefuehrt\n");

		for (ii = 0; ii < length; ++ii)
		{
			REAL(retVec)[ii] = vals[ii];
		}
		
//		Rprintf("umgeschrieben\n");

		UNPROTECT(nProtected);
		
//		Rprintf("UNProtect ausgefuehrt\n");

		delete vals;

//		Rprintf("delete vals ausgefuehrt\n");
		
		return retVec;
	}
} // extern C end
