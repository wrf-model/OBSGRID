#	Lower level Makefile for OBSGRID and plotting packages

include ../configure.oa
#include ../macros_little_r

.IGNORE:

#	obsgrid

SRC 		= 	driver.F error_handler.F linpak.F main.F \
			proc_final_analysis.F proc_first_guess.F proc_get_info_header.F  \
			proc_header.F proc_namelist.F proc_oa.F \
			proc_ob_access.F proc_obs_sort.F proc_qc.F sample.F yx2xy.F

OBJS		=	$(SRC:.F=.o)

MODULES 	=	date_pack_module.o error_module.o final_analysis_module.o  \
			first_guess_module.o \
			input_data_module.o map_utils_module.o map_utils_helper_module.o \
			namelist_module.o oa_module.o \
			ob_access_module.o obs_sort_module.o printout_module.o \
			qc0_module.o qc1_module.o qc2_module.o qc3_module.o \
			small_header_data_module.o tinterp_module.o

INCS		=	big_header.inc constants.inc error.inc first_guess.inc first_guess_size.inc \
			map.inc missing.inc namelist.inc proc_get_info_header.inc

INTS		=	error.int proc_final_analysis.int proc_first_guess.int \
			proc_get_info_header.int proc_header.int proc_oa.int proc_ob_access.int \
			proc_obs_sort.int proc_qc.int sample.int

NML		=	namelist.nml

COMMON		=	header.common namelist.common

#	plot_level

SRC_OBS		=	plot_level.F

OBJS_OBS	=	$(SRC_OBS:.F=.o)

MODULES_OBS 	=	module_header.o module_map_stuff.o module_ncarg.o module_read_station.o

#	plot_soundings

SRC_SND		=	plot_soundings.F

OBJS_SND	=	$(SRC_SND:.F=.o)

MODULES_SND 	=	module_mapinfo.o module_report.o module_skewt.o


all:			obsgrid.exe plot_level plot_soundings \
			$(OBJS) $(MODULES) $(INTS) $(INCS) $(NML) $(COMMON) \
			$(OBJS_OBS) $(MODULES_OBS) \
			$(OBJS_SND) $(MODULES_SND) \
			../configure.oa ./Makefile ../Makefile

obsgrid.exe: 		$(OBJS) $(MODULES) $(INTS) $(INCS) $(NML) $(COMMON) \
			../configure.oa ./Makefile ../Makefile 
			$(FC) -o $@ $(FCFLAGS) $(OBJS) $(MODULES) ${NETCDF_LIBS} ${NETCDF_INC} $(LDFLAGS)

plot_level: 		$(OBJS_OBS) $(MODULES_OBS) \
			../configure.oa ./Makefile ../Makefile 
			$(FC) -o $@ $(FCFLAGS) $(LDFLAGS) $(OBJS_OBS) $(MODULES_OBS) $(NCARG_LIBS)

plot_soundings:		$(OBJS_SND) $(MODULES_SND) \
			../configure.oa ./Makefile ../Makefile 
			$(FC) -o $@ $(FCFLAGS) $(LDFLAGS) $(OBJS_SND) $(MODULES_SND) $(NCARG_LIBS)

#	obsgrid files

date_pack_module.o:	../configure.oa ./Makefile ../Makefile

driver.o:		error.inc proc_get_info_header.inc first_guess_size.inc first_guess.inc \
			map.inc error.int map_utils_module.o map_utils_helper_module.o \
			proc_get_info_header.int proc_first_guess.int proc_qc.int \
			proc_obs_sort.int proc_final_analysis.int proc_oa.int \
			date_pack_module.o input_data_module.o namelist_module.o obs_sort_module.o \
			tinterp_module.o \
			../configure.oa ./Makefile ../Makefile

error_handler.o:	error_module.o \
			../configure.oa ./Makefile ../Makefile

error_module.o:		../configure.oa ./Makefile ../Makefile

final_analysis_module.o:	error.inc error.int big_header.inc sample.int \
			first_guess_size.inc first_guess.inc \
			first_guess_module.o input_data_module.o \
			date_pack_module.o \
			../configure.oa ./Makefile ../Makefile

first_guess_module.o:	error.inc error.int first_guess_size.inc first_guess.inc big_header.inc \
			input_data_module.o date_pack_module.o small_header_data_module.o \
			oa_module.o \
			../configure.oa ./Makefile ../Makefile


input_data_module.o:	small_header_data_module.o \
			../configure.oa ./Makefile ../Makefile

linpak.o:		../configure.oa ./Makefile ../Makefile

main.o:                 error.inc big_header.inc proc_header.int error.int \
			proc_get_info_header.int \
			namelist_module.o date_pack_module.o input_data_module.o \
			../configure.oa ./Makefile ../Makefile

map_utils_module.o:	../configure.oa ./Makefile ../Makefile

map_utils_helper_module.o:	map_utils_module.o \
			../configure.oa ./Makefile ../Makefile 

make_sfc_fdda.o:	diags_module.o \
			../configure.oa ./Makefile ../Makefile

namelist_module.o:	error.inc error.int namelist.inc namelist.common \
			../configure.oa ./Makefile ../Makefile

oa_module.o:		error.int error.inc first_guess_size.inc first_guess.inc \
			../configure.oa ./Makefile ../Makefile

ob_access_module.o:	error.int error.inc missing.inc \
			obs_sort_module.o date_pack_module.o  \
			../configure.oa ./Makefile ../Makefile

obs_sort_module.o:	missing.inc error.inc error.int proc_get_info_header.inc map.inc \
			constants.inc \
			date_pack_module.o map_utils_module.o map_utils_helper_module.o \
			../configure.oa ./Makefile ../Makefile

printout_module.o:	../configure.oa ./Makefile ../Makefile

proc_final_analysis.o:	first_guess_size.inc first_guess.inc proc_get_info_header.inc \
			proc_get_info_header.int sample.int header.common \
			final_analysis_module.o obs_sort_module.o \
			../configure.oa ./Makefile ../Makefile

proc_first_guess.o:	big_header.inc sample.int \
			input_data_module.o first_guess_module.o  \
			../configure.oa ./Makefile ../Makefile

proc_get_info_header.o:	header.common \
			../configure.oa ./Makefile ../Makefile

proc_header.o:		header.common \
			namelist_module.o \
			../configure.oa ./Makefile ../Makefile

proc_namelist.o:	error.inc namelist.inc namelist.nml namelist.common error.int \
			namelist_module.o \
			../configure.oa ./Makefile ../Makefile

proc_oa.o:		error.int error.inc first_guess_size.inc first_guess.inc \
			proc_ob_access.int \
			oa_module.o obs_sort_module.o \
			../configure.oa ./Makefile ../Makefile

proc_ob_access.o:	error.inc error.int \
			obs_sort_module.o ob_access_module.o  map_utils_module.o \
			map_utils_helper_module.o \
			../configure.oa ./Makefile ../Makefile

proc_obs_sort.o:	error.inc error.int \
			obs_sort_module.o \
			../configure.oa ./Makefile ../Makefile

proc_qc.o:		first_guess_size.inc first_guess.inc proc_ob_access.int \
			qc0_module.o qc1_module.o qc2_module.o qc3_module.o \
			obs_sort_module.o \
			../configure.oa ./Makefile ../Makefile

qc0_module.o:		obs_sort_module.o  \
			../configure.oa ./Makefile ../Makefile

qc1_module.o:		qc0_module.o \
			../configure.oa ./Makefile ../Makefile

qc2_module.o:		qc0_module.o \
			../configure.oa ./Makefile ../Makefile

qc3_module.o:		obs_sort_module.o \
			../configure.oa ./Makefile ../Makefile

sample.o:		printout_module.o \
			../configure.oa ./Makefile ../Makefile

small_header_data_module.o:	../configure.oa ./Makefile ../Makefile

tinterp_module.o:	small_header_data_module.o input_data_module.o \
			../configure.oa ./Makefile ../Makefile

temp_interp.o:		date_pack_module.o \
			../configure.oa ./Makefile ../Makefile

yx2xy.o:		../configure.oa ./Makefile ../Makefile

#	plot_level files

plot_level.o:		$(MODULES_OBS) \
			 ../configure.oa ./Makefile ../Makefile

module_header.o:	../configure.oa ./Makefile ../Makefile

module_map_stuff.o:	../configure.oa ./Makefile ../Makefile

module_ncarg.o:		../configure.oa ./Makefile ../Makefile

module_read_station.o:	../configure.oa ./Makefile ../Makefile

#	plot_soundings files

plot_soundings.o:	$(MODULES_SND) \
			 ../configure.oa ./Makefile ../Makefile

module_mapinfo.o:	../configure.oa ./Makefile ../Makefile

module_report.o:	../configure.oa ./Makefile ../Makefile

module_skewt.o:		../configure.oa ./Makefile ../Makefile

clean:
			$(RM) *.o *.f core obsgrid plot_level plot_soundings *.kmo ctmi* *.mod *_out_* output *.out
