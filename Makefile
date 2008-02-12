#	Top-level Makefile for LITTLE_R system - including objective analysis
#	and supporting ancillary plotting code

#	Macros, these should be generic for all machines

.IGNORE:

AR		=	ar ru
CD		=	cd
LN		=	ln -s
MAKE		=	make -i -f Makefile
RM		=	/bin/rm -f
RM_LIST		=	*.o *.M core *.kmo *.mod *_out_* *.cgm

#	If you don't have NCAR Graphics, then you can't look at the plots that could have
#	been produced by the plotting packages.  Do not build the two extra programs in this case.

I_HAVE_NCARG		=	all
I_DONT_HAVE_NCARG	=	little_r
PROGS			=	$(I_HAVE_NCARG)
#PROGS			=	$(I_DONT_HAVE_NCARG)

#	Targets for supported architectures

default:
	uname -a > .tmpfile
	@grep CRAY .tmpfile ; \
	if [ $$? = 0 ] ; then echo "Compiling for Cray"							; \
		echo "AR		=	$(AR)"				>  macros_little_r	; \
		echo "RM		=	$(RM)"				>> macros_little_r	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_little_r	; \
		echo "CD		=	$(CD)"				>> macros_little_r	; \
		echo "LN		=	$(LN)"				>> macros_little_r	; \
		echo "MAKE		=	$(MAKE)"			>> macros_little_r	; \
		echo "SHELL		=	/bin/sh"			>> macros_little_r	; \
		echo "TOUCH		=	touch"				>> macros_little_r	; \
		echo "CPP		=	cpp"				>> macros_little_r	; \
		echo "CPPFLAGS		=	-I. -C -P -DCRAY"		>> macros_little_r	; \
		echo "FC		=	f90"				>> macros_little_r	; \
		echo "FCFLAGS		=	-f free -x omp"			>> macros_little_r	; \
		echo "LDFLAGS		=	-Ca"				>> macros_little_r	; \
		echo "CCFLAGS		=	-DCRAY -I."			>> macros_little_r	; \
		echo "LOCAL_LIBRARIES	=	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm" >> macros_little_r	; \
		( $(CD) src ; $(MAKE) $(PROGS) )								; \
	else grep OSF .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for Compaq"						; \
		echo "AR		=	$(AR)"				>  macros_little_r	; \
		echo "RM		=	$(RM)"				>> macros_little_r	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_little_r	; \
		echo "CD		=	$(CD)"				>> macros_little_r	; \
		echo "LN		=	$(LN)"				>> macros_little_r	; \
		echo "MAKE		=	$(MAKE)"			>> macros_little_r	; \
		echo "SHELL		=	/bin/sh"			>> macros_little_r	; \
		echo "TOUCH		=	touch"				>> macros_little_r	; \
		echo "CPP		=	/usr/bin/cpp"			>> macros_little_r	; \
		echo "CPPFLAGS		=	-I. -C -P -DDEC"		>> macros_little_r	; \
		echo "FC		=	f90"				>> macros_little_r	; \
		echo "FCFLAGS		=	-automatic -convert big_endian -fast -fpe -free -pipeline -O4 -std " >> macros_little_r	; \
		echo "LDFLAGS		=	-fast -O4"			>> macros_little_r	; \
		echo "CCFLAGS		=	-DDEC -I."			>> macros_little_r	; \
		echo "LOCAL_LIBRARIES	=	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm"	>> macros_little_r	; \
		( $(CD) src ; $(MAKE) $(PROGS) )							; \
	else grep HP .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for HP"							; \
		echo "AR		=	$(AR)"				>  macros_little_r	; \
		echo "RM		=	$(RM)"				>> macros_little_r	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_little_r	; \
		echo "CD		=	$(CD)"				>> macros_little_r	; \
		echo "LN		=	$(LN)"				>> macros_little_r	; \
		echo "MAKE		=	$(MAKE)"			>> macros_little_r	; \
		echo "SHELL		=	/bin/sh"			>> macros_little_r	; \
		echo "TOUCH		=	touch"				>> macros_little_r	; \
		echo "CPP		=	/opt/langtools/lbin/cpp"	>> macros_little_r	; \
		echo "CPPFLAGS		=	-I. -C -P -DHP"			>> macros_little_r	; \
		echo "FC		=	f90"				>> macros_little_r	; \
		echo "FCFLAGS		=	+langlvl=90 +source=free"	>> macros_little_r	; \
		echo "LDFLAGS		=	" 				>> macros_little_r	; \
		echo "CCFLAGS		=	-DHP -I."			>> macros_little_r	; \
		echo "LOCAL_LIBRARIES	=	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm"	>> macros_little_r	; \
		( $(CD) src ; $(MAKE) $(PROGS) )							; \
	else grep AIX .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for IBM"							; \
		echo "AR		=	$(AR)"				>  macros_little_r	; \
		echo "RM		=	$(RM)"				>> macros_little_r	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_little_r	; \
		echo "CD		=	$(CD)"				>> macros_little_r	; \
		echo "LN		=	$(LN)"				>> macros_little_r	; \
		echo "MAKE		=	$(MAKE)"			>> macros_little_r	; \
		echo "SHELL		=	/bin/sh"			>> macros_little_r	; \
		echo "TOUCH		=	touch"				>> macros_little_r	; \
		echo "CPP		=	/usr/lib/cpp"			>> macros_little_r	; \
		echo "CPPFLAGS		=	-I. -C -P -DIBM"		>> macros_little_r	; \
		echo "FC		=	xlf90"				>> macros_little_r	; \
		echo "FCFLAGS		=	-qfree=f90 -qlanglvl=90pure -q64 -qarch=auto -qnosave -qmaxmem=-1 -qspillsize=20000 -Q"	>> macros_little_r	; \
		echo "LDFLAGS		=	-q64 " 			>> macros_little_r	; \
		echo "CCFLAGS		=	-DIBM -I."			>> macros_little_r	; \
		echo "LOCAL_LIBRARIES	=	-L/usr/local/lib64/r4i4 -lncarg -lncarg_gks -lncarg_c -lX11 -lm" >> macros_little_r	; \
		( $(CD) src ; $(MAKE) $(PROGS) )							; \
	else grep Darwin .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for MAC"							; \
		echo "AR		=	$(AR)"				>  macros_little_r	; \
		echo "RM		=	$(RM)"				>> macros_little_r	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_little_r	; \
		echo "CD		=	$(CD)"				>> macros_little_r	; \
		echo "LN		=	$(LN)"				>> macros_little_r	; \
		echo "MAKE		=	$(MAKE)"			>> macros_little_r	; \
		echo "SHELL		=	/bin/sh"			>> macros_little_r	; \
		echo "TOUCH		=	touch"				>> macros_little_r	; \
		echo "CPP		=	/usr/bin/cpp"			>> macros_little_r	; \
		echo "CPPFLAGS		=	-I. -C -P -DIBM -xassembler-with-cpp" >> macros_little_r	; \
		echo "FC		=	xlf90"				>> macros_little_r	; \
		echo "FCFLAGS		=	-qfree=f90 -qarch=auto "	>> macros_little_r	; \
		echo "LDFLAGS		=	-Wl,-stack_size,10000000,-stack_addr,0xc0000000" >> macros_little_r	; \
		echo "CCFLAGS		=	-DNOUNDERSCORE -DIBM -I."	>> macros_little_r	; \
		echo "LOCAL_LIBRARIES	=	-L/usr/local/ncarg/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm" >> macros_little_r	; \
		( $(CD) src ; $(MAKE) $(PROGS) )							; \
	else grep Linux .tmpfile									; \
	if [ $$? = 0 ] ; then echo "Compiling for Linux"						; \
		echo "AR		=	$(AR)"				>  macros_little_r	; \
		echo "RM		=	$(RM)"				>> macros_little_r	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_little_r	; \
		echo "CD		=	$(CD)"				>> macros_little_r	; \
		echo "LN		=	$(LN)"				>> macros_little_r	; \
		echo "MAKE		=	$(MAKE)"			>> macros_little_r	; \
		echo "SHELL		=	/bin/sh"			>> macros_little_r	; \
		echo "TOUCH		=	touch"				>> macros_little_r	; \
		echo "CPP		=	/lib/cpp"			>> macros_little_r	; \
		echo "CPPFLAGS		=	-I. -C -P -DDEC -traditional "	>> macros_little_r	; \
		echo "FC		=	pgf90"				>> macros_little_r	; \
		echo "FCFLAGS		=	-Mfreeform -pc 32 -byteswapio "	>> macros_little_r	; \
		echo "LDFLAGS		=	-L/usr/local/netcdf-3.6.0-p1-pgi/lib -lnetcdf -lm -I/usr/local/netcdf-3.6.0-p1-pgi/include" 				>> macros_little_r	; \
		echo "CCFLAGS		=	-DDEC -I. "			>> macros_little_r	; \
		echo "LOCAL_LIBRARIES	=	-L$(NCARG_ROOT)/lib -L/usr/X11R6/lib -lncarg -lncarg_gks -lncarg_c -lX11 -L$(PGI)/linux86/lib -L/usr/lib -lf2c" >> macros_little_r	; \
		( $(CD) src ; $(MAKE) $(PROGS) ) 							; \
	else grep IRIX .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for SGI"							; \
		echo "AR		=	$(AR)"				>  macros_little_r	; \
		echo "RM		=	$(RM)"				>> macros_little_r	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_little_r	; \
		echo "CD		=	$(CD)"				>> macros_little_r	; \
		echo "LN		=	$(LN)"				>> macros_little_r	; \
		echo "MAKE		=	$(MAKE)"			>> macros_little_r	; \
		echo "SHELL		=	/bin/sh"			>> macros_little_r	; \
		echo "TOUCH		=	touch"				>> macros_little_r	; \
		echo "CPP		=	/lib/cpp"			>> macros_little_r	; \
		echo "CPPFLAGS		=	-I. -C -P -DSGI"		>> macros_little_r	; \
		echo "FC		=	f90"				>> macros_little_r	; \
		echo "FCFLAGS		=	-freeform -n32 -O2 -I."		>> macros_little_r	; \
		echo "LDFLAGS		=	-n32 -O2 -lfastm"		>> macros_little_r	; \
		echo "CCFLAGS		=	-DSGI -I. -n32"			>> macros_little_r	; \
		echo "LOCAL_LIBRARIES	=	-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm"	>> macros_little_r	; \
		( $(CD) src ; $(MAKE) $(PROGS) )							; \
	else grep SUN .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for SUN"							; \
		echo "AR		=	$(AR)"				>  macros_little_r	; \
		echo "RM		=	$(RM)"				>> macros_little_r	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_little_r	; \
		echo "CD		=	$(CD)"				>> macros_little_r	; \
		echo "LN		=	$(LN)"				>> macros_little_r	; \
		echo "MAKE		=	$(MAKE)"			>> macros_little_r	; \
		echo "SHELL		=	/bin/sh"			>> macros_little_r	; \
		echo "TOUCH		=	touch"				>> macros_little_r	; \
		echo "CPP		=	/usr/ccs/lib/cpp"		>> macros_little_r	; \
		echo "CPPFLAGS		=	-I. -C -P -DSUN"		>> macros_little_r	; \
		echo "FC		=	f90"				>> macros_little_r	; \
		echo "FCFLAGS		=	-ansi -free"			>> macros_little_r	; \
		echo "LDFLAGS		=	" 				>> macros_little_r	; \
		echo "CCFLAGS		=	-DSUN -I."			>> macros_little_r	; \
		echo "LOCAL_LIBRARIES	=	-L$(NCARG_ROOT)/lib -L/usr/openwin/lib -L/usr/dt/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm" >> macros_little_r	; \
		( $(CD) src ; $(MAKE) $(PROGS) ) 							; \
	else grep UXP .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for Fujitsu"						; \
		echo "AR		=	$(AR)"				>  macros_little_r	; \
		echo "RM		=	$(RM)"				>> macros_little_r	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_little_r	; \
		echo "CD		=	$(CD)"				>> macros_little_r	; \
		echo "LN		=	$(LN)"				>> macros_little_r	; \
		echo "MAKE		=	$(MAKE)"			>> macros_little_r	; \
		echo "SHELL		=	/bin/sh"			>> macros_little_r	; \
		echo "TOUCH		=	touch"				>> macros_little_r	; \
		echo "CPP		=	/lib/cpp"			>> macros_little_r	; \
		echo "CPPFLAGS		=	-I. -C -P -DVPP -DBIT32"	>> macros_little_r	; \
		echo "FC		=	frt"				>> macros_little_r	; \
		echo "FCFLAGS		=	-Free -X9 -Am -sc -Kfast -Kfreealloc -Karraystack3"	>> macros_little_r	; \
		echo "LDFLAGS		=	-J" 				>> macros_little_r	; \
		echo "CCFLAGS		=	-DVPP -DBIT32 -I."		>> macros_little_r	; \
		echo "LOCAL_LIBRARIES	=	"				>> macros_little_r	; \
		( $(CD) src ; $(MAKE) $(PROGS) ) 							; \
	else echo "Do not know how to compile for the `cat .tmpfile` machine." 				; \
	fi ; \
	fi ; \
	fi ; \
	fi ; \
	fi ; \
	fi ; \
	fi ; \
	fi ; \
	fi ;

	$(RM) little_r ; $(LN) src/little_r .
	if [ $(PROGS) = $(I_HAVE_NCARG) ] ; then echo "Building plot programs"	; \
		$(RM) plot_level 	; 	$(LN) src/plot_level . 		; \
		$(RM) plot_soundings 	; 	$(LN) src/plot_soundings . 	; \
	fi

intel:
	echo "Compiling for Linux with INTEL compiler"						; \
		echo "AR		=	$(AR)"				>  macros_little_r	; \
		echo "RM		=	$(RM)"				>> macros_little_r	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_little_r	; \
		echo "CD		=	$(CD)"				>> macros_little_r	; \
		echo "LN		=	$(LN)"				>> macros_little_r	; \
		echo "MAKE		=	$(MAKE)"			>> macros_little_r	; \
		echo "SHELL		=	/bin/sh"			>> macros_little_r	; \
		echo "TOUCH		=	touch"				>> macros_little_r	; \
		echo "CPP		=	/lib/cpp"			>> macros_little_r	; \
		echo "CPPFLAGS		=	-I. -C -P -DDEC -traditional"	>> macros_little_r	; \
		echo "FC		=	ifort"				>> macros_little_r	; \
		echo "FCFLAGS		=	-FR -pc 32 -convert big_endian"	>> macros_little_r	; \
		echo "LDFLAGS		=	-i_dynamic" 				>> macros_little_r	; \
		echo "CCFLAGS		=	-DDEC -I."			>> macros_little_r	; \
		echo "LOCAL_LIBRARIES	=	-L$(NCARG_ROOT)/lib -L/usr/X11R6/lib -lncarg -lncarg_gks -lncarg_c -lX11 -L/usr/lib/gcc-lib/i386-redhat-linux/3.3.2  -lg2c" >> macros_little_r	; \
		( $(CD) src ; $(MAKE) $(PROGS) ) 							; \
	$(RM) little_r ; $(LN) src/little_r .
	if [ $(PROGS) = $(I_HAVE_NCARG) ] ; then echo "Building plot programs"	; \
		$(RM) plot_level 	; 	$(LN) src/plot_level . 		; \
		$(RM) plot_soundings 	; 	$(LN) src/plot_soundings . 	; \
	fi

clean:
	( $(CD) src ; $(MAKE) clean "CD = $(CD)" "RM = $(RM)" "RM_LIST = $(RM_LIST)" )
	$(RM) $(RM_LIST) .tmpfile
