#
#
FC=gfortran
FOPTS=-O
#
ALMA_OBJS=alma.o \
	  read_data_line.o \
	  to_uppercase.o \
	  config.o \
	  build_model.o \
	  write_log.o \
	  normalization.o \
	  direct_matrix.o \
	  inverse_matrix.o \
	  time_steps.o \
	  fluid_core_bc.o \
	  surface_bc.o \
	  love_numbers.o \
	  complex_rigidity.o \
	  lu.o \
	  salzer_weights.o \
          print_logo.o \
	  line_count.o
#
FMLIB_OBJS=fmsave.o fm.o fmzm90.o
#
#
all: ../alma.exe
#
clean:
	$(RM) ../alma.exe *.exe *.mod *.o
#
../alma.exe: $(FMLIB_OBJS) $(ALMA_OBJS) 
	$(FC) $(FOPTS) $(FMLIB_OBJS) $(ALMA_OBJS) -o ../alma.exe
#
#
%.o: %.f
	$(FC) $(FOPTS) -c $<
#
%.o: %.f90
	$(FC) $(FOPTS) -c $<
#
%.o: %.f95
	$(FC) $(FOPTS) -c $<
#
