include ../Makefile.inc

FPPFLAGS  = -fpp -mkl=parallel#-DMESSAGE_DETAIL -DDEBUG_RESTRICTION -DDEBUG_GRID #-DINPLACE #-DDEBUG1 #-DDEBUG_PROLONGATION #-DUSE_MKL

OBJDIR = ./obj
OBJS  = $(OBJDIR)/para_range.o \
	$(OBJDIR)/mpi_topology.o \
	$(OBJDIR)/geometry.o \
	$(OBJDIR)/matrix.o \
	$(OBJDIR)/poisson_matrix_operator.o \
	$(OBJDIR)/cg_poisson_matrix.o \
	$(OBJDIR)/rbgs_poisson_matrix.o \
	$(OBJDIR)/multigrid_common.o \
	$(OBJDIR)/multigrid.o \
	$(OBJDIR)/main_peps.o

LIBS	= -liomp5 -lpthread -lm

EXE	= peps_mg.ex
all	: $(OBJS)
	$(FC) -o $(EXE) $(FFLAGS) $(FPPFLAGS) $(OBJS) $(LIBS)
	mv $(EXE) ../run

clean	:
	rm -rf ../run/$(EXE) $(OBJDIR)


$(OBJDIR)/%.o: %.f90
	$(FC) -c $(FLAGS) $(FPPFLAGS) -module $(OBJDIR) $? -o $@ 

$(OBJDIR)/%.o: $(OBJDIR)/%.mod


