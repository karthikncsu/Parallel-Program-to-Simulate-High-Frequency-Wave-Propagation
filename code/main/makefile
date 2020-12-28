FX :=mpif90
FXFLAGS := -Ofast

BLASLOC := -L/usr/local/intel/mkl/lib/intel64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential

FILESMAIN := variables.o main.o readinput.o checkinput.o phymanager.o intsort.o findpos.o intrmvdup.o dotproduct.o
 	
ADD1D := ./../oned
FILES1D := $(ADD1D)/waveeqn1d_fd_expl.o 

ADD2D := ./../twod
FILES2D := $(ADD2D)/readgrid2d.o $(ADD2D)/planestrain_fem_expl.o $(ADD2D)/qudasmbly.o $(ADD2D)/pspqudcnt.o $(ADD2D)/pspqudll.o $(ADD2D)/planestrain_fem_expl_parallel.o $(ADD2D)/divisionofgrid2d.o $(ADD2D)/mapping.o 

ADD3D := ./../threed
FILES3D := $(ADD3D)/readgrid3d.o $(ADD3D)/triarea3d.o $(ADD3D)/psphexcnt.o $(ADD3D)/psphexll.o $(ADD3D)/invmat3.o $(ADD3D)/Bsubform.o $(ADD3D)/FEM_expl.o $(ADD3D)/hexasmbly.o $(ADD3D)/divisionofgrid3d.o $(ADD3D)/mapping3d.o $(ADD3D)/FEM_expl_parallel.o $(ADD3D)/divisionofgrid3d_1dx.o $(ADD3D)/FEM_expl_parallel_1dx.o $(ADD3D)/divisionofgrid3d_1dx_parallel.o $(ADD3D)/divisionofgrid3d_2d_parallel.o $(ADD3D)/FEM_expl_parallel_2d.o

TOTADD := -I $(ADD1D) $(ADD2D) $(ADD3D)

MODS :=variables.mod

OBJECTS := $(FILESMAIN) $(FILES1D) $(FILES2D) $(FILES3D)

all: project

project: $(OBJECTS)
	 $(FX) $(FXFLAGS) $(OBJECTS) -o project $(BLASLOC)

#variables.mod: variables.o variables.f
#		$(FX) $(FXFLAGS) -c variables.f

$(FILESMAIN): %.o: %.f
		$(FX) $(FXFLAGS) -c $< -o $@ $(TOTADD)

$(FILES1D): %.o: %.f
		$(FX) $(FXFLAGS) -c $< -o $@ -I.

$(FILES2D): %.o: %.f
		$(FX) $(FXFLAGS) -c $< -o $@ -I.

$(FILES3D): %.o: %.f
		$(FX) $(FXFLAGS) -c $< -o $@ -I.

clean:
	rm -rf $(OBJECTS) project *.mod

cleano:
	rm -rf $(OBJECTS)
