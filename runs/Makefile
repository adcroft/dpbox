ALL_EXPTS = $(wildcard [A-Z]*x[0-9]*)
RUN = srun

.PRECIOUS: %/INPUT/MOM.res.nc
.PRECIOUS: ICs/*

all: $(foreach e,$(ALL_EXPTS),$e/RESTART/MOM.res.nc)
64: $(foreach e,$(wildcard [A-Z]*x64),$e/RESTART/MOM.res.nc)
128: $(foreach e,$(wildcard [A-Z]*x128),$e/RESTART/MOM.res.nc)
256: $(foreach e,$(wildcard [A-Z]*x256),$e/RESTART/MOM.res.nc)
512: $(foreach e,$(wildcard [A-Z]*x512),$e/RESTART/MOM.res.nc)

# Interpolate initial conditions
ICs/MOM.res.nc.64: spinup.64x64/RESTART/MOM.res.nc 
	cp -p spinup.64x64/RESTART/MOM.res.nc ICs/MOM.res.nc.64
ICs/MOM.res.nc.%: ICs/MOM.res.nc.64
	python3 interpIC.py ICs/MOM.res.nc.64 $* $*
	mv MOM.res.nc.$* ICs/
ICs/md5sum.txt: $(foreach r,64 128 256 512,ICs/MOM.res.nc.$(r))

# Setup INPUT/MOM.res.nc
%/INPUT/MOM.res.nc: | ICs/md5sum.txt
	mkdir -p $(@D)
	cp ICs/MOM.res.nc.$(lastword $(subst x, ,$(suffix $*))) $@

# Run the model
%512x512/RESTART/MOM.res.nc: PE=1024
%256x256/RESTART/MOM.res.nc: PE=512
%128x128/RESTART/MOM.res.nc: PE=128
%64x64/RESTART/MOM.res.nc: PE=32
%32x32/RESTART/MOM.res.nc: PE=8
test.32x32/RESTART/MOM.res.nc spinup.64x64/RESTART/MOM.res.nc:
	mkdir -p $(@D)
	cd $(dir $(@D)) && $(RUN) -n $(PE) ../../build/MOM6 2>&1 | tee std.out
%/RESTART/MOM.res.nc: | %/INPUT/MOM.res.nc
	mkdir -p $(@D)
	cd $* && $(RUN) -n $(PE) ../../build/MOM6 2>&1 | tee std.out
