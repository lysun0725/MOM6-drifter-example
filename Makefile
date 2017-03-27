SHELL = bash -f
site = linux
platform = gnu

all: test_suite_gfortran

compile_libs:
	mkdir -p build/shared.$(platform).repro
	(cd build/shared.$(platform).repro/; \
	rm -f path_names; ../../src/mkmf/bin/list_paths ../../src/FMS; \
	../../src/mkmf/bin/mkmf -t ../../src/mkmf/templates/$(site)-$(platform).mk -p libfms.a -c "-Duse_libMPI -Duse_netCDF -DSPMD -DLAND_BND_TRACERS" path_names; \
	source ../../build/env/$(site)-$(platform); make  NETCDF=4 REPRO=1 libfms.a)

	mkdir -p build/MOM6.$(platform).repro
	(cd build/MOM6.$(platform).repro/; \
	rm -f path_names; ../../src/mkmf/bin/list_paths ../../src/MOM6/{config_src/dynamic,config_src/solo_driver/coupler_types*,src/*,src/*/*,pkg/*/src/*}  ;\
	../../src/mkmf/bin/mkmf -t ../../src/mkmf/templates/$(site)-$(platform).mk -o '-I../shared.$(platform).repro' -p libMOM6.a -l '-L../shared.$(platform).repro -lfms' -c "-Duse_libMPI -Duse_netCDF -DSPMD -DLAND_BND_TRACERS -Duse_AM3_physics" path_names; \
	source ../../build/env/$(site)-$(platform); make  NETCDF=4 REPRO=1 -j 8 libMOM6.a)

compile_MOM6_solo: compile_libs
	mkdir -p build/MOM6_solo.$(platform).repro
	(cd build/MOM6_solo.$(platform).repro/; \
	rm -f path_names; ../../src/mkmf/bin/list_paths ../../src/MOM6/config_src/solo_driver/*  ;\
	../../src/mkmf/bin/mkmf -t ../../src/mkmf/templates/$(site)-$(platform).mk -o  '-I../../src/MOM6/src/framework -I../../src/MOM6/config_src/dynamic -I../MOM6.$(platform).repro -I../shared.$(platform).repro' -l '-L../MOM6.$(platform).repro -lMOM6 -L../shared.$(platform).repro -lfms' -p MOM6 -c "-Duse_libMPI -Duse_netCDF -DSPMD -DLAND_BND_TRACERS -Duse_AM3_physics" path_names; \
	source ../../build/env/$(site)-$(platform); make  NETCDF=4 REPRO=1 -j 8 MOM6)

compile_drifters: compile_libs
	mkdir -p build/drifters.$(platform).repro
	(cd build/drifters.$(platform).repro/; \
	rm -f path_names; ../../src/mkmf/bin/list_paths ../../src/drifters/{particle_driver.F90}  ;\
	../../src/mkmf/bin/mkmf -t ../../src/mkmf/templates/$(site)-$(platform).mk -o  '-I../../src/MOM6/src/framework -I../../src/MOM6/config_src/dynamic -I../MOM6.$(platform).repro -I../shared.$(platform).repro' -l '-L../MOM6.$(platform).repro -lMOM6 -L../shared.$(platform).repro -lfms' -p drifters -c "-Duse_libMPI -Duse_netCDF -DSPMD -DLAND_BND_TRACERS -Duse_AM3_physics" path_names; \
	source ../../build/env/$(site)-$(platform); make  NETCDF=4 REPRO=1 -j 8 drifters)
	mkdir -p build/drifters.$(platform).debug
	(cd build/drifters.$(platform).debug/; \
	rm -f path_names; ../../src/mkmf/bin/list_paths ../../src/drifters/{particle_driver.F90}  ;\
	../../src/mkmf/bin/mkmf -t ../../src/mkmf/templates/$(site)-$(platform).mk -o  '-I../../src/MOM6/src/framework -I../../src/MOM6/config_src/dynamic -I../MOM6.$(platform).repro -I../shared.$(platform).repro' -l '-L../MOM6.$(platform).repro -lMOM6 -L../shared.$(platform).repro -lfms' -p drifters -c "-Duse_libMPI -Duse_netCDF -DSPMD -DLAND_BND_TRACERS -Duse_AM3_physics" path_names; \
	source ../../build/env/$(site)-$(platform); make  NETCDF=4 DEBUG=1 -j 8 drifters)


test_suite_gfortran: compile_MOM6_solo
	(cd ocean_only/double_gyre;mpirun -n 1 ../../build/MOM6_solo.$(platform).repro/MOM6 ; echo "double_gyre 1pe :" | tee  ../../.results_$(platform);diff -s ocean.stats ocean.stats.$(platform) | tee -a ../../.results_$(platform))
	(cd ocean_only/double_gyre;mpirun -n 1 ../../build/drifters.$(platform).debug/drifters ; echo "drifters :")
