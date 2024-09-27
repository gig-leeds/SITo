include Inc/SDR.top

gradients_mod.deps=$O/lsqrDataModule.o
warp2d_mod.deps=$O/sinc_mod.o $O/gradients2d_mod.o $O/lsqrModule.o $O/lsqrDataModule.o
lsqrModule.deps=$O/lsqrDataModule.o $O/lsqrblasInterface.o
acfd_born_mod_mpi.deps=$O/sdr_tmp_io_mpi.o $O/processing_mod.o

acfd_fwi_pre-idwc_mpi.deps=$O/acfd_born_mod_mpi.o $O/lsqrDataModule.o $O/lsqrModule.o $O/lsqrblas.o $O/lsqrblasInterface.o $O/processing_mod.o $O/warp2d_mod.o $O/sinc_mod.o $O/gradients2d_mod.o $O/sdr_tmp_io_mpi.o

bins: $B/acfd_fwi_pre-idwc_mpi.x

flags=-O3 -inline-level=2 -fpp -qopenmp
#flags=-g -check all -fpe0 -warn -traceback -debug extended -qopenmp
$O/%.o: $S/%.f90
	make ${$*.deps} none
	mpifort -o $@ $< -c ${flags} -I${RSFROOT}/include -I${FFTWROOT}/include -module $M -L${RSFROOT}/lib -L${FFTWROOT}/lib -lrsff90 -lrsf -lfftw3f

$B/%.x: $O/%.o
	make ${$*.deps} none
	mpifort -o $@ $O/$*.o ${flags} ${$*.deps} ${none} -L${RSFROOT}/lib -L${FFTWROOT}/lib -L${OPENMPROOT}/lib -lrsff90 -lrsf -lfftw3f

run: $B/acfd_fwi_pre-idwc_mpi.x
	mpirun -np 5 $B/acfd_fwi_pre-idwc_mpi.x init_model=vp_layers_start.rsf par=$P/acfd_layers_idwc.p par=$P/acfd_layers_warp.p output=rcvdat_layers_monitor_pre-idwc.rsf



include Inc/SDR.bot
