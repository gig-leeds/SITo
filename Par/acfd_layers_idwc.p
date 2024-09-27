 # file id
 fid=rcvdat_layers_baseline_
 fidm=rcvdat_layers_monitor_
 rid=rid_layers
 sid=1

 # number of damping elements
 npad=400

 # source geometry
 src_zx=src_layers_g1.rsf
 # source data
 src_dat=srcdat_layers_g1.rsf

 # receiver geometry
 rcv_zx=rcv_layers_g1.rsf


 
 # Names of some of the inversion output files
 wrp-file=prewrp-img.rsf
 grd-fwi=grd_fwi.rsf
 grd-wrp=grd_wrp.rsf
 grd-tot=grd_tot.rsf

 
 #### inversion parameters ???
 # number of inversion iterations
 niter=6

 perturbp=0.003       
 waterlevel=100000              

 # weight parameters required in minimising cost function
 # eps_rcv :receiver residual weight
 eps_rcv=0                          
 # image residual weight                      
 eps_wll=1

