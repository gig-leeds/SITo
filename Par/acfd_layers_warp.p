 # Here comes warp stuff
 # tapering values to taper rtm images
 ntaper=25
 ntop=50
 nbot=50

 mute_angle=45

 # 2nd order Tikhonov regularisation parameter for warp inversion
 eps_tik=1.0e10
 # LSQR regularisation
 eps_model=1.0e-13

 shiftguess=.0
 nguess=1
 dguess=1
 # warping iterations
 niter_warp=500
 # Maximum number of iterations in LSQR in absence of convergence
 itnlim_warp=1000

