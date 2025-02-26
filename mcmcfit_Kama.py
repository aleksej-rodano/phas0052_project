# Setup to use emcee for MCMC fitting -- example provided to the PHAS0052 24/25 Term 2 project team
# Mihkel Kama, 2017-2025
# mkama@ucl.ac.uk


# Import libraries
import os
import math
import numpy
import emcee
import matplotlib.pyplot as plt
import corner
import scipy
from scipy import integrate, special, stats
import time

#########################
Nwalkers 			= 100	# number of walkers
Nsteps				= 3000	# number of steps for walkers to take
Nburn 				= 20	# number of burn-in steps
width_parameter_a	= 3.0	# parameter controlling the width of distribution of parameter jumps
dosaveresults		= True
verbose 			= False
resultsdir 			= 'results/'
#########################

# Constants and unit conversion factors
auSI			= 1.496e+11 # [m] 1 astronomical unit
pcSI 			= 3.086e+16 # [m] 1 parsec

#### Parameter initialization, max ranges, labels (replace with your own)
# Initial guess values for all parameters
params_init 	= { 'T0':17.0, 'q1':0.5, 'q2':0.5, 'Sigma0':-0.2, 'p1':1.0, 'p2':3.0, 'rT0':90.0*auSI, 'rSbreak':110.0*auSI, 'dist':122.0*pcSI, 'rin':0.5*auSI, 'rout':800.0*auSI }
Npars 			= len( params_init.keys() )
# Create arrays of those initial values and their names
pos_ref 		= numpy.array( [ params_init['T0'], params_init['q1'], params_init['q2'], params_init['Sigma0'], params_init['p1'], params_init['p2'], params_init['rT0'], params_init['rSbreak'] ] )
labels 			= ['T0','q1','q2','Sigma0','p1','p2','rT0','rSbreak']
# Define the range of allowed variation for the initial parameter values around the initial guesses in pos_ref
parinit_minfrac	= numpy.array( [ 0.20, 0.1, 0.1, -3.0, 0.1, 0.1, 0.01, 0.01 ] )
parinit_maxfrac	= numpy.array( [ 10.0, 3.0, 1.9, 4.00, 1.9, 1.9, 1.9, 1.9 ] )


#### Plotting tools
# plot dimensions, font sizes, line thicknesses etc
px1 		= 14
py1 		= 5
px2 		= 7
py2 		= 9
fontsize 	= 22
fsize2		= 20
fanno		= 18
def plotmaker( xsize, ysize ):
	# Plot
	plt.figure( num=1, figsize=(xsize,ysize), dpi=80 )
	plt.rc('axes', linewidth=2)
	plt.rc('xtick.minor', size=4, width=2)
	plt.rc('xtick.major', size=8, width=2)
	plt.rc('ytick.minor', size=4, width=2)
	plt.rc('ytick.major', size=8, width=2)
	ax = plt.gca()
	for tick in ax.xaxis.get_major_ticks():
		tick.label1.set_fontsize(fontsize)
	for tick in ax.yaxis.get_major_ticks():
		tick.label1.set_fontsize(fontsize)	
	return


#### Fitting
def fitdata( pos_ref = pos_ref ):
	print ("-- running emcee MCMC...")
	time_init 						= time.time() # timer for our code runtime
	# Create a multidimensional array of initial parameter guesses for each walker by 
	# taking the first initial guesses in pos_ref as the center of a normal distribution 
	# whose width relates to parinit_minfrac and _maxfrac:
	pos 							= [ numpy.random.normal( pos_ref, numpy.abs( parinit_maxfrac - parinit_minfrac ), Npars ) for ii in range( Nwalkers ) ]
	# Initialise and run the sampler in emcee:
	sampler 						= emcee.EnsembleSampler( Nwalkers, Npars, lnprob, a=width_parameter_a )
	for i, result in enumerate( sampler.sample( pos, iterations=Nsteps ) ):
		if 100.0*i/Nsteps % 10 == 0:
			print("{0:5.1f}% with {1:5.1f}min elapsed".format( 100 * float( i ) / Nsteps, ( time.time() - time_init ) / 60.0 ) )
			print(" Max and mean acceptance fraction: {0:.3f}, {1:.3f}".format( numpy.max( sampler.acceptance_fraction ), numpy.mean( sampler.acceptance_fraction ) ) )
			Ninit 	 			= 10
			if i > Nburn: Ninit = Nburn
			fitparams	 		= getFitParams( sampler, Ninit )
			# Create a plot of the current results compared to the data:
			try:
				plotvis( data_cont, fitparams=fitparams, outfile=resultsdir + 'visfit_Nw%i_Ns%i_Nrad%i/fig_plotvis_step%i.pdf' % ( Nwalkers, Nsteps, Nrad, i ) )
			except:
				print ("!! problem generating plot using the mean parameter values :(")
				print "!! fitparams: ", fitparams
	return sampler


#### MCMC

def lnprob( pos_ref ):
		iT0, iq1, iq2, iSb, ip1, ip2, irT0, irSb	= pos_ref
		# set your priors here:
		if ( not 0.01*auSI < params_init['rin'] < irSb < params_init['rout'] ) \
		or ( not 0.01*auSI < params_init['rin'] < irT0 < params_init['rout'] ) \
		or ( not -1.0 <= iq1 < 3.0 ) \
		or ( not -1.0 <= iq2 < 10.0 ) \
		or ( not -1.0 <= ip1 < 3.0 ) \
		or ( not -1.0 <= ip2 < 10.0 ) \
		or ( not 2.7 < iT0 < 500.0 ) \
		or ( not -3.0 < iSb < 4.0 ):
			return -numpy.inf
		# define the parameters and assign values
		params 	= { 'T0':iT0, 'q1':iq1, 'q2':iq2, 'Sigma0':iSb, 'p1':ip1, 'p2':ip2, 'rT0':irT0, 'rSbreak':irSb }
		# calculate the figure of merit
		like_chi_sq			= -0.5 * numpy.sum( numpy.array( [ ( ( mod - obs ) / obserr )**2 + numpy.log(obserr**2) + numpy.log(2.0*numpy.pi) for obs, obserr, mod in zip( data_cont['visRe'], data_cont['visRe_std'], Visis ) ] ) )
		if verbose:
			print "-- "
			print "  params: ", params
			print "  lnprob_prefactor = %.4e		like_chi_sq = %.4e" % ( lnprob_prefactor, like_chi_sq )	# lnprob_prefactor is calculated later on
			print "-- "
		return lnprob_prefactor + like_chi_sq



def getFitParams( sampler1, Ninit ):
	fitparams 		= {	'T0':numpy.mean(sampler1.chain[:, Ninit:, 0]), \
							'q1':numpy.mean(sampler1.chain[:, Ninit:, 1]), \
							'q2':numpy.mean(sampler1.chain[:, Ninit:, 2]), \
							'Sigma0':numpy.mean(sampler1.chain[:, Ninit:, 3]), \
							'p1':numpy.mean(sampler1.chain[:, Ninit:, 4]), \
							'p2':numpy.mean(sampler1.chain[:, Ninit:, 5]), \
							'rT0':numpy.mean(sampler1.chain[:, Ninit:, 6]), \
							'rSbreak':numpy.mean(sampler1.chain[:, Ninit:, 7]), \
							'dist':params_init['dist'], \
							'rin':params_init['rin'], \
							'rout':params_init['rout'] }	# the last three parameters were fixed in the implementation I had for this, hence they're just read from the params_init structure
	return fitparams


#### Some plotting functions
def plotter( samplerdata ):
	plotmaker( 12, 10 )
	for i_var in range(5):
		plt.subplot( 5, 1, i_var+1 )
		for i in range(Nwalkers):
			plt.plot( range(Nsteps), samplerdata.chain[i,:,i_var], 'k-' )
	return

def plotter2( samplerdata, maindir=resultsdir, labels=labels ):
	for ii in range( Npars ):
		print "%s 	mean=%.6e 	stddev=%.6e" % ( labels[ii], numpy.median(samplerdata.chain[:, Nburn:, ii]), numpy.std(samplerdata.chain[:, Nburn:, ii]) )
	samples = samplerdata.chain[:, Nburn:, :].reshape((-1, Npars))
	fig = corner.corner( samples, labels=labels )
	fig.savefig( maindir + 'fig_latest_corner.pdf')
	plt.close('all')
	return

def plotter3( samplerdata, labels=labels ):
	for ii in range( Npars ):
		print "%s 	mean=%.2e 	stddev=%.2e" % ( labels[ii], numpy.median(samplerdata.chain[:, Nburn:, ii]), numpy.std(samplerdata.chain[:, Nburn:, ii]) )
	samples = samplerdata.chain[:, Nburn:, :].reshape((-1, Npars))
	fig = corner.corner( samples, labels=labels )
	fig.savefig('fig_latest_corner3.pdf')
	fig.close('all')
	return


#### Modelling

# I had some functions here that defined my model, you can have the lightcurve function (modelling_somethingsomething)



#### Support routines
# save the entire sampler to a binary file with a name assigned by the outFile string variable
def savesampler( samplerdata, outFile ):
	samples = samplerdata.chain[:, :, :].reshape((-1, Npars))
	numpy.save( outFile, samples )
	return
# check whether a directory exists; if not, then create it
def ensure_dir( dirName ):
	if not os.path.exists( dirName ):	os.makedirs( dirName )
	return

print "-- reading data..."

# I had functions to read the data from a file

# this only needs to be calculated once from the data so I've put it here rather than in the lnprob function which gets called a lot
lnprob_prefactor 		= - (1.0/2.0) * numpy.sum( [ numpy.log( 2.*math.pi * visRe_std**2 ) for visRe_std in data_cont['visRe_std'] ] )
print( "-- lnprob_prefactor = ", lnprob_prefactor )

if dovisfit:
	print "-- fitting..."
	sampler		= fitdata()
	try:
		print( "Sampler autocorrelation (max and mean): {0:.3f}, {1:.3f}".format( numpy.int(max(sampler.acor)), numpy.int(numpy.mean(sampler.acor)) ) )
		print( "All autocorrelations: ", sampler.acor )
	except:
		print "!! WARNING: failed to calculate the autocorrelation time(s). Maybe the chain is too short?"
	print( "Mean acceptance fractions: {0:.3f}".format( numpy.mean( sampler.acceptance_fraction ) ) )
	if verbose: print( "All acceptance fractions: ", sampler.acceptance_fraction)
	if dosaveresults:
		savesampler( sampler, resultsdir + 'visfit_Nw%i_Ns%i_Nrad%i/visfit_Nw%i_Ns%i_Nrad%i' % ( Nwalkers, Nsteps, Nrad, Nwalkers, Nsteps, Nrad ) )
		print "-- plotting visibilities..."
		fitparams 		= getFitParams( sampler, Nburn )
		try:
			plotvis( data_cont, fitparams=fitparams, outfile=resultsdir + 'visfit_Nw%i_Ns%i_Nrad%i/fig_plotvis_mean.pdf' % ( Nwalkers, Nsteps, Nrad )  )
		except:
				print "!! problem generating plot using the mean parameter values :("
				print "!! fitparams: ", fitparams
		for i in range( Npars ):
			plt.hist( sampler.flatchain[:,i], 100, color="k", histtype="step" )
			plt.title( "Dimension {0:d}".format( i ) )
			plt.savefig( '' + resultsdir + 'visfit_Nw%i_Ns%i_Nrad%i/fig_flatchain_par%i.pdf' % ( Nwalkers, Nsteps, Nrad, i ) )
			plt.close()
			plt.plot( sampler.chain[:,:,i].T, '-', color='k', alpha=0.3 )
			plt.savefig( '' + resultsdir + 'visfit_Nw%i_Ns%i_Nrad%i/fig_allwalkers_par%i.pdf' % ( Nwalkers, Nsteps, Nrad, i ) )
			plt.close()
		plotter2( sampler, maindir=resultsdir + 'visfit_Nw%i_Ns%i_Nrad%i/' % ( Nwalkers, Nsteps, Nrad ) )

print "-- fitdata completed."
