#!/usr/bin/env python

from subprocess import Popen, PIPE, STDOUT
import sys
import os
import re

# Use some backend
import matplotlib
matplotlib.use('Pdf')
extension = '.pdf'

from pylab import *
rc("font", family="sans-serif")
rc("font", size=16)

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from numpy import *
import scipy.integrate as integrate
import scipy.special as special

# Parameters for the md-simulation
N = 100
T = 1.8
rc = 2.5
dt = 0.01
nequil = 1000
nproduct = 5000
binwidth = 0.05

rhos = [0.84, 0.5, 0.1]

for rho in rhos:
	# File/Foldername where output is stored
	filename = 'N_%(N).d_rho_%(rho).2f_T_%(T).2f_rc_%(rc).2f_dt_%(dt).3f_nequil_%(nequil).d_nproduct_%(nproduct).d_binwidth_%(binwidth).3f' % \
	{'N': N, 'rho': rho, 'T': T, 'rc': rc, 'dt': dt, 'nequil': nequil, 'nproduct': nproduct, 'binwidth': binwidth}

	# Whether to simulate later on. Is set to false, if output already exists
	simulate = True

	# Create the output folder if it doesn't exist
	try:
		dirname = filename
		os.makedirs(dirname)
	except OSError, e:
		print >>sys.stderr, 'Execution failed:', e
		print 'We dont simulate, just create the Plots!'
		simulate = False

	# If output wasn't already there, start a simulation
	if simulate:
		# File for general simulation output
		output = open(dirname+'/'+filename+'_output', 'w+', 1)
		print 'Start simulation...'
		command = ['time', '../leapfrog', repr(N), repr(rho), repr(T), repr(rc), repr(dt), repr(nequil), repr(nproduct), repr(binwidth)]
		print '	   With command:', ' '.join(command)
		result = Popen(command, stdout=output, stderr=STDOUT, cwd=dirname).communicate()

# Now visualize the results
figNum = 0

# Create new figure
figNum += 1

# w, h = figaspect(2.)
# fig = Figure(figsize=(w,h))
# ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
# ax.imshow(A, **kwargs)

fig = plt.figure(figNum, frameon=False, figsize=(6,10))
fig.clf() # clear figure
fig.subplots_adjust(hspace=.3)
fig.suptitle(r'Pair-correlation-function $g(r)$')

subNum = 0

for rho in rhos:
	# File/Foldername where output is stored
	dirname = 'N_%(N).d_rho_%(rho).2f_T_%(T).2f_rc_%(rc).2f_dt_%(dt).3f_nequil_%(nequil).d_nproduct_%(nproduct).d_binwidth_%(binwidth).3f' % \
	{'N': N, 'rho': rho, 'T': T, 'rc': rc, 'dt': dt, 'nequil': nequil, 'nproduct': nproduct, 'binwidth': binwidth}

	averages = csv2rec(dirname+'/outAvgFinal.txt', delimiter="\t", names=['t', 'T','Tt','Etot','Etott'])

	print 'Plotting g(r) vs r'

	subNum += 1
	ax = fig.add_subplot(len(rhos),1,subNum)
	ax.text(2.5, 3, r'$T = %(T).3f$' % {'T': float(averages['T'][0])})
	ax.text(2.5, 2.3, r'$\rho = %(rho).3f$' % {'rho': rho})

	ax.set_xlabel(r'$r$')
	ax.set_ylabel(r'$g(r)$')

	data = csv2rec(dirname+'/outGr.txt', delimiter="\t", names=['x','y'])
	ax.plot(data['x'], data['y'], 'ro')

	ax.axis([0, 3.5, -1, 4])

plt.savefig('g_r'+extension)
plt.show()


	# pdf, bins, patches = ax.hist(x, 50, range=[xmin, xmax], normed=1, facecolor='green', alpha=0.75)
	# 
	# # Write all parameters into the plot title
	# if doExchanges:
	# 	title = r'$T = %(temp).3f,\, x_0 = %(x0).2f,\, N_\mathrm{exch} = %(nexchanges)d,\, n_\mathrm{step} = %(nsteps)d,\, \gamma = %(gamma)d,\, \mathrm{d}t= %(dt).3f$' %  \
	# 	{'temp': temp, 'gamma': gamma, 'x0': x0, 'nexchanges': nexchanges, 'nsteps': nsteps, 'dt': dt}
	# else:
	# 	title = r'$T = %(temp).3f,\, x_0 = %(x0).2f,\, N = %(nexchanges)d,\, \gamma = %(gamma)d,\, \mathrm{d}t= %(dt).3f$' %  \
	# 	{'temp': temp, 'gamma': gamma, 'x0': x0, 'nexchanges': nexchanges*nsteps, 'dt': dt}
	# plt.title(title)
	# 
	# # add expected distribution
	# # t = arange(bins[0], bins[-1], 0.01)
	# t = arange(xmin, xmax, 0.01)
	# w = (a*t**4 + b*t**2 + c*t);
	# # Calculate the normalization constant, i.e. integrate over the whole axis
	# CN,CNerr = integrate.quad(lambda t: exp(-(a*t**4 + b*t**2 + c*t)/temp), -inf, inf)
	# print "normalization constant "+repr(CN)
	# ax.plot(t, exp(-w/temp)/CN, 'b-', linewidth=lw)
	# 
	# ax.set_ylim(0, ylims[i])
	# i += 1
	# 
	# # ax.set_ylabel(r'$\exp[-\beta\, U(x)]$')
	# ax.set_ylabel(r'$1/Z\, \exp[-\beta\, U(x)]$')
	# ax.set_xlabel(r'Position $x$')
	# 
	# 
	# # also plot the potential
	# ax = plt.subplot(212)
	# ax.grid(True)
	# ax.plot(t, w, 'r-', linewidth=lw)
	# 
	# ax.set_ylabel(r'$g(r)$')
	# ax.set_xlabel(r'$r$')
	# 
	# # Save plot to file
	# plt.savefig(dirname+'/'+outfile+extension)


	# # Create another figure with the energy propability distributions
	# figNum += 1
	# fig = plt.figure(figNum, frameon=False)
	# fig.clf() # clear figure
	# 
	# # Title remains the same
	# plt.title(title)
	# 
	# # New plot
	# ax = fig.add_subplot(111)
	# 
	# # Get the replica index of each step from the data read above
	# x = data[:,0]
	# x = a*x**4 + b*x**2 + c*x
	# pdf, bins, patches = ax.hist(x, 100, range=[-4, 1], normed=1, facecolor='green', alpha=0.75, histtype='step')
	# 
	# ax.grid(True)
	# ax.set_yscale('log')
	# ax.set_ylim(0.1, 10)
	# ax.set_ylabel('Probability P')
	# ax.set_xlabel(r'Potential Energy U')
	# 
	# # Save plot to file
	# plt.savefig(dirname+'/energy_probability_'+outfile+extension)
	# 
	# 
	# # If replica exchange method wasn't applied, we're done
	# if not doExchanges:
	# 	continue
	# 
	# # Create another figure for the replica exchange visualization
	# figNum += 1
	# fig = plt.figure(figNum, frameon=False)
	# fig.clf() # clear figure
	# 
	# # Title remains the same
	# plt.title(title)
	# 
	# # Get the replica index of each step from the data read above
	# x = data[:,1]
	# t = range(0, alen(x))
	# 
	# # New plot
	# ax = fig.add_subplot(111)
	# ax.plot(t, x, 'r.', markevery=mev)
	# 
	# ax.grid(True)
	# ax.set_ylim(-.5, 7.5)
	# ax.set_ylabel('Replica')
	# ax.set_xlabel(r'Exchangestep $n$')
	# 
	# # Save plot to file
	# plt.savefig(dirname+'/replica_at_T_'+outfile+extension)


# # Now draw some plots for the temperature evolution of each replicum
# if doExchanges:
# 	for outfile in os.listdir(dirname):
# 		# Filter outputfiles
# 		match = re.search("^output_N_([0-9]+)$", outfile)
# 		if not match:
# 			continue
# 
# 		# Extract replica-id from filename
# 		N = int(match.group(1))
# 
# 		print "Plotting temperatures, that replicum "+repr(N)+" went through."
# 
# 		# Create another figure for the temperature exchange
# 		figNum += 1
# 		fig = plt.figure(figNum, frameon=False)
# 		fig.clf() # clear figure
# 		fig.subplots_adjust(hspace=.3)
# 		# fig.set_figheight(fig.get_figheight()*1.3)
# 
# 		# Plot title
# 		plt.title('Replicum '+repr(N))
# 
# 		# Get the temperatures for the current replicum
# 		data = loadtxt(dirname+'/'+outfile)
# 		x = data[:,0]
# 		t = range(0, alen(x))
# 
# 		# New plot with temperatures
# 		ax = plt.subplot(211)
# 		ax.plot(t, x, 'b-', markevery=mev)
# 		ax.set_yscale('log')
# 		ax.set_ylim(0.08, 4)
# 
# 		ax.grid(True)
# 		ax.set_ylabel(r'Temperature $T$')
# 		ax.set_xlabel(r'Exchangestep $n$')
# 
# 		# New plot with position
# 		x = data[:,1]
# 		t = range(0, alen(x))
# 
# 		ax = plt.subplot(212)
# 		ax.plot(t, x, 'r-', markevery=mev)
# 		ax.set_ylim(xmin, xmax)
# 
# 		ax.grid(True)
# 		ax.set_ylabel(r'Position $x$')
# 		ax.set_xlabel(r'Exchangestep $n$')
# 
# 		# Save plot to file
# 		plt.savefig(dirname+'/replica_N_'+outfile+extension)
# 
