#!/usr/bin/env python

from subprocess import Popen, PIPE, STDOUT
import sys
import os
import re

# Use some backend
import matplotlib

debug = False


if len(sys.argv) > 1 and sys.argv[1] == "1":
	debug = True

if not debug:
	matplotlib.use('Pdf')
	extension = '.pdf'

from pylab import *
rc("font", family="sans-serif")
rc("font", size=12)

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from numpy import *
import scipy.integrate as integrate
import scipy.special as special

# Parameters for the md-simulation
N = 100
dt = 0.001
nequil = 500000
nproduct = 0
# binwidth = 0.05
startconf = 1
me = 50

Ts = [0.1, 1.0, 10.0]
# Ts = [0.1]

rcs = [1.5]

rhos = [0.1]

for T in Ts:
	for rc in rcs:
		for rho in rhos:
			# File/Foldername where output is stored
			filename = 'N_%(N).d_rho_%(rho).2f_T_%(T).2f_rc_%(rc).2f_dt_%(dt).3f_nequil_%(nequil).d_nproduct_%(nproduct).d_startconf_%(startconf).3f' % \
			{'N': N, 'rho': rho, 'T': T, 'rc': rc, 'dt': dt, 'nequil': nequil, 'nproduct': nproduct, 'startconf': startconf}

			# Whether to simulate later on. Is set to false, if output already exists
			simulate = True

			# Create the output folder if it doesn't exist
			try:
				dirname = 'run/'+filename
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
				command = ['time', '../../src/leapfrog', repr(N), repr(rho), repr(T), repr(rc), repr(dt), repr(startconf), repr(nequil), repr(nproduct)]
				print '	   With command:', ' '.join(command)
				result = Popen(command, stdout=output, stderr=STDOUT, cwd=dirname).communicate()


# Now visualize the results
figNum = 0

# Create new figure


# =============================================
# Energy-evolution for all temperatures
figNum += 1

fig = plt.figure(figNum, frameon=False, figsize=(9,12)) # , figsize=(9,11)
fig.clf() # clear figure
fig.subplots_adjust(hspace=.4)
fig.suptitle(r'Time-evolution of Energies during Equilibration')

subNum = 0

locs = [3, 6, 6]
for T in Ts:
	subNum += 1
	# File/Foldername where output is stored
	dirname = 'run/'+'N_%(N).d_rho_%(rho).2f_T_%(T).2f_rc_%(rc).2f_dt_%(dt).3f_nequil_%(nequil).d_nproduct_%(nproduct).d_startconf_%(startconf).3f' % \
	{'N': N, 'rho': rho, 'T': T, 'rc': rc, 'dt': dt, 'nequil': nequil, 'nproduct': nproduct, 'startconf': startconf}

	averages = csv2rec(dirname+'/outAvgFinal.txt', delimiter="\t", names=['t', 'T', 'Etot'])

	print 'Plotting energies vs time'

	ax = fig.add_subplot(len(Ts),1,subNum)
	# ax.text(2, txtpos[subNum-1], r'$T = %(T).3f$' % {'T': float(averages['T'][0])})
	title(r'Temperature $T = %(T).3f$' % {'T': float(averages['T'][0])}, fontsize=12)
	# ax.text(2.5, 2.3, r'$\rho = %(rho).3f$' % {'rho': rho})
	# ax.text(2.5, 1.7, r'$rc = %(rc).3f$' % {'rc': rc})

	ax.set_xlabel(r'Time $t$')
	ax.set_ylabel(r'Energy $E(t)$')

	data = csv2rec(dirname+'/outAveragesEquil.txt', delimiter="\t", names=['t', 'T', 'Tavg', 'Etot', 'Etotavg', 'Ekin', 'Ekinavg', 'Epot', 'Epotavg'])
	ax.plot(data['t'], data['Ekin'], 'g+', markevery = me, label = r'$E_{kin}$')
	ax.plot(data['t'], data['Epot'], 'bx', markevery = me, label = r'$E_{pot}$')
	ax.plot(data['t'], data['Etot'], 'r.', markevery = me, label = r'$E_{tot}$')

	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles, labels, loc = locs[subNum-1])

	ax.set_xscale('log')
	ax.set_xlim(1, float(data['t'][-1]))

	# ax.axis([0, 3.5, -1, 4])

	if not debug:
		plt.savefig('energies'+extension)


# =============================================
# Momentum-evolution for all temperatures
figNum += 1

fig = plt.figure(figNum, frameon=False, figsize=(9,12)) # , figsize=(9,11)
fig.clf() # clear figure
fig.subplots_adjust(hspace=.4)
fig.suptitle(r'Evolution of Momentum-Components during Equilibration')

subNum = 0
locs = [2, 2, 2]
for T in Ts:
	subNum += 1
	# File/Foldername where output is stored
	dirname = 'run/'+'N_%(N).d_rho_%(rho).2f_T_%(T).2f_rc_%(rc).2f_dt_%(dt).3f_nequil_%(nequil).d_nproduct_%(nproduct).d_startconf_%(startconf).3f' % \
	{'N': N, 'rho': rho, 'T': T, 'rc': rc, 'dt': dt, 'nequil': nequil, 'nproduct': nproduct, 'startconf': startconf}

	averages = csv2rec(dirname+'/outAvgFinal.txt', delimiter="\t", names=['t', 'T', 'Etot'])

	print 'Plotting momentum vs time'

	ax = fig.add_subplot(len(Ts),1,subNum)
	# ax.text(2, txtpos[subNum-1], r'$T = %(T).3f$' % {'T': float(averages['T'][0])})
	title(r'Temperature $T = %(T).3f$' % {'T': float(averages['T'][0])}, fontsize=12)
	# ax.text(2.5, 2.3, r'$\rho = %(rho).3f$' % {'rho': rho})
	# ax.text(2.5, 1.7, r'$rc = %(rc).3f$' % {'rc': rc})

	ax.set_xlabel(r'Time $t$')
	ax.set_ylabel(r'Component of momentum $p_i$')

	data = csv2rec(dirname+'/outAveragesEquil.txt', delimiter="\t", names=['t', 'T', 'Tavg', 'Etot', 'Etotavg', 'Ekin', 'Ekinavg', 'Epot', 'Epotavg', 'vx', 'vxavg', 'vy', 'vyavg'])
	ax.plot(data['t'], data['vx'], 'r.', markevery = me, label = r'$p_x$')
	ax.plot(data['t'], data['vy'], 'bx', markevery = me, label = r'$p_y$')

	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles, labels, loc = locs[subNum-1])

	# ax.set_xscale('log')
	# ax.set_xlim(1, float(data['t'][-1]))

	# ax.axis([0, 3.5, -1, 4])

	if not debug:
		plt.savefig('momenta'+extension)


# =============================================
# Snapshots
locs = [2, 3, 3]
for T in Ts:
	figNum += 1

	fig = plt.figure(figNum, frameon=False, figsize=(7,12)) # , figsize=(9,11)
	fig.clf() # clear figure
	fig.subplots_adjust(hspace=.2)
	fig.suptitle(r'Snapshot of system configuration at $T=%(T).2f$' % {'T': T})
	subNum = 0;
	for config in ['start', 'equilibrated']:
		subNum += 1
		# File/Foldername where output is stored
		dirname = 'run/'+'N_%(N).d_rho_%(rho).2f_T_%(T).2f_rc_%(rc).2f_dt_%(dt).3f_nequil_%(nequil).d_nproduct_%(nproduct).d_startconf_%(startconf).3f' % \
		{'N': N, 'rho': rho, 'T': T, 'rc': rc, 'dt': dt, 'nequil': nequil, 'nproduct': nproduct, 'startconf': startconf}

		print 'Plotting snapshots'

		ax = fig.add_subplot(2,1,subNum)
		if config == 'start':
			title(r'Initial configuration $t = %(t).2f$' % {'t': 0}, fontsize=12)
		else:
			title(r'Final configuration $t = %(t).2f$' % {'t': nequil*dt}, fontsize=12)

		ax.set_xlabel(r'$x$')
		ax.set_ylabel(r'$y$')

		data = csv2rec(dirname+'/outCoords_'+config+'.txt', delimiter="\t", names=['N', 'x', 'y'])
		print data['x'][0]
		for i in range(0, len(data['x'])):
			ax.add_patch(Circle((data['x'][i], data['y'][i]), .5))

		ax.set_aspect(1)

		ax.axis([-1, max(data['x'])+1, -1, max(data['y'])+1])

		if not debug:
			filename = 'snapshot_T_%(T).2f' % {'T': T}
			filename = filename.replace('.', '_')+extension
			plt.savefig(filename)


if debug:
	plt.show()