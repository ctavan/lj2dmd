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
nequil = 200000
nproduct = 0
# binwidth = 0.05
startconf = 1

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
	ax.plot(data['t'], data['Ekin'], 'g+', markevery = 10, label = r'$E_{kin}$')
	ax.plot(data['t'], data['Epot'], 'bx', markevery = 10, label = r'$E_{pot}$')
	ax.plot(data['t'], data['Etot'], 'r.', markevery = 10, label = r'$E_{tot}$')

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
locs = [2, 3, 3]
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
	ax.plot(data['t'], data['vx'], 'r.', markevery = 10, label = r'$p_x$')
	ax.plot(data['t'], data['vy'], 'bx', markevery = 10, label = r'$p_y$')

	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles, labels, loc = locs[subNum-1])

	# ax.set_xscale('log')
	# ax.set_xlim(1, float(data['t'][-1]))

	# ax.axis([0, 3.5, -1, 4])

	if not debug:
		plt.savefig('momenta'+extension)

# # =============================================
# # g(r) for 3 denisities
# for rc in rcs:
# 	subNumC += 1
# 	subNumR = 0
# 	for rho in rhos_gr:
# 		# File/Foldername where output is stored
# 		dirname = 'run/'+'N_%(N).d_rho_%(rho).2f_T_%(T).2f_rc_%(rc).2f_dt_%(dt).3f_nequil_%(nequil).d_nproduct_%(nproduct).d_startconf_%(startconf).3f' % \
# 		{'N': N, 'rho': rho, 'T': T, 'rc': rc, 'dt': dt, 'nequil': nequil, 'nproduct': nproduct, 'startconf': startconf}
# 
# 		averages = csv2rec(dirname+'/outAvgFinal.txt', delimiter="\t", names=['t', 'T', 'Etot'])
# 
# 		print 'Plotting g(r) vs r'
# 
# 		subNumR += 1
# 		ax = fig.add_subplot(len(rhos_gr),len(rcs),(2*subNumR-1)+subNumC-1)
# 		ax.text(2.5, 3, r'$T = %(T).3f$' % {'T': float(averages['T'][0])})
# 		ax.text(2.5, 2.3, r'$\rho = %(rho).3f$' % {'rho': rho})
# 		ax.text(2.5, 1.7, r'$rc = %(rc).3f$' % {'rc': rc})
# 
# 		ax.set_xlabel(r'$r$')
# 		ax.set_ylabel(r'$g(r)$')
# 
# 		data = csv2rec(dirname+'/outGr.txt', delimiter="\t", names=['x','y'])
# 		ax.plot(data['x'], data['y'], 'ro')
# 
# 		ax.axis([0, 3.5, -1, 4])
# 
# 	if not debug:
# 		plt.savefig('g_r'+extension)
# 
# # =============================================
# # Plot time evolution of T
# figNum += 1
# 
# rc = 2.5
# 
# fig = plt.figure(figNum, frameon=False)
# fig.clf() # clear figure
# fig.subplots_adjust(hspace=.3)
# fig.suptitle(r'Instantaneous temperature and mean temperature')
# 
# ax = fig.add_subplot(1,1,1)
# 
# ax.set_xlabel(r'$t$')
# ax.set_ylabel(r'$T$')
# 
# rho = 0.5
# dirname = 'run/'+'N_%(N).d_rho_%(rho).2f_T_%(T).2f_rc_%(rc).2f_dt_%(dt).3f_nequil_%(nequil).d_nproduct_%(nproduct).d_startconf_%(startconf).3f' % \
# {'N': N, 'rho': rho, 'T': T, 'rc': rc, 'dt': dt, 'nequil': nequil, 'nproduct': nproduct, 'startconf': startconf}
# 
# data = csv2rec(dirname+'/outAverages.txt', delimiter="\t", names=['t', 'Tt','T','Etott','Etot'])
# 
# ax.plot(data['t'], data['Tt'], 'ro', data['t'], data['T'], 'b-', linewidth=3, markersize=4)
# 
# ax.set_ylim(T-0.7, T+0.7)
# 
# if not debug:
# 	plt.savefig('T_mean'+extension)
# 
# 
# # =============================================
# # g(r) for 3 denisities
# dp08 = []
# dp18 = []
# p08 = {}
# p18 = {}
# Tlist08 = {}
# Tlist18 = {}
# rhoss08 = []
# rhoss18 = []
# rho_invs08 = {}
# rho_invs18 = {}
# for T in Ts:
# 	p = {'r': [], 'lj': []}
# 	Tlist = {'r': [], 'lj': []}
# 	rho_invs = {'r': [], 'lj': []}
# 	for rc in rcs:
# 		rhoss = []
# 		for rho in rhos_p_inv:
# 			if rc == 2.5:
# 				rcname = 'lj'
# 			else:
# 				rcname = 'r'
# 
# 			# File/Foldername where output is stored
# 			dirname = 'run/'+'N_%(N).d_rho_%(rho).2f_T_%(T).2f_rc_%(rc).2f_dt_%(dt).3f_nequil_%(nequil).d_nproduct_%(nproduct).d_startconf_%(startconf).3f' % \
# 			{'N': N, 'rho': rho, 'T': T, 'rc': rc, 'dt': dt, 'nequil': nequil, 'nproduct': nproduct, 'startconf': startconf}
# 
# 			averages = csv2rec(dirname+'/outAvgFinal.txt', delimiter="\t", names=['t', 'T', 'Etot', 'rho', 'rhoinv', 'p'])
# 
# 			p[rcname].append(float(averages['p'][0]))
# 			rho_invs[rcname].append(float(averages['rhoinv'][0]))
# 			rhoss.append(float(averages['rho'][0])**2.0)
# 			Tlist[rcname].append(float(averages['T'][0]))
# 
# 
# 	dp = [x - y for (x,y) in zip(p['r'],p['lj'])]
# 	if T == 0.8:
# 		dp08 = dp
# 		p08 = p
# 		Tlist08 = Tlist
# 		rhoss08 = rhoss
# 		rho_invs08 = rho_invs
# 	else:
# 		dp18 = dp
# 		p18 = p
# 		Tlist18 = Tlist
# 		rhoss18 = rhoss
# 		rho_invs18 = rho_invs
# 
# print 'Plotting p(rho)'
# 
# Tavg = (float(sum(Tlist18['lj']))/float(len(Tlist18['lj']))+float(sum(Tlist18['r']))/float(len(Tlist18['r'])))/2.0
# 
# figNum += 1
# 
# fig = plt.figure(figNum, frameon=False, figsize=(6,6))
# fig.clf() # clear figure
# fig.subplots_adjust(hspace=.3)
# fig.suptitle(r'Pressure $P$ as a function of $1/\rho$')
# 
# ax = fig.add_subplot(1,1,1)
# ax.text(7, 10, r'$T = %(T).3f$' % {'T': Tavg})
# 
# ax.set_xlabel(r'$1/\rho$')
# ax.set_ylabel(r'$P$')
# 
# x = arange(0.5, 10, 0.01)
# y = 1/x*Tavg;
# 
# ax.plot(rho_invs18['lj'], p18['lj'], 'ro', label=r'$r_c = 2.5$')
# ax.plot(rho_invs18['r'], p18['r'], 'gs', label=r'$r_c = 2^{1/6}$')
# ax.plot(x, y, 'b-', label='Ideal gas')
# ax.legend()
# 
# ax.axis([0, 11, -1, 12])
# 
# if not debug:
# 	plt.savefig('p_rho_T18'+extension)
# 
# Tavg = (float(sum(Tlist08['lj']))/float(len(Tlist08['lj']))+float(sum(Tlist08['r']))/float(len(Tlist08['r'])))/2.0
# 
# figNum += 1
# 
# fig = plt.figure(figNum, frameon=False, figsize=(6,6))
# fig.clf() # clear figure
# fig.subplots_adjust(hspace=.3)
# fig.suptitle(r'Pressure $P$ as a function of $1/\rho$')
# 
# ax = fig.add_subplot(1,1,1)
# ax.text(7, max(p08['lj'])-2, r'$T = %(T).3f$' % {'T': Tavg})
# 
# ax.set_xlabel(r'$1/\rho$')
# ax.set_ylabel(r'$P$')
# 
# x = arange(0.5, 10, 0.01)
# y = 1/x*Tavg;
# 
# ax.plot(rho_invs08['lj'], p08['lj'], 'ro', label=r'$r_c = 2.5$')
# ax.plot(rho_invs08['r'], p08['r'], 'gs', label=r'$r_c = 2^{1/6}$')
# ax.plot(x, y, 'b-', label='Ideal gas')
# ax.legend()
# 
# ax.axis([0, 11, -0.5, 4])
# 
# if not debug:
# 	plt.savefig('p_rho_T08'+extension)
# 
# 
# 
# figNum += 1
# 
# fig = plt.figure(figNum, frameon=False, figsize=(6,6))
# fig.clf() # clear figure
# fig.subplots_adjust(hspace=.3)
# fig.suptitle(r'Pressure difference $P_R - P_{LJ}$')
# 
# ax = fig.add_subplot(1,1,1)
# 
# ax.set_xlabel(r'$\rho^2$')
# ax.set_ylabel(r'$\Delta P$')
# 
# # x = arange(0.5, 10, 0.01)
# # y = 1/x*Tavg;
# 
# ax.plot(rhoss08, dp08, 'ro', label=r'$T = 0.8$')
# ax.plot(rhoss18, dp18, 'gs', label=r'$T = 1.8$')
# # ax.plot(x, y, 'b-', label='Ideal gas')
# ax.legend()
# 
# # delp = []
# # for p in p08
# 
# 
# 
# if not debug:
# 	plt.savefig('vdw'+extension)


if debug:
	plt.show()