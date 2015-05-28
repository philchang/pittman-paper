import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as pl
import numpy as np
import matplotlib.backends.backend_pdf as mat_pdf
import argparse
import glob
import sys
import os

parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N3', type=int)
parser.add_argument('--bigsphere', action='store_true')
parser.add_argument('--smallsphere', action='store_true')
parser.add_argument('--particle', action='store_true')
parser.add_argument('--noparticle', action='store_true')
parser.add_argument('--shell', action='store_true')
parser.add_argument('--shellsphere', action='store_true')

args = parser.parse_args()

withBigSphere = args.bigsphere
withSmallSphere = args.smallsphere
withParticle = args.particle
withNoParticle = args.noparticle
withShell = args.shell
withShellSphere = args.shellsphere

mp = 1.6e-24
pi = np.pi
parsec = 3e18
sec_in_yr = np.pi* 1e7
Msun = 2e33
G = 6.67e-8

if (withBigSphere):
	compare_files = 'bigsphere_'
	vel_Plot_label = '3.0 pc'

if (withSmallSphere):
	compare_files = 'smallsphere_'
	vel_Plot_label = '0.01 pc'

if (withParticle):
	compare_files = 'part_'
	vel_Plot_label = 'Bulk by Particle Vel'	

if (withNoParticle):
	compare_files = 'nopart'
	vel_Plot_label = 'No Particles'	

if (withShell):
	compare_files = 'shell_'
	vel_Plot_label = 'Bulk by Shell'

if (withShellSphere):
	compare_files = 'shellsphere_'
	vel_Plot_label = 'Bulk by Sphere in Shell'

file_prefix = 'rad_profile'
quad = os.getcwd()[-5:]
# The file looks like this:
#'{file_prefix}{framestep}_{compare_files}_{particle_number}.out'
# i.e. rad_profile_0218_part_000.out

for framestep in range(args.start,args.end,args.step) :   

	file_exist = glob.glob('{0}_{1}_{2:04d}_{3}*.out'.format(file_prefix, quad, framestep, compare_files))
	if not file_exist:
		print 'File: "{0}_{1}_{2:04d}_{3}*.out" does not exist!'.format(file_prefix, quad, framestep, compare_files)
		sys.exit()

	particle_list = glob.glob('{0}_{1}_{2:04d}_{3}*.out'.format(file_prefix, quad, framestep, compare_files))
	if compare_files != 'nopart':
		print "Looking at each particle"
		last_particle = particle_list[-1]
		last_particle = int(last_particle[-7:-4])
		print 'Looking at: {0}_{1}_{2:04d}_{3}*.out'.format(file_prefix, quad, framestep, compare_files)
		print 'The total number of star particles in this timestep is: ', last_particle + 1
		for particle_number in range(last_particle + 1):
			print "Plotting particle", particle_number + 1, 'of', last_particle + 1
			rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin = np.loadtxt( "{0}_{1}_{2:04d}_{3}{4:03d}.out".format(file_prefix, quad, framestep, compare_files, particle_number ), unpack=True)

			pdf_timestep = mat_pdf.PdfPages("framestep_{0}_{1}_{2:04d}_{3}{4:03d}.pdf".format(file_prefix, quad, framestep, compare_files, particle_number))
#			print vrmsbin-vmagbin

			pl.clf()
			pl.loglog( rbin, -vrbin/1e5, 'b') 
			pl.loglog( rbin, vrmsbin/1e5, 'g.-') 
			pl.loglog( rbin, vmagbin/1e5)
			pl.loglog( rbin, vKbin/1e5, 'r--')
			pl.loglog( rbin, vphi_magbin/1e5, 'k+')

			pl.ylabel('Velocity ($km/s$)')
			pl.xlabel('Radius ($pc$)')
			pl.ylim(1e-2, 10e0)
			pl.xlim(1e-3, 3e0)
			pl.title('{0} framestep {1:04d}'.format(vel_Plot_label, framestep))
			pdf_timestep.savefig()
			
			pl.clf()
			# rhobin should be in g/cm**3
			pl.loglog(rbin, rhobin/mp) #plotting number density
			pl.ylabel('Density (${N}/{cm^3}$)')
			pl.xlabel('Radius ($pc$)')
			pdf_timestep.savefig()
			
			pl.clf()
			pl.loglog( rbin, (rbin*parsec)*vphi_magbin/np.sqrt(G*mTbin*(parsec*rbin)))
			pl.ylabel('Ratio of specific Ang momentum')
			pl.xlabel('Radius ($pc$)')
			pdf_timestep.savefig()

			pl.clf()
			pl.loglog( rbin, mTbin)
			pl.ylabel('Total Mass ($g$)')
			pl.xlabel('Radius ($pc$)')
			pdf_timestep.savefig()
			
			pl.clf()
			pl.loglog( rbin, -4.0*pi*rhobin*vrbin*(parsec*rbin)**2*sec_in_yr/Msun) 
			pl.ylabel('Mdot $M_{dot}/{yr}$')
			pl.xlabel('Radius $pc$')
			pdf_timestep.savefig()
			
			pdf_timestep.close()

	if compare_files == 'nopart':
		print "in nopart"
		rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin = np.loadtxt( "{0}_{1}_{2:04d}_{3}.out".format(file_prefix, quad, framestep, compare_files), unpack=True)
		
		pdf_timestep = mat_pdf.PdfPages("framestep_{0}_{1}_{2:04d}_{3}.pdf".format(file_prefix, quad, framestep, compare_files))
		
		pl.clf()
		pl.loglog( rbin, -vrbin/1e5, 'b') 
		pl.loglog( rbin, vrmsbin/1e5, 'g.-') 
		pl.loglog( rbin, vKbin/1e5, 'r--')
		#pl.loglog( rbin, norm/1e5, 'c+')
		#pl.loglog( rbin, vphi_magbin/1e5, 'k+')
		pl.ylabel('Velocity ($km/s$)')
		pl.xlabel('Radius ($pc$)')
		pl.ylim(1e-2, 10e0)
		pl.xlim(1e-3, 3e0)
		pl.title('{0} framestep {1:04d}'.format(vel_Plot_label, framestep))
		pdf_timestep.savefig()
		
		pl.clf()
		print rhobin[0], rhobin[-1]
		# rhobin should be in g/cm**3
		pl.loglog(rbin, rhobin/mp) #plotting number density
		pl.ylabel('Density (${N}/{cm^3}$)')
		pl.xlabel('Radius ($pc$)')
		pdf_timestep.savefig()
		
		pl.clf()
		pl.loglog( rbin, mTbin)
		pl.ylabel('Total Mass ($g$)')
		pl.xlabel('Radius ($pc$)')
		pdf_timestep.savefig()
		
		pl.clf()
		pl.loglog( rbin, -4.0*pi*rhobin*vrbin*(parsec*rbin)**2*sec_in_yr/Msun)
		pl.ylabel('Mdot $M_{dot}/{yr}$')
		pl.xlabel('Radius $pc$')
		pdf_timestep.savefig()
		
		pdf_timestep.close()
		
