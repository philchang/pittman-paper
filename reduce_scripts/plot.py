import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as pl
import numpy as np
import matplotlib.backends.backend_pdf as mat_pdf
import argparse
import glob
import sys
import os

def plotting_routine(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin):
	velocity_plot(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, vphi_magbin)
	density_plot(rbin, mTbin, rhobin)
	ang_moment_plot(rbin, mTbin, vphi_magbin)
	total_mass_plot(rbin, mTbin)
	Mdot_plot(rbin, vrbin, rhobin, mdotbin)
	pdf_close(pdf_timestep)

def velocity_plot(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, vphi_magbin):
	pl.clf()
	pl.loglog( rbin, -vrbin/1e5, 'bv-', label='Infall Velocity')
	pl.loglog( rbin, vrmsbin/1e5, 'g.-', label='Turbulent Velocity')  
	pl.loglog( rbin, vmagbin/1e5, 'm.', label='Total Velocity')  
	pl.loglog( rbin, vKbin/1e5, 'r--', label='Keplerian Velocity')  
	pl.loglog( rbin, vphi_magbin/1e5, 'k+', label='Angular Velocity')
	pl.loglog( rbin, test_all_together/1e5, 'c', label='Sum Velocity')  
	pl.axhline( c_s/1e5, color = 'c',label='cs')  

	pl.legend(loc=0, fontsize='medium', frameon=False)
	pl.ylabel('Velocity ($km/s$)')
	pl.xlabel('Radius ($pc$)')
	pl.ylim(1e-2, 10e0)
	pl.xlim(1e-3, 3e0)
	pl.title('{0} framestep {1:04d}'.format(vel_Plot_label, framestep))
	pdf_timestep.savefig()
	
def density_plot(rbin, mTbin, rhobin):
	pl.clf()
	# rhobin should be in g/cm**3
	pl.loglog(rbin, rhobin/mp, label='Number density') #plotting number density
	pl.loglog(rbin, 2.0/(rbin**3), label='-3')
	pl.loglog(rbin, 80./(rbin**2.5), label='-2.5')
	pl.loglog(rbin, 175./(rbin**2.0), label='-2')
	pl.legend(loc=0, fontsize='medium', frameon=False)
	pl.ylabel('Number Density (${N}/{cm^3}$)')
	pl.xlabel('Radius ($pc$)')
	pdf_timestep.savefig()
	
def ang_moment_plot(rbin, mTbin, vphi_magbin):
	pl.clf()
	pl.loglog( rbin, (rbin*parsec)*vphi_magbin/np.sqrt(G*mTbin*(parsec*rbin)))
	pl.ylabel('Ratio of specific Ang momentum')
	pl.xlabel('Radius ($pc$)')
	pdf_timestep.savefig()
	
def total_mass_plot(rbin, mTbin):
	pl.clf()
	pl.loglog( rbin, mTbin)
	pl.ylabel('Total Mass ($g$)')
	pl.xlabel('Radius ($pc$)')
	pdf_timestep.savefig()
	
def Mdot_plot(rbin, vrbin, rhobin, mdotbin):
	pl.clf()
	pl.loglog( rbin, -4.0*pi*rhobin*vrbin*(parsec*rbin)**2*sec_in_yr/Msun) 
	pl.ylabel('Mdot $M_{odot}/{yr}$')
	pl.xlabel('Radius $pc$')
	pdf_timestep.savefig()
	
def pdf_close(pdf_timestep):
	pdf_timestep.close()

def make_independent_plots(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin):
	print "in make independent plots"
	if (all_independent):
		velocity_alone = True
		density_alone = True
		angmv_alone = True
		mass_alone = True
		mdot_alone = True
		print "made it through conver all to true"
	if (velocity_alone):
		print "plotting just velocity now"
		pdf_timestep = mat_pdf.PdfPages("velocity_{0}_{1:04d}_{2}{3:03d}.pdf".format(quad, framestep, compare_files, particle_number))
		velocity_plot(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, vphi_magbin)
		pdf_close(pdf_timestep)

	if (density_alone):
		pdf_timestep = mat_pdf.PdfPages("density_{0}_{1:04d}_{2}{3:03d}.pdf".format(quad, framestep, compare_files, particle_number))
		density_plot(rbin, mTbin, rhobin)
		pdf_close(pdf_timestep)

	if (angmv_alone):
		pdf_timestep = mat_pdf.PdfPages("angmomentum_{0}_{1:04d}_{2}{3:03d}.pdf".format(quad, framestep, compare_files, particle_number))
		ang_moment_plot(rbin, mTbin, vphi_magbin)
		pdf_close(pdf_timestep)

	if (mass_alone):
		pdf_timestep = mat_pdf.PdfPages("mass_{0}_{1:04d}_{2}{3:03d}.pdf".format(quad, framestep, compare_files, particle_number))
		total_mass_plot(rbin, mTbin)
		pdf_close(pdf_timestep)

	if (mdot_alone):
		pdf_timestep = mat_pdf.PdfPages("mdot_{0}_{1:04d}_{2}{3:03d}.pdf".format(quad, framestep, compare_files, particle_number))
		Mdot_plot(rbin, vrbin, rhobin, mdotbin)
		pdf_close(pdf_timestep)


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
parser.add_argument('--velocity', action='store_true')
parser.add_argument('--density', action='store_true')
parser.add_argument('--mdot', action='store_true')
parser.add_argument('--mass', action='store_true')
parser.add_argument('--angmv', action='store_true')
parser.add_argument('--all', action='store_true')

args = parser.parse_args()

withBigSphere = args.bigsphere
withSmallSphere = args.smallsphere
withParticle = args.particle
withNoParticle = args.noparticle
withShell = args.shell
withShellSphere = args.shellsphere
velocity_alone = args.velocity
density_alone = args.density
angmv_alone = args.angmv
mdot_alone = args.mdot
mass_alone = args.mass
all_independent = args.all


mp = 1.6e-24
pi = np.pi
parsec = 3e18
sec_in_yr = np.pi* 1e7
Msun = 2e33
G = 6.67e-8
c_s = 2.64e4

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
	vel_Plot_label = 'Bulk velocity removed by Sphere in Shell'

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
			particle_file_exist =  glob.glob('{0}_{1}_{2:04d}_{3}{4:03d}.out'.format(file_prefix, quad, framestep, compare_files, particle_number))
			if not particle_file_exist:
				print 'File: "{0}_{1}_{2:04d}_{3}{4:03d}.out" does not exist!'.format(file_prefix, quad, framestep, compare_files, particle_number)
				continue
			# If the file in question exists, we then unpack
			rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin, test_all_together = np.loadtxt( "{0}_{1}_{2:04d}_{3}{4:03d}.out".format(file_prefix, quad, framestep, compare_files, particle_number ), unpack=True)
			# Set Output Filename
			pdf_timestep = mat_pdf.PdfPages("framestep_{0}_{1}_{2:04d}_{3}{4:03d}.pdf".format(file_prefix, quad, framestep, compare_files, particle_number))
			# Call to plotting routine
			plotting_routine(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin)
			#make_independent_plots(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin)
			if (all_independent):
				velocity_alone = True
				density_alone = True
				angmv_alone = True
				mass_alone = True
				mdot_alone = True
			if (velocity_alone):
				pdf_timestep = mat_pdf.PdfPages("velocity_{0}_{1:04d}_{2}{3:03d}.pdf".format(quad, framestep, compare_files, particle_number))
				velocity_plot(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, vphi_magbin)
				pdf_close(pdf_timestep)

			if (density_alone):
				pdf_timestep = mat_pdf.PdfPages("density_{0}_{1:04d}_{2}{3:03d}.pdf".format(quad, framestep, compare_files, particle_number))
				density_plot(rbin, mTbin, rhobin)
				pdf_close(pdf_timestep)
				
			if (angmv_alone):
				pdf_timestep = mat_pdf.PdfPages("angmomentum_{0}_{1:04d}_{2}{3:03d}.pdf".format(quad, framestep, compare_files, particle_number))
				ang_moment_plot(rbin, mTbin, vphi_magbin)
				pdf_close(pdf_timestep)
				
			if (mass_alone):
				pdf_timestep = mat_pdf.PdfPages("mass_{0}_{1:04d}_{2}{3:03d}.pdf".format(quad, framestep, compare_files, particle_number))
				total_mass_plot(rbin, mTbin)
				pdf_close(pdf_timestep)

			if (mdot_alone):
				pdf_timestep = mat_pdf.PdfPages("mdot_{0}_{1:04d}_{2}{3:03d}.pdf".format(quad, framestep, compare_files, particle_number))
				Mdot_plot(rbin, vrbin, rhobin, mdotbin)
				pdf_close(pdf_timestep)

	if compare_files == 'nopart':
		print "in nopart"
		rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin, test_all_together = np.loadtxt( "{0}_{1}_{2:04d}_{3}.out".format(file_prefix, quad, framestep, compare_files), unpack=True)
		# set file name
		pdf_timestep = mat_pdf.PdfPages("framestep_{0}_{1}_{2:04d}_{3}000.pdf".format(file_prefix, quad, framestep, compare_files))
		plotting_routine(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin)

		if (all_independent):
			velocity_alone = True
			density_alone = True
			angmv_alone = True
			mass_alone = True
			mdot_alone = True
		if (velocity_alone):
			pdf_timestep = mat_pdf.PdfPages("velocity_{0}_{1:04d}_{2}000.pdf".format(quad, framestep, compare_files))
			velocity_plot(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, vphi_magbin)
			pdf_close(pdf_timestep)
			
		if (density_alone):
			pdf_timestep = mat_pdf.PdfPages("density_{0}_{1:04d}_{2}000.pdf".format(quad, framestep, compare_files))
			density_plot(rbin, mTbin, rhobin)
			pdf_close(pdf_timestep)
			
		if (angmv_alone):
			pdf_timestep = mat_pdf.PdfPages("angmomentum_{0}_{1:04d}_{2}000.pdf".format(quad, framestep, compare_files))
			ang_moment_plot(rbin, mTbin, vphi_magbin)
			pdf_close(pdf_timestep)

		if (mass_alone):
			pdf_timestep = mat_pdf.PdfPages("mass_{0}_{1:04d}_{2}000.pdf".format(quad, framestep, compare_files))
			total_mass_plot(rbin, mTbin)
			pdf_close(pdf_timestep)

		if (mdot_alone):
			pdf_timestep = mat_pdf.PdfPages("mdot_{0}_{1:04d}_{2}000.pdf".format(quad, framestep, compare_files))
			Mdot_plot(rbin, vrbin, rhobin, mdotbin)
			pdf_close(pdf_timestep)
