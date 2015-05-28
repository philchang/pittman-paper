import matplotlib
matplotlib.use("Agg")

from yt.mods import *
import matplotlib.pyplot as p

import math
import numpy
import random
import time
import datetime
import os

bins = 70
meanDens = 3e-22
i = 0
Msun = 1.99e33
G = 6.67e-8
parsec = 3.09e18
logRhoMin = -3.0

def getRadialProfile_yt(pf, xc, yc, zc, vxc, vyc, vzc, fileout="rad_profile.out", radiusMin=1e-3, radiusSphere=3.0, particleMass = 0.0) : 
	if (Sphere_Bulk):
		sp = pf.h.sphere((xc, yc, zc), radiusSphere/pf['pc'])
		sp2 = pf.h.sphere((xc, yc, zc), Bulk_sphere_radius/pf['pc'])
		print Bulk_sphere_radius
		bulk = sp2.quantities["BulkVelocity"]()
	if (Particle_Bulk):
		sp = pf.h.sphere((xc, yc, zc), radiusSphere/pf['pc'])
		#sp2 = pf.h.sphere((xc, yc, zc), 0.01/pf['pc'])
		bulk = numpy.array([vxc, vyc, vzc])
	if (NoParticle_Bulk):
		sp = pf.h.sphere((xc, yc, zc), 4.0/pf['pc'])
		bulk = sp.quantities["BulkVelocity"]()
	sp.set_field_parameter("bulk_velocity", bulk)
	rp = BinnedProfile1D( sp, bins, "Radiuspc", 1e1**logRhoMin, radiusSphere, log_space=True)
	rp.add_fields("Density", weight="CellVolume")
	maxDens = rp["Density"][0]
	
	rp.add_fields("RadialVelocity")  
	rp.add_fields("VelocityMagnitude")
	
	rp.add_fields("CellMassMsun", accumulation=True,weight=None)

	rbin = rp["Radiuspc"]
	vrbin = rp["RadialVelocity"]
	mTbin = rp["CellMassMsun"] + particleMass
	vKbin = numpy.sqrt(G*mTbin*Msun/(parsec*rbin))
	rhobin = rp["Density"]
	mdotbin = 4.0*3.141*rbin**2*vrbin*rhobin
	rp.add_fields("TangentialVelocity")
	vrmsbin = rp["TangentialVelocity"]
	rp.add_fields("TangentialVelocity", weight="CellVolume")
	vrmsnbin = rp["TangentialVelocity"]
	rp.add_fields("VelocityMagnitude")
	vmagbin = rp["VelocityMagnitude"]	
	rp.add_fields("VelocityMagnitude", weight="CellVolume")
	vmagnbin = rp["VelocityMagnitude"]

	rp.add_fields(["AngularMomentumX", "AngularMomentumY", "AngularMomentumZ"])

	norm = numpy.sqrt(rp["AngularMomentumX"]**2 + rp["AngularMomentumY"]**2 + rp["AngularMomentumZ"]**2)
	angXbin = rp["AngularMomentumX"]/norm
	angYbin = rp["AngularMomentumY"]/norm
	angZbin = rp["AngularMomentumZ"]/norm
	
	numpy.savetxt(fileout, zip(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin), fmt="%15.9E")

def getRadialProfile_py(pf, xc, yc, zc, fileout="rad_profile.out", radiusMin=1e-3, radiusSphere=3.0, particleMass = 0.0) : 
	sp = pf.h.sphere((xc, yc, zc), radiusSphere/pf['pc'])

	x = sp["x"] - xc
	y = sp["y"] - yc
	z = sp["z"] - zc
	r = numpy.sqrt(x*x + y*y + z*z)
	# grams
	cellMass = sp["CellMassMsun"] * Msun
	# cgs
	dens = sp["Density"]
	# cm**3
	cellVolume = sp["CellVolume"]

	# These are in cm/s
	vx = sp["x-velocity"]
	vy = sp["y-velocity"]
	vz = sp["z-velocity"]

	# This is the total Mass
	Mtot = cellMass.sum()

	# loosely define the ang velocity
	# give proper size of lx
	# Its before bulk velocity subtraction
	lx = y*vz - vy*z
	ly = z*vx - x*vz
	lz = x*vy - y*vx

############################################################################
	if (ShellSphere_Bulk):
		# Find the Bulk Velocity of each shell, prior to removing from each cell in the shell
		print "Calculating the Bulk Velocity by sphere inside shell"
		ts = time.time()
		st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
		print st

		mbin = numpy.zeros(bins)
		menc = numpy.zeros(bins)
		vx_sphere_bulk_bin = numpy.zeros(bins)
		vy_sphere_bulk_bin = numpy.zeros(bins)
		vz_sphere_bulk_bin = numpy.zeros(bins)
		vx_bulkvelocity_bin = numpy.zeros(bins)
		vy_bulkvelocity_bin = numpy.zeros(bins)
		vz_bulkvelocity_bin = numpy.zeros(bins)
		
		lgradiusMin = math.log10( radiusMin)
		lgradiusSph = math.log10( radiusSphere)
		for i in range(r.size) : 
			if( r[i]/parsec < radiusMin) :
				continue
			index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
			if(index >= 0 and index < bins) :
				# Get the mass of each shell
				mbin[index] = mbin[index] + cellMass[i]
				
				# Find the Bulk velocity (Momentum) of each shell
				vx_sphere_bulk_bin[index] = vx_sphere_bulk_bin[index] + (vx[i]*cellMass[i])
				vy_sphere_bulk_bin[index] = vy_sphere_bulk_bin[index] + (vy[i]*cellMass[i])
				vz_sphere_bulk_bin[index] = vz_sphere_bulk_bin[index] + (vz[i]*cellMass[i])
		menc[0] = mbin[0]
		
		for shell in range(1,bins):
			menc[shell] = mbin[shell] + menc[shell-1]
			vx_sphere_bulk_bin[shell] = vx_sphere_bulk_bin[shell] + vx_sphere_bulk_bin[shell-1] 
			vy_sphere_bulk_bin[shell] = vy_sphere_bulk_bin[shell] + vy_sphere_bulk_bin[shell-1] 
			vz_sphere_bulk_bin[shell] = vz_sphere_bulk_bin[shell] + vz_sphere_bulk_bin[shell-1] 

		# Set the bulk velocity to be this bulk velocity		
		vx_bulkvelocity_bin = vx_sphere_bulk_bin/menc
		vy_bulkvelocity_bin = vy_sphere_bulk_bin/menc
		vz_bulkvelocity_bin = vz_sphere_bulk_bin/menc

		#vx_bulkvelocity_bin[:] = vx_bulkvelocity_bin[-1]
		#vy_bulkvelocity_bin[:] = vy_bulkvelocity_bin[-1]
		#vz_bulkvelocity_bin[:] = vz_bulkvelocity_bin[-1]

#		print vx_bulkvelocity_bin[-1], vy_bulkvelocity_bin[-1], vz_bulkvelocity_bin[-1]
#		sp2 = pf.h.sphere((xc, yc, zc), 3.0/pf['pc'])
#		bulk = sp2.quantities["BulkVelocity"]()
#		print bulk

###########################################################################
###########################################################################
	if (Shell_Bulk) or (NoParticle_Bulk):
		# Find the Bulk Velocity of each shell, prior to removing from each cell in the shell
		print 'Finding Bulk Velocity by Shell'
		ts = time.time()
		st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
		print st

		mbin = numpy.zeros(bins)
		vx_shell_bulk_bin = numpy.zeros(bins)
		vy_shell_bulk_bin = numpy.zeros(bins)
		vz_shell_bulk_bin = numpy.zeros(bins)
		vx_bulkvelocity_bin = numpy.zeros(bins)
		vy_bulkvelocity_bin = numpy.zeros(bins)
		vz_bulkvelocity_bin = numpy.zeros(bins)
		
		lgradiusMin = math.log10( radiusMin)
		lgradiusSph = math.log10( radiusSphere)
		for i in range(r.size) :
			if( r[i]/parsec < radiusMin) :
				continue
			index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
			if(index >= 0 and index < bins) :
				# Get the mass of each shell
				mbin[index] = mbin[index] + cellMass[i]
				# Find the Bulk Velocity of each shell
				vx_shell_bulk_bin[index] = vx_shell_bulk_bin[index] + (vx[i]*cellMass[i])
				vy_shell_bulk_bin[index] = vy_shell_bulk_bin[index] + (vy[i]*cellMass[i])
				vz_shell_bulk_bin[index] = vz_shell_bulk_bin[index] + (vz[i]*cellMass[i])


       		vx_shell_bulk_bin = vx_shell_bulk_bin/mbin
	       	vy_shell_bulk_bin = vy_shell_bulk_bin/mbin
		vz_shell_bulk_bin = vz_shell_bulk_bin/mbin

		# Set the bulk velocity to be this bulk velocity
		vx_bulkvelocity_bin = vx_shell_bulk_bin
		vy_bulkvelocity_bin = vy_shell_bulk_bin
		vz_bulkvelocity_bin = vz_shell_bulk_bin


############################################################################
	print 'Obtaining vr'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st

	#mbin = numpy.zeros(bins) # mbin calculated above
	rhobin = numpy.zeros(bins)
	volbin = numpy.zeros(bins)
	vrbin = numpy.zeros(bins)
	mdotbin = numpy.zeros(bins)

	vxbin = numpy.zeros(bins)
	vybin = numpy.zeros(bins)
	vzbin = numpy.zeros(bins)
	vr_nomass = numpy.zeros(vx.size)
	# get the radial velocity
	lgradiusMin = math.log10( radiusMin)
	lgradiusSph = math.log10( radiusSphere)
	for i in range(r.size) : 
		if( r[i]/parsec < radiusMin) :
			continue
		index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
		if(index >= 0 and index < bins) :
			#mbin[index] = mbin[index] + cellMass[i]
			# This is a volume weighted density. i.e. Calculate the mass
			# We'll then divide the mass by volbin
			rhobin[index] = rhobin[index] + cellMass[i]
			volbin[index] = volbin[index] + cellVolume[i]
			# getting a mass weighted speed out of a velocity (Mdot * distance). 
			#vr = ((vx[i] - vx_bulkvelocity_bin[index])*x[i] + (vy[i] - vy_bulkvelocity_bin[index])*y[i] + (vz[i] - vz_bulkvelocity_bin[index])*z[i])*cellMass[i]/r[i]

			# reset vx, vy, vz to remove bulk velocity
			vx[i] = vx[i] - vx_bulkvelocity_bin[index]
			vy[i] = vy[i] - vy_bulkvelocity_bin[index]
			vz[i] = vz[i] - vz_bulkvelocity_bin[index]
			vr_nomass[i] = (vx[i]*x[i] + vy[i]*y[i] + vz[i]*z[i])/r[i]
			vr = vr_nomass[i] * cellMass[i]
			#mdotbin[index] = mdotbin[index] + vr # vr is mdot right now
			mdotbin[index] = mdotbin[index] + vr/r[i]
			vrbin[index] = vrbin[index] + vr
		
	vrbin = vrbin/mbin
	# Check to see if these come out the same:
	#rhobin = mbin/volbin
	rhobin = rhobin/volbin

	# Find the middle radius of the bin
	lgrbin = lgradiusMin + (lgradiusSph-lgradiusMin)*numpy.arange(bins)/bins
	rbin = 1e1**lgrbin
	rbinparsec = rbin * parsec

#####################################################################
	print 'Obtaining the angular momentum'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st

	mbin = numpy.zeros(bins) # mbin calculated above
	menc = numpy.zeros(bins)
	angXbin = numpy.zeros(bins)
	angYbin = numpy.zeros(bins)
	angZbin = numpy.zeros(bins)
	vphi_magbin = numpy.zeros(bins)

	for i in range(r.size) : 
		if( r[i]/parsec < radiusMin) :
			continue
		index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
		if(index >= 0 and index < bins) :
			mbin[index] = mbin[index] + cellMass[i]
			# Now calculate the angular momentum (technically just r x v here)
			lx[i] = y[i]*vz[i] - vy[i]*z[i]
			ly[i] = z[i]*vx[i] - x[i]*vz[i]
			lz[i] = x[i]*vy[i] - y[i]*vx[i]
			
			# Don't need this as we have reset vx vy and vz to include this
			#lx[i] = y[i]*(vz[i] - vz_bulkvelocity_bin[index]) - (vy[i] - vy_bulkvelocity_bin[index])*z[i]
			#ly[i] = z[i]*(vx[i] - vx_bulkvelocity_bin[index]) - (vz[i] - vz_bulkvelocity_bin[index])*x[i]
			#lz[i] = x[i]*(vy[i] - vy_bulkvelocity_bin[index]) - (vx[i] - vx_bulkvelocity_bin[index])*y[i]
			# vphi is perp to this ang momentum vector
			angXbin[index] = lx[i] * cellMass[i] + angXbin[index]
			angYbin[index] = ly[i] * cellMass[i] + angYbin[index]
			angZbin[index] = lz[i] * cellMass[i] + angZbin[index]


	# This is the mass weighted ang momentum per bin not by sphere
	angXbin = angXbin/mbin
	angYbin = angYbin/mbin
	angZbin = angZbin/mbin
	L = angZbin
	L = numpy.sqrt(angXbin**2 + angYbin**2 + angZbin**2)
	norm = numpy.sqrt(angXbin**2 + angYbin**2 + angZbin**2)
	vphi_magbin = numpy.sqrt(angXbin**2 + angYbin**2 + angZbin**2)
	vphi_magbin = vphi_magbin / rbinparsec

## This next section calculates the avg angular momentum per shell
# Not the actual ang momentum
#
#	angX_sphere = angXbin
#	angY_sphere = angYbin
#	angZ_sphere = angZbin
#
#	menc[0] = mbin[0]
#	for shell in range(1,bins):
#		menc[shell] = mbin[shell] + menc[shell-1]
#		angX_sphere[shell] = angX_sphere[shell] + angX_sphere[shell - 1]
#		angY_sphere[shell] = angY_sphere[shell] + angY_sphere[shell - 1]
#		angZ_sphere[shell] = angZ_sphere[shell] + angZ_sphere[shell - 1]
#	
#	lx_sphere = angX_sphere/menc
#	ly_sphere = angY_sphere/menc
#	lz_sphere = angZ_sphere/menc
#	norm_avg = numpy.sqrt(lx_sphere**2 + ly_sphere**2 + lz_sphere**2)
#	# This is lmagbin (r*v)
#	vphi_magbin_avg = numpy.sqrt(lx_sphere**2 + ly_sphere**2 + lz_sphere**2)
#	# This is v = l/r
#	vphi_magbin_avg = vphi_magbin / rbinparsec

#####################################################################
	print 'Obtaining the rms velocity'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st


	# get the rms velocity
	mbin = numpy.zeros(bins) # Reset mbin as its used in calculating vkep
	vrmsbin = numpy.zeros(bins)
	vrmsnbin = numpy.zeros(bins)
	vmagbin = numpy.zeros(bins)
	vmagnbin = numpy.zeros(bins)
	nbin = numpy.zeros(bins)

#	angXbin = numpy.zeros(bins)
#	angYbin = numpy.zeros(bins)
#	angZbin = numpy.zeros(bins)
#
	for i in range(r.size) : 
		if( r[i]/parsec < radiusMin) :
			continue
		index = int((math.log10(r[i]/parsec)-lgradiusMin)*bins/(lgradiusSph - lgradiusMin))
		#index = int(r[i]*bins/(radiusSphere*parsec))
		if(index >= 0 and index < bins) :
			mbin[index] = mbin[index] + cellMass[i]
			nbin[index] = nbin[index] + 1
			rad = r[i]
			if( rad < 1e-5) : 
				 rad = 1e-5
			# This is the correct turbulent Velocity
			vrmsn = math.sqrt((vx[i]-vrbin[index]*x[i]/rad)**2. + (vy[i]-vrbin[index]*y[i]/rad)**2. + (vz[i]-vrbin[index]*z[i]/rad)**2.)
			vrms = vrmsn*cellMass[i]
			vmagn = math.sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])
			vmag = vmagn*cellMass[i]
			# Calculating this here to try to match up to what yt does
			# As yt gets the tangential velocity, and not the turbulent velocity.
			#vr = ((vx[i] - vx_bulkvelocity_bin[index])*x[i] + (vy[i] - vy_bulkvelocity_bin[index])*y[i] + (vz[i] - vz_bulkvelocity_bin[index])*z[i])/r[i]

			#vrmsn = math.sqrt( vmagn*vmagn - vr*vr)
			#vrms = vrmsn*cellMass[i]
			# End test to compare to yt.
			# Why do we set vr here?
			#vr = vr_nomass[i]
			vrmsbin[index] = vrmsbin[index] + vrms
			vrmsnbin[index] = vrmsnbin[index] + vrmsn
			vmagbin[index] = vmagbin[index] + vmag
			vmagnbin[index] = vmagnbin[index] + vmagn

#			lx[i] = y[i]*vz[i] - vy[i]*z[i]
#			ly[i] = z[i]*vx[i] - x[i]*vz[i]
#			lz[i] = x[i]*vy[i] - y[i]*vx[i]
#
#			# Don't need this as we have reset vx vy and vz to include this
#			#lx[i] = y[i]*(vz[i] - vz_bulkvelocity_bin[index]) - (vy[i] - vy_bulkvelocity_bin[index])*z[i]
#			#ly[i] = z[i]*(vx[i] - vx_bulkvelocity_bin[index]) - (vz[i] - vz_bulkvelocity_bin[index])*x[i]
#			#lz[i] = x[i]*(vy[i] - vy_bulkvelocity_bin[index]) - (vx[i] - vx_bulkvelocity_bin[index])*y[i]
#
#			angXbin[index] = lx[i] * cellMass[i] + angXbin[index]
#			angYbin[index] = ly[i] * cellMass[i] + angYbin[index]
#			angZbin[index] = lz[i] * cellMass[i] + angZbin[index]

	vrmsbin = vrmsbin/mbin
	vrmsnbin = vrmsnbin/nbin
	vmagbin = vmagbin/mbin
	vmagnbin = vmagnbin/nbin

###############################################################################

	# get the Kepler velocity
	print 'Obtaining the Kepler Velocity'
	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
	print st
	menc[0] = mbin[0]
	for shell in range(1,bins):
		menc[shell] = mbin[shell] + menc[shell-1]
	mTbin = menc
	# mTbin is in grams
	mTbin = mTbin + particleMass * Msun # Particle Mass was in solar masses
	# calculated rbin in the ang momentum section
	# and converted from being in parsecs
	#lgrbin = lgradiusMin + (lgradiusSph-lgradiusMin)*numpy.arange(bins)/bins
	#rbin = 1e1**lgrbin
	#vKbin = numpy.sqrt( G*mTbin*Msun/((rbin)*parsec))
	#vKbin = numpy.sqrt( G*mTbin/((rbin)*parsec))
	vKbin = numpy.sqrt( G*mTbin/rbinparsec)

	print vphi_magbin[0], vphi_magbin[1]
	print vKbin[0], vKbin[1]

	numpy.savetxt(fileout, zip(rbin, vrbin, vrmsbin, vrmsnbin, vKbin, vmagbin, vmagnbin, mTbin, rhobin, mdotbin, norm, angXbin, angYbin, angZbin, vphi_magbin), fmt="%15.9E")

################################################

def getDensityExclude( pf, xc, yc, zc, xexclude, yexclude, zexclude, fileoutrho_PDF="rho_PDF.out", radiusExclude=1.0, exclude=True) : 
	sp = pf.h.all_data()

	dens = sp["Density"]
	cellMass = sp["CellMassMsun"]
	x = sp["x"]
	y = sp["y"]
	z = sp["z"]
	
	# create a pdf
	volWeightedPDF = numpy.zeros(bins)
	massWeightedPDF = numpy.zeros(bins)
	volWeightedIPDF = numpy.zeros(bins)
	massWeightedIPDF = numpy.zeros(bins)

	rE2 = (radiusExclude*parsec)**2.

	for i in range(dens.size) :
		includeDens = True
		excludeDens = False
		if( exclude ) :
			for j in range(xexclude.size) :
				xpos = x[i] - xexclude[j]
				ypos = y[i] - yexclude[j]
				zpos = z[i] - zexclude[j]
				r2 = xpos*xpos + ypos*ypos + zpos*zpos
				if( r2 < rE2) :
					includeDens = False
					excludeDens = True
					break
		if( includeDens) : 	
			rho = dens[i]
			mcell = cellMass[i]
			lRho = math.log10(rho/meanDens)
			index = int((lRho - logRhoMin)*bins/(logRhoMax-logRhoMin))
			if(index >= 0 and index < bins) :
				volWeightedPDF[index] = volWeightedPDF[index] + 1
				massWeightedPDF[index] = massWeightedPDF[index] + mcell
		if( excludeDens) : 	
			rho = dens[i]
			mcell = cellMass[i]
			lRho = math.log10(rho/meanDens)
			index = int((lRho - logRhoMin)*bins/(logRhoMax-logRhoMin))
			if(index >= 0 and index < bins) :
				volWeightedIPDF[index] = volWeightedIPDF[index] + 1
				massWeightedIPDF[index] = massWeightedIPDF[index] + mcell
	# normalized 
	volWeightedPDF = volWeightedPDF/volWeightedPDF.sum()
	massWeightedPDF = massWeightedPDF/massWeightedPDF.sum()
	volWeightedIPDF = volWeightedIPDF/volWeightedIPDF.sum()
	massWeightedIPDF = massWeightedIPDF/massWeightedIPDF.sum()

	# write out
	lRhoArr = numpy.arange(bins)*(logRhoMax - logRhoMin)/bins + logRhoMin
	numpy.savetxt(fileoutrho_PDF, zip( lRhoArr, volWeightedPDF, massWeightedPDF,volWeightedIPDF, massWeightedIPDF), fmt="%15.9E")

def getDensityFunctions( pf, xc, yc, zc, fileoutrho_r="rho_r.out", fileoutrho_PDF="rho_PDF.out", radiusSphere=3.0, allData = False) : 
	sp = 0
	if( allData ) : 
		sp = pf.h.all_data()
	else :
		sp = pf.h.sphere([xc,yc,zc], radiusSphere/pf['pc'])

	dens = sp["Density"]
	cellMass = sp["CellMassMsun"]
	
	# create a pdf
	volWeightedPDF = numpy.zeros(bins)
	massWeightedPDF = numpy.zeros(bins)
	
	for i in range(dens.size) :
		rho = dens[i]
		mcell = cellMass[i]
		lRho = math.log10(rho/meanDens)
		index = int((lRho - logRhoMin)*bins/(logRhoMax-logRhoMin))
		if(index >= 0 and index < bins) :
			volWeightedPDF[index] = volWeightedPDF[index] + 1
			massWeightedPDF[index] = massWeightedPDF[index] + mcell
	
	# normalized 
	volWeightedPDF = volWeightedPDF/volWeightedPDF.sum()
	massWeightedPDF = massWeightedPDF/massWeightedPDF.sum()

	# write out
	lRhoArr = numpy.arange(bins)*(logRhoMax - logRhoMin)/bins + logRhoMin
	numpy.savetxt(fileoutrho_PDF, zip( lRhoArr, volWeightedPDF, massWeightedPDF), fmt="%15.9E")

	if( allData) :
		return

	# calculate rho as function of r
	x = sp["x"] - xc
	y = sp["y"] - yc
	z = sp["z"] - zc
	r = numpy.sqrt(x*x + y*y + z*z)

	volWeightedRho_r = numpy.zeros(bins)
	massWeightedRho_r = numpy.zeros(bins)
	vol_r = numpy.zeros(bins)
	mass_r = numpy.zeros(bins)
	
	for i in range(dens.size) :
		rho = dens[i]
		mcell = cellMass[i]
		dist = r[i] 
		index = int(bins*dist/(radiusSphere*parsec))
		if( index >= 0 and index < bins) : 
			volWeightedRho_r[index] = volWeightedRho_r[index] + rho
			massWeightedRho_r[index] = massWeightedRho_r[index] + rho*mcell
			vol_r[index] = vol_r[index] + 1
			mass_r[index] = mass_r[index] + mcell
	
	volWeightedRho_r = volWeightedRho_r/vol_r
	massWeightedRho_r = massWeightedRho_r/mass_r

	# write out
	rbin = radiusSphere*numpy.arange(bins)/bins
	numpy.savetxt(fileoutrho_r, zip( rbin, volWeightedRho_r, massWeightedRho_r), fmt="%15.9E")

def getLarsonsLaw(pf, xc, yc, zc, fileout="larson.out", radiusSphere=3.0,trials=10,randPoints = True) : 
	sp = pf.h.sphere([xc,yc,zc], radiusSphere/pf['pc'])

	x = sp["x"] - xc
	y = sp["y"] - yc
	z = sp["z"] - zc

	vx = sp["x-velocity"]
	vy = sp["y-velocity"]
	vz = sp["z-velocity"]
	vrmsbin = numpy.zeros(bins)
	nbin = numpy.zeros(bins)

	if ( not randPoints) : 
		trials = 1

	for trial in range(trials) : 
		ridx = 0
		if(not randPoints) :
			for i in range(x.size) :
				r2 = (x[i]*x[i] + y[i]*y[i] + z[i]*z[i])/(parsec*parsec)
				if ( r2 < 0.03*0.03) :
					ridx = i 
				
		else : 
			ridx = random.randint(0,x.size-1) 

		vxr = vx[ridx]
		vyr = vy[ridx]
		vzr = vz[ridx]
		xr = x[ridx]
		yr = y[ridx]
		zr = z[ridx]

		for i in range(x.size) : 
			if(i != ridx) : 
				dist = math.sqrt((x[i] - xr)**2. + (y[i] - yr)**2 + (z[i] - zr)**2)
				vrms = (vx[i] - vxr)**2. + (vy[i] - vyr)**2 + (vz[i] - vzr)**2

				index = int(bins*dist/(2.*radiusSphere*parsec))

				if(index < bins) :
					vrmsbin[index] = vrmsbin[index] + vrms
					nbin[index] = nbin[index] + 1
		

	vrmsbin = numpy.sqrt(vrmsbin/nbin)
	rbin = 2.*radiusSphere*numpy.arange(bins)/bins

	numpy.savetxt(fileout, zip(rbin, vrmsbin), fmt="%12.7e")



#pf = load("BB_hdf5_plt_cnt_0096")

# compressive case

import argparse
parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N2', type=int)
parser.add_argument('--smallsphere', action='store_true')
parser.add_argument('--bigsphere', action='store_true')
parser.add_argument('--particle', action='store_true')
parser.add_argument('--noparticle', action='store_true')
parser.add_argument('--shell', action='store_true')
parser.add_argument('--shellsphere', action='store_true')
parser.add_argument('--chk', action='store_true')

args = parser.parse_args()

YT_CALL = False
PY_CALL = False
# These settings eventually go to yt calls
withSmallSphere = args.smallsphere
withBigSphere = args.bigsphere
withParticle = args.particle
withNoParticle = args.noparticle
withCheckpoint = args.chk

Sphere_Bulk = False
Particle_Bulk = False
NoParticle_Bulk = False

if (withSmallSphere):
	compare_file = 'smallsphere'
	Sphere_Bulk = True
	Bulk_sphere_radius = 0.01
	YT_CALL = True

if (withBigSphere):
	compare_file = 'bigsphere'
	Sphere_Bulk = True
	Bulk_sphere_radius = 3.0
	YT_CALL = True

if (withParticle):
	compare_file = 'part'
	Particle_Bulk = True
	YT_CALL = True

if (withNoParticle):
	compare_file = 'nopart'
	NoParticle_Bulk = True
	YT_CALL = True

# These calls eventually go to the python written outputs
withShell = args.shell
withShellSphere = args.shellsphere

Shell_Bulk = False
ShellSphere_Bulk = False

if (withShell):
	compare_file = 'shell'
	Shell_Bulk = True
	PY_CALL = True


if (withShellSphere):
	compare_file = 'shellsphere'
	ShellSphere_Bulk = True
	PY_CALL = True

# These are universal
if (withCheckpoint):
	plt_prefix = "BB_hdf5_chk_"
else:
	plt_prefix = "BB_hdf5_plt_cnt_"

out_prefix = "rad_profile"
quad = os.getcwd()[-5:]
# The file looks like this:
#'{out_prefix}{framestep}_{compare_file}_{particle_number}.out'                                                                                               
# i.e. rad_profile_0218_part_000.out   
# compare file options area: bigsphere, smallsphere, part, nopart, shell, shellsphere
for i in range(args.start,args.end,args.step) :	
	print 'On File: {0}_{1:04d}'.format(plt_prefix, i)
	if not (withNoParticle):
		pf = load("{0}{1:04d}".format( plt_prefix, i))

		dd = pf.h.all_data()
		xp = dd["particle_posx"]
		yp = dd["particle_posy"]
		zp = dd["particle_posz"]
		vxp = dd["particle_velocity_x"]
		vyp = dd["particle_velocity_y"]
		vzp = dd["particle_velocity_z"]
		partMass = dd["ParticleMassMsun"]

		for j in range(xp.size) :
			xc = xp[j]	 
			yc = yp[j]	 
			zc = zp[j]
			
			vxc = vxp[j]
			vyc = vyp[j]
			vzc = vzp[j]	 
			if (withParticle) or (withSmallSphere) or (withBigSphere):
				print 'Using yt'
				print 'On particle:', j + 1, 'of:', xp.size
				getRadialProfile_yt(pf,xc,yc,zc, vxc, vyc, vzc, fileout="{0}_{1}_{2:04d}_{3}_{4:03d}.out".format( out_prefix, quad, i, compare_file, j), radiusSphere=3.0, particleMass = partMass[j])

			if (withShell) or (withShellSphere):
				print 'going to python script'
				print 'On particle:', j + 1, 'of:', xp.size
				getRadialProfile_py(pf, xc, yc, zc, fileout="{0}_{1}_{2:04d}_{3}_{4:03d}.out".format( out_prefix, quad, i, compare_file, j), radiusMin=1e-3, radiusSphere=3.0, particleMass = partMass[j])

	if (withNoParticle):
		# Use this for before star particle formation.
		# Make sure getRadialProfile_yt has bulk set to sphere and not particle Velocity
		pf = load("{0}{1:04d}".format( plt_prefix, i))
		dd = pf.h.all_data()
		max= dd.quantities["MaxLocation"]("Density")
		maxDens = max[0]/3e-22
		maxLoc = numpy.array(max[2:5])/3e18
		xc = max[2]
		yc = max[3]
		zc = max[4]
		vxc = 0
		vyc = 0
		vzc = 0
		getRadialProfile_py(pf, xc, yc, zc, fileout="{0}_{1}_{2:04d}_{3}.out".format( out_prefix, quad, i, compare_file), radiusMin=1e-3, radiusSphere=3.0, particleMass = 0.0)
