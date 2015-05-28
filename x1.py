import matplotlib
matplotlib.use("Agg")
import sys
import math

sys.path.append('./scripts/joishi/fourier_tools/fourier_tools/')

#import fourier_filter as ff

from yt.mods import * 
import numpy 
import argparse
import shutil
import matplotlib.pyplot as pl
import random
import Image
import os

parser = argparse.ArgumentParser(description = "start")
parser.add_argument('Start', metavar='N1', type=int)
parser.add_argument('End', metavar='N2', type=int)
parser.add_argument('Radius', metavar='N3', type=float)
args = parser.parse_args()
startIndex = args.Start
endIndex = args.End
sizeSphere = args.Radius

bins = 80
meanDens = 3e-22
i = 0
Msun = 1.99e33
G = 6.67e-8
parsec = 3.09e18
window = None
windowInitialized = False


def AccretionProfile(pf, xc, yc, zc, fileout="rad_profile.out", radiusSphere=sizeSphere, deltaR = 0.05, particleMass = 0.0) : 
	sp = pf.h.sphere([xc,yc,zc], (radiusSphere+deltaR)/pf['pc'])
	numPhi = 30
	numTheta = 15

	sp2 = pf.h.sphere([xc,yc,zc], 0.05/pf['pc'])
	com = sp2.quantities["CenterOfMass"]()
	xcom = com[0]
	ycom = com[1]
	zcom = com[2]

#	xc = xcom
#	yc = ycom
#	zc = zcom

	x = sp["x"] - xcom
	y = sp["y"] - ycom
	z = sp["z"] - zcom

	rho = sp["Density"]

	r = numpy.sqrt(x*x + y*y + z*z)

	vx = sp["x-velocity"]
	vy = sp["y-velocity"]
	vz = sp["z-velocity"]

	# now get the radial velocity

	# first get the central velocity
	sp2 = pf.h.sphere([xc,yc,zc], 0.05/pf['pc'])
	vcx = sp2["x-velocity"]
	vcy = sp2["y-velocity"]
	vcz = sp2["z-velocity"]
	cMass = sp2["CellMassMsun"]
	cMassTot = cMass.sum()

	vcx = (cMass*vcx).sum()/cMassTot
	vcy = (cMass*vcy).sum()/cMassTot
	vcz = (cMass*vcz).sum()/cMassTot
	
	# get velocities relative to center
	vx = vx - vcx
	vy = vy - vcy
	vz = vz - vcz

	pi = math.pi

	rhodotSph = numpy.zeros([numTheta, numPhi])
	rhoSph = numpy.zeros([numTheta, numPhi])
	iSph = numpy.zeros([numTheta, numPhi])
	aSph = numpy.zeros([numTheta, numPhi])

	dTheta = pi/numTheta*0.5
	dPhi = 2*pi/numTheta*0.5
	aSph[0,:] = (1.-math.cos(dTheta))*dPhi
	for i in range(1,numTheta) : 
		theta = pi * i/numTheta + dTheta
		theta1 = theta - dTheta
		theta2 = theta + dTheta
		aSph[i,:] = -(math.cos(theta2) - math.cos(theta1))*dPhi

	# get the radial velocity
	for i in range(r.size) : 
		if(abs(r[i] - radiusSphere*parsec)/parsec < deltaR) : 
			theta = math.acos( z[i]/r[i])
			rxy = math.sqrt(x[i]**2 + y[i]**2)
			phi = math.atan2(y[i]/rxy, x[i]/rxy)
			if(rxy < 0.1) : 
				phi = 0.
			rhodot = (vx[i]*x[i] + vy[i]*y[i] + vz[i]*z[i])*rho[i]/r[i]
			indexPhi = int((phi+pi)*numPhi/(2.*pi))
			indexTheta = int((theta)*numTheta/(pi))
			rhodotSph[indexTheta,indexPhi] = rhodotSph[indexTheta,indexPhi] + rhodot
			rhoSph[indexTheta,indexPhi] = rhoSph[indexTheta,indexPhi] + rho[i]
			iSph[indexTheta,indexPhi] = iSph[indexTheta,indexPhi] + 1
	
	rhodotSph = rhodotSph/numpy.maximum(iSph,numpy.ones([numTheta,numPhi]))
	rhoSph = rhoSph/numpy.maximum(iSph,numpy.ones([numTheta,numPhi]))
	mdotSph = rhodotSph * aSph
	numpy.savez( fileout, rhodotSph = rhodotSph, rhoSph=rhoSph, aSph=aSph)



# first get the masses of the most massive particles
prefix = "BB_hdf5_plt_cnt_"

initialized = False

# list of centers 
xcenters = numpy.array([])
ycenters = numpy.array([])
zcenters = numpy.array([])
ncenters = numpy.array([])
indcenters = numpy.array([])

starIndices = None
xstars = None
ystars = None
zstars = None

for fileIndex in range(endIndex,startIndex,-1) : 
	pf = load("{0}{1:04d}".format(prefix,fileIndex))

	partIndices = None
	xps = None
	yps = None
	zps = None
	mps = None
	
	# get all the star particles first 
	if("particle_index" in pf.h.derived_field_list) : 
		all_data = pf.h.all_data()
	
		partIndices = all_data["particle_index"]
		xps = all_data["particle_position_x"]
		yps = all_data["particle_position_y"]
		zps = all_data["particle_position_z"]
		mps = all_data["ParticleMassMsun"]
		print xps	


	if( initialized and partIndices != None) :
		# find missing star particles
		xfound = []
		yfound = []
		zfound = []
		indfound = []
	
		print "looking for star particles", starIndices.size, partIndices
		for i in range(starIndices.size) :
			index = starIndices[i] 
			if not index in partIndices : 	
				# check if something is close enough
				xp = xstars[i]
				yp = ystars[i]
				zp = zstars[i]

				found = False
				for partIndex,x,y,z in zip(partIndices,xps,yps,zps) : 
					r = math.sqrt( (x-xp)*(x-xp) + (y-yp)*(y-yp) + (z-zp)*(z-zp))
					if (r < 0.01*3e18) : 
						found = True
				if( not found) : 
					xfound.append(xp)
					yfound.append(yp)
					zfound.append(zp)
					indfound.append(index)
					print "Found density at ({0:5.2f}, {1:5.2f}, {2:5.2f})".format(xp/3e18, yp/3e18, zp/3e18)
		if( xfound) : # found something now append to list
			xcenters = numpy.append(xcenters,xfound)
			ycenters = numpy.append(ycenters,yfound)
			zcenters = numpy.append(zcenters,zfound)
			ncenters = numpy.append(ncenters,numpy.zeros(len(xfound)))
			indcenters = numpy.append(indcenters, indfound)
		# redefine star indices
		starIndices = partIndices
		xstars = xps
		ystars = yps
		zstars = zps

	elif( partIndices != None) : 
		starIndices = partIndices
		xstars = xps
		ystars = yps
		zstars = zps
		xcenters = numpy.array([])
		ycenters = numpy.array([])
		zcenters = numpy.array([])
		ncenters = numpy.array([])
		initialized = True 
		print "partIndices != None"

	xcenters = xps
	ycenters = yps
	zcenters = zps
	data = pf.h.all_data()
	ncenters = all_data["particle_index"]
	print ncenters
	print xcenters
	print "looking through ncenters" 
	for i in range(ncenters.size) :  
		xc = xcenters[i] 
		yc = ycenters[i] 
		zc = zcenters[i]  
		data = pf.h.sphere([xc,yc,zc], 0.01/pf['pc'])
		max= data.quantities["MaxLocation"]("Density")
		maxDens = max[0]/3e-22
		maxLoc = numpy.array(max[2:5])/3e18
		print fileIndex, maxDens, maxLoc
		xc = max[2]
		yc = max[3]
		zc = max[4]
		
		# redefine the new xcenters
		xcenters[i] = xc
		ycenters[i] = yc
		zcenters[i] = zc
		ncenters[i] = ncenters[i] + 1

		print "Working on {0} at position ({1:5.2f}, {2:5.2f}, {3:5.2f})".format(i,xc/3e18, yc/3e18, zc/3e18)
	
		pc = PlotCollection(pf, center = [xc,yc,zc])
		pc.add_profile_sphere( 5.0, "pc", ["Radiuspc", "Density"], weight="CellVolume",x_bins=150,x_log=False)
		pc.plots[-1].data.write_out("track_density_{0:03d}_{1:02d}.data".format(i,-int(ncenters[i])), format='%15.9E', bin_style='left')
		pc.save()

		AccretionProfile(pf, xc, yc, zc, fileout="Accretion_{0:03d}_{1:03}pc.npz".format(i, sizeSphere), radiusSphere=sizeSphere)

		#getRadialProfile(pf, xc, yc, zc, fileout="TrackDensity_{0:03d}_t={1:02d}.data".format(i,-int(ncenters[i])), radiusSphere=3.0)
	

#for i in range(indcenters.size) : 
#	print "{0} {1}".format(i,int(indcenters[i]))
 
#pf = load("{0}{1:04d}".format(prefix,84))
#for i in range(0,10): 
#	xc = random.random()*4.8e19
#	yc = random.random()*4.8e19
#	zc = random.random()*4.8e19
