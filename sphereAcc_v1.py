import matplotlib
matplotlib.use('Agg') #I tried a few other backends and got the same result
from matplotlib import rc
rc("text", usetex=True)
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

import matplotlib.pyplot as p 
import numpy
from math import pi
import Image
import os
import argparse

quad = os.getcwd()[-5:]
seed = os.getcwd()[-13:-6]
parser = argparse.ArgumentParser(description = "Radius")
parser.add_argument('Radius', metavar='N1', type=float)
args = parser.parse_args()
sizeSphere = args.Radius

r = 0.02*3e18
def plotSphAccretion( filename, fileout, histout) : 
	a = numpy.load(filename)
	rhodotSph = a['rhodotSph']
	rhoSph = a['rhoSph']
	aSph = a['aSph']
	print filename, rhodotSph.size
	img = rhodotSph
	file1 = "filename1_0_{0}.png".format(sizeSphere)
	file2 = "filename2_0_{0}.png".format(sizeSphere)

	p.figure(figsize=(15,10))
	p.axes(projection='mollweide',aspect=2)
	p.imshow(img,extent=(-pi,pi,-pi/2.,pi/2.),aspect=0.5)
	p.colorbar(orientation="horizontal")
	p.savefig(file1)

	img = numpy.log10(rhoSph)
	p.figure(figsize=(15,10))
	p.axes(projection='mollweide',aspect=2)
	p.imshow(img,extent=(-pi,pi,-pi/2.,pi/2.),aspect=0.5)
	p.colorbar(orientation="horizontal")
	p.savefig(file2)

	newImage = Image.new("RGB", (1500,2000))
	image1 = Image.open(file1) 
	image2 = Image.open(file2) 

	newImage.paste( image1, (0,0)) 
	newImage.paste( image2, (0,1000)) 
	newImage.save(fileout) 

	p.clf()
	num_bins = 100
	rhodotSphFlat = -rhodotSph.flatten()
	rhoSphFlat = rhoSph.flatten()
	rhodotSphFlat = numpy.maximum(numpy.zeros(rhodotSphFlat.size), rhodotSphFlat)
	aSphFlat = aSph.flatten()
	mdot = rhodotSphFlat * aSphFlat
	mTot = mdot.sum()
	mdot = mdot/mTot

	# Remove Zeros
	#isort = numpy.nonzero(rhoSphFlat)
	#rhoSort = rhoSphFlat[isort]
	#rhodotSort = rhodotSphFlat[isort]
	#aSphSort = aSphFlat[isort]
	#mdotSort = mdot[isort]
	#rhoSort = rhoSphFlat[isort]

	rhoSort = rhoSphFlat
	rhodotSort = rhodotSphFlat
	aSphSort = aSphFlat
	mdotSort = mdot
	rhoSort = rhoSphFlat

	isort = numpy.argsort( rhodotSort)
	isort = numpy.argsort( rhoSort)
#	isort = isort[::-1] # reverse

	rhodotSort = rhodotSort[isort]
	aSphSort = aSphSort[isort]
	mdotSort = mdotSort[isort]
	rhoSort = rhoSort[isort]
	aTotal= aSphSort.sum()

	mass = rhoSort*aSphSort
	avgRho = mass.sum()/aSphSort.sum()

#	print numpy.log10(rhoSort), numpy.cumsum(mdotSort)	
	p.semilogx( rhoSort, numpy.cumsum(mdotSort)/mdotSort.sum(), linewidth=3,color="black")
	#p.xlabel("Cumulative $\\dot{M}$",fontsize=19)
	#p.ylabel("Cumulative $f_{\\rm sky}$",fontsize=19)
	p.xlabel("$\\rho$",fontsize=60)
	p.ylabel("Cumulative $\\dot{M}$",fontsize=60)
	p.xticks(fontsize=60)
	p.yticks(fontsize=60)
	p.ylim(0., 1.)
	p.xlim(5e-23, 2e-19)
	p.hist(rhodotSphFlat*r*r/(2.99e33/3.15e7/1e5), num_bins, normed=False, cumulative=True, alpha=0.5,weights=aSph.flatten()/aTotal)
	p.xlabel("$\\dot{M}_{-5}\, [{\\rm M_{\\odot}\,yr^{-1}\\,str^{-1}}]$")
	p.ylabel("$f_{\\rm sky}$") 
	p.ylim(0., 1.)
	p.xlim(1e-19, 2e-17)
	p.savefig( histout)
	return numpy.cumsum(mdotSort)/mdotSort.sum(), numpy.cumsum(aSphSort/aTotal), rhoSort/avgRho

fileOutPrefix = "fileout"
histOutPrefix = "hist"
partIds = [0,1]
#partIds = [0,1,3]
#partIds = [2,4]
times = [0]

logRho = numpy.arange( -1.5, 1, 0.1)

cmdotAvg = None
caSphAvg = None
logRhoAvg = None
rhoNormAvg = None

# get the first one
i = 0
t = 0
time = times[t]

for i in range( len(partIds)) :
	filein = "Accretion_{0:03d}_{1}pc.npz".format(partIds[i],sizeSphere)
	fileIndex = i + t*len(partIds)
	fileout = "{0}_{1:02d}_{2}.png".format(fileOutPrefix, partIds[i],sizeSphere)
	histout = "{0}_{1:02d}_{2}.pdf".format(histOutPrefix, partIds[i],sizeSphere)
	cmdotSort, caSphSort, rhoNorm = plotSphAccretion( filein, fileout, histout)
	print rhoNorm.size, cmdotSort.size
	if( rhoNormAvg == None) : 
		cmdotAvg = cmdotSort
		rhoNormAvg = rhoNorm
	else :
		cmdotAvg = cmdotAvg + cmdotSort
		rhoNormAvg = rhoNormAvg + rhoNorm

	p.clf()

rhoNormAvg = rhoNormAvg/len(partIds)
cmdotAvg = cmdotAvg/len(partIds)

p.semilogx( rhoNormAvg, cmdotAvg, label="First star $t_6=0.1$", lw=3, ls="-.") 

#p.xlabel("Cumulative average normalized $\\dot{M}$",fontsize=19)
#p.ylabel("Cumulative average normalized $f_{\\rm sky}$",fontsize=19)
p.xlabel("Cumulative average normalized $\\rho/\\bar{\\rho}$",fontsize=40)
p.ylabel("Cumulative average normalized $\\dot{M}$",fontsize=40)
p.ylim(0., 1.)
#p.xlim(0., 1.)
p.xlim(0.03, 30.)
p.xticks(fontsize=45)
p.yticks(fontsize=45)


p.legend(loc="best",fontsize=40)
p.savefig( "Average_{0}_{1}_{2}.pdf".format(seed,quad,sizeSphere))

#for i in range(len(fileout)) : 
#	os.remove(fileout[i])
