import shutil
import argparse
import os
from yt.mods import *
import matplotlib.colorbar as cb
from PIL import Image

#def project( filein, xc, yc, zc, zoom_width, fileout = 'junk.png', withParticles=False) :
#	if(withProjectionFull):
#		print "Doing a projection plot of the Full Box."
#		p = ProjectionPlot(pf, "z", "Density")
#	if not (withProjectionFull):
#		p = ProjectionPlot(pf, "z", "Density", center = (xc, yc, zc), width = (zoom_width, 'pc'))
#        p.set_zlim("Density", 1e-3,1e0)
#	if(withParticles) :
#		p.annotate_particles(1.0/pf['unitary'], p_size=10.0)
#        pid = os.getpid()
#        p.save("junk{0:06d}".format(pid))
#        shutil.move("junk{0:06d}_Projection_z_Density.png".format(pid),fileout)
#
def Proj_plotting(pf, xc, yc, zc, zoom_width, quad, i, j):
        plot_out_prefix = 'movieframe'
	if(withProjectionFull):
		print "Doing a projection plot of the Full Box."
		p = ProjectionPlot(pf, "z", "Density")
		fileout="{0}_{1}_{2:04d}_fullprojection.png".format(plot_out_prefix, quad, i)
		if(withParticles) :
			p.annotate_particles(1.0/pf['unitary'], p_size=10.0)
			fileout="{0}_{1}_{2:04d}_fullproject_particle.png".format(plot_out_prefix, quad, i)
		pid = os.getpid()
		p.save("junk{0:06d}".format(pid))
		shutil.move("junk{0:06d}_Projection_{1}_Density.png".format(pid, plot_axis),fileout)
	if not (withProjectionFull):
		plot_field = plot_axis + '-velocity'
		fileout="{0}_{1}_{2:04d}_{3}_{4}_{5:03d}_{6}pc.png".format(plot_out_prefix, quad, i, plot_axis, plot_field, j, zoom_width)
		p_plot = ProjectionPlot(pf, plot_axis, plot_field, weight_field = 'Density', center = (xc, yc, zc), width = (zoom_width, 'pc'))
		pid = os.getpid()
		p_plot.save("junk{0:06d}".format(pid))
		shutil.move("junk{0:06d}_Projection_{1}_{2}_Density.png".format(pid, plot_axis, plot_field),fileout)


def slice(filein, xc=0.0, yc=0.0, zc=0.0, zoom_width = 2.0, fileout='junk.png', withParticles=False, field="Density") :
	pf = load(file)
	print moviefile
        plot_out_prefix = 'movieframe'
	# This is not the way/place to do this.
	#if not (withParticles): 
	#	print "No Particles use max density location"
	#	sp = pf.h.sphere("max", (1e-2, "pc"))
	
	sp = pf.h.sphere([xc, yc, zc], (zoom_width + 0.5, "pc"))

	# Get the angular momentum vector for the sphere.
	L = sp.quantities["AngularMomentumVector"]()

	# Create an OffAxisSlicePlot on the object with the L vector as its normal
	plot_field = 'Density'
	p = OffAxisSlicePlot(pf, L, plot_field, sp.center, (zoom_width, "pc"))
        #p = OffAxisSlicePlot(pf, L, "Density", sp.center, (1e-1, "pc"))
	p.set_zlim("Density", 1e-23,1e-14)
	#p.annotate_velocity(factor=16)
	pid = os.getpid()
        fileout="{0}_{1}_{2:04d}_{3}_{4:03d}_{5}pc.png".format(plot_out_prefix, quad, i, plot_field, j, zoom_width)
        p.save("junk{0:06d}".format(pid))
	shutil.move("junk{0:06d}_OffAxisSlice_{1}.png".format(pid,field),fileout)

def mslice( filein, fileout = 'junk.png',field="Density") :
	newImage = Image.new("RGB", (1070,2000))
	pf = load(file)
	print moviefile
	center = na.array([0.0, 0.0, -5.795e16])
#	center = na.array([0.0, 0.0, 0e12])
	#p = SlicePlot(pf, "z", "dens", center=center )
	p = SlicePlot(pf, "x", field, center="max", width=(4,"pc") )
	p.set_zlim("Density", 1e-23,1e-19)
	p.annotate_velocity(factor=16)
	
	pid = os.getpid()
        p.save("junk{0:06d}".format(pid))

	p = SlicePlot(pf, "y", field, center="max", width=(4,"pc") )
	p.set_zlim("Density", 1e-23,1e-19)
	p.annotate_velocity(factor=16)
        p.save("junk{0:06d}".format(pid))

	file1 = "junk{0:06d}_Slice_x_{1}.png".format(pid,field)
	file2 = "junk{0:06d}_Slice_y_{1}.png".format(pid,field)
	image1 = Image.open(file1)
	image2 = Image.open(file2)

	newImage.paste( image1, (0,0))
	newImage.paste( image2, (0,1000))
	newImage.save(fileout)

	os.remove(file1)
	os.remove(file2)

def multislice( filein, field="Density", fileout='junk.png') :
	pf = load(filein)
	orient = "horizontal"
	
	fig, axes, colorbars = get_multi_plot(2, 1, colorbar=orient, bw=4)
	#pc = PlotCollection(pf, [0., 0., 0.])
	pc = PlotCollection(pf, [2.4e19,2.4e19,2.4e19])

	p = pc.add_slice( field, "x", axes= axes[0][0], use_colorbar=False)
	p.set_cmap("bds_highcontrast")
	#p.modify["velocity"](factor=16)
	p.annotate_velocity(factor=16)

	p = pc.add_slice( field, "y", axes= axes[0][1], use_colorbar=False)
	p.set_cmap("bds_highcontrast")
	#p.show_velocity(factor=16)
	#p = pc.add_slice( field, "z", axes= axes[0][2], use_colorbar=False)
	#p.set_cmap("bds_highcontrast")
	#p.show_velocity(factor=16)

	for p, cax in zip(pc.plots, colorbars) : 
		cbar = cb.Colorbar(cax, p.image, orientation=orient)
		p.colorbar = cbar
		p._autoset_label()

	fig.savefig(fileout)
	

def ray_trace( filein, fileout='junk.png') : 
	pf = load( filein)
	pf.h
	field = "Density"
	pf.field_info[field].take_log = True
	
	dd = pf.h.all_data()
	mi, ma = dd.quantities["Extrema"](field)[0]

	mi, ma = na.log10(mi), na.log10(ma)

	tf = ColorTransferFunction((mi, ma))

	c = pf.domain_center
	L = na.array([0.0, 1.0, 0.0])
	W = 1.0/pf["unitary"]
	N = 512

	cam = pf.h.camera(c, L, W, N, tf)
	tf.add_layers(10, 0.01, colormap = "RdBu_r")
	cam.snapshot(fileout)

def histogram( filein, field="Density", fileout='junk.png', xbounds=[1e-24, 1e-18]) :
        pf = load(filein)
	pc = PlotCollection(pf)
	data = pf.h.all_data()
	if( xbounds == None) : 
		pc.add_profile_object(data, [field, "CellMassMsun"], weight=None)
	else :
		pc.add_profile_object(data, [field, "CellMassMsun"], weight=None, x_bounds=xbounds)
	
	pid = os.getpid()
        pc.save("junk{0:06d}".format(pid))
	shutil.move("junk{0:06d}_Profile1D_0_{1}_CellMassMsun.png".format(pid,field),fileout)



def find_particle_position( filein):
	pf = load("{0}{1:04d}".format( plt_prefix, i))

	dd = pf.h.all_data()
	xp = dd["particle_posx"]
	yp = dd["particle_posy"]
	zp = dd["particle_posz"]
	partMass = dd["ParticleMassMsun"]


parser = argparse.ArgumentParser(description = "start number to end number")

parser.add_argument('start', metavar='N1', type=int)
parser.add_argument('end', metavar='N2', type=int)
parser.add_argument('step', metavar='N3', type=int)
parser.add_argument('zoom_width', metavar='N4', type=float)
parser.add_argument('plot_axis', metavar='axis_', type=str)
parser.add_argument('--particle', action='store_true')
parser.add_argument('--chk', action='store_true')
parser.add_argument('--noparticle', action='store_true')
# What kind of visualization do we want?
parser.add_argument('--project', action='store_true')
parser.add_argument('--projectfull', action='store_true')
parser.add_argument('--slice', action='store_true')


#parser.add_argument('--parsec', metavar='N4', action='store_true')
args = parser.parse_args()
zoom_width = args.zoom_width
withParticles=args.particle
withNoParticles = args.noparticle
withProjection=args.project
withProjectionFull=args.projectfull
withSlice=args.slice
plot_axis = args.plot_axis
withCheckpoint = args.chk
if (withCheckpoint) :
	prefix = "BB_hdf5_chk_"	
else:
	prefix = "BB_hdf5_plt_cnt_"
quad = os.getcwd()[-5:]

if zoom_width > 5.0:
	print "The width of the box extends beyond the data values"
	print "that this quadrant calculated with AMR."
	print "System now exiting."
	sys.exit()

for i in range(args.start,args.end,args.step) :
	pf = load("{0}{1:04d}".format(prefix, i))
	dd = pf.h.all_data()
	if (withNoParticles):
		max= dd.quantities["MaxLocation"]("Density")
		maxDens = max[0]/3e-22
		maxLoc = numpy.array(max[2:5])/3e18
		xc = max[2]
		yc = max[3]
		zc = max[4]
		file = prefix+"{0:04d}".format(i)
		print file
		j = 0
		if (withProjection) or (withProjectionFull):
			moviefile = "movie_{0}_frame_{1:04d}_prjct_{2}pc_{3:04d}.png".format(quad, i, zoom_width, j)
			#project( file, xc, yc, zc, zoom_width, moviefile, withParticles)
			Proj_plotting(pf, xc, yc, zc, zoom_width, quad, i, j)
		if (withSlice):
			moviefile = "movie_{0}_frame_{1:04d}_slice_{2}pc_{3:04d}.png".format(quad, i, zoom_width, j)
			slice( file, xc, yc, zc, zoom_width, moviefile, withParticles)
			
	else:
		xp = dd["particle_posx"]
		yp = dd["particle_posy"]
		zp = dd["particle_posz"]
		#partMass = dd["ParticleMassMsun"]
		for j in range(xp.size) :
			xc = xp[j]
			yc = yp[j]
			zc = zp[j]
			file = prefix+"{0:04d}".format(i)
			print file
			print "Currently on particle ", j + 1, "of ", xp.size
			if (withProjection) or (withProjectionFull):
				moviefile = "movie_{0}_frame_{1:04d}_prjct_{2}pc_{3:04d}.png".format(quad, i, zoom_width, j)
			#project( file, xc, yc, zc, zoom_width, moviefile, withParticles)
				Proj_plotting(pf, xc, yc, zc, zoom_width, quad, i, j)
			if (withSlice):
				moviefile = "movie_{0}_frame_{1:04d}_slice_{2}pc_{3:04d}.png".format(quad, i, zoom_width, j)
				slice( file, xc, yc, zc, zoom_width, moviefile, withParticles)

print "Finished, closing up shop"
	#project( file, moviefile, withParticles=args.particle)
	#mslice( file, moviefile)
        #slice( file, moviefile)
	#slice( file, moviefile,field="accx")
	#multislice( file, field="Density", fileout=moviefile)
	#multislice( "one_proc.hdf5", field="accx", fileout=moviefile)
	#ray_trace(file, moviefile)
	#histogram( file, field="Density", fileout=moviefile)
