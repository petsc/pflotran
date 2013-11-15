from __future__ import division
from matplotlib.tri import Triangulation
from optparse import OptionParser
import numpy as np

def ParseArcGIS(infile,skip,verbose=False):
    """
    Parses elevation data into a dictionary.

    Parameters:
    -----------
    infile : string
        a ASCII file containing the elevation data in the following format.

    Sample format:
    --------------
    north: 182399.5
    south: 181642.5
    east: 755265.39600081
    west: 754356.38649747
    rows: 1514
    cols: 1818
    285.056305 285.054688 285.051331 ...

    OR

    ncols 361
    nrows 337
    xllcorner 585865
    yllcorner 7910248
    cellsize 0.500000
    NODATA_value -9999
    4.57731485 4.57674599 4.58373213 ...
    """
    f = open(infile,'r')
    if verbose: print "Reading data from %s..." % (infile)
    lines = f.readlines()
    data = {}

    style = 0
    if lines[0].count(":") > 0: style = 1

    if style == 1:
        row  = 0
        for line in lines:
            if line.count(":") > 0:
                line=line.split(":")
                data[line[0]] = float(line[1])
            else:
                if data.has_key("z") == False:
                    data["z"] = np.zeros((int(data["rows"]),int(data["cols"])))
                line=line.split()
                data["z"][row,:] = np.asarray(line)
                row += 1

    if style == 0:
        row = 0
        for line in lines:
            line = line.split()
            if len(line) == 2:
                data[line[0]] = float(line[1])
            else:
                if data.has_key("z") == False:
                    data["z"] = np.zeros((int(data["nrows"]),int(data["ncols"])))
                data["z"][row,:] = np.asarray(line)
                row += 1
        data["rows"]  = data["nrows"]
        data["cols"]  = data["ncols"]
        data["south"] = data["xllcorner"]
        data["north"] = data["xllcorner"] + data["cellsize"]*data["nrows"]
        data["west"]  = data["yllcorner"]
        data["east"]  = data["yllcorner"] + data["cellsize"]*data["ncols"]

    f.close()
    if verbose:
        print "\tFound %s x %s = %s datapoints" % ("{:,}".format(int(data["rows"])),
                                                   "{:,}".format(int(data["cols"])),
                                                   "{:,}".format(int(data["rows"]*data["cols"])))
        print "\tMin/Max surface elevation = (%.6f,%.6f)" % (data["z"].min(),data["z"].max())

    X,Y = np.meshgrid(np.linspace(data["south"],data["north"],data["rows"]),
                      np.linspace(data["west"], data["east"], data["cols"]))
    data["X"] = X[::skip,::skip]
    data["Y"] = Y[::skip,::skip]
    data["Z"] = (data["z"].T)[::skip,::skip]
    if verbose and skip > 1:
        print "\tThinning dataset (skip = %d)..." % skip
        print "\t\t%s surface points reduced to %s" % ("{:,}".format(int(data["rows"]*data["cols"])),
                                                       "{:,}".format(int(data["X"].shape[0]*data["X"].shape[1])))
    return data

def CheckHandedness(vertices,cells,verbose=False):
    """
    Given the vertices and cells of a triangular mesh, check that the
    cells conform to the right-hand rule.

    Parameters:
    -----------
    vertices : float array
        array of mesh vertices (nvertices,2)
    cells : int array
        array of mesh connectivity (ncells,3)
    """
    if verbose: print "Checking surface mesh handedness..."
    nv   = vertices.shape[0]
    nc   = cells.shape[0]
    flip = 0
    for c in range(nc):
        v1 = vertices[cells[c,1],:]-vertices[cells[c,0],:]
        v2 = vertices[cells[c,2],:]-vertices[cells[c,0],:]
        if np.cross(v1,v2) < 0:
            rem        = cells[c,1]
            cells[c,1] = cells[c,2]
            cells[c,2] = rem
            flip      += 1
    if verbose: print "\tFlipped %d of %d triangular cells" % (flip,cells.shape[0])

def CreatePrismaticMesh(data,depths,verbose=False):
    """
    Given a surface dataset, creates a prismatic mesh (0-based) and
    adds the data to the input dictionary.

    Parameters:
    -----------
    data : dictionary
        output from ParseArcGIS
    depths : array of float
        depth of each layer of cells (measure to bottom of cell)
    """
    if verbose: print "Creating surface mesh..."
    X = data["X"].flatten()
    Y = data["Y"].flatten()
    Z = data["Z"].flatten()
    tri = Triangulation(X,Y)
    nsimplex  = tri.triangles.shape[0]
    simplices = tri.triangles
    if verbose: print "\tContains %s triangular elements" % ("{:,}".format(nsimplex))
    if verbose: print "Creating prismatic mesh..."
    ncells = depths.shape[0]-1
    n   = X.shape[0]       # number of points in each layer
    N   = n*(ncells+1)     # number of total points
    xyz = np.zeros((N,3))
    for i in range(ncells+1):
        xyz[i*n:(i+1)*n,0] = X
        xyz[i*n:(i+1)*n,1] = Y
        xyz[i*n:(i+1)*n,2] = Z - depths[i]
    ne  = nsimplex       # number of elements in each layer
    Ne  = ne*ncells      # number of total elements
    cells = np.zeros((Ne,6),dtype='int')
    for i in range(ncells):
        cells[i*ne:(i+1)*ne,:3] = simplices + (i+1)*n
        cells[i*ne:(i+1)*ne,3:] = simplices + i*n
    if verbose:
        print "\tTotal vertices: %s" % ("{:,}".format(N))
        print "\tTotal elements: %s" % ("{:,}".format(Ne))

    # add mesh data to the dictionary
    data["n"]     = n
    data["N"]     = N
    data["xyz"]   = xyz
    data["ne"]    = ne
    data["Ne"]    = Ne
    data["cells"] = cells # 0-based indexing
    return data

def WriteVTK(outfile,data,verbose=False):
    """
    Writes a VTK file of the mesh found in the data dictionary.

    Parameters:
    -----------
    outfile : string
        name of the data into which to write the VTK data
    data : dictionary
        generated by CreatePrismaticMesh
    """
    if verbose: print "Writing vtk file %s..." % outfile
    f = open(outfile,'w')
    f.write('# vtk DataFile Version 2.0\n')
    f.write('prismatic mesh\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS %d float\n' % (data["N"]))
    data["xyz"].tofile(f,sep=' ',format='%s')
    f.write('\nCELLS %d %d\n' % (data["Ne"],data["Ne"]*7))
    cell = np.zeros((data["Ne"],7),dtype='int')
    cell[:,0] = 6
    cell[:,1:] = data["cells"]
    cell.tofile(f,sep=' ',format='%s')
    cell = np.ones(data["Ne"],dtype='int')*13
    f.write('\nCELL_TYPES %d\n' % (data["Ne"]))
    cell.tofile(f,sep=' ',format='%s')

def WriteJou(outfile,data,max_depth,solid,verbose):
    """
    Writes a CUBIT journal file which reconstructs a smooth surface
    from the input data. Optionally can include a solid extrustion to
    max_depth.

    Parameters:
    -----------
    outfile : string
        name of the output journal file
    data : dictionary
        output from ParseArcGIS
    max_depth : float
        maximum depth under the surface to model
    solid : boolean
        flag to include a solid
    """
    def _writeSpline(points):
        crv = 'create curve spline'
        for i in range(points.shape[1]):
            crv += ' location %.16e %.16e %.16e' % (points[0,i],points[1,i],points[2,i])
        crv += '\n'
        return crv
    x = data["X"]; y = data["Y"]; z = data["Z"];
    xyz = np.asarray([x,y,z])
    if verbose: print "Writing CUBIT journal file %s..." % outfile
    f = open(outfile,'w')
    for r in range(x.shape[0]):
        crv = _writeSpline(xyz[:,r,:])
        f.write('%s' % crv)
    for c in range(x.shape[1]):
        crv = _writeSpline(xyz[:,:,c])
        f.write('%s' % crv)
    f.write('create surface net u curve 1 to %d v curve %d to %d\n' % (x.shape[0],x.shape[0]+1,x.shape[0]+x.shape[1]))
    f.write('delete curve 1 to %d\n' % (x.shape[0]+x.shape[1]))
    if solid:
        dx  = (x.max()-x.min())
        dy  = (y.max()-y.min())
        dz  = (z.max()-z.min())
        pad = 0.05
        shiftx = x.min()+dx/2
        shifty = y.min()+dy/2
        shiftz = z.min()-(((1+pad)*dz+max_depth)/2-(1+pad)*dz)
        f.write('brick x %.15e y %.15e z %.15e\n' % ((1-pad)*dx,(1-pad)*dy,(1+pad)*dz+max_depth))
        f.write('move Volume 2 location %.15e %.15e %.15e include_merged\n' % (shiftx,shifty,shiftz))
        f.write('webcut volume 2 with sheet body 1\n')
        f.write('delete volume 3\n')
        f.write('delete surface 1\n')
    f.close()

def WriteHDF5(outfile,data,ncells,materials,verbose=False):
    """
    Writes a HDF5 file compatible with PFLOTRAN. Defines regions: all,
    top, Layer[1-ncells].

    Parameters:
    -----------
    outfile : string
        name of the output h5 file
    data : dictionary
        output from ParseArcGIS
    ncells : int
        number of cells used in the depth-direction
    materials : array of int
        material ids for each cell layer (size of ncells)
    """
    if verbose: print "Writing h5 file %s..." % outfile
    import h5py as h5
    nside = 3
    f = h5.File(outfile,'w')
    V         = f.create_dataset('Domain/Vertices',(data["N"],3),dtype='f8')
    V[...]    = data["xyz"]
    C         = f.create_dataset('Domain/Cells',(data["Ne"],2*nside+1),dtype='u8')
    C[...]    = np.ones((data["Ne"],2*nside+1))*2*nside
    C[...,1:] = data["cells"]+1
    A         = f.create_dataset('Regions/all',(data["Ne"],),dtype='u8')
    A[...]    = np.asarray(range(data["Ne"]),dtype='u8')+1 # 1-based
    T         = f.create_dataset('Regions/top',(data["ne"],nside+1),dtype='u8')
    T[...]    = np.ones((data["ne"],nside+1))*nside
    T[...,1:] = data["cells"][:data["ne"],3:]+1 # 1-based
    M         = f.create_dataset('Regions/Material ID',(data["Ne"],),dtype='u8')
    for i in range(ncells):
        L     = f.create_dataset('Regions/Layer%d' % (i+1),(data["ne"],),dtype='u8')
        L[...]= np.asarray(range(i*data["ne"],(i+1)*data["ne"]),dtype='u8')+1 # 1-based
        M[i*data["ne"]:(i+1)*data["ne"]] = materials[i]
    f.close()

def Preview(data):
    import pylab as plt
    cb = plt.imshow(data["Z"],cmap='terrain')
    plt.colorbar(cb)
    plt.show()

def SplitList(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

desc = """This script reads in elevation data in ArcGIS format to create various
intermediate representations of that data. It can be used to create a
prismatic mesh that can be viewed as a VTK file or output in a
PFLOTRAN-ready hdf5 format. It can also be used to generate CUBIT
journal files that can be used to create a surface/solid
representation that CUBIT can then mesh.
"""
epil = """Send questions or problems to:
Nathan Collier,
nathaniel.collier@gmail.com"""
parser = OptionParser(usage="python %prog [options] datafile.txt",
                      description=desc,
                      epilog=epil)
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="suppress verbose information")
parser.add_option("-p", "--preview",
                  action="store_true", dest="preview", default=False,
                  help="show color contour plot of thinned out dataset")
parser.add_option("--h5",
                  action="store_true", dest="h5", default=False,
                  help="output a PFLOTRAN hdf5 file")
parser.add_option("--vtk",
                  action="store_true", dest="vtk", default=False,
                  help="output a vtk file")
parser.add_option("--jou",
                  action="store_true", dest="jou", default=False,
                  help="output a CUBIT jou file which reconstructs the surface")
parser.add_option("--solid",
                  action="store_true", dest="solid", default=False,
                  help="appends a solid extrustion to the file generated with the --jou option")
parser.add_option("-s","--skip",
                  action="store", type="int", dest="skip", default=1,
                  help="increase to thin out the dataset")
parser.add_option("-n","--ncells",
                  action="store", type="int", dest="ncells", default=-1,
                  help="number of cells in the depth direction")
parser.add_option("-d","--depth",
                  action="callback", type="string", dest="depth",default=["40"],
                  callback=SplitList,
                  help="depth below surface to model")
parser.add_option("-m","--material",
                  action="callback", type="string", dest="material",default=["1"],
                  callback=SplitList,
                  help="specify material ids for each depth")

(options,args) = parser.parse_args()
if options.ncells < 0:
    options.depth.insert(0,0)
    options.depth = np.asarray(options.depth,dtype=float)
else:
    options.depth = np.linspace(0,float(options.depth[-1]),options.ncells+1)
options.ncells = options.depth.shape[0]-1
if len(options.material) == 1:
    options.material = np.ones(options.ncells,dtype=int)*int(options.material[0])
else:
    options.material = np.asarray(options.material,dtype=int)
    if options.material.shape[0] != options.depth.shape[0]-1:
        raise ValueError("You must specify a material for each layer! [len(depth)=%d] != [len(material)=%d]" % 
                         (options.material.shape[0],options.depth.shape[0]-1))
    
if len(args) == 0 or not(options.vtk or options.h5 or options.jou or options.preview):
    parser.print_help()
    from sys import exit
    exit(1)
fhead = args[0].split(".")[0]
data = ParseArcGIS(args[0],options.skip,options.verbose)
if options.preview: Preview(data)
if (options.vtk or options.h5):
    xyz = CreatePrismaticMesh(data,options.depth,options.verbose)
    if options.vtk: WriteVTK("%s.vtk" % fhead,data,options.verbose)
    if options.h5:  WriteHDF5("%s.h5" % fhead,data,options.ncells,options.material,options.verbose)
if options.jou or options.solid:
    WriteJou("%s.jou" % fhead,data,options.depth.max(),options.solid,options.verbose)
