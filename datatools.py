# -*- coding: utf-8 -*-
"""
Front Matter
=============

Created on Fri Jun 21 12:42:46 2013

Author: Robie Hennigar

A compilation of functions that are used for analysing the data 
contained in .nc files.

Requirements
===================================
Absolutely Necessary:

* Numpy
* SciPy
* Matplotlib version 1.3.0
* Numexpr

Optional, but recommended:

* Numba

.. Note:: A convenient package that is both easy to install and contains nearly all the required packages for the python code herein is Anaconda (available here http://www.continuum.io/downloads).  Note that if you choose to use  Anaconda, it is necessary to install matplotlib v1.3.0  manually, as Anaconda comes with an older version.


Functions
=========
"""
#load modules
#Numerical modules
import numpy as np
import matplotlib.tri as mplt
import bisect
import numexpr as ne


#I/O modules
import glob
from scipy.io import netcdf
import scipy.io as sio


def loadnc(datadir, singlename=None, dim='2D'):
    """Loads a .nc  data file
    
    :Parameters: 
    	**datadir** -- The path to the directory where the data is stored. 
    	
    	**singlename (optional)** -- the name of the .nc file of interest,
            not needed if there is only one .nc file in datadir
        
        **dim = {'2D', '3D'} (optional)** -- The dimension of the data of
            interest.  Default is 2D for a 2D data file.  
    """
    #identify the file to load
    if singlename == None:
        files = glob.glob(datadir + "*.nc")
        filepath = files[0]
    else: 
        filepath = datadir + singlename
   
    #initialize a dictionary for the data.
    data = {}
    #load data
    ncid = netcdf.netcdf_file(filepath, 'r')
    data['x'] = ncid.variables['x'].data
    data['y'] = ncid.variables['y'].data
    data['ua'] = ncid.variables['ua'].data
    data['va'] = ncid.variables['va'].data
    data['nv'] = np.transpose(ncid.variables['nv'].data.astype(int))-[1] #python index
    data['zeta'] = ncid.variables['zeta'].data
    data['h'] = ncid.variables['h'].data
    data['time'] = ncid.variables['time'].data
    data['siglev'] = ncid.variables['siglev'].data
    data['siglay'] = ncid.variables['siglay'].data
    data['nbe'] = np.transpose(ncid.variables['nbe'].data.astype(int))-[1] #python index
    data['a1u'] = np.transpose(ncid.variables['a1u'].data)
    data['a2u'] = np.transpose(ncid.variables['a2u'].data)
    data['nele'] = ncid.dimensions['nele']
    data['node'] = ncid.dimensions['node']

    if dim == '3D':
        data['u'] = ncid.variables['u'].data
        data['v'] = ncid.variables['v'].data
        data['ww'] = ncid.variables['ww'].data
    elif (dim != '2D' and dim != '3D'):
        raise Exception("Dim must be '2D', '3D', or absent.")
    ncid.close()
    
    #Now we get the long/lat data.  Note that this method will search 
    #for long/lat files in the datadir and up to two levels above 
    #the datadir.  
    
    long_matches = []
    lat_matches = []
    
    if glob.glob(datadir + "*_long.dat"):
    	long_matches.append(glob.glob(datadir + "*_long.dat")[0])
    if glob.glob(datadir + "*_lat.dat"):
    	lat_matches.append(glob.glob(datadir + "*_lat.dat")[0])
    if glob.glob(datadir + "../*_long.dat"):
    	long_matches.append(glob.glob(datadir + "../*_long.dat")[0])
    if glob.glob(datadir + "../*_lat.dat"):
    	lat_matches.append(glob.glob(datadir + "*../_lat.dat")[0])
    if glob.glob(datadir + "../../*_long.dat"):
    	long_matches.append(glob.glob(datadir + "../../*_long.dat")[0])
    if glob.glob(datadir + "../../*_lat.dat"):
    	lat_matches.append(glob.glob(datadir + "../../*_lat.dat")[0])
          
    #let's make sure that long/lat files were found.
    if (len(long_matches) > 0 and len(lat_matches) > 0):
        data['lon'] = np.loadtxt(long_matches[0])
        data['lat'] = np.loadtxt(lat_matches[0])
    else:
        print "No long/lat files found. Long/lat set to x/y"
        data['lon'] = data['x']
        data['lat'] = data['y']

    data['trigrid'] = mplt.Triangulation(data['lon'], data['lat'], \
        data['nv'])
    return data
        
def ncdatasort(data):
    """From the nc data provided, common variables are produced.
    
    :Parameters: **data** -- a data dictionary of data from a .nc file
         
	| 
	       
    :Returns: **data** -- Python data dictionary updated to include uvnode and uvnodell
    """

    nodexy = np.zeros((len(data['x']),2))
    nodexy[:,0]
    
    x = data['x']
    y = data['y']
    nv = data['nv']
    lon = data['lon']
    lat = data['lat']
    
    #make uvnodes by averaging the values of ua/va over the nodes of
    #an element
    uvnode = np.empty((len(nv[:,0]),2))
    uvnode[:,0] = (x[nv[:,0]] + x[nv[:,1]] + x[nv[:,2]]) / 3.0
    uvnode[:,1] = (y[nv[:,0]] + y[nv[:,1]] + y[nv[:,2]]) / 3.0

    nodell = np.empty((len(lon),2))
    nodell[:,0] = lon
    nodell[:,1] = lat
    uvnodell = np.empty((len(nv[:,0]),2))
    uvnodell[:,0] = (lon[nv[:,0]] + lon[nv[:,1]] + lon[nv[:,2]]) / 3.0
    uvnodell[:,1] = (lat[nv[:,0]] + lat[nv[:,1]] + lat[nv[:,2]]) / 3.0
    
    data['uvnode'] = uvnode
    data['uvnodell'] = uvnodell
    data['nodell'] = nodell
    
    return data
def tri_finder(XY, data):
    """Determines which element a given X, Y value is in.
    
    :Parameters: **XY** -- an N x 2 dimensional array of points to locate in triangles. The first column consists of the x coordinate (long), while the second column consists of the y coordinate (lat).
    
	**data** -- dictionary containing all the data from the .nc file.
	
    :Returns: **indicies** -- a list of indices, in order, giving the 			triangle that each x,y pair of XY belongs to.
        
    .. Note:: if indices contains -1, this means that the cooresponding XY 
        value did not belong to any of the elements.
    """
    
    #get the relevant variables
    trigrid = data['trigrid']
    trifinder = mplt.TrapezoidMapTriFinder(trigrid)
    indices = [] #indices of the triangle where the points are found.
    #find the polygon that X, Y is in.
        
    l = len(XY[:,0])
    
    for i in xrange(l): 
        element = trifinder.__call__(XY[i,0],XY[i,1])
        indices.append(element.item())
        
    #test for correctness
    badInds = np.where(indices == -1)[0]
    if len(badInds) > 0:
        print "No triangle found for some of the points given.  Please \
            ensure that the following entries in XY are correct: "
        print badInds
    return indices

                     
def interp_vel(XY, triInds, data, numjit=False):
    """Interpolates velocity data from the known locations to the x, y
    coordinates given in the array XY.
    
    :Parameters:
        **XY** -- Nx2 array of x, y values (in long/lat) where the velocity should be interpolated to.
        
        **TriInds** -- indices of the triangles that each XY point belongs to.  This is returned by tri_finder.
        
        **data** -- Typical data dictionary returned by loadnc or ncdatasort.
        
        **jit = {True, False} (optional)** -- Optional just in time compilation of code.  This is very worthwhile for larger datasets as the compiled code runs much faster.  Requires Numba.
    
    .. Note:: Since this contains a (potentially large) for loop, the first 
        call of the code will be slow.  For faster interpolation 
        (10-100 times faster) install Numba and make use of the jit (just in
        time) compilation feature.
    """
    if numjit:
        #this assumes numba is installed on the users computer.
        #we must first cast the data in float form to ensure it will work.
        #This method is incredibly fast compared to the numba-free method,
        #so it is highly recommended for interpolating large numbers of points.
        from numba import autojit
        u = np.transpose(data['ua']).astype(float)
        v = np.transpose(data['va']).astype(float)  
        uvnodell = data['uvnodell']
        a1u = data['a1u'].astype(float)
        a2u = data['a2u'].astype(float)    
        nbe = data['nbe']
        nv = data['nv']
        
        inds = np.array(triInds)
        numba_interp = autojit()(jit_interp_vel)
        UI, VI = numba_interp(XY, inds, u, v, uvnodell, a1u, a2u, nbe, nv)        
    else:  
        #name the required data pieces
        u = np.transpose(data['ua'])
        v = np.transpose(data['va'])  
        uvnodell = data['uvnodell']
        a1u = data['a1u']
        a2u = data['a2u']    
        nbe = data['nbe']
        nv = data['nv']
         
        i = 0
        lnv = len(nv[:,0])
        UI = np.empty((len(XY[:,0]), len(u[0,:])))
        VI = np.empty(UI.shape)
        for j in triInds:
            if j != 0 and j != -1:
                e = nbe[j,:]
                e[e == 0] = lnv - 1
                xy0c = XY[i,:] - uvnodell[j,:]
                
                #Now we calculate derivatives.  Note that we attempt to do this
                #in a vectorized fashion, as much as possible, while 
                #avoiding the large matrix products that occur in 
                #totally vectorized code
                u_tmp = np.vstack((u[j,:], u[e,:]))
                v_tmp = np.vstack((v[j,:], v[e,:]))
                
                a1 = a1u[j,:]
                a2 = a2u[j,:]
                
                dudx = np.dot(a1, u_tmp)
                dudy = np.dot(a2, u_tmp)
                
                dvdx = np.dot(a1, v_tmp)
                dvdy = np.dot(a1, v_tmp)
                
                UI[i,:] = u[j,:] + xy0c[0] * dudx + xy0c[1] * dudy
                VI[i,:] = v[j,:] + xy0c[1] * dvdx + xy0c[1] * dvdy
                
            
            else:
                UI[i,:] = np.nan 
                VI[i,:] = np.nan
            i += 1
    return UI, VI

def jit_interp_vel(XY, triInds, u, v, uvnodell, a1u, a2u, nbe, nv):
    """
    .. warning:: This function is not meant to be called by the user. Use "interp_vel".
    """
    i = 0
    lnv = len(nv[:,0])
    UI = np.empty((len(XY[:,0]), len(u[0,:])))
    VI = np.empty(UI.shape)
    for j in triInds:
        if j != 0 and j != -1:
            e = nbe[j,:]
            e[e == 0] = lnv - 1
            xy0c = XY[i,:] - uvnodell[j,:]
            
            #Now we calculate derivatives.  Note that we attempt to do this
            #in a vectorized fashion, as much as possible, while 
            #avoiding the large matrix products that occur in 
            #totally vectorized code
            u_tmp = np.vstack((u[j,:], u[e,:]))
            v_tmp = np.vstack((v[j,:], v[e,:]))
            
            a1 = a1u[j,:]
            a2 = a2u[j,:]
            
            dudx = np.array([np.dot(a1, u_tmp)])
            dudy = np.array([np.dot(a2, u_tmp)])
            
            dvdx = np.array([np.dot(a1, v_tmp)])
            dvdy = np.array([np.dot(a1, v_tmp)])
            
            UI[i,:] = u[j,:] + xy0c[0] * dudx + xy0c[1] * dudy
            VI[i,:] = v[j,:] + xy0c[1] * dvdx + xy0c[1] * dvdy
            
        
        else:
            UI[i,:] = np.nan 
            VI[i,:] = np.nan
            i += 1
    return UI, VI
    
def get_elements(data, region):
    """Takes uvnodes and a  region (specified by the corners of a
    rectangle) and determines the elements of uvnode that lie within the
    region
    """
    uvnode = data['uvnodell']
    elements = np.where((uvnode[:,0] >= region[0]) & \
        (uvnode[:,0] <= region[1]) & \
        (uvnode[:,1] >= region[2]) & \
        (uvnode[:,1] <= region[3]))[0]
        
    return elements
    
def get_nodes(data, region):
    """Takes nodexy and a region (specified by the corners of a rectangle)
    and determines the nodes that lie in the region
    """
    nodexy = data['nodell']
    region = np.abs(region)
    nodexy = np.abs(nodexy)
    nodes = np.where((nodexy[:,0] >= region[0]) & \
        (nodexy[:,0] <= region[1]) & \
        (nodexy[:,1] >= region[2]) & \
        (nodexy[:,1] <= region[3]))[0]
    return nodes 
    
def regioner(region, data, name=None, savedir=None, dim='2D'):
    """
    Takes as input a region (given by a four elemenTakes as input a region 
    (given by a four element NumPy array),
    and some standard data output by ncdatasort and loadnc2d_python
    and returns only the data that lies within the region specified
    in the region arrayt NumPy array),
    and some standard data output by ncdatasort and loadnc2d_python
    and returns only the data that lies within the region specified
    in the region array
    
    :Parameters:
        **region** -- four element array containing the four corners of the 
            region box.  Entires should be in the following form:
            [long1, long2, lat1, lat2] with the following property:
            abs(long1) < abs(long2), etc.
        
        **data** -- standard python data dictionary for these files
        
        **name** -- what should  the region be called in the output file?
        
        **savedir** -- where should the resultant data be saved? Default is 			none, i.e. the data will not be saved, only returned.
        
        **dim = {'2D', '3D'}** the dimension of the data to use regioner 				on.  Default is 2D.
    """
    #short name for relevant variables
    
    nv = data['nv']
    nbe = data['nbe']
    a1u = data['a1u']
    a2u = data['a2u']
    
    
    l = nv.shape[0]
    if savedir == None:
        regionData = "/not/a/real/path"
    else:
        regionData = savedir + name + "_region_" + str(region[0]) \
    		+ "_" + str(region[1]) + "_" + str(region[2]) + "_" + str(region[3]) + ".mat"
    files = glob.glob(regionData)
    
    #find the nodes that lie in the region
    idx = get_nodes(data, region)
    if len(files)==0:
        #There is currently no file with this particular region data
        #build a new data set for this region
        
        #first, reindex elements in the region
        element_index_tmp = np.zeros(l, int)
        nv_rs = nv.reshape(l*3, order='F')
        #find indices that sort nv_rs
        nv_sortedind = nv_rs.argsort()
        #sort the array
        nv_sortd = nv_rs[nv_sortedind]                
        #pick out the values in the region
        for i in xrange(len(idx)):
            i1 = bisect.bisect_left(nv_sortd, idx[i])
            i2 = bisect.bisect_right(nv_sortd, idx[i])
            inds = nv_sortedind[i1:i2]
            element_index_tmp[inds % l] = 1
        element_index = np.where(element_index_tmp == 1)[0]
        node_index = np.unique(nv[element_index,:])
        #create new linkage arrays
        nv_tmp = nv[element_index,:]
        L = len(nv_tmp[:,0])
        nv_tmp2 = np.empty((1, L*3.0))
        
        #make a new array of the node labellings for the tri's in the region
        
        nv2 = nv_tmp.reshape(L * 3, order='F')
        nv2_sortedind = nv2.argsort()
        nv2_sortd = nv2[nv2_sortedind]
        
        for i in xrange(len(node_index)):
            i1 = bisect.bisect_left(nv2_sortd, node_index[i])
            i2 = bisect.bisect_right(nv2_sortd, node_index[i])
            inds = nv2_sortedind[i1:i2]
            nv_tmp2[0, inds] = i 
        
        nv_new = np.reshape(nv_tmp2, (L, 3), 'F')
        #now do the same for nbe
        nbe_index = np.unique(nbe[element_index, :])
        nbe_tmp = nbe[element_index,:]
        lnbe = len(nbe_tmp[:,0])
        nbe_tmp2 = np.empty((1, lnbe*3))

        nbe2 = nbe_tmp.reshape(lnbe*3, order='F')
        nbe_sortedind = nbe2.argsort()
        nbe_sortd = nbe2[nbe_sortedind]
        
        for i in xrange(len(nbe_index)):
            i1 = bisect.bisect_left(nbe_sortd, nbe_index[i])
            i2 = bisect.bisect_right(nbe_sortd, nbe_index[i])
            inds = nbe_sortedind[i1:i2]
            nbe_tmp2[0, inds] = i
        
        nbe_new = np.reshape(nbe_tmp2, (lnbe,3), 'F')
        nbe_new[nbe_new > len(nv_new[:,0]), :] = 0
        
        #create new variables for the region
        data['node_index'] = node_index
        data['element_index'] = element_index
        data['nbe'] = nbe_new
        data['nv'] = nv_new
        data['a1u'] = a1u[element_index, :]
        data['a2u'] = a2u[element_index, :]
        data['h'] = data['h'][node_index]
        data['uvnodell'] = data['uvnodell'][element_index,:]
        data['x'] = data['x'][node_index]
        data['y'] = data['y'][node_index]
        data['zeta'] = data['zeta'][:,node_index]
        data['ua'] = data['ua'][:,element_index]
        data['va'] = data['va'][:,element_index]
        data['lon'] = data['lon'][node_index]
        data['lat'] = data['lat'][node_index]
        data['trigrid'] = mplt.Triangulation(data['lon'], data['lat'], \
            data['nv'])
        #take care of extra variables if data file is 3D    
        if dim=='3D':
            
            data['u'] = data['u'][:,:,element_index]
            data['v'] = data['v'][:,:,element_index]
            data['ww'] = data['ww'][:,:,element_index]
        #save the data if that was requested.
        if savedir != None and name != None:
            mat_save(data, regionData)
    return data

def mat_save(data, saveDirName, dim='2D'):
    """
    Save .nc data to a mat file.
    
    :Parameters:
        **data** -- the standard data dictionary
        
        **saveDirName** -- the path to where the data should be saved,
        along with the name. Ex: "/home/user/Desktop/data.mat"
        
        **dim ={'2D', '3D'} (optional)** -- the dimension of the data file.
        
    """
    dtype = float  
    rdata={}
    if dim=='3D': 
        rdata['a1u'] = data['a1u'].astype(dtype)
        rdata['a2u'] = data['a2u'].astype(dtype)
        rdata['h'] = data['h'].astype(dtype)
        rdata['uvnodell'] = data['uvnodell']
        rdata['x'] = data['x'].astype(dtype)
        rdata['y'] = data['y'].astype(dtype)
        rdata['u'] = data['u'].astype(dtype)
        rdata['v'] = data['v'].astype(dtype)
        rdata['ww'] = data['ww'].astype(dtype)
        rdata['zeta'] = data['zeta'].astype(dtype)
        rdata['h'] = data['h'].astype(dtype)
        rdata['ua'] = data['ua'].astype(dtype)
        rdata['va'] = data['va'].astype(dtype)

    elif dim=='2D':
        rdata['a1u'] = data['a1u'].astype(dtype)
        rdata['a2u'] = data['a2u'].astype(dtype)
        rdata['h'] = data['h'].astype(dtype)
        rdata['uvnodell'] = data['uvnodell']
        rdata['x'] = data['x'].astype(dtype)
        rdata['y'] = data['y'].astype(dtype)
        rdata['zeta'] = data['zeta'].astype(dtype)
        rdata['ua'] = data['ua'].astype(dtype)
        rdata['va'] = data['va'].astype(dtype)
    else:
        raise Exception("Dim must be '2D', '3D', or absent.")
    sio.savemat(saveDirName, rdata, oned_as='column')
    
def h5_save(data, savedir, filename, cast=False):
	"""Save data into an htf5 format.  This is faster than saving to 
	.mat files.  Note this this code assumes that the data is already cast,
	unless it explicitly told to cast  the data by setting cast=True.
	
	:Parameters:
		**data** -- Standard python data dictionary with data to be saved.
		
		**savedir** -- The directory where the resultant hft5 file should be 				saved.  Include a '/' at the end, like this: /path/to/save/
		
		**filename** -- The name for the file.
		
		**cast = {True, False}** -- Optional. If True the data will be cast 		to float form before saving.  This is necessary if the data has not 		yet been cast.
	"""
	#a list of data entries that will need to be casted.
	toCast = ['x', 'y', 'lon', 'lat', 'h', 'a1u', 'a2u']
	#a list of data entries that would not need to be cast to floats
	noCast = ['nv', 'nbe']
	#make sure there are not any 'unsavable' structures in the dictionary
	key = data.keys()
	#initialize list of items not to be saved.
	noWrite = []
	for i in key:
		try:
			data[i].shape
		except:
			noWrite.append(i)
	#cast the data to float form.
	if cast:
		for i in toCast:
			data[i] = data[i].astype('float32')
	#get the dictionary that will be saved ready.
	sdict = {i:data[i] for i in key if i not in noWrite}
	#initalize a dictionary of random variables.  This serves 
	#only to make it possible to do the following in a condensed form
	rdata = {}
	for i in sdict.keys():
		rdata[i] = 4
	
	#save to h5 file.
	f = h5py.File(savedir + filename + '.h5py', 'w')
	for i in sdict.keys():
		if i not in noCast:
			rdata[i] = f.create_dataset(i, sdict[i].shape, 'f')
			rdata[i][...] = sdict[i]
		else:
			rdata[i] = f.create_dataset(i, sdict[i].shape, 'i')
			rdata[i][...] = sdict[i]
	f.close()
def calc_speed(data):
    """
    Calculates the speed from ua and va
    
    :Parameters:
        **data** -- the standard python data dictionary

    .. note:: We use numexpr here because, with multiple threads, it is
    about 30 times faster than direct calculation.        
    """
    #name required variables
    ua = data['ua']
    va = data['va']

    #we can take advantage of multiple cores to do this calculation
    ne.set_num_threads(ne.detect_number_of_cores())
    #calculate the speed at each point. 
    data['speed'] = ne.evaluate("sqrt(ua*ua + va*va)")
    return data
    
def calc_energy(data):
    """Calculate the energy of the entire system.

    :Parameters: **data** -- the standard python data dictionary
    """
    #name relevant variables
    x = data['x']
    y = data['y']
    nv = data['nv']
    rho = 1000 #density of water
    #initialize EK to zero
    l = nv.shape[0]
    area = np.zeros(l)
    for i in xrange(l):
        #first, calculate the area of the triangle
        xCoords = x[nv[i,:]]
        yCoords = y[nv[i,:]]
        #Compute two vectors for the area calculation.  
        v1x = xCoords[1] - xCoords[0]
        v2x = xCoords[2] - xCoords[0]
        v1y = yCoords[1] - yCoords[0]
        v2y = yCoords[2] - yCoords[0]
        #calculate the area as the determinant
        area[i] = abs(v1x*v2y - v2x*v1y)
    #get a vector of speeds.
    sdata = calc_speed(data)
    speed = sdata['speed']   
        
    #calculate EK, use numexpr for speed (saves ~15 seconds)
    ne.set_num_threads(ne.detect_number_of_cores())
    ek = ne.evaluate("sum(rho * area * speed * speed, 1)")
    ek = ek / 4     
    data['ek'] = ek
    print ek.shape
    return data
        
  
    
