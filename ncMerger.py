"""
@author: Robie Hennigar

This is optimized python code to merge and save multiple .nc files.
"""

from datatools import *
import scipy.io as sio
from scipy.io import netcdf
import numpy as np
import glob
import bisect
import os
from datetime import datetime, timedelta
import h5py

#enter the data directory
datadir = '/home/robie/Documents/python/'
#enter the number of days of data per .mat file
numDays = 0.5
#enter the directory where the .mat files should be saved
savedir = '/home/robie/Desktop/Python/'
#enter name for .mat files
filename = "merged"

def size_check(datadir):
	"""Determine how large the merged arrays will be, for preallocation.
	"""
	files = glob.glob(datadir + '*.nc')
	numFiles = len(files)
	ncid = netcdf.netcdf_file(files[0])
	nele = ncid.dimensions['nele']
	node = ncid.dimensions['node']
	ncid.close()
	timeDim = 0
	for i in xrange(numFiles):
		ncid = netcdf.netcdf_file(files[i])
		timeDim += len(ncid.variables['time'].data)
		ncid.close()
	return nele, node, timeDim
	
def nan_index(data):
	"""Determines, for a given data set, what time series
	contain nans.  Also calculates the fraction of
	a data set that is nans (a measure of the 'goodness' of
	the data set).
	"""
	#name necessary variables
	time = data['time']
	ua = data['ua']
	va = data['va']
	#initialize list of nan-containing time series
	nanInd = []
	for i in xrange(time.shape[0]):
		checkArray = []
		key = ['ua', 'va']
		for j in key:
			checkArray.append(np.isnan(np.sum(u[i,:,:])))
		if True in checkArray:
			nanInd.append(i)
	#calculate percentage of nans
	nanFrac = len(nanInd)/len(time)
	return nanInd, nanFrac
	
def time_sorter(datadir):
	"""Sorts the output of glob so that earlier time series are loaded
	first.
	"""
	files = glob.glob(datadir + "*.nc")
	first_time = np.zeros(len(files))
	for i in xrange(len(files)):
		ncid = netcdf.netcdf_file(files[i])
		first_time[i] = ncid.variables['time'][0]
		ncid.close()
	ordered_time = first_time.argsort()
	return ordered_time
	
def merge_nc(datadir,savedir, intelligent=False):
	"""Merge code for N .nc files.  Currently works for 2D files only.
	
	::Parameters:
		datadir: the directory containing the .nc files
		intelligent={True,False}: Default False.  This option is useful 
			if there are several overlapping data files, some of which
			are expected to contain many nans.  If this option is set to 
			True, the code will run slower, but the output should be cleaner.
	"""
	#generate a list of all files in the directory
	files = glob.glob(datadir + '*.nc')
	numFiles = len(files)
	print numFiles
	#get the proper indices to load the files.
	ordered_time = time_sorter(datadir)
	#load the first file
	singlename = os.path.basename(files[ordered_time[0]])
	data = loadnc(datadir, singlename=singlename, dim='2D')
	#name variables
	time = data['time']
	ua = data['ua']
	va=data['va']
	zeta = data['zeta']
	
	#get the dimensions we will need
	nele, node, timeDim = size_check(datadir)
	
	#create arrays that the data will be put into.  For the large arrays
	#use memory mapping for speed.  We will memory map into the datadir.
	#this will free up RAM for calculations.
	
	Time = np.zeros(timeDim)
#	UA = np.memmap(savedir + 'ua', dtype='float32', mode='w+', shape=(timeDim,nele))
#	VA = np.memmap(savedir + 'va', dtype='float32', mode='w+', shape=(timeDim,nele))
#	ZETA = np.memmap(savedir + 'zeta', dtype='float32', mode='w+', shape=(timeDim,node))
	UA = np.zeros((timeDim, nele))
	VA = np.zeros((timeDim, nele))
	ZETA = np.zeros((timeDim, node))	
	#enter the elements of the first time series
	l = len(time)
	Time[0:l] = time
	UA[0:l,:] = ua
	VA[0:l,:] = va
	ZETA[0:l,:] = zeta
	nt = len(time) - 1
	for i in xrange(1,numFiles):
		#load a file
		ncid = netcdf.netcdf_file(files[ordered_time[i]],'r')
		ua_temp = ncid.variables['ua'].data
		va_temp = ncid.variables['va'].data
		zeta_temp = ncid.variables['zeta'].data
		time_temp = ncid.variables['time'].data
		ncid.close()
		ltt = len(time_temp)
		#discover is there is any overlap in time series
		start = bisect.bisect_left(Time,time_temp[0])
		
		#determine if there was a match
		if start != len(Time):
			#there is a matching index, i.e. the data overlaps.
			start_ind = start
			NT = ltt + start_ind
			
			Time[start_ind:NT] = time_temp
			UA[start_ind:NT,:] = ua_temp
			VA[start_ind:NT,:] = va_temp
			ZETA[start_ind:NT,:] = zeta_temp
			nt = NT
			print "had some overlap"
			
		else:
			#there was no matching index, i.e. the data does not overlap
			numT = ltt + nt
			Time[nt:numT] = time_temp
			UA[nt:numT,:] = ua_temp
			VA[nt:numT,:] = va_temp
			ZETA[nt:numT,:] = zeta_temp
			nt += ltt
		print "loaded file " + str(i+1) + " of " + str(numFiles) + "."
	data['ua'] = UA
	data['va'] = VA
	data['zeta'] = ZETA
	data['time'] = Time
	return data
if __name__ == '__main__':
	a = datetime.now()
	data = merge_nc(datadir,savedir)
	#due to the way python reads .nc files, we will need to cast some
	#of the data before saving it, otherwise it will not save correctly.
	toCast = ['x', 'y', 'lon', 'lat', 'h', 'a1u', 'a2u']
	noCast = ['nv', 'nbe']

	dtype = 'float32'
	for i in toCast:
		data[i] = data[i].astype(dtype)
	noWrite = ['trigrid', 'siglay', 'siglev']
	for i in data.keys():
		try:
			data[i].shape
		except:
			noWrite.append(i)			
	mdict = {i:data[i] for i in data.keys() if i not in noWrite}
	b = datetime.now()
	print "saving"
	sio.savemat(savedir + filename, mdict, oned_as='column')
	#f = h5py.File(savedir+filename+'.h5py', 'w')
	rdata = {}
	for i in mdict.keys():
		rdata[i] = 4
#	for i in mdict.keys():
#		if i not  in noCast:
#			rdata[i] = f.create_dataset(i, mdict[i].shape, 'f') 
#			rdata[i][...] = mdict[i]
#		else:
#			rdata[i] = f.create_dataset(i, mdict[i].shape, 'i')
#			rdata[i][...] = mdict[i]
#	f.close()
	d = datetime.now()
	c = b - a
	e = d - b
	print c.seconds, c.microseconds
	print e.seconds, e.microseconds
