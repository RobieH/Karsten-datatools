from mpi4py import MPI
import scipy.io as sio
from datatools import *
from scipy.io import netcdf as nc

datadir = '/home/robie/Desktop/Practice/'
saveName = '/home/robie/Desktop/parSave.nc'

if __name__ == '__main__':
	#open up communication
	comm = MPI.COMM_WORLD
	size = comm.size
	rank = comm.rank
	
	#load data to all processors
	data = loadnc(datadir)
	
	#the goal here is to perform the casting of data on multiple cores,
	#as this is the part of the saving that has previously taken 
	#the most time.
	
	#initialize the dictionary on all processors
	rdata = {}
	#initialize the dictionary key as a list on all processors
	key = data.keys()
	
	#we want to have some processors cast some data, and some processors 
	#cast others.  To do this, we will make a list on each processor
	#containing the indices of certain dictionary elements. 
	
	#figure out how many elements are in the dictionary
	ldict = len(key)
	comm.Barrier()
		
	#make separate lists on each processor for the elements of the 
	#dictionary each processor is responsible for
	entries = [rank + size * i for i in range(int(ldict/size) + 1) if rank + size * i < ldict]
	
	#now we start the casting procedure.  Note that certain matrices do
	#not need to be casted, so we insert exceptions.  Furthermore, 
	#certain entries cannot be saved, so we check for these as well. 
	#To make this process a bit easier, we make a list with the invalid 
	#entries.
	
	if rank == 0:
		#create a list of invalid keys
		invalid = []
		for i in key:
			try:
				data[i].dtype
			except:
				invalid.append(i)
	else:
		invalid = None
	#make a list of the entries that will not require casting.
	noCast = ['nv', 'nbe', 'lon', 'lat', 'uvnodell', 'nodell']
	
	#broadcast the result to all processors.
	invalid = comm.bcast(invalid, root=0)
	comm.Barrier()
	
	#all processors now have a copy of the dictionary indices, and
	#those  which are invalid, so let us start casting the data.
	
	#first have each processor generate its key
	toCast = [key[i] for i in entries if (key[i] not in invalid and key[i] not in noCast)]
	comm.Barrier()
	#have each processor cast its data
	dtype = float
	for i in toCast:
		rdata[i] = data[i].astype(dtype)
	#the data that needed to be casted has been casted.  Now we begin
	#saving the data.  We let the head node create the file, save its
	#structures, along with the no cast structures.
	if rank == 0:
		for i in noCast:
			try:
				rdata[i] = data[i]
			except:
				pass
"""THE CODE WORKS UP TO THIS POINT.  NEED TO FIGURE OUT A WAY TO APPEND DATA TO AN ALREADY CREATED FILE... PERHAPS USING .NC FILES AS AN INTERMEDIATE WOULD BE A GOOD IDEA"""

		#initialize netcdf file
		ncid = nc.netcdf_file(saveName, 'w')
		for j in rdata.keys():
			ncid.createVariable(j, 
			
	comm.Barrier()	
	#gather the data
#	if rank != 0:
#		sio.savemat(saveName, rdata)
#	toSave = [key[i] for i in entries if key[i] not in invalid]
#	comm.Barrier()
	
