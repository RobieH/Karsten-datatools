from mpi4py import MPI
import scipy.io as sio
from datatools import *

datadir = '/home/robie/Desktop/Practice/'

if __name__ == '__main__':
	#open up communication
	comm = MPI.COMM_WORLD
	size = comm.size
	rank = comm.rank
	
	#load data to all processors
	data = loadnc('datadir')
	
	#the goal here 
