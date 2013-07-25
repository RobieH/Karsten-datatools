from datatools import *
from os import path
from time import sleep

datadir = '/output/'
dim = '2D'

while True:
	try:
		files = glob.glob(datadir + "*.nc")

		#only look at the actual output files
		for i in files:
			if "restart" in i:
				files.remove(i)
		
		#get the latest run
		num = []
		for i in files:
			num.append(int(i[-7:-3]))
	
		index = num.index(max(num))

		#load the relevant file
		singlename = path.basename(files[index])
		data = loadnc(datadir, singlename, dim)

		#check for nans
		nanInd, nanFrac = nan_index(data, dim)

		if nanFrac > 0:
			break
	except:
		pass
