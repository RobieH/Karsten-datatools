from plottools import *
from mpi4py import MPI


def save_plot(i):
	plt.tripcolor(trigrid, speed[i,:])
	plt.colorbar()
	plt.title('speed' + ' Time Series' + str(i))
	plt.savefig(savedir + 'plot' + '_' + str(i) + '.png')
	plt.clf()

if __name__ == '__main__':
	comm = MPI.COMM_WORLD
	root = 0
	if comm.rank == 0:
		print "loading data"
		data = loadnc('/home/robie/Desktop/Practice/')
		print "data loaded.  Calculating speed."
		sdata = calc_speed(data)
		speed = sdata['speed']
		print "about to start main process."
		trigrid = data['trigrid']
		savedir = '/home/robie/Desktop/plots/'
	comm.Barrier()
	rows = [comm.rank + comm.size * i for i in range(int(20/comm.size)+1) if comm.rank + comm.size*i < 20]
	if comm.rank == root:
		if not os.path.exists(savedir):
			os.makedirs(savedir)
	comm.Barrier()
	for i in rows:
		save_plot(i)
