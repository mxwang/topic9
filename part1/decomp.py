from mpi4py import MPI
import sys

comm = MPI.COMM_WORLD()
psize = comm.Get_size()    # process size
prank = comm.Get_rank()    # process id

sys.stdout.write("Netcdf! I am process %d of %d.\n" % (prank, psize)) 
