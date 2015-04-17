# Domain Decoposito
# Muxuan Wang

#!/usr/bin/env python
from netCDF4 import Dataset
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as cols
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from struct import *
from mpi4py import MPI

class Netcdf_Reader:

    #return number of dimensions
    def ndims(self,dims):
        
        #print(dims)
        numDims = len(dims)
        print("number of dimensions = " +str(numDims))
        print("---------")
        return dims

    #return the name of the dimensions
    def get_global_dim_names(self, dimensions):
        for key in dimensions:
            print("the name of the dimensions: "+ key)
        print("---------")

    #return the length of the dimensions
    def get_global_dim_length(self,dimensions):
        for key in dimensions:
            print("the length of dimension["+ key + "] = "+ str(len(dimensions[key])))
        print("---------")
        

        
    #return number of variables
    def nvars(self,dict):        
        vars = dict.variables
        nvar = len(vars)
        print("number of variables = " + str(nvar))
        print("---------")
        return vars

    #return the names of variables stored in the data set
    def get_var_names(self, variables):
        for key in variables:
            print("the name of variables: " + key )
        print("---------")

    #return the dimensions of the input variable 'var'
    def get_var_dim_lens(self,dims,vars,var, dim_name):
        print("---------variable '" + var + "'---------")
        print("shape = " +  str(vars[var].shape))
        vdims = vars[var].dimensions
        # for vd in vdims:
        #     print("dimension [" + vd + "] =" + str(len(dims[vd])))
        return len(dims[dim_name])
        
      
    #return the names of the dimensions for the variable 'var' 
    def get_var_dim_names(self,vars,var):
        vdims = vars[var].dimensions
        print("---------variable '" + var + "'---------")
        for key in vdims:
            print("dimension name: [" + key + "]" )
    
    #return the number of grid points that have the variable 'var' defined
    def get_num_of_points(self, dims, vars, var):
        grids = 1
        vdims = vars[var].dimensions
        for vd in vdims:
            grids = grids * len(dims[vd])
        print("number of grid points : " + str(grids))
        return grids
        
    #return the number of cells that have the variable 'var' defined 
    def get_num_of_Cells(self, ngrid, vars, var):
        vdims = vars[var].dimensions
        nvdims = len(vdims)
        cell = math.pow(2, nvdims)
        print("number of cells :" + str(ngrid//cell))
        
    #return the data array for the variable 'var'.
    def get_data(self, vars,var):
        #print("---------variable '" + var + "'---------")
        var_val = vars[var][:]
        #print var_val
        return var_val
            
        
    # plot the histogram of the variable 'var'.  
    def plot_histogram(self, vars, var):
        print("---------variable '" + var + "'---------")
        x = np.array(vars[var])
        #plt.hist(x, 10)
        #plt.show()

    #return data at computational point
    def data_at_comp_pos(self, temp_var, x, y):
        
        xlen = len(temp_var)
        ylen = len(temp_var[0][:])
        x0 = math.floor(x)
        y0 = math.floor(y)
        dx = x - x0
        dy = y - y0

        v1 = temp_var[x0][y0]
        v2 = temp_var[x0][y0]
        v3 = temp_var[x0][y0]
        v4 = temp_var[x0][y0]
        v = temp_var[x0][y0]
        
        
        if(x0 < xlen-1 and y0 < xlen-1):
            v1 = temp_var[x0][y0]
            v2 = temp_var[x0 + 1][y0]
            v3 = temp_var[x0 + 1][y0 + 1]
            v4 = temp_var[x0][y0 + 1]
            v12 = self.lerp(v1, v2, dx)
            v34 = self.lerp(v4, v3, dy)
            v = self.lerp(v12, v34,dy)
        elif (x0 == xlen - 1 and y0 != ylen -1):
            v1 = temp_var[x0][y0]
            v4 = temp_var[x0][y0 + 1]
            v = self.lerp(v1, v4, 0.5)
        elif (y0 == ylen - 1 and x0 != xlen -1):
            v1 = temp_var[x0][y0]
            v2 = temp_var[x0 + 1][y0]
            v = self.lerp(v1,v2,0.5)
        elif(x0 == xlen -1 and y0 == ylen -1):
            v = temp_var[x0][y0]
        return v

    #upsample: 360(lon) * 340(lat)
    #downsmaple:90(lon) * 85(lat)
    #input new lat and new lon
    def resample_tos(self, ratio, vars, time, resampled_var):
        temp_var = vars["tos"][:][time]
        len_lat = len(resampled_var[0])
        len_lon = len(resampled_var[0][0])
        #print len_lat, len_lon
        for y in xrange(0, len_lon):
            for x in xrange(0, len_lat):
                #calculate the source coordinates(u,v)
                u = x * (1.0/ratio)
                v = y * (1.0/ratio)
                value = self.data_at_comp_pos(temp_var, u, v)
                resampled_var[time][x][y] = value;
        
                
    #returen the data at a point with given lat and lon 
    def data_at_phys_pos(self, vars, var, lat, lon):
        lat_idx = int(lat - vars["lat"][0])  #calculate idx
        lon_idx = int((lon - vars["lon"][0])/2)
        #data_at_phy_pos(self, vars,  lon_idx, lat_idx)
        #print lat_idx, lon_idx, var[lon_idx][lat_idx]
        data = var[lon_idx][lat_idx]
        return data
            
    def lerp(self, v1, v2, ratio):
        p = v1 * (1 - ratio) + v2 * ratio
        return p

    #composite all 2d slices of 'var' between min and max along input 'dim'
    #using back_to_front formular
    def volume_composite(self,var, dim, dic, min_idx, max_idx, val, col_data):

        var_comp_len = max_idx - min_idx + 1
        var_len = len(var)
        var_width = len(var[0])
        var_height = len(var[0][0])

        var_comp = np.empty([var_comp_len, var_width, var_height])
        i = 0
        for idx in xrange(min_idx, max_idx + 1):
            var_comp[:][i] = var[:][idx]
            i += 1 
           
        val_max = np.nanmax(var_comp)
        val_min = np.nanmin(var_comp)

        colorConverter = cols.ColorConverter()
        if (dim == 'height'):
            
            for i in xrange(0, var_width):
                for j in xrange(0, var_height):
                    c = np.zeros(3)
                    
                    
                    if(dic < 0):                #along negative direction
                        for idx in xrange(0, var_comp_len):
                            xalpha = self.get_alpha_x(var_comp[idx][i][j], val_min, val_max, val) 
                            xnor = self.get_color_x(var_comp[idx][i][j], val_min, val_max)
                            
                            xcolor = colorConverter.to_rgb(str(xnor))
                            c[0] = c[0] * (1 - xalpha) + xcolor[0] * xalpha
                            c[1] = c[1] * (1 - xalpha) + xcolor[1] * xalpha
                            c[2] = c[2] * (1 - xalpha) + xcolor[2] * xalpha
                    else:                         #along positive direction
                        alpha = 0                           #front-to-back
                        for idx in xrange(0, var_comp_len):
                            xalpha = self.get_alpha_x(var_comp[idx][i][j], val_min, val_max, val) 
                            #float, gray 0-1
                            xnor = self.get_color_x(var_comp[idx][i][j], val_min, val_max)
                            
                            xcolor = colorConverter.to_rgb(repr(xnor))
                            #xcolor = colorConverter.to_rgb(repr(var_comp[idx][i][j]))
                            # c[0] = c[0] * (1 - xalpha) + xcolor[0] * xalpha
                            # c[1] = c[1] * (1 - xalpha) + xcolor[1] * xalpha
                            # c[2] = c[2] * (1 - xalpha) + xcolor[2] * xalpha

                            
                            c[0] = c[0] + xcolor[0] * xalpha * ( 1- alpha)
                            c[1] = c[1] + xcolor[1] * xalpha * ( 1- alpha)
                            c[2] = c[2] + xcolor[2] * xalpha * ( 1- alpha)
                            alpha = alpha + xalpha * ( 1 - alpha)
                            if(alpha >=0.9):
                                break

                            
                    col_data.append(c[0])
                    col_data.append(c[1])
                    col_data.append(c[2])

            #find the min and max value for transfer function
           
            
            #set matplotlib color map
            # jet = cm = plt.get_cmap('jet') 
            # cNorm  = colors.Normalize(vmin = val_min, vmax = val_max, clip=True)
            # scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
            
        return var_comp

    def get_alpha_x(self, var_val, min, max, opa):
        if(np.isnan(var_val)):
            return 1
        else:
            val = (var_val - min)/(max - min)
            if(val > opa):
                return (1-val)/(1-opa)
            else:
                return val/opa

    #return a normalized rgba array corresponding to var_val.
    def get_color_x(self, var_val, min, max):
        if(np.isnan(var_val)):
            return 1
        else:
            xnor = (var_val - min)/(max - min)
            #xcolor = cols.ColorConverter.to_rgb(xnor)
        return xnor;

    #return a 2D array for ploting 2D contour
    def contour2D(self, dimensions, var, dim, idx, iso_val):
        #get the correct plane for contouring
        height_len = len(dimensions['height'])
        lat_len = len(dimensions['lat'])
        lon_len = len(dimensions['lon'])
        
        if(dim == 'height'):
            temp_2d = np.empty([lat_len, lon_len])
            for i in xrange(0, lat_len):
                for j in xrange(0, lon_len):
                    temp_2d[i][j] = var[idx][i][j]
        elif(dim == 'lat'):
            temp_2d = np.empty([height_len, lon_len])
            for i in xrange(0, height_len):
                for j in xrange(0, lon_len):
                    temp_2d[i][j] = var[i][idx][j]
        elif(dim == 'lon'):
            temp_2d = np.empty([height_len, lat_len])
            for i in xrange(0, height_len):
                for j in xrange(0, lat_len):
                    temp_2d[i][j] = var[i][j][idx]
            
        #compute isocontour for the selected plane
        x_len = len(temp_2d[:])
        y_len = len(temp_2d[0][:])
               
        plt.figure()
        m = plt.imshow(temp_2d)
        plt.colorbar(m)
        x_data_total = []
        y_data_total = []
       
        for x in xrange(0, (x_len-1)-1):
            for y in xrange(0, (y_len-1)-1):
                x_data = []
                y_data = []
                if((temp_2d[x][y] < iso_val and temp_2d[x+1][y] > iso_val) or (temp_2d[x][y] > iso_val and temp_2d[x+1][y] < iso_val)):
                    p_x = (temp_2d[x][y] - iso_val)/(temp_2d[x][y] - temp_2d[x+1][y]) * (x +1 - x) + x
                    x_data.append(p_x)
                    y_data.append(y)
                                
                if((temp_2d[x+1][y] < iso_val and temp_2d[x+1][y+1] > iso_val) or (temp_2d[x+1][y] > iso_val and temp_2d[x+1][y+1] < iso_val)):
                    p_y = (temp_2d[x+1][y] - iso_val)/(temp_2d[x+1][y] - temp_2d[x+1][y+1]) * (y +1 - y) + y
                    x_data.append(x+1)
                    y_data.append(p_y)
                                
                if((temp_2d[x+1][y+1] < iso_val and temp_2d[x][y+1] > iso_val) or (temp_2d[x+1][y+1] > iso_val and temp_2d[x][y+1] < iso_val)):
                    p_x = (temp_2d[x][y+1] - iso_val)/(temp_2d[x][y+1] - temp_2d[x+1][y+1]) * (x +1 - x) + x
                    x_data.append(p_x)
                    y_data.append(y+1)
                                
                if((temp_2d[x][y] < iso_val and temp_2d[x][y+1] > iso_val) or ( temp_2d[x][y] >iso_val and temp_2d[x][y+1] < iso_val)):
                    p_y = (temp_2d[x][y] - iso_val)/(temp_2d[x][y] - temp_2d[x][y+1]) * (y +1 - y) + y
                    x_data.append(x)
                    y_data.append(p_y)
                
    
                x_data_total.append(x_data)
                y_data_total.append(y_data)
                
        
        
        plt.axis([0, x_len, 0, y_len])
        

        for i in xrange(0, len(x_data_total)):
            plt.plot(x_data_total[i], y_data_total[i])
        
        plt.savefig('out.png')
    
        
    def RK4(self, x, y, z, step_size, num_steps, u_data, v_data, w_data):
        x_trace = []
        y_trace = []
        z_trace = []
        x_trace.append(x)
        y_trace.append(y)
        z_trace.append(z)
        steps = 2 * step_size
        for nstep in xrange(0, num_steps):
            if(x >= 47 or y>= 47 or z>= 47):
                break
            ax = steps * u_data[x][y][z]
            newv = self.lerp_in_cube(x + ax/2, y, z, u_data)
            bx = steps * newv
            newv = self.lerp_in_cube(x + bx/2, y, z, u_data)
            cx  = steps * newv
            newv = self.lerp_in_cube(x + cx/2, y, z, u_data)
            dx = steps * newv
            pn_x = x + (ax + 2*bx + 2*cx + dx)/6 
            x_trace.append(pn_x)

            ay = steps * v_data[x][y][z]
            newv = self.lerp_in_cube(x, y + ay/2, z, v_data)
            by = steps * newv
            newv = self.lerp_in_cube(x, y + by/2, z, v_data)
            cy  = steps * newv
            newv = self.lerp_in_cube(x, y + cy/2, z, v_data)
            dy = steps * newv
            pn_y = y + (ay + 2*by + 2*cy + dy)/6 
            y_trace.append(pn_y)

            az = steps * w_data[x][y][z]
            newv = self.lerp_in_cube(x, y, z + az/2, w_data)
            bz = steps * newv
            newv = self.lerp_in_cube(x, y, z + bz/2, w_data)
            cz  = steps * newv
            newv = self.lerp_in_cube(x, y, z + cz/2, w_data)
            dz = steps * newv
            pn_z = z + (az + 2*bz + 2*cz + dz)/6 
            z_trace.append(pn_z)

            x = pn_x
            y = pn_y
            z = pn_z
            nstep += step_size 
        #print x_trace, y_trace, z_trace
        
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot(x_trace, y_trace, z_trace)
        #plt.savefig('rk42.png')
        plt.show()

    #interpolation 
    def lerp_in_cube(self,x,y,z,data):
            #find the most left point of the cell
            bx = math.trunc(x)
            by = math.trunc(y)
            bz = math.trunc(z)
            #alph, beta, gema
            a = x - bx
            b = y - by
            g = z - bz
            #finding the cell containing interpolated position
            v0 = data[bx][by][bz]
            v1 = data[bx+1][by][bz]
            v2 = data[bx][by][bz+1]
            v3 = data[bx+1][by][bz+1]
            v4 = data[bx][by+1][bz]
            v5 = data[bx+1][by+1][bz]
            v6 = data[bx][by+1][bz+1]
            v7 = data[bx+1][by+1][bz+1]
            
            p_val = (1-a)*(1-b)*(1-g)*v0 + a*(1-b)*(1-g)*v1 + (1-a)*b*(1-g)*v2 + a*b*(1-g)*v3 + (1-a)*(1-b)*g*v4 + a*(1-b)*g*v5 + (1-a)*b*g*v6 + a*b*g*v7
            
            return p_val

    #aprun -n num_procs python <program> <data file>  <xpartitions> <ypartitions> <zpartitions>
    def decompose(self, xpart, ypart, zpart, raw_data):
        size = MPI.COMM_WORLD.Get_size()
        rank = MPI.COMM_WORLD.Get_rank()
        #sys.stdout.write("Helloworld! I am process %d of %d\n" % (rank, size))

        xdim, ydim, zdim = unpack('3i', raw_data[0:12])
            

        data_size = xdim * ydim * zdim                
        if (rank == 0):
                    
            #x, y, z , x+dimX*(y+dimY*z)
            data = unpack('25000000f', raw_data[12:16 * data_size])
            data_array = np.asarray(data)

            print "mean from data_array: ", np.mean(data_array)

            xsub = int(math.ceil(float(xdim)/xpart))
            ysub = int(math.ceil(float(ydim)/ypart))
            zsub = int(math.ceil(float(zdim)/zpart))

            bound = np.zeros((size, 6), dtype = np.int)
           
            
            for ridx in xrange(0, size-1):
                
                zidx = (ridx) % zpart
                yidx = ((ridx) / zpart) % ypart
                xidx = (ridx) / (ypart * zpart)
              
                xmin = xsub * xidx
                xmax = xsub + xmin
                if(xmax > xdim - 1 ):
                    xmax = xdim -1
                ymin = ysub * yidx
                ymax = ysub + ymin
                if(ymax > ydim - 1):
                    ymax = ydim -1
                zmin = zsub * zidx
                zmax = zsub + zmin
                if(zmax > zdim - 1):
                    zmax = zdim -1

                if(xmin == 0):
                    bound[ridx+1][0] = xmin
                    bound[ridx+1][3] = xmax
                else:
                    bound[ridx+1][0] = xmin + 1
                    bound[ridx+1][3] = xmax
                if(ymin == 0):
                    bound[ridx+1][1] = ymin
                    bound[ridx+1][4] = ymax
                else:
                    bound[ridx+1][1] = ymin + 1
                    bound[ridx+1][4] = ymax
                if(zmin == 0):
                    bound[ridx+1][2] = zmin
                    bound[ridx+1][5] = zmax
                else:
                    bound[ridx+1][2] = zmin + 1
                    bound[ridx+1][5] = zmax

                #print "rank", ridx, "boundary:", bound[ridx][:]
           
        else:
            bound = None
            
        
            
        bound = MPI.COMM_WORLD.scatter(bound,root = 0)
        
        #bound[0],[3] -> x
        #bound[1],[4] -> y
        #boudn[2],[5] -> z
        if(rank != 0):
            print "Subvolume <", bound[0], bound[3], "> <" ,bound[1], bound[4], "> <", bound[2], bound[5],"> is assigned to process <", str(rank), ">"

        #--------------------------
        #broad cast to all other processes slice by slice along z
        local_buffer = np.empty(shape = (0), dtype = np.int)
        
        for i in xrange(0, zdim):
            
            if(rank == 0):
                
                sliced_data = np.zeros(250000)
                
                #x + dimX * (y + dimY * z)
                sidx_min = 0 + 500 * (0 + 500 * i)
                sidx_max = 499 + 500 * (499 + 500 * i)
                sliced_data = data_array[sidx_min:sidx_max]
                
            else:
                sliced_data = None

            #broadcasting.....
            sliced_data = MPI.COMM_WORLD.bcast(sliced_data,root = 0)
                      
            if(i >= bound[2] and i <= bound[5] and i != 0):
                local_buffer = np.append(local_buffer, self.copy_local_data(i,bound,sliced_data,rank))
        #total = np.zeros(1)
        total = np.zeros(size-1)
        #integral = np.zeros(1)
         
        if(rank != 0):
            #integral[0] = np.mean(local_buffer)
            #integral = np.mean(local_buffer)
            print "Process <", rank , "> has data <", bound[0], bound[3], "> <" ,bound[1], bound[4], "> <", bound[2], bound[5],">,  mean = <",  np.mean(local_buffer), ">"

        #reduce node receives results with a collective "reduce"
        #MPI.COMM_WORLD.Reduce(integral, total, op = MPI.SUM, root = 0)

        #gathering
        if(rank != 0)
            total = MPI.COMM_WORLD.gather(np.mean(local_buffer), root = 0)
        
        if rank == 0:
            print "Process", rank, "receives local means <", total[:][0], "> and the overall mean = <", np.mean(total), ">", 
            print "sum", np.sum(total), "mean / size -1", np.sum(total)/(size -1)

            
    def copy_local_data(self, z, bound, sliced_data, rank):
            
        min = bound[1] * 500 + bound[0]
        max = bound[4] * 500 + bound[3]

        return sliced_data[min:max]
        
            

            
            
            

                



    
            
                        
def main(argv):
    
    # rootgrp = Dataset(argv[1], 'r')
    # dims =rootgrp.dimensions
    reader = Netcdf_Reader()
    # dimensions = reader.ndims(dims)
    # reader.get_global_dim_names(dimensions)
    # reader.get_global_dim_length(dimensions)
    # variables = reader.nvars(rootgrp)
     
    #-------------2d contour--------------------------
    #var_temp = reader.get_data(variables, var)
    
    #reader.contour2D(dims,var_temp, dim, idx, val)
    #------------------------------------------------

    #------------RK_4-----------------------
    # u_data = reader.get_data(variables, 'U')
    # v_data = reader.get_data(variables, 'V')
    # w_data = reader.get_data(variables, 'W')
    

    # if(x >= 47 or y >= 47 or z >= 47):
    #     print "ERROR: the seed is out of bound"
    #     sys.exit()
    # else:
    #     reader.RK4(x, y, z, step_size, num_steps, u_data, v_data, w_data)
    #------------------------------------------------

    #------- domain decoposition & prallel mean caculation-------
    #aprun -n num_procs python <program> <data file>  <xpartitions> <ypartitions> <zpartitions>
    raw_data = open(sys.argv[1]).read()
    xpart = (int)(sys.argv[2])
    ypart = (int)(sys.argv[3])
    zpart = (int)(sys.argv[4])

    # #decoding binary data from file
    # xdim, ydim, zdim = unpack('3i', raw_data[0:12])
    # #print xdim, ydim, zdim

    # data_size = xdim * ydim * zdim
    # #x, y, z , x+dimX*(y+dimY*z)
    # data = unpack('25000000f', raw_data[12:16 * data_size])
    # data_array = np.asarray(data)

    #print "mean: ", np.mean(data_array)

    reader.decompose(xpart, ypart, zpart, raw_data)
    
    #-------round up
    # xsub = int(math.ceil(float(xdim)/xpart))
    # ysub = int(math.ceil(float(ydim)/ypart))
    # zsub = int(math.ceil(float(zdim)/zpart))
    
    # pro_size = xpart * ypart * zpart
    
    
    # data = np.zeros((pro_size, 6), dtype = np.int)
    # for rank in xrange(0, pro_size):
    #     print "rank: ", rank, "-----"
    #     zidx = rank % zpart
    #     yidx = (rank / zpart) % ypart
    #     xidx = rank / (ypart * zpart)
    #     print "(", xidx, yidx, zidx, ")"
    #     xmin = xsub * xidx
    #     xmax = xsub + xmin
    #     if(xmax > xdim):
    #         xmax = xdim
       
    #     ymin = ysub * yidx
    #     ymax = ysub + ymin
    #     if(ymax > ydim):
    #         ymax = ydim
       
    #     zmin = zsub * zidx
    #     zmax = zsub + zmin
    #     if(zmax > zdim):
    #         zmax = zdim
       
    #     if(xmin == 0):
    #         data[rank][0] = xmin
    #         data[rank][3] = xmax
    #     else:
    #         data[rank][0] = xmin + 1
    #         data[rank][3] = xmax
    #     if(ymin == 0):
    #         data[rank][1] = ymin
    #         data[rank][4] = ymax
    #     else:
    #         data[rank][1] = ymin + 1
    #         data[rank][4] = ymax
    #     if(zmin == 0):
    #         data[rank][2] = zmin
    #         data[rank][5] = zmax
    #     else:
    #         data[rank][2] = zmin + 1
    #         data[rank][5] = zmax
    #     print "x:", data[rank][0], data[rank][3]
    #     print "y:", data[rank][1], data[rank][4]
    #     print "z:", data[rank][2], data[rank][5]
       
if __name__ == '__main__':
    main(sys.argv)


