import argparse
import netCDF4
import numpy
import scipy.interpolate

parser = argparse.ArgumentParser()
parser.add_argument('input_restart_file', type=str, help='Input restart file')
parser.add_argument('ni', type=int, help='Size in i-direction')
parser.add_argument('nj', type=int, help='Size in j-direction')

args = parser.parse_args()

def copy_var(old, new, varname):
    oldvar = old.variables[varname]
    var = new.createVariable(varname, oldvar.dtype, oldvar.dimensions)
    for att in oldvar.ncattrs():
        if not att in ['checksum']:
            var.setncattr(att, oldvar.__getattr__(att))
    return var
def interp(old, new, xh, yh, lon, lat, varname, var):
    nj,ni = old.variables[varname][0,0].shape
    X,Y,Z = numpy.zeros((ni+2)),numpy.zeros((nj+2)),numpy.zeros((nj+2,ni+2))
    X[1:-1] = xh
    X[0] = 2.*xh[0] - xh[1]
    X[-1] = 2.*xh[-1] - xh[-2]
    Y[1:-1] = yh
    Y[0] = 2.*yh[0] - yh[1]
    Y[-1] = 2.*yh[-1] - yh[-2]
    for k in range(old.variables[varname].shape[1]):
        Z[1:-1,1:-1] = old.variables[varname][0,k][:]
        Z[:,0] = Z[:,-2]
        Z[:,-1] = Z[:,1]
        Z[0,:] = Z[-2,:]
        Z[-1,:] = Z[1,:]
        f = scipy.interpolate.interp2d(X, Y, Z)
        var[0,k] = f( lon[:], lat[:] )

old = netCDF4.Dataset(args.input_restart_file, 'r')

ni,nj = args.ni,args.nj
with netCDF4.Dataset('MOM.res.nc','w') as new:
    yq = old.variables['latq'][:]
    xq = old.variables['lonq'][:]
    yh = old.variables['lath'][:]
    xh = old.variables['lonh'][:]
    lenlon = xq[-1] - xq[0]
    lenlat = yq[-1] - yq[0]

    new.createDimension('Time',None)
    new.createDimension('lath',nj)
    new.createDimension('lonh',ni)
    new.createDimension('latq',nj+1)
    new.createDimension('lonq',ni+1)
    new.createDimension('Layer',2)
    Time = copy_var(old, new, 'Time')
    lath = copy_var(old, new, 'lath')
    lonh = copy_var(old, new, 'lonh')
    latq = copy_var(old, new, 'latq')
    lonq = copy_var(old, new, 'lonq')
    layer = copy_var(old, new, 'Layer')
    latq[:] = yq[0] + (numpy.arange(nj+1)/nj)*lenlat
    lonq[:] = xq[0] + (numpy.arange(ni+1)/ni)*lenlon
    lath[:] = yq[0] + ((numpy.arange(nj)+0.5)/nj)*lenlat
    lonh[:] = xq[0] + ((numpy.arange(ni)+0.5)/ni)*lenlon
    layer[:] = old.variables['Layer'][:]
    Time[0] = old.variables['Time'][0]
    h = copy_var(old, new, 'h')
    interp(old, new, xh, yh, lonh[:], lath[:], 'h', h)
    u = copy_var(old, new, 'u')
    interp(old, new, xq, yh, lonq[:], lath[:], 'u', u)
    v = copy_var(old, new, 'v')
    interp(old, new, xh, yq, lonh[:], latq[:], 'v', v)

