# IPython log file

import netCDF4 as nc
f=nc.Dataset('drifters.res.nc','w')
f.createDimension('i')
iv=f.createVariable('i','i4',('i'))
lonv=f.createVariable('lon','f8',('i'))
latv=f.createVariable('lat','f8',('i'))
depthv=f.createVariable('depth','f8',('i'))
timev=f.createVariable('time','f8',('i'))
idv=f.createVariable('id','i4',('i'))

lonv.long_name='longitude'
lonv.units='degrees_E'
latv.long_name='latitude'
latv.units='degrees_N'
depthv.long_name='depth'
depthv.units='meters'
timev.units='days since 1900-01-01 00:00:00'


iv[:]=[1,2,3,4]
lonv[:]=[5.,10.,15.,20.]
latv[:]=[40.,40.,40.,40.]
depthv[:]=[15.,15.,15.,15.]
idv[:]=[1,2,3,4]
timev[:]=[0.,0.,0.,0.]

f.sync()
f.close()

