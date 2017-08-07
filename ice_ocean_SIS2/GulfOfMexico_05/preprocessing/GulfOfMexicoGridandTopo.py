#This script requires the MIDAS package, in addition to netCDF4 and numpy
#MIDAS is located at https://github.com/mjharriso/MIDAS
#Anaconda installation guidelines:
#https://github.com/mjharriso/MIDAS/blob/master/.travis.yml


from midas.rectgrid import *
from midas.rectgrid_gen import *
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


# Generate approximately 1/24 degree supergrid (combined tracer and C-grid velocity grids)
#x=np.linspace(262.,279.5,441)
#y=np.linspace(18.,30.5,301)
ni=70
nj=50
x=np.linspace(262.,279.5,ni+1)
y=np.linspace(18.,30.5,nj+1)

X,Y=np.meshgrid(x,y)

# sgrid is a MIDAS supergrid class structure
sgrid=supergrid(xdat=X,ydat=Y,axis_units='degrees')
sgrid.grid_metrics()
model_grid=quadmesh(supergrid=sgrid)


# Read 30 sec GEBCO data (from http://www.gebco.net/data_and_products/gridded_bathymetry_data)

topo_path='GEBCO_2014_2D.nc'
#quadmesh is a MIDAS grid class
topo_grid=quadmesh(topo_path,var='elevation',simple_grid=True,cyclic=True)
ybnds=[sgrid.y.min()-0.1,sgrid.y.max()+0.1]
xbnds=[sgrid.x.min()-0.1,sgrid.x.max()+0.1]
ystart=np.where(topo_grid.latq>ybnds[0])[0][0]
yend=np.where(topo_grid.latq>ybnds[1])[0][0]
xstart=np.where(topo_grid.lonq>xbnds[0]-360.)[0][0]
xend=np.where(topo_grid.lonq>xbnds[1]-360.)[0][0]
region=topo_grid.indexed_region(j=(ystart,yend),i=(xstart,xend))
TOPO=state(topo_path,grid=topo_grid,geo_region=region,fields=['elevation'])


# Calculate a mean elevation of the 30-sec data on the model grid

R=TOPO.subtile('elevation',target=model_grid)



def show_depth(grid,depth):
    fig=plt.figure(figsize=(12,8))
    plt.pcolormesh(grid.x_T_bounds,grid.y_T_bounds,np.ma.masked_where(depth<=0.,depth),\
                   vmin=0,vmax=5000);
    plt.colorbar();
    plt.title('Model Bathymetry from GEBCO (m)')


#show_depth(model_grid,-sq(R.mean))
#plt.show()


def ice9it(i,j,depth,shallow=0.0):
  # Iterative implementation of "ice 9"
  wetMask = 0*depth
  (nj,ni) = wetMask.shape
  stack = set()
  stack.add( (j,i) )
  while stack:
    (j,i) = stack.pop()
    if wetMask[j,i] or depth[j,i] <= shallow: continue
    wetMask[j,i] = 1
    if i>0: stack.add( (j,i-1) )
    else: stack.add( (j,ni-1) )
    if i<ni-1: stack.add( (j,i+1) )
    else: stack.add( (0,j) )
    if j>0: stack.add( (j-1,i) )
    if j<nj-1: stack.add( (j+1,i) )
    else: stack.add( (j,ni-1-i) )
  return wetMask

def applyIce9(depth, i0, j0):

  numNotLand = np.count_nonzero(depth)
  print '# of wet points before Ice 9 = %i'%(numNotLand)
  notLand = ice9it(i0,j0,depth)
  newDepth = depth*notLand
  return newDepth


# Remove isolated cells

model_grid.D=-sq(R.mean)
model_grid.D[model_grid.D<0.0]=0.0
Depth=applyIce9(model_grid.D,np.round(ni/4),np.round(nj/4))

#show_depth(model_grid,Depth)
#plt.show()


model_grid.D=Depth
model_grid.wet[model_grid.D==0.]=0
model_grid.wet[model_grid.D>0.]=1

# Save interpolated topography to file

S=state(grid=model_grid)
vdict=R.var_dict['mean']
S.add_field_from_array(model_grid.D,'depth',var_dict=vdict)
S.write_nc('topog.nc',fields=['depth'])
f=nc.Dataset('topog.nc','a')
f.createDimension('ntiles',1)
f.close()


# Save supergrid to file

sgrid.write_nc('ocean_hgrid.nc')
f=nc.Dataset('ocean_hgrid.nc','a')
str=f.createDimension('string',255)
tile=f.createVariable('tile','S1',('string'))
tile[0:4]='tile1'
var=f.variables['tile']
dat = np.empty(1,'S'+repr(len(var)))
dat[0]='tile1'
dc=nc.stringtochar(dat)
var[:]=dc
f.close()


# Generate mosaic files for FMS exchange grid

name = 'ocean_mosaic'
rg = nc.Dataset(name+'.nc','w')
rg.createDimension('ntiles',1)
rg.createDimension('string',255)
mosaic = rg.createVariable('mosaic','c',('string',))
mosaic.standard_name = 'grid_mosaic_spec'
mosaic.children = 'contacts'
mosaic.grid_descriptor = ''
gridlocation = rg.createVariable('gridlocation','c',('string',))
gridlocation.standard_name = 'grid_file_location'
gridfiles = rg.createVariable('gridfiles','c',('ntiles','string',))
gridtiles = rg.createVariable('gridtiles','c',('ntiles','string',))
rg.grid_version = '0.2'
# Fill in data
mosaic[:] = '\000' * 255
mosaic[:12] = 'ocean_mosaic'
gridlocation[:] = '\000' * 255
gridlocation[:2] = './'
gridfiles[:] = '\000' * 255
gridfiles[0,:14] = 'ocean_hgrid.nc'
gridtiles[:] = '\000' * 255
gridtiles[0,:5] = 'tile1'
rg.close()

def set_string(variable, value):
    """Sets "variable" to "value" padded with blanks where
    "variable" is a netcdf variable object and "value" is a string."""
    variable[:] = '\000' * variable.shape[0]
    variable[:len(value)] = value

dx=nc.Dataset('ocean_hgrid.nc').variables['dx'][:]
dy=nc.Dataset('ocean_hgrid.nc').variables['dy'][:]
d2x=dx+numpy.roll(dx,shift=-1,axis=1)
d2x=d2x[:,::2]
DX=0.5*(d2x+numpy.roll(d2x,shift=-1,axis=0))
DX=DX[:-1:2,:]
d2y=dy+numpy.roll(dy,shift=-1,axis=0)
d2y=d2y[::2,:]
DY=0.5*(d2y+numpy.roll(d2y,shift=-1,axis=1))
DY=DY[:,:-1:2]


nj=model_grid.jm;ni=model_grid.im
snj=2*nj;sni=2*ni
Ocean_Depth=nc.Dataset('topog.nc').variables['depth'][:]
nl=len(numpy.where(Ocean_Depth==0.)[0])
print 'Number of land points= ',nl
AREA=DX*DY
rg = nc.Dataset('atmos_mosaic_tile1Xland_mosaic_tile1.nc','w',format='NETCDF3_CLASSIC') # atmos_mosaic_tile1Xland_mosaic_tile1.nc
rg2 = nc.Dataset('land_mask.nc','w',format='NETCDF3_CLASSIC') # atmos_mosaic_tile1Xland_mosaic_tile1.nc
rg.createDimension('string',255)
rg.createDimension('ncells',nl)  #It is unclear whether this works when nl=0. It does work for nl>0
rg.createDimension('two',2)
contact = rg.createVariable('contact','c',('string',))
contact.standard_name = 'grid_contact_spec'
contact.contact_type = 'exchange'
contact.parent1_cell = 'tile1_cell'
contact.parent2_cell = 'tile2_cell'
contact.xgrid_area_field = 'xgrid_area'
contact.distant_to_parent1_centroid = 'tile1_distance'
contact.distant_to_parent2_centroid = 'tile2_distance'
tile1_cell = rg.createVariable('tile1_cell','i4',('ncells','two',))
tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
tile2_cell = rg.createVariable('tile2_cell','i4',('ncells','two',))
tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
xgrid_area = rg.createVariable('xgrid_area','f8',('ncells',))
xgrid_area.standard_name = 'exchange_grid_area'
xgrid_area.units = 'm2'
tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
rg.grid_version = '0.2'

rg2.createDimension('nx',ni)
rg2.createDimension('ny',nj)
mask=rg2.createVariable('mask','f8',('ny','nx'))
mask.standard_name  = 'land fraction at T-cell centers'
mask.units = 'none'
mask[:,:]=0.0
rg2.grid_version = '0.2'

contact[:] = '\000' * 255
contact[:37] = 'atmos_mosaic:tile1::land_mosaic:tile1'
count=0
for j in range(nj):
    for i in range(ni):
        if Ocean_Depth[j,i]==0.:
            tile1_cell[count] = [i+1,j+1]
            tile2_cell[count] = [i+1,j+1]
            tile1_distance[count] = [0.,0.]
            tile2_distance[count] = [0.,0.]
            xgrid_area[count] = AREA[j,i]
            count=count+1
            mask[j,i]=1.0
rg.close()
rg2.close()

rg = nc.Dataset('atmos_mosaic_tile1Xocean_mosaic_tile1.nc','w',format='NETCDF3_CLASSIC') # atmos_mosaic_tile1Xocean_mosaic_tile1.nc
rg2 = nc.Dataset('ocean_mask.nc','w',format='NETCDF3_CLASSIC') # atmos_mosaic_tile1Xland_mosaic_tile1.nc
rg.createDimension('string',255)
rg.createDimension('ncells',ni*nj-nl) # -1 is for a single land point
print 'ncells= ',ni*nj-nl
rg.createDimension('two',2)
contact = rg.createVariable('contact','c',('string',))
contact.standard_name = 'grid_contact_spec'
contact.contact_type = 'exchange'
contact.parent1_cell = 'tile1_cell'
contact.parent2_cell = 'tile2_cell'
contact.xgrid_area_field = 'xgrid_area'
contact.distant_to_parent1_centroid = 'tile1_distance'
contact.distant_to_parent2_centroid = 'tile2_distance'
tile1_cell = rg.createVariable('tile1_cell','i4',('ncells','two',))
tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
tile2_cell = rg.createVariable('tile2_cell','i4',('ncells','two',))
tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
xgrid_area = rg.createVariable('xgrid_area','f8',('ncells',))
xgrid_area.standard_name = 'exchange_grid_area'
xgrid_area.units = 'm2'
tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
rg.grid_version = '0.2'
# Fill in data
contact[:] = '\000' * 255
contact[:38] = 'atmos_mosaic:tile1::ocean_mosaic:tile1'

rg2.createDimension('nx',ni)
rg2.createDimension('ny',nj)
mask=rg2.createVariable('mask','f8',('ny','nx'))
mask.standard_name  = 'ocean fraction at T-cell centers'
mask.units = 'none'
mask[:,:]=0.0
rg2.grid_version = '0.2'

count=0
for j in range(nj):
    for i in range(ni):
        if Ocean_Depth[j,i]!=0:
            tile1_cell[count] = [i+1,j+1]
            tile2_cell[count] = [i+1,j+1]
            tile1_distance[count] = [0.,0.]
            tile2_distance[count] = [0.,0.]
            xgrid_area[count] = AREA[j,i]
            count=count+1
            mask[j,i]=1.0
rg.close()
rg2.close()

rg = nc.Dataset('land_mosaic_tile1Xocean_mosaic_tile1.nc','w',format='NETCDF3_CLASSIC') # land_mosaic_tile1Xocean_mosaic_tile1.nc
rg.createDimension('string',255)
rg.createDimension('ncells',ni*nj-nl) # -1 is for a single land point
rg.createDimension('two',2)
contact = rg.createVariable('contact','c',('string',))
contact.standard_name = 'grid_contact_spec'
contact.contact_type = 'exchange'
contact.parent1_cell = 'tile1_cell'
contact.parent2_cell = 'tile2_cell'
contact.xgrid_area_field = 'xgrid_area'
contact.distant_to_parent1_centroid = 'tile1_distance'
contact.distant_to_parent2_centroid = 'tile2_distance'
tile1_cell = rg.createVariable('tile1_cell','i4',('ncells','two',))
tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
tile2_cell = rg.createVariable('tile2_cell','i4',('ncells','two',))
tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
xgrid_area = rg.createVariable('xgrid_area','f8',('ncells',))
xgrid_area.standard_name = 'exchange_grid_area'
xgrid_area.units = 'm2'
tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
rg.grid_version = '0.2'
# Fill in data
contact[:] = '\000' * 255
contact[:37] = 'land_mosaic:tile1::ocean_mosaic:tile1'
count=0
for j in range(nj):
    for i in range(ni):
        if Ocean_Depth[j,i]>0.:
            tile1_cell[count] = [i+1,j+1]
            tile2_cell[count] = [i+1,j+1]
            tile1_distance[count] = [0.,0.]
            tile2_distance[count] = [0.,0.]
            xgrid_area[count] = AREA[j,i]
            count=count+1
            
            
rg.close()

g = nc.Dataset('grid_spec.nc','w',format='NETCDF3_CLASSIC') # land_mosaic_tile1Xocean_mosaic_tile1.nc
rg.createDimension('string',255)
rg.createDimension('nfile_aXo',1) # -1 is for a single land point
rg.createDimension('nfile_aXl',1) # -1 is for a single land point
rg.createDimension('nfile_lXo',1) # -1 is for a single land point
atm_mosaic_dir = rg.createVariable('atm_mosaic_dir','c',('string',))
atm_mosaic_dir.standard_name = 'directory_storing_atmosphere_mosaic'
atm_mosaic_file = rg.createVariable('atm_mosaic_file','c',('string',))
atm_mosaic_file.standard_name = 'atmosphere_mosaic_file_name'
atm_mosaic = rg.createVariable('atm_mosaic','c',('string',))
atm_mosaic.standard_name = 'atmosphere_mosaic_name'
lnd_mosaic_dir = rg.createVariable('lnd_mosaic_dir','c',('string',))
lnd_mosaic_dir.standard_name = 'directory_storing_land_mosaic'
lnd_mosaic_file = rg.createVariable('lnd_mosaic_file','c',('string',))
lnd_mosaic_file.standard_name = 'land_mosaic_file_name'
lnd_mosaic = rg.createVariable('lnd_mosaic','c',('string',))
lnd_mosaic.standard_name = 'land_mosaic_name'
ocn_mosaic_dir = rg.createVariable('ocn_mosaic_dir','c',('string',))
ocn_mosaic_dir.standard_name = 'directory_storing_ocean_mosaic'
ocn_mosaic_file = rg.createVariable('ocn_mosaic_file','c',('string',))
ocn_mosaic_file.standard_name = 'ocean_mosaic_file_name'
ocn_mosaic = rg.createVariable('ocn_mosaic','c',('string',))
ocn_mosaic.standard_name = 'ocean_mosaic_name'
ocn_topog_dir = rg.createVariable('ocn_topog_dir','c',('string',))
ocn_mosaic_dir.standard_name = 'directory_storing_ocean_topog'
ocn_topog_file = rg.createVariable('ocn_topog_file','c',('string',))
ocn_topog_file.standard_name = 'ocean_topog_file_name'
aXo_file = rg.createVariable('aXo_file','c',('nfile_aXo','string',))
aXo_file.standard_name = 'atmXocn_exchange_grid_file'
aXl_file = rg.createVariable('aXl_file','c',('nfile_aXl','string',))
aXl_file.standard_name = 'atmXlnd_exchange_grid_file'
lXo_file = rg.createVariable('lXo_file','c',('nfile_lXo','string',))
lXo_file.standard_name = 'lndXocn_exchange_grid_file'
#Global attributes
rg.grid_version = '0.2'
rg.code_version = "$Name:  $"
rg.history = " "

atm_mosaic_dir[:] = '\000' * 255
atm_mosaic_dir[:2] = './'
atm_mosaic_file[:] = '\000' * 255
atm_mosaic_file[:15] = 'ocean_mosaic.nc'
atm_mosaic[:] = '\000' * 255
atm_mosaic[:12] = 'atmos_mosaic'
lnd_mosaic_dir[:] = '\000' * 255
lnd_mosaic_dir[:2] = './'
lnd_mosaic_file[:] = '\000' * 255
lnd_mosaic_file[:15] = 'ocean_mosaic.nc'
lnd_mosaic[:] = '\000' * 255
lnd_mosaic[:11] = 'land_mosaic'
ocn_mosaic_dir[:] = '\000' * 255
ocn_mosaic_dir[:2] = './'
ocn_mosaic_file[:] = '\000' * 255
ocn_mosaic_file[:15] = 'ocean_mosaic.nc'
ocn_mosaic[:] = '\000' * 255
ocn_mosaic[:12] = 'ocean_mosaic'
ocn_topog_dir[:] = '\000' * 255
ocn_topog_dir[:2] = './'
ocn_topog_file[:] = '\000' * 255
ocn_topog_file[:8] = 'topog.nc'
aXo_file[:,:] = '\000' * 255
aXo_file[:,:40] = 'atmos_mosaic_tile1Xocean_mosaic_tile1.nc'
aXl_file[:,:] = '\000' * 255
aXl_file[:,:39] = 'atmos_mosaic_tile1Xland_mosaic_tile1.nc'
lXo_file[:,:] = '\000' * 255
lXo_file[:,:39] = 'land_mosaic_tile1Xocean_mosaic_tile1.nc'

rg.close()
                                                                                                
