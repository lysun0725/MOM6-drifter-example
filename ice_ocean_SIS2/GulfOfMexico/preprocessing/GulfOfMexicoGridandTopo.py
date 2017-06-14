
# coding: utf-8

# In[48]:

#get_ipython().magic(u'pylab inline')
from midas.rectgrid import *
from midas.rectgrid_gen import *
import numpy as np
from IPython.display import Image
from IPython.core.display import HTML 


# In[49]:

import matplotlib.pyplot as plt


# In[60]:

x=np.linspace(262.,279.5,441)
y=np.linspace(18.,30.5,301)
X,Y=np.meshgrid(x,y)


# In[61]:

sgrid=supergrid(xdat=X,ydat=Y,axis_units='degrees')
sgrid.grid_metrics()
model_grid=quadmesh(supergrid=sgrid)


# In[62]:

topo_path='GEBCO_2014_2D.nc'
topo_grid=quadmesh(topo_path,var='elevation',simple_grid=True,cyclic=True)
ybnds=[sgrid.y.min()-0.1,sgrid.y.max()+0.1]
xbnds=[sgrid.x.min()-0.1,sgrid.x.max()+0.1]
ystart=np.where(topo_grid.latq>ybnds[0])[0][0]
yend=np.where(topo_grid.latq>ybnds[1])[0][0]
xstart=np.where(topo_grid.lonq>xbnds[0]-360.)[0][0]
xend=np.where(topo_grid.lonq>xbnds[1]-360.)[0][0]
region=topo_grid.indexed_region(j=(ystart,yend),i=(xstart,xend))
TOPO=state(topo_path,grid=topo_grid,geo_region=region,fields=['elevation'])


# In[63]:

R=TOPO.subtile('elevation',target=model_grid)


# In[76]:

def show_depth(grid,depth):
    fig=plt.figure(figsize=(12,8))
    plt.pcolormesh(grid.x_T_bounds,grid.y_T_bounds,np.ma.masked_where(depth<=0.,depth),vmin=0,vmax=5000);
    plt.colorbar();
    plt.title('Model Bathymetry from GEBCO (m)')


# In[77]:

show_depth(model_grid,-sq(R.mean))


# In[129]:

def ice9it(i,j,depth,shallow=0.0):
  # Iterative implementation of "ice 9"
  wetMask = 0*depth
  (nj,ni) = wetMask.shape
  stack = set()
  stack.add( (j,i) )
  while stack:
    (j,i) = stack.pop()
#    print i,j,depth[j,i]
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

def applyIce9(depth, i0,j0):

  numNotLand = np.count_nonzero(depth)
  print '# of wet points before Ice 9 = %i'%(numNotLand)
  notLand = ice9it(i0,j0,depth)
  newDepth = depth*notLand
  return newDepth


# In[127]:

model_grid.D=-sq(R.mean)
model_grid.D[model_grid.D<0.0]=0.0
Depth=applyIce9(model_grid.D,110,75)


# In[128]:

show_depth(model_grid,Depth)


# In[130]:

model_grid.D=Depth
model_grid.wet[model_grid.D==0.]=0
model_grid.wet[model_grid.D>0.]=1

S=state(grid=model_grid)
vdict=R.var_dict['mean']
S.add_field_from_array(model_grid.D,'depth',var_dict=vdict)

S.write_nc('topog.nc',fields=['depth'])

sgrid.write_nc('ocean_hgrid.nc')
