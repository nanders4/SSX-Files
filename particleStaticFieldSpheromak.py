import numpy as np
import netCDF4 as n4
import math
import h5py
from scipy.special import j0, j1, jn_zeros

'''This script generates a magnetic field superimposed on an orbit in VAPOR.
Though specific to a spheromak field here, it is easy to adapt for any function that
returns magnetic field. Simply replace 'getSpheromakFieldAtPosition' with the desired
function (must take x,y,z coordinates and return Bx,By,Bz), and also change the call.'''

def getSpheromakFieldAtPosition(x, y, z, center=(0,0,1), B0=1, R=1, L=1):
   """The spheromak center in z must be L.
   """

   # parameters
   j1_zero1 = jn_zeros(1,1)[0]
   kr = j1_zero1/R
   kz = np.pi/L

   lam = np.sqrt(kr**2 + kz**2)

   # construct cylindrical coordinates centered on center
   r = np.sqrt((x- center[0])**2 + (y- center[1])**2)
   theta = np.arctan2(y,x)
   centZ = z - center[2]


   # calculate cylindrical fields
   Br = -B0 * kz/kr * j1(kr*r) * np.cos(kz*centZ)
   Bt = B0 * lam/kr * j1(kr*r) * np.sin(kz*centZ)

   # convert back to cartesian, place on grid.
   Bx = Br*np.cos(theta) - Bt*np.sin(theta)
   By = Br*np.sin(theta) + Bt*np.cos(theta)
   Bz = B0 * j0(kr*r) * np.sin(kz*centZ)
   return Bx, By, Bz

def getExtrema(dataArray):
    max = np.NINF
    min = np.Inf
    for i in range(len(dataArray)):
        if(dataArray[i]>max):
            max = dataArray[i]
        elif(dataArray[i]<min):
            min = dataArray[i]
    return max, min

def dataToMesh(dataArray, min, max, meshPoints):
    dimSize = max-min
    outArray = np.zeros(len(dataArray))
    for i in range(len(dataArray)):
        outArray[i] = int(round((meshPoints-1)*(dataArray[i]-min)/dimSize))
    return outArray


file = h5py.File('analytical_data.h5', 'r')
data = file['r'][10]
print("Data shape: ", data.shape)
#get data into proper format
#x y z t

xs = data[:,0]
ys = data[:,1]
zs = data[:,2]

#make 3D grid from data
timeCompresssion = 1
totalTimes = data.shape[0]

nx = 250
ny = 250
nz = 250
nt = int(round(totalTimes*timeCompresssion))

xMin = -100
yMin = -100
zMin = 0
xMax = 100
yMax = 100
zMax = 100

Bx = np.zeros((nx, ny, nz))
By = np.zeros((nx, ny, nz))
Bz = np.zeros((nx, ny, nz))
for i in range(nx):
    x=(xMax-xMin)*(i/nx)+xMin
    for j in range(ny):
        y=(yMax-yMin)*(j/ny)+yMin
        if(np.abs(y)<.000001):
            y=.000001
        for k in range(nz):
            z=(zMax-zMin)*(k/nz)+zMin
            if(np.abs(z)<.000001):
                z=.000001
            r = np.array([x,y,z])
            Barray = getSpheromakFieldAtPosition(x,y,z, center=(0,0,100), B0=60.256574486577478, R=100, L=100)
            #print(Barray)
            if(math.isnan(Barray[0])==False):
                Bx[i][j][k] = Barray[0]
            if(math.isnan(Barray[1])==False):
                By[i][j][k] = Barray[1]
            if(math.isnan(Barray[2])==False):
                Bz[i][j][k] = Barray[2]

xMaxP, xMinP = getExtrema(xs)
yMaxP, yMinP = getExtrema(ys)
zMaxP, zMinP = getExtrema(zs)

xInds = dataToMesh(xs, xMinP, xMaxP, nx)
yInds = dataToMesh(ys, yMinP, yMaxP, ny)
zInds = dataToMesh(zs, zMinP, zMaxP, nz)
partMesh = np.zeros((nx, ny, nz))

for i in range(nt):
    index = int(round((totalTimes)*(i/nt)))
    partMesh[int(xInds[index])-1][int(yInds[index])-1][int(zInds[index])-1] = 10

#for netCDF
x = np.linspace(-1,1,nx)
y = np.linspace(-1,1,ny)
z = np.linspace(0,1,nz)

dataset = n4.Dataset('particle.nc', 'w', format='NETCDF4')

xset = dataset.createDimension('x', nx)
yset = dataset.createDimension('y', ny)
zset = dataset.createDimension('z', nz)

xs = dataset.createVariable('x', np.float64, ('x',))
ys = dataset.createVariable('y', np.float64, ('y',))
zs = dataset.createVariable('z', np.float64, ('z',))

ps = dataset.createVariable('p', np.float64, ('x','y','z'))
Bxs = dataset.createVariable('Bx', np.float64, ('x','y','z'))
Bys = dataset.createVariable('By', np.float64, ('x','y','z'))
Bzs = dataset.createVariable('Bz', np.float64, ('x','y','z'))

xs[:] = x
ys[:] = y
zs[:] = z
ps[:] = partMesh
Bxs[:] = Bx
Bys[:] = By
Bzs[:] = Bz

dataset.close()
