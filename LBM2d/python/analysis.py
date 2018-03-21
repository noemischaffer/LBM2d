#!/opt/local/bin/python2.7
# Filename: analyse.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#------------------------------------------------#
def get_input_param():
  inputf=open("data/input.param","r")
  data=inputf.readlines()
  for line in data:
    words=line.split()
  Nxs=words[0]
  Nx=int(Nxs)
  Nys=words[1]
  Ny=int(Nys)
  Lxs=words[2]
  Lx=float(Lxs)
  Lys=words[3]
  Ly=float(Lys)
  vunits=words[4]
  vunit=float(vunits)
  taus=words[5]
  tau=float(taus)
  iTmaxs=words[6]
  iTMAX=int(iTmaxs)
  return Nx,Ny,Lx,Ly,vunit,tau,iTMAX
#------------------------------------------------#
def read_uu():
  Nx,Ny,Lx,Ly,vunit,tau,iTMAX=get_input_param()
  ux=np.zeros([Nx,Ny])
  uy=np.zeros([Nx,Ny])
  xx=np.linspace(0,Lx,Nx)
  yy=np.linspace(0,Ly,Ny)
#  xx=np.zeros([Nx,Ny])
#  yy=np.zeros([Nx,Ny])
  dx=Lx/(Nx-1)
  dy=Ly/(Ny-1)
  uu=np.loadtxt('data/usnap')
  #print uu.shape
  for iy in range(0,Ny):
    #print iy
    for ix in range(0,Nx):
      ieven=2*ix
      iodd=2*ix+1
      #print ix,ieven,iodd
      ux[ix,iy]=uu[iy,ieven]
      uy[ix,iy]=uu[iy,iodd]
#      xx[ix,iy]=ix*dx
#      yy[ix,iy]=iy*dy
#      ux[ix,iy]=xx[ix]*xx[ix]+yy[iy]
#      uy[ix,iy]=xx[ix]-yy[iy]*yy[iy]
  return Nx,Ny,Lx,Ly,ux,uy,xx,yy
#------------------------------------------------#
def read_domain():
  Nx,Ny,Lx,Ly,vunit,tau,iTMAX=get_input_param()
  solid=np.loadtxt('data/solid.now')
  xx=np.linspace(0,Lx,Nx)
  yy=np.linspace(0,Ly,Ny)
  plt.interactive(False)
  plt.imshow(solid)
  plt.axis('equal')
  plt.show()
  #for iy in range(0,Ny):
 #	for ix in range(0,Nx):
#	  plt.plot(xx[ix],yy[iy],'r.')
  return solid
#------------------------------------------------#
def streamline():
  Nx,Ny,Lx,Ly,ux,uy,xx,yy=read_uu()
  print(Nx,Ny)
  print(xx.shape,ux.shape)
  print(xx.shape,yy.shape)
# plt.interactive(False)
  plt.streamplot(yy,xx,ux,uy)
  plt.axis('equal')
  plt.show()
#  plt.plot(xx[:,4],ux[:,4],'.-')
#------------------------------------------------#
if __name__ == '__main__':
  print('standalone routine')
  #enkf(1,1)
else:
  print('Importing analysis module')
#------------------------------------------------#

