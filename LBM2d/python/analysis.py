#!/opt/local/bin/python2.7
# Filename: analyse.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
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
  taus=words[4]
  tau=float(taus)
  iTmaxs=words[5]
  iTMAX=int(iTmaxs)
  return Nx,Ny,Lx,Ly,tau,iTMAX
#------------------------------------------------#
def read_rho():
  Nx,Ny,Lx,Ly,tau,iTMAX=get_input_param()
  xx=np.linspace(0,Lx,Nx)
  yy=np.linspace(0,Ly,Ny)
  dx=1.
  dy=1.
  rho=np.loadtxt('data/rho_snap')
  return Nx,Ny,Lx,Ly,xx,yy,rho
#------------------------------------------------#
def read_uu():
  Nx,Ny,Lx,Ly,tau,iTMAX=get_input_param()
  ux=np.zeros([Nx+1,Ny+1])
  uy=np.zeros([Nx+1,Ny+1])
  xx=np.linspace(0,Lx,Nx+1)
  yy=np.linspace(0,Ly,Ny+1)
  dx=1.
  dy=1.
  uu=np.loadtxt('data/usnap')
  #print uu.shape
  for iy in range(0,Ny+1):
    #print iy
    for ix in range(0,Nx+1):
      ieven=2*iy
      iodd=2*iy+1
      #print ix,ieven,iodd
      ux[ix,iy]=uu[ix,ieven]
      uy[ix,iy]=uu[ix,iodd]
  return Nx,Ny,Lx,Ly,ux,uy,xx,yy
#------------------------------------------------#
def read_domain(lshow):
  Nx,Ny,Lx,Ly,tau,iTMAX=get_input_param()
  solid=np.loadtxt('data/solid.now')
  xx=np.linspace(0,Lx,Nx+1)
  yy=np.linspace(0,Ly,Ny+1)
  plt.interactive(False)
  im=plt.imshow(solid,interpolation='none',cmap=cm.Greys)
  plt.axis('equal')
  if lshow == 1 :
    plt.show()
  return solid
#------------------------------------------------#
def streamline():
  Nx,Ny,Lx,Ly,ux,uy,xx,yy=read_uu()
  plt.close()
  plt.interactive(False)
  solid=read_domain(0)
  plt.streamplot(xx,yy,np.transpose(ux),np.transpose(uy))
  plt.axis('tight')
  plt.axis('equal')
  plt.show()
  return ux,uy
#  plt.plot(xx[:,4],ux[:,4],'.-')
#------------------------------------------------#
def pts():
  ts=np.loadtxt('data/ts')
  plt.close()
  plt.interactive(False)
  plt.plot(ts[:,0],ts[:,1],'.')
  plt.plot(ts[:,0],ts[:,2],'s')
  plt.show()
  return ts
#------------------------------------------------#
if __name__ == '__main__':
  print('standalone routine')
  #enkf(1,1)
else:
  print('Importing analysis module')
#------------------------------------------------#

