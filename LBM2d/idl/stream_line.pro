pro stream_line


device, decomposed=0
loadct, 5
!P.Color = '000000'xL ;Plot with black and have white background
!P.Background = 'FFFFFF'xL


Nx=0
Ny=0
  openr,1,'/home/noemi/Documents/LBM/LBM2d/samples/shear_test/velocity.txt'
  readf,1,Nx
  readf,1,Ny
  x=indgen(Nx)
  y=indgen(Ny)
  uu=fltarr(Nx,Ny)
  help, uu
  vv=fltarr(Nx,Ny)
  readf,1,uu
  readf,1,vv
 ; xaxis=x[1:Nx]
 ; yaxis=y[1:Ny]

set_plot, 'ps'

device, filename='/home/noemi/Documents/LBM/LBM2d/samples/shear_test/stream_lines.eps', /encapsulated, bits_per_pixel=8, /color, xsize = 32, ysize = 32, /TIMES

close, 1

stream=streamline(uu,vv, xtitle='x', ytitle='y', xrange=[0,128], yrange=[0,128])
; Change some properties.
stream.THICK = 2
stream.AUTO_COLOR = 1
stream.RGB_TABLE = 33

bar=colorbar(target=stream,  POSITION=[0.35,0.04,0.7,0.06])

device, /close


set_plot, 'x' 


end
