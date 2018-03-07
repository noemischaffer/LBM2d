pro vorticity


device, decomposed=0
loadct, 5
!P.Color = '000000'xL ;Plot with black and have white background
!P.Background = 'FFFFFF'xL


Nx=0
Ny=0

ct = COLORTABLE(72, /reverse)


  openr,1,'/home/noemi/Documents/LBM/LBM2d/samples/shear_test/vorticity_idl.txt'
  readf,1,Nx
  readf,1,Ny
  x=indgen(Nx+1)
  y=indgen(Ny+1)
  uu_curl=fltarr(Nx,Ny)
  readf,1,uu_curl
  xaxis=x[1:Nx]
  yaxis=y[1:Ny]

set_plot, 'ps'

device, filename='/home/noemi/Documents/LBM/LBM2d/samples/shear_test/vorticity.eps', /encapsulated, bits_per_pixel=8, /color, xsize = 32, ysize = 32, /TIMES

close, 1

;plotimage, uu_curl[*,*], range= [min(uu_curl[*,*]), max(uu_curl[*,*])], xtitle='!8x', ytitle='!8y', /interp, ncolors= 100, xthick = 2, ythick = 2, font_name = 'times', font_size = 20, /interpolate, xtickformat = '(I3)', ytickformat = '(I3)', charsize=1.5
vort = contour(uu_curl[*,*],  xtitle='!8x', ytitle='!8y',  xthick = 2, ythick = 2, font_name = 'times', font_size = 20, xtickformat = '(I3)', ytickformat = '(I3)', RGB_TABLE=ct)


bar=colorbar(target=vort,  POSITION=[0.35,0.04,0.7,0.06])

device, /close


set_plot, 'x' 


end
