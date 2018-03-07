pro solid_points


device, decomposed=0
loadct, 5
!P.Color = '000000'xL ;Plot with black and have white background
!P.Background = 'FFFFFF'xL


Nx=0
Ny=0
openr,1,'/home/noemi/Documents/LBM/LBM2d/samples/shear_test/solid_points.txt'

readf,1,Nx
readf,1,Ny
i=indgen(1,Nx*Ny)
j=indgen(1,Nx*Ny)
is_solid=indgen(1, Nx*Ny)
readf,1, is_solid
readf,1, i
readf,1, j
set_plot, 'ps'

device, filename='/home/noemi/Documents/LBM/LBM2d/samples/shear_test/solid_points.eps', /encapsulated, bits_per_pixel=8, /color, xsize = 32, ysize = 32, /TIMES

close, 1

plot, i,j, psym=2, xrange=[1,Nx], yrange=[1,Ny]
;plot, i,j, psym=3



save, Nx, Ny,i, j, filename='solid_cooordinates.dat'

device, /close

set_plot, 'x'

end
