pro surface_points


device, decomposed=0
loadct, 5
!P.Color = '000000'xL ;Plot with black and have white background
!P.Background = 'FFFFFF'xL


openr,1,'/home/noemi/Documents/LBM/LBM2d/samples/shear_test/surface.txt'



Nsurf=0
readf,1, Nsurf


x=indgen(Nsurf)
y=indgen(Nsurf)


  readf,1,x
  readf,1,y

set_plot, 'ps'

device, filename='/home/noemi/Documents/LBM/LBM2d/samples/shear_test/surface.eps', /encapsulated, bits_per_pixel=8, /color, xsize = 32, ysize = 32, /TIMES

restore, filename='solid_cooordinates.dat'

print, Nx

plot,x,y,psym=2, xrange=[1,Nx], yrange=[1,Ny]
oplot, i, j, color=100, psym=5


device, /close

close, 1

set_plot, 'x'

end
