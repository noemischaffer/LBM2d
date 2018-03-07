pro rho_uu_data

u0 = 0.799154


;rho=fltarr(7,7)
Nx=0
Ny=0
openr,1,'/home/noemi/Documents/LBM/LBM2d/samples/shear_test/velocity.txt'
readf,1,Nx
readf,1,Ny
x=indgen(Nx+1)
y=indgen(Ny+1)
x_nat = indgen(12+1)
uux=fltarr(Nx,Ny,1)
uuy=fltarr(Nx,Ny,1)
readf,1,uux
readf,1,uuy
xaxis=x[1:Nx]
yaxis=y[1:Ny]
xaxis_nat=x_nat[1:12]

;u_theo = u0 * (1.0-(xaxis^2/6^2))

close, 1
parabola=plot(xaxis[*], uux[10,*,0], xtitle='!8x', ytitle='!8u!D!8x', xrange=[1,Nx], yrange = [0,1])
;parab_theo = plot(xaxis[*], u_theo, color = 'red', overplot = 1)
end

