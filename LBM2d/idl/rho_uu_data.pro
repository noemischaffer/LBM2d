pro rho_uu_data


;rho=fltarr(7,7)
Nx=0
Ny=0
openr,1,'density_velocity.txt'
readf,1,Nx
readf,1,Ny
x=indgen(Nx+1)
y=indgen(Ny+1)
uux=fltarr(Nx,Ny,1)
uuy=fltarr(Nx,Ny,1)
readf,1,uux
readf,1,uuy
xaxis=x[1:Nx]
yaxis=y[1:Ny]
print,uux[Nx/2,*,0]
print, xaxis[*]
parabola=plot(xaxis[*], uux[Nx-3,*,0], xtitle='!8x', ytitle='!8u!D!8x', xrange=[1,Nx])
close,1
end

