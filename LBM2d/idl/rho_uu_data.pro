pro rho_uu_data

Nx=0
Ny=0
openr,1,'/home/noemi/Documents/LBM_Runs/poiseuille/velocity.txt'
readf,1,Nx
readf,1,Ny
print, Nx
print, Ny
x=indgen(Nx+3)
y=indgen(Ny+3)
uux=fltarr(Nx+2,Ny+2,1)
uuy=fltarr(Nx+2,Ny+2,1)
readf,1,uux
readf,1,uuy
xaxis=x[0:Nx+2]
;print, xaxis
yaxis=y[0:Ny+2]
print, uux
;u_theo = u0 * (1.0-(xaxis^2/6^2))

close, 1
parabola=plot(xaxis[*], uux[10,*,0], xtitle='!8x', ytitle='!8u!D!8x', yrange = [0,max(uux)])
;parab_theo = plot(xaxis[*], u_theo, color = 'red', overplot = 1)
end

