pro stream_line


device, decomposed=0
loadct, 5
!P.Color = '000000'xL ;Plot with black and have white background
!P.Background = 'FFFFFF'xL


Nx=0
Ny=0

;for n=0,3001,200 do begin
;  if (n eq 0) then begin
;  endif else begin
    ;openr,Unit,'/home/noemi/Documents/LBM_Runs/DiskInGravity/velocity_Re'+strcompress(string(n,  Format='(I05)'), /remove_all) +'.txt', /GET_LUN
    openr,Unit,'/home/noemi/Documents/LBM_Runs/DiskInGravity/velocity.txt', /GET_LUN
;    print, 'opened file '+n 
    readf,Unit,Nx
    readf,Unit,Ny
    x=indgen(Nx+2)
    y=indgen(Ny+2)
    uu=fltarr(Nx+2,Ny+2)
    vv=fltarr(Nx+2,Ny+2)
    readf,Unit,uu
    readf,Unit,vv

    close, Unit
 
    stream=streamline(uu,vv, xtitle='x', ytitle='y', xrange=[0,Nx], yrange=[0,Ny]) ;,title = 'timestep= ' + string(n, f='(I04)'))
    ; Change some properties.
    stream.THICK = 2
    stream.AUTO_COLOR = 1
    stream.RGB_TABLE = 33

    bar=colorbar(target=stream,  POSITION=[0.35,0.04,0.7,0.06])
    stream.save, strcompress(string(n,  Format='(I05)'), /remove_all) + '_velocity_Re.png'
    stream.close
 ; endelse
;endfor
; make movie 
;mencoder "mf://*.png" -mf fps=10 -o movie.mp4 -ovc lavc -lavcopts vcodec=mpeg4

end
