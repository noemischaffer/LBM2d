! Contains routines needed by most other modules
!***************************************************************
module Sub
  implicit none
  public
  double precision :: tini = 5.0d0*TINY(1.0d0) 
contains
!***************************************************************
subroutine normalize2d(vec2d)
  double precision,dimension(2) :: vec2d
  double precision :: norm
  norm=sqrt(dot2d(vec2d,vec2d))
  if (norm .gt. tini) vec2d=vec2d/norm
endsubroutine normalize2d
!***************************************************************
real function dot2d(A,B)
  double precision, dimension(2) :: A,B
    dot2d = A(1)*B(1)+A(2)*B(2)
endfunction dot2d
!***************************************************************
endmodule Sub
