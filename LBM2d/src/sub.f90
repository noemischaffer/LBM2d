! Contains routines needed by most other modules
!***************************************************************
module Sub
  implicit none
  public
  double precision :: tini = 5.0d0*TINY(1.0d0) 
contains
!***************************************************************
subroutine normalize2d(vec2d, norm)
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
subroutine dotmv(A,B,C)
  double precision, dimension(2,2) :: A
  double precision, dimension(2) :: B
  double precision, dimension(2) :: C
    C(1) = A(1,1)*B(1)+A(1,2)*B(2)
    C(2) = A(2,1)*B(1)+A(2,2)*B(2)
endsubroutine dotmv
!***************************************************************
endmodule Sub
