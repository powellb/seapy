      subroutine bfort(x,y,z,nin,nout)
      implicit none
      integer, intent(in) :: nin, nout
      !f2py integer, intent(in) :: nin, nout
      real(8), intent(in), dimension(nin) :: x
      !f2py  real, intent(in), dimension(nin) :: x
      real(8), intent(inout), dimension(nout) :: y
      !f2py  real, intent(inout), dimension(nout) :: y
      real(8), intent(out) :: z
      !f2py real, intent(out) :: z

      integer :: i
      
      z=0
      DO i=1,nout
        y(i) = x(i) + y(i)
        z = z + y(i)
      END DO

      end
      
      