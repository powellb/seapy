!=======================================================================
!  Copyright (c) 2002-2020 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This module contains several all purpuse generic routines:          !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    hindices      Finds model grid cell for any datum.                !
!    try_range     Binary search of model grid cell for any datum.     !
!    inside        Closed polygon datum search.                        !
!                                                                      !
!=======================================================================
!
      SUBROUTINE hindices (angler, Xgrd, Ygrd,                          &
     &                     Xpos, Ypos, Ipos, Jpos,                      &
     &                     Lm, Mm, Npos)
!
!=======================================================================
!                                                                      !
!  Given any geographical locations Xpos and Ypos, this routine finds  !
!  the corresponding array cell indices (Ipos, Jpos) of gridded  data  !
!  Xgrd and Ygrd containing each requested location. This indices are  !
!  usually used for interpolation.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     angler      Gridded data angle between X-axis and true EAST      !
!                   (radians).                                         !
!     Xgrd        Gridded data X-locations (usually, longitude).       !
!     Ygrd        Gridded data Y-locations (usually, latitude).        !
!     Xpos        Requested X-locations to process (usually longitude).!
!     Ypos        Requested Y-locations to process (usually latitude). !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Ipos       Fractional I-cell index containing locations in data. !
!     Jpos       Fractional J-cell index containing locations in data. !
!                                                                      !
!  Calls:    try_Range                                                 !
!                                                                      !
!=======================================================================
!
!
!  Imported variable declarations.
!
      implicit none

      integer, intent(in) :: Lm, Mm, Npos
      !f2py integer, intent(in) :: Lm,Mm,Npos

      real(8), intent(in), dimension(Mm,Lm) :: angler, Xgrd, Ygrd
      !f2py real, intent(in), dimension(Mm,Lm) :: angler,Xgrd,Ygrd

      real(8), intent(in), dimension(Npos) :: Xpos, Ypos
      !f2py real, intent(in), dimension(Npos) :: Xpos, Ypos

      real(8), intent(out), dimension(Npos) :: Ipos, Jpos
      !f2py real, intent(out), dimension(Npos) :: Ipos, Jpos
!
!  Local variable declarations.
!
      logical :: found, foundi, foundj

      integer :: Imax, Imin, Jmax, Jmin, i, i0, j, j0, mp, np

      real(8) :: aa2, ang, bb2, diag2, dx, dy, phi
      real(8) :: xfac, xpp, yfac, ypp
      real(8), parameter :: IJspv = -999.0
      real(8), parameter :: Eradius = 6371315.0
      real(8), parameter :: pi = 3.14159265358979323846
      real(8), parameter :: deg2rad = pi / 180.0

      integer :: Is, Ie, Js, Je

      logical :: try_range

      Is = 1
      Ie = Mm
      Js = 1
      Je = Lm
!
!-----------------------------------------------------------------------
!  Determine grid cell indices containing requested position points.
!  Then, interpolate to fractional cell position.
!-----------------------------------------------------------------------
!
      DO np=1, Npos
         Ipos(np)=IJspv
         Jpos(np)=IJspv
!
!     Check each position to find if it falls inside the whole domain.
!     Once it is stablished that it inside, find the exact cell to which
!     it belongs by successively dividing the domain by a half (binary
!     search).
!
         found=try_range(Xgrd, Ygrd, Lm, Mm,                            &
     &        Is, Ie, Js, Je,                                           &
     &        Xpos(np), Ypos(np))
         IF (found) THEN
            Imin=Is
            Imax=Ie
            Jmin=Js
            Jmax=Je
            DO while (((Imax-Imin).gt.1).or.((Jmax-Jmin).gt.1))
               IF ((Imax-Imin).gt.1) THEN
                  i0=(Imin+Imax)/2
                  found=try_range(Xgrd, Ygrd, Lm, Mm,                   &
     &                 Imin, i0, Jmin, Jmax,                            &
     &                 Xpos(np), Ypos(np))
                  IF (found) THEN
                     Imax=i0
                  ELSE
                     Imin=i0
                  END IF
               END IF
               IF ((Jmax-Jmin).gt.1) THEN
                  j0=(Jmin+Jmax)/2
                  found=try_range(Xgrd, Ygrd, Lm, Mm,                   &
     &                 Imin, Imax, Jmin, j0,                            &
     &                 Xpos(np), Ypos(np))
                  IF (found) THEN
                     Jmax=j0
                  ELSE
                     Jmin=j0
                  END IF
               END IF
            END DO
            found=(Is.le.Imin).and.(Imin.le.Ie).and.                    &
     &           (Is.le.Imax).and.(Imax.le.Ie).and.                     &
     &           (Js.le.Jmin).and.(Jmin.le.Je).and.                     &
     &           (Js.le.Jmax).and.(Jmax.le.Je)
!
!     Knowing the correct cell, calculate the exact indices, accounting
!     for a possibly rotated grid.  If spherical, convert all positions
!     to meters first.
!
            IF (found) THEN
               yfac=Eradius*deg2rad
               xfac=yfac*COS(Ypos(np)*deg2rad)
               xpp=(Xpos(np)-Xgrd(Imin,Jmin))*xfac
               ypp=(Ypos(np)-Ygrd(Imin,Jmin))*yfac
!
!     Use Law of Cosines to get cell parallelogram "shear" angle.
!
               diag2=((Xgrd(Imin+1,Jmin)-Xgrd(Imin,Jmin+1))*xfac)**2+      &
     &              ((Ygrd(Imin+1,Jmin)-Ygrd(Imin,Jmin+1))*yfac)**2
               aa2=((Xgrd(Imin,Jmin)-Xgrd(Imin+1,Jmin))*xfac)**2+          &
     &              ((Ygrd(Imin,Jmin)-Ygrd(Imin+1,Jmin))*yfac)**2
               bb2=((Xgrd(Imin,Jmin)-Xgrd(Imin,Jmin+1))*xfac)**2+          &
     &              ((Ygrd(Imin,Jmin)-Ygrd(Imin,Jmin+1))*yfac)**2
               phi=ASIN((diag2-aa2-bb2)/(2.0*SQRT(aa2*bb2)))
!
!     Transform float position into curvilinear coordinates. Assume the
!     cell is rectanglar, for now.
!
               ang=angler(Imin,Jmin)
               dx=xpp*COS(ang)+ypp*SIN(ang)
               dy=ypp*COS(ang)-xpp*SIN(ang)
!
!     Correct for parallelogram.
!
               dx=dx+dy*TAN(phi)
               dy=dy/COS(phi)
!
!     Scale with cell side lengths to translate into cell indices.
!
               dx=MIN(MAX(0.0,dx/SQRT(aa2)),1.0)
               dy=MIN(MAX(0.0,dy/SQRT(bb2)),1.0)
               Ipos(np)=REAL(Imin-1,8)+dx
               Jpos(np)=REAL(Jmin-1,8)+dy
            END IF
         END IF
      END DO

      RETURN
      END SUBROUTINE hindices

      FUNCTION try_range (Xgrd, Ygrd,  Lm, Mm,                          &
     &                            Imin, Imax, Jmin, Jmax, Xo, Yo)
!
!=======================================================================
!                                                                      !
!  Given a grided domain with matrix coordinates Xgrd and Ygrd, this   !
!  function finds if the point (Xo,Yo)  is inside the box defined by   !
!  the requested corners (Imin,Jmin) and (Imax,Jmax). It will return   !
!  logical switch  try_range=.TRUE.  if (Xo,Yo) is inside, otherwise   !
!  it will return false.                                               !
!                                                                      !
!  Calls:   inside                                                     !
!                                                                      !
!=======================================================================
!
!
!  Imported variable declarations.
!
      integer, intent(in) :: Imin, Imax, Jmin, Jmax

      real(8), intent(in), dimension(Mm, Lm) :: Xgrd
      real(8), intent(in), dimension(Mm, Lm) :: Ygrd

      real(8), intent(in) :: Xo, Yo

!
!  Local variable declarations.
!
      integer ::  Nb, i, j, shft, ic

      real(8), dimension(2*(Jmax-Jmin+Imax-Imin)+1) :: Xb, Yb
      logical :: try_range, inside

!
!-----------------------------------------------------------------------
!  Define closed polygon.
!-----------------------------------------------------------------------
!
!  Note that the last point (Xb(Nb),Yb(Nb)) does not repeat first
!  point (Xb(1),Yb(1)).  Instead, in function inside, it is implied
!  that the closing segment is (Xb(Nb),Yb(Nb))-->(Xb(1),Yb(1)). In
!  fact, function inside sets Xb(Nb+1)=Xb(1) and Yb(Nb+1)=Yb(1).
!
      Nb=2*(Jmax-Jmin+Imax-Imin)
      shft=1-Imin
      DO i=Imin,Imax-1
        Xb(i+shft)=Xgrd(i,Jmin)
        Yb(i+shft)=Ygrd(i,Jmin)
      END DO
      shft=1-Jmin+Imax-Imin
      DO j=Jmin,Jmax-1
        Xb(j+shft)=Xgrd(Imax,j)
        Yb(j+shft)=Ygrd(Imax,j)
      END DO
      shft=1+Jmax-Jmin+2*Imax-Imin
      DO i=Imax,Imin+1,-1
        Xb(shft-i)=Xgrd(i,Jmax)
        Yb(shft-i)=Ygrd(i,Jmax)
      END DO
      shft=1+2*Jmax-Jmin+2*(Imax-Imin)
      DO j=Jmax,Jmin+1,-1
        Xb(shft-j)=Xgrd(Imin,j)
        Yb(shft-j)=Ygrd(Imin,j)
      END DO
!
!-----------------------------------------------------------------------
!  Check if point (Xo,Yo) is inside of the defined polygon.
!-----------------------------------------------------------------------
!
      try_range=inside(Nb, Xb, Yb, Xo, Yo)
      RETURN
      END FUNCTION try_range

      FUNCTION inside (Nb, Xb, Yb, Xo, Yo)
!
!=======================================================================
!                                                                      !
!  Given the vectors Xb and Yb of size Nb, defining the coordinates    !
!  of a closed polygon,  this function find if the point (Xo,Yo) is    !
!  inside the polygon.  If the point  (Xo,Yo)  falls exactly on the    !
!  boundary of the polygon, it still considered inside.                !
!                                                                      !
!  This algorithm does not rely on the setting of  Xb(Nb)=Xb(1) and    !
!  Yb(Nb)=Yb(1).  Instead, it assumes that the last closing segment    !
!  is (Xb(Nb),Yb(Nb)) --> (Xb(1),Yb(1)).                               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Reid, C., 1969: A long way from Euclid. Oceanography EMR,         !
!      page 174.                                                       !
!                                                                      !
!  Algorithm:                                                          !
!                                                                      !
!  The decision whether the point is  inside or outside the polygon    !
!  is done by counting the number of crossings from the ray (Xo,Yo)    !
!  to (Xo,-infinity), hereafter called meridian, by the boundary of    !
!  the polygon.  In this counting procedure,  a crossing is counted    !
!  as +2 if the crossing happens from "left to right" or -2 if from    !
!  "right to left". If the counting adds up to zero, then the point    !
!  is outside.  Otherwise,  it is either inside or on the boundary.    !
!                                                                      !
!  This routine is a modified version of the Reid (1969) algorithm,    !
!  where all crossings were counted as positive and the decision is    !
!  made  based on  whether the  number of crossings is even or odd.    !
!  This new algorithm may produce different results  in cases where    !
!  Xo accidentally coinsides with one of the (Xb(k),k=1:Nb) points.    !
!  In this case, the crossing is counted here as +1 or -1 depending    !
!  of the sign of (Xb(k+1)-Xb(k)).  Crossings  are  not  counted if    !
!  Xo=Xb(k)=Xb(k+1).  Therefore, if Xo=Xb(k0) and Yo>Yb(k0), and if    !
!  Xb(k0-1) < Xb(k0) < Xb(k0+1),  the crossing is counted twice but    !
!  with weight +1 (for segments with k=k0-1 and k=k0). Similarly if    !
!  Xb(k0-1) > Xb(k0) > Xb(k0+1), the crossing is counted twice with    !
!  weight -1 each time.  If,  on the other hand,  the meridian only    !
!  touches the boundary, that is, for example, Xb(k0-1) < Xb(k0)=Xo    !
!  and Xb(k0+1) < Xb(k0)=Xo, then the crossing is counted as +1 for    !
!  segment k=k0-1 and -1 for segment k=k0, resulting in no crossing.   !
!                                                                      !
!  Note 1: (Explanation of the logical condition)                      !
!                                                                      !
!  Suppose  that there exist two points  (x1,y1)=(Xb(k),Yb(k))  and    !
!  (x2,y2)=(Xb(k+1),Yb(k+1)),  such that,  either (x1 < Xo < x2) or    !
!  (x1 > Xo > x2).  Therefore, meridian x=Xo intersects the segment    !
!  (x1,y1) -> (x2,x2) and the ordinate of the point of intersection    !
!  is:                                                                 !
!                                                                      !
!                 y1*(x2-Xo) + y2*(Xo-x1)                              !
!             y = -----------------------                              !
!                          x2-x1                                       !
!                                                                      !
!  The mathematical statement that point  (Xo,Yo)  either coinsides    !
!  with the point of intersection or lies to the north (Yo>=y) from    !
!  it is, therefore, equivalent to the statement:                      !
!                                                                      !
!         Yo*(x2-x1) >= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 > 0      !
!  or                                                                  !
!         Yo*(x2-x1) <= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 < 0      !
!                                                                      !
!  which, after noting that  Yo*(x2-x1) = Yo*(x2-Xo + Xo-x1) may be    !
!  rewritten as:                                                       !
!                                                                      !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) >= 0,   if   x2-x1 > 0      !
!  or                                                                  !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) <= 0,   if   x2-x1 < 0      !
!                                                                      !
!  and both versions can be merged into  essentially  the condition    !
!  that (Yo-y1)*(x2-Xo)+(Yo-y2)*(Xo-x1) has the same sign as x2-x1.    !
!  That is, the product of these two must be positive or zero.         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: Nb

      real(8), intent(in) :: Xo, Yo

      real(8), intent(inout) :: Xb(Nb+1), Yb(Nb+1)
!
!  Local variable declarations.
!
      integer, parameter :: Nstep =128

      integer :: crossings, i, inc, k, kk, nc

      integer, dimension(Nstep) :: Sindex

      real(8) :: dx1, dx2, dxy

      logical :: inside
!
!-----------------------------------------------------------------------
!  Find intersections.
!-----------------------------------------------------------------------
!
!  Set crossings counter and close the contour of the polygon.
!
      crossings=0
      Xb(Nb+1)=Xb(1)
      Yb(Nb+1)=Yb(1)
!
!  The search is optimized.  First select the indices of segments
!  where Xb(k) is different from Xb(k+1) and Xo falls between them.
!  Then, further investigate these segments in a separate loop.
!  Doing it in two stages takes less time because the first loop is
!  pipelined.
!
      DO kk=0,Nb-1,Nstep
        nc=0
        DO k=kk+1,MIN(kk+Nstep,Nb)
           IF (((Xb(k+1)-Xo)*(Xo-Xb(k)).ge.0.0).and.                     &
     &          (Xb(k).ne.Xb(k+1))) THEN
            nc=nc+1
            Sindex(nc)=k
          END IF
        END DO
        DO i=1,nc
          k=Sindex(i)
          IF (Xb(k).ne.Xb(k+1)) THEN
            dx1=Xo-Xb(k)
            dx2=Xb(k+1)-Xo
            dxy=dx2*(Yo-Yb(k))-dx1*(Yb(k+1)-Yo)
            inc=0
            IF ((Xb(k).eq.Xo).and.(Yb(k).eq.Yo)) THEN
              crossings=1
              goto 10
            ELSE IF (((dx1.eq.0.0).and.(Yo.ge.Yb(k  ))).or.             &
     &              ((dx2.eq.0.0).and.(Yo.ge.Yb(k+1)))) THEN
              inc=1
            ELSE IF ((dx1*dx2.gt.0.0).and.                              &
     &              ((Xb(k+1)-Xb(k))*dxy.ge.0.0)) THEN  ! see note 1
              inc=2
            END IF
            IF (Xb(k+1).gt.Xb(k)) THEN
              crossings=crossings+inc
            ELSE
              crossings=crossings-inc
            END IF
          END IF
        END DO
      END DO
!
!  Determine if point (Xo,Yo) is inside of closed polygon.
!
  10  IF (crossings.eq.0) THEN
        inside=.FALSE.
      ELSE
        inside=.TRUE.
      END IF
      RETURN
      END FUNCTION inside
