c This routine was written by Bruce Cornuelle at Scripps Institute
c of Oceanography. Because it uses LAPACK routines, it is under
c the BSD-License consistent with LAPACK
      
      
      subroutine oa2d(x,y,datain,xhat,yhat,a,b,pmap,dataout,
     &                   errout,PtIn,PtOut,PtIns,debug)
    	implicit none
      integer, intent(in) :: PtIn,PtOut,PtIns
      !f2py integer, intent(in) :: PtIn,PtOut,PtIns
      real(8), intent(in), dimension(PtIn) :: x, y, datain
      !f2py  real, intent(in), dimension(PtIn) :: x, y, datain
      real(8), intent(in), dimension(PtOut) :: xhat, yhat
      !f2py  real, intent(in), dimension(PtOut) :: xhat, yhat
      real(8), intent(inout), dimension(PtOut,PtIns) :: pmap
      !f2py  real, intent(in,out), dimension(PtOut,PtIns) :: pmap
      real(8), intent(in) :: a, b
      !f2py  real, intent(in) :: a, b
      real(8), intent(out), dimension(PtOut) :: dataout, errout
      !f2py  real, intent(out), dimension(PtOut) :: dataout, errout
      logical, intent(in) :: debug
      
    	integer :: ii, jj,i,j,Isel
    	real(8)    :: distb, dista, dist_ref,dist
      real(8)    :: xs(PtIn),ys(PtIn),datains(PtIn)



c init to zero matricies
    	do i=1, PtOut
    		dataout(i)=0
        errout(i)=0
      enddo

      if (pmap(1,1).eq.0) then
        if (debug .eqv. .true.) then
          print *, 'PtIns = ', PtIns
          print *, '  building map '
        endif
    	 call makeMap(x,y,datain,PtIn,
     &             pmap,PtIns,xhat,yhat,PtOut)

        if (debug .eqv. .true.) then
          print *, '  building map DONE.'
        endif
      endif

c do the objective map
    	do i=1,PtOut
            do ii=1,PtIns
                 Isel           =      pmap( i    ,ii )
    	     xs(ii)=         x( Isel )
    	     ys(ii)=         y( Isel )
    	     datains(ii)=    datain( Isel )
    	  enddo

	    call oa(xs,ys,datains,xhat(i),yhat(i),PtIns,1,
     &                         dataout(i),errout(i),a,b)
      enddo
	return
	end

c========================================================
c	select closest points
c========================================================

      subroutine makeMap(x,y,datain,PtIn,
     &                    pmap,PtIns,xhat,yhat,PtOut)
    	integer PtIn,PtOut, ii, jj,i,j,PtIns,jout
    	real*8 x(PtIn),y(PtIn),datain(PtIn)
    	real*8 xhat(PtOut),yhat(PtOut)
    	real*8 pmap(PtOut,PtIns)
    	real*8 distb, dista, dist_ref,dist(PtIns)

      real*8 xs(PtIn),ys(PtIn),datains(PtIn)

c     make sure the number of points to select is
c     not bigger then input points
	if (PtIn.lt.PtIns) then
	   print*, 'PtIn lower than PtIns '
	   return
	endif

	do jout=1,PtOut
c     now find PtIns closest points
      do ii=1,PtIns
        dist(ii)=30000
      enddo

	do i=1,PtIn
	  if (datain(i).eq.-999999.0) goto 68
          dista=abs(x(i)-xhat(jout))
          distb=abs(y(i)-yhat(jout))
          dist_ref=dista + distb

          do ii=1,PtIns
             if (dist_ref .le. dist(ii)) then
              do jj=PtIns,ii+1,-1
                 dist(jj)=dist(jj-1)
                 xs(jj)=xs(jj-1)
                 ys(jj)=ys(jj-1)
		     pmap(jout,jj) = pmap(jout,jj-1)
              enddo
                dist(ii)=dist_ref
                xs(ii)=x(i)
		    ys(ii)=y(i)
		    pmap(jout,ii) = i
              goto 68
             endif
          enddo


68      continue
        enddo
        enddo



c	do i=1,PtIns
c	  xs(i)=x(i)
c	  ys(i)=y(i)
c	  datains(i)=datain(i)
c	enddo

	return
	end

      subroutine oa3d(x,y,z,datain,xhat,yhat,zhat,a,b,pmap,dataout,
     &                   errout,PtIn,PtOut,PtZin,PtZout,PtIns,debug)
      implicit none
      integer, intent(in) :: PtIn,PtOut,PtIns,PtZin,PtZout
      !f2py integer, intent(in) :: PtIn,PtOut,PtIns,PtZin,PtZout
      real(8), intent(in), dimension(PtIn) :: x, y
      !f2py  real, intent(in), dimension(PtIn) :: x, y
      real(8), intent(in), dimension(PtIn, PtZin) :: z, datain
      !f2py real(8), intent(in), dimension(PtIn, PtZin) :: z, datain
      real(8), intent(in), dimension(PtOut) :: xhat, yhat
      !f2py  real, intent(in), dimension(PtOut) :: xhat, yhat
      real(8), intent(in), dimension(PtOut,PtZout) :: zhat
      !f2py real(8), intent(in), dimension(PtOut,PtZout) :: zhat
      real(8), intent(inout), dimension(PtOut,PtIns) :: pmap
      !f2py  real, intent(in,out), dimension(PtOut,PtIns) :: pmap
      real(8), intent(in) :: a, b
      !f2py  real, intent(in) :: a, b
      real(8), intent(out), dimension(PtOut,PtZout) :: dataout
      !f2py real(8), intent(out), dimension(PtOut,PtZout) :: dataout
      real(8), intent(out), dimension(PtOut) :: errout
      !f2py  real, intent(out), dimension(PtOut) :: errout
      logical, intent(in) :: debug

      integer :: ii, jj,i,j,Isel,iz
      real(8), dimension(1,PtZout) :: dataout_tmp
      real(8)    :: zrOut(PtOut,PtZout),zrIn(PtIn,PtZin)
      real(8)    :: xs(PtIns),ys(PtIns),zs(PtZin),out_y(PtZout)
      real(8)    :: datainsZ(PtZin),datains(PtIns,PtZout)
      real(8)    :: zhat_tmp(PtZout)

c added by bdc:
      integer    :: nextrap,izfirst(PtIns),iii,ngood,izfirstmax
      real(8)    :: dataouts(PtIns)
      real(8)    :: dstart,dcrit
      real(8)    :: distb, dista, dist_ref,dist(PtIns)
      integer    :: ilook,pmap2(PtIns),iiz
      integer    :: idobdc,nptmin,idbbdc
      real(8)    :: epsz
      real(8)    :: zsbot(PtIns)

c init to zero matrices
      do i=1, PtOut
        do iz=1,PtZout
          dataout(i,iz)=0
        enddo
        errout(i)=0
      enddo

c bdc: hardwired parameters for enhanced mapping
c idobdc : 0 = do it old way: just use the vertical extrapolation
c          so points always exist.
c    1 = new way: don't allow vertical extrapolation beyond epsz,
c        and re-map using only the good points.  (Takes a long time!)
      idobdc = 1
c      print *, 'idobdc = ',idobdc
c dcrit = critical value of abs sum distance:
c   for choosing nearest points when looking for more points.
c   values bigger than this are too far away, and those points aren't used
      dcrit = 3
c      print *, 'max absolute sum distance = ',dcrit
c epsz = tolerance allowed in depth extrapolation.  If a point
c    is less than epsz*local depth too shallow, it still can be used
      epsz = 0.05
c      print *, 'fractional depth tolerance = ',epsz
c nptmin = minimum number of points to have to use the objective map
      nptmin = 5
c      print *, 'minimum number of points required for a map= ',nptmin
c if nptmin = 2, and there is one point, it will be used for the output
c idbbdc : debug parameter for output.  0 = none, 1 = more, 2 = more, etc.
      idbbdc = 0

      if (idobdc .eq. 1 .and. idbbdc .gt. 0) then
        open(10,file='bdc.txt')
      endif

      if (pmap(1,1).eq.0) then
        print *, 'PtIns = ', PtIns
        print *, '  need to build pmap for subroutine oa_ex '
        return
      endif
c  pmap(1:PtOut,1:PtIns) tells which PtIns selected input points go
c    with each output point.  They were chosen outside on the basis
c    of nearest neighbors.

c      print*, PtIn,PtOut,PtZin,PtZout, PtIns
c do the objective map
c loop over output points (to be interpolated to)
      do i=1,PtOut

c the original program does vertical interpolation first, before horizontal,
c meaning that in areas of strong topography, surface values
c are extrapolated downward, which is scary.
c This should be modified, so that only ANOMALIES with respect to
c horizontal averages are extrapolated.

c     vertical linear interpolation
c z points at i'th output point
        do iz=1,PtZout
          zhat_tmp (iz) = zhat(i, iz)
        enddo

c count number of points where extrapolation happened.
        nextrap = 0
c loop over the PtIns nearest (in x,y) input points
c (ii is input point index, i is output index)
        do ii=1,PtIns
          Isel           =      pmap( i    ,ii )
          xs      (ii   )=         x( Isel     )
          ys      (ii   )=         y( Isel     )
c loop over the layers of the input points
          do iz=1,PtZin
            datainsZ(iz)=    datain( Isel, iz )
            zs      (iz)=         z( isel ,iz )
          enddo
c so zs(1) is the deepest depth at input point isel
c save it for diagnostics
          zsbot(ii) = zs(1)

c interpolate from zs, datainsZ, to out_y at depths zhat_tmp
c CAUTION!  This extrapolates linearly to missing depths!
          call lintrp (PtZin,zs,datainsZ,PtZout,zhat_tmp,out_y)
c now out_y is at the OUTPUT depths, although still at the input
c  point locations.
         ! assignment loop: put the 1-d array into a 2-d array
          do iz=1,PtZout
            datains(ii, iz ) = out_y(iz)
          enddo

          if (idobdc .eq. 1) then
c check for vertical extrapolation

c note: depths are negative, and point 1 is at the bottom
c          if (-zhat_tmp(PtZout) .gt. zs (PtZin)) then
            if (-zhat_tmp(1)*(1-epsz) .gt. -zs (1)) then
c the output max depth is larger than the input max depth at this input point
c some of the input values were extrapolated to deeper output depths
c add to count of input points that were extrapolated
              nextrap = nextrap + 1
c find the first output depth index that was NOT extrapolated
c (deepest is 1)
c            iz = PtZout
              iz = 1
              do while (-zhat_tmp(iz)*(1-epsz) .gt. -zs(1))
                iz = iz + 1
              enddo
c now iz holds the first (deepest) good output depth for input station ii.
              izfirst(nextrap) = iz
            endif
c end of:          if (idobdc .eq. 1) then
          endif
c end of:        do ii=1,PtIns
c end of loop over input points for this output point
        enddo
c now datains(1:PtIns,1:Ptzout) holds the input data interpolated to the
c  output depths.

c	  print*, ' OA - Estimated COV ', i
c now do the objective mapping to the output point locations
        call oa3(xs,ys,datains,xhat(i),yhat(i),PtIns,1,PtZout,
     c    dataout_tmp,errout(i),a,b)

        ! assignment loop
        do iz = 1,PtZout
          dataout(i,iz) = dataout_tmp(1,iz)
        enddo

c now dataout holds the values made using vertical extrapolation

c bdc: new step added: fix up the points that were vertically extrapolated
c  using a plane fit
c        if (nextrap .gt. 0 .and. idobdc .eq. 1) then
c nextrap is only > 0 IF idobdc=1
        if (nextrap .gt. 0) then
c some points were extrapolated
          if (idbbdc .gt. 0) then
!           write(10,*) 'Output pt = ',i, ' nextrap = ',nextrap
!           write(10,*) 'xhat(i)= ',xhat(i),' yhat(i)= ',yhat(i)
!           write(10,*) 'bottom depth =  ',-zhat_tmp(1)
!           write(10,*) 'izfirst= ',(izfirst(iii),iii=1,nextrap)
          if (idbbdc .gt. 1) then
          if (i .eq. 64) then
!             write(10,*) 'special output for point 64'
c depths
!             write(10,1002)(zhat_tmp(iz),dataout(i,iz),iz=1,PtZout)
 1002       format(2g13.5)
c now put out all the associated points
c lon, lat, depth
!             write(10,1001)(xs(ii),ys(ii),zsbot(ii),ii=1,PtIns)
 1001         format(2f13.7,g13.5)
c values of the stations
            do iz=1,PtZout
!               write(10,1003)(datains(ii,iz),ii=1,PtIns)
            enddo
 1003       format(30g13.5)
          endif
          endif
          endif
c          if (nextrap .ge. PtIns) then
c            print*, 'Stopping: NO pts are good'
c            stop 'no good points in oa'
c          endif
c          if (nextrap .gt. PtIns/2) then
c            print*, 'Warning: fewer than 1/2 selected pts are good'
c          endif

c find the shallowest first good depth
          izfirstmax = izfirst(1)
          if (nextrap .gt. 1) then
            do iii = 2,nextrap
              izfirstmax = max(izfirstmax,izfirst(iii))
            enddo
          endif
          if (idbbdc .gt. 0) then
!           write(10,*) 'deepest all good depth level = ',izfirstmax
!           write(10,*) ' =  ',-zhat_tmp(izfirstmax)
          endif
c-----------------------------
c this was already done above
c-----------------------------
c do oa over all depths, just like before:
c (it's too hard to just do it from the first good depth and up)
c          call oa(xs,ys,datains,xhat(i),yhat(i),PtIns,1,PtZout,
c     c      dataout_tmp,errout(i),a,b)
c
cc put in the output array
c          do iz = 1,PtZout
c            dataout(i,iz) = dataout_tmp(1,iz)
c	  enddo
c-----------------------------

c now loop over output depths where some points were extrapolated
c work down from shallow
          do iz = izfirstmax-1,1,-1
c find all good points at this level
            ngood = 0
            do ii=1,PtIns
              Isel           =      pmap( i    ,ii )

c              write(10,*) 'iz= ',iz,', testing input point ',Isel
c              write(10,*) 'x(isel)= ',x(isel),', y(isel)= ',y(isel)
c              write(10,*) 'max depth here is ', -z( isel, 1)
c              write(10,*) 'desired depth is ',-zhat_tmp(iz)

              if (-zhat_tmp(iz)*(1-epsz) .le. -z( isel, 1)) then
c this output depth is less than or equal to the input max depth:
c it's a good station
                ngood = ngood+1
                xs(ngood)=         x( Isel     )
                ys(ngood)=         y( Isel     )
c data values (already have been interpolated to this depth)
                dataouts(ngood) = datains(ii, iz )
c                write(10,*) 'value saved= ',dataouts(ii)
c              else
c                write(10,*) 'iz= ',iz,', rejected input point ',Isel
c                write(10,*) 'x(isel)= ',x(isel),', y(isel)= ',y(isel)
c                write(10,*) 'max depth here is ', -z( isel, 1)
c                write(10,*) 'desired depth is ',-zhat_tmp(iz)
              endif
            enddo
c now have a list of the ngood good stations.
c use the horizontal objective map (with a least-squares plane fit)
c to extrapolate to the bad stations

            if (idbbdc .gt. 0) then
!             write(10,*) ngood,' good points found'
!             write(10,*) iz,'=iz; desired depth is ',-zhat_tmp(iz)
c            if (izfirstmax .gt. 6 .and. iz .lt. 3 .and.
c     c          -zhat_tmp(iz) .gt. 500) then
c              write(10,*) 'map input values, lat, lon: '
c              write(10,1001)(xs(ii),ys(ii),dataouts(ii),ii=1,ngood)
c            endif
            endif
            if (ngood .ge. nptmin) then
c have enough good points to just do the mapping
              call oa3(xs,ys,dataouts,xhat(i),yhat(i),ngood,1,1,
     c          dataout_tmp,errout(i),a,b)

              if (idbbdc .gt. 0) then
!               write(10,*) ngood,' good points; doing objective map'
c              if (izfirstmax .gt. 5 .and. iz .lt. 3) then
!                 write(10,*) 'original value = ',dataout(i,iz)
!                 write(10,*) 'new value = ',dataout_tmp(1,1)
c                if (i .ge. 200 .and. iz .eq. 1) then
c                  write(10,*) 'at i = 200; stop the program'
c                  close(10)
c                  stop 'testing subroutine oa_ex'
c                endif
c                if (-zhat_tmp(iz) .gt. 500 .and.
c     c               abs(dataout(i,iz)) .gt.  0.8) then
c                  write(10,*) 'extrapolation gap > 5; stop the program'
c                  stop 'testing subroutine oa_ex'
c                endif
              endif
c save the output for this depth
              dataout(i,iz) = dataout_tmp(1,1)

            else
c not enough points for an objective map

c decide whether or not to increase the search radius
c              ilook = 0
              ilook = 1
              if (ilook .eq. 1) then
c                print *, 'looking for more points'
!                 write(10,*) 'looking for more points'
c critical value: values bigger than this are too far away
                dcrit = 3
c                print *, 'max absolute sum distance = ',dcrit
!                 write(10,*) 'max absolute sum distance = ',dcrit

c have to look for more nearby points.
c (take code from makemap.f)

c     now find PtIns closest points with good data
c start count of good points (within max distance)
                ngood = 0
c starting distance: should be huge
                dstart = 300000
                do ii=1,PtIns
                  dist(ii)=dstart
                enddo
                do ii=1,PtIn
c check for good point
c check for deep enough
c check to see if this input point is as deep as the desired level
                  if (-zhat_tmp(iz)*(1-epsz) .gt. -z( ii, 1)) goto 68
c take deepest point for datain
                  if (datain(ii,1).eq.-999999.0) goto 68
c stay here if point is deep enough and no bad flag
                  dista=abs(x(ii)-xhat(i))
                  distb=abs(y(ii)-yhat(i))
                  dist_ref=dista + distb
                  if (dist_ref .le. dcrit .and.
     c                dist_ref .lt. dist(PtIns)) then
c this is a good point, and can get into the list
                    ngood = ngood + 1

c find the place for this distance in the sorted list of distances
c (smallest = 1st)
c if two points tie, the first one gets a lower number
                    do iii=1,PtIns
                      if (dist_ref .lt. dist(iii)) then
c adjust the points after it in the list
c go from the end to iii+1
                        do jj=PtIns,iii+1,-1
                          dist(jj)=dist(jj-1)
c                         xs(jj)=xs(jj-1)
c                         ys(jj)=ys(jj-1)
                          pmap2(jj) = pmap2(jj-1)
                        enddo
c add the point at iii
                        dist(iii)=dist_ref
c                        xs(iii)=x(ii)
c                        ys(iii)=y(ii)
                        pmap2(iii) = ii
                        goto 68
                      endif
                    enddo
                  endif
 68               continue
c end of:    do ii=1,PtIn
                enddo

c limit the number of good points to the max allowed
                ngood = min(ngood,PtIns)
!                 write(10,*) ngood,' good points found'
c                if (dist(PtIns) .lt. dstart) then
cc now have PtIns good points nearby.
                if (ngood .gt. 0) then
!                   write(10,*) 'distances: ',(dist(ii),ii=1,ngood)

c do the vertical interpolation
c loop over the ngood nearest (in x,y) input points
c (ii is input point index, i is output index)
                  do ii=1,ngood
                    Isel           =      pmap2(ii )
                    xs      (ii   )=         x( Isel     )
                    ys      (ii   )=         y( Isel     )

c load the array to be interpolated from
c loop over the layers of the input points
c (caution: iz is the level in the OUTPUT depths)
                    do iiz=1,PtZin
                      datainsZ(iiz)=    datain( Isel, iiz )
                      zs      (iiz)=         z( isel ,iiz )
                    enddo
c only need the value at one depth: zhat_tmp(iz)
                    call lintrp (PtZin,zs,datainsZ,1,zhat_tmp(iz),out_y)
                    dataouts(ii) = out_y(1)
c end of loop over good points; now array dataouts() is loaded
                  enddo
                  if (ngood .ge. nptmin) then
c do an objective map on the good points
!                   write(10,*) 'original value = ',dataout(i,iz)

                    call oa3(xs,ys,dataouts,xhat(i),yhat(i),ngood,1,1,
     c                dataout_tmp,errout(i),a,b)
!                     write(10,*) 'mapped value = ',dataout_tmp(1,1)
c save the output for this depth
                    dataout(i,iz) = dataout_tmp(1,1)
                  else
!                     write(10,*) 'original value = ',dataout(i,iz)
                    if (nptmin .eq. 2 .and. ngood .eq. 1) then
!                       write(10,*) 'only one good point, use that value'
!                       write(10,*) 'new value = ',dataouts(1)
                      dataout(i,iz) = dataouts(1)
                    else
!                       write(10,*) 'too few good points; take pt above'
                      dataout(i,iz) = dataout(i,iz+1)
!                       write(10,*) 'new value = ',dataout(i,iz+1)
                    endif
                  endif
                else
!                   write(10,*) 'no good points; take pt above'
                  dataout(i,iz) = dataout(i,iz+1)
                endif
              else
!                 write(10,*) 'no new search'
!                 write(10,*) 'original value = ',dataout(i,iz)
                if (nptmin .eq. 2 .and. ngood .eq. 1) then
!                   write(10,*) 'only one good point, use that value'
!                   write(10,*) 'new value = ',dataouts(1)
                  dataout(i,iz) = dataouts(1)
                else
!                   write(10,*) 'no good points; just take pt above'
!                   write(10,*) 'new value = ',dataout(i,iz+1)
                  dataout(i,iz) = dataout(i,iz+1)
                endif
c               stop 'testing subroutine oa_ex'
c end of:  if (ilook .eq. 1) then
              endif
c end of:  if (ngood .gt. nptmin) then
            endif
c end of:   do iz = izfirstmax-1,1,-1
          enddo

c end of:     if (nextrap .gt. 0 .and. idobdc .eq. 1) then
        endif

      enddo
      return
      end


      subroutine lintrp (n,x,y,ni,xi,yi)
c
c=======================================================================
c  Copyright (c) 1996 Rutgers University                             ===
c=======================================================================
c                                                                    ===
c  Given arrays X and Y of length N, which tabulate a function,      ===
c  Y = F(X),  with the Xs  in ascending order, and given array       ===
c  XI of lenght NI, this routine returns a linear interpolated       ===
c  array YI.                                                         ===
c                                                                    ===
c=======================================================================
c
c-----------------------------------------------------------------------
c  Define local variable.
c-----------------------------------------------------------------------
c
      implicit none
      integer i, ii, j, n, ni
      real*8 d1, d2, dydx1, dydx2
      real*8 x(n), y(n), xi(ni), yi(ni)
c
c-----------------------------------------------------------------------
c  Begin executable code.
c-----------------------------------------------------------------------
c
c slope at start
      dydx1=(y(2)-y(1))/ abs(x(2) - x(1))
c slope at end
	dydx2=(y(n)-y(n-1))/ abs(x(n) - x(n-1))

c loop over output points
      do 30 j=1,ni
        if (xi(j).le.x(1)) then
c output abscissa is before x(1); have to extrapolate!
          ii=1
c use the starting slope computed above.
c bdc: commented out extrapolation
c	  print*, ' lintrp: extrapolating above first point!'
c          yi(j)=y(1) - dydx1 * abs(xi(j) - x(1))
          yi(j)=y(1)

        elseif (xi(j).ge.x(n))then
c output abscissa is after x(n); have to extrapolate!!
c extrapolate with a constant
          yi(j)=y(n)
c bdc: commented out extrapolation
c	  print*, ' lintrp: extrapolating below last point!'
c	    yi(j)=y(n) + dydx2 * abs(xi(j) - x(n))
        else
c are within the range; find bracketing points
          do 10 i=1,n-1
            if ((x(i).lt.xi(j)).and.(xi(j).le.x(i+1))) then
c now xi(j) > x(i) and xi(j) <= x(i+1)  (so xi(j) is between i and i+1)
              ii=i
              goto 20
            endif
 10       continue
c compute linear weights
 20       d1=xi(j)-x(ii)
          d2=x(ii+1)-xi(j)
          yi(j)=(d1*y(ii+1)+d2*y(ii))/(d1+d2)
        endif
 30   continue
      return
      end

c========================================================
c	Objective analysis
c  added variable k for exponent on cov & noise eqns
c  set k=2 for gaussian, 1 for markov expl  3/28/06
c
c========================================================
      subroutine oa3(x,y,datain,xhat,yhat,PtIn,PtOut,PtZout,
     c                                      dataout,errout,a,b)
	implicit none
	integer PtIn,PtOut, ii, jj,i,j,iz,PtZout,k
	real*8 x(PtIn),y(PtIn),xhat(PtOut),yhat(PtOut),datap(PtIn),
     C     dataout(PtOut,PtZout),errout(PtOut),datain(PtIn,PtZout)
c        real*8 issing,agd,datain_tmp(PtIn),dataout_tmp(PtOut)
        real*8 agd,datain_tmp(PtIn),dataout_tmp(PtOut)
        integer issing
c OA variables
        real*8 a,b
	real*8 GD(PtOut,PtIn), DD(PtIn,PtIn), GG(PtOut),
     c      DDinv(PtIn,PtIn) ,m(3),GD_DD(PtOut,PtIn),
     c      GDt(PtIn,PtOut)
c extra parameters added for bdc revised RemoveMean subroutine
	real*8 xmean,ymean,xl,yl
      integer ipiv(PtIn),INFO

c	print*, '               PtIn      PtOut'
c	print*, 'Index ',       PtIn,     PtOut

c set parameter
c      a=1.05
c	b=1.15
c        print *,'Using ',a,b

c init to zero matricies
	do i=1, PtOut
	   do iz=1,PtZout
		dataout(i,iz)=0
         enddo
            errout(i)=0
      enddo

c BELOW: changed to k instead of 2 for expt
c set k = 1 (markov) or 2 (gaussian)

c Data - Grid Covariance = GD

      k=1
c      print *,'k = ',k

      do i=1,PtOut
	do j=1,PtIn
	  GD(i,j)=exp(-(  abs((x(j)-xhat(i))/a)**k +
     c                    abs((y(j)-yhat(i))/b)**k   ))
      enddo
	enddo

c Data - Data Covariance = DD
      do i=1,PtIn
	do j=1,PtIn
	  DD(i,j)=exp(-(  abs((x(j)-x(i))/a)**k +
     c                    abs((y(j)-y(i))/b)**k     ))
      enddo
	enddo

c Add noise to Diagonal	of data-data Covariance
	do i=1,PtIn
	  DD(i,i)=DD(i,i) + DD(i,i)*0.1
	enddo

c grid-grid covariance: just used for error
      do i=1,PtOut
        GG(i)=exp(-(  abs((xhat(i)-xhat(i))/a)**k +
     c                  abs((yhat(i)-yhat(i))/b)**k  ))
      enddo

c Compute estimate at grid point
c      print *,'Inverse'
c Manu's routine ..slow!
c      call inverse(DD,PtIn,Ptin,DDinv,issing)
c      call dot(GD,DDinv,PtOut,PtIn,PtIn,GD_DD)
c      call dot(GD_DD,datap,PtOut,PtIn,1,dataout)

c LAPACK routine ...:)
      call ftranspose(GD,PtOut,PtIn,GDt)
c now GDt is the transpose of the grid-data covariance
c solve the problem DD * x = GDt
	call DGESV( PtIn, PtOut, DD, PtIn, ipiv, GDt, PtIn, INFO )
c      print *,'Transpose', INFO
c transpose back to get GD_DD = GD*(DD)^-1
      call ftranspose(GDt,PtIn,PtOut,GD_DD)
c now GD_DD is the operator to map in 2-d

c
c loop on PtZout
       do iz=1,PtZout
       do i=1,PtIn
          datain_tmp(i)=datain(i,iz)
	 enddo

c remove mean (and slope) using a plane fit
c       call removeMean(x,y,datain_tmp,PtIn,m,datap)
c revised version from bdc
       call removeMean3(x,y,datain_tmp,PtIn,m,datap,xmean,ymean,xl,yl)
c datap are the output de-meaned values, m is the parameter vector
c do the estimation:
       call fdot(GD_DD,datap,PtOut,PtIn,1,dataout_tmp)
c  A=GD/DD;
c  dhat=A*d;
c put the mean back in:
c      call addMean3(xhat,yhat,dataout_tmp,PtOut,m)
      call addMean3(xhat,yhat,dataout_tmp,PtOut,m,xmean,ymean,xl,yl)

	  do i=1,PtOut
	     dataout(i,iz) = dataout_tmp(i)
	  enddo
	enddo

c compute errormap
	do i=1,PtOut
	  agd=0
	  do j=1,PtIn
	     agd=agd + GD_DD(i,j) * GD(i,j)
	  enddo
	  errout(i)=GG(i) - agd
	enddo

	return
	end




c========================================================
c	remove mean
c========================================================
      subroutine removeMean3(x,y,datain,Pt,m,datap,xmean,ymean,xl,yl)
	implicit none
	integer Pt, ii, jj,i,j,r,c
	real*8 x(Pt),y(Pt),datain(Pt),m(3),G(Pt,3),datap(Pt),
     c      dmean(Pt),Gt(3,Pt), b(3,1),Gp(3,3),A(3,3)
         integer issing
c extra output pars added by bdc
	real*8 xmean,ymean,xl,yl
c extra working parameters
	real*8 xmin,xmax,ymin,ymax

   ! setup matrix G
      r=Pt
	c=3

c compute forward problem matrix for plane fit
c compute min and max, and mean to non-dimensionalize
        xmin = x(1)
        ymin = y(1)
        xmax = x(1)
        ymax = y(1)
        xmean = 0;
        ymean = 0;
      do i=1,r
        xmean = xmean + x(i)
        ymean = ymean + y(i)
        xmin = min(xmin,x(i))
        ymin = min(ymin,y(i))
        xmax = max(xmax,x(i))
        ymax = max(ymax,y(i))
      enddo
      xmean = xmean/r
      ymean = ymean/r
c these parameters are output so the forward model in addmean is the same
      if (xmax .gt. xmin) then
        xl = xmax-xmin
      else
        xl = 1.0
      endif
      if (ymax .gt. ymin) then
        yl = ymax-ymin
      else
        yl = 1.0
      endif
c compute non-dimensionalized forward problem matrix
      do i=1,r
        G(i,1)=(x(i)-xmean)/xl
        G(i,2)=(y(i)-ymean)/yl
        G(i,3)=1.d0
      enddo

      call ftranspose(G,r,c,Gt)

      call fdot(Gt,datain ,c,r,1, b)
      call fdot(Gt,G,c,r,c,Gp)
c bdc: add error to diagonal of Gp (3x3)
c assume errors are all 1, so add noise/signal ratio
c this n/s level was picked arbitrarily!
c remember that Gp (1,1) is r, and Gp(2,2) and Gp(3,3) are probably about .25r,
c so if r is 30, 0.1 is a small number to add.
      do j=1,3
        Gp(j,j) = Gp(j,j) + 0.1
      enddo

      call inverse(Gp,c,c,A,issing)

      if (issing .eq. 1) then
          print*, 'Remove Mean: error in computing inverse'
      else
          call fdot(A,b,c,c,1,m)
	    call fdot(G,m,r,c,1,dmean)
      endif

	do i=1,Pt
	  datap(i)=datain(i) - dmean(i)
	enddo


	return
	end

c========================================================
c	add mean
c========================================================
      subroutine addMean3(x,y,datain,Pt,m,xmean,ymean,xl,yl)
	implicit none
	integer Pt, ii, jj,i,j,r,c
	real*8 xmean,ymean,xl,yl
	real*8 x(Pt),y(Pt),datain(Pt),m(3),G(Pt,3),
     c      dmean(Pt),Gt(3,Pt), b(3,1),Gp(3,3),A(3,3)
        integer issing
c extra input parameters added by bdc

   ! setup matrix G
      r=Pt
	c=3
c old: dimensional, no de-meaned:
c      do i=1,r
c        G(i,1)=x(i)
c        G(i,2)=y(i)
c        G(i,3)=1.d0
c      enddo
c compute forward problem matrix for plane fit
c compute non-dimensionalized forward problem matrix
      do i=1,r
        G(i,1)=(x(i)-xmean)/xl
        G(i,2)=(y(i)-ymean)/yl
        G(i,3)=1.d0
      enddo

      call fdot(G,m,r,c,1,dmean)
	do i=1,Pt
	  datain(i)=datain(i) + dmean(i)
	enddo

	return
	end

c========================================================
c	Objective analysis
c========================================================
      subroutine oa(x,y,datain,xhat,yhat,PtIn,PtOut,dataout,errout,a,b)
    	implicit none
    	integer PtIn,PtOut, ii, jj,i,j
    	real*8 x(PtIn),y(PtIn),xhat(PtOut),yhat(PtOut),datap(PtIn),
     &                 dataout(PtOut),errout(PtOut),datain(PtIn)
    	real*8 distb, dista, dist_ref,dist
      real*8 issing,agd
c OA variables
      real*8 a,b
    	real*8 GD(PtOut,PtIn), DD(PtIn,PtIn), GG(PtOut),
     &      DDinv(PtIn,PtIn) ,m(3),GD_DD(PtOut,PtIn),
     &      GDt(PtIn,PtOut)
      integer ipiv(PtIn),INFO

c	print*, '               PtIn      PtOut'
c	print*, 'Index ',       PtIn,     PtOut

c set parameter
c      a=1.05
c	b=1.15
c        print *,'Using ',a,b

c init to zero matricies
    	do i=1, PtOut
    		dataout(i)=0
                errout(i)=0
          enddo

c remove mean as a plane
       call removeMean(x,y,datain,PtIn,m,datap)

c Data - Grid Covariance = GD
            do i=1,PtOut
      	do j=1,PtIn
      	  GD(i,j)=exp(-( ((x(j)-xhat(i))**2)/a**2 +
     &                     ((y(j)-yhat(i))**2)/b**2 ))
            enddo
      	enddo

c Data - Data Covariance = DD
            do i=1,PtIn
      	do j=1,PtIn
      	  DD(i,j)=exp(-( ((x(j)-x(i))**2)/a**2 +
     &                    ((y(j)-y(i))**2)/b**2 ))
            enddo
      	enddo

c Add noise to Diagonal	and Grid -Grid Covariance
      	do i=1,PtIn
      	  DD(i,i)=DD(i,i) + DD(i,i)*0.1
c	  DD(i,i)=DD(i,i) + 0.1
      	enddo

      do i=1,PtOut
        GG(i)=exp(-( ((xhat(i)-xhat(i))**2)/a**2 +
     &                     ((yhat(i)-yhat(i))**2)/b**2 ))
	enddo

c Compute estimate at grid point
c      print *,'Inverse'
c Manu's routine ..slow!
c      call inverse(DD,PtIn,Ptin,DDinv,issing)
c      call dot(GD,DDinv,PtOut,PtIn,PtIn,GD_DD)
c      call dot(GD_DD,datap,PtOut,PtIn,1,dataout)

c LAPACK routine ...:)
      call ftranspose(GD,PtOut,PtIn,GDt)
	call DGESV( PtIn, PtOut, DD, PtIn, ipiv, GDt, PtIn, INFO )
c      print *,'Transpose', INFO
      call ftranspose(GDt,PtIn,PtOut,GD_DD)
      call fdot(GD_DD,datap,PtOut,PtIn,1,dataout)
c  A=GD/DD;
c  dhat=A*d;

      call addMean(xhat,yhat,dataout,PtOut,m)

c compute errormap
	do i=1,PtOut
	  agd=0
	  do j=1,PtIn
	     agd=agd + GD_DD(i,j) * GD(i,j)
	  enddo
	  errout(i)=GG(i) - agd
	enddo

	return
	end




c========================================================
c	remove mean
c========================================================
      subroutine removeMean(x,y,datain,Pt,m,datap)
    	implicit none
    	integer Pt, ii, jj,i,j,r,c
    	real*8 x(Pt),y(Pt),datain(Pt),m(3),G(Pt,3),datap(Pt),
     &       dmean(Pt),Gt(3,Pt), b(3,1),Gp(3,3),A(3,3),issing

       ! setup matrix G
          r=Pt
    	    c=3

          do i=1,r
            G(i,1)=x(i)
            G(i,2)=y(i)
            G(i,3)=1.d0
          enddo

          ! transpose Gt
          call ftranspose(G,r,c,Gt)

          ! b=G'*d
          call fdot(Gt,datain ,c,r,1, b)
          ! Gp=G'*G
          call fdot(Gt,G,c,r,c,Gp)
c  MANU
          do i=1,c
             Gp(i,i) = Gp(i,i) + 0.01
          enddo
          ! A=inv(Gp)
          call inverse(Gp,c,c,A,issing)

          if (issing .eq. 1) then
              print*, 'Remove Mean: error in computing inverse'
          else
              !m=A*b
              call fdot(A,b,c,c,1,m)
          !dmean=G*m
    	      call fdot(G,m,r,c,1,dmean)
          endif

      	do i=1,Pt
      	  datap(i)=datain(i) - dmean(i)
      	enddo


      	return
      	end

c========================================================
c	add mean
c========================================================
      subroutine addMean(x,y,datain,Pt,m)
    	implicit none
    	integer Pt, ii, jj,i,j,r,c
    	real*8 x(Pt),y(Pt),datain(Pt),m(3),G(Pt,3)
      real*8 dmean(Pt),Gt(3,Pt), b(3,1),Gp(3,3),A(3,3),issing

! setup matrix G
      r=Pt
      c=3

      do i=1,r
        G(i,1)=x(i)
        G(i,2)=y(i)
        G(i,3)=1.0d0
      enddo

! compute dmean=G*m
      call fdot(G,m,r,c,1,dmean)
      do ii=1,Pt
        datain(ii)=datain(ii) + dmean(ii)
      enddo

      return
      end























! Interpolation routines

! - FIT PLANE

      subroutine fitplane(X,Y,V_,rx,ry,dhat)
      implicit none
      integer r,c,in,i,j,issing
      parameter (r=4,c=3,in=1)
      real*8 X(r),Y(r),G(r,c),Gt(c,r)
      real*8 rx,ry,dmean, dhat,V_(4)
      real*8 d(r), b(c),Gp(c,c), A(c,c),m(c)

    ! setup matrix G
      dmean = 0.d0
      do i=1,r
        G(i,1)=X(i)
        G(i,2)=Y(i)
        G(i,3)=1.d0
        dmean = dmean + V_(i)
      enddo
      dmean = dmean/4.
      do i=1,r
        d(i)=V_(i) - dmean
      enddo

      call ftranspose(G,r,c,Gt)
      call fdot(Gt,d,c,r,1,b)
      call fdot(Gt,G,c,r,c,Gp)

    ! add noise to the diagonal
    !  do i=1,c
    !    Gp(i,i) = Gp(i,i) + 0.01
    !  enddo

      call inverse(Gp,c,c,A,issing)
      if (issing .eq. 1) then
          dhat=999999.0
      else
        call fdot(A,b,c,c,in,m)
        dhat = dmean + rx*m(1) + ry*m(2) + m(3)
      endif

      return
      end



!-------------------------------------------------------------------------
!
! Manu's linear algebra selection!
!
!
! - INVERSE (using row-echelon form)
!
      subroutine inverse(A,r,c,Inv,issing)
      integer r,c
      real*8 A(r,c),AI(r,c+c),Inv(r,c)
      real*8 singtest,issing

    ! generate AI
      do i=1,r
       do j=1,c
         AI(i,j)=A(i,j)
         AI(i,c+j)=0.
        enddo
        AI(i,c+i)=1.
      enddo

      issing=0

    ! reduce echelon form for AI
       call rref(AI,r,c+c)

    ! assign result to I
       do i=1,r
        singtest=0
        do j=1,c
          Inv(i,j)=AI(i,c+j)
          singtest=singtest+AI(i,j)
        enddo
        if (abs(singtest) .lt. 1.0e-15) then
          print*, 'matrix is singular'
          issing=1
        endif
       enddo
       return
       end
!
! - PRINT MATRIX
!
      subroutine printm(A,r,c)
      integer r,c
      real*8 A(r,c)
      write(*,*) '  '
      do i=1,r
         write(*,'(10f12.4)') (A(i,j),j=1,c)
       enddo
      write(*,*) '  '
      return
      end
!
! - ROW-ECHELON REDUCTION
!
      subroutine rref(A,r,c)
      integer r,c,i,j,m,n,k,ii,jj
      real*8 A(r,c),p,rtmp,tol

    ! Loop over the entire matrix.
      i = 1
      j = 1
      m=r
      n=c
      tol=1.0e-14

10    continue

      if ((i .le. m) .and. (j .le. n)) then

    ! Find value and index of largest element in the remainder of column j
    ! (find the mximum absolute value in the first row
    ! and assign to k the row number)

      p=0
      k=0
      do ii=i,m
        if (abs(A(ii,j)) .gt. p) then
           p=abs(A(ii,j))
           k=ii
        endif
      enddo
      k=k+1-1
      if (p .le. tol) then
    ! The column is negligible, zero it out.
        do ii=i,m
          A(ii,j) = 0
          j=j+1
         enddo
      else

    ! Swap i-th and k-th rows.
       do jj=j,n
       rtmp=A(i,jj)
         A(i,jj)=A(k,jj)
         A(k,jj)=rtmp
       enddo

    ! Divide the pivot row by the pivot element.
       rtmp=A(i,j)
       do jj=j,n
         A(i,jj) = A(i,jj)/rtmp
       enddo

    ! Subtract multiples of the pivot row from all the other rows.
       do k = 1,i-1
         rtmp=A(k,j)
         do jj=j,n
           A(k,jj) = A(k,jj) - rtmp*A(i,jj)
         enddo
       enddo
       do k = i+1,m
         rtmp=A(k,j)
         do jj=j,n
           A(k,jj) = A(k,jj) - rtmp*A(i,jj)
         enddo
       enddo

       i=i+1
       j=j+1
      endif
      goto 10
      endif

      return
      end

!
! - DOT PRODUCT  a(r,ri)*m(ri,c) = b(r,c)
!
      subroutine fdot(a,m,r,ri,c,b)
      integer r,ri,c
      real*8 a(r,ri), m(ri,c), b(r,c)
      integer i,j

      do i=1,r
       do j=1,c
        b(i,j)=0.0
        do ii=1,ri
          b(i,j)=b(i,j)+a(i,ii)*m(ii,j)
        enddo
       enddo
      enddo
      return
      end
!
! - TRANSPOSE  a=matrix, r,c=rows, columns t=transpose
!
      subroutine ftranspose(a,r,c,t)
      integer r,c
      real*8 a(r,c), t(c,r)
      integer i,j
      do i=1,r
       do j=1,c
        t(j,i)=a(i,j)
       enddo
      enddo
      return
      end
      SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGESV computes the solution to a real system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
*  The LU decomposition with partial pivoting and row interchanges is
*  used to factor A as
*     A = P * L * U,
*  where P is a permutation matrix, L is unit lower triangular, and U is
*  upper triangular.  The factored form of A is then used to solve the
*  system of equations A * X = B.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N coefficient matrix A.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          The pivot indices that define the permutation matrix P;
*          row i of the matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS matrix of right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
*                has been completed, but the factor U is exactly
*                singular, so the solution could not be computed.
*
*  =====================================================================
*
*     .. External Subroutines ..
      EXTERNAL           DGETRF, DGETRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGESV ', -INFO )
         RETURN
      END IF
*
*     Compute the LU factorization of A.
*
      CALL DGETRF( N, N, A, LDA, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL DGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,
     $                INFO )
      END IF
      RETURN
*
*     End of DGESV
*
      END
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, JP
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M )
     $         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGETF2
*
      END
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A' * X = B.
*
*        Solve U'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve L'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of DGETRS
*
      END
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
*  Further Details
*  ===============
*
*  Modified by
*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  Purpose
*  =======
*
*  IEEECK is called from the ILAENV to verify that Infinity and
*  possibly NaN arithmetic is safe (i.e. will not trap).
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies whether to test just for inifinity arithmetic
*          or whether to test for infinity and NaN arithmetic.
*          = 0: Verify infinity arithmetic only.
*          = 1: Verify infinity and NaN arithmetic.
*
*  ZERO    (input) REAL
*          Must contain the value 0.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  ONE     (input) REAL
*          Must contain the value 1.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  RETURN VALUE:  INTEGER
*          = 0:  Arithmetic failed to produce the correct answers
*          = 1:  Arithmetic produced the correct answers
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. Executable Statements ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      RETURN
      END
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR and QZ methods
*               for nonsymmetric eigenvalue problems.
*          = 9: maximum size of the subproblems at the bottom of the
*               computation tree in the divide-and-conquer algorithm
*               (used by xGELSD and xGESDD)
*          =10: ieee NaN arithmetic can be trusted not to trap
*          =11: infinity arithmetic can be trusted not to trap
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER            IEEECK
      EXTERNAL           IEEECK
*     ..
*     .. Executable Statements ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000,
     $        1100 ) ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  900 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
 1000 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
*
 1100 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
*
*     End of ILAENV
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END

      SUBROUTINE DSWAP (N, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DSWAP
C***PURPOSE  Interchange two vectors.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A5
C***TYPE      DOUBLE PRECISION (SSWAP-S, DSWAP-D, CSWAP-C, ISWAP-I)
C***KEYWORDS  BLAS, INTERCHANGE, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DX  input vector DY (unchanged if N .LE. 0)
C       DY  input vector DX (unchanged if N .LE. 0)
C
C     Interchange double precision DX and double precision DY.
C     For I = 0 to N-1, interchange  DX(LX+I*INCX) and DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DSWAP
      INTEGER N, INCX, INCY
      DOUBLE PRECISION DX(*), DY(*), DTEMP1, DTEMP2, DTEMP3
C***FIRST EXECUTABLE STATEMENT  DSWAP
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP1 = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 3.
C
   20 M = MOD(N,3)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   30 CONTINUE
      IF (N .LT. 3) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP1 = DX(I)
        DTEMP2 = DX(I+1)
        DTEMP3 = DX(I+2)
        DX(I) = DY(I)
        DX(I+1) = DY(I+1)
        DX(I+2) = DY(I+2)
        DY(I) = DTEMP1
        DY(I+1) = DTEMP2
        DY(I+2) = DTEMP3
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   70 CONTINUE
      RETURN
      END
      SUBROUTINE DSCAL (N, DA, DX, INCX)
C***BEGIN PROLOGUE  DSCAL
C***PURPOSE  Multiply a vector by a constant.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A6
C***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DSCAL
      DOUBLE PRECISION DA, DX(*)
      INTEGER I, INCX, IX, M, MP1, N
C***FIRST EXECUTABLE STATEMENT  DSCAL
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE DGER (M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
C***BEGIN PROLOGUE  DGER
C***PURPOSE  Perform the rank 1 operation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B4
C***TYPE      DOUBLE PRECISION (DGER-D)
C***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
C***AUTHOR  Dongarra, J. J., (ANL)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C  DGER   performs the rank 1 operation
C
C     A := alpha*x*y' + A,
C
C  where alpha is a scalar, x is an m element vector, y is an n element
C  vector and A is an m by n matrix.
C
C  Parameters
C  ==========
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( m - 1)*abs( INCX)).
C           Before entry, the incremented array X must contain the m
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients. On exit, A is
C           overwritten by the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C***REFERENCES  Dongarra, J. J., Du Croz, J., Hammarling, S., and
C                 Hanson, R. J.  An extended set of Fortran basic linear
C                 algebra subprograms.  ACM TOMS, Vol. 14, No. 1,
C                 pp. 1-17, March 1988.
C***ROUTINES CALLED  XERBLA
C***REVISION HISTORY  (YYMMDD)
C   861022  DATE WRITTEN
C   910605  Modified to meet SLATEC prologue standards.  Only comment
C           lines were modified.  (BKS)
C***END PROLOGUE  DGER
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C***FIRST EXECUTABLE STATEMENT  DGER
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of DGER  .
C
      END



      SUBROUTINE DTRSM (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $   B, LDB)
C***BEGIN PROLOGUE  DTRSM
C***PURPOSE  Solve one of the matrix equations.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B6
C***TYPE      DOUBLE PRECISION (STRSM-S, DTRSM-D, CTRSM-C)
C***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
C***AUTHOR  Dongarra, J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S. (NAG)
C***DESCRIPTION
C
C  DTRSM  solves one of the matrix equations
C
C     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
C
C  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
C  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
C
C     op( A ) = A   or   op( A ) = A'.
C
C  The matrix X is overwritten on B.
C
C  Parameters
C  ==========
C
C  SIDE   - CHARACTER*1.
C           On entry, SIDE specifies whether op( A ) appears on the left
C           or right of X as follows:
C
C              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
C
C              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
C
C           Unchanged on exit.
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix A is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANSA - CHARACTER*1.
C           On entry, TRANSA specifies the form of op( A ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSA = 'N' or 'n'   op( A ) = A.
C
C              TRANSA = 'T' or 't'   op( A ) = A'.
C
C              TRANSA = 'C' or 'c'   op( A ) = A'.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit triangular
C           as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of B. M must be at
C           least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of B.  N must be
C           at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
C           zero then  A is not referenced and  B need not be set before
C           entry.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
C           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
C           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
C           upper triangular part of the array  A must contain the upper
C           triangular matrix  and the strictly lower triangular part of
C           A is not referenced.
C           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
C           lower triangular part of the array  A must contain the lower
C           triangular matrix  and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
C           A  are not referenced either,  but are assumed to be  unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
C           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
C           then LDA must be at least max( 1, n ).
C           Unchanged on exit.
C
C  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
C           Before entry,  the leading  m by n part of the array  B must
C           contain  the  right-hand  side  matrix  B,  and  on exit  is
C           overwritten by the solution matrix  X.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in  the  calling  (sub)  program.   LDB  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C***REFERENCES  Dongarra, J., Du Croz, J., Duff, I., and Hammarling, S.
C                 A set of level 3 basic linear algebra subprograms.
C                 ACM TOMS, Vol. 16, No. 1, pp. 1-17, March 1990.
C***ROUTINES CALLED  LSAME, XERBLA
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910605  Modified to meet SLATEC prologue standards.  Only comment
C           lines were modified.  (BKS)
C***END PROLOGUE  DTRSM
C     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
C
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C***FIRST EXECUTABLE STATEMENT  DTRSM
C
C     Test the input parameters.
C
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
C
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
C     And when  alpha.eq.zero.
C
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
C
C     Start the operations.
C
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
C
C           Form  B := alpha*inv( A )*B.
C
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
C
C           Form  B := alpha*inv( A' )*B.
C
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
C
C           Form  B := alpha*B*inv( A ).
C
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
C
C           Form  B := alpha*B*inv( A' ).
C
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of DTRSM .
C
      END



      INTEGER FUNCTION IDAMAX (N, DX, INCX)
C***BEGIN PROLOGUE  IDAMAX
C***PURPOSE  Find the smallest index of that component of a vector
C            having the maximum magnitude.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A2
C***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C   IDAMAX  smallest index (zero if N .LE. 0)
C
C     Find smallest index of maximum magnitude of double precision DX.
C     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  IDAMAX
      DOUBLE PRECISION DX(*), DMAX, XMAG
      INTEGER I, INCX, IX, N
C***FIRST EXECUTABLE STATEMENT  IDAMAX
      IDAMAX = 0
      IF (N .LE. 0) RETURN
      IDAMAX = 1
      IF (N .EQ. 1) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increments not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DMAX = ABS(DX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
        XMAG = ABS(DX(IX))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increments equal to 1.
C
   20 DMAX = ABS(DX(1))
      DO 30 I = 2,N
        XMAG = ABS(DX(I))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
   30 CONTINUE
      RETURN
      END



      SUBROUTINE DGEMM (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $   BETA, C, LDC)
C***BEGIN PROLOGUE  DGEMM
C***PURPOSE  Perform one of the matrix-matrix operations.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B6
C***TYPE      DOUBLE PRECISION (SGEMM-S, DGEMM-D, CGEMM-C)
C***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
C***AUTHOR  Dongarra, J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S. (NAG)
C***DESCRIPTION
C
C  DGEMM  performs one of the matrix-matrix operations
C
C     C := alpha*op( A )*op( B ) + beta*C,
C
C  where  op( X ) is one of
C
C     op( X ) = X   or   op( X ) = X',
C
C  alpha and beta are scalars, and A, B and C are matrices, with op( A )
C  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANSA - CHARACTER*1.
C           On entry, TRANSA specifies the form of op( A ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSA = 'N' or 'n',  op( A ) = A.
C
C              TRANSA = 'T' or 't',  op( A ) = A'.
C
C              TRANSA = 'C' or 'c',  op( A ) = A'.
C
C           Unchanged on exit.
C
C  TRANSB - CHARACTER*1.
C           On entry, TRANSB specifies the form of op( B ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSB = 'N' or 'n',  op( B ) = B.
C
C              TRANSB = 'T' or 't',  op( B ) = B'.
C
C              TRANSB = 'C' or 'c',  op( B ) = B'.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry,  M  specifies  the number  of rows  of the  matrix
C           op( A )  and of the  matrix  C.  M  must  be at least  zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry,  N  specifies the number  of columns of the matrix
C           op( B ) and the number of columns of the matrix C. N must be
C           at least zero.
C           Unchanged on exit.
C
C  K      - INTEGER.
C           On entry,  K  specifies  the number of columns of the matrix
C           op( A ) and the number of rows of the matrix op( B ). K must
C           be at least  zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
C           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
C           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
C           part of the array  A  must contain the matrix  A,  otherwise
C           the leading  k by m  part of the array  A  must contain  the
C           matrix A.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
C           LDA must be at least  max( 1, m ), otherwise  LDA must be at
C           least  max( 1, k ).
C           Unchanged on exit.
C
C  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
C           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
C           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
C           part of the array  B  must contain the matrix  B,  otherwise
C           the leading  n by k  part of the array  B  must contain  the
C           matrix B.
C           Unchanged on exit.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
C           LDB must be at least  max( 1, k ), otherwise  LDB must be at
C           least  max( 1, n ).
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
C           supplied as zero then C need not be set on input.
C           Unchanged on exit.
C
C  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
C           Before entry, the leading  m by n  part of the array  C must
C           contain the matrix  C,  except when  beta  is zero, in which
C           case C need not be set on entry.
C           On exit, the array  C  is overwritten by the  m by n  matrix
C           ( alpha*op( A )*op( B ) + beta*C ).
C
C  LDC    - INTEGER.
C           On entry, LDC specifies the first dimension of C as declared
C           in  the  calling  (sub)  program.   LDC  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C***REFERENCES  Dongarra, J., Du Croz, J., Duff, I., and Hammarling, S.
C                 A set of level 3 basic linear algebra subprograms.
C                 ACM TOMS, Vol. 16, No. 1, pp. 1-17, March 1990.
C***ROUTINES CALLED  LSAME, XERBLA
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910605  Modified to meet SLATEC prologue standards.  Only comment
C           lines were modified.  (BKS)
C***END PROLOGUE  DGEMM
C     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C***FIRST EXECUTABLE STATEMENT  DGEMM
C
C     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
C     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
C     and  columns of  A  and the  number of  rows  of  B  respectively.
C
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
C
C     Test the input parameters.
C
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     And if  alpha.eq.zero.
C
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
C
C     Start the operations.
C
      IF( NOTB )THEN
         IF( NOTA )THEN
C
C           Form  C := alpha*A*B + beta*C.
C
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B + beta*C
C
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
C
C           Form  C := alpha*A*B' + beta*C
C
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B' + beta*C
C
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of DGEMM .
C
      END





