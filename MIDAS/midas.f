      PROGRAM MIDAS
c
c   Author: Geoff Blewitt.  Copyright (C) 2015.
c
c   (1) take all possible pairs of data seperated by precisely one year.
c   (2) if a corresponding pair cannot be found, match with the next 
c       closest epoch that has not yet been used. 
c       2015-07-17: can have unused epochs; particularly bad for campaigns.
c       Solution: make algorithm symmetric, forward and back in time.
c       2015-10-10: MIDAS4 has option to use a step epoch file.  Pairs are
c       not allowed to cross epochs. Daily solution of step epoch is ignored.
c       Modify 
c   (3) find the median, v50, and median of absolute deviation (MAD)
c   (4) trim distribution for data > 2 std dev
c   (5) recompute median and MAD
c   (6) estimate standard error in median using MAD
c
c   Read data from MIDAS.TENV in tenv2 format, 
c    or from MIDAS.TENU in free format: station, time, east, north, up
c   Optionally read step epochs from MIDAS.STEPS
c   Write one line of statistics out to stdout, and write
c   to MIDAS.RENV a new tenv2 file containing residuals to the fit. 
c
      implicit none

c   following is the label and version number for this fit
      character*6, parameter :: label = 'MIDAS4'
      integer, parameter :: maxm=9999, minm=6, minn=10, maxn=19999
      integer, parameter :: maxstep=99 
c     setting tol=0.001d+0 forces pairs separated by precisely 1 year
      real*8,  parameter :: tol=0.001d+0, tolneg = -tol
      real*8,  parameter :: deltmin = 1.d0+real(minm)/365.25+tol
      real*8,  parameter :: svmax=9.999999d+0, rmax=99.999999d+0

      integer i, j, k, m, n1, n2, n, mgood, ndata, ios, nstep, istep
      logical freeformat
      integer ip(2,maxn) 
      real*8 delt, fdt, dt
      integer ne, nn, nu
      real*8 ve50, vn50, vu50, de50, dn50, du50, ce, cn, cu
      real*8 sve, svn, svu, sde, sdn, sdu
      real*8 fe, fn, fu
      real*8 ts, tstep(maxstep+1), tbstep(maxstep+1)
      real*8 t(maxm), tb(maxm)
      real*8 ve(maxn), vn(maxn), vu(maxn)
      real*8 de(maxn), dn(maxn), du(maxn)
      logical good(maxm)
      real*8 xe50, xn50, xu50 
      character*80 line
      character*4 sta, stepsta
      character*7 date(maxm)
      real*8 xe(maxm), xn(maxm), xu(maxm)
      real*8 se(maxm), sn(maxm), su(maxm), ah(maxm)
      real*8 c1(maxm), c2(maxm), c3(maxm) 
      real*8 re(maxm), rn(maxm), ru(maxm)
      integer mjd(maxm), mgpsw(maxm), mwday(maxm)
      real*8 qmedian


c-----DATA INPUT SECTION---------------------------------------------

      open(unit=91,action='write',file='MIDAS.ERR',form='formatted',
     +     iostat=ios, status='new')
      open(unit=51,action='read',file='MIDAS.TENV',form='formatted',
     +     iostat=ios, status='old') 
      if(ios.eq.0) then
         freeformat = .FALSE.
      else
         open(unit=51,action='read',file='MIDAS.TENU',form='formatted',
     +     iostat=ios, status='old') 
         if(ios.ne.0) then
           write(91,*) 'FATAL: Input file not found'
           call exit(1)
         endif
         freeformat = .TRUE.
      endif


      m = 1
      if(freeformat) then
        do
          read(51,*,iostat=ios) sta, t(m), xe(m), xn(m), xu(m)
          if(ios.ne.0) exit
          m = m + 1  
          if(m > maxm) then
             write(91,'(2a,i4)') sta,' WARNING: m epochs > ',maxm
             exit
          endif
        enddo
      else
        do
          read(51,*,iostat=ios) sta, date(m), t(m), mjd(m), mgpsw(m), 
     +      mwday(m), xe(m), xn(m), xu(m), ah(m), se(m), sn(m), su(m),
     +      c1(m), c2(m), c3(m)
          if(ios.ne.0) exit
          if(xe(m)>99.9) cycle
          if(xn(m)>99.9) cycle
          if(xu(m)>99.9) cycle
          if(se(m)>0.03) cycle
          if(sn(m)>0.03) cycle
          if(su(m)>0.1) cycle
          m = m + 1  
          if(m > maxm) then
             write(91,'(2a,i4)') sta,' WARNING: m epochs >= ',maxm
             exit
          endif
        enddo
      endif
      m = m - 1
      close(51)

      if ( m < 1 )  then
         write(91,'(2a)') sta,' FATAL: No data'
         call exit(1)
      endif

      if ( m < 2 )  then
         write(91,'(2a)') sta,' FATAL: No pairs'
         call exit(1)
      endif

c   Compute time span
      delt = t(m) - t(1)
      if ( delt < deltmin )  then
         write(91,'(2a,f5.3,a,f5.3)') sta,' FATAL: Time span(yr) ',
     +      delt,' < ',deltmin
         call exit(1)
      endif

      if ( m < minm ) then
         write(91,'(2a,i1)') sta,' FATAL: m epochs < ',minm
         call exit(1)
      endif

      nstep = 0
      open(unit=52,action='read',file='MIDAS.STEPIN',form='formatted',
     +     iostat=ios, status='old')
      if(ios.eq.0) then
         open(unit=62,action='write',file='MIDAS.STEPOUT',
     +      form='formatted',iostat=ios, status='new')
         do 
           read(52,'(a80)',iostat=ios) line
           if(ios.ne.0) exit
           read(line,*,iostat=ios) stepsta, ts
           if(ios.ne.0) exit

c          only use steps from this station within the data span
           if(stepsta.ne.sta) cycle
           if(ts<t(1) .or. ts>t(m) ) cycle

           if(nstep>=maxstep) then
             write(91,'(2a,i4)') sta,' WARNING: nstep > ',maxstep
             exit
           endif
           nstep = nstep + 1
           tstep(nstep) = ts
           write(62,'(a)') trim(line)
         enddo
      endif

c-----DATA SELECTION SECTION-----------------------------------------

c   select pairs forward in positive time
      call selectpair(m,maxn,tol,t,n1,ip,nstep,tstep)

c   now select pairs backward in negative time to ensure 
c     the algorithm is time symmetric
      call tback(m,t,tb,nstep,tstep,tbstep)
      call selectpair(m,maxn,tol,tb,n2,ip(1,n1+1),nstep,tbstep)
c
c     correct pair indices for forward positive time
      n = n1+n2
      do k = n1+1, n
         i = ip(1,k)
         j = ip(2,k)
         ip(1,k) = m+1-j
         ip(2,k) = m+1-i
      enddo

c   initialize data points contributing to median
      do i = 1, m
         good(i) = .FALSE.
      enddo

c   compute velocity data
      do k = 1, n
         i = ip(1,k)
         j = ip(2,k)
         good(i) = .TRUE.
         good(j) = .TRUE.
         dt = t(j) - t(i)

c         debug/scrutinize pair selection
c         print'(2i6,2a8,3f10.4)',i,j,date(i),date(j),t(i),t(j),dt

         ve(k) = (xe(j)-xe(i))/dt
         vn(k) = (xn(j)-xn(i))/dt
         vu(k) = (xu(j)-xu(i))/dt
      enddo

      if ( n >= maxn ) then
         write(91,'(2a,i5)') sta,' WARNING: n pairs >= ',maxn
      endif
      if ( n < minn )  then
         write(91,'(2a,i2)') sta,' FATAL: n pairs < ',minn
         call exit(1)
      endif

c   Number of data points used
      mgood = 0
      do i = 1, m
         if ( good(i) ) then
            mgood = mgood + 1
         endif
      enddo

c-----ESTIMATION SECTION---------------------------------------------

c   Median of the distribution
      ve50 = qmedian(n,ve)
      vn50 = qmedian(n,vn)
      vu50 = qmedian(n,vu)

c   Absolute velocity deviation from the median
      do i = 1, n
         de(i) = abs(ve(i)-ve50)
         dn(i) = abs(vn(i)-vn50)
         du(i) = abs(vu(i)-vu50)
      enddo
      
c   Median Absolute Deviation (MAD)
      de50 = qmedian(n,de)
      dn50 = qmedian(n,dn)
      du50 = qmedian(n,du)

c   Estimated standard deviation of velocities
      sde = 1.4826*de50
      sdn = 1.4826*dn50
      sdu = 1.4826*du50

c   Delete velocities more than 2 std deviations from median
      ce = 2.0*sde
      cn = 2.0*sdn
      cu = 2.0*sdu
      ne = 0
      nn = 0
      nu = 0
      do i = 1, n
         if (de(i) < ce) then
            ne = ne + 1
            ve(ne) = ve(i)
         endif
         if (dn(i) < cn) then
            nn = nn + 1
            vn(nn) = vn(i)
         endif
         if (du(i) < cu) then
            nu = nu + 1
            vu(nu) = vu(i)
         endif
      enddo

c     Recompute median
      ve50 = qmedian(ne,ve)
      vn50 = qmedian(nn,vn)
      vu50 = qmedian(nu,vu)

c     Recompute absolute deviations
      do i = 1, ne
         de(i) = abs(ve(i)-ve50)
      enddo
      do i = 1, nn
         dn(i) = abs(vn(i)-vn50)
      enddo
      do i = 1, nu
         du(i) = abs(vu(i)-vu50)
      enddo
      
c     Recompute MAD     
      de50 = qmedian(ne,de)
      dn50 = qmedian(nn,dn)
      du50 = qmedian(nu,du)

c     Estimated standard deviation of velocities
c     Multiply by theoretical factor of 1.4826 
      sde = 1.4826*de50
      sdn = 1.4826*dn50
      sdu = 1.4826*du50

c     Standard errors for the median velocity
c     Multiply by theoretical factor of sqrt(pi/2) = 1.2533
c     Divide number of data by 4 since we use coordinate data a nominal 4 times
      sve = 1.2533*sde/sqrt(real(ne)/4.0)
      svn = 1.2533*sdn/sqrt(real(nn)/4.0)
      svu = 1.2533*sdu/sqrt(real(nu)/4.0)

c     Scale standard errors by ad hoc factor of 3 to be realistic
      sve = 3.0*sve
      svn = 3.0*svn
      svu = 3.0*svu

c   Compute intercept at first epoch
c   Intercepts: xe50, xn50, xu50
      do i = 1, m
         dt = t(i) - t(1)
         re(i) = xe(i) - ve50*dt
         rn(i) = xn(i) - vn50*dt
         ru(i) = xu(i) - vu50*dt
      enddo
      xe50 = qmedian(m,re)
      xn50 = qmedian(m,rn)
      xu50 = qmedian(m,ru)
c
c   Compute residuals
      do i = 1, m
         re(i) = re(i) - xe50
         rn(i) = rn(i) - xn50
         ru(i) = ru(i) - xu50
      enddo

c-----OUTPUT RESULTS SECTION-----------------------------------------

c   Protect against format errors and identify gross problems
      if ( ve50 > svmax .or. ve50 < -svmax ) ve50 = svmax
      if ( vn50 > svmax .or. vn50 < -svmax ) vn50 = svmax
      if ( vu50 > svmax .or. vu50 < -svmax ) vu50 = svmax
      if ( sve > svmax ) sve = svmax
      if ( svn > svmax ) svn = svmax
      if ( svu > svmax ) svu = svmax
      if ( sde > svmax ) sde = svmax
      if ( sdn > svmax ) sdn = svmax
      if ( sdu > svmax ) sdu = svmax
      do i = 1, m
        if ( re(i) > rmax .or. re(i) < -rmax ) re(i) = rmax
        if ( rn(i) > rmax .or. rn(i) < -rmax ) rn(i) = rmax
        if ( ru(i) > rmax .or. ru(i) < -rmax ) ru(i) = rmax
      enddo

c   Fraction of pairs removed
      fe = real(n-ne)/real(n)
      fn = real(n-nn)/real(n)
      fu = real(n-nu)/real(n)

c   Write out solution  
      write(6,'(a4,1x,a6,2f10.4,f8.4,i5,i5,i7,
     +        3f10.6,1x,3f9.6,
     +        3f11.6, 3f6.3, 3f9.6, i3)')
     +        sta, label, t(1), t(m), delt, m, mgood, n,
     +        ve50, vn50, vu50, sve, svn, svu,
     +        xe50, xn50, xu50, fe, fn, fu, sde, sdn, sdu, nstep
c
c   Write residuals
      if(freeformat) then
        open(unit=61,action='write',file='MIDAS.RENU',form='formatted',
     +     iostat=ios, status='new')
        if(ios.ne.0) then
           write(91,'(2a)') sta,' FATAL: Cannot open output file'
           call exit(1) 
        endif
        do i = 1, m
          write(61,'(a4,x,f9.4,3f11.6)')
     +       sta, t(i), re(i), rn(i), ru(i)
        enddo
      else
        open(unit=61,action='write',file='MIDAS.RENV',form='formatted',
     +     iostat=ios, status='new')
        if(ios.ne.0) then
           write(91,'(2a)') sta,' FATAL: Cannot open output file'
           call exit(1) 
        endif
        do i = 1, m
          write(61,'(a4,x,a7,x,f9.4,i6,i5,i2,3f11.6,f8.4,3f9.6,3f10.6)')
     +       sta, date(i), t(i), mjd(i), mgpsw(i), mwday(i), 
     +       re(i), rn(i), ru(i), ah(i), se(i), sn(i), su(i), 
     +       c1(i), c2(i), c3(i)
        enddo
      endif
      close(61)

      end
c--------------------------------------------------------------
      subroutine selectpair(m,maxn,tol,t,n,ip,nstep,tstep)
c
c     Given a time tag array t(m), select pairs ip(2,n)
c
c     Moves forward in time: for each time tag, pair it with only
c     one future time tag.
c     First attempt to form a pair within tolerance tol of 1 year.
c     If this fails, then find next unused partner.
c     If this fails, cycle through all possible future partners again.
c
c     MIDAS calls this twice -- firstly forward in time, and
c     secondly backward in time with negative tags and data.
c     This ensures a time symmetric solution.

c     2010-10-12: now allow for apriori list of step epochs
c     - do not select pairs that span or include the step epoch

c     input
      integer m, maxn, nstep
      real*8 tol,t(m), tstep(nstep+1)

c     output
      integer n,ip(2,maxn)

c     local
      integer i,j,k,i2,istep
      real*8 dt, fdt
       
      k = 0
      n = 0
      istep = 1
      do i = 1, m
         if( n >= maxn ) exit
         if(t(i) > (t(m)+tol-1.0)) exit

c        scroll through steps until next step time is later than epoch 1
         do
            if(istep>nstep) exit
            if(t(i) < tstep(istep)+tol) exit
            istep = istep + 1
         enddo 
         if(istep<=nstep) then
            if(t(i) > (tstep(istep)+tol-1.0)) cycle
         endif 

         do j = i+1, m
            if(k<j) k = j
            if(istep<=nstep) then
              if(t(j) > (tstep(istep)-tol)) exit
            endif

            dt = t(j) - t(i)

c         time difference from 1 year
            fdt = (dt - 1.0)

c         keep searching if pair less than one year
            if(fdt<-tol) cycle
c
c         try to find a matching pair within tolerance of 1 year
            if(fdt<tol) then
              i2 = j
            else
c             otherwise, if greater than 1 year, cycle through remaining data
              i2 = k
              dt = t(i2) - t(i)
              if(istep<=nstep) then
                if(t(i2) > (tstep(istep)-tol)) then
                   k = 0
                   cycle
                endif
              endif
              if(k==m) k = 0
              k = k + 1
            endif
c
c         data pair has been found
            n = n + 1
            ip(1,n) = i
            ip(2,n) = i2
            exit
         enddo
      enddo

      end
c
c-----------------------------------------
      subroutine tback(m,t,tb,nstep,tstep,tbstep)
 
c     reverses the order of time array t(m)
c     and multiplies result by -1
c     
     
c     input
      integer m, nstep

c     input
      real*8 t(m), tstep(nstep+1), tbstep(nstep+1)

c     output
      real*8 tb(m)

c     local 
      integer i
      
      do i = 1, m
        tb(m+1-i) = -t(i)
      enddo

      do i = 1, nstep
        tbstep(nstep+1-i) = -tstep(i)
      enddo

      end

c---------------------------------------
      real*8 function qmedian(n,arr)
c
c     Wrapper around "qselect" function to find median.
c     This wrapper ensures array "arr" is not overwritten
c     If n is even, qmedian is the mean of the 2 middle numbers

      integer n, i, k
      integer, parameter :: maxn=19999 
      real*8 arr(n), work(maxn)
      real*8 qselect

      do i = 1, n
        work(i) = arr(i)
      enddo

      k = (n+1)/2
      qmedian = qselect(k,n,work)

      if (k==n/2) then
         qmedian = 0.5d+0*(qmedian + qselect(k+1,n,work))
      endif

      end

c---------------------------------------
      real*8 function qselect(k,n,arr)
c
c   Quickselect algorithm by Hoare [1961].
c   Rewritten in modular form to get rid of confusing goto statements.

c   To find the median, set k = (n+1)/2
c   No account is made here as to whether n is odd or even
c
c   Caution: array "arr" gets overwritten (is "in place")
c   This makes if faster for subsequent selects, but care
c   should be taken if "arr" is to be used externally

      integer k, n
      real*8 arr(n)
      integer i, ir, j, l, mid
      real*8 a, temp

      l = 1
      ir = n

      do while (ir > l+1)
         mid = (l+ir)/2
         temp = arr(mid)
         arr(mid) = arr(l+1)
         arr(l+1) = temp

         if (arr(l) > arr(ir)) then
            temp = arr(l)
            arr(l) = arr(ir)
            arr(ir) = temp
         endif
         if (arr(l+1) > arr(ir)) then
            temp = arr(l+1)
            arr(l+1) = arr(ir)
            arr(ir) = temp
         endif
         if (arr(l) > arr(l+1)) then
            temp = arr(l)
            arr(l) = arr(l+1)
            arr(l+1) = temp
         endif

         i = l+1
         j = ir
         a = arr(i)

         do
            i = i+1
            if (arr(i) < a) cycle
            do
               j = j-1
               if (arr(j) <= a) exit
            enddo
            if (j < i) exit
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
         enddo

         arr(l+1) = arr(j)
         arr(j) = a
         if (j >= k) ir = j-1
         if (j <= k) l = i
      enddo

      if (ir == 2) then
         if (arr(ir) < arr(l)) then
            temp = arr(l)
            arr(l) = arr(ir)
            arr(ir) = temp
         endif
      endif
      qselect = arr(k)

      end
