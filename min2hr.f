c23456
c     convert OHP magnetometer 1min values to 1hr values
c
      character*50 inpfile, inpfile0, oupfile, oupfile2
      character*50 oupfile3
      real*8 hx0(45000), hy0(45000), hz0(45000), fp0(45000)
      real*8 tts0(45000), ttc0(45000)
      real*8 hxm(45000), hym(45000), hzm(45000), fpm(45000)
      real*8 ttsm(45000), ttcm(45000)
      real*8 xhx(45000), xhy(45000), xhz(45000), xfp(45000)
      real*8 xtts(45000), xttc(45000)
      real*8 hxom(750), hyom(750), hzom(750) 
      real*8 ttsom(750), ttcom(750), fpom(750)
c
      real*8 cfil(0:45)
c     
      real*8 sp, critic
      integer clen,nnn,istart(1000),iend(0:1000)
      integer mm, istsp(1000), iensp(1000)
c
      isec=0
      fff=40700.
c
      write(*,*) 'start day'
      write(*,*) 'year?'
      read(*,*) nyrs
      write(*,*) 'month?'
      read(*,*) nmons
c
      nyr=nyrs-2000
      nmon=nmons
c
      if(nmon.eq.1) then
        nmonm=12
        nyrm=nyr-1
        else
         nmonm=nmon-1
         nyrm=nyr
        endif
c
      write(inpfile,'(a3,2i2.2,a4)') 'mut',nyr,nmon,'.min'
      write(*,*) inpfile
c      
      write(inpfile0,'(a3,2i2.2,a4)') 'mut',nyrm,nmonm,'.min'
      write(*,*) inpfile0      
c
      write(oupfile,'(a3,2i2.2,a3)') 'mut',nyr,nmon,'.hr'
      write(*,*) oupfile
      write(oupfile3,'(a3,2i2.2,a3)') 'hdz',nyr,nmon,'.hr'
      write(*,*) oupfile3

c
      write(oupfile2,'(a1,2i2.2,a3)') 'p',nyr,nmon,'.hr'
      write(*,*) oupfile2
c
      open(2,file=inpfile,form='formatted',status='old')
      open(12,file=inpfile0,form='formatted',status='old')
c      
      open(3,file=oupfile,form='formatted',status='new')
      open(13,file=oupfile2,form='formatted',status='new')
      open(23,file=oupfile3,form='formatted',status='new')
c
      pi=acos(-1.0)
      r2m=180.*60./pi
c
c-- set filter coefficients
c
       cfil(0) =  0.02519580
       cfil(1) =  0.02514602
       cfil(2) =  0.02499727
       cfil(3) =  0.02475132
       cfil(4) =  0.02441104
       cfil(5) =  0.02398040
       cfil(6) =  0.02346437
       cfil(7) =  0.02286881
       cfil(8) =  0.02220039
       cfil(9) =  0.02146643
       cfil(10) =  0.02067480
       cfil(11) =  0.01983377
       cfil(12) =  0.01895183
       cfil(13) =  0.01803763
       cfil(14) =  0.01709976
        cfil(15) =  0.01614667
        cfil(16) =  0.01518651
        cfil(17) =  0.01422707
        cfil(18) =  0.01327563
        cfil(19) =  0.01233892
        cfil(20) =  0.01142303
        cfil(21) =  0.01053338
        cfil(22) =  0.00967467
        cfil(23) =  0.00885090
        cfil(24) =  0.00806530
        cfil(25) =  0.00732042
        cfil(26) =  0.00661811
        cfil(27) =  0.00595955
        cfil(28) =  0.00534535
        cfil(29) =  0.00477552
        cfil(30) =  0.00424959
        cfil(31) =  0.00376666
        cfil(32) =  0.00332543
        cfil(33) =  0.00292430
        cfil(34) =  0.00256140
        cfil(35) =  0.00223468
        cfil(36) =  0.00194194
        cfil(37) =  0.00168089
        cfil(38) =  0.00144918
        cfil(39) =  0.00124449
        cfil(40) =  0.00106449
        cfil(41) =  0.00090693
        cfil(42) =  0.00076964
        cfil(43) =  0.00065055
        cfil(44) =  0.00054772
        cfil(45) =  0.00045933     
c
c-- set parameters
c
      clen=15
c
c
c  month -1 
c
      call calender2(nyrm,nmonm,ndaym)
      ndatam=ndaym*1440
c
      do itm = 1, ndatam
      read(12,100) nyr1,nmon1,nday1,ihr,imin,isec,
     @    hxm(itm),hym(itm),hzm(itm),fpm(itm),
     @    ttsm(itm),ttcm(itm), ibhx, ibhy, ibhz
      fpm(itm)=fpm(itm)-fff
      enddo
c
c month 0
c      
      call calender2(nyr,nmon,nday)
      ndata=nday*1440      
      ndatx=ndata+100
      ndathh=nday*24
c
      do itm = 1, ndata
      read(2,100) nyr1,nmon1,nday1,ihr,imin,isec,
     @    hx0(itm),hy0(itm),hz0(itm),fp0(itm),
     @    tts0(itm),ttc0(itm), ibhx, ibhy, ibhz
      fp0(itm)=fp0(itm)-fff
      enddo
c
      call msdata(hxm,hx0,xhx,ndata,ndatam)
      call msdata(hym,hy0,xhy,ndata,ndatam)
      call msdata(hzm,hz0,xhz,ndata,ndatam)
      call msdata(fpm,fp0,xfp,ndata,ndatam)      
      call msdata(ttsm,tts0,xtts,ndata,ndatam)
      call msdata(ttcm,ttc0,xttc,ndata,ndatam)
c
c
c-- spike check
c
      write(*,*) 'Hx'
      sp=99999.9
      critic=50.
      call spike(xhx,ndatx,sp,critic,clen,istart,iend,nnn,
     @   istsp,iensp,mm)
      call tscorr(xhx,ndatx,istsp,iensp,mm)
c
      write(*,*) 'Hy'
      critic=50.
      call spike(xhy,ndatx,sp,critic,clen,istart,iend,nnn,
     @   istsp,iensp,mm)
      call tscorr(xhy,ndatx,istsp,iensp,mm)
c
      write(*,*) 'Hz'
      critic=50.
      call spike(xhz,ndatx,sp,critic,clen,istart,iend,nnn,
     @ istsp,iensp,mm)
      call tscorr(xhz,ndata,istsp,iensp,mm)
c      
      write(*,*) 'F'
      critic=50.
      call spike(xfp,ndatx,sp,critic,clen,istart,iend,nnn,
     @ istsp,iensp,mm)
      call tscorr(xfp,ndata,istsp,iensp,mm)      
c
      sp=999.9
      write(*,*) 'Ts'
      critic=2.
      call spike(xtts,ndatx,sp,critic,clen,istart,iend,nnn,
     @   istsp,iensp,mm)
      call tscorr(xtts,ndatx,istsp,iensp,mm)
      write(*,*) 'Tc'
      critic=2.
      call spike(xttc,ndatx,sp,critic,clen,istart,iend,nnn,
     @   istsp,iensp,mm)
      call tscorr(xttc,ndatx,istsp,iensp,mm)
c
c
c-- calculate one-min values w/ interamagnet filter
c
      call omfilt(xhx,hxom,cfil,ndata)
      call omfilt(xhy,hyom,cfil,ndata)
      call omfilt(xhz,hzom,cfil,ndata)
      call omfilt(xfp,fpom,cfil,ndata)
      call omfilt(xtts,ttsom,cfil,ndata)
      call omfilt(xttc,ttcom,cfil,ndata)
c
      imin=0
      isec=0
      write(*,*) ndathh
c
      do itm=1, ndathh
        itm0=itm-1
        call i2t(itm0,idd,ihh)
         if(abs(hxom(itm)).gt.1000.) hxom(itm)=999.99
         if(abs(hyom(itm)).gt.1000.) hyom(itm)=999.99
         if(abs(hzom(itm)).gt.1000.) hzom(itm)=999.99
        fpom(itm)=fpom(itm)+fff
        idd2=idd+1
      write(3,100) nyr,nmon,idd2,ihh,imin,isec,
     @    hxom(itm),hyom(itm),hzom(itm),fpom(itm),
     @    ttsom(itm),ttcom(itm), ibhx, ibhy, ibhz
c
        shx=hxom(itm)+ibhx
        shy=hyom(itm)+ibhy
        shz=hzom(itm)+ibhz
        shh=sqrt(shx**2+shy**2)
        sdd=r2m*atan(shy/shx)
      write(23,115) nyr,nmon,idd2,ihh,imin,isec,
     @    shh, sdd, shz, fpom(itm)
c
         timp=float(idd2)+float(ihh)/24.
        write(13, 200) timp, shx, shy, shz, shh, sdd, fpom(itm),
     @    ttsom(itm), ttcom(itm) 
       enddo
c
c      enddo
c
c
      close(3)
      close(13)
      close(2)
      close(12)      
100     format(I2.2,x,I2.2,x,I2.2,x,I2.2,x,I2.2,x,I2.2,x,f8.3,x,
     @         f8.3,x,f8.3,x,f9.1,x,f8.3,x,f8.3,x,i7,x,i7,x,i7)
115     format(I2.2,x,I2.2,x,I2.2,x,I2.2,x,I2.2,x,I2.2,x,f10.3,x,
     @         f10.3,x,f10.3,x,f9.1)

200     format(f15.8,x,f10.3,x,f10.3,x,f10.3,x,f10.3,x,f10.3,x,
     @     f9.1,x,f8.3,x,f8.3)
111     format(i3,f12.9)
c23456
c
      stop
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine i2t(in,iday,ihr)
c
       iday=in/24
       ihr=in-iday*24
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine tscorr(xx,ndata,istsp,iensp,mm)
      real*8 xx(ndata)
      integer istsp(1000), iensp(1000) 
c
      if(mm.ge.1) then 
       do iii=1, mm 
         write(*,*) 'istart',istsp(iii), 'iend',iensp(iii)
          if(istsp(iii).eq.1.or.iensp(iii).eq.ndata) write(*,*) 'nya-!' 
           ngintv=iensp(iii)-istsp(iii)+1
            dxx=xx(iensp(iii)+1)-xx(istsp(iii)-1)   
            dxx=dxx/float(ngintv+1) 
             do jt=istsp(iii), iensp(iii)
              xx(jt)=dxx*(jt-float(istsp(iii)-1))+xx(istsp(iii)-1)
              write(*,*) jt, xx(jt)
             enddo
             write(*,*) 'data gap at', istsp(iii), 'corrected'
       enddo
      else
      endif
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine spike(x,ndata,sp,critic,clen,istart,iend,n,
     @                 istsp,iensp,m)
c -------1---------2---------3---------4---------5---------6---------7--
c *******************************************************************
c *  Input                                                          *
c *     x(ndata) : data array                                       *
c *     sp       : special value for data gaps                      *
c *     critic   : maxmum value for a normal difference             *
c *     clen     : maxmum length for a short gap                    *
c *  Output                                                         *
c *     istart(i): increment # of the beginning of i-th big gap     *
c *     iend(i)  : increment # of the end of i-th big gap           *
c *     n        : # of big gaps                                    *
c *******************************************************************
c
      Implicit none
c * Input *
      Integer ndata, clen, m
      Real*8  x(ndata), sp, critic 
c * Variables *
      Integer i,idf,j,flag
      Real*8  df, dummy
      Real*8 sum
c * Output *
      Real  dfmean,dfdev
      Integer n,istart(1000),iend(0:1000),istsp(1000),iensp(1000)
c ***************************************************************
c ************************* Main ********************************
      iend(0) = 0
c --------------------- Compute a mean value of deferences
      sum = 0.D0
      dfmean = 0.D0
      n = 0
      Do 300 i=2, ndata 
        If ( (x(i).ne.sp).and.(x(i-1).ne.sp) ) then
          df = abs( x(i) - x(i-1) )
          sum = sum + df
          n = n + 1
        Endif
  300 Continue
      If ( n .gt. 0 ) dfmean = sum / dble(n)
      Write(6,*) 'mean of deferences=', dfmean, 
     &                       ' ( # of efficient data=',n,')'
c --------------------- Compute a standard deviation of deferences
      sum = 0.D0
      dfdev = 0.D0
      n = 0
      Do 400 i=2, ndata 
        If ( (x(i).ne.sp).and.(x(i-1).ne.sp) ) then
          df = abs( x(i) - x(i-1) )
          sum = sum + (dble(df)-dble(dfmean))**2
          n = n + 1
        Endif
  400 Continue
      If ( n .gt. 0 ) dfdev = sum / dble(n-1)
      Write(6,*) 'dev of deferences=', dfdev, 
     &                       ' ( # of efficient data=',n,')'
c --------------------- Set the thresh hold of a spike
      If ((n.gt.0).and.(critic.gt.dfmean+20*sqrt(dfdev))) then
        critic = dfmean + 20*sqrt(dfdev)
      Endif
c --------------------- Search spikes
      Write(6,*)  'critical value for spike=',critic
      flag = 0
      n = 2
      m=0
      If ( x(1) .eq. sp ) then
 1000   n = n + 1
        If ( x(n-1) .eq. sp ) go to 1000
      Endif
      Do 700 i=n, ndata
        If (flag .eq. 0) df = abs( x(i) - x(i-1) )
        If (flag .eq. 1) df = abs( x(i) - x(istart(1)) )
        If ( (x(i).eq.sp) .or. (df.gt.critic) ) then
          If ( (flag .eq. 1) .and. (i-istart(1) .gt. clen) ) then
            dummy = (dfmean+3.D0*sqrt(dfmean))*dble(i-istart(1))
            If ( df .le. dummy ) flag = 0
          Elseif (flag .eq. 0) then
            flag = 1
            istart(1) = i - 1
          Endif
        Elseif ( flag .eq. 1 ) then
          flag = 0
          iend(1) = i - 1
          idf = iend(1) - istart(1)
          If ( ( idf .le. clen ).and.( istart(1) .ge. 1 ) ) then
            Do 100 j = 1, idf-1
              x(istart(1)+j) = sp
  100       Continue
            Write(6, '(3I7)') istart(1)+1,iend(1),iend(1)-istart(1)
            m=m+1
            istsp(m)=istart(1)+1
            iensp(m)=iend(1)
          Endif
        Endif
  700 Continue
c --------------------- Search sequences of the special value
      Write(6,*)  
      Write(6,*)  'special value for gap=',sp
      flag = 0
      n = 0
      Do 900 i=1, ndata
        If ( x(i).eq.sp ) then
          If ( flag .eq. 0 ) then
            flag = 1
            n = n + 1
            istart(n) = i
          Endif
        Elseif ( flag .eq. 1 ) then
          flag = 0
          iend(n) = i - 1 
          If ( iend(n)-istart(n)+1 .le. clen ) n = n - 1
        Endif
  900 Continue
      If ( flag .eq. 1 ) iend(n) = ndata
      Write(6,*) 'GAP', ' ( # of gap=', n, ')' 
      Do 200 i=1, n
        Write(6, '(3I7)') istart(i),iend(i),iend(i)-istart(i)+1
  200 Continue
c --------------------- End of program
      return
      End

c


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      Subroutine calender(year,month,day)
c
      Implicit none
c
      Integer year,month,day,monthday(12)
      Data monthday / 31,28,31,30,31,30,31,31,30,31,30,31/
c
      If ( mod(year,4) .eq. 0 ) then
        monthday(2) = 29
      Else
        monthday(2) = 28
      Endif
c
      If ( day .lt. 1 ) then
        month = month - 1
        If ( month .lt. 1 ) then
           month = 12
           year = year - 1
        Endif
        day = monthday(month)
      Endif
c
      If ( day .gt. monthday(month) ) then
        day = 1
        month = month + 1
        If ( month .gt. 12 ) then
           month = 1
           year = year + 1
           If ( year .ge. 100 ) year = year - 100
        Endif
      Endif
c
      Return
      End
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      Subroutine calender2(year,month,nday)
c
      Implicit none
c
      Integer year,month,nday,monthday(12)
      Data monthday / 31,28,31,30,31,30,31,31,30,31,30,31/
c
      If ( mod(year,4) .eq. 0 ) then
        monthday(2) = 29
      Else
        monthday(2) = 28
      Endif
c
      nday=monthday(month)
c
      return
      end     
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine msdata(hxm,hx0,xhx,ndata,ndatam)
      real*8 hxm(45000), hx0(45000),xhx(45000)
c
      do i=1,99
       nst=ndatam-99
       xhx(i)=hxm(nst+i)
      enddo
      do i=1, ndata
       xhx(i+99)=hx0(i)
      enddo
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine omfilt(xhx,hxom,cfil,ndata)
      real*8 xhx(45000),hxom(750),cfil(0:45)
c
      ndat=ndata/60
      do ii=1, ndat
        itpoint=100+(ii-1)*60       
        xom=cfil(0)*xhx(itpoint)
         do is=1,45
          xom=xom+cfil(is)*(xhx(itpoint-is)+xhx(itpoint+is))
         enddo 
        hxom(ii)=xom
      enddo
      return
      end
