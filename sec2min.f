c23456
c     convert OHP magnetometer binary file to ascii file
c
c      character*50 fmt, inpfile, oupfile
c
      character*50 oupfile, oupfile2, oupfile3
      real*8 hx0(86400), hy0(86400), hz0(86400), fp0(1440)
      real*8 tts0(86400), ttc0(86400)
      real*8 hxm(86400), hym(86400), hzm(86400), fpm(86400)
      real*8 ttsm(86400), ttcm(86400)
      real*8 hxp(86400), hyp(86400), hzp(86400), fpp(1440)
      real*8 ttsp(86400), ttcp(86400)
      real*8 xhx(86600), xhy(86600), xhz(86600)
      real*8 xtts(86600), xttc(86600)
c
      real*8 cfil(0:45)
c     
      real*8 sp, critic
      integer clen,nnn,istart(1000),iend(0:1000)
      integer mm, istsp(1000), iensp(1000)
c
      real*8 hxom(1440), hyom(1440), hzom(1440) 
      real*8 ttsom(1440), ttcom(1440)
c
c      real*8 tts(86400), tc(86400)
c
      isec=0
c
      write(*,*) 'start day'
      write(*,*) 'year?'
      read(*,*) nyrs
      write(*,*) 'month?'
      read(*,*) nmons
      write(*,*) 'day?'
      read(*,*) ndays
      write(*,*) 'number of days?'
      read(*,*) idays
c
      nyr=nyrs-2000
      nmon=nmons
c
      write(oupfile,'(a3,2i2.2,a4)') 'mut',nyr,nmon,'.min'
      write(*,*) oupfile
c
      write(oupfile2,'(a1,2i2.2,a4)') 'p',nyr,nmon,'.min'
      write(*,*) oupfile2
c
      write(oupfile3,'(a3,2i2.2,a4)') 'hdz',nyr,nmon,'.min'
      write(*,*) oupfile3
      
c
      open(3,file=oupfile,form='formatted',status='new')
      open(13,file=oupfile2,form='formatted',status='new')
      open(23,file=oupfile3,form='formatted',status='new')
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
      ibhx=39000
      ibhy=0
      ibhz=11000
c
      ndatx=86600
      pi=acos(-1.0)
      r2m=180.*60./pi
c
c  day -1 
c
      nday=ndays-1
      nmon=nmons
      nyr=nyrs-2000
c23456
      call mmm(nyr,nmon,nday,fpm,hxm,hym,hzm,
     @  ttsm,ttcm,ibhx,ibhy,ibhz)
c
c   day 0
c
       nday0=ndays
       nmon0=nmons
       nyr0=nyrs-2000
      call mmm(nyr0,nmon0,nday0,fp0,hx0,hy0,hz0,
     @  tts0,ttc0,ibhx,ibhy,ibhz)
c       
      do id=2, idays+1
       nday=ndays+id-1
       nmon=nmons
       nyr=nyrs-2000
      call mmm(nyr,nmon,nday,fpp,hxp,hyp,hzp,
     @  ttsp,ttcp,ibhx,ibhy,ibhz)
c
      call msdata(hxm,hx0,hxp,xhx)
      call msdata(hym,hy0,hyp,xhy)
      call msdata(hzm,hz0,hzp,xhz)
      call msdata(ttsm,tts0,ttsp,xtts)
      call msdata(ttcm,ttc0,ttcp,xttc)
c
c-- spike check
c
      write(*,*) 'Hx'
      sp=99999.9
      critic=150.
      call spike(xhx,ndatx,sp,critic,clen,istart,iend,nnn,
     @   istsp,iensp,mm)
      call tscorr(xhx,ndatx,istsp,iensp,mm)
c
      write(*,*) 'Hy'
      critic=150.
      call spike(xhy,ndatx,sp,critic,clen,istart,iend,nnn,
     @   istsp,iensp,mm)
      call tscorr(xhy,ndatx,istsp,iensp,mm)
c
      write(*,*) 'Hz'
      critic=150.
      call spike(xhz,ndatx,sp,critic,clen,istart,iend,nnn,
     @ istsp,iensp,mm)
      call tscorr(xhz,ndata,istsp,iensp,mm)
c
      sp=999.9
      write(*,*) 'Ts'
      critic=10.
      call spike(xtts,ndatx,sp,critic,clen,istart,iend,nnn,
     @   istsp,iensp,mm)
      call tscorr(xtts,ndatx,istsp,iensp,mm)
      write(*,*) 'Tc'
      critic=10.
      call spike(xttc,ndatx,sp,critic,clen,istart,iend,nnn,
     @   istsp,iensp,mm)
      call tscorr(xttc,ndatx,istsp,iensp,mm)
c
c-- calculate one-min values w/ interamagnet filter
c
      call omfilt(xhx,ndatx,hxom,cfil)
      call omfilt(xhy,ndatx,hyom,cfil)
      call omfilt(xhz,ndatx,hzom,cfil)
      call omfilt(xtts,ndatx,ttsom,cfil)
      call omfilt(xttc,ndatx,ttcom,cfil)
c
      do itm=1, 1440
        itm0=itm-1
        call i2t(itm0,ihr,imin)
        if(abs(hxom(itm)).gt.1000.) hxom(itm)=999.99
        if(abs(hyom(itm)).gt.1000.) hyom(itm)=999.99
        if(abs(hzom(itm)).gt.1000.) hzom(itm)=999.99
      write(3,100) nyr0,nmon0,nday0,ihr,imin,isec,
     @    hxom(itm),hyom(itm),hzom(itm),fp0(itm),
     @    ttsom(itm),ttcom(itm), ibhx, ibhy, ibhz
c
        shx=hxom(itm)+ibhx
        shy=hyom(itm)+ibhy
        shz=hzom(itm)+ibhz
        sdd=r2m*atan(shy/shx)
        shh=sqrt(shx**2+shy**2)
      write(23,115) nyr0,nmon0,nday0,ihr,imin,isec,
     @    shh, sdd, shz, fp0(itm)
c
         timp=float(nday0)+float(ihr)/24.+float(imin)/1440.
        write(13, 200) timp, shx, shy, shz, shh, sdd, fp0(itm),
     @    ttsom(itm), ttcom(itm) 
       enddo
c
      nyr0=nyr
      nmon0=nmon
      nday0=nday
      do kk=1,86400
       hxm(kk)=hx0(kk)
       hym(kk)=hy0(kk)
       hzm(kk)=hz0(kk)
       ttsm(kk)=tts0(kk)
       ttcm(kk)=ttc0(kk)
       hx0(kk)=hxp(kk)
       hy0(kk)=hyp(kk)
       hz0(kk)=hzp(kk)
       tts0(kk)=ttsp(kk)
       ttc0(kk)=ttcp(kk)
      enddo
      do kk=1, 1440
       fp0(kk)=fpp(kk)
      enddo
c
      enddo
c
c
      close(3)
      close(13)
100     format(I2.2,x,I2.2,x,I2.2,x,I2.2,x,I2.2,x,I2.2,x,f8.3,x,
     @         f8.3,x,f8.3,x,f9.1,x,f8.3,x,f8.3,x,i7,x,i7,x,i7)
115     format(I2.2,x,I2.2,x,I2.2,x,I2.2,x,I2.2,x,I2.2,x,f10.3,x,
     @         f10.3,x,f10.3,x,f9.1)
200     format(f15.8,x,f10.3,x,f10.3,x,f10.3,x,f10.3,x,f10.3,x,
     @   f9.1,x,f8.3,x,f8.3)
111     format(i3,f12.9)
c23456
c
      stop
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c23456
      subroutine mmm(nyr,nmon,nday,fp,hx,hy,hz,
     @     tts,ttc,ibhx,ibhy,ibhz)
c
      character*50 fmt, inpfile
      integer*1 iyr,imn,idy,ihr,imi,isc
      integer*2 itemp,iox,ioy,ioz,ihx,ihy,ihz,itc
      integer*4 iff
      integer inums, ihour, imin, isec 
c
      real*8 hx(86400), hy(86400), hz(86400), fp(1440)
      real*8 tts(86400),ttc(86400)
c
      parameter (fmt='(6a1,7a2,1a4,a2)')
      ipac=26
c
      call calender(nyr,nmon,nday)
c
      write(inpfile,'(a1,3i2.2)') 'd',nyr,nmon,nday
      write(*,*) inpfile
c
      do ii=1, 86400
       hx(ii)=99999.99
       hy(ii)=99999.99
       hz(ii)=99999.99
       tts(ii)=999.99
       ttc(ii)=999.99
      enddo
      do ii=1, 1440
        fp(ii)=99999.9
      enddo
c
      open(2,file=inpfile,access='direct',form='formatted',
     @   recl=ipac,status='old')
c
      do j=1, 86420
        irecn=j
       read(2,fmt,rec=irecn,err=99) iyr,imn,idy,ihr,imi,isc,
     @   itemp,iox,ioy,ioz,ihx,ihy,ihz,iff,itc
c
         ihour=ihr
         imin=imi
         isec=isc
c
         inums=ihour*3600+imin*60+isec+1
c
         hx(inums)=float(iox*10)+float(ihx)/100.-ibhx
         hy(inums)=float(ioy*10)+float(ihy)/100.-ibhy
         hz(inums)=float(ioz*10)+float(ihz)/100.-ibhz
         tts(inums)=float(itemp)/100.
         ttc(inums)=float(itc)/100.
       if(iff.gt.100) then
         inumm=ihour*60+imin+1
         fp(inumm)=float(iff)/10.
       else
       endif
c
      enddo
99    close(2)
c
        write(*,*) iyr, imn, idy 
       return
       end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine i2t(in,ihr,imin)
      ihr=in/60
      imin=in-ihr*60
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
      Real*8  sum
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
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine msdata(hxm,hx0,hxp,xhx)
      real*8 hxm(86400), hx0(86400), hxp(86400), xhx(86600)
c
      do i=1,99
       xhx(i)=hxm(86301+i)
      enddo
      do i=1, 86400
       xhx(i+99)=hx0(i)
      enddo
      do i=1,100
       xhx(i+86499)=hxp(i)
      enddo
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine omfilt(xhx,ndatx,hxom,cfil)
      real*8 xhx(ndatx),hxom(1440),cfil(0:45)
c
      do ii=1, 1440
        itpoint=100+(ii-1)*60       
        xom=cfil(0)*xhx(itpoint)
         do is=1,45
          xom=xom+cfil(is)*(xhx(itpoint-is)+xhx(itpoint+is))
         enddo 
        hxom(ii)=xom
      enddo
      return
      end
