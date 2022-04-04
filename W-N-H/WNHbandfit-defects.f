c **********************************************************************
c **********************************************************************
c *                   This is the main program                         *
c **********************************************************************
      program bandfitting
      implicit double precision (a-h,o-z)
      parameter(mnum=100,neinum=150)
      parameter(mnp=30,no1=9,no2=1,no3=4)
C WW, HH, WH
      parameter(ngu_WW=77,ngu_HH=6,ngu_WH=19)
c NN, WN, NH
      parameter(ngu_NN=18,ngu_NH=9,ngu_WN=67,ngu_tot=86)
c
      parameter(mnkp=186)
      parameter(nz=41)
      common/con1/vk(3,mnkp,nz),nkp(nz)
      common/con2/alat(3,3,nz),ralat(3,3,nz),tou(mnp,3,nz)
     &      ,itou(mnp,nz),acon(nz),mmm(nz),npa1(nz),npa2(nz),npa3(nz),mm
      common/con3/pi,r_cut,rc_W,rc_H,rc_N,rc_WH,rc_NH,rc_WN,efe
c WW, HH:
c     common/data1/guess_WW(ngu_WW),guess_HH(ngu_HH)
c WH:
c     common/data2/guess_WH(ngu_WH)
c NN, NH:
c     common/data3/guess_NN(ngu_NN),guess_NH(ngu_NH)
c WN
      common/data4/guess_WN(ngu_tot),bund(ngu_tot,2)
c
      common/location/al(3,mnp,mnum,nz),nt(mnp,mnum,nz)
      common/nei1/dxij(mnp,mnum,mnp,nz),dyij(mnp,mnum,mnp,nz)
     &           ,dzij(mnp,mnum,mnp,nz),rij0(mnp,mnum,mnp,nz)
     &           ,nn(mnp,nz),nb(neinum,mnp,nz),nbb(neinum,mnp,nz)
      common/nei2/rjk0(neinum,neinum,mnp,nz),nnc(neinum,mnp,nz)
     &           ,nbc(neinum,neinum,mnp,nz),nbbc(neinum,neinum,mnp,nz)
      common/hmat/hmatr(180,180,2,mnkp)
      common/eig/eigen(180,mnkp,nz)
      common/dewe/dee(84,mnkp,nz),wee(84,mnkp,nz)
      common/ndewe/ndee(nz)
      common/minadata/del,ndiv
      external fe
      dimension x(ngu_tot)
c
       efe=1E10
       call para 
      do nzz=1,nz
       call build(nzz)
       call neighbor(nzz)
       enddo
c
       open(4,file="fitpara.out")
       open(5,file="fitngu.out")
       write(4,*) "Start fitting"
      call mina(fe,ngu_tot,ndiv,del,bund,guess_WN,x,FOFX,IERR)
       last=1
       error=fe(x)
       close(4)
       write(5,142) (guess_WN(i),i=1,ngu_WN)
       write(5,143) (guess_WN(i),i=1+ngu_WN,ngu_tot)
124   format(A4,10f8.4)
142   format(7f10.6)
143   format(4f10.6)
c
      end
c **********************************************************************

c **********************************************************************
c *      Use data to find the hamiltonian matrix elements              *
c **********************************************************************
      double precision function fe(guess)
      implicit double precision (a-h,o-z)
      parameter(mnum=100,neinum=150)
      parameter(mnp=30,no1=9,no2=1,no3=4)
      parameter(ngu_WW=77,ngu_HH=6,ngu_WH=19)
c NN, WN, NH
      parameter(ngu_NN=18,ngu_NH=9,ngu_WN=67,ngu_tot=86)
      parameter(mnkp=186)
      parameter(nz=41)
      common/con1/vk(3,mnkp,nz),nkp(nz)
      common/con2/alat(3,3,nz),ralat(3,3,nz),tou(mnp,3,nz)
     &      ,itou(mnp,nz),acon(nz),mmm(nz),npa1(nz),npa2(nz),npa3(nz),mm
      common/con3/pi,r_cut,rc_W,rc_H,rc_N,rc_WH,rc_NH,rc_WN,efe
c WW, HH:
      common/data1/guess_WW(ngu_WW),guess_HH(ngu_HH)
c WH:
      common/data2/guess_WH(ngu_WH)
c NN, NH:
      common/data3/guess_NN(ngu_NN),guess_NH(ngu_NH)
c WN
c     common/data4/guess_WN(ngu_WN)
c
      common/location/al(3,mnp,mnum,nz),nt(mnp,mnum,nz)
      common/nei1/dxij(mnp,mnum,mnp,nz),dyij(mnp,mnum,mnp,nz)
     &           ,dzij(mnp,mnum,mnp,nz),rij0(mnp,mnum,mnp,nz)
     &           ,nn(mnp,nz),nb(neinum,mnp,nz),nbb(neinum,mnp,nz)
      common/nei2/rjk0(neinum,neinum,mnp,nz),nnc(neinum,mnp,nz)
     &           ,nbc(neinum,neinum,mnp,nz),nbbc(neinum,neinum,mnp,nz)
      common/hmat/hmatr(180,180,2,mnkp)
      common/eig/eigen(180,mnkp,nz)
      common/dewe/dee(84,mnkp,nz),wee(84,mnkp,nz)
      common/ndewe/ndee(nz)

      dimension guess(ngu_tot)
      dimension disc(nz)
      dimension ddc(200,mnkp,nz)
c
c BCPan:
c       guess_NN(21)=guess(29)
c       guess_NN(22)=guess(30)
c
      do nzz=1,nz
         call hmatrix(guess,nzz)
         call c_band(hmatr,nzz)
      enddo
c
      do i=0,nz
       disc(i)=0.0d0
         write(19,*)'test 1 ',i,nkp(i),ndee(i)
       enddo
      fe=0.0d0
c
c****************** For W1N1 1.00 *****************************
       do i=2, ndee(1)+1
          ii=i-1
          do k=1,nkp(1)
             ddc(ii,k,1)=eigen(i,k,1)-eigen(1,k,1)
          enddo
       enddo
 
      do i=1,ndee(1)
         do k=1,nkp(1)
          disc(1)=disc(1)+((ddc(i,k,1)-dee(i,k,1))/0.008)**2*wee(i,k,1)
         enddo
      enddo
       adding=adding+(((dee(1,1,2)-dee(1,1,1))-0.36d0)/0.1)**2
       adding=adding+(((dee(1,1,3)-dee(1,1,1))-0.36d0)/0.1)**2
1212   format(6f10.5)
c
c****************** For W1N1 0.96 *****************************
       do i=2, ndee(2)+1
          ii=i-1
          do k=1,nkp(2)
             ddc(ii,k,2)=eigen(i,k,2)-eigen(1,k,2)
          enddo
       enddo

      do i=1,ndee(2)
         do k=1,nkp(2)
          disc(2)=disc(2)+((ddc(i,k,2)-dee(i,k,2))/0.008)**2*wee(i,k,2)
         enddo
      enddo
c****************** For W1N1 1.04 *****************************
       do i=2, ndee(3)+1
          ii=i-1
          do k=1,nkp(3)
             ddc(ii,k,3)=eigen(i,k,3)-eigen(1,k,3)
          enddo
       enddo

      do i=1,ndee(3)
         do k=1,nkp(3)
          disc(3)=disc(3)+((ddc(i,k,3)-dee(i,k,3))/0.008)**2*wee(i,k,3)
         enddo
      enddo
c
c****************** For W1N2 1.00 *****************************
       dede=0.1d0
       do i=2,ndee(4)+1
          ii=i-1
          do k=1,nkp(4)
             ddc(ii,k,4)=eigen(i,k,4)-eigen(1,k,4)
          enddo
       enddo
 
      do i=1,ndee(4)
         do k=1,nkp(4)
            dede=0.008d0
c           if(k.eq.1.or.k.eq.2.or.k.eq.3)dede=0.01d0
          disc(4)=disc(4)+((ddc(i,k,4)-dee(i,k,4))/dede)**2*wee(i,k,4)
         enddo
      enddo
c****************** For W1N2 0.94 *****************************
       dede=0.1d0
       do i=2,ndee(5)+1
          ii=i-1
          do k=1,nkp(5)
             ddc(ii,k,5)=eigen(i,k,5)-eigen(1,k,5)
          enddo
       enddo

      do i=1,ndee(5)
         do k=1,nkp(5)
            dede=0.008d0
c           if(k.eq.1.or.k.eq.2.or.k.eq.3)dede=0.01d0
          disc(5)=disc(5)+((ddc(i,k,5)-dee(i,k,5))/dede)**2*wee(i,k,5)
         enddo
      enddo
c****************** For W1N2 1.04 *****************************
       dede=0.1d0
       do i=2,ndee(6)+1
          ii=i-1
          do k=1,nkp(6)
             ddc(ii,k,6)=eigen(i,k,6)-eigen(1,k,6)
          enddo
       enddo

      do i=1,ndee(6)
         do k=1,nkp(6)
            dede=0.008d0
c           if(k.eq.1.or.k.eq.2.or.k.eq.3)dede=0.01d0
          disc(6)=disc(6)+((ddc(i,k,6)-dee(i,k,6))/dede)**2*wee(i,k,6)
         enddo
      enddo
c****************** For W2N1 1.00 *****************************
       do i=2,ndee(7)+1
          ii=i-1
          do k=1,nkp(7)
             ddc(ii,k,7)=eigen(i,k,7)-eigen(1,k,7)
          enddo
       enddo
 
      do i=1,ndee(7)
         do k=1,nkp(7)
          disc(7)=disc(7)+((ddc(i,k,7)-dee(i,k,7))/0.01d0)**2*wee(i,k,7)
         enddo
      enddo
c****************** For W2N1 0.96 *****************************
       do i=2,ndee(8)+1
          ii=i-1
          do k=1,nkp(8)
             ddc(ii,k,8)=eigen(i,k,8)-eigen(1,k,8)
          enddo
       enddo

      do i=1,ndee(8)
         do k=1,nkp(8)
          disc(8)=disc(8)+((ddc(i,k,8)-dee(i,k,8))/0.01d0)**2*wee(i,k,8)
         enddo
      enddo
c****************** For W2N1 1.04 *****************************
       do i=2,ndee(9)+1
          ii=i-1
          do k=1,nkp(9)
             ddc(ii,k,9)=eigen(i,k,9)-eigen(1,k,9)
          enddo
       enddo

      do i=1,ndee(9)
         do k=1,nkp(9)
          disc(9)=disc(9)+((ddc(i,k,9)-dee(i,k,9))/0.01d0)**2*wee(i,k,9)
         enddo
      enddo
c****************** For W3N3 1.00 *******************************
       do i=2,ndee(10)+1
          ii=i-1
          do k=1,nkp(10)
             ddc(ii,k,10)=eigen(i,k,10)-eigen(1,k,10)
          enddo
       enddo
 
      do i=1,ndee(10)
         do k=1,nkp(10)
         disc(10)=disc(10)+
     &            ((ddc(i,k,10)-dee(i,k,4))/0.001d0)**2*wee(i,k,10)
         enddo
      enddo
c****************** For W3N3 0.96 *******************************
       do i=2,ndee(11)+1
          ii=i-1
          do k=1,nkp(11)
             ddc(ii,k,11)=eigen(i,k,11)-eigen(1,k,11)
          enddo
       enddo

      do i=1,ndee(11)
         do k=1,nkp(11)
         disc(11)=disc(11)+
     &            ((ddc(i,k,11)-dee(i,k,11))/0.001d0)**2*wee(i,k,11)
         enddo
      enddo
c****************** For W3N3 1.04 *******************************
       do i=2,ndee(12)+1
          ii=i-1
          do k=1,nkp(12)
             ddc(ii,k,12)=eigen(i,k,12)-eigen(1,k,12)
          enddo
       enddo

      do i=1,ndee(12)
         do k=1,nkp(12)
         disc(12)=disc(12)+
     &           ((ddc(i,k,12)-dee(i,k,12))/0.001d0)**2*wee(i,k,12)
         enddo
      enddo
c****************** For W4N2 1.00 *******************************
       do i=2,ndee(13)+1
          ii=i-1
          do k=1,nkp(13)
             ddc(ii,k,13)=eigen(i,k,13)-eigen(1,k,13)
          enddo
       enddo
 
      do i=1,ndee(13)
         do k=1,nkp(13)
        disc(13)=disc(13)+
     &           ((ddc(i,k,13)-dee(i,k,13))/0.0005d0)**2*wee(i,k,13)
         enddo
      enddo
c****************** For W4N2 0.96 *******************************
       do i=2,ndee(14)+1
          ii=i-1
          do k=1,nkp(14)
             ddc(ii,k,14)=eigen(i,k,14)-eigen(1,k,14)
          enddo
       enddo

      do i=1,ndee(14)
         do k=1,nkp(14)
        disc(14)=disc(14)+
     &           ((ddc(i,k,14)-dee(i,k,14))/0.0005d0)**2*wee(i,k,14)
         enddo
      enddo
c****************** For W4N2 1.04 *******************************
       do i=2,ndee(15)+1
          ii=i-1
          do k=1,nkp(15)
             ddc(ii,k,15)=eigen(i,k,15)-eigen(1,k,15)
          enddo
       enddo

      do i=1,ndee(15)
         do k=1,nkp(15)
        disc(15)=disc(15)+
     &           ((ddc(i,k,15)-dee(i,k,15))/0.0005d0)**2*wee(i,k,15)
         enddo
      enddo
c****************** For W2N1H1 *****************************
       do i=2,ndee(16)+1
          ii=i-1
          do k=1,nkp(16)
             ddc(ii,k,16)=eigen(i,k,16)-eigen(1,k,16)
          enddo
       enddo
 
      do i=1,ndee(16)
         do k=1,nkp(16)
          disc(16)=disc(16)+
     &             ((ddc(i,k,16)-dee(i,k,16))/0.01)**2*wee(i,k,16)
         enddo
      enddo
c****************** For W2N1H2 *****************************
       do i=2,ndee(17)+1
          ii=i-1
          do k=1,nkp(17)
             ddc(ii,k,17)=eigen(i,k,17)-eigen(1,k,17)
          enddo
       enddo
 
      do i=1,ndee(17)
         do k=1,nkp(17)
         disc(17)=disc(17)+
     &            ((ddc(i,k,17)-dee(i,k,17))/0.01)**2*wee(i,k,17)
         enddo
      enddo
c****************** For NH3 *******************************
       do i=2,ndee(18)+1
          ii=i-1
          do k=1,nkp(18)
             ddc(ii,k,18)=eigen(i,k,18)-eigen(1,k,18)
          enddo
       enddo
 
      do i=1,ndee(18)
         do k=1,nkp(18)
           disc(18)=disc(18)+(ddc(i,k,18)-dee(i,k,18))**2*wee(i,k,18)
         enddo
      enddo
c        dTB=eigen(2,1,8)-eigen(1,1,8)
c        adding=adding+((dTB-9.90d0)/0.01d0)**2
c        dTB=eigen(4,1,8)-eigen(1,1,8)
c        adding=adding+((dTB-14.99d0)/0.01d0)**2
c****************** For N2 *******************************
       do i=2,ndee(19)+1
          ii=i-1
          do k=1,nkp(19)
             ddc(ii,k,19)=eigen(i,k,19)-eigen(1,k,19)
          enddo
       enddo
 
       do i=1,ndee(19)
         do k=1,nkp(19)
          disc(19)=disc(19)+
     &             ((ddc(i,k,19)-dee(i,k,19))/0.003)**2*wee(i,k,19)
         enddo
      enddo

         dTB=eigen(2,1,19)-eigen(1,1,19)
         adding=adding+((dTB-14.64d0)/0.003d0)**2
         dTB=eigen(3,1,19)-eigen(1,1,19)
         adding=adding+((dTB-16.64d0)/0.003d0)**2
         dTB=eigen(4,1,19)-eigen(1,1,19)
         adding=adding+((dTB-16.64d0)/0.003d0)**2
         dTB=eigen(5,1,19)-eigen(1,1,19)
         adding=adding+((dTB-18.01d0)/0.003d0)**2
c****************** For WH5 ******************************
       do i=2,ndee(20)+1
          ii=i-1
          do k=1,nkp(20)
             ddc(ii,k,20)=eigen(i,k,20)-eigen(1,k,20)
          enddo
       enddo
c
c****************** For W1N1-cluster *********************
       do i=2,ndee(21)+1
          ii=i-1
          do k=1,nkp(21)
             ddc(ii,k,21)=eigen(i,k,21)-eigen(1,k,21)
          enddo
       enddo
 
       do i=1,ndee(21)
         do k=1,nkp(21)
            disc(21)=disc(21)+(ddc(i,k,21)-dee(i,k,21))**2*wee(i,k,21)
         enddo
      enddo
c****************** For W1N2-cluster**********************
       do i=2,ndee(22)+1
          ii=i-1
          do k=1,nkp(22)
             ddc(ii,k,22)=eigen(i,k,22)-eigen(1,k,22)
          enddo
       enddo
 
       do i=1,ndee(22)
         do k=1,nkp(22)
            disc(22)=disc(22)+(ddc(i,k,22)-dee(i,k,22))**2*wee(i,k,22)
         enddo
      enddo
c****************** For W2N1-cluster**********************
       do i=2,ndee(23)+1
          ii=i-1
          do k=1,nkp(23)
             ddc(ii,k,23)=eigen(i,k,23)-eigen(1,k,23)
          enddo
       enddo
 
       do i=1,ndee(23)
         do k=1,nkp(23)
            disc(23)=disc(23)+(ddc(i,k,23)-dee(i,k,23))**2*wee(i,k,23)
         enddo
      enddo
c****************** For W2N1-L 1.00 **********************
       do i=2,ndee(24)+1
          ii=i-1
          do k=1,nkp(24)
             ddc(ii,k,24)=eigen(i,k,24)-eigen(1,k,24)
          enddo
       enddo
 
       do i=1,ndee(24)
         do k=1,nkp(24)
            disc(24)=disc(24)+
     &              ((ddc(i,k,24)-dee(i,k,24))/0.005d0)**2*wee(i,k,24)
         enddo
      enddo
c****************** For W2N1-L 0.96 **********************
       do i=2,ndee(25)+1
          ii=i-1
          do k=1,nkp(25)
             ddc(ii,k,25)=eigen(i,k,25)-eigen(1,k,25)
          enddo
       enddo

       do i=1,ndee(25)
         do k=1,nkp(25)
            disc(25)=disc(25)+
     &              ((ddc(i,k,25)-dee(i,k,25))/0.005d0)**2*wee(i,k,25)
         enddo
      enddo
c****************** For W2N1-L 1.04 **********************
       do i=2,ndee(26)+1
          ii=i-1
          do k=1,nkp(26)
             ddc(ii,k,26)=eigen(i,k,26)-eigen(1,k,26)
          enddo
       enddo

       do i=1,ndee(26)
         do k=1,nkp(26)
            disc(26)=disc(26)+
     &              ((ddc(i,k,26)-dee(i,k,26))/0.005d0)**2*wee(i,k,26)
         enddo
      enddo
c****************** For W6N12 1.00 **********************
       do i=2,ndee(27)+1
          ii=i-1
          do k=1,nkp(27)
             ddc(ii,k,27)=eigen(i,k,27)-eigen(1,k,27)
          enddo
       enddo
 
       do i=1,ndee(27)
         do k=1,nkp(27)
         disc(27)=disc(27)+
     &           ((ddc(i,k,27)-dee(i,k,27))/0.001)**2*wee(i,k,27)
         enddo
      enddo
c****************** For W6N12 0.96 **********************
       do i=2,ndee(28)+1
          ii=i-1
          do k=1,nkp(28)
             ddc(ii,k,28)=eigen(i,k,28)-eigen(1,k,28)
          enddo
       enddo

       do i=1,ndee(28)
         do k=1,nkp(28)
         disc(28)=disc(28)+
     &           ((ddc(i,k,28)-dee(i,k,28))/0.001)**2*wee(i,k,28)
         enddo
      enddo
c****************** For W6N12 1.04 **********************
       do i=2,ndee(29)+1
          ii=i-1
          do k=1,nkp(29)
             ddc(ii,k,29)=eigen(i,k,29)-eigen(1,k,29)
          enddo
       enddo

       do i=1,ndee(29)
         do k=1,nkp(29)
         disc(29)=disc(29)+
     &           ((ddc(i,k,29)-dee(i,k,29))/0.001)**2*wee(i,k,29)
         enddo
      enddo
c****************** For W12N18 1.00 **********************
       do i=2,ndee(30)+1
          ii=i-1
          do k=1,nkp(30)
             ddc(ii,k,30)=eigen(i,k,30)-eigen(1,k,30)
          enddo
       enddo
 
       do i=1,ndee(30)
         do k=1,nkp(30)
         disc(30)=disc(30)+
     &        ((ddc(i,k,30)-dee(i,k,30))/0.005)**2*wee(i,k,30)
         enddo
      enddo
c****************** For W12N18 0.96 **********************
       do i=2,ndee(31)+1
          ii=i-1
          do k=1,nkp(31)
             ddc(ii,k,31)=eigen(i,k,31)-eigen(1,k,31)
          enddo
       enddo

       do i=1,ndee(31)
         do k=1,nkp(31)
         disc(31)=disc(31)+
     &        ((ddc(i,k,31)-dee(i,k,31))/0.005)**2*wee(i,k,31)
         enddo
      enddo
c****************** For W12N18 1.04 **********************
       do i=2,ndee(32)+1
          ii=i-1
          do k=1,nkp(32)
             ddc(ii,k,32)=eigen(i,k,32)-eigen(1,k,32)
          enddo
       enddo

       do i=1,ndee(32)
         do k=1,nkp(32)
         disc(32)=disc(32)+
     &        ((ddc(i,k,32)-dee(i,k,32))/0.005)**2*wee(i,k,32)
         enddo
      enddo
c****************** For W4N6 1.00 **********************
       do i=2,ndee(33)+1
          ii=i-1
          do k=1,nkp(33)
             ddc(ii,k,33)=eigen(i,k,33)-eigen(1,k,33)
          enddo
       enddo

       do i=1,ndee(33)
         do k=1,nkp(33)
         disc(33)=disc(33)+
     &        ((ddc(i,k,33)-dee(i,k,33))/0.005)**2*wee(i,k,33)
         enddo
      enddo
c****************** For W4N6 0.96 **********************
       do i=2,ndee(34)+1
          ii=i-1
          do k=1,nkp(34)
             ddc(ii,k,34)=eigen(i,k,34)-eigen(1,k,34)
          enddo
       enddo

       do i=1,ndee(34)
         do k=1,nkp(34)
         disc(34)=disc(34)+
     &        ((ddc(i,k,34)-dee(i,k,34))/0.005)**2*wee(i,k,34)
         enddo
      enddo
c****************** For W4N6 1.04 **********************
       do i=2,ndee(35)+1
          ii=i-1
          do k=1,nkp(35)
             ddc(ii,k,35)=eigen(i,k,35)-eigen(1,k,35)
          enddo
       enddo

       do i=1,ndee(35)
         do k=1,nkp(35)
         disc(35)=disc(35)+
     &        ((ddc(i,k,35)-dee(i,k,35))/0.005)**2*wee(i,k,35)
         enddo
      enddo
c****************** For W4N2_1*********************
       do i=2,ndee(36)+1
          ii=i-1
          do k=1,nkp(36)
             ddc(ii,k,36)=eigen(i,k,36)-eigen(1,k,36)
          enddo
       enddo

       do i=1,ndee(36)
         do k=1,nkp(36)
         disc(36)=disc(36)+
     &        ((ddc(i,k,36)-dee(i,k,36))/0.005)**2*wee(i,k,36)
         enddo
      enddo
c****************** For W4N2_2*********************
       do i=2,ndee(37)+1
          ii=i-1
          do k=1,nkp(37)
             ddc(ii,k,37)=eigen(i,k,37)-eigen(1,k,37)
          enddo
       enddo

       do i=1,ndee(37)
         do k=1,nkp(37)
         disc(37)=disc(37)+
     &        ((ddc(i,k,37)-dee(i,k,37))/0.005)**2*wee(i,k,37)
         enddo
      enddo
c*******20*********** For W12N18_V_N1*********************
       do i=2,ndee(38)+1
          ii=i-1
          do k=1,nkp(38)
             ddc(ii,k,38)=eigen(i,k,38)-eigen(1,k,38)
          enddo
       enddo

       do i=1,ndee(38)
         do k=1,nkp(38)
         disc(38)=disc(38)+
     &        ((ddc(i,k,38)-dee(i,k,38))/0.005)**2*wee(i,k,38)
         enddo
      enddo
c*******21*********** For W12N18_V_N2*********************
       do i=2,ndee(39)+1
          ii=i-1
          do k=1,nkp(39)
             ddc(ii,k,39)=eigen(i,k,39)-eigen(1,k,39)
          enddo
       enddo

       do i=1,ndee(39)
         do k=1,nkp(39)
         disc(39)=disc(39)+
     &        ((ddc(i,k,39)-dee(i,k,39))/0.005)**2*wee(i,k,39)
         enddo
      enddo
c*******22*********** For W12N18_V_W1*********************
       do i=2,ndee(40)+1
          ii=i-1
          do k=1,nkp(40)
             ddc(ii,k,40)=eigen(i,k,40)-eigen(1,k,40)
          enddo
       enddo

       do i=1,ndee(40)
         do k=1,nkp(40)
         disc(40)=disc(40)+
     &        ((ddc(i,k,40)-dee(i,k,40))/0.005)**2*wee(i,k,40)
         enddo
      enddo

c*******23*********** For W12N18_V_W2*********************
       do i=2,ndee(41)+1
          ii=i-1
          do k=1,nkp(41)
             ddc(ii,k,41)=eigen(i,k,41)-eigen(1,k,41)
          enddo
       enddo

       do i=1,ndee(41)
         do k=1,nkp(41)
         disc(41)=disc(41)+
     &        ((ddc(i,k,41)-dee(i,k,41))/0.005)**2*wee(i,k,41)
         enddo
      enddo

c adding:
c
c        dTB=eigen(2,1,15)-eigen(1,1,15)
c        disnew=disnew+((dTB-0.49d0)/0.001d0)**2

c        dTB=eigen(3,1,15)-eigen(2,1,15)
c        disnew=disnew+((dTB-0.49d0)/0.01d0)**2

c        dTB=eigen(4,1,15)-eigen(3,1,15)
c        disnew=disnew+((dTB-0.49d0)/0.01d0)**2

c        dTB=eigen(15,1,15)-eigen(14,1,15)
c        disnew=disnew+((dTB-0.09d0)/0.001d0)**2
c
c        dTB=eigen(15,1,15)-eigen(13,1,15)
c        disnew=disnew+((dTB-0.09d0)/0.001d0)**2
c        dTB=eigen(15,1,15)-eigen(12,1,15)
c        disnew=disnew+((dTB-0.09d0)/0.001d0)**2
c
c        dTB=eigen(3,2,15)-eigen(2,2,15)
c        disnew=disnew+((dTB-0.04d0)/0.0001d0)**2
c        dTB=eigen(4,2,15)-eigen(2,2,15)
c        disnew=disnew+((dTB-0.04d0)/0.0001d0)**2

c     do k=2,nkp(27)
c        dLDA=dee(47,k,27)-dee(47,1,27)
c        dTB =eigen(48,k,27)-eigen(48,1,27)
c        disnew=disnew+((dTB-0.02d0)/0.002d0)**2
c     enddo
c
c        dTB=eigen(48,1,27)-eigen(47,1,27)
c        adding=adding+((dTB-0.82d0)/0.01d0)**2
c !!!
c        dTB=eigen(1,1,27)-eigen(1,1,5)
c        adding=adding+((dTB-1.223d0)/0.05d0)**2
c
c        dTB=eigen(1,1,15)-eigen(1,1,4)
c        adding=adding+((dTB-2.648d0)/0.05d0)**2

c
c        dTB=eigen(49,1,15)-eigen(48,1,15)
c        disnew=disnew+((dTB-0.72d0)/0.01d0)**2
c
c        dTB=eigen(50,1,15)-eigen(49,1,15)
c        disnew=disnew+((dTB-0.33d0)/0.01d0)**2
 
c        dTB=eigen(52,1,15)-eigen(51,1,15)
c        disnew=disnew+((dTB-0.60d0)/0.01d0)**2
 
c        dTB=eigen(48,29,15)-eigen(47,29,15)
c        disnew=disnew+((dTB-0.10d0)/0.01d0)**2
c
c******************* END **************************************
c 2020-11:
c      write(4,134)eigen(1,1,5),eigen(1,1,3)
c    &                           ,eigen(1,1,4),eigen(1,1,1),eigen(1,1,2)
134    format(5f12.6)
c      adderror1=(((eigen(1,1,3)-eigen(1,1,5))-0.5d0)/0.31)**2
c      adderror2=(((eigen(1,1,4)-eigen(1,1,5))-0.8d0)/0.31)**2
c      adderror3=(((eigen(1,1,1)-eigen(1,1,5))-1.0d0)/0.31)**2
c      adderror4=(((eigen(1,1,2)-eigen(1,1,5))-1.4d0)/0.31)**2
c 2020-11
c      disc(0)=adderror1+adderror2+adderror3+adderror4
       disc(0)=0.0d0+adding 
c      disc(0)=0.0d0
       do i=1,nz
       disc(0)=disc(0)+disc(i)
      enddo
       errorr=disc(0)
       fe=dsqrt(errorr)
c      write(19,*)'fe efe ',fe,efe
      if(fe.lt.efe) then
         efe=fe
        write(4,*)'fitting parameters'
c
        write(4,*)" W1N1 1.0 ",disc(1),eigen(1,1,1)
        do k=1,nkp(1)
           write(4,124)'dee',(dee(i,k,1),i=1,ndee(1))
           write(4,124)'ddc',(ddc(i,k,1),i=1,ndee(1))
           write(4,*)
        enddo
c
        write(4,*)" W1N1 0.9",disc(2),eigen(1,1,2)
        do k=1,nkp(2)
           write(4,124)'dee',(dee(i,k,2),i=1,ndee(2))
           write(4,124)'ddc',(ddc(i,k,2),i=1,ndee(2))
           write(4,*)
        enddo
c
        write(4,*)" W1N1 1.04",disc(3),eigen(1,1,3)
        do k=1,nkp(3)
           write(4,124)'dee',(dee(i,k,3),i=1,ndee(3))
           write(4,124)'ddc',(ddc(i,k,3),i=1,ndee(3))
           write(4,*)
        enddo
c
        write(4,*)" W1N2 1.0 ",disc(4),eigen(1,1,4)
        do k=1,nkp(4)
           write(4,124)'dee',(dee(i,k,4),i=1,ndee(4))
           write(4,124)'ddc',(ddc(i,k,4),i=1,ndee(4))
           write(4,*)
        enddo
c
        write(4,*)" W1N2 0.96",disc(5),eigen(1,1,5)
        do k=1,nkp(5)
           write(4,124)'dee',(dee(i,k,5),i=1,ndee(5))
           write(4,124)'ddc',(ddc(i,k,5),i=1,ndee(5))
           write(4,*)
        enddo
c
        write(4,*)" W1N2 1.04",disc(6),eigen(1,1,6)
        do k=1,nkp(6)
           write(4,124)'dee',(dee(i,k,6),i=1,ndee(6))
           write(4,124)'ddc',(ddc(i,k,6),i=1,ndee(6))
           write(4,*)
        enddo
c
        write(4,*)" W2N1 1.00",disc(7),eigen(1,1,7)
        do k=1,nkp(7)
           write(4,124)'dee',(dee(i,k,7),i=1,ndee(7))
           write(4,124)'ddc',(ddc(i,k,7),i=1,ndee(7))
           write(4,*)
        enddo
c
        write(4,*)" W2N1 0.96",disc(8),eigen(1,1,8)
        do k=1,nkp(8)
           write(4,124)'dee',(dee(i,k,8),i=1,ndee(8))
           write(4,124)'ddc',(ddc(i,k,8),i=1,ndee(8))
           write(4,*)
        enddo
c
        write(4,*)" W2N1 1.04",disc(9),eigen(1,1,9)
        do k=1,nkp(9)
           write(4,124)'dee',(dee(i,k,9),i=1,ndee(9))
           write(4,124)'ddc',(ddc(i,k,9),i=1,ndee(9))
           write(4,*)
        enddo
c
        write(4,*)" W3N3 1.00",disc(10),eigen(1,1,10)
        do k=1,nkp(10)
           write(4,124)'dee',(dee(i,k,10),i=1,ndee(10))
           write(4,124)'ddc',(ddc(i,k,10),i=1,ndee(10))
           write(4,*)
        enddo
c
        write(4,*)" W3N3 0.96",disc(11),eigen(1,1,11)
        do k=1,nkp(11)
           write(4,124)'dee',(dee(i,k,11),i=1,ndee(11))
           write(4,124)'ddc',(ddc(i,k,11),i=1,ndee(11))
           write(4,*)
        enddo
c
        write(4,*)" W3N3 1.04",disc(12),eigen(1,1,12)
        do k=1,nkp(12)
           write(4,124)'dee',(dee(i,k,12),i=1,ndee(12))
           write(4,124)'ddc',(ddc(i,k,12),i=1,ndee(12))
           write(4,*)
        enddo
c
        write(4,*)" W4N2 1.00",disc(13),eigen(1,1,13)
        do k=1,nkp(13)
           write(4,124)'dee',(dee(i,k,13),i=1,ndee(13))
           write(4,124)'ddc',(ddc(i,k,13),i=1,ndee(13))
           write(4,*)
        enddo
c
        write(4,*)" W4N2 0.96",disc(14),eigen(1,1,14)
        do k=1,nkp(14)
           write(4,124)'dee',(dee(i,k,14),i=1,ndee(14))
           write(4,124)'ddc',(ddc(i,k,14),i=1,ndee(14))
           write(4,*)
        enddo
c
        write(4,*)" W4N2 1.04",disc(15),eigen(1,1,15)
        do k=1,nkp(15)
           write(4,124)'dee',(dee(i,k,15),i=1,ndee(15))
           write(4,124)'ddc',(ddc(i,k,15),i=1,ndee(15))
           write(4,*)
        enddo
c
        write(4,*)" W2N1H1 ",disc(16),eigen(1,1,16)
        do k=1,nkp(16)
           write(4,124)'dee',(dee(i,k,16),i=1,ndee(16))
           write(4,124)'ddc',(ddc(i,k,16),i=1,ndee(16))
           write(4,*)
        enddo
c
        write(4,*)" W2N1H2 ",disc(17),eigen(1,1,17)
        do k=1,nkp(17)
           write(4,124)'dee',(dee(i,k,17),i=1,ndee(17))
           write(4,124)'ddc',(ddc(i,k,17),i=1,ndee(17))
           write(4,*)
        enddo
c
        write(4,*)" NH3 ",disc(18)
        do k=1,nkp(18)
           write(4,124)'dee',(dee(i,k,18),i=1,ndee(18))
           write(4,124)'ddc',(ddc(i,k,18),i=1,ndee(18))
           write(4,*)
        enddo
 
        write(4,*)" N2 ",disc(19)
        do k=1,nkp(19)
           write(4,124)'dee',(dee(i,k,19),i=1,ndee(19))
           write(4,124)'ddc',(ddc(i,k,19),i=1,ndee(19))
           write(4,*)
        enddo
 
        write(4,*)" WH5 ",disc(20)
        do k=1,nkp(20)
           write(4,124)'dee',(dee(i,k,20),i=1,ndee(20))
           write(4,124)'ddc',(ddc(i,k,20),i=1,ndee(20))
           write(4,*)
        enddo
c
        write(4,*)" W1N1-c ",disc(21)
        do k=1,nkp(21)
           write(4,124)'dee',(dee(i,k,21),i=1,ndee(21))
           write(4,124)'ddc',(ddc(i,k,21),i=1,ndee(21))
           write(4,*)
        enddo
c
        write(4,*)" W1N2-c ",disc(22)
        do k=1,nkp(22)
           write(4,124)'dee',(dee(i,k,22),i=1,ndee(22))
           write(4,124)'ddc',(ddc(i,k,22),i=1,ndee(22))
           write(4,*)
        enddo
c
        write(4,*)" W2N1-c ",disc(23)
        do k=1,nkp(23)
           write(4,124)'dee',(dee(i,k,23),i=1,ndee(23))
           write(4,124)'ddc',(ddc(i,k,23),i=1,ndee(23))
           write(4,*)
        enddo
c
        write(4,*)" W2N1-L 1.00",disc(24),eigen(1,1,24)
        do k=1,nkp(24)
           write(4,124)'dee',(dee(i,k,24),i=1,20)
           write(4,124)'ddc',(ddc(i,k,24),i=1,20)
           write(4,*)
           write(4,124)'dee',(dee(i,k,24),i=21, ndee(24))
           write(4,124)'ddc',(ddc(i,k,24),i=21, ndee(24))
           write(4,*)
        enddo
c
        write(4,*)" W2N1-L 0.96",disc(25),eigen(1,1,25)
        do k=1,nkp(25)
           write(4,124)'dee',(dee(i,k,25),i=1,20)
           write(4,124)'ddc',(ddc(i,k,25),i=1,20)
           write(4,*)
           write(4,124)'dee',(dee(i,k,25),i=21, ndee(25))
           write(4,124)'ddc',(ddc(i,k,25),i=21, ndee(25))
           write(4,*)
        enddo
c
        write(4,*)" W2N1-L 1.04",disc(26),eigen(1,1,26)
        do k=1,nkp(26)
           write(4,124)'dee',(dee(i,k,26),i=1,20)
           write(4,124)'ddc',(ddc(i,k,26),i=1,20)
           write(4,*)
           write(4,124)'dee',(dee(i,k,26),i=21, ndee(26))
           write(4,124)'ddc',(ddc(i,k,26),i=21, ndee(26))
           write(4,*)
        enddo
c
        write(4,*)" W6N12 1.00",disc(27),eigen(1,1,27)
        do k=1,nkp(27)
           write(4,124)'dee',(dee(i,k,27),i=1,20)
           write(4,124)'ddc',(ddc(i,k,27),i=1,20)
           write(4,*)
           write(4,124)'dee',(dee(i,k,27),i=21,40)
           write(4,124)'ddc',(ddc(i,k,27),i=21,40)
           write(4,*)
           write(4,124)'dee',(dee(i,k,27),i=41,ndee(27))
           write(4,124)'ddc',(ddc(i,k,27),i=41,ndee(27))
           write(4,*)
        enddo
c
        write(4,*)" W6N12 0.96",disc(28),eigen(1,1,28)
        do k=1,nkp(28)
           write(4,124)'dee',(dee(i,k,28),i=1,20)
           write(4,124)'ddc',(ddc(i,k,28),i=1,20)
           write(4,*)
           write(4,124)'dee',(dee(i,k,28),i=21,40)
           write(4,124)'ddc',(ddc(i,k,28),i=21,40)
           write(4,*)
           write(4,124)'dee',(dee(i,k,28),i=41,ndee(28))
           write(4,124)'ddc',(ddc(i,k,28),i=41,ndee(28))
           write(4,*)
        enddo
c
        write(4,*)" W6N12 1.04",disc(29),eigen(1,1,29)
        do k=1,nkp(29)
           write(4,124)'dee',(dee(i,k,29),i=1,20)
           write(4,124)'ddc',(ddc(i,k,29),i=1,20)
           write(4,*)
           write(4,124)'dee',(dee(i,k,29),i=21,40)
           write(4,124)'ddc',(ddc(i,k,29),i=21,40)
           write(4,*)
           write(4,124)'dee',(dee(i,k,29),i=41,ndee(29))
           write(4,124)'ddc',(ddc(i,k,29),i=41,ndee(29))
           write(4,*)
        enddo
c
        write(4,*)" W12N18 1.00",disc(30)
        do k=1,nkp(30)
           write(4,124)'dee',(dee(i,k,30),i=1,20)
           write(4,124)'ddc',(ddc(i,k,30),i=1,20)
           write(4,*)
           write(4,124)'dee',(dee(i,k,30),i=21,40)
           write(4,124)'ddc',(ddc(i,k,30),i=21,40)
           write(4,*)
           write(4,124)'dee',(dee(i,k,30),i=41,60)
           write(4,124)'ddc',(ddc(i,k,30),i=41,60)
           write(4,*)
           write(4,124)'dee',(dee(i,k,30),i=61,80)
           write(4,124)'ddc',(ddc(i,k,30),i=61,80)
           write(4,*)
           write(4,124)'dee',(dee(i,k,30),i=81,ndee(30))
           write(4,124)'ddc',(ddc(i,k,30),i=81,ndee(30))
           write(4,*)
        enddo
c
        write(4,*)" W12N18 0.96",disc(31)
        do k=1,nkp(31)
           write(4,124)'dee',(dee(i,k,31),i=1,20)
           write(4,124)'ddc',(ddc(i,k,31),i=1,20)
           write(4,*)
           write(4,124)'dee',(dee(i,k,31),i=21,40)
           write(4,124)'ddc',(ddc(i,k,31),i=21,40)
           write(4,*)
           write(4,124)'dee',(dee(i,k,31),i=41,60)
           write(4,124)'ddc',(ddc(i,k,31),i=41,60)
           write(4,*)
           write(4,124)'dee',(dee(i,k,31),i=61,80)
           write(4,124)'ddc',(ddc(i,k,31),i=61,80)
           write(4,*)
           write(4,124)'dee',(dee(i,k,31),i=81,ndee(31))
           write(4,124)'ddc',(ddc(i,k,31),i=81,ndee(31))
           write(4,*)
        enddo
c
        write(4,*)" W12N18 1.04",disc(32)
        do k=1,nkp(32)
           write(4,124)'dee',(dee(i,k,32),i=1,20)
           write(4,124)'ddc',(ddc(i,k,32),i=1,20)
           write(4,*)
           write(4,124)'dee',(dee(i,k,32),i=21,40)
           write(4,124)'ddc',(ddc(i,k,32),i=21,40)
           write(4,*)
           write(4,124)'dee',(dee(i,k,32),i=41,60)
           write(4,124)'ddc',(ddc(i,k,32),i=41,60)
           write(4,*)
           write(4,124)'dee',(dee(i,k,32),i=61,80)
           write(4,124)'ddc',(ddc(i,k,32),i=61,80)
           write(4,*)
           write(4,124)'dee',(dee(i,k,32),i=81,ndee(32))
           write(4,124)'ddc',(ddc(i,k,32),i=81,ndee(32))
           write(4,*)
        enddo
c
           write(4,*)" W4N6 1.00",disc(33)
        do k=1,nkp(33)
           write(4,124)'dee',(dee(i,k,33),i=1,20)
           write(4,124)'ddc',(ddc(i,k,33),i=1,20)
           write(4,*)
           write(4,124)'dee',(dee(i,k,33),i=21,ndee(33))
           write(4,124)'ddc',(ddc(i,k,33),i=21,ndee(33))
           write(4,*)
        enddo
c
           write(4,*)" W4N6 0.96",disc(34)
        do k=1,nkp(34)
           write(4,124)'dee',(dee(i,k,34),i=1,20)
           write(4,124)'ddc',(ddc(i,k,34),i=1,20)
           write(4,*)
           write(4,124)'dee',(dee(i,k,34),i=21,ndee(34))
           write(4,124)'ddc',(ddc(i,k,34),i=21,ndee(34))
           write(4,*)
        enddo
c
           write(4,*)" W4N6 1.04",disc(35)
        do k=1,nkp(35)
           write(4,124)'dee',(dee(i,k,35),i=1,20)
           write(4,124)'ddc',(ddc(i,k,35),i=1,20)
           write(4,*)
           write(4,124)'dee',(dee(i,k,35),i=21,ndee(35))
           write(4,124)'ddc',(ddc(i,k,35),i=21,ndee(35))
           write(4,*)
        enddo
c
           write(4,*)"W4N2_1",disc(36)
        do k=1,nkp(36)
           write(4,124)'dee',(dee(i,k,36),i=1,ndee(36))
           write(4,124)'ddc',(ddc(i,k,36),i=1,ndee(36))
c          write(4,*)
c          write(4,124)'dee',(dee(i,k,18),i=21,ndee(18))
c          write(4,124)'ddc',(ddc(i,k,18),i=21,ndee(18))
           write(4,*)
        enddo
        
c
           write(4,*)" W4N2_2",disc(37)
        do k=1,nkp(37)
           write(4,124)'dee',(dee(i,k,37),i=1,ndee(37))
           write(4,124)'ddc',(ddc(i,k,37),i=1,ndee(37))
c          write(4,*)
c          write(4,124)'dee',(dee(i,k,19),i=21,ndee(19))
c          write(4,124)'ddc',(ddc(i,k,19),i=21,ndee(19))
           write(4,*)
        enddo
c
           write(4,*)" V_N1  ",disc(38)
        do k=1,nkp(38)
           write(4,124)'dee',(dee(i,k,38),i=1,20)
           write(4,124)'ddc',(ddc(i,k,38),i=1,20)
           write(4,124)'dee',(dee(i,k,38),i=21,40)
           write(4,124)'ddc',(ddc(i,k,38),i=21,40)
           write(4,124)'dee',(dee(i,k,38),i=41,60)
           write(4,124)'ddc',(ddc(i,k,38),i=41,60)
           write(4,124)'dee',(dee(i,k,38),i=61,ndee(38))
           write(4,124)'ddc',(ddc(i,k,38),i=61,ndee(38))
c          write(4,*)
c          write(4,124)'dee',(dee(i,k,19),i=21,ndee(19))
c          write(4,124)'ddc',(ddc(i,k,19),i=21,ndee(19))
           write(4,*)
        enddo
c
           write(4,*)" V_N2  ",disc(39)
        do k=1,nkp(39)
           write(4,124)'dee',(dee(i,k,39),i=1,20)
           write(4,124)'ddc',(ddc(i,k,39),i=1,20)
           write(4,124)'dee',(dee(i,k,39),i=21,40)
           write(4,124)'ddc',(ddc(i,k,39),i=21,40)
           write(4,124)'dee',(dee(i,k,39),i=41,60)
           write(4,124)'ddc',(ddc(i,k,39),i=41,60)
           write(4,124)'dee',(dee(i,k,39),i=61,ndee(39))
           write(4,124)'ddc',(ddc(i,k,39),i=61,ndee(39))
c          write(4,*)
c          write(4,124)'dee',(dee(i,k,19),i=21,ndee(19))
c          write(4,124)'ddc',(ddc(i,k,19),i=21,ndee(19))
           write(4,*)
        enddo
c
           write(4,*)" V_W1  ",disc(40)
        do k=1,nkp(40)
           write(4,124)'dee',(dee(i,k,40),i=1,20)
           write(4,124)'ddc',(ddc(i,k,40),i=1,20)
           write(4,124)'dee',(dee(i,k,40),i=21,40)
           write(4,124)'ddc',(ddc(i,k,40),i=21,40)
           write(4,124)'dee',(dee(i,k,40),i=41,60)
           write(4,124)'ddc',(ddc(i,k,40),i=41,60)
           write(4,124)'dee',(dee(i,k,40),i=61,ndee(40))
           write(4,124)'ddc',(ddc(i,k,40),i=61,ndee(40))
c          write(4,*)
c          write(4,124)'dee',(dee(i,k,19),i=21,ndee(19))
c          write(4,124)'ddc',(ddc(i,k,19),i=21,ndee(19))
           write(4,*)
        enddo
c
           write(4,*)" V_W2  ",disc(41)
        do k=1,nkp(41)
           write(4,124)'dee',(dee(i,k,41),i=1,20)
           write(4,124)'ddc',(ddc(i,k,41),i=1,20)
           write(4,124)'dee',(dee(i,k,41),i=21,40)
           write(4,124)'ddc',(ddc(i,k,41),i=21,40)
           write(4,124)'dee',(dee(i,k,41),i=41,60)
           write(4,124)'ddc',(ddc(i,k,41),i=41,60)
           write(4,124)'dee',(dee(i,k,41),i=61,ndee(41))
           write(4,124)'ddc',(ddc(i,k,41),i=61,ndee(41))
c          write(4,*)
c          write(4,124)'dee',(dee(i,k,19),i=21,ndee(19))
c          write(4,124)'ddc',(ddc(i,k,19),i=21,ndee(19))
           write(4,*)
        enddo

c        write(4,124)'dee',((dee(48,k,15)-dee(48,1,15)),k=1,nkp(15))
c
         dTB1=eigen(48,1,27)-eigen(47,1,27)
         dTB2=eigen(49,1,27)-eigen(48,1,27)
         dTB3=eigen(50,1,27)-eigen(49,1,27)
         dTB4=eigen(52,1,27)-eigen(51,1,27)

      write(4,125)'ddc ',((eigen(48,k,27)-eigen(48,1,27)),k=1,25)
      write(4,125)'ddc ',((eigen(48,k,27)-eigen(48,1,27)),k=26,nkp(27))
      write(4,124)'27: ',dTB1,dTB2,dTB3,dTB4
c
      write(4,*)  'fe=',fe
c     write(4,142) (guess(i),i=1,ngu_WN)
      write(4,142) (guess(i),i=1,56)
      write(4,142) (guess(i),i=57,60)
      write(4,142) (guess(i),i=61,67)
      write(4,143) (guess(i),i=68,ngu_tot)
c
      endif
124   format(A4,20f6.2)
125   format(A4,30f6.2)
142   format(7f10.6)
143   format(4f10.6)
c
      end
c **********************************************************************

c **********************************************************************
c *      Choose a strucure, and also give some const for compute       *
c **********************************************************************
      subroutine para
      implicit double precision (a-h,o-z)
      parameter(mnum=100,neinum=150)
      parameter(mnp=30,no1=9,no2=1,no3=4)
C WW, HH, WH
      parameter(ngu_WW=77,ngu_HH=6,ngu_WH=19)
c NN, WN, NH
      parameter(ngu_NN=18,ngu_NH=9,ngu_WN=67,ngu_tot=86)
c
      parameter(mnkp=186)
      parameter(nz=41)
      common/con1/vk(3,mnkp,nz),nkp(nz)
      common/con2/alat(3,3,nz),ralat(3,3,nz),tou(mnp,3,nz)
     &      ,itou(mnp,nz),acon(nz),mmm(nz),npa1(nz),npa2(nz),npa3(nz),mm
      common/con3/pi,r_cut,rc_W,rc_H,rc_N,rc_WH,rc_NH,rc_WN,efe
c WW, HH:
      common/data1/guess_WW(ngu_WW),guess_HH(ngu_HH)
c WH:
      common/data2/guess_WH(ngu_WH)
c NN, NH:
      common/data3/guess_NN(ngu_NN),guess_NH(ngu_NH)
c WN
      common/data4/guess_WN(ngu_tot),bund(ngu_tot,2)
c
      common/location/al(3,mnp,mnum,nz),nt(mnp,mnum,nz)
      common/nei1/dxij(mnp,mnum,mnp,nz),dyij(mnp,mnum,mnp,nz)
     &           ,dzij(mnp,mnum,mnp,nz),rij0(mnp,mnum,mnp,nz)
     &           ,nn(mnp,nz),nb(neinum,mnp,nz),nbb(neinum,mnp,nz)
      common/nei2/rjk0(neinum,neinum,mnp,nz),nnc(neinum,mnp,nz)
     &           ,nbc(neinum,neinum,mnp,nz),nbbc(neinum,neinum,mnp,nz)
      common/hmat/hmatr(180,180,2,mnkp)
      common/eig/eigen(180,mnkp,nz)
      common/dewe/dee(84,mnkp,nz),wee(84,mnkp,nz)
      common/ndewe/ndee(nz)
      common/minadata/del,ndiv

      dimension alat1(3,3),aij(3,3),nkp1(nz)
c
c   According to the value of npa(k),we get the number of different kind of atoms 
      open(1,file="para.out")
      open(2,file="kpts.in")
      open(3,file="bandpara.txt")
c
        pi=4.0d0*datan(1.0d0)
        read(3,*)
      do k=1,nz
        read(3,*)  npa1(k),npa2(k),npa3(k),acon(k)
       do i=1,3
        read(3,*)  (alat(i,j,k),j=1,3)
       enddo
c 
       nnpa=npa1(k)+npa2(k)+npa3(k)
c
       do i=1,nnpa
        read(3,*)  (tou(i,j,k),j=1,3),itou(i,k)
       enddo
      enddo
c
        read(3,*)
        read(3,*)   mm,r_cut,rc_W,rc_H,rc_WH
        read(3,*)
        read(3,*)   rc_N,rc_NH,rc_WN
        read(3,*)
        read(3,*)  (guess_WW(i),i=1,ngu_WW)
        read(3,*)
        read(3,*)  (guess_HH(i),i=1,ngu_HH)
        read(3,*)
        read(3,*)  (guess_WH(i),i=1,ngu_WH)
        read(3,*)
        read(3,*)  (guess_NN(i),i=1,ngu_NN)
        read(3,*)
        read(3,*)  (guess_NH(i),i=1,ngu_NH)
        read(3,*)
        read(3,*)  (guess_WN(i),i=1,ngu_tot)
        read(3,*)
c
        do i=1,ngu_tot
          read(3,*) bund(i,1),bund(i,2),ibund
         if(ibund.eq.0) then
          bund(i,1)=guess_WN(i)
          bund(i,2)=guess_WN(i)
         endif
        enddo
        read(3,*)  del,ndiv
c
        do k=1,nz
         read(3,*) ndee(k),nkp(k)
          write(19,*)'test read3',k,ndee(k),nkp(k)
          do kp=1,nkp(k)
              read(3,*)  (dee(i,kp,k),i=1,ndee(k))
              read(3,*)  (wee(i,kp,k),i=1,ndee(k))
              read(3,*)
c             if(k.eq.27)write(19,*)kp,(dee(i,kp,k),i=1,ndee(k))
          enddo
        enddo
c       stop
c
         s=0.0d0
       do k=1,nz
         do i=1,ndee(k)
            do kp=1,nkp(k)
               s=s+wee(i,kp,k)
            enddo
         enddo
       enddo
c     
       do k=1,nz
         do i=1,ndee(k)
            do kp=1,nkp(k)
               wee(i,kp,k)=wee(i,kp,k)/s
            enddo
         enddo
       enddo
c
      do k=1,nz
        do i=1,3
         do j=1,3
            alat(i,j,k)=alat(i,j,k)*acon(k)
         enddo
        enddo
      enddo
c
      do k=1,nz
c
          do i=1,3
           do j=1,3
            alat1(i,j)=alat(i,j,k)
           enddo
          enddo
c
          detA=+alat1(1,1)*(alat1(2,2)*alat1(3,3)-alat1(2,3)*alat1(3,2))
     &         -alat1(1,2)*(alat1(2,1)*alat1(3,3)-alat1(2,3)*alat1(3,1))
     &         +alat1(1,3)*(alat1(2,1)*alat1(3,2)-alat1(2,2)*alat1(3,1))
c
          aij(1,1)=  alat1(2,2)*alat1(3,3)-alat1(2,3)*alat1(3,2)
          aij(1,2)=-(alat1(2,1)*alat1(3,3)-alat1(2,3)*alat1(3,1))
          aij(1,3)=  alat1(2,1)*alat1(3,2)-alat1(2,2)*alat1(3,1)
          aij(2,1)=  alat1(3,2)*alat1(1,3)-alat1(3,3)*alat1(1,2)
          aij(2,2)=-(alat1(3,1)*alat1(1,3)-alat1(3,3)*alat1(1,1))
          aij(2,3)=  alat1(3,1)*alat1(1,2)-alat1(3,2)*alat1(1,1)
          aij(3,1)=  alat1(1,2)*alat1(2,3)-alat1(1,3)*alat1(2,2)
          aij(3,2)=-(alat1(1,1)*alat1(2,3)-alat1(1,3)*alat1(2,1))
          aij(3,3)=  alat1(1,1)*alat1(2,2)-alat1(1,2)*alat1(2,1)
c
          do i=1,3
           do j=1,3
            ralat(i,j,k)=2.0d0*3.14159265358979d0*aij(i,j)/detA
           enddo
          enddo
c
      enddo
c
      do k=1,nz
c 
       nnpa=npa1(k)+npa2(k)+npa3(k)
c
c     if(k.le.12) then
       do ii=1,nnpa
       rx=tou(ii,1,k)
       ry=tou(ii,2,k)
       rz=tou(ii,3,k)
       tou(ii,1,k)=rx*alat(1,1,k)+ry*alat(2,1,k)+rz*alat(3,1,k)
       tou(ii,2,k)=rx*alat(1,2,k)+ry*alat(2,2,k)+rz*alat(3,2,k)
       tou(ii,3,k)=rx*alat(1,3,k)+ry*alat(2,3,k)+rz*alat(3,3,k)
       enddo
c     endif
      enddo
c
      do k=1,nz
        read(2,*)nkp1(k)
        write(19,*)'test read2 ',k,nkp1(k)
        if(nkp1(k).ne.nkp(k))stop 'nkp1'
        
        do i=1,nkp(k)
           read(2,*)(vk(j,i,k),j=1,3)
c          write(19,*)i,vk(1,i,k),vk(2,i,k),vk(3,i,k)
        enddo
      enddo
c
        write(1,*)  "The parameter of energyband program"
      do k=1,nz
        write(1,9)  npa1(k),npa2(k),npa3(k),1.000
       do i=1,3
        write(1,10)  (alat(i,j,k),j=1,3)
       enddo
          nnpa=npa1(k)+npa2(k)+npa3(k)
       do i=1,nnpa
        write(1,11)  (tou(i,j,k),j=1,3),itou(i,k)
       enddo
      enddo
9       format(3I4,f13.6)
10      format(3f13.7)
11      format(3f13.7,4x,I1)
c
        write(1,*)  "The number of repeat cell, r_cut and tree rc"
        write(1,12)   mm, r_cut,rc_W,rc_H,rc_WH
12      format(I3,4x,4f13.6)
        write(1,*)  "The guess_WW value "
        write(1,22)  (guess_WW(i),i=1,ngu_WW)
        write(1,*)  "The guess_HH value "
        write(1,22)  (guess_HH(i),i=1,ngu_HH)
        write(1,*)  "The guess_WH value "
        write(1,22)  (guess_WH(i),i=1,ngu_WH)
        write(1,*)  "The guess_NN value "
        write(1,22)  (guess_NN(i),i=1,ngu_NN)
        write(1,*)  "The guess_NH value "
        write(1,22)  (guess_NH(i),i=1,ngu_NH)
        write(1,*)  "The guess_WN value "
        write(1,22)  (guess_WN(i),i=1,ngu_tot)
        write(1,*)
22      format(7f13.6)
c
      do j=1,nz
        write(1,*)"The structure is",j
        write(1,*)"The number of kpoints is",nkp(j)
        write(1,*)"The kpoints is as follows"
       do i=1,nkp(j)
        write(1,23)  (vk(i,k,j),k=1,3)
       enddo
      enddo
23     format(3f13.6)
c
      close(1)
      close(2)
      close(3)
      end
c **********************************************************************

c **********************************************************************
c *      Use data in program para to build the structure               *
c **********************************************************************
      subroutine build(nzz) 
      implicit double precision (a-h,o-z)
      parameter(mnum=100,neinum=150)
      parameter(mnp=30,no1=9,no2=1,no3=4)
C WW, HH, WH
      parameter(ngu_WW=77,ngu_HH=6,ngu_WH=19)
c NN, WN, NH
      parameter(ngu_NN=18,ngu_NH=9,ngu_WN=67,ngu_tot=86)
c
      parameter(mnkp=186)
      parameter(nz=41)
      common/con1/vk(3,mnkp,nz),nkp(nz)
      common/con2/alat(3,3,nz),ralat(3,3,nz),tou(mnp,3,nz)
     &     ,itou(mnp,nz),acon(nz),mmm(nz),npa1(nz),npa2(nz),npa3(nz),mm
      common/con3/pi,r_cut,rc_W,rc_H,rc_N,rc_WH,rc_NH,rc_WN,efe
c WW, HH:
      common/data1/guess_WW(ngu_WW),guess_HH(ngu_HH)
c WH:
      common/data2/guess_WH(ngu_WH)
c NN, NH:
      common/data3/guess_NN(ngu_NN),guess_NH(ngu_NH)
c WN
      common/data4/guess_WN(ngu_tot),bund(ngu_tot,2)
c
      common/location/al(3,mnp,mnum,nz),nt(mnp,mnum,nz)
      common/nei1/dxij(mnp,mnum,mnp,nz),dyij(mnp,mnum,mnp,nz)
     &           ,dzij(mnp,mnum,mnp,nz),rij0(mnp,mnum,mnp,nz)
     &           ,nn(mnp,nz),nb(neinum,mnp,nz),nbb(neinum,mnp,nz)
      common/nei2/rjk0(neinum,neinum,mnp,nz),nnc(neinum,mnp,nz)
     &           ,nbc(neinum,neinum,mnp,nz),nbbc(neinum,neinum,mnp,nz)
      common/hmat/hmatr(180,180,2,mnkp)
      common/eig/eigen(180,mnkp,nz)
      common/dewe/dee(84,mnkp,nz),wee(84,mnkp,nz)
      common/ndewe/ndee(nz)
      common/minadata/del,ndiv
      dimension alx(mnp),aly(mnp),alz(mnp)
c
       mmm(nzz)=1
       nnpa=npa1(nzz)+npa2(nzz)+npa3(nzz)
      do ii=1,nnpa
       al(1,ii,1,nzz)=tou(ii,1,nzz)
       al(2,ii,1,nzz)=tou(ii,2,nzz)
       al(3,ii,1,nzz)=tou(ii,3,nzz)
       nt(ii,1,nzz)=itou(ii,nzz)
      enddo
c
      do i=-mm,mm
        do j=-mm,mm
          do k=-mm,mm
           if(i.ne.0.or.j.ne.0.or.k.ne.0) then
             jnpa=0
            do ii=1,nnpa
             alx(ii)=i*alat(1,1,nzz)+j*alat(2,1,nzz)
     &               +k*alat(3,1,nzz)+al(1,ii,1,nzz) 
             aly(ii)=i*alat(1,2,nzz)+j*alat(2,2,nzz)
     &               +k*alat(3,2,nzz)+al(2,ii,1,nzz)
             alz(ii)=i*alat(1,3,nzz)+j*alat(2,3,nzz)
     &               +k*alat(3,3,nzz)+al(3,ii,1,nzz)
             do jj=1,nnpa
             dx=alx(ii)-al(1,jj,1,nzz)
             dy=aly(ii)-al(2,jj,1,nzz)
             dz=alz(ii)-al(3,jj,1,nzz)
             rij=dsqrt(dx**2+dy**2+dz**2)
             if(rij.lt.r_cut)  jnpa=jnpa+1
             enddo
            enddo
c
            if(jnpa.ne.0) then
              mmm(nzz)=mmm(nzz)+1
             do jj=1,nnpa
              al(1,jj,mmm(nzz),nzz)=alx(jj)
              al(2,jj,mmm(nzz),nzz)=aly(jj)
              al(3,jj,mmm(nzz),nzz)=alz(jj)
                nt(jj,mmm(nzz),nzz)=itou(jj,nzz)
             enddo
            endif
c
           endif
          enddo
        enddo
      enddo
c
      end
c **********************************************************************

c **********************************************************************
c *      Use data in program build to find the neighborlist of atom 1  *
c **********************************************************************
      subroutine neighbor(nzz) 
      implicit double precision (a-h,o-z)
      parameter(mnum=100,neinum=150)
      parameter(mnp=30,no1=9,no2=1,no3=4)
C WW, HH, WH
      parameter(ngu_WW=77,ngu_HH=6,ngu_WH=19)
c NN, WN, NH
      parameter(ngu_NN=18,ngu_NH=9,ngu_WN=67,ngu_tot=86)
c
      parameter(mnkp=186)
      parameter(nz=41)
      common/con1/vk(3,mnkp,nz),nkp(nz)
      common/con2/alat(3,3,nz),ralat(3,3,nz),tou(mnp,3,nz)
     &      ,itou(mnp,nz),acon(nz),mmm(nz),npa1(nz),npa2(nz),npa3(nz),mm
      common/con3/pi,r_cut,rc_W,rc_H,rc_N,rc_WH,rc_NH,rc_WN,efe
c WW, HH:
      common/data1/guess_WW(ngu_WW),guess_HH(ngu_HH)
c WH:
      common/data2/guess_WH(ngu_WH)
c NN, NH:
      common/data3/guess_NN(ngu_NN),guess_NH(ngu_NH)
c WN
      common/data4/guess_WN(ngu_tot),bund(ngu_tot,2)
c
      common/location/al(3,mnp,mnum,nz),nt(mnp,mnum,nz)
      common/nei1/dxij(mnp,mnum,mnp,nz),dyij(mnp,mnum,mnp,nz)
     &           ,dzij(mnp,mnum,mnp,nz),rij0(mnp,mnum,mnp,nz)
     &           ,nn(mnp,nz),nb(neinum,mnp,nz),nbb(neinum,mnp,nz)
      common/nei2/rjk0(neinum,neinum,mnp,nz),nnc(neinum,mnp,nz)
     &           ,nbc(neinum,neinum,mnp,nz),nbbc(neinum,neinum,mnp,nz)
c     common/hmat/hmatr(no1*mnp,no1*mnp,2,mnkp,nz)
c     common/eig/eigen(no1*mnp,mnkp,nz)
      common/dewe/dee(84,mnkp,nz),wee(84,mnkp,nz)
      common/ndewe/ndee(nz)
c
      nnpa=npa1(nzz)+npa2(nzz)+npa3(nzz)
      do ii=1,nnpa
        nn(ii,nzz)=0
      enddo
c
      do ii=1,nnpa
       do jt=1,neinum
        nnc(jt,ii,nzz)=0
       enddo
      enddo
c
      do ii=1,nnpa
       do j=1,mmm(nzz)
        do jj=1,nnpa
         if(j.ne.1.or.ii.ne.jj) then
          dx=al(1,jj,j,nzz)-al(1,ii,1,nzz)
          dy=al(2,jj,j,nzz)-al(2,ii,1,nzz)
          dz=al(3,jj,j,nzz)-al(3,ii,1,nzz)
          rij=dsqrt(dx**2+dy**2+dz**2)
          dxij(jj,j,ii,nzz)=dx
          dyij(jj,j,ii,nzz)=dy
          dzij(jj,j,ii,nzz)=dz
          rij0(jj,j,ii,nzz)=rij
c
          if(rij.lt.r_cut) then
          nn(ii,nzz)=nn(ii,nzz)+1
          nb(nn(ii,nzz),ii,nzz)=j
          nbb(nn(ii,nzz),ii,nzz)=jj
          endif
         endif
        enddo
       enddo
      enddo

c     do ii=1,nnpa
c      write(19,*)'nz, ii, nn ',nzz,ii, nn(ii,nzz)
c     enddo
c
      do ii=1,nnpa
       do jt=1,nn(ii,nzz)
           j=nb(jt,ii,nzz)
           jj=nbb(jt,ii,nzz)
         do kc=1,nn(ii,nzz)
           k=nb(kc,ii,nzz)
           kk=nbb(kc,ii,nzz)
          if(jt.ne.kc) then
           rjk=dsqrt((al(1,jj,j,nzz)-al(1,kk,k,nzz))**2
     &              +(al(2,jj,j,nzz)-al(2,kk,k,nzz))**2
     &              +(al(3,jj,j,nzz)-al(3,kk,k,nzz))**2)
           if(rjk.lt.r_cut) then
            nnc(jt,ii,nzz)=nnc(jt,ii,nzz)+1
             nbc(nnc(jt,ii,nzz),jt,ii,nzz)=k
            nbbc(nnc(jt,ii,nzz),jt,ii,nzz)=kk
            rjk0(nnc(jt,ii,nzz),jt,ii,nzz)=rjk
c            write(19,*)'test ',nzz,ii,j,rjk
           endif
          endif
         enddo
       enddo
      enddo
c
c     if(nzz.eq.1) then
c       if(nt(1,1,nzz).eq.1) then
c       write(*,*) "atom W",(al(i,1,1,nzz),i=1,3)
c       else
c       write(*,*) "atom Cu",(al(i,1,1,nzz),i=1,3)
c       endif
c     do jt=1,nn(1,nzz)
c       j=nb(jt,1,nzz)
c       jj=nbb(jt,1,nzz)
c       if(nt(jj,j,nzz).eq.1) then
c       write(*,*) "atom W", (al(i,jj,j,nzz),i=1,3)
c       else
c       write(*,*) "atom Cu",(al(i,jj,j,nzz),i=1,3)
c       endif
c     enddo
c     endif
c
      end
c **********************************************************************

c **********************************************************************
c *      Use data to find the hamiltonian matrix elements              *
c **********************************************************************
      subroutine hmatrix(guess,nzz)
      implicit double precision (a-h,o-z)
      parameter(mnum=100,neinum=150)
      parameter(mnp=30,no1=9,no2=1,no3=4)
C WW, HH, WH
      parameter(ngu_WW=77,ngu_HH=6,ngu_WH=19)
c NN, WN, NH
      parameter(ngu_NN=18,ngu_NH=9,ngu_WN=67,ngu_tot=86)
c
      parameter(mnkp=186)
      parameter(nz=41)
      common/con1/vk(3,mnkp,nz),nkp(nz)
      common/con2/alat(3,3,nz),ralat(3,3,nz),tou(mnp,3,nz)
     &      ,itou(mnp,nz),acon(nz),mmm(nz),npa1(nz),npa2(nz),npa3(nz),mm
      common/con3/pi,r_cut,rc_W,rc_H,rc_N,rc_WH,rc_NH,rc_WN,efe
c WW, HH:
      common/data1/guess_WW(ngu_WW),guess_HH(ngu_HH)
c WH:
      common/data2/guess_WH(ngu_WH)
c NN, NH:
      common/data3/guess_NN(ngu_NN),guess_NH(ngu_NH)
c WN
c     common/data4/guess_WN(ngu_WN)
c
      common/location/al(3,mnp,mnum,nz),nt(mnp,mnum,nz)
      common/nei1/dxij(mnp,mnum,mnp,nz),dyij(mnp,mnum,mnp,nz)
     &           ,dzij(mnp,mnum,mnp,nz),rij0(mnp,mnum,mnp,nz)
     &           ,nn(mnp,nz),nb(neinum,mnp,nz),nbb(neinum,mnp,nz)
      common/nei2/rjk0(neinum,neinum,mnp,nz),nnc(neinum,mnp,nz)
     &           ,nbc(neinum,neinum,mnp,nz),nbbc(neinum,neinum,mnp,nz)
      common/hmat/hmatr(180,180,2,mnkp)
      common/eig/eigen(180,mnkp,nz)
      common/dewe/dee(84,mnkp,nz),wee(84,mnkp,nz)
      common/ndewe/ndee(nz)
      dimension guess(ngu_tot),onsite_1(3)
c
         ED_W=guess_WW(75)
         ES_W=guess_WW(75)+guess_WW(76)
         EP_W=guess_WW(75)+guess_WW(77)
c
c        dvalue=guess_WH(19)
c        ES_H=guess_HH(6)+dvalue
c
c      modified here for adjust the ES_N and EP_N:
c        guess_NH(9)=guess(86)
c         
         ES_H=guess_HH(6)+guess_WH(19)
c
         shift=guess_NH(9)-guess_WH(19)
c
         ES_N=guess(84)-shift          !!! Ns=Ns0+shift 
         EP_N=guess(84)+guess(85)-shift         !!! Np=Np0+shift 
c        ES_N=guess_NN(17)-shift          !!! Ns=Ns0+shift 
c        EP_N=guess_NN(17)+guess_NN(18)-shift         !!! Np=Np0+shift 

c        EP_N=guess_NN(21)+guess(29)-shift         !!! Np=Np0+shift 
c
c        write(19,*)',ED_W,ES_W,EP_W ',ED_W,ES_W,EP_W
c        write(19,*)',Es_N,Ep_N, Es_H',ES_N,EP_N,Es_H
c        stop
c
        nnpam=no1*npa1(nzz)+no2*npa2(nzz)+no3*npa3(nzz)
        nnpa=npa1(nzz)+npa2(nzz)+npa3(nzz)
c      DO L=1,mnkp
       DO L=1,NKP(NZZ)
        DO K=1,2
          DO J=1,nnpam
c         DO J=1,no1*mnp
           DO I=1,nnpam
c          DO I=1,no1*mnp
            HMATR(I,J,K,L)=0.0D0
           ENDDO
          ENDDO
        ENDDO
       ENDDO
c
      DO L=1,NKP(NZZ)
        DO II=1,nnpa
         if(itou(ii,nzz).eq.1) then  !!! 1 for w, and 2 for H .
c          if(nzz.eq.8)write(19,*)'pass W !'
         es=es_W
         ep=ep_W
         ed=ed_W
c
         IL=(II-1)*no1+1
         hmatr(IL,IL,1,L)=es
         do j=1,3
         IL=IL+1
         hmatr(IL,IL,1,L)=ep
         enddo
         do j=1,5
         IL=IL+1
         hmatr(IL,IL,1,L)=ed
         enddo
         endif
c
         if(itou(ii,nzz).eq.2) then   !!! 1 for w, and 2 for H .
         es=es_H
         IL=npa1(nzz)*no1+(II-npa1(nzz))
         hmatr(IL,IL,1,L)=es
c         if(nzz.eq.8)write(19,*)'pass H !',L,IL,IL,hmatr(IL,IL,1,L,nzz)
         endif
c
         if(itou(ii,nzz).eq.3) then   !!! 1 for w, and 2 for H, 3 for N.
         epij=0.0d0
         es=es_N
         delta_eP=0.0d0
c
c        if(npa1(nzz).eq.0.or.npa2(nzz).eq.0)goto 1212
c        if(nzz.eq.8.or.nzz.eq.9.or.nzz.eq.10)goto 1212
c        do jt=1,nn(ii,nzz)
c          j=nb(jt,ii,nzz)
c          jj=nbb(jt,ii,nzz)
c BCPan:
c          ntjj=nt(jj,j,nzz)
c          rij=rij0(jj,j,ii,nzz)
c          if(ntjj.eq.1.or.ntjj.eq.3)then
c          epij=epij+0.005d0/(1.0d0+
c    &          dexp(-onsite_1(ntjj)*rij-onsite_2))
c          endif
c       enddo
c212     continue
c
c          
c        delta_eP=onsite_3*epij+onsite_4*epij**2+onsite_5*epij**3
c        ep=ep_N+delta_eP          
c        es=es_N+delta_eP  
          ep=ep_N
          es=es_N

c          write(19,1)nzz,ii,epij                
1         format(2I5,2f12.5)
          IL=npa1(nzz)*no1+npa2(nzz)*no2+4*(II-npa1(nzz)-npa2(nzz)-1)+1
         hmatr(IL,IL,1,L)=es
c         if(nzz.eq.8)write(19,*)'pass N !',L,IL,IL,hmatr(IL,IL,1,L,nzz)
         do j=1,3
         IL=IL+1
         hmatr(IL,IL,1,L)=ep
c         if(nzz.eq.8)write(19,*)'pass N !',L,IL,IL,hmatr(IL,IL,1,L,nzz)
         enddo
         endif

        ENDDO
      ENDDO
c BCPan
c
      if(npa1(nzz).gt.0) call c_hmat11(guess_WW,dcu,nzz) ! W-W
      if(npa2(nzz).gt.0) call c_hmat22(guess_HH,dw,nzz)  ! H-H
c     if(npa3(nzz).gt.0) call c_hmat33(guess_NN,nzz)     ! N-N 
      if(npa3(nzz).gt.0) call c_hmat33(guess,nzz)     ! N-N 
      if(npa1(nzz).gt.0.and.npa2(nzz).gt.0) call c_hmat12(guess_WH,nzz)
      if(npa2(nzz).gt.0.and.npa3(nzz).gt.0) call c_hmat23(guess_NH,nzz)
      if(npa1(nzz).gt.0.and.npa3(nzz).gt.0) call c_hmat13(guess,nzz) 
1213   format(23f6.2)
c
      end
************************************************************************

c **********************************************************************
c *      Use data to find the hamiltonian matrix elements              *
c **********************************************************************
      subroutine c_hmat11(guess,dcu,nzz)
      implicit double precision (a-h,o-z)
      parameter(mnum=100,neinum=150)
      parameter(mnp=30,no1=9)
C WW, HH, WH
      parameter(ngu_WW=77,ngu_HH=6,ngu_WH=19)
c NN, WN, NH
      parameter(ngu_NN=18,ngu_NH=9,ngu_WN=67)
c
      parameter(mnkp=186)
      parameter(nz=41)
      common/con1/vk(3,mnkp,nz),nkp(nz)
      common/con2/alat(3,3,nz),ralat(3,3,nz),tou(mnp,3,nz)
     &      ,itou(mnp,nz),acon(nz),mmm(nz),npa1(nz),npa2(nz),npa3(nz),mm
      common/con3/pi,r_cut,rc_W,rc_H,rc_N,rc_WH,rc_NH,rc_WN,efe
c     common/data1/guess_WW(ngu_WW),guess_HH(ngu_HH)
c     common/data2/guess_WH(ngu_WH)
      common/location/al(3,mnp,mnum,nz),nt(mnp,mnum,nz)
      common/nei1/dxij(mnp,mnum,mnp,nz),dyij(mnp,mnum,mnp,nz)
     &           ,dzij(mnp,mnum,mnp,nz),rij0(mnp,mnum,mnp,nz)
     &           ,nn(mnp,nz),nb(neinum,mnp,nz),nbb(neinum,mnp,nz)
      common/nei2/rjk0(neinum,neinum,mnp,nz),nnc(neinum,mnp,nz)
     &           ,nbc(neinum,neinum,mnp,nz),nbbc(neinum,neinum,mnp,nz)
      common/hmat/hmatr(180,180,2,mnkp)
      common/eig/eigen(180,mnkp,nz)
c     common/hmat/hmatr(no1*mnp,no1*mnp,2,mnkp,nz)
c     common/eig/eigen(no1*mnp,mnkp,nz)
      common/dewe/dee(84,mnkp,nz),wee(84,mnkp,nz)
      common/ndewe/ndee(nz)
c
      dimension guess(ngu_WW)
      dimension vbi(10)
      dimension arf1(10),arf2(10),arf3(10),arf4(10)
     &         ,arf6(10),arf7(10),sor_W(10)
      dimension ehmat(no1,no1)
c     dimension bta1(2)
c
         arf5=guess(71)
c        bta1(1)=guess(72)
         bta1=guess(72)
         bta2=guess(73)
         bta3=guess(74)
         sq3=dsqrt(3.0d0)
c
      do m=1,10
         n=(m-1)*7
         arf1(m)=guess(n+1)
         arf2(m)=guess(n+2)
         arf3(m)=guess(n+3)
         arf4(m)=guess(n+4)
         arf6(m)=guess(n+5)
         arf7(m)=guess(n+6)
         sor_W(m)=guess(n+7)
      enddo
c
c        bta1(2)=dcu
        nnpam=no1*npa1(nzz)+no2*npa2(nzz)+no3*npa3(nzz)
        nnpa=npa1(nzz)+npa2(nzz)+npa3(nzz)
c     For W-W Block.
      do ii=1,nnpa
          ntii=nt(ii,1,nzz)
        if(ntii.eq.1) then  ! 1 for W.
        do jt=1,nn(ii,nzz)
           j=nb(jt,ii,nzz)
           jj=nbb(jt,ii,nzz)
c BCPan:
           ntjj=nt(jj,j,nzz)
c
          if(ntjj.eq.1) then ! 1 for W
           rij=rij0(jj,j,ii,nzz)
c          write(19,*)'W-W ',nzz,ii,jj,rij
c
           cosl=dxij(jj,j,ii,nzz)/rij
           cosm=dyij(jj,j,ii,nzz)/rij
           cosn=dzij(jj,j,ii,nzz)/rij
c
           esij=0.0d0
c
          do kt=1,nnc(jt,ii,nzz)
           k=nbc(kt,jt,ii,nzz) !Count in common neighbors
           kk=nbbc(kt,jt,ii,nzz)
           ntkk=nt(kk,k,nzz)
c
           if(ntkk.eq.1)then
           rjk=rjk0(kt,jt,ii,nzz)
           rik=rij0(kk,k,ii,nzz)
c          esij=esij+rij/(rik+rjk)*bta1(ntkk)
           esij=esij+rij/(rik+rjk)*bta1
     *               *dexp(-bta2*((rik+rjk)/rij)**bta3)
           endif
           enddo
           delt=2.0d0/pi*datan(esij)
          do L=1,10
           rrij=rij*(1.0d0+arf5*delt+arf6(L)*delt**2+arf7(L)*delt**3)
           vbi(L)=arf1(L)*rrij**arf2(L)*dexp(-arf3(L)*rrij**arf4(L))
     &                   /(1+dexp((rrij-rc_W)/sor_W(L)))
          enddo
c
            ehmat(1,1)=vbi(1)
            ehmat(1,2)=cosl*vbi(2)
            ehmat(1,3)=cosm*vbi(2)
            ehmat(1,4)=cosn*vbi(2)
            ehmat(1,5)=sq3*cosl*cosm*vbi(5)
            ehmat(1,6)=sq3*cosm*cosn*vbi(5)
            ehmat(1,7)=sq3*cosl*cosn*vbi(5)
            ehmat(1,8)=0.5d0*sq3*(cosl**2-cosm**2)*vbi(5)
            ehmat(1,9)=(cosn**2-0.5d0*(cosl**2+cosm**2))*vbi(5)
            ehmat(2,2)=cosl**2*vbi(3)+(1.0d0-cosl**2)*vbi(4)
            ehmat(2,3)=cosl*cosm*(vbi(3)-vbi(4))
            ehmat(2,4)=cosl*cosn*(vbi(3)-vbi(4))
            ehmat(2,5)=sq3*cosl**2*cosm*vbi(6)
     &                 +cosm*(1.0d0-2.0d0*cosl**2)*vbi(7)
            ehmat(2,6)=sq3*cosl*cosm*cosn*vbi(6)
     &                 -2.0d0*cosl*cosm*cosn*vbi(7)
            ehmat(2,7)=sq3*cosl**2*cosn*vbi(6)
     &                 +cosn*(1.0d0-2*cosl**2)*vbi(7)
            ehmat(2,8)=0.5d0*sq3*cosl*(cosl**2-cosm**2)*vbi(6)
     &                 +cosl*(1.0d0-cosl**2+cosm**2)*vbi(7)
            ehmat(2,9)=cosl*(cosn**2-0.5d0*(cosl**2+cosm**2))*vbi(6)
     &                 -sq3*cosl*cosn**2*vbi(7)
            ehmat(3,3)=cosm**2*vbi(3)+(1.0d0-cosm**2)*vbi(4)
            ehmat(3,4)=cosm*cosn*(vbi(3)-vbi(4))
            ehmat(3,5)=sq3*cosl*cosm**2*vbi(6)
     &                 +cosl*(1.0d0-2*cosm**2)*vbi(7)
            ehmat(3,6)=sq3*cosm**2*cosn*vbi(6)
     &                 +cosn*(1.0d0-2*cosm**2)*vbi(7)
            ehmat(3,7)=sq3*cosl*cosm*cosn*vbi(6)
     &                 -2.0d0*cosl*cosm*cosn*vbi(7)
            ehmat(3,8)=0.5d0*sq3*cosm*(cosl**2-cosm**2)*vbi(6)
     &                 -cosm*(1.0d0+cosl**2-cosm**2)*vbi(7)
            ehmat(3,9)=cosm*(cosn**2-0.5d0*(cosl**2+cosm**2))*vbi(6)
     &                 -sq3*cosm*cosn**2*vbi(7)
            ehmat(4,4)=cosn**2*vbi(3)+(1.0d0-cosn**2)*vbi(4)
            ehmat(4,5)=sq3*cosl*cosm*cosn*vbi(6)
     &                 -2.0d0*cosl*cosm*cosn*vbi(7)
            ehmat(4,6)=sq3*cosm*cosn**2*vbi(6)
     &                 +cosm*(1.0d0-2.0d0*cosn**2)*vbi(7)
            ehmat(4,7)=sq3*cosl*cosn**2*vbi(6)
     &                 +cosl*(1.0d0-2.0d0*cosn**2)*vbi(7)
            ehmat(4,8)=0.5d0*sq3*cosn*(cosl**2-cosm**2)*vbi(6)
     &                 -cosn*(cosl**2-cosm**2)*vbi(7)
            ehmat(4,9)=cosn*(cosn**2-0.5d0*(cosl**2+cosm**2))*vbi(6)
     &                 +sq3*cosn*(cosl**2+cosm**2)*vbi(7)
            ehmat(5,5)=3.0d0*cosl**2*cosm**2*vbi(8)
     &                 +(cosl**2+cosm**2-4.0d0*cosl**2*cosm**2)*vbi(9)
     &                 +(cosn**2+cosl**2*cosm**2)*vbi(10)
            ehmat(5,6)=3.0d0*cosl*cosm**2*cosn*vbi(8)
     &                 +cosl*cosn*(1.0d0-4.0d0*cosm**2)*vbi(9)
     &                 +cosl*cosn*(cosm**2-1.0d0)*vbi(10)
            ehmat(5,7)=3.0d0*cosl**2*cosm*cosn*vbi(8)
     &                 +cosm*cosn*(1.0d0-4.0d0*cosl**2)*vbi(9)
     &                 +cosm*cosn*(cosl**2-1.0d0)*vbi(10)
            ehmat(5,8)=1.5d0*cosl*cosm*(cosl**2-cosm**2)*vbi(8)
     &                 +2.0d0*cosl*cosm*(cosm**2-cosl**2)*vbi(9)
     &                 +0.5d0*cosl*cosm*(cosl**2-cosm**2)*vbi(10)
            ehmat(5,9)=sq3*cosl*cosm*(cosn**2-0.5d0*(cosl**2+cosm**2))
     &                    *vbi(8)-2.0d0*sq3*cosl*cosm*cosn**2*vbi(9)
     &                 +0.5d0*sq3*cosl*cosm*(1.0d0+cosn**2)*vbi(10)
            ehmat(6,6)=3.0d0*cosm**2*cosn**2*vbi(8)
     &                 +(cosm**2+cosn**2-4.0d0*cosm**2*cosn**2)*vbi(9)
     &                 +(cosl**2+cosm**2*cosn**2)*vbi(10)
            ehmat(6,7)=3.0d0*cosl*cosm*cosn**2*vbi(8)
     &                 +cosl*cosm*(1.0d0-4.0d0*cosn**2)*vbi(9)
     &                 +cosl*cosm*(cosn**2-1.0d0)*vbi(10)
            ehmat(6,8)=1.5d0*cosm*cosn*(cosl**2-cosm**2)*vbi(8)
     &                -cosm*cosn*(1.0d0+2.0d0*(cosl**2-cosm**2))*vbi(9)
     &                +cosm*cosn*(1.0d0+0.5d0*(cosl**2-cosm**2))*vbi(10)
            ehmat(6,9)=sq3*cosm*cosn*(cosn**2-0.5d0*(cosl**2+cosm**2))
     &                    *vbi(8)
     &                 +sq3*cosm*cosn*(cosl**2+cosm**2-cosn**2)*vbi(9)
     &                 -0.5d0*sq3*cosm*cosn*(cosl**2+cosm**2)*vbi(10)
            ehmat(7,7)=3.0d0*cosl**2*cosn**2*vbi(8)
     &                 +(cosl**2+cosn**2-4.0d0*cosl**2*cosn**2)*vbi(9)
     &                 +(cosm**2+cosl**2*cosn**2)*vbi(10)
            ehmat(7,8)=1.5d0*cosl*cosn*(cosl**2-cosm**2)*vbi(8)
     &                +cosl*cosn*(1.0d0-2.0d0*(cosl**2-cosm**2))*vbi(9)
     &                -cosl*cosn*(1.0d0-0.5d0*(cosl**2-cosm**2))*vbi(10)
            ehmat(7,9)=sq3*cosl*cosn*(cosn**2-0.5d0*(cosl**2+cosm**2))
     &                    *vbi(8)
     &                 +sq3*cosl*cosn*(cosl**2+cosm**2-cosn**2)*vbi(9)
     &                 -0.5d0*sq3*cosl*cosn*(cosl**2+cosm**2)*vbi(10)
            ehmat(8,8)=0.75d0*(cosl**2-cosm**2)**2*vbi(8)
     &                 +(cosl**2+cosm**2-(cosl**2-cosm**2)**2)*vbi(9)
     &                 +(cosn**2+0.25d0*(cosl**2-cosm**2)**2)*vbi(10)
            ehmat(8,9)=0.5d0*sq3*(cosl**2-cosm**2)
     &                       *(cosn**2-0.5d0*(cosl**2+cosm**2))*vbi(8)
     &              +sq3*cosn**2*(cosm**2-cosl**2)*vbi(9)
     &              +sq3/4.0d0*(1.0d0+cosn**2)*(cosl**2-cosm**2)*vbi(10)
            ehmat(9,9)=(cosn**2-0.5d0*(cosl**2+cosm**2))**2*vbi(8)
     &                  +3.0d0*cosn**2*(cosl**2+cosm**2)*vbi(9)
     &                  +0.75d0*(cosl**2+cosm**2)**2*vbi(10)
c
          do ia=1,no1-1
           do ib=ia+1,no1
            if( (ia.eq.1.and.(ib.ge.2.and.ib.le.4)).or.
     &        ((ib.ge.5.and.ib.le.9).and.(ia.ge.2.and.ia.le.4)) ) then
              ehmat(ib,ia)=-ehmat(ia,ib)
            else
              ehmat(ib,ia)=ehmat(ia,ib)
            endif
           enddo
          enddo
c
               dx=dxij(jj,j,ii,nzz)
               dy=dyij(jj,j,ii,nzz)
               dz=dzij(jj,j,ii,nzz)
            do mk=1,nkp(nzz)
               v1=vk(1,mk,nzz)
               v2=vk(2,mk,nzz)
               v3=vk(3,mk,nzz)
              x1=+v1*ralat(1,1,nzz)+v2*ralat(2,1,nzz)
     &           +v3*ralat(3,1,nzz)
              y1=+v1*ralat(1,2,nzz)+v2*ralat(2,2,nzz)
     &           +v3*ralat(3,2,nzz)
              z1=+v1*ralat(1,3,nzz)+v2*ralat(2,3,nzz)
     &           +v3*ralat(3,3,nzz)
              pha=(x1*dx+y1*dy+z1*dz)
              prel=dcos(pha)
              pimg=dsin(pha)
c
             do iit=(ii-1)*no1+1,ii*no1
              do jjt=(jj-1)*no1+1,jj*no1
                 kkt=mod(iit-1,no1)+1
                 mmt=mod(jjt-1,no1)+1
c BCPan:
c               write(19,*)'W-W ',ii,jj,iit,jjt,kkt,mmt

                 hmatr(iit,jjt,1,mk)=hmatr(iit,jjt,1,mk)
     &                                  +ehmat(kkt,mmt)*prel
                 hmatr(iit,jjt,2,mk)=hmatr(iit,jjt,2,mk)
     &                                  +ehmat(kkt,mmt)*pimg
              enddo
             enddo
c
            enddo  ! end K recycle
            endif 
        enddo
        endif  ! end of W-W
      enddo
c
      end
c **********************************************************************

c **********************************************************************
c *      Use data to find the hamiltonian matrix elements              *
c **********************************************************************
      subroutine c_hmat22(guess,dw,nzz)
      implicit double precision (a-h,o-z)
      parameter(mnum=100,neinum=150)
      parameter(mnp=30,no1=9,no2=1,no3=4)
C WW, HH, WH
      parameter(ngu_WW=77,ngu_HH=6,ngu_WH=19)
c NN, WN, NH
      parameter(ngu_NN=18,ngu_NH=9,ngu_WN=67)
c
      parameter(mnkp=186)
      parameter(nz=41)
      common/con1/vk(3,mnkp,nz),nkp(nz)
      common/con2/alat(3,3,nz),ralat(3,3,nz),tou(mnp,3,nz)
     &      ,itou(mnp,nz),acon(nz),mmm(nz),npa1(nz),npa2(nz),npa3(nz),mm
      common/con3/pi,r_cut,rc_W,rc_H,rc_N,rc_WH,rc_NH,rc_WN,efe
c     common/data1/guess_WW(ngu_WW),guess_HH(ngu_HH)
c     common/data2/guess_WH(ngu_WH)
      common/location/al(3,mnp,mnum,nz),nt(mnp,mnum,nz)
      common/nei1/dxij(mnp,mnum,mnp,nz),dyij(mnp,mnum,mnp,nz)
     &           ,dzij(mnp,mnum,mnp,nz),rij0(mnp,mnum,mnp,nz)
     &           ,nn(mnp,nz),nb(neinum,mnp,nz),nbb(neinum,mnp,nz)
      common/nei2/rjk0(neinum,neinum,mnp,nz),nnc(neinum,mnp,nz)
     &           ,nbc(neinum,neinum,mnp,nz),nbbc(neinum,neinum,mnp,nz)
      common/hmat/hmatr(180,180,2,mnkp)
      common/eig/eigen(180,mnkp,nz)
c     common/hmat/hmatr(no1*mnp,no1*mnp,2,mnkp,nz)
c     common/eig/eigen(no1*mnp,mnkp,nz)
      common/dewe/dee(84,mnkp,nz),wee(84,mnkp,nz)
      common/ndewe/ndee(nz)
c
      dimension guess(ngu_HH)
      dimension vbi(10)
      dimension arf1(10),arf2(10),arf3(10),arf4(10)
     &         ,arf6(10),arf7(10),sor_H(10)
      dimension ehmat(no2,no2)
      dimension bta1(2)
c
c        arf5=guess(71)
c        bta1(2)=guess(72)
c        bta2=guess(73)
c        bta3=guess(74)
         sq3=dsqrt(3.0d0)
c      do m=1,10
       do m=1,1 
         n=(m-1)*7
         arf1(m)=guess(n+1)
         arf2(m)=guess(n+2)
         arf3(m)=guess(n+3)
         arf4(m)=guess(n+4)
         sor_H(m)=guess(n+5)
       enddo
c
c     For  H-H  Block.
        nnpam=no1*npa1(nzz)+no2*npa2(nzz)+no3*npa3(nzz)
        nnpa=npa1(nzz)+npa2(nzz)+npa3(nzz)
      do ii=1,nnpa
          ntii=nt(ii,1,nzz)
        if(ntii.eq.2) then  ! 2 for H .
        do jt=1,nn(ii,nzz)
           j=nb(jt,ii,nzz)
           jj=nbb(jt,ii,nzz)
           ntjj=nt(jj,j,nzz)
c
          if(ntjj.eq.2) then ! 2 for H .
           rij=rij0(jj,j,ii,nzz)
c
           cosl=dxij(jj,j,ii,nzz)/rij
           cosm=dyij(jj,j,ii,nzz)/rij
           cosn=dzij(jj,j,ii,nzz)/rij
c
           esij=0.0d0
c
c         do kt=1,nnc(jt,ii,nzz)
c           k=nbc(kt,jt,ii,nzz) !Count incommon neighbors.
c           kk=nbbc(kt,jt,ii,nzz)
c           ntkk=nt(kk,k,nzz)
c
c           rjk=rjk0(kt,jt,ii,nzz)
c           rik=rij0(kk,k,ii,nzz)
c           esij=esij+rij/(rik+rjk)*bta1(ntkk)
c    *                *dexp(-bta2*((rik+rjk)/rij)**bta3)
c         enddo
c           delt=2.0d0/pi*datan(esij)
c         do L=1,10
          do L=1,1
c          rrij=rij*(1.0d0+arf5*delt+arf6(L)*delt**2+arf7(L)*delt**3)
           rrij=rij
           vbi(L)=arf1(L)*rrij**arf2(L)*dexp(-arf3(L)*rrij**arf4(L))
     &                   /(1+dexp((rrij-rc_H)/sor_H(L)))
          enddo
c
            ehmat(1,1)=vbi(1)
c BCPan
c           if(nzz.eq.8)then
c          write(19,*)'test H-H',arf1(1),arf2(1),arf3(1),
c    1         arf4(1),rc_H,sor_H(1)
c          write(19,*)'ehmat ',ii,jj,rrij,vbi(1)
c           endif
c
c           ehmat(1,2)=cosl*vbi(2)
c           ehmat(1,3)=cosm*vbi(2)
c           ehmat(1,4)=cosn*vbi(2)
c           ehmat(1,5)=sq3*cosl*cosm*vbi(5)
c           ehmat(1,6)=sq3*cosm*cosn*vbi(5)
c           ehmat(1,7)=sq3*cosl*cosn*vbi(5)
c           ehmat(1,8)=0.5d0*sq3*(cosl**2-cosm**2)*vbi(5)
c           ehmat(1,9)=(cosn**2-0.5d0*(cosl**2+cosm**2))*vbi(5)
c           ehmat(2,2)=cosl**2*vbi(3)+(1.0d0-cosl**2)*vbi(4)
c           ehmat(2,3)=cosl*cosm*(vbi(3)-vbi(4))
c           ehmat(2,4)=cosl*cosn*(vbi(3)-vbi(4))
c           ehmat(2,5)=sq3*cosl**2*cosm*vbi(6)
c    &                 +cosm*(1.0d0-2.0d0*cosl**2)*vbi(7)
c           ehmat(2,6)=sq3*cosl*cosm*cosn*vbi(6)
c    &                 -2.0d0*cosl*cosm*cosn*vbi(7)
c           ehmat(2,7)=sq3*cosl**2*cosn*vbi(6)
c    &                 +cosn*(1.0d0-2*cosl**2)*vbi(7)
c           ehmat(2,8)=0.5d0*sq3*cosl*(cosl**2-cosm**2)*vbi(6)
c    &                 +cosl*(1.0d0-cosl**2+cosm**2)*vbi(7)
c           ehmat(2,9)=cosl*(cosn**2-0.5d0*(cosl**2+cosm**2))*vbi(6)
c    &                 -sq3*cosl*cosn**2*vbi(7)
c           ehmat(3,3)=cosm**2*vbi(3)+(1.0d0-cosm**2)*vbi(4)
c           ehmat(3,4)=cosm*cosn*(vbi(3)-vbi(4))
c           ehmat(3,5)=sq3*cosl*cosm**2*vbi(6)
c    &                 +cosl*(1.0d0-2*cosm**2)*vbi(7)
c           ehmat(3,6)=sq3*cosm**2*cosn*vbi(6)
c    &                 +cosn*(1.0d0-2*cosm**2)*vbi(7)
c           ehmat(3,7)=sq3*cosl*cosm*cosn*vbi(6)
c    &                 -2.0d0*cosl*cosm*cosn*vbi(7)
c           ehmat(3,8)=0.5d0*sq3*cosm*(cosl**2-cosm**2)*vbi(6)
c    &                 -cosm*(1.0d0+cosl**2-cosm**2)*vbi(7)
c           ehmat(3,9)=cosm*(cosn**2-0.5d0*(cosl**2+cosm**2))*vbi(6)
c    &                 -sq3*cosm*cosn**2*vbi(7)
c           ehmat(4,4)=cosn**2*vbi(3)+(1.0d0-cosn**2)*vbi(4)
c           ehmat(4,5)=sq3*cosl*cosm*cosn*vbi(6)
c    &                 -2.0d0*cosl*cosm*cosn*vbi(7)
c           ehmat(4,6)=sq3*cosm*cosn**2*vbi(6)
c    &                 +cosm*(1.0d0-2.0d0*cosn**2)*vbi(7)
c           ehmat(4,7)=sq3*cosl*cosn**2*vbi(6)
c    &                 +cosl*(1.0d0-2.0d0*cosn**2)*vbi(7)
c           ehmat(4,8)=0.5d0*sq3*cosn*(cosl**2-cosm**2)*vbi(6)
c    &                 -cosn*(cosl**2-cosm**2)*vbi(7)
c           ehmat(4,9)=cosn*(cosn**2-0.5d0*(cosl**2+cosm**2))*vbi(6)
c    &                 +sq3*cosn*(cosl**2+cosm**2)*vbi(7)
c           ehmat(5,5)=3.0d0*cosl**2*cosm**2*vbi(8)
c    &                 +(cosl**2+cosm**2-4.0d0*cosl**2*cosm**2)*vbi(9)
c    &                 +(cosn**2+cosl**2*cosm**2)*vbi(10)
c           ehmat(5,6)=3.0d0*cosl*cosm**2*cosn*vbi(8)
c    &                 +cosl*cosn*(1.0d0-4.0d0*cosm**2)*vbi(9)
c    &                 +cosl*cosn*(cosm**2-1.0d0)*vbi(10)
c           ehmat(5,7)=3.0d0*cosl**2*cosm*cosn*vbi(8)
c    &                 +cosm*cosn*(1.0d0-4.0d0*cosl**2)*vbi(9)
c    &                 +cosm*cosn*(cosl**2-1.0d0)*vbi(10)
c           ehmat(5,8)=1.5d0*cosl*cosm*(cosl**2-cosm**2)*vbi(8)
c    &                 +2.0d0*cosl*cosm*(cosm**2-cosl**2)*vbi(9)
c    &                 +0.5d0*cosl*cosm*(cosl**2-cosm**2)*vbi(10)
c           ehmat(5,9)=sq3*cosl*cosm*(cosn**2-0.5d0*(cosl**2+cosm**2))
c    &                    *vbi(8)-2.0d0*sq3*cosl*cosm*cosn**2*vbi(9)
c    &                 +0.5d0*sq3*cosl*cosm*(1.0d0+cosn**2)*vbi(10)
c           ehmat(6,6)=3.0d0*cosm**2*cosn**2*vbi(8)
c    &                 +(cosm**2+cosn**2-4.0d0*cosm**2*cosn**2)*vbi(9)
c    &                 +(cosl**2+cosm**2*cosn**2)*vbi(10)
c           ehmat(6,7)=3.0d0*cosl*cosm*cosn**2*vbi(8)
c    &                 +cosl*cosm*(1.0d0-4.0d0*cosn**2)*vbi(9)
c    &                 +cosl*cosm*(cosn**2-1.0d0)*vbi(10)
c           ehmat(6,8)=1.5d0*cosm*cosn*(cosl**2-cosm**2)*vbi(8)
c    &                -cosm*cosn*(1.0d0+2.0d0*(cosl**2-cosm**2))*vbi(9)
c    &                +cosm*cosn*(1.0d0+0.5d0*(cosl**2-cosm**2))*vbi(10)
c           ehmat(6,9)=sq3*cosm*cosn*(cosn**2-0.5d0*(cosl**2+cosm**2))
c    &                    *vbi(8)
c    &                 +sq3*cosm*cosn*(cosl**2+cosm**2-cosn**2)*vbi(9)
c    &                 -0.5d0*sq3*cosm*cosn*(cosl**2+cosm**2)*vbi(10)
c           ehmat(7,7)=3.0d0*cosl**2*cosn**2*vbi(8)
c    &                 +(cosl**2+cosn**2-4.0d0*cosl**2*cosn**2)*vbi(9)
c    &                 +(cosm**2+cosl**2*cosn**2)*vbi(10)
c           ehmat(7,8)=1.5d0*cosl*cosn*(cosl**2-cosm**2)*vbi(8)
c    &                +cosl*cosn*(1.0d0-2.0d0*(cosl**2-cosm**2))*vbi(9)
c    &                -cosl*cosn*(1.0d0-0.5d0*(cosl**2-cosm**2))*vbi(10)
c           ehmat(7,9)=sq3*cosl*cosn*(cosn**2-0.5d0*(cosl**2+cosm**2))
c    &                    *vbi(8)
c    &                 +sq3*cosl*cosn*(cosl**2+cosm**2-cosn**2)*vbi(9)
c    &                 -0.5d0*sq3*cosl*cosn*(cosl**2+cosm**2)*vbi(10)
c           ehmat(8,8)=0.75d0*(cosl**2-cosm**2)**2*vbi(8)
c    &                 +(cosl**2+cosm**2-(cosl**2-cosm**2)**2)*vbi(9)
c    &                 +(cosn**2+0.25d0*(cosl**2-cosm**2)**2)*vbi(10)
c           ehmat(8,9)=0.5d0*sq3*(cosl**2-cosm**2)
c    &                       *(cosn**2-0.5d0*(cosl**2+cosm**2))*vbi(8)
c    &              +sq3*cosn**2*(cosm**2-cosl**2)*vbi(9)
c    &              +sq3/4.0d0*(1.0d0+cosn**2)*(cosl**2-cosm**2)*vbi(10)
c           ehmat(9,9)=(cosn**2-0.5d0*(cosl**2+cosm**2))**2*vbi(8)
c    &                  +3.0d0*cosn**2*(cosl**2+cosm**2)*vbi(9)
c    &                  +0.75d0*(cosl**2+cosm**2)**2*vbi(10)
c BCPan
c         do ia=1,1
c           write(19,19)(ehmat(ia,ib),ib=1,1)
c19          format(9f12.5)
c         enddo
c
c          do ia=1,no2
c           do ib=ia+1,no2
c            if( (ia.eq.1.and.(ib.ge.2.and.ib.le.4)).or.
c     &        ((ib.ge.5.and.ib.le.9).and.(ia.ge.2.and.ia.le.4)) ) then
c              ehmat(ib,ia)=-ehmat(ia,ib)
c            else
c             ehmat(ib,ia)=ehmat(ia,ib)
c            endif
c           enddo
c          enddo
c   ??????????  BCPan
c
               dx=dxij(jj,j,ii,nzz)
               dy=dyij(jj,j,ii,nzz)
               dz=dzij(jj,j,ii,nzz)
            do mk=1,nkp(nzz)
               v1=vk(1,mk,nzz)
               v2=vk(2,mk,nzz)
               v3=vk(3,mk,nzz)
              x1=+v1*ralat(1,1,nzz)+v2*ralat(2,1,nzz)
     &           +v3*ralat(3,1,nzz)
              y1=+v1*ralat(1,2,nzz)+v2*ralat(2,2,nzz)
     &           +v3*ralat(3,2,nzz)
              z1=+v1*ralat(1,3,nzz)+v2*ralat(2,3,nzz)
     &           +v3*ralat(3,3,nzz)
              pha=(x1*dx+y1*dy+z1*dz)
              prel=dcos(pha)
              pimg=dsin(pha)
c
             ITH=npa1(nzz)*no1
             do iit=ITH+(ii-npa1(nzz))*no2,ITH+(ii-npa1(nzz))*no2
               do jjt=ITH+(jj-npa1(nzz))*no2,ITH+(jj-npa1(nzz))*no2
c               mmt=mod(jjt-1,no2)+1
c               kkt=mod(iit-1,no2)+1
c               write(19,*)'H-H',ii,jj,iit,jjt
               kkt=1
               mmt=1

               hmatr(iit,jjt,1,mk)=hmatr(iit,jjt,1,mk)
     &                                +ehmat(kkt,mmt)*prel
               hmatr(iit,jjt,2,mk)=hmatr(iit,jjt,2,mk)
     &                                +ehmat(kkt,mmt)*pimg
c BCPan
c               if(nzz.eq.1)write(19,*)hmatr(iit,jjt,1,mk,nzz)
               enddo
             enddo
c
            enddo  ! end K recycle
            endif 
        enddo
        endif  ! end of  H-H .
      enddo
c
      end
c **********************************************************************
c **********************************************************************
c *      Use data to find the hamiltonian matrix elements              *
c **********************************************************************
      subroutine c_hmat33(guess,nzz)
      implicit double precision (a-h,o-z)
      parameter(mnum=100,neinum=150)
      parameter(mnp=30,no1=9,no2=1,no3=4)
C WW, HH, WH
      parameter(ngu_WW=77,ngu_HH=6,ngu_WH=19)
c NN, WN, NH
      parameter(ngu_NN=18,ngu_NH=9,ngu_WN=67,ngu_tot=86)
c
      parameter(mnkp=186)
      parameter(nz=41)
      common/con1/vk(3,mnkp,nz),nkp(nz)
      common/con2/alat(3,3,nz),ralat(3,3,nz),tou(mnp,3,nz)
     &      ,itou(mnp,nz),acon(nz),mmm(nz),npa1(nz),npa2(nz),npa3(nz),mm
      common/con3/pi,r_cut,rc_W,rc_H,rc_N,rc_WH,rc_NH,rc_WN,efe
c     common/data1/guess_WW(ngu_WW),guess_HH(ngu_HH)
c     common/data2/guess_WH(ngu_WH)
      common/location/al(3,mnp,mnum,nz),nt(mnp,mnum,nz)
      common/nei1/dxij(mnp,mnum,mnp,nz),dyij(mnp,mnum,mnp,nz)
     &           ,dzij(mnp,mnum,mnp,nz),rij0(mnp,mnum,mnp,nz)
     &           ,nn(mnp,nz),nb(neinum,mnp,nz),nbb(neinum,mnp,nz)
      common/nei2/rjk0(neinum,neinum,mnp,nz),nnc(neinum,mnp,nz)
     &           ,nbc(neinum,neinum,mnp,nz),nbbc(neinum,neinum,mnp,nz)
      common/hmat/hmatr(180,180,2,mnkp)
      common/eig/eigen(180,mnkp,nz)
c     common/hmat/hmatr(no1*mnp,no1*mnp,2,mnkp,nz)
c     common/eig/eigen(no1*mnp,mnkp,nz)
      common/dewe/dee(84,mnkp,nz),wee(84,mnkp,nz)
      common/ndewe/ndee(nz)
c
      dimension guess(ngu_tot)
      dimension vbi(10)
      dimension arf1(10),arf2(10),arf3(10),arf4(10)
     &         ,arf6(10),arf7(10),sor_N(10)
      dimension ehmat(no3,no3)
c     dimension bta1(2)
c
c        arf5=guess(71)
c        bta1(1)=guess(72)
c        bta1=guess(72)
c        bta2=guess(73)
c        bta3=guess(74)
         sq3=dsqrt(3.0d0)
c
      do m=1,4 
         n=(m-1)*4+67
         arf1(m)=guess(n+1)
         arf2(m)=guess(n+2)
         arf3(m)=guess(n+3)
         arf4(m)=guess(n+4)
         sor_N(m)=0.05d0
c        write(19,*)guess(n+1),guess(n+2),guess(n+3),guess(n+4)
      enddo
c     do m=1,4 
c        n=(m-1)*4   
c        arf1(m)=guess(n+1)
c        arf2(m)=guess(n+2)
c        arf3(m)=guess(n+3)
c        arf4(m)=guess(n+4)
c        sor_N(m)=0.05d0
c     enddo
c
c        bta1(2)=dcu
        nnpam=no1*npa1(nzz)+no2*npa2(nzz)+no3*npa3(nzz)
        nnpa=npa1(nzz)+npa2(nzz)+npa3(nzz)
c     For N-N Block.
c 
      do ii=1,nnpa
          ntii=nt(ii,1,nzz)
        if(ntii.eq.3) then  ! 1 for N.
        do jt=1,nn(ii,nzz)
           j=nb(jt,ii,nzz)
           jj=nbb(jt,ii,nzz)
c BCPan:
           ntjj=nt(jj,j,nzz)
c
          if(ntjj.eq.3) then ! 1 for N
           rij=rij0(jj,j,ii,nzz)
c          write(*,*) rij
c
c           cosl=dxij(jj,j,ii,nzz)/rij
c           cosm=dyij(jj,j,ii,nzz)/rij
c           cosn=dzij(jj,j,ii,nzz)/rij
c
           CL=dxij(jj,j,ii,nzz)/rij
           CM=dyij(jj,j,ii,nzz)/rij
           CN=dzij(jj,j,ii,nzz)/rij
           CL2=CL*CL
           CM2=CM*CM
           CN2=CN*CN
c
           esij=0.0d0
c
c         do kt=1,nnc(jt,ii,nzz)
c          k=nbc(kt,jt,ii,nzz) !Count in common neighbors
c          kk=nbbc(kt,jt,ii,nzz)
c          ntkk=nt(kk,k,nzz)
c
c          rjk=rjk0(kt,jt,ii,nzz)
c          rik=rij0(kk,k,ii,nzz)
c          esij=esij+rij/(rik+rjk)*bta1(ntkk)
c          esij=esij+rij/(rik+rjk)*bta1
c    *               *dexp(-bta2*((rik+rjk)/rij)**bta3)
c         enddo
c          delt=2.0d0/pi*datan(esij)
          do L=1,4
c          rrij=rij*(1.0d0+arf5*delt+arf6(L)*delt**2+arf7(L)*delt**3)
           rrij=rij
           vbi(L)=arf1(L)*rrij**arf2(L)*dexp(-arf3(L)*rrij**arf4(L))
     &                   /(1+dexp((rrij-rc_N)/sor_N(L)))
          enddo
c
          SSS=VBI(1)
          SPS=VBI(2)
          PSS=VBI(2)
          PPS=VBI(3)
          PPP=VBI(4)
c          
           EHMAT(1,1)=+SSS
           EHMAT(1,2)=CL*SPS
           EHMAT(1,3)=CM*SPS
           EHMAT(1,4)=CN*SPS
           EHMAT(2,1)=-CL*PSS
           EHMAT(2,2)=+CL2*PPS+(1.0D0-CL2)*PPP
           EHMAT(2,3)=+CL*CM*(PPS-PPP)
           EHMAT(2,4)=+CL*CN*(PPS-PPP)
           EHMAT(3,1)=-CM*PSS
           EHMAT(3,2)=+EHMAT(2,3)
           EHMAT(3,3)=+CM2*PPS+(1.0D0-CM2)*PPP
           EHMAT(3,4)=+CM*CN*(PPS-PPP)
           EHMAT(4,1)=-CN*PSS
           EHMAT(4,2)=+EHMAT(2,4)
           EHMAT(4,3)=+EHMAT(3,4)
           EHMAT(4,4)=+CN**2*PPS+(1.0D0-CN2)*PPP

c           ehmat(1,1)=vbi(1)
c           ehmat(1,2)=cosl*vbi(2)
c           ehmat(1,3)=cosm*vbi(2)
c           ehmat(1,4)=cosn*vbi(2)
c           ehmat(2,2)=cosl**2*vbi(3)+(1.0d0-cosl**2)*vbi(4)
c           ehmat(2,3)=cosl*cosm*(vbi(3)-vbi(4))
c           ehmat(2,4)=cosl*cosn*(vbi(3)-vbi(4))
c           ehmat(3,3)=cosm**2*vbi(3)+(1.0d0-cosm**2)*vbi(4)
c           ehmat(3,4)=cosm*cosn*(vbi(3)-vbi(4))
c           ehmat(4,4)=cosn**2*vbi(3)+(1.0d0-cosn**2)*vbi(4)
c
c         do ia=1,no3-1
c          do ib=ia+1,no3
c           if( (ia.eq.1.and.(ib.ge.2.and.ib.le.4)).or.
c    &        ((ib.ge.5.and.ib.le.9).and.(ia.ge.2.and.ia.le.4)) ) then
c             ehmat(ib,ia)=-ehmat(ia,ib)
c           else
c             ehmat(ib,ia)=ehmat(ia,ib)
c           endif
c          enddo
c         enddo
c
               dx=dxij(jj,j,ii,nzz)
               dy=dyij(jj,j,ii,nzz)
               dz=dzij(jj,j,ii,nzz)
            do mk=1,nkp(nzz)
               v1=vk(1,mk,nzz)
               v2=vk(2,mk,nzz)
               v3=vk(3,mk,nzz)
              x1=+v1*ralat(1,1,nzz)+v2*ralat(2,1,nzz)
     &           +v3*ralat(3,1,nzz)
              y1=+v1*ralat(1,2,nzz)+v2*ralat(2,2,nzz)
     &           +v3*ralat(3,2,nzz)
              z1=+v1*ralat(1,3,nzz)+v2*ralat(2,3,nzz)
     &           +v3*ralat(3,3,nzz)
              pha=(x1*dx+y1*dy+z1*dz)
              prel=dcos(pha)
              pimg=dsin(pha)
c
             ITN=npa1(nzz)*no1+npa2(nzz)*no2
             do iit=ITN+(ii-npa1(nzz)-npa2(nzz)-1)*no3+1,
     &              ITN+(ii-npa1(nzz)-npa2(nzz)-1)*no3+4
                  kkt=iit-(ITN+(ii-npa1(nzz)-npa2(nzz)-1)*no3)
               do jjt=ITN+(jj-npa1(nzz)-npa2(nzz)-1)*no3+1,
     &                ITN+(jj-npa1(nzz)-npa2(nzz)-1)*no3+4

                 mmt=jjt-(ITN+(jj-npa1(nzz)-npa2(nzz)-1)*no3)
c                kkt=mod(iit-1,no3)+1
c                mmt=mod(jjt-1,no3)+1
c BCPan
c                write(19,*)'N-N',ii,jj,iit,jjt,kkt,mmt

                 hmatr(iit,jjt,1,mk)=hmatr(iit,jjt,1,mk)
     &                                  +ehmat(kkt,mmt)*prel
                 hmatr(iit,jjt,2,mk)=hmatr(iit,jjt,2,mk)
     &                                  +ehmat(kkt,mmt)*pimg
c BCPan:
c        if(NZZ.eq.34)then
c     write(19,*)'N-N ',ii,jj,iit,jjt,kkt,mmt,hmatr(iit,jjt,1,mk,nzz)
c        endif
c
c
              enddo
             enddo
c
            enddo  ! end K recycle
            endif 
        enddo
        endif  ! end of N-N
      enddo
c
      end

c **********************************************************************
c *      Use data to find the hamiltonian matrix elements              *
c **********************************************************************
      subroutine c_hmat12(guess,nzz)
      implicit double precision (a-h,o-z)
      parameter(mnum=100,neinum=150)
      parameter(mnp=30,no1=9,no2=1)
C WW, HH, WH
      parameter(ngu_WW=77,ngu_HH=6,ngu_WH=19)
c NN, WN, NH
      parameter(ngu_NN=18,ngu_NH=9,ngu_WN=67)
c
      parameter(mnkp=186)
      parameter(nz=41)
      common/con1/vk(3,mnkp,nz),nkp(nz)
      common/con2/alat(3,3,nz),ralat(3,3,nz),tou(mnp,3,nz)
     &      ,itou(mnp,nz),acon(nz),mmm(nz),npa1(nz),npa2(nz),npa3(nz),mm
      common/con3/pi,r_cut,rc_W,rc_H,rc_N,rc_WH,rc_NH,rc_WN,efe
c     common/data1/guess_WW(ngu_WW),guess_HH(ngu_HH)
c     common/data2/guess_WH(ngu_WH)
      common/location/al(3,mnp,mnum,nz),nt(mnp,mnum,nz)
      common/nei1/dxij(mnp,mnum,mnp,nz),dyij(mnp,mnum,mnp,nz)
     &           ,dzij(mnp,mnum,mnp,nz),rij0(mnp,mnum,mnp,nz)
     &           ,nn(mnp,nz),nb(neinum,mnp,nz),nbb(neinum,mnp,nz)
      common/nei2/rjk0(neinum,neinum,mnp,nz),nnc(neinum,mnp,nz)
     &           ,nbc(neinum,neinum,mnp,nz),nbbc(neinum,neinum,mnp,nz)
      common/hmat/hmatr(180,180,2,mnkp)
      common/eig/eigen(180,mnkp,nz)
c     common/hmat/hmatr(no1*mnp,no1*mnp,2,mnkp,nz)
c     common/eig/eigen(no1*mnp,mnkp,nz)
      common/dewe/dee(84,mnkp,nz),wee(84,mnkp,nz)
      common/ndewe/ndee(nz)
c
      dimension guess(ngu_WH)
      dimension vb(14)
      dimension arf1(14),arf2(14),arf3(14),arf4(14)
     &         ,arf6(14),arf7(14)
      dimension emat(no1,no1)
      dimension bta1(2)
c
c         arf5=guess(85)
c         bta1(1)=guess(86)
c         bta1(2)=guess(87)
c         bta2=guess(88)
c         bta3=guess(89)
 
         SQ3=dsqrt(3.0d0)
c     do m=1,14
      do m=1,3
         n=(m-1)*6
         arf1(m)=guess(n+1)
         arf2(m)=guess(n+2)
         arf3(m)=guess(n+3)
         arf4(m)=guess(n+4)
         arf6(m)=guess(n+5)
         arf7(m)=guess(n+6)
c         sor(m)=guess(n+7)
      enddo
         sor_WH=0.1d0
c
c     For W-H  Block.
        nnpam=no1*npa1(nzz)+no2*npa2(nzz)+no3*npa3(nzz)
        nnpa=npa1(nzz)+npa2(nzz)+npa3(nzz)
      do ii=1,nnpa
          ntii=nt(ii,1,nzz)
        if(ntii.eq.1) then  ! 1 for W.
        do jt=1,nn(ii,nzz)
           j=nb(jt,ii,nzz)
           jj=nbb(jt,ii,nzz)
           ntjj=nt(jj,j,nzz)
c
          if(ntjj.eq.2) then ! 2 for H.
           rij=rij0(jj,j,ii,nzz)
c BCPan
           CL=dxij(jj,j,ii,nzz)/rij
           CM=dyij(jj,j,ii,nzz)/rij
           CN=dzij(jj,j,ii,nzz)/rij
           CL2=CL*CL
           CM2=CM*CM
           CN2=CN*CN
c
c          esij=0.0d0
c BCPan:
c          do kt=1,nnc(jt,ii,nzz)
c            k=nbc(kt,jt,ii,nzz) !Count in neighbors of num i,ii to find the same neighbors
c            kk=nbbc(kt,jt,ii,nzz)
c            ntkk=nt(kk,k,nzz)
c
c            rjk=rjk0(kt,jt,ii,nzz)
c            rik=rij0(kk,k,ii,nzz)
c BCPan
c            esij=esij+rij/(rik+rjk)*bta1(ntkk)
c    *                *dexp(-bta2*((rik+rjk)/rij)**bta3)
c          enddo
c            delt=2.0d0/pi*datan(esij)
c         do L=1,14
          do L=1,3
c          rrij=rij*(1.0d0+arf5*delt+arf6(L)*delt**2+arf7(L)*delt**3)
           rrij=rij
           vb(L)=arf1(L)*rrij**arf2(L)*dexp(-arf3(L)*rrij**arf4(L))
     &                   /(1+dexp((rrij-rc_WH)/sor_WH))
c          if(nzz.eq.7)write(19,*)'rij= ',rij,vb(l)
          enddo
C BOND INTEGRELS SSS,SPS,PSS,PPS,PPP,SDS,DSS,PDS,DPS,PDP,DPP,DDS,DDP,DDD
          SSS=VB( 1)
          PSS=VB( 2)
          DSS=VB( 3)
c         SPS=VB( 2)
c         SDS=VB( 3)

c         PSS=VB( 3)
c         PPS=VB( 4)
c         PPP=VB( 5)
c         SDS=VB( 6)
c         DSS=VB( 7)
c         PDS=VB( 8)
c         DPS=VB( 9)
c         PDP=VB(10)
c         DPP=VB(11)
c         DDS=VB(12)
c         DDP=VB(13)
c         DDD=VB(14)
C
c
           EMAT(1,1)=+SSS
           EMAT(2,1)=-CL*PSS
           EMAT(3,1)=-CM*PSS
           EMAT(4,1)=-CN*PSS
           EMAT(5,1)=+SQ3*CL*CM*DSS
           EMAT(6,1)=+SQ3*CM*CN*DSS
           EMAT(7,1)=+SQ3*CN*CL*DSS
           EMAT(8,1)=+0.5D0*SQ3*(CL2-CM2)*DSS
           EMAT(9,1)=+(CN2-0.5D0*(CL2+CM2))*DSS
c          EMAT(1,1)=+SSS
c          EMAT(1,2)=+CL*SPS
c          EMAT(1,3)=+CM*SPS
c          EMAT(1,4)=+CN*SPS
c          EMAT(1,5)=+SQ3*CL*CM*SDS
c          EMAT(1,6)=+SQ3*CM*CN*SDS
c          EMAT(1,7)=+SQ3*CN*CL*SDS
c          EMAT(1,8)=+0.5D0*SQ3*(CL2-CM2)*SDS
c          EMAT(1,9)=+(CN2-0.5D0*(CL2+CM2))*SDS

c BCPan
c         if(nzz.eq.1)write(19,*)(vb(iii),iii=1,3),CL,CM,CN
c         if(nzz.eq.1)write(19,*)(emat(1,iii), iii=1,9)
c
c         EMAT(2,1)=-CL*PSS
c          EMAT(2,2)=+CL2*PPS+(1.0D0-CL2)*PPP
c          EMAT(2,3)=+CL*CM*(PPS-PPP)
c          EMAT(2,4)=+CL*CN*(PPS-PPP)
c          EMAT(2,5)=+SQ3*CL2*CM*PDS
c     &              +CM*(1.0D0-2.0D0*CL2)*PDP
c          EMAT(2,6)=+CL*CM*CN*(SQ3*PDS-2.0*PDP)
c          EMAT(2,7)=+SQ3*CL2*CN*PDS
c     &              +CN*(1.0D0-2.0D0*CL2)*PDP
c          EMAT(2,8)=+0.5D0*SQ3*CL*(CL2-CM2)*PDS
c     &              +CL*(1.0D0-CL2+CM2)*PDP
c          EMAT(2,9)=+CL*(CN2-0.5D0*(CL2+CM2))*PDS
c     &              -SQ3*CL*CN2*PDP
c
c         EMAT(3,1)=-CM*PSS
c          EMAT(3,2)=+EMAT(2,3)
c          EMAT(3,3)=+CM2*PPS+(1.0D0-CM2)*PPP
c          EMAT(3,4)=+CM*CN*(PPS-PPP)
c          EMAT(3,5)=+SQ3*CL*CM2*PDS
c     &              +CL*(1.0D0-2*CM2)*PDP
c          EMAT(3,6)=+SQ3*CM2*CN*PDS
c     &              +CN*(1.0D0-2*CM2)*PDP
c          EMAT(3,7)=+CL*CM*CN*(SQ3*PDS-2.0D0*PDP)
c          EMAT(3,8)=+0.5D0*SQ3*CM*(CL2-CM2)*PDS
c     &              -CM*(1.0D0+CL2-CM2)*PDP
c          EMAT(3,9)=+CM*(CN2-0.5D0*(CL2+CM2))*PDS
c     &              -SQ3*CM*CN2*PDP
c
c         EMAT(4,1)=-CN*PSS
c          EMAT(4,2)=+EMAT(2,4)
c          EMAT(4,3)=+EMAT(3,4)
c          EMAT(4,4)=+CN**2*PPS+(1.0D0-CN2)*PPP
c          EMAT(4,5)=+CL*CM*CN*(SQ3*PDS-2.0D0*PDP)
c          EMAT(4,6)=+SQ3*CM*CN2*PDS
c     &              +CM*(1.0D0-2.0D0*CN2)*PDP
c          EMAT(4,7)=+SQ3*CL*CN2*PDS
c     &              +CL*(1.0D0-2.0D0*CN2)*PDP
c          EMAT(4,8)=+0.5D0*SQ3*CN*(CL2-CM2)*PDS
c     &              -CN*(CL2-CM2)*PDP
c          EMAT(4,9)=+CN*(CN2-0.5D0*(CL2+CM2))*PDS
c     &              +SQ3*CN*(CL2+CM2)*PDP
c
c         EMAT(5,1)=+SQ3*CL*CM*DSS
c          EMAT(5,2)=-SQ3*CL2*CM*DPS
c     &              -CM*(1.0D0-2.0D0*CL2)*DPP
c          EMAT(5,3)=-SQ3*CL*CM2*DPS
c     &              -CL*(1.0D0-2*CM2)*DPP
c          EMAT(5,4)=-CL*CM*CN*(SQ3*DPS-2.0D0*DPP)
c          EMAT(5,5)=+3.0D0*CL2*CM2*DDS
c     &              +(CL2+CM2-4.0D0*CL2*CM2)*DDP
c     &              +(CN2+CL2*CM2)*DDD
c          EMAT(5,6)=+3.0D0*CL*CM2*CN*DDS
c     &              +CL*CN*(1.0D0-4.0D0*CM2)*DDP
c     &              +CL*CN*(CM2-1.0D0)*DDD
c          EMAT(5,7)=+3.0D0*CL2*CM*CN*DDS
c     &              +CM*CN*(1.0D0-4.0D0*CL2)*DDP
c     &              +CM*CN*(CL2-1.0D0)*DDD
c          EMAT(5,8)=+CL*CM*(CL2-CM2)*(1.5d0*DDS-2.0d0*DDP+0.5D0*DDD)
c          EMAT(5,9)=+SQ3*CL*CM*(CN2-0.5D0*(CL2+CM2))*DDS
c     &              -2.0D0*SQ3*CL*CM*CN2*DDP
c     &              +0.5D0*SQ3*CL*CM*(1.0D0+CN2)*DDD
c
c         EMAT(6,1)=+SQ3*CM*CN*DSS
c          EMAT(6,2)=-CL*CM*CN*(SQ3*DPS-2.0*DPP)
c          EMAT(6,3)=-SQ3*CM2*CN*DPS
c     &              -CN*(1.0D0-2*CM2)*DPP
c          EMAT(6,4)=-SQ3*CM*CN2*DPS
c     &              -CM*(1.0D0-2.0D0*CN2)*DPP
c          EMAT(6,5)=EMAT(5,6)
c          EMAT(6,6)=+3.0D0*CM2*CN2*DDS
c     &              +(CM2+CN2-4.0D0*CM2*CN2)*DDP
c     &              +(CL2+CM2*CN2)*DDD
c          EMAT(6,7)=+3.0D0*CL*CM*CN2*DDS
c     &              +CL*CM*(1.0D0-4.0D0*CN2)*DDP
c     &              +CL*CM*(CN2-1.0D0)*DDD
c          EMAT(6,8)=+1.5D0*CM*CN*(CL2-CM2)*DDS
c     &              -CM*CN*(1.0D0+2.0D0*(CL2-CM2))*DDP
c     &              +CM*CN*(1.0D0+0.5D0*(CL2-CM2))*DDD
c          EMAT(6,9)=+SQ3*CM*CN*(CN2-0.5D0*(CL2+CM2))*DDS
c     &              +SQ3*CM*CN*(CL2+CM2-CN2)*DDP
c     &              -0.5D0*SQ3*CM*CN*(CL2+CM2)*DDD
c
c         EMAT(7,1)=+SQ3*CN*CL*DSS
c          EMAT(7,2)=-SQ3*CL2*CN*DPS
c     &              -CN*(1.0D0-2.0D0*CL2)*DPP
c          EMAT(7,3)=-CL*CM*CN*(SQ3*DPS-2.0D0*DPP)
c          EMAT(7,4)=-SQ3*CL*CN2*DPS
c     &              -CL*(1.0D0-2.0D0*CN2)*DPP
c          EMAT(7,5)=EMAT(5,7)
c          EMAT(7,6)=EMAT(6,7)
c          EMAT(7,7)=+3.0D0*CL2*CN2*DDS
c     &              +(CL2+CN2-4.0D0*CL2*CN2)*DDP
c     &              +(CM2+CL2*CN2)*DDD
c          EMAT(7,8)=+1.5D0*CL*CN*(CL2-CM2)*DDS
c     &              +CL*CN*(1.0D0-2.0D0*(CL2-CM2))*DDP
c     &              -CL*CN*(1.0D0-0.5D0*(CL2-CM2))*DDD
c          EMAT(7,9)=+SQ3*CL*CN*(CN2-0.5D0*(CL2+CM2))*DDS
c     &              +SQ3*CL*CN*(CL2+CM2-CN2)*DDP
c     &              -0.5D0*SQ3*CL*CN*(CL2+CM2)*DDD
c
c         EMAT(8,1)=+0.5D0*SQ3*(CL2-CM2)*DSS
c          EMAT(8,2)=-0.5D0*SQ3*CL*(CL2-CM2)*DPS
c     &              -CL*(1.0D0-CL2+CM2)*DPP
c          EMAT(8,3)=-0.5D0*SQ3*CM*(CL2-CM2)*DPS
c     &              +CM*(1.0D0+CL2-CM2)*DPP
c          EMAT(8,4)=-0.5D0*SQ3*CN*(CL2-CM2)*DPS
c     &              +CN*(CL2-CM2)*DPP
c          EMAT(8,5)=EMAT(5,8)
c          EMAT(8,6)=EMAT(6,8)
c          EMAT(8,7)=EMAT(7,8)
c          EMAT(8,8)=+0.75D0*(CL2-CM2)**2*DDS
c     &              +(CL2+CM2-(CL2-CM2)**2)*DDP
c     &              +(CN2+0.25D0*(CL2-CM2)**2)*DDD
c          EMAT(8,9)=+0.5D0*SQ3*(CL2-CM2)*(CN2-0.5D0*(CL2+CM2))*DDS
c     &              +SQ3*CN2*(CM2-CL2)*DDP
c     &              +0.25d0*SQ3*(1.0D0+CN**2)*(CL2-CM2)*DDD
c
c         EMAT(9,1)=+(CN2-0.5D0*(CL2+CM2))*DSS
c          EMAT(9,2)=-CL*(CN2-0.5D0*(CL2+CM2))*DPS
c     &              +SQ3*CL*CN2*DPP
c          EMAT(9,3)=-CM*(CN2-0.5D0*(CL2+CM2))*DPS
c     &              +SQ3*CM*CN2*DPP
c          EMAT(9,4)=-CN*(CN2-0.5D0*(CL2+CM2))*DPS
c     &              -SQ3*CN*(CL2+CM2)*DPP
c          EMAT(9,5)=EMAT(5,9)
c          EMAT(9,6)=EMAT(6,9)
c          EMAT(9,7)=EMAT(7,9)
c          EMAT(9,8)=EMAT(8,9)
c          EMAT(9,9)=+(CN2-0.5D0*(CL2+CM2))**2*DDS
c     &              +3.0D0*CN**2*(CL2+CM2)*DDP
c     &              +0.75D0*(CL2+CM2)**2*DDD

          dx=dxij(jj,j,ii,nzz)
          dy=dyij(jj,j,ii,nzz)
          dz=dzij(jj,j,ii,nzz)
       do mk=1,nkp(nzz)
          v1=vk(1,mk,nzz)
          v2=vk(2,mk,nzz)
          v3=vk(3,mk,nzz)
         x1=+v1*ralat(1,1,nzz)+v2*ralat(2,1,nzz)
     &      +v3*ralat(3,1,nzz)
         y1=+v1*ralat(1,2,nzz)+v2*ralat(2,2,nzz)
     &      +v3*ralat(3,2,nzz)
         z1=+v1*ralat(1,3,nzz)+v2*ralat(2,3,nzz)
     &      +v3*ralat(3,3,nzz)
         pha=(x1*dx+y1*dy+z1*dz)
         prel=dcos(pha)
         pimg=dsin(pha)
c
            IHS = (ii-1)*no1 ! Row start point
            JHS = npa1(nzz)*no1+(jj-npa1(nzz))*no2 ! Column start point.
        do iit= IHS+1,IHS+no1
              kkt=iit-IHS
           do jjt=JHS,JHS
              mmt=jjt-JHS+1
c              kkt=mod(iit-1,no1)+1
c              mmt=mod(jjt-npa1(nzz)*no1-1,no2)+1
c BCPan:
c            write(19,*)'W-H ',ii,jj,iit,jjt,kkt,mmt
c
             hmatr(iit,jjt,1,mk)=hmatr(iit,jjt,1,mk)
     &                              +emat(kkt,mmt)*prel
             hmatr(iit,jjt,2,mk)=hmatr(iit,jjt,2,mk)
     &                              +emat(kkt,mmt)*pimg
C BCPan:
c            if(nzz.eq.5.and.mk.eq.1) write(19,19)nzz,ii,jj,iit,jjt,kkt,
c     1           mmt,hmatr(iit,jjt,1,mk,nzz)
 19         format(7i5,f12.6)
C
            enddo
          enddo
c
        enddo ! end K recycle
        endif 
       enddo
       endif  ! end of W-H.
      enddo
c     stop
c
      end
c **********************************************************************

c **********************************************************************
c *      Use data to find the hamiltonian matrix elements              *
c **********************************************************************
      subroutine c_hmat23(guess,nzz)
      implicit double precision (a-h,o-z)
      parameter(mnum=100,neinum=150)
      parameter(mnp=30,no1=9,no2=1,no3=4)
C WW, HH, WH
      parameter(ngu_WW=77,ngu_HH=6,ngu_WH=19)
c NN, WN, NH
      parameter(ngu_NN=18,ngu_NH=9,ngu_WN=67)
c
      parameter(mnkp=186)
      parameter(nz=41)
      common/con1/vk(3,mnkp,nz),nkp(nz)
      common/con2/alat(3,3,nz),ralat(3,3,nz),tou(mnp,3,nz),itou(mnp,nz),
     &      acon(nz),mmm(nz),npa1(nz),npa2(nz),npa3(nz),mm
      common/con3/pi,r_cut,rc_W,rc_H,rc_N,rc_WH,rc_NH,rc_WN,efe
c     common/data1/guess1(ngu1),guess2(ngu2)
c     common/data2/guess12(ngu3)
      common/location/al(3,mnp,mnum,nz),nt(mnp,mnum,nz)
      common/nei1/dxij(mnp,mnum,mnp,nz),dyij(mnp,mnum,mnp,nz)
     &           ,dzij(mnp,mnum,mnp,nz),rij0(mnp,mnum,mnp,nz)
     &           ,nn(mnp,nz),nb(neinum,mnp,nz),nbb(neinum,mnp,nz)
      common/nei2/rjk0(neinum,neinum,mnp,nz),nnc(neinum,mnp,nz)
     &           ,nbc(neinum,neinum,mnp,nz),nbbc(neinum,neinum,mnp,nz)
      common/hmat/hmatr(180,180,2,mnkp)
      common/eig/eigen(180,mnkp,nz)
c     common/hmat/hmatr(no1*mnp,no1*mnp,2,mnkp,nz)
c     common/eig/eigen(no1*mnp,mnkp,nz)
c     common/dewe/dee(20,mnkp,nz),wee(20,mnkp,nz),bund(ngu_tot,2),
c    &            del,ndee(nz),ndiv
c
      dimension guess(ngu_NH)
      dimension vb(14)
      dimension arf1(14),arf2(14),arf3(14),arf4(14)
     &         ,arf6(14),arf7(14)
      dimension emat(no1,no1)
      dimension bta1(2)
c
!         arf5=guess(85)
!         bta1(1)=guess(86)
!         bta1(2)=guess(87)
!         bta2=guess(88)
!         bta3=guess(89)
c
         SQ3=dsqrt(3.0d0)
      do m=1,2
         n=(m-1)*4
         arf1(m)=guess(n+1)
         arf2(m)=guess(n+2)
         arf3(m)=guess(n+3)
         arf4(m)=guess(n+4)
      enddo
         sor_NH=0.05d0
c
c     For N-H Block.
      nnpa=npa1(nzz)+npa2(nzz)+npa3(nzz)
      do ii=1,nnpa
          ntii=nt(ii,1,nzz)
        if(ntii.eq.2) then  ! 1 for H.
        do jt=1,nn(ii,nzz)
           j=nb(jt,ii,nzz)
           jj=nbb(jt,ii,nzz)
           ntjj=nt(jj,j,nzz)
c
          if(ntjj.eq.3) then ! 2 for N.
           rij=rij0(jj,j,ii,nzz)
c
           CL=dxij(jj,j,ii,nzz)/rij
           CM=dyij(jj,j,ii,nzz)/rij
           CN=dzij(jj,j,ii,nzz)/rij
           CL2=CL*CL
           CM2=CM*CM
           CN2=CN*CN
c
           esij=0.0d0
!          do kt=1,nnc(jt,ii,nzz)
!            k=nbc(kt,jt,ii,nzz) !Count in neighbors of num i,ii to find the same neighbors
!            kk=nbbc(kt,jt,ii,nzz)
!            ntkk=nt(kk,k,nzz)
c
!            rjk=rjk0(kt,jt,ii,nzz)
!            rik=rij0(kk,k,ii,nzz)
!            esij=esij+rij/(rik+rjk)*bta1(ntkk)
!     *                *dexp(-bta2*((rik+rjk)/rij)**bta3)
!          enddo
!            delt=2.0d0/pi*datan(esij)
          do L=1,2
!           rrij=rij*(1.0d0+arf5*delt+arf6(L)*delt**2+arf7(L)*delt**3)
           rrij=rij
           vb(L)=arf1(L)*rrij**arf2(L)*dexp(-arf3(L)*rrij**arf4(L))
     &                   /(1+dexp((rrij-rc_NH)/sor_NH))
c BCPan:
c        if(NZZ.eq.8)then
c        write(19,*)'test N-H ',nzz,L,rrij,arf1(L),arf2(L),arf3(L),a
c    &                         rf4(L),rc_NH,sor_NH
c        endif
c
          enddo
C BOND INTEGRELS SSS,SPS,PSS,PPS,PPP,SDS,DSS,PDS,DPS,PDP,DPP,DDS,DDP,DDD
          SSS=VB( 1)
          PSS=VB( 2)
C
          EMAT(1,1)=+SSS
c
          EMAT(2,1)=-CL*PSS
c
          EMAT(3,1)=-CM*PSS
c
          EMAT(4,1)=-CN*PSS
c
          EMAT(1,2)=CL*PSS
          EMAT(1,3)=CM*PSS
          EMAT(1,4)=CN*PSS
c
!          if(nzz.eq.1) then
!             write(10086,*) "emat(1,1)=", emat(1,1), ii,jj
!             write(10086,*) "emat(2,1)=", emat(2,1), ii,jj
!             write(10086,*) "emat(3,1)=", emat(3,1), ii,jj
!             write(10086,*) "emat(4,1)=", emat(4,1), ii,jj
!          endif
c
c         do ia=1,no3-1
c          do ib=ia+1,no3
c           if( (ia.eq.1.and.(ib.ge.2.and.ib.le.4)).or.
c    &        ((ib.ge.5.and.ib.le.9).and.(ia.ge.2.and.ia.le.4)) ) then
c             emat(ib,ia)=-emat(ia,ib)
c           else
c             emat(ib,ia)=emat(ia,ib)
c           endif
c          enddo
c         enddo
c
           dx=dxij(jj,j,ii,nzz)
           dy=dyij(jj,j,ii,nzz)
           dz=dzij(jj,j,ii,nzz)
       do mk=1,nkp(nzz)
           v1=vk(1,mk,nzz)
           v2=vk(2,mk,nzz)
           v3=vk(3,mk,nzz)
          x1=+v1*ralat(1,1,nzz)+v2*ralat(2,1,nzz)
     &      +v3*ralat(3,1,nzz)
          y1=+v1*ralat(1,2,nzz)+v2*ralat(2,2,nzz)
     &      +v3*ralat(3,2,nzz)
          z1=+v1*ralat(1,3,nzz)+v2*ralat(2,3,nzz)
     &      +v3*ralat(3,3,nzz)
          pha=(x1*dx+y1*dy+z1*dz)
          prel=dcos(pha)
          pimg=dsin(pha)
c        prel=1.0
c        pimg=0.0
c
         ITN=npa1(nzz)*no1+npa2(nzz)*no2
         ITH=npa1(nzz)*no1
       do iit=ITH+(ii-npa1(nzz))*no2,ITH+(ii-npa1(nzz))*no2
          kkt=1
          do jjt=ITN+(jj-npa1(nzz)-npa2(nzz)-1)*no3+1,
     &           ITN+(jj-npa1(nzz)-npa2(nzz)-1)*no3+4
             mmt=jjt-(ITN+(jj-npa1(nzz)-npa2(nzz)-1)*no3)
c
c             kkt=mod(iit-1,no2)+1
c             mmt=mod(jjt-npa1(nzz)*no1,no3)+1
          hmatr(iit,jjt,1,mk)=hmatr(iit,jjt,1,mk)
     &                              +emat(kkt,mmt)*prel
          hmatr(iit,jjt,2,mk)=hmatr(iit,jjt,2,mk)
     &                              +emat(kkt,mmt)*pimg
c BCPan:
c        if(NZZ.eq.5)then
c     write(19,*)'H-N ',ii,jj,iit,jjt,kkt,mmt,hmatr(iit,jjt,1,mk,nzz)
c        endif
c
c         hmatr(jjt,iit,1,mk,nzz)= hmatr(iit,jjt,1,mk,nzz)
c         hmatr(jjt,iit,2,mk,nzz)=-hmatr(iit,jjt,2,mk,nzz)
           enddo
         enddo
c
        enddo ! end K recycle
        endif 
       enddo
       endif  ! end of N-H.
      enddo
c
      end
c **********************************************************************
c *      Use data to find the hamiltonian matrix elements              *
c **********************************************************************
      subroutine c_hmat13(guess,nzz)
      implicit double precision (a-h,o-z)
      parameter(mnum=100,neinum=150)
      parameter(mnp=30,no1=9,no2=1,no3=4)
C WW, HH, WH
      parameter(ngu_WW=77,ngu_HH=6,ngu_WH=19)
c NN, WN, NH
      parameter(ngu_NN=18,ngu_NH=9,ngu_WN=67,ngu_tot=86)
c
      parameter(mnkp=186)
      parameter(nz=41)
      common/con1/vk(3,mnkp,nz),nkp(nz)
      common/con2/alat(3,3,nz),ralat(3,3,nz),tou(mnp,3,nz)
     &      ,itou(mnp,nz),acon(nz),mmm(nz),npa1(nz),npa2(nz),npa3(nz),mm
      common/con3/pi,r_cut,rc_W,rc_H,rc_N,rc_WH,rc_NH,rc_WN,efe
c     common/data1/guess_WW(ngu_WW),guess_HH(ngu_HH)
c     common/data2/guess_WH(ngu_WH)
      common/location/al(3,mnp,mnum,nz),nt(mnp,mnum,nz)
      common/nei1/dxij(mnp,mnum,mnp,nz),dyij(mnp,mnum,mnp,nz)
     &           ,dzij(mnp,mnum,mnp,nz),rij0(mnp,mnum,mnp,nz)
     &           ,nn(mnp,nz),nb(neinum,mnp,nz),nbb(neinum,mnp,nz)
      common/nei2/rjk0(neinum,neinum,mnp,nz),nnc(neinum,mnp,nz)
     &           ,nbc(neinum,neinum,mnp,nz),nbbc(neinum,neinum,mnp,nz)
      common/hmat/hmatr(180,180,2,mnkp)
      common/eig/eigen(180,mnkp,nz)
      common/dewe/dee(84,mnkp,nz),wee(84,mnkp,nz)
      common/ndewe/ndee(nz)
c     common/dewe/dee(20,mnkp,nz),wee(20,mnkp,nz),bund(ngu_tot,2),
c    &            del,ndee(nz),ndiv
c
      dimension guess(ngu_tot)
      dimension vb(14)
      dimension arf1(14),arf2(14),arf3(14),arf4(14),arf5(14)
     &         ,arf6(14),arf7(14)
      dimension emat(no1,no1)
      dimension bta1_WN(3)
c
c         arf5=guess(85)
          bta1_WN(1)=guess(57)
          bta1_WN(3)=guess(58)
          bta2_WN=guess(59)
          bta3_WN=guess(60)
 
         SQ3=dsqrt(3.0d0)
c     do m=1,14
      do m=1,8
         n=(m-1)*7
         arf1(m)=guess(n+1)
         arf2(m)=guess(n+2)
         arf3(m)=guess(n+3)
         arf4(m)=guess(n+4)
         arf5(m)=guess(n+5)
         arf6(m)=guess(n+6)
         arf7(m)=guess(n+7)
c        sor_WN(m)=guess(n+7)
      enddo
         sor_WN=0.1d0
c
c     For W-N  Block.
        nnpam=no1*npa1(nzz)+no2*npa2(nzz)+no3*npa3(nzz)
        nnpa=npa1(nzz)+npa2(nzz)+npa3(nzz)
      do ii=1,nnpa
          ntii=nt(ii,1,nzz)
        if(ntii.eq.1) then  ! 1 for W.
        do jt=1,nn(ii,nzz)
           j=nb(jt,ii,nzz)
           jj=nbb(jt,ii,nzz)
           ntjj=nt(jj,j,nzz)
c
          if(ntjj.eq.3) then ! 2 for N.
           rij=rij0(jj,j,ii,nzz)
c BCPan
           CL=dxij(jj,j,ii,nzz)/rij
           CM=dyij(jj,j,ii,nzz)/rij
           CN=dzij(jj,j,ii,nzz)/rij
           CL2=CL*CL
           CM2=CM*CM
           CN2=CN*CN
c
          esij=0.0d0
          do kt=1,nnc(jt,ii,nzz)
           k=nbc(kt,jt,ii,nzz) !Count in common neighbors
           kk=nbbc(kt,jt,ii,nzz)
           ntkk=nt(kk,k,nzz)
c
           if(ntkk.eq.1.or.ntkk.eq.3)then
           rjk=rjk0(kt,jt,ii,nzz)
           rik=rij0(kk,k,ii,nzz)
c          esij=esij+rij/(rik+rjk)*bta1(ntkk)
           esij=esij+rij/(rik+rjk)*bta1_WN(ntkk)
     *               *dexp(-bta2_WN*((rik+rjk)/rij)**bta3_WN)
           endif
           enddo
           delt=2.0d0/pi*datan(esij)
c
          do L=1,8
c          rrij=rij*(1.0d0+arf5(L)*delt+arf6(L)*delt**2+arf7(L)*delt**3)
           rrij=rij*(1.0d0+arf5(L)*delt+arf6(L)*delt**2+arf7(L)*delt**3+
     &               guess(61)*delt**4+guess(62)*delt**5)
c          rrij=rij
           vb(L)=arf1(L)*rrij**arf2(L)*dexp(-arf3(L)*rrij**arf4(L))
     &                   /(1+dexp((rrij-rc_WN)/sor_WN))
c          if(nzz.eq.1)write(19,*)'rij= ',ii,jj,rij,L,vb(L),
c    &        arf1(L),arf2(L),arf3(L),arf4(L),rc_WN,sor_WN
          enddo
C BOND INTEGRELS SSS,SPS,PSS,PPS,PPP,SDS,DSS,PDS,DPS,PDP,DPP,DDS,DDP,DDD
          SSS=VB(1)     ! used
          SPS=VB(2)     ! used
          PSS=VB(3)     ! no use W p is too high in energy
          PPS=VB(4)     ! no use
          PPP=VB(5)     ! no use
          DSS=VB(6)     ! used
          DPS=VB(7)     ! used
          DPP=VB(8)     ! used

c         SSS=VB( 1)
c         PSS=VB( 2)
c         DSS=VB( 3)
c         SPS=VB( 2)
c         SDS=VB( 3)

c         PSS=VB( 3)
c         PPS=VB( 4)
c         PPP=VB( 5)
c         SDS=VB( 6)
c         DSS=VB( 7)
c         PDS=VB( 8)
c         DPS=VB( 9)
c         PDP=VB(10)
c         DPP=VB(11)
c         DDS=VB(12)
c         DDP=VB(13)
c         DDD=VB(14)
C
c
           EMAT(1,1)=+SSS
c          EMAT(1,2)=CL*PSS
c          EMAT(1,3)=CM*PSS
c          EMAT(1,4)=CN*PSS
           EMAT(1,2)=CL*SPS
           EMAT(1,3)=CM*SPS
           EMAT(1,4)=CN*SPS
c          EMAT(1,1)=+SSS
c          EMAT(1,2)=CL*SPS
c          EMAT(1,3)=CM*SPS
c          EMAT(1,4)=CN*SPS
c          EMAT(1,5)=+SQ3*CL*CM*DSS
c          EMAT(1,6)=+SQ3*CM*CN*DSS
c          EMAT(1,7)=+SQ3*CN*CL*DSS
c          EMAT(1,8)=+0.5D0*SQ3*(CL2-CM2)*DSS
c          EMAT(1,9)=+(CN2-0.5D0*(CL2+CM2))*DSS
c          EMAT(1,1)=+SSS
c          EMAT(1,2)=+CL*SPS
c          EMAT(1,3)=+CM*SPS
c          EMAT(1,4)=+CN*SPS
c          EMAT(1,5)=+SQ3*CL*CM*SDS
c          EMAT(1,6)=+SQ3*CM*CN*SDS
c          EMAT(1,7)=+SQ3*CN*CL*SDS
c          EMAT(1,8)=+0.5D0*SQ3*(CL2-CM2)*SDS
c          EMAT(1,9)=+(CN2-0.5D0*(CL2+CM2))*SDS

c BCPan
c
           EMAT(2,1)=-CL*PSS
           EMAT(2,2)=+CL2*PPS+(1.0D0-CL2)*PPP
           EMAT(2,3)=+CL*CM*(PPS-PPP)
           EMAT(2,4)=+CL*CN*(PPS-PPP)
c          EMAT(2,5)=+SQ3*CL2*CM*PDS
c     &              +CM*(1.0D0-2.0D0*CL2)*PDP
c          EMAT(2,6)=+CL*CM*CN*(SQ3*PDS-2.0*PDP)
c          EMAT(2,7)=+SQ3*CL2*CN*PDS
c     &              +CN*(1.0D0-2.0D0*CL2)*PDP
c          EMAT(2,8)=+0.5D0*SQ3*CL*(CL2-CM2)*PDS
c     &              +CL*(1.0D0-CL2+CM2)*PDP
c          EMAT(2,9)=+CL*(CN2-0.5D0*(CL2+CM2))*PDS
c     &              -SQ3*CL*CN2*PDP
c
           EMAT(3,1)=-CM*PSS
           EMAT(3,2)=+EMAT(2,3)
           EMAT(3,3)=+CM2*PPS+(1.0D0-CM2)*PPP
           EMAT(3,4)=+CM*CN*(PPS-PPP)
c          EMAT(3,5)=+SQ3*CL*CM2*PDS
c     &              +CL*(1.0D0-2*CM2)*PDP
c          EMAT(3,6)=+SQ3*CM2*CN*PDS
c     &              +CN*(1.0D0-2*CM2)*PDP
c          EMAT(3,7)=+CL*CM*CN*(SQ3*PDS-2.0D0*PDP)
c          EMAT(3,8)=+0.5D0*SQ3*CM*(CL2-CM2)*PDS
c     &              -CM*(1.0D0+CL2-CM2)*PDP
c          EMAT(3,9)=+CM*(CN2-0.5D0*(CL2+CM2))*PDS
c     &              -SQ3*CM*CN2*PDP
c
           EMAT(4,1)=-CN*PSS
           EMAT(4,2)=+EMAT(2,4)
           EMAT(4,3)=+EMAT(3,4)
           EMAT(4,4)=+CN**2*PPS+(1.0D0-CN2)*PPP
c          EMAT(4,5)=+CL*CM*CN*(SQ3*PDS-2.0D0*PDP)
c          EMAT(4,6)=+SQ3*CM*CN2*PDS
c     &              +CM*(1.0D0-2.0D0*CN2)*PDP
c          EMAT(4,7)=+SQ3*CL*CN2*PDS
c     &              +CL*(1.0D0-2.0D0*CN2)*PDP
c          EMAT(4,8)=+0.5D0*SQ3*CN*(CL2-CM2)*PDS
c     &              -CN*(CL2-CM2)*PDP
c          EMAT(4,9)=+CN*(CN2-0.5D0*(CL2+CM2))*PDS
c     &              +SQ3*CN*(CL2+CM2)*PDP
c
           EMAT(5,1)=+SQ3*CL*CM*DSS
           EMAT(5,2)=-SQ3*CL2*CM*DPS
     &              -CM*(1.0D0-2.0D0*CL2)*DPP
           EMAT(5,3)=-SQ3*CL*CM2*DPS
     &              -CL*(1.0D0-2*CM2)*DPP
           EMAT(5,4)=-CL*CM*CN*(SQ3*DPS-2.0D0*DPP)
c          EMAT(5,5)=+3.0D0*CL2*CM2*DDS
c     &              +(CL2+CM2-4.0D0*CL2*CM2)*DDP
c     &              +(CN2+CL2*CM2)*DDD
c          EMAT(5,6)=+3.0D0*CL*CM2*CN*DDS
c     &              +CL*CN*(1.0D0-4.0D0*CM2)*DDP
c     &              +CL*CN*(CM2-1.0D0)*DDD
c          EMAT(5,7)=+3.0D0*CL2*CM*CN*DDS
c     &              +CM*CN*(1.0D0-4.0D0*CL2)*DDP
c     &              +CM*CN*(CL2-1.0D0)*DDD
c          EMAT(5,8)=+CL*CM*(CL2-CM2)*(1.5d0*DDS-2.0d0*DDP+0.5D0*DDD)
c          EMAT(5,9)=+SQ3*CL*CM*(CN2-0.5D0*(CL2+CM2))*DDS
c     &              -2.0D0*SQ3*CL*CM*CN2*DDP
c     &              +0.5D0*SQ3*CL*CM*(1.0D0+CN2)*DDD
c
           EMAT(6,1)=+SQ3*CM*CN*DSS
           EMAT(6,2)=-CL*CM*CN*(SQ3*DPS-2.0*DPP)
           EMAT(6,3)=-SQ3*CM2*CN*DPS
     &              -CN*(1.0D0-2*CM2)*DPP
           EMAT(6,4)=-SQ3*CM*CN2*DPS
     &              -CM*(1.0D0-2.0D0*CN2)*DPP
c          EMAT(6,5)=EMAT(5,6)
c          EMAT(6,6)=+3.0D0*CM2*CN2*DDS
c     &              +(CM2+CN2-4.0D0*CM2*CN2)*DDP
c     &              +(CL2+CM2*CN2)*DDD
c          EMAT(6,7)=+3.0D0*CL*CM*CN2*DDS
c     &              +CL*CM*(1.0D0-4.0D0*CN2)*DDP
c     &              +CL*CM*(CN2-1.0D0)*DDD
c          EMAT(6,8)=+1.5D0*CM*CN*(CL2-CM2)*DDS
c     &              -CM*CN*(1.0D0+2.0D0*(CL2-CM2))*DDP
c     &              +CM*CN*(1.0D0+0.5D0*(CL2-CM2))*DDD
c          EMAT(6,9)=+SQ3*CM*CN*(CN2-0.5D0*(CL2+CM2))*DDS
c     &              +SQ3*CM*CN*(CL2+CM2-CN2)*DDP
c     &              -0.5D0*SQ3*CM*CN*(CL2+CM2)*DDD
c
           EMAT(7,1)=+SQ3*CN*CL*DSS
           EMAT(7,2)=-SQ3*CL2*CN*DPS
     &              -CN*(1.0D0-2.0D0*CL2)*DPP
           EMAT(7,3)=-CL*CM*CN*(SQ3*DPS-2.0D0*DPP)
           EMAT(7,4)=-SQ3*CL*CN2*DPS
     &              -CL*(1.0D0-2.0D0*CN2)*DPP
c          EMAT(7,5)=EMAT(5,7)
c          EMAT(7,6)=EMAT(6,7)
c          EMAT(7,7)=+3.0D0*CL2*CN2*DDS
c     &              +(CL2+CN2-4.0D0*CL2*CN2)*DDP
c     &              +(CM2+CL2*CN2)*DDD
c          EMAT(7,8)=+1.5D0*CL*CN*(CL2-CM2)*DDS
c     &              +CL*CN*(1.0D0-2.0D0*(CL2-CM2))*DDP
c     &              -CL*CN*(1.0D0-0.5D0*(CL2-CM2))*DDD
c          EMAT(7,9)=+SQ3*CL*CN*(CN2-0.5D0*(CL2+CM2))*DDS
c     &              +SQ3*CL*CN*(CL2+CM2-CN2)*DDP
c     &              -0.5D0*SQ3*CL*CN*(CL2+CM2)*DDD
c
           EMAT(8,1)=+0.5D0*SQ3*(CL2-CM2)*DSS
           EMAT(8,2)=-0.5D0*SQ3*CL*(CL2-CM2)*DPS
     &              -CL*(1.0D0-CL2+CM2)*DPP
           EMAT(8,3)=-0.5D0*SQ3*CM*(CL2-CM2)*DPS
     &              +CM*(1.0D0+CL2-CM2)*DPP
           EMAT(8,4)=-0.5D0*SQ3*CN*(CL2-CM2)*DPS
     &              +CN*(CL2-CM2)*DPP
c          EMAT(8,5)=EMAT(5,8)
c          EMAT(8,6)=EMAT(6,8)
c          EMAT(8,7)=EMAT(7,8)
c          EMAT(8,8)=+0.75D0*(CL2-CM2)**2*DDS
c     &              +(CL2+CM2-(CL2-CM2)**2)*DDP
c     &              +(CN2+0.25D0*(CL2-CM2)**2)*DDD
c          EMAT(8,9)=+0.5D0*SQ3*(CL2-CM2)*(CN2-0.5D0*(CL2+CM2))*DDS
c     &              +SQ3*CN2*(CM2-CL2)*DDP
c     &              +0.25d0*SQ3*(1.0D0+CN**2)*(CL2-CM2)*DDD
c
           EMAT(9,1)=+(CN2-0.5D0*(CL2+CM2))*DSS
           EMAT(9,2)=-CL*(CN2-0.5D0*(CL2+CM2))*DPS
     &              +SQ3*CL*CN2*DPP
           EMAT(9,3)=-CM*(CN2-0.5D0*(CL2+CM2))*DPS
     &              +SQ3*CM*CN2*DPP
           EMAT(9,4)=-CN*(CN2-0.5D0*(CL2+CM2))*DPS
     &              -SQ3*CN*(CL2+CM2)*DPP
c          EMAT(9,5)=EMAT(5,9)
c          EMAT(9,6)=EMAT(6,9)
c          EMAT(9,7)=EMAT(7,9)
c          EMAT(9,8)=EMAT(8,9)
c          EMAT(9,9)=+(CN2-0.5D0*(CL2+CM2))**2*DDS
c     &              +3.0D0*CN**2*(CL2+CM2)*DDP
c     &              +0.75D0*(CL2+CM2)**2*DDD

          dx=dxij(jj,j,ii,nzz)
          dy=dyij(jj,j,ii,nzz)
          dz=dzij(jj,j,ii,nzz)
          do mk=1,nkp(nzz)
          v1=vk(1,mk,nzz)
          v2=vk(2,mk,nzz)
          v3=vk(3,mk,nzz)
c                write(6,*)'test Kpoints ',nzz,mk,v1,v2,v3,dx,dy,dz
         x1=+v1*ralat(1,1,nzz)+v2*ralat(2,1,nzz)
     &      +v3*ralat(3,1,nzz)
         y1=+v1*ralat(1,2,nzz)+v2*ralat(2,2,nzz)
     &      +v3*ralat(3,2,nzz)
         z1=+v1*ralat(1,3,nzz)+v2*ralat(2,3,nzz)
     &      +v3*ralat(3,3,nzz)
         pha=(x1*dx+y1*dy+z1*dz)
         prel=dcos(pha)
         pimg=dsin(pha)
c
           IHS = (ii-1)*no1 ! Row start point
           ITN=npa1(nzz)*no1+npa2(nzz)*no2

        do iit= IHS+1,IHS+no1
           kkt=iit-IHS
           do jjt=ITN+(jj-npa1(nzz)-npa2(nzz)-1)*no3+1,
     &            ITN+(jj-npa1(nzz)-npa2(nzz)-1)*no3+4
             mmt=jjt-(ITN+(jj-npa1(nzz)-npa2(nzz)-1)*no3)
c             kkt=mod(iit-1,no1)+1
c             mmt=mod(jjt-2,no3)+1
             hmatr(iit,jjt,1,mk)=hmatr(iit,jjt,1,mk)
     &                              +emat(kkt,mmt)*prel
             hmatr(iit,jjt,2,mk)=hmatr(iit,jjt,2,mk)
     &                              +emat(kkt,mmt)*pimg
c
             hmatr(jjt,iit,1,mk)=hmatr(jjt,iit,1,mk)
     &                              +emat(mmt,kkt)*prel
             hmatr(jjt,iit,2,mk)=hmatr(jjt,iit,2,mk)
     &                              -emat(mmt,kkt)*pimg

c BCPan:
c        if(NZZ.eq.5)then
c     write(19,*)'W-N ',ii,jj,iit,jjt,kkt,mmt,hmatr(iit,jjt,1,mk,nzz)
c        endif
c

 19         format(4i5,f12.6)
C
            enddo
          enddo
c
        enddo ! end K recycle
        endif 
       enddo
       endif  ! end of W-N.
      enddo
c     stop
c
      end
c **********************************************************************
c *                   This is the  program band                        *
c **********************************************************************
      subroutine c_band(hmatr,nzz) 
      implicit double precision (a-h,o-z)
      parameter(mnum=100,neinum=150)
      parameter(mnp=30,no1=9,no2=1,no3=4)
C WW, HH, WH
      parameter(ngu_WW=77,ngu_HH=6,ngu_WH=19)
c NN, WN, NH
      parameter(ngu_NN=18,ngu_NH=9,ngu_WN=67)
c
      parameter(mnkp=186)
      parameter(nz=41)
      common/con1/vk(3,mnkp,nz),nkp(nz)
      common/con2/alat(3,3,nz),ralat(3,3,nz),tou(mnp,3,nz)
     &      ,itou(mnp,nz),acon(nz),mmm(nz),npa1(nz),npa2(nz),npa3(nz),mm
      common/con3/pi,r_cut,rc_W,rc_H,rc_N,rc_WH,rc_NH,rc_WN,efe
c WW, HH:
c     common/data1/guess_WW(ngu_WW),guess_HH(ngu_HH)
c WH:
c     common/data2/guess_WH(ngu_WH)
c NN, NH:
c     common/data3/guess_NN(ngu_NN),guess_NH(ngu_NH)
c WN
c     common/data4/guess_WN(ngu_WN)
c
      common/location/al(3,mnp,mnum,nz),nt(mnp,mnum,nz)
      common/nei1/dxij(mnp,mnum,mnp,nz),dyij(mnp,mnum,mnp,nz)
     &           ,dzij(mnp,mnum,mnp,nz),rij0(mnp,mnum,mnp,nz)
     &           ,nn(mnp,nz),nb(neinum,mnp,nz),nbb(neinum,mnp,nz)
      common/nei2/rjk0(neinum,neinum,mnp,nz),nnc(neinum,mnp,nz)
     &           ,nbc(neinum,neinum,mnp,nz),nbbc(neinum,neinum,mnp,nz)
c     common/hmat/hmatr(180,180,2,mnkp,nz)
      common/eig/eigen(180,mnkp,nz)
c     common/hmat/hmatr(no*mnp,no*mnp,2,mnkp,nz)
c     common/eig/eigen(no1*mnp,mnkp,nz)
c     common/dewe/dee(20,mnkp,nz),wee(20,mnkp,nz),bund(ngu_tot,2),
c    &            del,ndee(nz),ndiv
c
      dimension hmatr(180,180,2,mnkp)
c     complex*16 ap(no1*mnp*(no1*mnp+1)/2)
c     dimension RWORK(3*no1*mnp-2),W(no1*mnp)
      complex*16 ap(180*(180+1)/2)
      dimension RWORK(3*180-2),W(180)
      character V,U
c     complex*16 WORK(2*no1*mnp-1)
c     complex*16 Z(no1*mnp,no1*mnp)
      complex*16 WORK(2*180-1)
      complex*16 Z(180,180)
c
       Mspace=no1*npa1(nzz)+no2*npa2(nzz)+no3*npa3(nzz)
c BCPan
c      if(nzz.eq.7)write(19,*)'Mspace= ',Mspace
         LDZ=Mspace
      do mk=1,nkp(nzz)
          nap=0
c        if(NZZ.eq.16.and.mk.eq.1)write(19,*)'test 2'
c        if(NZZ.eq.16.and.mk.eq.1)then
c     write(19,1212)((hmatr(mmt,kkt,1,1,nzz),mmt=1,Mspace),kkt=1,Mspace)
c        endif
1212   format(23f6.2)
c
        do kkt=1,Mspace
         do mmt=1,kkt
          nap=nap+1
          ap(nap)=cmplx(hmatr(mmt,kkt,1,mk),hmatr(mmt,kkt,2,mk))
c BCPan:
c         if(nzz.eq.8)then
c           write(19,*)nzz,mk,mmt,kkt,nap,ap(nap)
c         endif
c
         enddo
        enddo
c BCPan
c      if(NZZ.eq.8)write(19,1212)(ap(ik),ik=1,nap,2)
c      write(*,*)'test zhpev',nzz
       call zhpev("V","U",Mspace,ap,W,Z,LDZ,WORK,RWORK,INFO)
c       write(19,1212)'test zhpev-2',nzz
c
        do ik=1,Mspace
         do im=ik+1,Mspace
           t1=W(ik)
           t2=W(im)
           if(t1.gt.t2) then
           T=W(ik)
           W(ik)=W(im)
           W(im)=T
           else
           endif
         enddo
        enddo
c
        do ii=1,Mspace
           eigen(ii,mk,nzz)=W(ii)
        enddo
c     if(NZZ.eq.5)write(19,*)(w(ii),ii=1,Mspace)
c     if(nzz.eq.7)stop 
      enddo
103    format(3f12.8)
c
      end
c **********************************************************************


c***********************************************************************
c*        call the subprogram to calculate the eigenvalue              *
c**************************of the tbmatrix******************************
      SUBROUTINE ZHPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK,
     $                  INFO )
*
*  -- LAPACK driver routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK(*), W(*)
      COMPLEX*16         AP(*), WORK(*), Z( LDZ,*)
*     ..
*
*  Purpose
*  =======
*
*  ZHPEV computes all the eigenvalues and, optionally, eigenvectors of a
*  complex Hermitian matrix in packed storage.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)
*          On entry, the upper or lower triangle of the Hermitian matrix
*          A, packed columnwise in a linear array.  The j-th column of A
*          is stored in the array AP as follows:
*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
*
*          On exit, AP is overwritten by values generated during the
*          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
*          and first superdiagonal of the tridiagonal matrix T overwrite
*          the corresponding elements of A, and if UPLO = 'L', the
*          diagonal and first subdiagonal of T overwrite the
*          corresponding elements of A.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  Z       (output) COMPLEX*16 array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
*          eigenvectors of the matrix A, with the i-th column of Z
*          holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (max(1, 2*N-1))
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(1, 3*N-2))
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            WANTZ
      INTEGER            IINFO, IMAX, INDE, INDRWK, INDTAU, INDWRK,
     $                   ISCALE
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANHP
      EXTERNAL           LSAME, DLAMCH, ZLANHP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSTERF, XERBLA, ZDSCAL, ZHPTRD, ZSTEQR,
     $                   ZUPGTR
*     ..
*     .. Intrinsic Functions ..
*PAN  INTRINSIC          SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( UPLO, 'L' ) .OR. LSAME( UPLO, 'U' ) ) )
     $          THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -7
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHPEV ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = AP( 1 )
         RWORK( 1 ) = 1
         IF( WANTZ )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = ZLANHP( 'M', UPLO, N, AP, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         CALL ZDSCAL( ( N*( N+1 ) ) / 2, SIGMA, AP, 1 )
      END IF
*
*     Call ZHPTRD to reduce Hermitian packed matrix to tridiagonal form.
*
      INDE = 1
      INDTAU = 1
      CALL ZHPTRD( UPLO, N, AP, W, RWORK( INDE ), WORK( INDTAU ),
     $             IINFO )
*
*     For eigenvalues only, call DSTERF.  For eigenvectors, first call
*     ZUPGTR to generate the orthogonal matrix, then call ZSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         INDWRK = INDTAU + N
         CALL ZUPGTR( UPLO, N, AP, WORK( INDTAU ), Z, LDZ,
     $                WORK( INDWRK ), IINFO )
         INDRWK = INDE + N
         CALL ZSTEQR( JOBZ, N, W, RWORK( INDE ), Z, LDZ,
     $                RWORK( INDRWK ), INFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
      RETURN
*
*     End of ZHPEV
*
      END
C

      SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, RT1, RT2
*     ..
*
*  Purpose
*  =======
*
*  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, and RT2
*  is the eigenvalue of smaller absolute value.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  B       (input) DOUBLE PRECISION
*          The (1,2) and (2,1) elements of the 2-by-2 matrix.
*
*  C       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  RT1     (output) DOUBLE PRECISION
*          The eigenvalue of larger absolute value.
*
*  RT2     (output) DOUBLE PRECISION
*          The eigenvalue of smaller absolute value.
*
*  Further Details
*  ===============
*
*  RT1 is accurate to a few ulps barring over/underflow.
*
*  RT2 may be inaccurate if there is massive cancellation in the
*  determinant A*C-B*B; higher precision or correctly rounded or
*  correctly truncated arithmetic would be needed to compute RT2
*  accurately in all cases.
*
*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*  Underflow is harmless if the input data is 0 or exceeds
*     underflow_threshold / macheps.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
*     ..
*     .. Intrinsic Functions ..
*PAN     INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
*
*     End of DLAE2
*
      END
      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA,*)
*     ..
*
*  Purpose
*  =======
*
*  DLASCL multiplies the M by N real matrix A by the real scalar
*  CTO/CFROM.  This is done without over/underflow as long as the final
*  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
*  A may be full, upper triangular, lower triangular, upper Hessenberg,
*  or banded.
*
*  Arguments
*  =========
*
*  TYPE    (input) CHARACTER*1
*          TYPE indices the storage type of the input matrix.
*          = 'G':  A is a full matrix.
*          = 'L':  A is a lower triangular matrix.
*          = 'U':  A is an upper triangular matrix.
*          = 'H':  A is an upper Hessenberg matrix.
*          = 'B':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the lower
*                  half stored.
*          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the upper
*                  half stored.
*          = 'Z':  A is a band matrix with lower bandwidth KL and upper
*                  bandwidth KU.
*
*  KL      (input) INTEGER
*          The lower bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  KU      (input) INTEGER
*          The upper bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  CFROM   (input) DOUBLE PRECISION
*  CTO     (input) DOUBLE PRECISION
*          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
*          without over/underflow if the final result CTO*A(I,J)/CFROM
*          can be represented without over/underflow.  CFROM must be
*          nonzero.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)
*          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
*          storage type.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  INFO    (output) INTEGER
*          0  - successful exit
*          <0 - if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. Intrinsic Functions ..
* PAN INTRINSIC          ABS, MAX, MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
*
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
*
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
      CFROMC = CFROM
      CTOC = CTO
*
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
      END IF
*
      IF( ITYPE.EQ.0 ) THEN
*
*        Full matrix
*
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( ITYPE.EQ.1 ) THEN
*
*        Lower triangular matrix
*
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Upper triangular matrix
*
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Upper Hessenberg matrix
*
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
*
      ELSE IF( ITYPE.EQ.4 ) THEN
*
*        Lower half of a symmetric band matrix
*
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
*
      ELSE IF( ITYPE.EQ.5 ) THEN
*
*        Upper half of a symmetric band matrix
*
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
*
      ELSE IF( ITYPE.EQ.6 ) THEN
*
*        Band matrix
*
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
*
      END IF
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of DLASCL
*
      END
      SUBROUTINE DSTERF( N, D, E, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D(*), E(*)
*     ..
*
*  Purpose
*  =======
*
*  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
*  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm failed to find all of the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDM1, LENDP1,
     $                   LENDSV, LM1, LSV, M, MM1, NM1, NMAXIT
      DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC,
     $                   OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN,
     $                   SIGMA, SSFMAX, SSFMIN, TST
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           DLAMCH, DLANST, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
*     ..
*     .. Intrinsic Functions ..
*PAN     INTRINSIC          ABS, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DSTERF', -INFO )
         RETURN
      END IF
      IF( N.LE.1 )
     $   RETURN
*
*     Determine the unit roundoff for this environment.
*
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues of the tridiagonal matrix.
*
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 170
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
      DO 40 I = L, LEND - 1
         E( I ) = E( I )**2
   40 CONTINUE
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GE.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   50    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 60 M = L, LENDM1
               TST = ABS( E( M ) )
               IF( TST.LE.EPS2*ABS( D( M )*D( M+1 ) ) )
     $            GO TO 70
   60       CONTINUE
         END IF
*
         M = LEND
*
   70    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 90
*
*        If remaining matrix is 2 by 2, use DLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L+1 ) THEN
            RTE = SQRT( E( L ) )
            CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 50
            GO TO 150
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         MM1 = M - 1
         DO 80 I = MM1, L, -1
            BB = E( I )
            R = P + BB
            IF( I.NE.M-1 )
     $         E( I+1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
   80    CONTINUE
*
         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50
*
*        Eigenvalue found.
*
   90    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 50
         GO TO 150
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
  100    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 110 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )
               IF( TST.LE.EPS2*ABS( D( M )*D( M-1 ) ) )
     $            GO TO 120
  110       CONTINUE
         END IF
*
         M = LEND
*
  120    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 140
*
*        If remaining matrix is 2 by 2, use DLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L-1 ) THEN
            RTE = SQRT( E( L-1 ) )
            CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 100
            GO TO 150
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         LM1 = L - 1
         DO 130 I = M, LM1
            BB = E( I )
            R = P + BB
            IF( I.NE.M )
     $         E( I-1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
  130    CONTINUE
*
         E( LM1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100
*
*        Eigenvalue found.
*
  140    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 100
         GO TO 150
*
      END IF
*
*     Undo scaling if necessary
*
  150 CONTINUE
      IF( ISCALE.EQ.1 )
     $   CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
      IF( ISCALE.EQ.2 )
     $   CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.EQ.NMAXIT ) THEN
         DO 160 I = 1, N - 1
            IF( E( I ).NE.ZERO )
     $         INFO = INFO + 1
  160    CONTINUE
         RETURN
      END IF
      GO TO 10
*
*     Sort eigenvalues in increasing order.
*
  170 CONTINUE
      CALL DLASRT( 'I', N, D, INFO )
*
      RETURN
*
*     End of DSTERF
*
      END
      subroutine  zdscal(n,da,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*)
      double precision da
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = dcmplx(da,0.0d0)*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = dcmplx(da,0.0d0)*zx(i)
   30 continue
      return
      end
      SUBROUTINE ZLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
c
c  -- LAPACK auxiliary routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c     September 30, 1994
c
c     .. Scalar Arguments ..
*PAN   INTRINSIC          ABS
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*16         TAU
c     ..
c     .. Array Arguments ..
      COMPLEX*16  C(LDC,*), V(*), WORK(*)
c     ..
c
c  Purpose
c  =======
c
c  ZLARF applies a complex elementary reflector H to a complex M-by-N
c  matrix C, from either the left or the right. H is represented in the
c  form
c
c        H = I - tau * v * v'
c
c  where tau is a complex scalar and v is a complex vector.
c
c  If tau = 0, then H is taken to be the unit matrix.
c
c  To apply H' (the conjugate transpose of H), supply conjg(tau) instead
c  tau.
c
c  Arguments
c  =========
c
c  SIDE    (input) CHARACTER*1
c          = 'L': form  H * C
c          = 'R': form  C * H
c
c  M       (input) INTEGER
c          The number of rows of the matrix C.
c
c  N       (input) INTEGER
c          The number of columns of the matrix C.
c
c  V       (input) COMPLEX*16 array, dimension
c                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
c                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
c          The vector v in the representation of H. V is not used if
c          TAU = 0.
c
c  INCV    (input) INTEGER
c          The increment between elements of v. INCV <> 0.
c
c  TAU     (input) COMPLEX*16
c          The value tau in the representation of H.
c
c  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
c          On entry, the M-by-N matrix C.
c          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
c          or C * H if SIDE = 'R'.
c
c  LDC     (input) INTEGER
c          The leading dimension of the array C. LDC >= max(1,M).
c
c  WORK    (workspace) COMPLEX*16 array, dimension
c                         (N) if SIDE = 'L'
c                      or (M) if SIDE = 'R'
c
c  =====================================================================
c
c     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
c     ..
c     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZGERC
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. Executable Statements ..
c
      IF( LSAME( SIDE, 'L' ) ) THEN
c
c        Form  H * C
c
         IF( TAU.NE.ZERO ) THEN
c
c           w := C' * v
c
            CALL ZGEMV( 'Conjugate transpose', M, N, ONE, C, LDC, V,
     $                  INCV, ZERO, WORK, 1 )
c
c           C := C - v * w'
c
            CALL ZGERC( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
c
c        Form  C * H
c
         IF( TAU.NE.ZERO ) THEN
c
c           w := C * v
c
            CALL ZGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV,
     $                  ZERO, WORK, 1 )
c
c           C := C - w * v'
c
            CALL ZGERC( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
c
c     End of ZLARF
c
      END
      subroutine  zscal(n,za,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex za,zx(*)
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = za*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = za*zx(i)
   30 continue
      return
      end
      SUBROUTINE ZUNG2L( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, *), TAU(*), WORK(*)
*     ..
*
*  Purpose
*  =======
*
*  ZUNG2L generates an m by n complex matrix Q with orthonormal columns,
*  which is defined as the last n columns of a product of k elementary
*  reflectors of order m
*
*        Q  =  H(k) . . . H(2) H(1)
*
*  as returned by ZGEQLF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the (n-k+i)-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by ZGEQLF in the last k columns of its array
*          argument A.
*          On exit, the m-by-n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGEQLF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZSCAL
*     ..
*     .. Intrinsic Functions ..
* PAN INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNG2L', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns 1:n-k to columns of the unit matrix
*
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = 1, K
         II = N - K + I
*
*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
*
         A( M-N+II, II ) = ONE
         CALL ZLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A,
     $               LDA, WORK )
         CALL ZSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
*
*        Set A(m-k+i+1:m,n-k+i) to zero
*
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of ZUNG2L
*
      END
      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
*     ..
*
*  Purpose
*  =======
*
*  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
*  eigenvector for RT1, giving the decomposition
*
*     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
*     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  B       (input) DOUBLE PRECISION
*          The (1,2) element and the conjugate of the (2,1) element of
*          the 2-by-2 matrix.
*
*  C       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  RT1     (output) DOUBLE PRECISION
*          The eigenvalue of larger absolute value.
*
*  RT2     (output) DOUBLE PRECISION
*          The eigenvalue of smaller absolute value.
*
*  CS1     (output) DOUBLE PRECISION
*  SN1     (output) DOUBLE PRECISION
*          The vector (CS1, SN1) is a unit right eigenvector for RT1.
*
*  Further Details
*  ===============
*
*  RT1 is accurate to a few ulps barring over/underflow.
*
*  RT2 may be inaccurate if there is massive cancellation in the
*  determinant A*C-B*B; higher precision or correctly rounded or
*  correctly truncated arithmetic would be needed to compute RT2
*  accurately in all cases.
*
*  CS1 and SN1 are accurate to a few ulps barring over/underflow.
*
*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*  Underflow is harmless if the input data is 0 or exceeds
*     underflow_threshold / macheps.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN
*     ..
*     .. Intrinsic Functions ..
*PAN   INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
*
*     Compute the eigenvector
*
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
*
*     End of DLAEV2
*
      END
      SUBROUTINE DLASRT( ID, N, D, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D(*)
*     ..
*
*  Purpose
*  =======
*
*  Sort the numbers in D in increasing order (if ID = 'I') or
*  in decreasing order (if ID = 'D' ).
*
*  Use Quick Sort, reverting to Insertion sort on arrays of
*  size <= 20. Dimension of STACK limits N to about 2**32.
*
*  Arguments
*  =========
*
*  ID      (input) CHARACTER*1
*          = 'I': sort D in increasing order;
*          = 'D': sort D in decreasing order.
*
*  N       (input) INTEGER
*          The length of the array D.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the array to be sorted.
*          On exit, D has been sorted into increasing order
*          (D(1) <= ... <= D(N) ) or into decreasing order
*          (D(1) >= ... >= D(N) ), depending on ID.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
*     ..
*     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
*     ..
*     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input paramters.
*
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
*
*        Do Insertion sort on D( START:ENDD )
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
*
         ELSE
*
*           Sort into increasing order
*
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
*
         END IF
*
      ELSE IF( ENDD-START.GT.SELECT ) THEN
*
*        Partition D( START:ENDD ) and stack parts, largest one first
*
*        Choose partition entry as median of 3
*
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX )
     $         GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX )
     $         GO TO 80
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
*
*           Sort into increasing order
*
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX )
     $         GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX )
     $         GO TO 110
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 )
     $   GO TO 10
      RETURN
*
*     End of DLASRT
*
      END
      SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA,*)
*     ..
*
*  Purpose
*  =======
*
*  ZLASET initializes a 2-D array A to BETA on the diagonal and
*  ALPHA on the offdiagonals.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be set.
*          = 'U':      Upper triangular part is set. The lower triangle
*                      is unchanged.
*          = 'L':      Lower triangular part is set. The upper triangle
*                      is unchanged.
*          Otherwise:  All of the matrix A is set.
*
*  M       (input) INTEGER
*          On entry, M specifies the number of rows of A.
*
*  N       (input) INTEGER
*          On entry, N specifies the number of columns of A.
*
*  ALPHA   (input) COMPLEX*16
*          All the offdiagonal array elements are set to ALPHA.
*
*  BETA    (input) COMPLEX*16
*          All the diagonal array elements are set to BETA.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
*                   A(i,i) = BETA , 1 <= i <= min(m,n)
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
*PAN     INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Set the diagonal to BETA and the strictly upper triangular
*        part of the array to ALPHA.
*
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
         DO 30 I = 1, MIN( N, M )
            A( I, I ) = BETA
   30    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Set the diagonal to BETA and the strictly lower triangular
*        part of the array to ALPHA.
*
         DO 50 J = 1, MIN( M, N )
            DO 40 I = J + 1, M
               A( I, J ) = ALPHA
   40       CONTINUE
   50    CONTINUE
         DO 60 I = 1, MIN( N, M )
            A( I, I ) = BETA
   60    CONTINUE
*
      ELSE
*
*        Set the array to BETA on the diagonal and ALPHA on the
*        offdiagonal.
*
         DO 80 J = 1, N
            DO 70 I = 1, M
               A( I, J ) = ALPHA
   70       CONTINUE
   80    CONTINUE
         DO 90 I = 1, MIN( M, N )
            A( I, I ) = BETA
   90    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASET
*
      END
      SUBROUTINE ZSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D(*), E(*), WORK(*)
      COMPLEX*16         Z( LDZ,*)
*     ..
*
*  Purpose
*  =======
*
*  ZSTEQR computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the implicit QL or QR method.
*  The eigenvectors of a full or band complex Hermitian matrix can also
*  be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this
*  matrix to tridiagonal form.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'V':  Compute eigenvalues and eigenvectors of the original
*                  Hermitian matrix.  On entry, Z must contain the
*                  unitary matrix used to reduce the original matrix
*                  to tridiagonal form.
*          = 'I':  Compute eigenvalues and eigenvectors of the
*                  tridiagonal matrix.  Z is initialized to the identity
*                  matrix.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
*          On entry, if  COMPZ = 'V', then Z contains the unitary
*          matrix used in the reduction to tridiagonal form.
*          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
*          orthonormal eigenvectors of the original Hermitian matrix,
*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*          of the symmetric tridiagonal matrix.
*          If COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          eigenvectors are desired, then  LDZ >= max(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
*          If COMPZ = 'N', then WORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm has failed to find all the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero; on exit, D
*                and E contain the elements of a symmetric tridiagonal
*                matrix which is unitarily similar to the original
*                matrix.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ),
     $                   CONE = ( 1.0D0, 0.0D0 ) )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,
     $                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,
     $                   NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,
     $                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASRT, XERBLA,
     $                   ZLASET, ZLASR, ZSWAP
*     ..
*     .. Intrinsic Functions ..
*PAN     INTRINSIC          ABS, MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSTEQR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 )
     $      Z( 1, 1 ) = CONE
         RETURN
      END IF
*
*     Determine the unit roundoff and over/underflow thresholds.
*
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues and eigenvectors of the tridiagonal
*     matrix.
*
      IF( ICOMPZ.EQ.2 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )
*
      NMAXIT = N*MAXIT
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GT.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+
     $             SAFMIN )GO TO 60
   50       CONTINUE
         END IF
*
         M = LEND
*
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80
*
*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL ZLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
*
   70    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL ZLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
*
*        Eigenvalue found.
*
   80    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 140
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+
     $             SAFMIN )GO TO 110
  100       CONTINUE
         END IF
*
         M = LEND
*
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130
*
*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL ZLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
*
  120    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL ZLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
*
*        Eigenvalue found.
*
  130    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 140
*
      END IF
*
*     Undo scaling if necessary
*
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      END IF
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.EQ.NMAXIT ) THEN
         DO 150 I = 1, N - 1
            IF( E( I ).NE.ZERO )
     $         INFO = INFO + 1
  150    CONTINUE
         RETURN
      END IF
      GO TO 10
*
*     Order eigenvalues and eigenvectors.
*
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
*
*        Use Quick Sort
*
         CALL DLASRT( 'I', N, D, INFO )
*
      ELSE
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
      RETURN
*
*     End of ZSTEQR
*
      END
      SUBROUTINE ZUNG2R( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA,*), TAU(*), WORK(*)
*     ..
*
*  Purpose
*  =======
*
*  ZUNG2R generates an m by n complex matrix Q with orthonormal columns,
*  which is defined as the first n columns of a product of k elementary
*  reflectors of order m
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by ZGEQRF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by ZGEQRF in the first k columns of its array
*          argument A.
*          On exit, the m by n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGEQRF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZSCAL
*     ..
*     .. Intrinsic Functions ..
*PAN     INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNG2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns k+1:n to columns of the unit matrix
*
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = K, 1, -1
*
*        Apply H(i) to A(i:m,i:n) from the left
*
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL ZLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )
     $      CALL ZSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
*
*        Set A(1:i-1,i) to zero
*
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of ZUNG2R
*
      END
      SUBROUTINE DLARTG( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
*     ..
*
*  Purpose
*  =======
*
*  DLARTG generate a plane rotation so that
*
*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a slower, more accurate version of the BLAS1 routine DROTG,
*  with the following other differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*        floating point operations (saves work in DBDSQR when
*        there are zeros on the diagonal).
*
*  If F exceeds G in magnitude, CS will be positive.
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The first component of vector to be rotated.
*
*  G       (input) DOUBLE PRECISION
*          The second component of vector to be rotated.
*
*  CS      (output) DOUBLE PRECISION
*          The cosine of the rotation.
*
*  SN      (output) DOUBLE PRECISION
*          The sine of the rotation.
*
*  R       (output) DOUBLE PRECISION
*          The nonzero component of the rotated vector.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
*PAN     INTRINSIC          ABS, INT, LOG, MAX, SQRT
*     ..
*     .. Save statement ..
      SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
      END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 )
     $         GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 )
     $         GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
*
*     End of DLARTG
*
      END
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (version 2.0) --
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
      SUBROUTINE ZHPTRD( UPLO, N, AP, D, E, TAU, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D(*), E(* )
      COMPLEX*16         AP(*), TAU(*)
*     ..
*
*  Purpose
*  =======
*
*  ZHPTRD reduces a complex Hermitian matrix A stored in packed form to
*  real symmetric tridiagonal form T by a unitary similarity
*  transformation: Q**H * A * Q = T.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)
*          On entry, the upper or lower triangle of the Hermitian matrix
*          A, packed columnwise in a linear array.  The j-th column of A
*          is stored in the array AP as follows:
*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*          of A are overwritten by the corresponding elements of the
*          tridiagonal matrix T, and the elements above the first
*          superdiagonal, with the array TAU, represent the unitary
*          matrix Q as a product of elementary reflectors; if UPLO
*          = 'L', the diagonal and first subdiagonal of A are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements below the first subdiagonal, with
*          the array TAU, represent the unitary matrix Q as a product
*          of elementary reflectors. See Further Details.
*
*  D       (output) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i).
*
*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*
*  TAU     (output) COMPLEX*16 array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in AP,
*  overwriting A(1:i-1,i+1), and tau is stored in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,
*  overwriting A(i+2:n,i), and tau is stored in TAU(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO, HALF
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   HALF = ( 0.5D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, I1, I1I1, II
      COMPLEX*16         ALPHA, TAUI
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZHPMV, ZHPR2, ZLARFG
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHPTRD', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A.
*        I1 is the index in AP of A(1,I+1).
*
         I1 = N*( N-1 ) / 2 + 1
         AP( I1+N-1 ) = DBLE( AP( I1+N-1 ) )
         DO 10 I = N - 1, 1, -1
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(1:i-1,i+1)
*
            ALPHA = AP( I1+I-1 )
            CALL ZLARFG( I, ALPHA, AP( I1 ), 1, TAUI )
            E( I ) = ALPHA
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(1:i,1:i)
*
               AP( I1+I-1 ) = ONE
*
*              Compute  y := tau * A * v  storing y in TAU(1:i)
*
               CALL ZHPMV( UPLO, I, TAUI, AP, AP( I1 ), 1, ZERO, TAU,
     $                     1 )
*
*              Compute  w := y - 1/2 * tau * (y'*v) * v
*
               ALPHA = -HALF*TAUI*ZDOTC( I, TAU, 1, AP( I1 ), 1 )
               CALL ZAXPY( I, ALPHA, AP( I1 ), 1, TAU, 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w' - w * v'
*
               CALL ZHPR2( UPLO, I, -ONE, AP( I1 ), 1, TAU, 1, AP )
*
            END IF
            AP( I1+I-1 ) = E( I )
            D( I+1 ) = AP( I1+I )
            TAU( I ) = TAUI
            I1 = I1 - I
   10    CONTINUE
         D( 1 ) = AP( 1 )
      ELSE
*
*        Reduce the lower triangle of A. II is the index in AP of
*        A(i,i) and I1I1 is the index of A(i+1,i+1).
*
         II = 1
         AP( 1 ) = DBLE( AP( 1 ) )
         DO 20 I = 1, N - 1
            I1I1 = II + N - I + 1
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(i+2:n,i)
*
            ALPHA = AP( II+1 )
            CALL ZLARFG( N-I, ALPHA, AP( II+2 ), 1, TAUI )
            E( I ) = ALPHA
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(i+1:n,i+1:n)
*
               AP( II+1 ) = ONE
*
*              Compute  y := tau * A * v  storing y in TAU(i:n-1)
*
               CALL ZHPMV( UPLO, N-I, TAUI, AP( I1I1 ), AP( II+1 ), 1,
     $                     ZERO, TAU( I ), 1 )
*
*              Compute  w := y - 1/2 * tau * (y'*v) * v
*
               ALPHA = -HALF*TAUI*ZDOTC( N-I, TAU( I ), 1, AP( II+1 ),
     $                 1 )
               CALL ZAXPY( N-I, ALPHA, AP( II+1 ), 1, TAU( I ), 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w' - w * v'
*
               CALL ZHPR2( UPLO, N-I, -ONE, AP( II+1 ), 1, TAU( I ), 1,
     $                     AP( I1I1 ) )
*
            END IF
            AP( II+1 ) = E( I )
            D( I ) = AP( II )
            TAU( I ) = TAUI
            II = I1I1
   20    CONTINUE
         D( N ) = AP( II )
      END IF
*
      RETURN
*
*     End of ZHPTRD
*
      END
      SUBROUTINE ZLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C(*), S(*)
      COMPLEX*16         A( LDA, *)
*     ..
*
*  Purpose
*  =======
*
*  ZLASR   performs the transformation
*
*     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
*
*     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
*
*  where A is an m by n complex matrix and P is an orthogonal matrix,
*  consisting of a sequence of plane rotations determined by the
*  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
*  and z = n when SIDE = 'R' or 'r' ):
*
*  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
*
*     P = P( z - 1 )*...*P( 2 )*P( 1 ),
*
*  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
*
*     P = P( 1 )*P( 2 )*...*P( z - 1 ),
*
*  where  P( k ) is a plane rotation matrix for the following planes:
*
*     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
*        the plane ( k, k + 1 )
*
*     when  PIVOT = 'T' or 't'  ( Top pivot ),
*        the plane ( 1, k + 1 )
*
*     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
*        the plane ( k, z )
*
*  c( k ) and s( k )  must contain the  cosine and sine that define the
*  matrix  P( k ).  The two by two plane rotation part of the matrix
*  P( k ), R( k ), is assumed to be of the form
*
*     R( k ) = (  c( k )  s( k ) ).
*              ( -s( k )  c( k ) )
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          Specifies whether the plane rotation matrix P is applied to
*          A on the left or the right.
*          = 'L':  Left, compute A := P*A
*          = 'R':  Right, compute A:= A*P'
*
*  DIRECT  (input) CHARACTER*1
*          Specifies whether P is a forward or backward sequence of
*          plane rotations.
*          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
*          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
*
*  PIVOT   (input) CHARACTER*1
*          Specifies the plane for which P(k) is a plane rotation
*          matrix.
*          = 'V':  Variable pivot, the plane (k,k+1)
*          = 'T':  Top pivot, the plane (1,k+1)
*          = 'B':  Bottom pivot, the plane (k,z)
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  If m <= 1, an immediate
*          return is effected.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  If n <= 1, an
*          immediate return is effected.
*
*  C, S    (input) DOUBLE PRECISION arrays, dimension
*                  (M-1) if SIDE = 'L'
*                  (N-1) if SIDE = 'R'
*          c(k) and s(k) contain the cosine and sine that define the
*          matrix P(k).  The two by two plane rotation part of the
*          matrix P(k), R(k), is assumed to be of the form
*          R( k ) = (  c( k )  s( k ) ).
*                   ( -s( k )  c( k ) )
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          The m by n matrix A.  On exit, A is overwritten by P*A if
*          SIDE = 'R' or by A*P' if SIDE = 'L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP
      COMPLEX*16         TEMP
*     ..
*     .. Intrinsic Functions ..
*PAN     INTRINSIC          MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLASR ', INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  P * A
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        Form A * P'
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of ZLASR
*
      END
      subroutine  zswap (n,zx,incx,zy,incy)
c
c     interchanges two vectors.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(N),zy(N),ztemp
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = zx(ix)
        zx(ix) = zy(iy)
        zy(iy) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
   20 do 30 i = 1,n
        ztemp = zx(i)
        zx(i) = zy(i)
        zy(i) = ztemp
   30 continue
      return
      end
      SUBROUTINE ZUPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDQ, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         AP(*), Q( LDQ,*), TAU(*), WORK(*)
*     ..
*
*  Purpose
*  =======
*
*  ZUPGTR generates a complex unitary matrix Q which is defined as the
*  product of n-1 elementary reflectors H(i) of order n, as returned by
*  ZHPTRD using packed storage:
*
*  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
*
*  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U': Upper triangular packed storage used in previous
*                 call to ZHPTRD;
*          = 'L': Lower triangular packed storage used in previous
*                 call to ZHPTRD.
*
*  N       (input) INTEGER
*          The order of the matrix Q. N >= 0.
*
*  AP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)
*          The vectors which define the elementary reflectors, as
*          returned by ZHPTRD.
*
*  TAU     (input) COMPLEX*16 array, dimension (N-1)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZHPTRD.
*
*  Q       (output) COMPLEX*16 array, dimension (LDQ,N)
*          The N-by-N unitary matrix Q.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= max(1,N).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N-1)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IINFO, IJ, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNG2L, ZUNG2R
*     ..
*     .. Intrinsic Functions ..
*PAN     INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUPGTR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Q was determined by a call to ZHPTRD with UPLO = 'U'
*
*        Unpack the vectors which define the elementary reflectors and
*        set the last row and column of Q equal to those of the unit
*        matrix
*
         IJ = 2
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
   10       CONTINUE
            IJ = IJ + 2
            Q( N, J ) = CZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            Q( I, N ) = CZERO
   30    CONTINUE
         Q( N, N ) = CONE
*
*        Generate Q(1:n-1,1:n-1)
*
         CALL ZUNG2L( N-1, N-1, N-1, Q, LDQ, TAU, WORK, IINFO )
*
      ELSE
*
*        Q was determined by a call to ZHPTRD with UPLO = 'L'.
*
*        Unpack the vectors which define the elementary reflectors and
*        set the first row and column of Q equal to those of the unit
*        matrix
*
         Q( 1, 1 ) = CONE
         DO 40 I = 2, N
            Q( I, 1 ) = CZERO
   40    CONTINUE
         IJ = 3
         DO 60 J = 2, N
            Q( 1, J ) = CZERO
            DO 50 I = J + 1, N
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
   50       CONTINUE
            IJ = IJ + 2
   60    CONTINUE
         IF( N.GT.1 ) THEN
*
*           Generate Q(2:n,2:n)
*
            CALL ZUNG2R( N-1, N-1, N-1, Q( 2, 2 ), LDQ, TAU, WORK,
     $                   IINFO )
         END IF
      END IF
      RETURN
*
*     End of ZUPGTR
*
      END
      SUBROUTINE ZGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
*PAN     INTRINSIC          ABS
      COMPLEX*16         ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*), X(*), 
     1            Y(*)
*     ..
*
*  Purpose
*  =======
*
*  ZGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
*
*     y := alpha*conjg( A' )*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
      LOGICAL            NOCONJ
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
*     INTRINSIC          DCONJG, MAX
*     INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
      NOCONJ = LSAME( TRANS, 'T' )
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 110, J = 1, N
               TEMP = ZERO
               IF( NOCONJ )THEN
                  DO 90, I = 1, M
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
               ELSE
                  DO 100, I = 1, M
                     TEMP = TEMP + DCONJG( A( I, J ) )*X( I )
  100             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  110       CONTINUE
         ELSE
            DO 140, J = 1, N
               TEMP = ZERO
               IX   = KX
               IF( NOCONJ )THEN
                  DO 120, I = 1, M
                     TEMP = TEMP + A( I, J )*X( IX )
                     IX   = IX   + INCX
  120             CONTINUE
               ELSE
                  DO 130, I = 1, M
                     TEMP = TEMP + DCONJG( A( I, J ) )*X( IX )
                     IX   = IX   + INCX
  130             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  140       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZGEMV .
*
      END
      SUBROUTINE ZGERC ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
*PAN     INTRINSIC          ABS
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*), X(*), Y(*)
*     ..
*
*  Purpose
*  =======
*
*  ZGERC  performs the rank 1 operation
*
*     A := alpha*x*conjg( y' ) + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
*PAN     INTRINSIC          DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
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
         CALL XERBLA( 'ZGERC ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*DCONJG( Y( JY ) )
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
               TEMP = ALPHA*DCONJG( Y( JY ) )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of ZGERC .
*
      END
      SUBROUTINE ZHPMV ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
*     .. Scalar Arguments ..
*PAN     INTRINSIC          ABS
      COMPLEX*16         ALPHA, BETA
      INTEGER            INCX, INCY, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      COMPLEX*16 AP(*),X(*), 
     1            Y(*)
*     ..
*
*  Purpose
*  =======
*
*  ZHPMV  performs the matrix-vector operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n hermitian matrix, supplied in packed form.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the matrix A is supplied in the packed
*           array AP as follows:
*
*              UPLO = 'U' or 'u'   The upper triangular part of A is
*                                  supplied in AP.
*
*              UPLO = 'L' or 'l'   The lower triangular part of A is
*                                  supplied in AP.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  AP     - COMPLEX*16       array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with UPLO = 'U' or 'u', the array AP must
*           contain the upper triangular part of the hermitian matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
*           and a( 2, 2 ) respectively, and so on.
*           Before entry with UPLO = 'L' or 'l', the array AP must
*           contain the lower triangular part of the hermitian matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
*           and a( 3, 1 ) respectively, and so on.
*           Note that the imaginary parts of the diagonal elements need
*           not be set and are assumed to be zero.
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y. On exit, Y is overwritten by the updated
*           vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     .. Local Scalars ..
      COMPLEX*16         TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
*PAN     INTRINSIC          DCONJG, DBLE
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 6
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZHPMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set up the start points in  X  and  Y.
*
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of the array AP
*     are accessed sequentially with one pass through AP.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      KK = 1
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  y  when AP contains the upper triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               K     = KK
               DO 50, I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*AP( K )
                  TEMP2  = TEMP2  + DCONJG( AP( K ) )*X( I )
                  K      = K      + 1
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*DBLE( AP( KK + J - 1 ) )
     $                         + ALPHA*TEMP2
               KK     = KK     + J
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO 70, K = KK, KK + J - 2
                  Y( IY ) = Y( IY ) + TEMP1*AP( K )
                  TEMP2   = TEMP2   + DCONJG( AP( K ) )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*DBLE( AP( KK + J - 1 ) )
     $                           + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
               KK      = KK      + J
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y  when AP contains the lower triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J ) + TEMP1*DBLE( AP( KK ) )
               K      = KK     + 1
               DO 90, I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*AP( K )
                  TEMP2  = TEMP2  + DCONJG( AP( K ) )*X( I )
                  K      = K      + 1
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
               KK     = KK     + ( N - J + 1 )
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY ) + TEMP1*DBLE( AP( KK ) )
               IX      = JX
               IY      = JY
               DO 110, K = KK + 1, KK + N - J
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*AP( K )
                  TEMP2   = TEMP2   + DCONJG( AP( K ) )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
               KK      = KK      + ( N - J + 1 )
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZHPMV .
*
      END
      subroutine zaxpy(n,za,zx,incx,zy,incy)
c
c     constant times a vector plus a vector.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(N),zy(N),za
      integer i,incx,incy,ix,iy,n
      double precision dcabs1
      if(n.le.0)return
      if (dcabs1(za) .eq. 0.0d0) return
      if (incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zy(iy) + za*zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zy(i) = zy(i) + za*zx(i)
   30 continue
      return
      end
      SUBROUTINE ZHPR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
*     .. Scalar Arguments ..
*PAN  INTRINSIC          ABS
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      COMPLEX*16 AP(*),X(*),Y(*)
*     ..
*
*  Purpose
*  =======
*
*  ZHPR2  performs the hermitian rank 2 operation
*
*     A := alpha*x*conjg( y' ) + conjg( alpha )*y*conjg( x' ) + A,
*
*  where alpha is a scalar, x and y are n element vectors and A is an
*  n by n hermitian matrix, supplied in packed form.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the matrix A is supplied in the packed
*           array AP as follows:
*
*              UPLO = 'U' or 'u'   The upper triangular part of A is
*                                  supplied in AP.
*
*              UPLO = 'L' or 'l'   The lower triangular part of A is
*                                  supplied in AP.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  AP     - COMPLEX*16       array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with  UPLO = 'U' or 'u', the array AP must
*           contain the upper triangular part of the hermitian matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
*           and a( 2, 2 ) respectively, and so on. On exit, the array
*           AP is overwritten by the upper triangular part of the
*           updated matrix.
*           Before entry with UPLO = 'L' or 'l', the array AP must
*           contain the lower triangular part of the hermitian matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
*           and a( 3, 1 ) respectively, and so on. On exit, the array
*           AP is overwritten by the lower triangular part of the
*           updated matrix.
*           Note that the imaginary parts of the diagonal elements need
*           not be set, they are assumed to be zero, and on exit they
*           are set to zero.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     .. Local Scalars ..
      COMPLEX*16         TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
*PAN     INTRINSIC          DCONJG, DBLE
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZHPR2 ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Set up the start points in X and Y if the increments are not both
*     unity.
*
      IF( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
         JX = KX
         JY = KY
      END IF
*
*     Start the operations. In this version the elements of the array AP
*     are accessed sequentially with one pass through AP.
*
      KK = 1
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  A  when upper triangle is stored in AP.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*DCONJG( Y( J ) )
                  TEMP2 = DCONJG( ALPHA*X( J ) )
                  K     = KK
                  DO 10, I = 1, J - 1
                     AP( K ) = AP( K ) + X( I )*TEMP1 + Y( I )*TEMP2
                     K       = K       + 1
   10             CONTINUE
                  AP( KK + J - 1 ) = DBLE( AP( KK + J - 1 ) ) +
     $                               DBLE( X( J )*TEMP1 + Y( J )*TEMP2 )
               ELSE
                  AP( KK + J - 1 ) = DBLE( AP( KK + J - 1 ) )
               END IF
               KK = KK + J
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*DCONJG( Y( JY ) )
                  TEMP2 = DCONJG( ALPHA*X( JX ) )
                  IX    = KX
                  IY    = KY
                  DO 30, K = KK, KK + J - 2
                     AP( K ) = AP( K ) + X( IX )*TEMP1 + Y( IY )*TEMP2
                     IX      = IX      + INCX
                     IY      = IY      + INCY
   30             CONTINUE
                  AP( KK + J - 1 ) = DBLE( AP( KK + J - 1 ) ) +
     $                               DBLE( X( JX )*TEMP1 +
     $                                     Y( JY )*TEMP2 )
               ELSE
                  AP( KK + J - 1 ) = DBLE( AP( KK + J - 1 ) )
               END IF
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + J
   40       CONTINUE
         END IF
      ELSE
*
*        Form  A  when lower triangle is stored in AP.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1   = ALPHA*DCONJG( Y( J ) )
                  TEMP2   = DCONJG( ALPHA*X( J ) )
                  AP( KK ) = DBLE( AP( KK ) ) +
     $                       DBLE( X( J )*TEMP1 + Y( J )*TEMP2 )
                  K        = KK               + 1
                  DO 50, I = J + 1, N
                     AP( K ) = AP( K ) + X( I )*TEMP1 + Y( I )*TEMP2
                     K       = K       + 1
   50             CONTINUE
               ELSE
                  AP( KK ) = DBLE( AP( KK ) )
               END IF
               KK = KK + N - J + 1
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1    = ALPHA*DCONJG( Y( JY ) )
                  TEMP2    = DCONJG( ALPHA*X( JX ) )
                  AP( KK ) = DBLE( AP( KK ) ) +
     $                       DBLE( X( JX )*TEMP1 + Y( JY )*TEMP2 )
                  IX       = JX
                  IY       = JY
                  DO 70, K = KK + 1, KK + N - J
                     IX      = IX      + INCX
                     IY      = IY      + INCY
                     AP( K ) = AP( K ) + X( IX )*TEMP1 + Y( IY )*TEMP2
   70             CONTINUE
               ELSE
                  AP( KK ) = DBLE( AP( KK ) )
               END IF
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + N - J + 1
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZHPR2 .
*
      END
      double precision function dcabs1(z)
      double complex z,zz
      double precision t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = dabs(t(1)) + dabs(t(2))
      return
      end
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      DLAMCH = RMACH
      RETURN
*
*     End of DLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
*PAN  INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A, B    (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
*PAN     INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5
*
      END
      DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  DLANST  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real symmetric tridiagonal matrix A.
*
*  Description
*  ===========
*
*  DLANST returns the value
*
*     DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANST as described
*          above.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, DLANST is
*          set to zero.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of A.
*
*  E       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) sub-diagonal or super-diagonal elements of A.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, SCALE, SUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. Intrinsic Functions ..
*PAN      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         ANORM = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         ANORM = ABS( D( N ) )
         DO 10 I = 1, N - 1
            ANORM = MAX( ANORM, ABS( D( I ) ) )
            ANORM = MAX( ANORM, ABS( E( I ) ) )
   10    CONTINUE
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.
     $         LSAME( NORM, 'I' ) ) THEN
*
*        Find norm1(A).
*
         IF( N.EQ.1 ) THEN
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ),
     $              ABS( E( N-1 ) )+ABS( D( N ) ) )
            DO 20 I = 2, N - 1
               ANORM = MAX( ANORM, ABS( D( I ) )+ABS( E( I ) )+
     $                 ABS( E( I-1 ) ) )
   20       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( N.GT.1 ) THEN
            CALL DLASSQ( N-1, E, 1, SCALE, SUM )
            SUM = 2*SUM
         END IF
         CALL DLASSQ( N, D, 1, SCALE, SUM )
         ANORM = SCALE*SQRT( SUM )
      END IF
*
      DLANST = ANORM
      RETURN
*
*     End of DLANST
*
      END
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
*     ..
*
*  Purpose
*  =======
*
*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*          X and Y specify the values x and y.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
*     ..
*     .. Intrinsic Functions ..
*PAN      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
*
*     End of DLAPY2
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     January 31, 1994
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
      double complex function zdotc(n,zx,incx,zy,incy)
c
c     forms the dot product of a vector.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*),zy(*),ztemp
      integer i,incx,incy,ix,iy,n
      ztemp = (0.0d0,0.0d0)
      zdotc = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = ztemp + dconjg(zx(ix))*zy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      zdotc = ztemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = ztemp + dconjg(zx(i))*zy(i)
   30 continue
      zdotc = ztemp
      return
      end
      DOUBLE PRECISION FUNCTION ZLANHP( NORM, UPLO, N, AP, WORK )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         AP( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLANHP  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  complex hermitian matrix A,  supplied in packed form.
*
*  Description
*  ===========
*
*  ZLANHP returns the value
*
*     ZLANHP = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in ZLANHP as described
*          above.
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          hermitian matrix A is supplied.
*          = 'U':  Upper triangular part of A is supplied
*          = 'L':  Lower triangular part of A is supplied
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, ZLANHP is
*          set to zero.
*
*  AP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)
*          The upper or lower triangle of the hermitian matrix A, packed
*          columnwise in a linear array.  The j-th column of A is stored
*          in the array AP as follows:
*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
*          Note that the  imaginary parts of the diagonal elements need
*          not be set and are assumed to be zero.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
*          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
*          WORK is not referenced.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASSQ
*     ..
*     .. Intrinsic Functions ..
*PAN      INTRINSIC          ABS, DBLE, MAX, SQRT
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            K = 0
            DO 20 J = 1, N
               DO 10 I = K + 1, K + J - 1
                  VALUE = MAX( VALUE, ABS( AP( I ) ) )
   10          CONTINUE
               K = K + J
               VALUE = MAX( VALUE, ABS( DBLE( AP( K ) ) ) )
   20       CONTINUE
         ELSE
            K = 1
            DO 40 J = 1, N
               VALUE = MAX( VALUE, ABS( DBLE( AP( K ) ) ) )
               DO 30 I = K + 1, K + N - J
                  VALUE = MAX( VALUE, ABS( AP( I ) ) )
   30          CONTINUE
               K = K + N - J + 1
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.
     $         ( NORM.EQ.'1' ) ) THEN
*
*        Find normI(A) ( = norm1(A), since A is hermitian).
*
         VALUE = ZERO
         K = 1
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
   50          CONTINUE
               WORK( J ) = SUM + ABS( DBLE( AP( K ) ) )
               K = K + 1
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( DBLE( AP( K ) ) )
               K = K + 1
               DO 90 I = J + 1, N
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         K = 2
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL ZLASSQ( J-1, AP( K ), 1, SCALE, SUM )
               K = K + J
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL ZLASSQ( N-J, AP( K ), 1, SCALE, SUM )
               K = K + N - J + 1
  120       CONTINUE
         END IF
         SUM = 2*SUM
         K = 1
         DO 130 I = 1, N
            IF( DBLE( AP( K ) ).NE.ZERO ) THEN
               ABSA = ABS( DBLE( AP( K ) ) )
               IF( SCALE.LT.ABSA ) THEN
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               ELSE
                  SUM = SUM + ( ABSA / SCALE )**2
               END IF
            END IF
            IF( LSAME( UPLO, 'U' ) ) THEN
               K = K + I + 1
            ELSE
               K = K + N - I + 1
            END IF
  130    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      ZLANHP = VALUE
      RETURN
*
*     End of ZLANHP
*
      END
      SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX*16         ALPHA, TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARFG generates a complex elementary reflector H of order n, such
*  that
*
*        H' * ( alpha ) = ( beta ),   H' * H = I.
*             (   x   )   (   0  )
*
*  where alpha and beta are scalars, with beta real, and x is an
*  (n-1)-element complex vector. H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a complex scalar and v is a complex (n-1)-element
*  vector. Note that H is not hermitian.
*
*  If the elements of x are all zero and alpha is real, then tau = 0
*  and H is taken to be the unit matrix.
*
*  Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) COMPLEX*16
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) COMPLEX*16 array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  TAU     (output) COMPLEX*16
*          The value tau.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY3, DZNRM2
      COMPLEX*16         ZLADIV
      EXTERNAL           DLAMCH, DLAPY3, DZNRM2, ZLADIV
*     ..
*     .. Intrinsic Functions ..
* BCPAN INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, SIGN
        INTRINSIC               DBLE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZDSCAL, ZSCAL
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = DZNRM2( N-1, X, INCX )
      ALPHR = DBLE( ALPHA )
      ALPHI = DIMAG( ALPHA )
*
      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        general case
*
         BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         RSAFMN = ONE / SAFMIN
*
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL ZDSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = DZNRM2( N-1, X, INCX )
            ALPHA = DCMPLX( ALPHR, ALPHI )
            BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
            TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
            ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA )
            CALL ZSCAL( N-1, ALPHA, X, INCX )
*
*           If ALPHA is subnormal, it may lose relative accuracy
*
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
            ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA )
            CALL ZSCAL( N-1, ALPHA, X, INCX )
            ALPHA = BETA
         END IF
      END IF
*
      RETURN
*
*     End of ZLARFG
*
      END
      DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y, Z
*     ..
*
*  Purpose
*  =======
*
*  DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
*  unnecessary overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*  Z       (input) DOUBLE PRECISION
*          X, Y and Z specify the values x, y and z.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, ZABS
*     ..
*     .. Intrinsic Functions ..
*PAN      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO ) THEN
         DLAPY3 = ZERO
      ELSE
         DLAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+
     $            ( ZABS / W )**2 )
      END IF
      RETURN
*
*     End of DLAPY3
*
      END
      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASSQ  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in SCALE and SUMSQ and
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector x.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) DOUBLE PRECISION
*          The vector for which a scaled sum of squares is computed.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
*     ..
*     .. Intrinsic Functions ..
*PAN      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of DLASSQ
*
      END
      DOUBLE PRECISION FUNCTION DZNRM2( N, X, INCX )
*     .. Scalar Arguments ..
      INTEGER                           INCX, N
*     .. Array Arguments ..
      COMPLEX*16                        X( * )
*     ..
*
*  DZNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DZNRM2 := sqrt( conjg( x' )*x )
*
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to ZLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
*     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      INTEGER               IX
      DOUBLE PRECISION      NORM, SCALE, SSQ, TEMP
*     .. Intrinsic Functions ..
*PAN      INTRINSIC             ABS, DIMAG, DBLE, SQRT
      INTRINSIC             DBLE
*     ..
*     .. Executable Statements ..
      IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE
         SCALE = ZERO
         SSQ   = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
*
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( DBLE( X( IX ) ).NE.ZERO )THEN
               TEMP = ABS( DBLE( X( IX ) ) )
               IF( SCALE.LT.TEMP )THEN
                  SSQ   = ONE   + SSQ*( SCALE/TEMP )**2
                  SCALE = TEMP
               ELSE
                  SSQ   = SSQ   +     ( TEMP/SCALE )**2
               END IF
            END IF
            IF( DIMAG( X( IX ) ).NE.ZERO )THEN
               TEMP = ABS( DIMAG( X( IX ) ) )
               IF( SCALE.LT.TEMP )THEN
                  SSQ   = ONE   + SSQ*( SCALE/TEMP )**2
                  SCALE = TEMP
               ELSE
                  SSQ   = SSQ   +     ( TEMP/SCALE )**2
               END IF
            END IF
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      END IF
*
      DZNRM2 = NORM
      RETURN
*
*     End of DZNRM2.
*
      END
      DOUBLE COMPLEX   FUNCTION ZLADIV( X, Y )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      COMPLEX*16         X, Y
*     ..
*
*  Purpose
*  =======
*
*  ZLADIV := X / Y, where X and Y are complex.  The computation of X / Y
*  will not overflow on an intermediary step unless the results
*  overflows.
*
*  Arguments
*  =========
*
*  X       (input) COMPLEX*16
*  Y       (input) COMPLEX*16
*          The complex scalars X and Y.
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   ZI, ZR
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLADIV
*     ..
*     .. Intrinsic Functions ..
*PAN  INTRINSIC          DBLE, DCMPLX, DIMAG
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
      CALL DLADIV( DBLE( X ), DIMAG( X ), DBLE( Y ), DIMAG( Y ), ZR,
     $             ZI )
      ZLADIV = DCMPLX( ZR, ZI )
*
      RETURN
*
*     End of ZLADIV
*
      END
      SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLASSQ returns the values scl and ssq such that
*
*     ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is
*  assumed to be at least unity and the value of ssq will then satisfy
*
*     1.0 .le. ssq .le. ( sumsq + 2*n ).
*
*  scale is assumed to be non-negative and scl returns the value
*
*     scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
*            i
*
*  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
*  SCALE and SUMSQ are overwritten by scl and ssq respectively.
*
*  The routine makes only one pass through the vector X.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) DOUBLE PRECISION
*          The vector x as described above.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with the value  scl .
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with the value  ssq .
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   TEMP1
*     ..
*     .. Intrinsic Functions ..
*PAN      INTRINSIC          ABS, DBLE, DIMAG
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( DBLE( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( DBLE( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
            IF( DIMAG( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( DIMAG( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASSQ
*
      END
      SUBROUTINE DLADIV( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
*     ..
*
*  Purpose
*  =======
*
*  DLADIV performs complex division in  real arithmetic
*
*                        a + i*b
*             p + i*q = ---------
*                        c + i*d
*
*  The algorithm is due to Robert L. Smith and can be found
*  in D. Knuth, The art of Computer Programming, Vol.2, p.195
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*  C       (input) DOUBLE PRECISION
*  D       (input) DOUBLE PRECISION
*          The scalars a, b, c, and d in the above expression.
*
*  P       (output) DOUBLE PRECISION
*  Q       (output) DOUBLE PRECISION
*          The scalars p and q in the above expression.
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   E, F
*     ..
*     .. Intrinsic Functions ..
*PAN      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( ABS( D ).LT.ABS( C ) ) THEN
         E = D / C
         F = C + D*E
         P = ( A+B*E ) / F
         Q = ( B-A*E ) / F
      ELSE
         E = C / D
         F = D + C*E
         P = ( B+A*E ) / F
         Q = ( -A+B*E ) / F
      END IF
*
      RETURN
*
*     End of DLADIV
*
      END
      SUBROUTINE MINA(FN,NV,NDIV,DEL,A,GUESS,X,FOFX,IERR)               MINA.2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NOP = 86)
c                                                                       ADDRESS.2
c     SANDIA MATHEMATICAL PROGRAM LIBRARY                               ADDRESS.3
c     APPLIED MATHEMATICS DIVISION 2613                                 ADDRESS.
c     SANDIA LABORATORIES                                               ADDRESS.5
c     ALBUQUERQUE, NEW MEXICO  87185                                    JUN0278.1
c     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978                     JUN0278.2
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.1
c                    ISSUED BY SANDIA LABORATORIES                     *MAR1378.2
c  *                   A PRIME CONTRACTOR TO THE                       *MAR1378.3
c  *                UNITED STATES DEPARTMENT OF ENERGY                 *MAR1378.4
c  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * * *MAR1378.5
c  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE   *MAR1378.6
c  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE      *MAR1378.7
c  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,    *MAR1378.8
c  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES  *MAR1378.9
c  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL      *MAR1378.10
c  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR     *MAR1378.11
c  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS      *MAR1378.12
c  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE          *MAR1378.13
c  * OWNED RIGHTS.                                                     *MAR1378.14
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.15
c  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS     *MAR1378.16
c  * PART IS SAND77-1441.                                              *MAR1378.17
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.18
c                                                                       ADDRESS.9
c     ORIGINAL ROUTINE WAS H2 SAND MIN, BY Z. BEISINGER AND S. BELL     MINA.4
c     PRESENT VERSION BY R E JONES                                      MINA.5
c                                                                       MINA.6
c     ABSTRACT                                                          MINA.7
c        MINA FINDS AN APPROXIMATE MINIMUM OF A REAL FUNCTION OF        MINA.8
c        NV VARIABLES, GIVEN AN INITIAL ESTIMATE OF THE POSITION OF     MINA.9
c        THE MINIMUM AND RANGES FOR EACH OF THE VARIABLES.              MINA.10
c        MINA USES A SELECTIVE DIRECTED SEARCH OF A SURROUNDING         MINA.11
c        NV-DIMENSIONAL GRID OF POINTS TO FIND A DIRECTION IN WHICH     MINA.12
c        THE FUNCTION DECREASES.  IT THEN PROCEEDS IN THIS DIRECTION    MINA.13
c        AS FAR AS THE FUNCTION DECREASES, THEN DETERMINES A NEW        MINA.14
c        DIRECTION TO TRAVEL.  WHEN NO SUCH DIRECTION IS FOUND THE      MINA.15
c        SEARCH INCREMENT FACTOR IS DECREASED AND THE PROCESS           MINA.16
c        IS REPEATED.                                                   MINA.17
c                                                                       MINA.18
c     DESCRIPTION OF ARGUMENTS                                          MINA.19
c        THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST  MINA.20
c              A(NV,2), GUESS(NV), X(NV)                                MINA.21
c                                                                       MINA.22
c        INPUT--                                                        MINA.23
c        FN   - NAME OF FUNCTION OF NV VARIABLES TO BE MINIMIZED.       MINA.24
c               (THIS NAME MUST APPEAR IN AN EXTERNAL STATEMENT.)       MINA.25
c               FORM OF THE CALLING SEQUENCE MUST BE FUNCTION FN(X),    MINA.26
c               WHERE X IS AN ARRAY OF NV VARIABLE VALUES. THE          MINA.27
c               ORDERING OF THE VARIABLES IS ARBITRARY, EXCEPT          MINA.28
c               THAT IT MUST AGREE WITH THE ORDERING USED IN            MINA.29
c               ARRAYS A AND GUESS.                                     MINA.30
c        NV   - NUMBER OF VARIABLES.  (NV .GE. 1)                       MINA.31
c        NDIV - NUMBER OF REFINEMENTS OF THE SEARCH INCREMENTS TO USE.  MINA.32
c               AT EACH REFINEMENT, THE INCREMENT IN EACH DIMENSION     MINA.33
c               IS DIVIDED BY 10.  (USUALLY NDIV IS ABOUT 3 OR 4.)      MINA.34
c        DEL  - FRACTION OF VARIABLE RANGE (IN EACH DIMENSION) TO USE   MINA.35
c               AS THE INITIAL INCREMENT (IN THAT DIMENSION)            MINA.36
c        A    - ARRAY OF SEARCH BOUNDS, DIMENSIONED NV BY 2.            MINA.37
c               A(I,1) SHOULD BE THE LOWER BOUND OF THE I-TH VARIABLE.  MINA.38
c               A(I,2) SHOULD BE THE UPPER BOUND OF THE I-TH VARIABLE.  MINA.39
c        GUESS- ARRAY OF NV INITIAL VALUES.  GUESS(I) SHOULD BE THE     MINA.40
c               INITIAL VALUE TO USE FOR THE I-TH VARIABLE.             MINA.41
c                                                                       MINA.42
c        OUTPUT--                                                       MINA.43
c        X    - ARRAY (DIMENSIONED NV) GIVING THE VALUES OF THE         MINA.44
c               VARIABLES AT THE MINIMUM.  X(I) WILL BE THE VALUE       MINA.45
c               OF THE I-TH VARIABLE.                                   MINA.46
c        FOFX - FUNCTION VALUE AT THE MINIMUM                           MINA.47
c        IERR - A STATUS CODE                                           MINA.48
c              -NORMAL CODE                                             MINA.49
c               =1 MEANS THE SEARCH FOR A MINIMUM PROCEEDED FOR THE     MINA.50
c                  SPECIFIED NUMBER OF REFINEMENTS.                     MINA.51
c              -ABNORMAL CODES                                          MINA.52
c               =2 MEANS NV IS GREATER THAN 50                          MINA.53
c               =3 MEANS A RANGE MINIMUM IS GREATER THAN THE            MINA.54
c                  CORRESPONDING MAXIMUM                                MINA.55
c                                                                       MINA.57
      DIMENSION A(NOP,2),GUESS(NOP),X(NOP)                              MINA.58
      DIMENSION XNOW(NOP),XNEW(NOP),R(NOP)                                 MINA.59
c                                                                       MESS.2
c     CALL MLMESS(10HSANDIA7.2 )                                        MESS3.1
c                                                                       MESS.4
c                                                                       MINA.60
c     INITIALIZE                                                        MINA.61
c                                                                       MINA.62
      IERR = 1                                                          MINA.64
      IF (NV.LE.NOP) GO TO 2                                             MINA.65
cc      CALL ERRCHK(33,33HIN MINA  , NV IS GREATER THAN 100.)            MINA.66
      write(6,*)'nv is greater than 100'
      IERR = 2                                                          MINA.67
      RETURN                                                            MINA.68
    2 NX = NV                                                           MINA.69
      IDIV = 0                                                          MINA.70
      DO 5 I=1,NX                                                       MINA.71
      XNOW(I) = GUESS(I)                                                MINA.72
      IF (XNOW(I).LT.A(I,1)) XNOW(I) = A(I,1)                           MINA.73
      IF (XNOW(I).GT.A(I,2)) XNOW(I) = A(I,2)                           MINA.74
      IF (A(I,1)-A(I,2)) 5,5,4                                          MINA.75
    4  write(6,*) 'range minimum greater than maximum' 
       print *,A(i,1),A(i,2)
c    4 CALL ERRCHK(46,46HIN MINA  , RANGE MINIMUM GREATER THAN MAXIMUM.) MINA.76
      IERR = 3                                                          MINA.77
      RETURN                                                            MINA.78
    5 R(I) = A(I,2)-A(I,1)                                              MINA.79
      DELTA = DEL                                                       MINA.80
      IF (DELTA.LE.0.0) DELTA = 0.1                                     MINA.81
      FNOW = FN(XNOW)                                                   MINA.82
c                                                                       MINA.83
c     FIND NEW DIRECTION                                                MINA.84
c                                                                       MINA.85
    7 DO 8 I=1,NX                                                       MINA.86
    8 XNEW(I) = XNOW(I)                                                 MINA.87
      FOLD = FNOW                                                       MINA.88
   10 DO 40 I=1,NX                                                      MINA.89
      IF (XNOW(I).GE.A(I,2)) GO TO 20                                   MINA.90
      XNEW(I) = MIN(XNOW(I)+DELTA*R(I),A(I,2))
      FNEW = FN(XNEW)                                                   MINA.92
      IF (FNEW.LT.FNOW) GO TO 30                                        MINA.93
   20 IF (XNOW(I) .LE. A(I,1)) GO TO 25                                 MINA.94
      XNEW(I) = MAX(XNOW(I)-DELTA*R(I),A(I,1))
      FNEW = FN(XNEW)                                                   MINA.96
      IF (FNEW.LT.FNOW) GO TO 30                                        MINA.97
   25 XNEW(I) = XNOW(I)                                                 MINA.98
      GO TO 40                                                          MINA.99
   30 FNOW = FNEW                                                       MINA.100
   40 CONTINUE                                                          MINA.101
      ISTEP = 1                                                         MINA.102
c                                                                       MINA.103
c     REFINE IF NEEDED                                                  MINA.104
c                                                                       MINA.105
      IF (FNOW.LT.FOLD) GO TO 50                                        MINA.106
      IF (IDIV.GE.NDIV) GO TO 100                                       MINA.107
      DELTA = DELTA*0.1                                                 MINA.108
      IDIV = IDIV+1                                                     MINA.109
      GO TO 10                                                          MINA.110
c                                                                       MINA.111
c     TRY TO CONTINUE IN CHOSEN DIRECTION                               MINA.112
c                                                                       MINA.113
   50 ICHNG = 0                                                         MINA.114
      FAc = 1.0                                                         MINA.115
      IF ((ISTEP/10)*10.EQ.ISTEP) FAc = 2.0                             MINA.116
      DO 60 I=1,NX                                                      MINA.117
      DX = (XNEW(I)-XNOW(I))*FAc                                        MINA.118
      XNOW(I) = XNEW(I)                                                 MINA.119
      IF (DX) 52,54,56                                                  MINA.120
   52 XNEW(I) = MAX(XNOW(I)+DX,A(I,1))
      IF (XNEW(I).LT.XNOW(I)) ICHNG = 1                                 MINA.122
      GO TO 60                                                          MINA.123
   54 XNEW(I) = XNOW(I)                                                 MINA.124
      GO TO 60                                                          MINA.125
   56 XNEW(I) = MIN(XNOW(I)+DX,A(I,2))
      IF (XNEW(I).GT.XNOW(I)) ICHNG = 1                                 MINA.127
   60 CONTINUE                                                          MINA.128
      IF (ICHNG.EQ.0) GO TO 7                                           MINA.129
      FNEW = FN(XNEW)                                                   MINA.130
      IF (FNEW.GE.FNOW) GO TO 7                                         MINA.131
      FNOW = FNEW                                                       MINA.132
      ISTEP = ISTEP+1                                                   MINA.133
      GO TO 50                                                          MINA.134
c                                                                       MINA.135
c     RETURN ANSWERS                                                    MINA.136
c                                                                       MINA.137
  100 FOFX = FOLD                                                       MINA.138
      DO 110 I=1,NX                                                     MINA.139
  110 X(I) = XNOW(I)                                                    MINA.140
      RETURN                                                            MINA.141
      END                                                               MINA.142
      SUBROUTINE SIMIN (F,K,EPS,ANS,S,NEV,ICONT,Y)                      SIMIN.2
      IMPLICIT double precision (A-H,O-Z)
c                                                                       ADDRESS.2
c     SANDIA MATHEMATICAL PROGRAM LIBRARY                               ADDRESS.3
c     APPLIED MATHEMATICS DIVISION 2613                                 ADDRESS.4
c     SANDIA LABORATORIES                                               ADDRESS.5
c     ALBUQUERQUE, NEW MEXICO  87185                                    JUN0278.1
c     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978                     JUN0278.2
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.1
c                    ISSUED BY SANDIA LABORATORIES                     *MAR1378.2
c  *                   A PRIME CONTRACTOR TO THE                       *MAR1378.3
c  *                UNITED STATES DEPARTMENT OF ENERGY                 *MAR1378.4
c  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * * *MAR1378.5
c  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE   *MAR1378.6
c  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE      *MAR1378.7
c  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,    *MAR1378.8
c  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES  *MAR1378.9
c  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL      *MAR1378.10
c  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR     *MAR1378.11
c  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS      *MAR1378.12
c  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE          *MAR1378.13
c  * OWNED RIGHTS.                                                     *MAR1378.14
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.15
c  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS     *MAR1378.16
c  * PART IS SAND77-1441.                                              *MAR1378.17
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.18
c                                                                       ADDRESS.9
c     ORIGINAL ROUTINE BY L F SHAMPINE, AS DESCRIBED IN REF.1 BELOW.    SIMIN.4
c     PREPARATION FOR MATH LIBRARY BY R E JONES.                        SIMIN.5
c                                                                       SIMIN.6
c     ABSTRACT                                                          SIMIN.7
c        SIMIN FINDS AN APPROXIMATE MINIMUM OF A REAL FUNCTION OF K     SIMIN.8
c        VARIABLES, GIVEN AN INITIAL ESTIMATE OF THE POSITION OF THE    SIMIN.9
c        MINIMUM. THE SIMPLEX METHOD IS USED. SEE REFERENCE 1 BELOW     SIMIN.10
c        FOR A FULL EXPLANATION OF THIS METHOD.  BRIEFLY, A SET OF      SIMIN.11
c        K+1 POINTS IN K-DIMENSIONAL SPACE IS CALLED A SIMPLEX.         SIMIN.12
c        THE MINIMIZATION PROCESS ITERATES BY REPLACING THE POINT       SIMIN.13
c        WITH THE LARGEST FUNCTION VALUE BY A NEW POINT WITH A          SIMIN.14
c        SMALLER FUNCTION VALUE.  ITERATION CONTINUES UNTIL ALL THE     SIMIN.15
c        POINTS CLUSTER SUFFICIENTLY CLOSE TO A MINIMUM.                SIMIN.16
c                                                                       SIMIN.17
c     REFERENCES                                                        SIMIN.18
c      1. L F SHAMPINE, A ROUTINE FOR UNCONSTRAINED OPTIMIZATION,       SIMIN.19
c         SC-TM-72130  OR  SC-RR-720657                                 SIMIN.20
c      2. J A NELDER AND R MEAD, A SIMPLEX METHOD FOR FUNCTION          SIMIN.21
c         MINIMIZATION, COMPUTER JOURNAL, 7(1965) 308-313               SIMIN.22
c                                                                       SIMIN.23
c     DESCRIPTION OF PARAMETERS                                         SIMIN.24
c      --INPUT--                                                        SIMIN.25
c        F  - NAME OF FUNCTION OF K VARIABLES TO BE MINIMIZED.          SIMIN.26
c             (THIS NAME MUST APPEAR IN AN EXTERNAL STATEMENT.)         SIMIN.27
c             FORM OF THE CALLING SEQUENCE MUST BE FUNCTION F(X),       SIMIN.28
c             WHERE X IS AN ARRAY OF K VARIABLES.                       SIMIN.29
c        K  - THE NUMBER OF VARIABLES.  K MUST BE AT LEAST 2.           SIMIN.30
c             NORMALLY K SHOULD BE LESS THAN ABOUT 10, AS SIMIN         SIMIN.31
c             BECOMES LESS EFFECTIVE FOR LARGER VALUES OF K.            SIMIN.32
c        EPS- THE CONVERGENCE CRITERION.  LET YAVG BE THE AVERAGE       SIMIN.33
c             VALUE OF THE FUNCTION F AT THE K+1 POINTS OF THE          SIMIN.34
c             SIMPLEX, AND LET R BE THEIR STANDARD ERROR.  (THAT IS,    SIMIN.35
c             THE ROOT-MEAN-SQUARE OF THE SET OF VALUES (Y(I)-YAVG),    SIMIN.36
c             WHERE Y(I) IS THE FUNCTION VALUE AT THE I-TH POINT OF     SIMIN.37
c             THE SIMPLEX.)  THEN--                                     SIMIN.38
c             IF EPS.GT.0, CONVERGENCE IS OBTAINED IF  R.LE.EPS.        SIMIN.39
c             IF EPS.LT.0, CONVERGENCE IS IF  R.LE.ABS(EPS*YAVG).       SIMIN.40
c             IF EPS=0, THE PROCESS WILL NOT CONVERGE BUT INSTEAD WILL  SIMIN.41
c             QUIT WHEN NEV FUNCTION EVALUATIONS HAVE BEEN USED.        SIMIN.42
c        ANS- AN ARRAY OF LENGTH K CONTAINING A GUESS FOR THE LOCATION  SIMIN.43
c             OF A MINIMUM OF F.                                        SIMIN.44
c        S  - A SCALE PARAMETER, WHICH MAY BE A SIMPLE VARIABLE OR AN   SIMIN.45
c             ARRAY OF LENGTH K.  USE OF AN ARRAY IS SIGNALLED BY       SIMIN.46
c             SETTING S(1) NEGATIVE.                                    SIMIN.47
c             -SIMPLE VARIABLE CASE.  HERE S IS THE LENGTH OF EACH      SIMIN.48
c             SIDE OF THE INITIAL SIMPLEX.  THUS, THE INITIAL SEARCH    SIMIN.49
c             RANGE IS THE SAME FOR ALL THE VARIABLES.                  SIMIN.50
c             -ARRAY CASE.  HERE THE LENGTH OF SIDE I OF THE INITIAL    SIMIN.51
c             SIMPLEX IS ABS(S(I)).  THUS, THE INITIAL SEARCH RANGE     SIMIN.52
c             MAY BE DIFFERENT FOR DIFFERENT VARIABLES.                 SIMIN.53
c             NOTE-- THE VALUE(S) USED FOR S ARE NOT VERY CRITICAL.     SIMIN.54
c             ANY REASONABLE GUESS SHOULD DO O.K.                       SIMIN.55
c        NEV- THE MAXIMUM NUMBER OF FUNCTION EVALUATIONS TO BE USED.    SIMIN.56
c             (THE ACTUAL NUMBER USED MAY EXCEED THIS SLIGHTLY SO THE   SIMIN.57
c             LAST SEARCH ITERATION MAY BE COMPLETED.)                  SIMIN.58
c        ICONT - ICONT SHOULD BE ZERO ON ANY CALL TO SIMIN WHICH        SIMIN.59
c             IS NOT A CONTINUATION OF A PREVIOUS CALL.                 SIMIN.60
c             IF ICONT=1 THE PROBLEM WILL BE CONTINUED.  IN THIS        SIMIN.61
c             CASE THE WORK ARRAY Y MUST BE THE SAME ARRAY THAT WAS     SIMIN.62
c             USED IN THE CALL THAT IS BEING CONTINUED (AND THE VALUES  SIMIN.63
c             IN IT MUST BE UNCHANGED).  THE REASON FOR THIS IS THAT    SIMIN.64
c             IF ICONT=1 THEN THE ARGUMENT S IS IGNORED AND THE SIMPLEX SIMIN.65
c             AND RELATED FUNCTION VALUES THAT WERE STORED IN ARRAY Y   SIMIN.66
c             DURING A PREVIOUS EXECUTION ARE USED TO CONTINUE THAT     SIMIN.67
c             PREVIOUS PROBLEM.                                         SIMIN.68
c        Y  - A WORK ARRAY CONTAINING AT LEAST K*K + 5*K + 1 WORDS.     SIMIN.69
c             IF ICONT=1 THIS MUST BE THE SAME ARRAY USED IN THE CALL   SIMIN.70
c             THAT IS BEING CONTINUED.                                  SIMIN.71
c      --OUTPUT--                                                       SIMIN.72
c        ANS- ANS WILL CONTAIN THE LOCATION OF THE POINT WITH THE       SIMIN.73
c             SMALLEST VALUE OF THE FUNCTION THAT WAS FOUND.            SIMIN.74
c        S  - IN THE SIMPLE VARIABLE CASE S WILL BE RETURNED AS THE     SIMIN.75
c             AVERAGE DISTANCE FROM THE VERTICES TO THE CENTROID OF     SIMIN.76
c             THE SIMPLEX.                                              SIMIN.77
c             IN THE ARRAY CASE S(I) WILL BE RETURNED AS THE AVERAGE    SIMIN.78
c             DISTANCE IN THE I-TH DIMENSION OF VERTICES FROM           SIMIN.79
c             THE CENTROID.  (S(1) WILL BE NEGATED.)                    SIMIN.80
c             NOTE-- THE VALUE(S) RETURNED IN S ARE USEFUL FOR          SIMIN.81
c             ASSESSING THE FLATNESS OF THE FUNCTION NEAR THE           SIMIN.82
c             MINIMUM.  THE LARGER THE VALUE OF S (FOR A GIVEN          SIMIN.83
c             VALUE OF EPS), THE FLATTER THE FUNCTION.                  SIMIN.84
c        NEV- NEV WILL BE THE COUNT OF THE ACTUAL NUMBER OF FUNCTION    SIMIN.85
c             EVALUATIONS USED.                                         SIMIN.86
c        Y  - WILL CONTAIN ALL DATA NEEDED TO CONTINUE THE MINIMIZATION SIMIN.87
c             SEARCH EFFICIENTLY IN A SUBSEQUENT CALL.                  SIMIN.88
c             NOTE -- THE FIRST K+1 ELEMENTS OF Y WILL CONTAIN THE      SIMIN.89
c             FUNCTION VALUES AT THE K+1 POINTS OF THE LATEST SIMPLEX.  SIMIN.90
c             THE NEXT K*(K+1) ELEMENTS OF Y WILL BE THE K+1 POINTS     SIMIN.91
c             OF THE SIMPLEX (IN EXACT CORRESPONDENSE TO THE ARRAY      SIMIN.92
c             P DISCUSSED IN REFERENCE 1 ABOVE).  THE REMAINING 3*K     SIMIN.93
c             WORDS ARE TEMPORARY WORKING STORAGE ONLY.                 SIMIN.94
c                                                                       SIMIN.96
      DIMENSION ANS(K),S(K),Y(1)                                        SIMIN.97
      EXTERNAL F                                                        SIMIN.98
c                                                                       MESS.2
c     CALL MLMESS(10HSANDIA7.2 )                                        MESS3.1
c                                                                       MESS.4
c                                                                       SIMIN.99
      IF(K.GE.2 .AND. S(1).NE.0.) GO TO 10                              SIMIN.101
      write(6,*) 's(1)=0 or k is less than 2'
c      CALL ERRCHK(39,39HIN SIMIN , S(1)=0. OR K IS LESS THAN 2.)        SIMIN.102
      RETURN                                                            SIMIN.103
c   10 IF (K.GT.100) CALL ERRCHK(31,31HIN SIMIN , K IS LARGER THAN 100.)   SIMIN.104
   10 write(6,*) 'k is large than 100'
c                                                                       SIMIN.105
      IP = K+2                                                          SIMIN.106
      Ic = IP+K*(K+1)                                                   SIMIN.107
      IR = IC+K                                                         SIMIN.108
      IRR = IR+K                                                        SIMIN.109
      CALL SIMINA(F,K,EPS,ANS,S,NEV,ICONT,Y,Y(IP),Y(IC),Y(IR),Y(IRR))   SIMIN.110
      RETURN                                                            SIMIN.111
      END                                                               SIMIN.112
      SUBROUTINE SIMINA(F,K,EPS,ANS,S,NEV,ICONT,Y,P,PC,PR,PRR)          SIMIN.113
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION ANS(K),S(K),Y(3),P(K,3),PC(K),PR(K),PRR(K)              SIMIN.114
      DATA ALPHA,BETA,GAMMA/1.0,0.5,2.0/                                SIMIN.115
c                                                                       SIMIN.116
c     SIMINA IS A CORE MINIMIZATION ROUTINE CALLED ONLY BY SIMIN.       SIMIN.117
c                                                                       SIMIN.118
c     INITIALIZATION                                                    SIMIN.119
c                                                                       SIMIN.120
      KOUNT = 0                                                         SIMIN.122
      KK = K                                                            SIMIN.123
      IF (KK.LE.1) GO TO 99                                             SIMIN.124
      ONEK = 1.0/FLOAT(KK)                                              SIMIN.125
      KP1 = KK+1                                                        SIMIN.126
      ONEKP1 = 1.0/FLOAT(KP1)                                           SIMIN.127
      TOL = FLOAT(KP1)*EPS**2                                           SIMIN.128
c                                                                       SIMIN.129
c     INITIAL SIMPLEX                                                   SIMIN.130
c                                                                       SIMIN.131
      IF (ICONT.GE.1) GO TO 10                                          SIMIN.132
      IF (S(1)) 4,99,1                                                  SIMIN.133
    1 SKP1 = S(1)*ONEKP1                                                SIMIN.134
      DO 2 I=1,KP1                                                      SIMIN.135
      DO 2 J=1,KK                                                       SIMIN.136
    2 P(J,I) = ANS(J) - SKP1                                            SIMIN.137
      DO 3 J=1,KK                                                       SIMIN.138
    3 P(J,J+1) = P(J,J+1) + S(1)                                        SIMIN.139
      GO TO 7                                                           SIMIN.140
c                                                                       SIMIN.141
    4 DO 5 I=1,KP1                                                      SIMIN.142
      DO 5 J=1,KK                                                       SIMIN.143
    5 P(J,I) = ANS(J) - ABS(S(J))*ONEKP1                                SIMIN.144
      DO 6 J=1,KK                                                       SIMIN.145
    6 P(J,J+1) = P(J,J+1) + ABS(S(J))                                   SIMIN.146
c                                                                       SIMIN.147
c     FUNCTION VALUES FOR INITIAL SIMPLEX                               SIMIN.148
c                                                                       SIMIN.149
    7 I1 = 1                                                            SIMIN.150
      DO 8 I=1,KP1                                                      SIMIN.151
      Y(I) = F(P(1,I))                                                  SIMIN.152
      IF (Y(I).GT.Y(I1)) I1 = I                                         SIMIN.153
    8 CONTINUE                                                          SIMIN.154
      YANS = F(ANS)                                                     SIMIN.155
      KOUNT = KP1+1                                                     SIMIN.156
      IF (YANS.GE.Y(I1)) GO TO 10                                       SIMIN.157
      Y(I1) = YANS                                                      SIMIN.158
      DO 9 J=1,KK                                                       SIMIN.159
    9 P(J,I1) = ANS(J)                                                  SIMIN.160
c                                                                       SIMIN.161
c     RE-START / NEXT ITERATION                                         SIMIN.162
c     IF K.LT.0 VALUES IN THE P AND Y ARRAYS (AND ONLY THESE VALUES)    SIMIN.163
c     WILL NOT HAVE BEEN DEFINED IN THIS CALL.  THIS IS NON-ANSI USAGE. SIMIN.164
c                                                                       SIMIN.165
c     FIRST FIND LARGEST, SECOND LARGEST, AND SMALLEST FUNCTION VALUES. SIMIN.166
c                                                                       SIMIN.167
   10 I1 = 1                                                            SIMIN.168
      IL = 1                                                            SIMIN.169
      DO 12 I=2,KP1                                                     SIMIN.170
      IF (Y(I).LT.Y(IL)) IL = I                                         SIMIN.171
      IF (Y(I).GT.Y(I1)) I1 = I                                         SIMIN.172
   12 CONTINUE                                                          SIMIN.173
      I2 = IL                                                           SIMIN.174
      DO 13 I=1,KP1                                                     SIMIN.175
      IF (I.EQ.I1) GO TO 13                                             SIMIN.176
      IF (Y(I).GT.Y(I2)) I2 = I                                         SIMIN.177
   13 CONTINUE                                                          SIMIN.178
c                                                                       SIMIN.179
c     COMPUTE CENTROID, LEAVING OUT P(*,I1)                             SIMIN.180
c                                                                       SIMIN.181
      DO 15 J=1,KK                                                      SIMIN.182
      SUM = 0.0                                                         SIMIN.183
      DO 14 I=1,KP1                                                     SIMIN.184
      IF (I.EQ.I1) GO TO 14                                             SIMIN.185
      SUM = SUM + P(J,I)                                                SIMIN.186
   14 CONTINUE                                                          SIMIN.187
   15 PC(J) = SUM*ONEK                                                  SIMIN.188
c                                                                       SIMIN.189
c     FORM REFLECTED POINT AND TEST                                     SIMIN.190
c                                                                       SIMIN.191
      DO 20 J=1,KK                                                      SIMIN.192
   20 PR(J) = PC(J) + ALPHA*(PC(J)-P(J,I1))                             SIMIN.193
      YR = F(PR)                                                        SIMIN.194
      KOUNT = KOUNT+1                                                   SIMIN.195
      IF (YR.LT.Y(IL)) GO TO 30                                         SIMIN.196
      IF (YR.GE.Y(I2)) GO TO 40                                         SIMIN.197
c                                                                       SIMIN.198
c     ACCEPT REFLECTED POINT                                            SIMIN.199
c                                                                       SIMIN.200
   21 Y(I1) = YR                                                        SIMIN.201
      DO 22 J=1,KK                                                      SIMIN.202
   22 P(J,I1) = PR(J)                                                   SIMIN.203
      GO TO 60                                                          SIMIN.204
c                                                                       SIMIN.205
c     EXPAND IN FAVORABLE DIRECTION AND TEST                            SIMIN.206
c                                                                       SIMIN.207
   30 DO 31 J=1,KK                                                      SIMIN.208
   31 PRR(J) = PR(J) + GAMMA*(PR(J)-PC(J))                              SIMIN.209
      YRR = F(PRR)                                                      SIMIN.210
      KOUNT = KOUNT+1                                                   SIMIN.211
      IF (YRR.GE.YR) GO TO 21                                           SIMIN.212
c                                                                       SIMIN.213
c     ACCEPT EXPANDED POINT                                             SIMIN.214
c                                                                       SIMIN.215
      Y(I1) = YRR                                                       SIMIN.216
      DO 32 J=1,KK                                                      SIMIN.217
   32 P(J,I1) = PRR(J)                                                  SIMIN.218
      GO TO 60                                                          SIMIN.219
c                                                                       SIMIN.220
c     DECIDE WHETHER TO ACCEPT REFLECTED POINT.                         SIMIN.221
c                                                                       SIMIN.222
   40 IF (YR.GE.Y(I1)) GO TO 42                                         SIMIN.223
      Y(I1) = YR                                                        SIMIN.224
      DO 41 J=1,KK                                                      SIMIN.225
   41 P(J,I1) = PR(J)                                                   SIMIN.226
c                                                                       SIMIN.227
c     TRY CONTRACTION.                                                  SIMIN.228
c                                                                       SIMIN.229
   42 DO 43 J=1,KK                                                      SIMIN.230
   43 PR(J) = PC(J) + BETA*(P(J,I1)-PC(J))                              SIMIN.231
      YCT = F(PR)                                                       SIMIN.232
      KOUNT = KOUNT+1                                                   SIMIN.233
      IF (YCT.GT.Y(I1)) GO TO 50                                        SIMIN.234
      Y(I1) = YCT                                                       SIMIN.235
      DO 44 J=1,KK                                                      SIMIN.236
   44 P(J,I1) = PR(J)                                                   SIMIN.237
      GO TO 60                                                          SIMIN.238
c                                                                       SIMIN.239
c     ALL EFFORTS FAILED.  SHRINK THE SIMPLEX ABOUT BEST POINT.         SIMIN.240
c                                                                       SIMIN.241
   50 DO 52 I=1,KP1                                                     SIMIN.242
      IF (I.EQ.IL) GO TO 52                                             SIMIN.243
      DO 51 J=1,KK                                                      SIMIN.244
   51 P(J,I) = 0.5*(P(J,I)+P(J,IL))                                     SIMIN.245
      Y(I) = F(P(1,I))                                                  SIMIN.246
   52 CONTINUE                                                          SIMIN.247
      KOUNT = KOUNT+KP1                                                 SIMIN.248
c                                                                       SIMIN.249
c     CHECK FOR CONVERGENCE                                             SIMIN.250
c                                                                       SIMIN.251
   60 IF (KOUNT.GE.NEV) GO TO 65                                        SIMIN.252
      IF (EPS.EQ.0.0) GO TO 10                                          SIMIN.253
      SUM = 0.0                                                         SIMIN.254
      DO 61 I=1,KP1                                                     SIMIN.255
   61 SUM = SUM + Y(I)                                                  SIMIN.256
      YAVG = SUM*ONEKP1                                                 SIMIN.257
      SUM = 0.0                                                         SIMIN.258
      DO 62 I=1,KP1                                                     SIMIN.259
   62 SUM = SUM + (Y(I)-YAVG)**2                                        SIMIN.260
      IF (EPS) 64,63,63                                                 SIMIN.261
   63 IF (SUM-TOL) 65,65,10                                             SIMIN.262
   64 IF (SUM-TOL*ABS(YAVG)) 65,65,10                                   SIMIN.263
c                                                                       SIMIN.264
c     CONVERGENCE OBTAINED.                                             SIMIN.265
c                                                                       SIMIN.266
c     COMPUTE CENTROID                                                  SIMIN.267
c                                                                       SIMIN.268
   65 DO 68 J=1,KK                                                      SIMIN.269
      SUM = 0.0                                                         SIMIN.270
      DO 67 I=1,KP1                                                     SIMIN.271
   67 SUM = SUM+P(J,I)                                                  SIMIN.272
   68 PC(J) = SUM*ONEKP1                                                SIMIN.273
      IF (S(1)) 73,69,69                                                SIMIN.274
c                                                                       SIMIN.275
c     COMPUTE S(1) AS AVERAGE DISTANCE OF VERTICES FROM CENTROID.       SIMIN.276
c                                                                       SIMIN.277
   69 DIST = 0.0                                                        SIMIN.278
      DO 71 I=1,KP1                                                     SIMIN.279
      SUM = 0.0                                                         SIMIN.280
      DO 70 J=1,KK                                                      SIMIN.281
   70 SUM = SUM + (P(J,I)-PC(J))**2                                     SIMIN.282
   71 DIST = DIST + SQRT(SUM)                                           SIMIN.283
      S(1) = DIST*ONEKP1                                                SIMIN.284
      GO TO 80                                                          SIMIN.285
c                                                                       SIMIN.286
c     COMPUTE S(J) AS AVERAGE DISTANCE IN J-TH DIMENSION OF             SIMIN.287
c     VERTICES FROM THE CENTROID.                                       SIMIN.288
c                                                                       SIMIN.289
   73 DO 75 J=1,KK                                                      SIMIN.290
      SUM = 0.0                                                         SIMIN.291
      DO 74 I=1,KP1                                                     SIMIN.292
   74 SUM = SUM + ABS(P(J,I)-PC(J))                                     SIMIN.293
   75 S(J) = SUM*ONEKP1                                                 SIMIN.294
      S(1) = -S(1)                                                      SIMIN.295
c                                                                       SIMIN.296
c     RETURN P(*,IL) AS ANSWER                                          SIMIN.297
c                                                                       SIMIN.298
   80 IL = 1                                                            SIMIN.299
      DO 82 I=2,KP1                                                     SIMIN.300
      IF (Y(I).LT.Y(IL)) IL = I                                         SIMIN.301
   82 CONTINUE                                                          SIMIN.302
      DO 84 J=1,KK                                                      SIMIN.303
   84 ANS(J) = P(J,IL)                                                  SIMIN.304
      NEV = KOUNT                                                       SIMIN.305
      RETURN                                                            SIMIN.306
c                                                                       SIMIN.307
c     ERROR MESSAGE                                                     SIMIN.308
c                                                                       SIMIN.309
c     99 CALL ERRCHK(39,39HIN SIMINA, S(1)=0. OR K IS LESS THAN 2.)        SIMIN.310
   99 write(6,*) 's(1)=0. or k is less than 2.'
      RETURN                                                            SIMIN.311
      END                                                               SIMIN.312



