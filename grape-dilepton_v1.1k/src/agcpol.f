      subroutine AgcPol(agc)
      implicit NONE
*----------- Arguments -----------
      complex*16 agc(2,2,2,2,2,2)
*---------------------------------
*------------ COMMONs ------------
      include './inc/graepia.h'
*---------------------------------
*-------- Local variables --------
      integer  ih1,ih3,ih4,ih5,ih6
      real*8   thetah,sin_thetah,cos_thetah, phi, rad, Pol
      parameter ( rad=1.74532925199D-2 )
      complex*16  zu,zd,zun,zdn,zzz,zi
      integer  icnt
      save  sin_thetah,cos_thetah, phi, Pol, icnt
*---------------------------------
*------------ DATA ------------
      data  icnt/0/
*------------------------------
      if (.NOT.Lpol_Ebeam)  RETURN
      if (icnt .EQ. 0) then
         Pol = Ebeam_pol(1)
          if (Pol .GT. +1D0)  Pol = +1D0
          if (Pol .LT. -1D0)  Pol = -1D0
         thetah = Ebeam_pol(2)/2D0 *rad
          sin_thetah = sin(thetah)
          cos_thetah = cos(thetah)
         phi    = Ebeam_pol(3)     *rad
         icnt = icnt + 1
      endif
      zi=(0D0,1D0)
      do 1000 ih6 = 1, 2
      do 1000 ih5 = 1, 2
      do 1000 ih4 = 1, 2
      do 1000 ih3 = 1, 2
      do 1000 ih1 = 1, 2
         if (thetah .EQ. 0D0) then
           zun = agc(ih1,2,ih3,ih4,ih5,ih6)
           zdn = agc(ih1,1,ih3,ih4,ih5,ih6)
         else
           zu = agc(ih1,2,ih3,ih4,ih5,ih6)
           zd = agc(ih1,1,ih3,ih4,ih5,ih6)
           zun =  cos_thetah             *zu + sin_thetah*exp(zi*phi)*zd
           zdn = -sin_thetah*exp(-zi*phi)*zu + cos_thetah            *zd
         endif
         zzz = sqrt(1D0+Pol)*zun + sqrt(1D0-Pol)*zdn
         agc(ih1,1,ih3,ih4,ih5,ih6) = zzz
         agc(ih1,2,ih3,ih4,ih5,ih6) = (0D0,0D0)
 1000 continue
      return
      end
