      subroutine eBpol(agc)
      implicit NONE
*----------- Arguments -----------
      complex*16  agc(2,2,2,2,2,2)
*---------------------------------
*------------ COMMONs ------------
      include './inc/graepia.h'
*---------------------------------
*----------- Functions -----------
*---------------------------------
*-------- Local variables --------
      integer  ih1,ih2,ih3,ih4,ih5,ih6
*---------------------------------
*------------- DATA -------------
*--------------------------------
      if (.NOT.Lpol_Ebeam)  RETURN
      do 1000 ih6 = 1, 2
      do 1000 ih5 = 1, 2
      do 1000 ih4 = 1, 2
      do 1000 ih3 = 1, 2
      do 1000 ih1 = 1, 2
 1000    continue
      return
      end
