*(Written by T.Abe on 16/03/2000)
      double precision function  Alpha_Run(Q2,alpha0)
      implicit NONE
*----------- Arguments -----------
      double precision  Q2,alpha0
*---------------------------------
*----------- Functions -----------
      double precision  Ph_Alpha_Run
      external          Ph_Alpha_Run
*---------------------------------
*-------- Local variables --------
      integer  i
      double precision  am(3), CR, am2, pi, Q2_loc, Q2_sav,Alpha_Run_sav
      parameter ( pi=3.14159265359D0 )
      save  am, Q2_sav,Alpha_Run_sav
*---------------------------------
*------------- DATA -------------
      data  am/ 0.51099907D-3, 105.658389D-3, 1.77705D0 /
      data  Q2_sav/0D0/, Alpha_Run_sav/1D0/
*--------------------------------
      if (Q2 .EQ. Q2_sav) then
        Alpha_Run = Alpha_Run_sav
        RETURN
      else
        Q2_loc = abs(Q2)
        Q2_sav = Q2
      endif
      CR = 0D0
      do 100 i = 1, 3
        am2 = am(i)*am(i)
        if (Q2_loc .LE. am2)  GOTO 100
        CR = CR + log(Q2_loc/am2)
 100  continue
      CR = alpha0/3D0/pi*CR + Ph_Alpha_Run(Q2_loc)
      Alpha_Run = 1D0 /(1D0-CR)
      Alpha_Run_sav = Alpha_Run
      return
      end
*##########################################################################
*(Ref: Phys.Lett.B356 (1995) 398-403)
      double precision function  Ph_Alpha_Run(Q2)
      implicit NONE
*----------- Arguments -----------
      double precision  Q2
*---------------------------------
*---------------------------------
*-------- Local variables --------
      integer      Num_reg, i
       parameter ( Num_reg=5 )
      double precision  A(Num_reg),B(Num_reg),C(Num_reg)
     &                 ,Q2max(Num_reg)
      save  A,B,C,Q2max
*---------------------------------
*------------- DATA -------------
      data  Q2max/ 4D0, 16D0, 100D0, 8317.44D0, 1D10 /
      data  A/ 0D0, 0D0, 0D0, 0.00122270D0, 0.00164178D0 /
      data  B/ 0.00228770D0, 0.00251507D0, 0.00279328D0, 0.00296694D0
     &        ,0.00292051D0 /
      data  C/ 4.08041425D0, 3.09624477D0, 2.07463133D0, 1D0, 1D0 /
*--------------------------------
      do 1000 i = 1, Num_reg
        if (Q2 .LT. Q2max(i)) then
          Ph_Alpha_Run = A(i) + B(i)*log(1D0+C(i)*Q2)
          RETURN
        endif
 1000 continue
      Ph_Alpha_Run = 0D0
      return
      end
*##########################################################################
