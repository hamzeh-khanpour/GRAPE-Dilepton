*##################################################################
      double precision function PYR(IDUMMY)
      implicit NONE
*------------ Argument ------------
      integer IDUMMY
*----------------------------------
*----------- Function -----------
      double precision DRN
      external         DRN
*--------------------------------
*------- Local variables -------
      integer  Iseed
*-------------------------------
      PYR = DRN(Iseed)
      return
      end
*##################################################################
      double precision function RLU(IDUMMY)
      implicit NONE
*------------ Argument ------------
      integer IDUMMY
*----------------------------------
*----------- Function -----------
      double precision DRN
      external         DRN
*--------------------------------
*------- Local variables -------
      integer  Iseed
*-------------------------------
      RLU = DRN(Iseed)
      return
      end
*##################################################################
      double precision function RNDM(IDUMMY)
      implicit NONE
*------------ Argument ------------
      integer IDUMMY
*----------------------------------
*----------- Function -----------
      double precision DRN
      external         DRN
*--------------------------------
*------- Local variables -------
      integer  Iseed
*-------------------------------
      RNDM = DRN(Iseed)
      return
      end
*##################################################################
