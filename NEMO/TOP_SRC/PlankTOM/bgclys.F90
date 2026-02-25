       SUBROUTINE bgclys
#if defined key_planktom && defined key_top
!!======================================================================
!!                         ***  ROUTINE bgclys  ***
!! TOP : PlankTOM ecosystem model 
!!======================================================================
!! History : 
!!      original h3clys : 1988-07 E. MAIER-REIMER      MPI HAMBURG
!!      additions       : 1998    O. Aumont
!!      modifications   : 1999    C. Le Quere
!!   
!!----------------------------------------------------------------------
!!    'key_trc_piic'      :    additional tracer with pre-industrial DIC 
!!----------------------------------------------------------------------
!! ** Purpose : CALCULATES DEGREE OF CACO3 SATURATION IN THE WATER
!!              COLUMN, DISSOLUTION/PRECIPITATION OF CACO3 AND LOSS
!!              OF CACO3 TO THE CACO3 SEDIMENT POOL.
!! ** Action  : 
!!
!! References : 
!!          SCARBOROUGH, J. (1958) NUMERICAL MATHEMATICAL ANALYSIS.
!!          OXFORD UNIVERSITY PRESS, LONDON, 4TH ED., 576 PP..
!!          PlankTOM12 manual, Buitenhuis et al. 2023
!!                                 https://zenodo.org/records/8388158 
!!
!!    Method : 
!!          [H+] AND [CO3--] FOR THE ACTUAL TIME STEP ARE CALCULATED
!!     BY NEWTON-RAPHSON ITERATION (E.G. SCARBOROUGH, 1958).
!!
!!    External:  NONE
!!
!!      VARIABLE           TYPE    PURPOSE.
!!      --------           ----    --------
!! 
!!      *NYEAR*            INTEGER COUNTS TIMESTEPS (YEARS) OF INTEGRATION
!!                           (INTEGER, INPUT)
!!      *KITTER*           INTEGER SETS UPPER LIMIT FOR NUMBER OF ITERATIONS
!!                                 TO DETERMINE [CO3--] AND [H+]
!!      *AKW*              REAL    APPROXIMATE VALUE OF IONIC PRODUCT OF
!!                                 WATER
!!      *H*                REAL    [H+], DUMMY VARIABLE
!!      *ALKA*             REAL    GIVEN ALKALINITY [EQV/L], DUMMY VARIABLE
!!      *C*                REAL    GIVEN [SUM(12C)O2] [MOLE/L], DUMMY VARIABLE
!!      *A*                REAL    ALKALINITY [EQV/L] AS FUNCTION OF [CO3--]
!!                                 AND [H+], DUMMY VARIABLE
!!      *DELCO3*           REAL    DEVIATION OF ACTUAL CACO3 CONCENTRATION FROM
!!                                 SATURATION VALUE
!!      *UNDSAT*           REAL    UNDERSATURATION OF CACO3 (E.G. 3.=THREEFOLD)
!!      *EXCE14*           REAL    SUPERSATURATION IN CA(14C)O3 (E.G. 3.=
!!                                 THREEFOLD)
!!      *DISP14*           REAL    FRACTION CACO3 (14C) THAT IS DISSOLVED
!!      *EXCE13*           REAL    SUPERSATURATION IN CA(13C)O3 (E.G. 3.=
!!                                 THREEFOLD)
!!      *DISP13*           REAL    FRACTION CACO3 (13C) THAT IS DISSOLVED
!!      *SEDLOS*           REAL    FRACTION OF CACO3 IN THE BOTTOM WATER LAYER
!!                                 LOST TO THE CACO3 SEDIMENT POOL
!!      *SEDLOI*           REAL    FRACTION OF CACO3 IN THE BOTTOM WATER LAYER
!!                                 WHICH REMAINS IN THE WATER COLUMN
!! ---------------------------------------------------------------------------
! parameters and commons
      USE trp_trc
      USE sms_planktom
      USE oce_trc
      IMPLICIT NONE
! 
! local variables
      INTEGER ji, jj, jk, jn
      INTEGER kitter
      REAL bot, alka, a, c, delco3, h, omecal, omeara, remara, remco3, ah2
!
!! -----------------------------------------------------------
!!* 2. ITERATION TO DETERMINE [CO3--] AND [H+]
!!     (NEWTON-RAPHSON METHOD:
!!     THE VALUES OF [SUM(CO2)] AND [ALK] ARE GIVEN,
!!     DESIRED ROOTS OF [CO3--] AND [H+] FOR THAT PAIR
!!     ARE DETERMINED BY SOLVING NUMERICALLY THE SYSTEM
!!     OF THE TWO NONLINEAR EQUATIONS
!!     1) [ALK]GIVEN      - [ALK]([CO3--],[H+])      = 0 (=F)
!!     2) [SUM(CO2)]GIVEN - [SUM(CO2)]([CO3--],[H+]) = 0 (=GG)
!! -----------------------------------------------------------
!
!* 2.1  SET MAX. NUMBER OF ITERATIONS
! -----------------------------------
!
      kitter = 10
!
!* 2.4  BEGIN OF ITERATION
! ------------------------
!
      DO jn = 1, kitter
!
!* 2.5  COMPUTE [H+] CONCENTRATIONS
! ---------------------------------
!
        DO jk = 1, jpkm1
          DO jj = 2, nlcj-1
            DO ji = 2, nlci-1
!
!* 2.6  SET DUMMY VARIABLE FOR TOTAL BORATE
! -----------------------------------------
!
              bot = borat(ji,jj,jk)
!
!* 2.7  SET DUMMY VARIABLE FOR [H+]
! ---------------------------------
!
              h = amax1(hi(ji,jj,jk),1.E-10)
!
!* 2.8  SET DUMMY VARIABLE FOR [ALK]GIVEN AND [SUM(CO2)]GIVEN
! -----------------------------------------------------------
!
              alka = trn(ji,jj,jk,jptal) 
              c = trn(ji,jj,jk,jpdic) 
!
!* 2.9 CALCULATE Carbonate Alkalinity
! ------------------------------------
!
              A = alka- &
     &            (akw3(ji,jj,jk)/h-h+bot/(1.+h/akb3(ji,jj,jk)))
!
!* 2.10 CALCULATE [H+]
! --------------------
!
              ah2 = sqrt((c-a)**2+4.*(a*ak23(ji,jj,jk)/ak13(ji,jj,jk)) &
     &            *(2*c-a))
              ah2 = 0.5*ak13(ji,jj,jk)/a*((c-a)+ah2)
              hi(ji,jj,jk) = ah2
#  if defined key_trc_piic
!
! Calculate the pre-industrial component of H+ concentration
! ----------------------------------------------------------
!
              h = amax1(pihi(ji,jj,jk),1.E-10)
              c = trn(ji,jj,jk,jppiic)
              A = alka- &
     &            (akw3(ji,jj,jk)/h-h+bot/(1.+h/akb3(ji,jj,jk)))
              ah2 = sqrt((c-a)**2+4.*(a*ak23(ji,jj,jk)/ak13(ji,jj,jk)) &
     &            *(2*c-a))
              ah2 = 0.5*ak13(ji,jj,jk)/a*((c-a)+ah2)
              pihi(ji,jj,jk) = ah2
#  endif
            END DO
          END DO
        END DO
      END DO 
!
!* 2.11 COMPUTE [CO3--] AND [HCO3-] CONCENTRATIONS
! ------------------------------------------------
!
      DO jk = 1, jpkm1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
            bot = borat(ji,jj,jk)
            A = trn(ji,jj,jk,jptal)-(akw3(ji,jj,jk)/hi(ji,jj,jk) &
     &        -hi(ji,jj,jk)+bot/(1.+hi(ji,jj,jk)/akb3(ji,jj,jk)))
            co3(ji,jj,jk) = a/(2.+hi(ji,jj,jk)/ak23(ji,jj,jk))
            hco3(ji,jj,jk) = a/(1.+2.*ak23(ji,jj,jk)/hi(ji,jj,jk))
!
!     ---------------------------------------------------------
!*    3. CALCULATE DEGREE OF CACO3 SATURATION AND CORRESPONDING
!        DISSOLOUTION OF CACO3 (BE AWARE OF MGCO3)
!     ---------------------------------------------------------
            omecal = co3(ji,jj,jk)*calcon/aksp(ji,jj,jk)
!
!* 3.3  AMOUNT CACO3 (12C) THAT RE-ENTERS SOLUTION
!       (ACCORDING TO THIS FORMULATION ALSO SOME PARTICULATE
!       CACO3 GETS DISSOLVED EVEN IN THE CASE OF OVERSATURATION)
! --------------------------------------------------------------
!       WHEN YOU CHANGE THIS, CHANGE THE SEDIMENT MODEL IN bgcbio AS WELL
!
            remco3 = max(trn(ji,jj,jk,jpcal)*rn_lyscal/rjjss*(1.-omecal),0.)**rn_lyoco3
            remco3 = min(trn(ji,jj,jk,jpcal)*rfactr,remco3)
!
!* 3.8  CHANGE OF [ALK], PARTICULATE [CACO3], AND DIC DUE TO CACO3 
!       DISSOLUTION
!       note the absense of rfact/raajj
! ----------------------------------------------------------------
!
            tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal)-remco3
            omeara = co3(ji,jj,jk)*calcon/aksara(ji,jj,jk)
            remara = max(trn(ji,jj,jk,jpara)*rn_lysara/rjjss*(1.-omeara),0.)**rn_lyoco3
            remara = min(trn(ji,jj,jk,jpara)*rfactr,remara)
            tra(ji,jj,jk,jpara) = tra(ji,jj,jk,jpara)-remara
            tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal)+ &
     &          2.*remco3+2.*remara
            tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic)+                  &
     &          remco3+remara
#  if defined key_trc_piic
!
! Ignore the fact that preindustrial dissolution is slightly lower
! because including it would necessitate including preindustrial 
! state variables for jptal, jpcal and jpara
!
            tra(ji,jj,jk,jppiic) = tra(ji,jj,jk,jppiic)+                  &
     &          remco3+remara
#  endif
            discarb(ji,jj,jk) = (remco3+remara)*1e3
          END DO
        END DO
      END DO
#endif
      RETURN
      END
