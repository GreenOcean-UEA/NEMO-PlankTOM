      SUBROUTINE bgcflx
#if defined key_planktom && defined key_top
!!======================================================================
!!                         ***  ROUTINE bgcflx  ***
!! TOP : PlankTOM ecosystem model 
!!======================================================================
!! History : original      : 1988-07 E. Maier-Reimer
!!           additions     : 1998    O. Aumont
!!           modifications : 1999&2004    C. Le Quere
!!           PlankTOM see bgcbio
!!----------------------------------------------------------------------
!!    'key_trc_n2o'       :                        calculates N2O fluxes
!!    'key_trc_ch4'       :                        calculates CH4 fluxes
!!    'key_trc_dms'       :            includes the DMS cycle and fluxes
!!    'key_trc_cfc11'     :            includes the CFC11 fluxes
!!    'key_trc_piic'      :    additional tracer with pre-industrial DIC
!!    'key_trc_diaadd'    :       save tracer diagnostic files 2D and 3D
!!                                         (to use with key_iomput only)
!!    'key_iomput'        :                              use IOM library
!!----------------------------------------------------------------------
!! ** Purpose  : calculates gas exchange and chemistry at sea surface 
!! ** Action   :  
!! ** External : none
!! ** Method   : solving system of two non-linear simultaneous equations 
!!               for [H2CO3]i and [H+] 
!!
!! References : 
!!              Broecker, W.S., and T.-H. Pend (1982) Tracers in the sea.
!!              Eldigio press, Lamont-Doherty Geological Observatory,
!!              Palisades, N.Y., 690 PP..
!!
!!              Mook, W.G. (1986), 13C in atmospheric CO2.
!!              Netherlands Journal of Sea Research, 20(2/3): 211-223.
!!
!!              Scarborough, J. (1958) Numerical Mathematical Analysis.
!!              Oxford Univertiry Press, London, 4TH ED., 576 PP..
!!
!!              Wanninkhof, R. (2014) doi:10.4319/lom.2014.12.351
!!              Limnology and Oceanography: Methods, 12: 351-362
!!
!!              PlankTOM12 manual, Buitenhuis et al. 2023
!!                                 https://zenodo.org/records/8388158 
!! 
!!      VARIABLE           TYPE    PURPOSE.
!!      --------           ----    --------
!!
!!      alka              REAL    actual [ALK] [EQV/L]
!!      akb               REAL    1. dissoc. constant of boric acid
!!      krorr          INTEGER    counts iterations in Newton-Raphson algorithm
!!      ct1               REAL    actual valua of inorg. C
!!      x1                REAL    [H+] expressed as sqrt(AK1*AK2)/[H+]
!!      zpa               REAL    atmospheric pco2 [ppm]
!!      fld               REAL    [CO2] flux atmosphere -> ocean in [mol/m2/s]
!!      flu               REAL    [CO2] flux ocean -> atmosphere in [mol/m2/s]
!!      relw13            REAL    ratio [sum((13C)O2)]/[sum((12C)O2)] in surface layer
!!      rela13            REAL    ratio [sum((13C)O2)]/[sum((12C)O2)] in the atmosphere
!!
!!----------------------------------------------------------------------
!
! parameters and commons
      USE trc
      USE trp_trc
      USE sms_planktom
      USE oce_trc
#  if defined key_iomput
      USE iom
#  endif
      USE wrk_nemo
      IMPLICIT NONE
! 
! local variables
      INTEGER ji, jj, krorr,yearrc
      CHARACTER (len=20) :: diaadd
      REAL ttc,ttc2,ttc3,ttc4, ws
      REAL ttcdms
      REAL za, zpa, zvapor, ztkel
      REAL oxy16
      REAL zph,ah2,bot
      REAL ct1,alka
      REAL schmico2, schmio2
      REAL schmidms
#  if defined key_trc_piic
      REAL piza,piph,piah2,piic
#  endif
#  if defined key_trc_n2o
      REAL schmin2o
#  endif
#  if defined key_trc_ch4
      REAL schmich4
#  endif
#  if defined key_trc_cfc11
      REAL schmicfc11,pcfc11,hemisp
#  endif
#  include "domzgr_substitute.h90"
!
! Initialisation of variables
! ----------------------------
! 
      dpco2  = 0.0
      pco2   = 0.0
      cflx   = 0.0
      flu16  = 0.0
#  if defined key_trc_dms
      fludms = 0.0
#  endif
      fld    = 0.0
      flu    = 0.0
!
!! ------------------------------------------------------------
!! Calculate surface chemistry (pCO2 AND [H+] in surface layer)
!! ------------------------------------------------------------
!
      yearrc = ndastp/10000-int(yrcfc(1))+1 !1 record per year in atmco2cfc.dat
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
!
! total borate
            bot  = borat(ji,jj,1)
            ct1  = trn(ji,jj,1,jpdic)
            alka = trn(ji,jj,1,jptal)
            zph  = amax1(hi(ji,jj,1),1.E-10)
#  if defined key_trc_piic
            piic = trn(ji,jj,1,jppiic)
            piph = amax1(pihi(ji,jj,1),1.E-10)
#  endif
!
! calculate [alk]([CO3--], [HCO3-])
          DO krorr = 1, 15
            za = alka- &
     &          (akw3(ji,jj,1)/zph-zph+bot/(1.+zph/akb3(ji,jj,1)))
!
! calculate [H+] AND [H2CO3]
            ah2 = sqrt((ct1-za)**2+4*(za*ak23(ji,jj,1)/ak13(ji,jj,1)) &
     &           * (2*ct1-za)) 
            zph = 0.5*ak13(ji,jj,1)/za*((ct1-za)+ah2)
#  if defined key_trc_piic
            piza = alka- &
     &          (akw3(ji,jj,1)/piph-piph+bot/(1.+piph/akb3(ji,jj,1)))
            piah2 = sqrt((piic-piza)**2+4*(piza*ak23(ji,jj,1)/ak13(ji,jj,1)) &
     &           * (2*piic-piza))
            piph = 0.5*ak13(ji,jj,1)/piza*((piic-piza)+piah2)
#  endif
          END DO
            hi(ji,jj,1)  = zph
            h2co3(ji,jj) = (2*ct1-za)/(2.+ak13(ji,jj,1)/zph)
#  if defined key_trc_piic
            pihi(ji,jj,1)  = piph
            pih2co3(ji,jj) = (2*piic-piza)/(2.+ak13(ji,jj,1)/piph)
#  endif
        END DO
      END DO
!
!! --------------
!! compute fluxes
!! --------------
!
! Compute gas exchange coefficients
! ---------------------------------
!
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
!
          ttc    = min(39.,tsnche(ji,jj,1,1))
          ttc2   = ttc**2
          ttc3   = ttc**3
          ttc4   = ttc**4
#  if defined key_trc_dms
          ttcdms = min(36.,tsnche(ji,jj,1,1))
#  endif
!
! use wind speed in m/s
          ws = wndm(ji,jj)
!
! this is Wanninkhof (2014) equation 12, in cm/h
          kgwanin(ji,jj) = 0.251*ws*ws 
!
! convert from cm/h to m/s and apply ice cover
          kgwanin(ji,jj) = kgwanin(ji,jj) /100./3600.  &
     &                   * (1-fr_i(ji,jj))
!
! Compute Schmidt number for CO2 (Wanninkhof 2014)
! ------------------------------------------------
!
          schmico2 = 2116.8-136.25*ttc+4.7353*ttc2-0.092307*ttc3+0.0007555*ttc4
          zkgco2(ji,jj) = kgwanin(ji,jj)*sqrt(660./schmico2)
!
! Compute Schmidt number for O2 (Wanninkhof 2014)
! -----------------------------------------------
!
          schmio2 = 1920.4-135.6*ttc+5.2122*ttc2-0.10939*ttc3+0.00093777*ttc4
          zkgo2(ji,jj) = kgwanin(ji,jj)*sqrt(660./schmio2)
!
#  if defined key_trc_dms
!
! Compute Schmidt number for DMS (Wanninkhof 2014)
!-------------------------------------------------
!
          schmidms = 2855.7-177.63*ttcdms+6.0438*ttcdms**2 &
                                 -0.11645*ttcdms**3+0.00094743*ttcdms**4
          zkgdms(ji,jj) = kgwanin(ji,jj)*sqrt(660./schmidms)
#  endif
#  if defined key_trc_n2o
!
! Compute Schmidt number for N2O (Wanninkhof 2014)
! ------------------------------------------------
!
          schmin2o = 2356.2-166.38*ttc+6.3952*ttc2-0.13422*ttc3+0.0011506*ttc4
          zkgn2o(ji,jj) = kgwanin(ji,jj)*sqrt(660./schmin2o)
#  endif
#  if defined key_trc_ch4
!
! Compute Schmidt number for CH4 (Wanninkhof 2014)
! ------------------------------------------------
!
          schmich4 = 2101.2-131.54*ttc+4.4931*ttc**2-0.08676*ttc**3+0.00070663*ttc**4
          zkgch4(ji,jj) = kgwanin(ji,jj)*sqrt(660./schmich4)
#  endif
#  if defined key_trc_cfc11
! 
! Compute Schmidt number for CFC-11 (Wanninkhof 2014)
! ---------------------------------------------------
!
          schmicfc11 = 3579.2-222.63*ttc+7.5749*ttc2-0.14595*ttc3+0.0011874*ttc4
          zkgcfc11(ji,jj) = kgwanin(ji,jj)*sqrt(660./schmicfc11)
#  endif
!
! Correct atmospheric CO2 (pco2at in ppm) for 100% water vapor 
! following Sarmiento et al., JGR 1992
! ------------------------------------------------------------
!
           ztkel  = tsnche(ji,jj,1,1)+temzer
           zvapor = exp( 20.1050 - 0.0097982 * ztkel - 6163.10/ztkel)
           zpa    = pco2at * (1. - zvapor)
!
! compute CO2 flux for the sea and air, 
! fld and flu are in mol/m2/s
! -------------------------------------
!
           fld(ji,jj)   = zpa*chemc(ji,jj,1)*1.e3*zkgco2(ji,jj)*tmask(ji,jj,1)
           flu(ji,jj)   = h2co3(ji,jj)      *1.e3*zkgco2(ji,jj)*tmask(ji,jj,1)
           dpco2(ji,jj) = min(h2co3(ji,jj)/chemc(ji,jj,1)-zpa,4000.)
!
! add tendency in mol/L/s
! -----------------------
!
           tra(ji,jj,1,jpdic) = tra(ji,jj,1,jpdic)+(fld(ji,jj)-flu(ji,jj))/1000./fse3t(ji,jj,1)
#  if defined key_trc_piic
           cflx(ji,jj) = (278.*(1.-zvapor)*chemc(ji,jj,1)               &
             -pih2co3(ji,jj))*zkgco2(ji,jj)*tmask(ji,jj,1)
           tra(ji,jj,1,jppiic) = tra(ji,jj,1,jppiic)+cflx(ji,jj)/fse3t(ji,jj,1)
           cflx(ji,jj) = cflx(ji,jj)*1.e3
#  endif
!
! calculate cumulative flux mol/timestep
           qcumul(jpdic)=qcumul(jpdic)+(fld(ji,jj)-flu(ji,jj))*rfact*e1t(ji,jj) &
     &         *e2t(ji,jj)
!
! Compute O2 flux 
! ---------------
!
          oxy16 = trn(ji,jj,1,jpoxy)
          flu16(ji,jj) = (atcox*chemc(ji,jj,3)*(1.-zvapor) - oxy16) &
     &      * zkgo2(ji,jj) * tmask(ji,jj,1)
          tra(ji,jj,1,jpoxy) = tra(ji,jj,1,jpoxy)+flu16(ji,jj)/fse3t(ji,jj,1)
          qcumul(jpoxy) = qcumul(jpoxy)+flu16(ji,jj)*rfact*e1t(ji,jj)*e2t(ji,jj)
#  if defined key_trc_dms
!
! Compute DMS flux
! ----------------
!
          fludms(ji,jj) = -trn(ji,jj,1,jpdms)*zkgdms(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpdms) = tra(ji,jj,1,jpdms)+fludms(ji,jj)/fse3t(ji,jj,1)
!
! calculate cumulative flux 
          qcumul(jpdms) = qcumul(jpdms)+fludms(ji,jj)*1.e3*rfact*e1t(ji,jj)      &
     &      *e2t(ji,jj)*tmask(ji,jj,1)
#  endif
#  if defined key_trc_n2o
!
! Compute N2O flux
! ----------------
!
          flun2s(ji,jj) = (pn2oat*1e-9*chemc(ji,jj,4)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpn2s))*zkgn2o(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpn2s) = tra(ji,jj,1,jpn2s)+flun2s(ji,jj)/fse3t(ji,jj,1)
          dpn2s(ji,jj) = trn(ji,jj,1,jpn2s)/chemc(ji,jj,4)-pn2oat*1e-9*(1.-zvapor)
#  endif
#  if defined key_trc_ch4
          fluch1(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch1))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch1) = tra(ji,jj,1,jpch1)+fluch1(ji,jj)/fse3t(ji,jj,1)
          fluch2(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch2))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch2) = tra(ji,jj,1,jpch2)+fluch2(ji,jj)/fse3t(ji,jj,1)
          fluch3(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch3))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch3) = tra(ji,jj,1,jpch3)+fluch3(ji,jj)/fse3t(ji,jj,1)
          fluch4(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch4))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch4) = tra(ji,jj,1,jpch4)+fluch4(ji,jj)/fse3t(ji,jj,1)
          fluch5(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch5))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch5) = tra(ji,jj,1,jpch5)+fluch5(ji,jj)/fse3t(ji,jj,1)
          dpch4(ji,jj)=trn(ji,jj,1,jpch2)/chemc(ji,jj,5)-pch4at*1e-9*(1.-zvapor)
          fluch6(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch6))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch6) = tra(ji,jj,1,jpch6)+fluch6(ji,jj)/fse3t(ji,jj,1)
          fluch7(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch7))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch7) = tra(ji,jj,1,jpch7)+fluch7(ji,jj)/fse3t(ji,jj,1)
          fluch8(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch8))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch8) = tra(ji,jj,1,jpch8)+fluch8(ji,jj)/fse3t(ji,jj,1)
          fluch9(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch9))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch9) = tra(ji,jj,1,jpch9)+fluch9(ji,jj)/fse3t(ji,jj,1)
          fluch10(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch10))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch10) = tra(ji,jj,1,jpch10)+fluch10(ji,jj)/fse3t(ji,jj,1)
          fluch11(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch11))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch11) = tra(ji,jj,1,jpch11)+fluch11(ji,jj)/fse3t(ji,jj,1)
          fluch12(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch12))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch12) = tra(ji,jj,1,jpch12)+fluch12(ji,jj)/fse3t(ji,jj,1)
          fluch13(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch13))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch13) = tra(ji,jj,1,jpch13)+fluch13(ji,jj)/fse3t(ji,jj,1)
          fluch14(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch14))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch14) = tra(ji,jj,1,jpch14)+fluch14(ji,jj)/fse3t(ji,jj,1)
          fluch15(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch15))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch15) = tra(ji,jj,1,jpch15)+fluch15(ji,jj)/fse3t(ji,jj,1)
          fluch16(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch16))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch16) = tra(ji,jj,1,jpch16)+fluch16(ji,jj)/fse3t(ji,jj,1)
          fluch17(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch17))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch17) = tra(ji,jj,1,jpch17)+fluch17(ji,jj)/fse3t(ji,jj,1)
          fluch18(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch18))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch18) = tra(ji,jj,1,jpch18)+fluch18(ji,jj)/fse3t(ji,jj,1)
          fluch19(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch19))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch19) = tra(ji,jj,1,jpch19)+fluch19(ji,jj)/fse3t(ji,jj,1)
          fluch20(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch20))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch20) = tra(ji,jj,1,jpch20)+fluch20(ji,jj)/fse3t(ji,jj,1)
          fluch21(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch21))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch21) = tra(ji,jj,1,jpch21)+fluch21(ji,jj)/fse3t(ji,jj,1)
          fluch22(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch22))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch22) = tra(ji,jj,1,jpch22)+fluch22(ji,jj)/fse3t(ji,jj,1)
          fluch23(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch23))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch23) = tra(ji,jj,1,jpch23)+fluch23(ji,jj)/fse3t(ji,jj,1)
          fluch24(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch24))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch24) = tra(ji,jj,1,jpch24)+fluch24(ji,jj)/fse3t(ji,jj,1)
          fluch25(ji,jj) = (pch4at*1e-9*chemc(ji,jj,5)*(1.-zvapor)-        &
     &      trn(ji,jj,1,jpch25))*zkgch4(ji,jj)*tmask(ji,jj,1)
          tra(ji,jj,1,jpch25) = tra(ji,jj,1,jpch25)+fluch25(ji,jj)/fse3t(ji,jj,1)
#  endif
#  if defined key_trc_cfc11
          hemisp = min(max((gphit(ji,jj)+10.)/20.,0.),1.)
          pcfc11 = (sipcfc(yearrc,1,1)*hemisp+sipcfc(yearrc,1,2)*(1.-hemisp))*1e-12
           tra(ji,jj,1,jpc11) = tra(ji,jj,1,jpc11)+(pcfc11*chemc(ji,jj,6)&
     &       *(1.-zvapor)-trn(ji,jj,1,jpc11))*zkgcfc11(ji,jj)*tmask(ji,jj,1)/fse3t(ji,jj,1)
#  endif
        END DO
      END DO
!
!! ----------------
!! Save diagnostics
!! ----------------
!
#  if defined key_trc_diaadd
#   if defined key_iomput
#    if defined key_trc_piic
       CALL iom_put("PICflx", cflx )
#    endif
       cflx(2:nlci-1,2:nlcj-1)    = (fld(2:nlci-1,2:nlcj-1)-flu(2:nlci-1,2:nlcj-1))
       pco2(2:nlci-1,2:nlcj-1)    = min(h2co3(2:nlci-1,2:nlcj-1)/chemc(2:nlci-1,2:nlcj-1,1),4000.)
       flu16(2:nlci-1,2:nlcj-1)   = flu16(2:nlci-1,2:nlcj-1)*1000.
#    if defined key_trc_dms
       fludms(2:nlci-1,2:nlcj-1)  = fludms(2:nlci-1,2:nlcj-1)*1000.
#    endif
#    if defined key_trc_n2o
       flun2s(2:nlci-1,2:nlcj-1)  = flun2s(2:nlci-1,2:nlcj-1)*1000.
#    endif
#    if defined key_trc_ch4
       fluch1(2:nlci-1,2:nlcj-1)  = fluch1(2:nlci-1,2:nlcj-1)*1000.
       fluch2(2:nlci-1,2:nlcj-1)  = fluch2(2:nlci-1,2:nlcj-1)*1000.
       fluch3(2:nlci-1,2:nlcj-1)  = fluch3(2:nlci-1,2:nlcj-1)*1000.
       fluch4(2:nlci-1,2:nlcj-1)  = fluch4(2:nlci-1,2:nlcj-1)*1000.
       fluch5(2:nlci-1,2:nlcj-1)  = fluch5(2:nlci-1,2:nlcj-1)*1000.
       fluch6(2:nlci-1,2:nlcj-1)  = fluch6(2:nlci-1,2:nlcj-1)*1000.
       fluch7(2:nlci-1,2:nlcj-1)  = fluch7(2:nlci-1,2:nlcj-1)*1000.
       fluch8(2:nlci-1,2:nlcj-1)  = fluch8(2:nlci-1,2:nlcj-1)*1000.
       fluch9(2:nlci-1,2:nlcj-1)  = fluch9(2:nlci-1,2:nlcj-1)*1000.
       fluch10(2:nlci-1,2:nlcj-1) = fluch10(2:nlci-1,2:nlcj-1)*1000.
       fluch11(2:nlci-1,2:nlcj-1) = fluch11(2:nlci-1,2:nlcj-1)*1000.
       fluch12(2:nlci-1,2:nlcj-1) = fluch12(2:nlci-1,2:nlcj-1)*1000.
       fluch13(2:nlci-1,2:nlcj-1) = fluch13(2:nlci-1,2:nlcj-1)*1000.
       fluch14(2:nlci-1,2:nlcj-1) = fluch14(2:nlci-1,2:nlcj-1)*1000.
       fluch15(2:nlci-1,2:nlcj-1) = fluch15(2:nlci-1,2:nlcj-1)*1000.
       fluch16(2:nlci-1,2:nlcj-1) = fluch16(2:nlci-1,2:nlcj-1)*1000.
       fluch17(2:nlci-1,2:nlcj-1) = fluch17(2:nlci-1,2:nlcj-1)*1000.
       fluch18(2:nlci-1,2:nlcj-1) = fluch18(2:nlci-1,2:nlcj-1)*1000.
       fluch19(2:nlci-1,2:nlcj-1) = fluch19(2:nlci-1,2:nlcj-1)*1000.
       fluch20(2:nlci-1,2:nlcj-1) = fluch20(2:nlci-1,2:nlcj-1)*1000.
       fluch21(2:nlci-1,2:nlcj-1) = fluch21(2:nlci-1,2:nlcj-1)*1000.
       fluch22(2:nlci-1,2:nlcj-1) = fluch22(2:nlci-1,2:nlcj-1)*1000.
       fluch23(2:nlci-1,2:nlcj-1) = fluch23(2:nlci-1,2:nlcj-1)*1000.
       fluch24(2:nlci-1,2:nlcj-1) = fluch24(2:nlci-1,2:nlcj-1)*1000.
       fluch25(2:nlci-1,2:nlcj-1) = fluch25(2:nlci-1,2:nlcj-1)*1000.
#    endif
!
       WHERE (tmask(:,:,1) .EQ. 0 )
         cflx = ncf_fill
         dpco2 = ncf_fill
         flu16 = ncf_fill
#    if defined key_trc_dms
         fludms = ncf_fill
#    endif
#    if defined key_trc_n2o
         flun2s = ncf_fill
#    endif
#    if defined key_trc_ch4
         fluch4 = ncf_fill
#    endif
       END WHERE
!
       diaadd="Cflx"
       CALL iom_put(diaadd, cflx(:,:) )
       diaadd="dpCO2"
       CALL iom_put(diaadd, dpco2(:,:) )
       CALL iom_put("pCO2", pco2(:,:) )
       diaadd="Oflx"
       CALL iom_put(diaadd, flu16(:,:) )
#    if defined key_trc_dms
       diaadd="DMSflux"
       CALL iom_put(diaadd, fludms(:,:) )
#    endif
#    if defined key_trc_n2o
       CALL iom_put("N2Sflux", flun2s(:,:) )
       CALL iom_put("dpN2S", dpn2s)
#    endif
#    if defined key_trc_ch4
       CALL iom_put("CH1flux", fluch1(:,:) )
       CALL iom_put("CH2flux", fluch2(:,:) )
       CALL iom_put("CH3flux", fluch3(:,:) )
       CALL iom_put("CH4flux", fluch4(:,:) )
       CALL iom_put("CH5flux", fluch5(:,:) )
       CALL iom_put("dpCH2",   dpch4)
       CALL iom_put("CH6flux", fluch6(:,:) )
       CALL iom_put("CH7flux", fluch7(:,:) )
       CALL iom_put("CH8flux", fluch8(:,:) )
       CALL iom_put("CH9flux", fluch9(:,:) )
       CALL iom_put("CH10flux",fluch10(:,:) )
       CALL iom_put("CH11flux",fluch11(:,:) )
       CALL iom_put("CH12flux",fluch12(:,:) )
       CALL iom_put("CH13flux",fluch13(:,:) )
       CALL iom_put("CH14flux",fluch14(:,:) )
       CALL iom_put("CH15flux",fluch15(:,:) )
       CALL iom_put("CH16flux",fluch16(:,:) )
       CALL iom_put("CH17flux",fluch17(:,:) )
       CALL iom_put("CH18flux",fluch18(:,:) )
       CALL iom_put("CH19flux",fluch19(:,:) )
       CALL iom_put("CH20flux",fluch20(:,:) )
       CALL iom_put("CH21flux",fluch21(:,:) )
       CALL iom_put("CH22flux",fluch22(:,:) )
       CALL iom_put("CH23flux",fluch23(:,:) )
       CALL iom_put("CH24flux",fluch24(:,:) )
       CALL iom_put("CH25flux",fluch25(:,:) )
#    endif
#   endif
#  endif
#endif
      RETURN
      END

