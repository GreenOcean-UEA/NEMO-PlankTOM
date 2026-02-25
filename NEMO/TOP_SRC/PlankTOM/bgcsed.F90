      SUBROUTINE bgcsed(kt)
#if defined key_planktom && defined key_top
!!======================================================================
!!                         ***  ROUTINE bgcsed  ***
!! TOP : PlankTOM ecosystem model 
!!======================================================================
!! History : PISCES       ! 2001 (O. Aumont) Original code 
!!           PlankTOM see bgcbio
!!           
!!----------------------------------------------------------------------
!! ** Purpose : compute loss of organic matter in the sediments necessary to balance
!!          input of elements from rivers and dust. 
!!          This is by no way a sediment model. Sediment biogeochemistry is
!!          computed at the end of subroutine bgcsnk.F90 
!!
!! ** Action  : 
!! ** External : NONE
!!
!! References : PlankTOM12 manual, Buitenhuis et al. 2023
!!                                 https://zenodo.org/records/8388158 
!! ---------------------------------------------------------------------
!
! parameters and commons
      USE trp_trc
      USE sms_planktom
      USE oce_trc
      USE lib_mpp
      USE lbclnk
      USE dom_oce , ONLY :   mbathy     =>   mbathy     !: number of ocean level (=0,  & 1, ... , jpk-1) 
      IMPLICIT NONE
! 
! local variables
      INTEGER ji, jj, jn, ikt, kt
      REAL dummyv(5),sumsed(12),botcor(10),litres,ztemp
#  include "domzgr_substitute.h90"
!!
!! ------------------------------------------------------------
!! Loss of biogenic silicon and organic carbon in the sediments 
!! ------------------------------------------------------------
!!
!    First, the total inventory is computed
!    --------------------------------------
!
      sumsed=0.
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
          ikt = max(mbathy(ji,jj)-1,1)
          litres = 1000.*e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,ikt)
!
!         Si budget
          sumsed(1) = sumsed(1)+trnsed(ji,jj,jpdsi)*litres
!
!         P budget
          sumsed(2) = sumsed(2)+trnsed(ji,jj,jpgoc)*litres
          sumsed(7) = sumsed(7)+(trnsed(ji,jj,jppoc)+trnsed(ji,jj,jphoc))*litres
!
!         Fe budget
          sumsed(3) = sumsed(3)+(trnsed(ji,jj,jpbfe)+trnsed(ji,jj,jpufe))*litres
          sumsed(8) = sumsed(8)+trnsed(ji,jj,jpsfe)*litres
!
!         DOC budget
          sumsed(4) = sumsed(4)+trn(ji,jj,ikt,jpdoc)*litres
!
!         CaCO3 budget
          sumsed(5)  = sumsed(5)+trnsed(ji,jj,jpcal)*litres
          sumsed(10) = sumsed(10)+trnsed(ji,jj,jpara)*litres
!
!         N budget
          sumsed(11) = sumsed(11)+trnsed(ji,jj,jpgon)*litres
!
!         C budget
          sumsed(12) = sumsed(12)+trn(ji,jj,ikt,jpdic)*litres
        END DO
      END DO
      IF( lk_mpp ) CALL mpp_sum( sumsed,12)
!
!    Then the surface inputs of each element are scaled against the 
!    sediment inventories, and the bottom concentrations are multiplied
!    by this scaling factor. Thus, the amount lost in the sediments equals
!    the supply at the surface (dust+rivers)
!    -------------------------------------------------------------------
!
!    loop over dSi, totalPOC, totalPOFe, DOC, totalCaCO3
      DO jn = 1, 5
!       
        dummyv(jn) = sumsed(jn)+sumsed(jn+5)
!
!    if the quantity it wants to remove to balance surface inputs (extinp)
!    is less than what is present in the sediment (dummyv) 
        IF (extinp(jn) .LE. dummyv(jn)) THEN
!
!    then remove extinp from the sediment
          botcor(jn)  = 1.-extinp(jn)/dummyv(jn)
        ELSE
!
!    else remove everything there is in the sediment
          botcor(jn)  = 0.
!
!    and keep track of how much it was unable to remove
          sedcor(jn)  = sedcor(jn)+extinp(jn)-dummyv(jn)
        ENDIF
      END DO
!
!    if the removal of P (variable 2) is less than what you need to
!    remove for N, then remove additional N
      ztemp = extinp(6)-sumsed(7)*(1.-botcor(2))*ratn2c
      IF (ztemp .LE. sumsed(11)) THEN
        botcor(6)  = 1.-max(ztemp,0.)/sumsed(11)
        sedcor(6)  = sedcor(6)+min(ztemp,0.)
      ELSE
        botcor(6)  = 0.
        sedcor(6)  = sedcor(6)+ztemp-sumsed(11)
      ENDIF
!
!    if the removal of P (variable 2) is less than what you need to
!    remove for C, then remove additional C
      IF (extinp(7) .LE. sumsed(12)) THEN
        botcor(7)  = 1.-extinp(7)/sumsed(12)
      ELSE
        botcor(7)  = 0.
        sedcor(7)  = sedcor(7)+extinp(7)-sumsed(12)
      ENDIF
!
!     only output the sediment inventory at the beginning and end of the run
      IF ( kt.EQ.nit000 .OR. kt.EQ.nitend ) THEN
        IF (lwp) THEN
          WRITE(numout,31) "sediment inventory (Tmol):", (dummyv(jn)*1.e-12, jn=1,5),sumsed(11)*1.e-12,sumsed(12)*1.e-12
31        FORMAT(a26,7f11.3)
        END IF
      END IF
!
!    Remove material from the sediment
!    ---------------------------------
!
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
          ikt = max(mbathy(ji,jj)-1,1)
          trnsed(ji,jj,jpdsi)  = trnsed(ji,jj,jpdsi)*botcor(1)
          trnsed(ji,jj,jpgoc)  = trnsed(ji,jj,jpgoc)*botcor(2)
          trnsed(ji,jj,jphoc)  = trnsed(ji,jj,jphoc)*botcor(2)
          trnsed(ji,jj,jppoc)  = trnsed(ji,jj,jppoc)*botcor(2)
          trnsed(ji,jj,jpbfe)  = trnsed(ji,jj,jpbfe)*botcor(3)
          trnsed(ji,jj,jpufe)  = trnsed(ji,jj,jpufe)*botcor(3)
          trnsed(ji,jj,jpsfe)  = trnsed(ji,jj,jpsfe)*botcor(3)
          trn(ji,jj,ikt,jpdoc) = trn(ji,jj,ikt,jpdoc)*botcor(4)
          trnsed(ji,jj,jpcal)  = trnsed(ji,jj,jpcal)*botcor(5)
          trnsed(ji,jj,jpara)  = trnsed(ji,jj,jpara)*botcor(5)
          trnsed(ji,jj,jpgon)  = trnsed(ji,jj,jpgon)*botcor(6)
          trn(ji,jj,ikt,jpdic) = trn(ji,jj,ikt,jpdic)*botcor(7)
        END DO
      END DO
#endif
      RETURN
      END
