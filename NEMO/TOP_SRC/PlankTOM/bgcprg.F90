       SUBROUTINE bgcprg(kt)
#if defined key_planktom &&  defined key_top
!!======================================================================
!!                         ***  ROUTINE bgcprg  ***
!! TOP : PlankTOM ecosystem model
!!======================================================================
!! History : original  : E. Maier-Reimer (GBC 1993): h3cprg
!!           additions : O. Aumont (1998)
!!           additions : C. Le Quere (1999)
!!           additions : O. Aumont (2001)
!!           additions : O. Aumont , EK (05/2001): add h3cadj 
!!           see also bgcbio 
!!----------------------------------------------------------------------
!! ** Purpose : Call Biological sources and sinks subroutines
!! ** Action  : - calls : bgcche
!!                        bgcint
!!                        bgclys
!!                        bgcbio
!!                        bgcsed
!!                        bgcflx
!!
!! References : PlankTOM12 manual, Buitenhuis et al. 2023
!!                                 https://zenodo.org/records/8388158
!! ---------------------------------------------------------------------
!
! parameters and commons
      USE trp_trc
      USE sms_planktom
      USE oce_trc
      USE lbclnk
      USE lib_mpp
      USE iom
      IMPLICIT NONE
! 
! local variables
      INTEGER ktask, kt, jn
      INTEGER ji, jj, jk
      REAL(wp) total(jptra), totaf(jptra), totic, totif, totec, totef, units
!
! Initialisation of variables and calculation of inventories
! ----------------------------------------------------------
!
! If no zpft adds to HOC, set HOC to 0 so that the results are the same as the model without HOC
      IF ((kt .EQ. nit000) .and. (ngochoc .LE. 1)) THEN
        trn(:,:,:,jphoc) = 0.
        trb(:,:,:,jphoc) = 0.
        trn(:,:,:,jpufe) = 0.
        trb(:,:,:,jpufe) = 0.
      ENDIF
!
      IF (kt .EQ. nit000) THEN
        qcumul = 0.
      END IF
!
      total = 0.
      DO jn = 1, jptra
        DO jk = 1, jpk
          DO jj = 2, nlcj-1
            DO ji = 2, nlci-1
              total(jn) = total(jn)+trn(ji,jj,jk,jn)*volumt(ji,jj,jk)
            END DO
          END DO
        END DO
      END DO
#  if defined key_mpp_mpi
      CALL mpp_sum(total,jptra)
#  endif
!
      totic = total(jpdoc)
      DO jn = jppoc, jpdia+jppft-1
        totic = totic+total(jn)
      END DO
!
      totif = total(jpsfe)+total(jpbfe)
      DO jn = jpbac, jpbac+jpzft
        totif = totif+total(jn)*ferat3
      END DO
!
      DO jn = jpdfe, jpdfe+jppft-1
        totif = totif+total(jn)
      END DO
!
      IF (kt .EQ. nit000) THEN
        IF (lwp) THEN
          WRITE(numout,*) ''    
          WRITE(numout,*) 'total inventory by element'
          WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~'
          WRITE(numout,*) 'bgcprg start C  ',totic+total(jpdic)
          WRITE(numout,*) 'bgcprg start P  ',totic+total(jppo4) 
          WRITE(numout,*) 'bgcprg start Fe ',totif+total(jpfer)
          WRITE(numout,*) 'bgcprg start Si ',total(jpsil)+total(jpbsi)+total(jpdsi) 
          WRITE(numout,*) 'bgcprg start O2 ',total(jpoxy)-rato2c*totic

          WRITE(numout,*) '' 
          WRITE(numout,*) 'total inventory by tracer'
          WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~'
          WRITE(numout,*) 'bgcprg start jptra'
          DO jn = 1, jptra
            WRITE(numout,*) '      tracer start ',ctrcnm(jn),total(jn)
          END DO
        END IF
      END IF
!
! Interpolate chemical variables
! ------------------------------
!
      CALL bgcint(kt)
!
! Compute chemical variables
! --------------------------
!
      CALL bgcche
!
! Compute CaCO3 saturation
! ------------------------
!
      CALL bgclys
!
! Compute biology (ecosystem sources and sinks)
! ---------------------------------------------
!
      CALL bgcbio(kt)
!
! Close budgets
! -------------
!
      CALL bgcsed(kt)
!
      IF ( kt == nitend ) THEN
#  if defined key_mpp_mpi
        CALL mpp_sum(sedcor,6)
#  endif
!
! print the applied correction 
        units = rfact/raass*100.
        IF(lwp) WRITE(numout,*)' bottom water correction lacked ',     &
     &        sedcor,' mol Si, C, Fe, DOC, Alk, N, which is ',         &
     &      (sedcor(ji)/extinp(ji)*units,ji=1,6),'% of external inputs.'
#    if defined key_trc_dms
        IF(lwp) WRITE(numout,*)' dms sinks converted to sources ',     &
     &          dms_snkcount, ' times'
        IF(lwp) WRITE(numout,*)' dmd sinks converted to sources ',     &
     &          dmd_snkcount, ' times'
#    endif
      ENDIF
!
! Compute surface fluxes
! ----------------------
!
      CALL bgcflx
!
! calculate changes in total passive tracer inventories, to detect coding errors
      totaf = 0.
      DO jn = 1, jptra
        DO jk = 1, jpk-1
          DO jj = 2, nlcj-1
            DO ji = 2, nlci-1
              totaf(jn) = totaf(jn)+trn(ji,jj,jk,jn)*volumt(ji,jj,jk)
            END DO
          END DO
        END DO
      END DO
!
#  if defined key_mpp_mpi
      CALL mpp_sum(totaf,jptra)
#  endif
!
      DO jn = 1, jptra
        qcumul(jn) = qcumul(jn)+totaf(jn)-total(jn)
      END DO
!
      IF (kt .EQ. nitend) THEN
!
        totec = totaf(jpdoc)
        DO jn = jppoc, jpdia+jppft-1
          totec = totec+totaf(jn)
        END DO
!
        totef = total(jpsfe)+total(jpbfe)
        DO jn = jpbac, jpbac+jpzft
          totef = totef+totaf(jn)*ferat3
        END DO
!
        DO jn = jpdfe, jpdfe+jppft-1
          totef = totef+totaf(jn)
        END DO
!
        IF (lwp) THEN
          WRITE(numout,*) '' 
          WRITE(numout,*) 'total inventory by element'
          WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~'
          WRITE(numout,*) 'bgcprg end C  ',totec+totaf(jpdic)
          WRITE(numout,*) 'bgcprg end P  ',totec+totaf(jppo4)
          WRITE(numout,*) 'bgcprg end Fe ',totef+totaf(jpfer)
          WRITE(numout,*) 'bgcprg end Si ',totaf(jpsil)+totaf(jpbsi)+totaf(jpdsi)
          WRITE(numout,*) 'bgcprg end O2 ',totaf(jpoxy)-rato2c*totec

          WRITE(numout,*) '' 
          WRITE(numout,*) 'total inventory by tracer'
          WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~'
          WRITE(numout,*) 'bgcprg end jptra'
          DO jn = 1, jptra
            WRITE(numout,*) '      tracer end ',ctrcnm(jn),totaf(jn)
          END DO
        END IF
!
        totec = qcumul(jpdoc)
        DO jn = jppoc, jpdia+jppft-1
          totec = totec+qcumul(jn)
        END DO
!
        totef = total(jpsfe)+total(jpbfe)
        DO jn = jpbac, jpbac+jpzft
          totef = totef+qcumul(jn)*ferat3
        END DO
!
        DO jn = jpdfe, jpdfe+jppft-1
          totef = totef+qcumul(jn)
        END DO
!
        IF (lwp) THEN
          WRITE(numout,*) '' 
          WRITE(numout,*) 'total budget by element'
          WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~'
          WRITE(numout,*) 'bgcprg budget C  ',totec+qcumul(jpdic)
          WRITE(numout,*) 'bgcprg budget P  ',totec+qcumul(jppo4)
          WRITE(numout,*) 'bgcprg budget Fe ',totef+qcumul(jpfer)
          WRITE(numout,*) 'bgcprg budget Si ',qcumul(jpsil)+qcumul(jpbsi)+qcumul(jpdsi)
          WRITE(numout,*) 'bgcprg budget O2 ',qcumul(jpoxy)-rato2c*totec

          WRITE(numout,*) '' 
          WRITE(numout,*) 'total budget by tracer'
          WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~'
          WRITE(numout,*) 'bgcprg budget jptra'
          DO jn = 1, jptra
            WRITE(numout,*) '      tracer budget ',ctrcnm(jn),qcumul(jn)
          END DO
        END IF
!
      END IF
!
      DO jn = 1, jptra
        CALL lbc_lnk(trn(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk(tra(:,:,:,jn), 'T', 1. )
      END DO
      trb(:,:,:,:) = trn(:,:,:,:)
#endif
!
      RETURN
      END

