#if defined key_planktom && defined key_top
      SUBROUTINE bgcpro(kt)
!!======================================================================
!!                         ***  ROUTINE bgcpro  ***
!! TOP : PlankTOM ecosystem model 
!!======================================================================
!! History : original p4zprod       ! 2002 (O. Aumont) 
!!           PlankTOM see bgcbio
!!----------------------------------------------------------------------
!! ** Purpose : compute the phytoplankton production depending on
!!              light, temperature and nutrient availability for all
!!              Plankton Functional Types (PFTs)
!! ** Action  :  
!! ** External  : etot used here is calculated in traqsr.F90
!!
!! References : Manizza et al. 2005 doi:10.1029/2004GL020778
!!              Buitenhuis and Geider, L&O 2010
!!              Le Quéré et al., Biogeosciences, 2016
!!              PlankTOM12 manual, Buitenhuis et al. 2023
!!                                 https://zenodo.org/records/8388158 
!!----------------------------------------------------------------------
!
! parameters and commons
      USE trc
      USE trp_trc
      USE sms_planktom
      USE oce_trc
      USE traqsr  , ONLY : ln_qsr_sms
      USE iom
      IMPLICIT NONE
! 
! local variables
      INTEGER ji, jj, jk, jl, kt
      REAL silpot1,silpot2,silfac,pislopen,ysopt
      REAL xlimb1 ! PO4 limitation of bacteria
      REAL xlimb2 ! iron limitation of bacteria
      REAL xlimb3 ! DOC limitation of bacteria
!
      REAL xlim1                      ! iron quota limitation of maximumm iron uptake rate
      REAL xlim2(jpdia:jpdia+jppft-1) ! iron uptake limitation 
      REAL xlim3                      ! iron quota limitation of photosynthesis
      REAL xlim4(jpdia:jpdia+jppft-1) ! PO4 limitation
      REAL xlim5(jpdia:jpdia+jppft-1) ! silicate limitation
      REAL xlim6(jpdia:jpdia+jppft-1) ! NO3 limitation
      REAL xlim8                      ! light limitation
!
      REAL parlux,xchl,ekg,ekr
      REAL dinlim
! 
      REAL pcphot(jpdia:jpdia+jppft-1),quopfe(jpdia:jpdia+jppft-1)
      REAL rhochl(jpdia:jpdia+jppft-1),vcfer(jpdia:jpdia+jppft-1)
      REAL perfrm(jpdia:jpdia+jppft-1)
      REAL pctnut,docpro
      REAL zmask,zkrdphy
!
      LOGICAL ln_lop  !  calculate limits of phytoplankton
      NAMELIST/namtrc_lop/ln_lop
!
#  include "domzgr_substitute.h90"
!
      REWIND( numnat_ref )  !  namtrc_lop in reference namelist : ln_lop
      READ  ( numnat_ref, namtrc_lop )

      ! output on first timestep only
      IF ( kt .EQ. nit000 ) THEN
        IF (lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) 'limits of phytoplankton'
          WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~'

          IF(ln_lop) THEN
            WRITE(numout,*) 'ln_lop = TRUE : calculating limits of phytoplankton'
          ELSE
            WRITE(numout,*) 'ln_lop = FALSE : do not calculate limits of phytoplankton'
          END IF

          WRITE(numout,*)
        END IF
      END IF
!
! Initialisation of variables
! ----------------------------
!
! calculate the maximum photosynthesis rate
      DO jl = jpdia, jpdia+jppft-1
        pcmax(jl) = rn_mumpft(jl)*(1.+rn_resphy(jl))/rjjss
      END DO
!
! set silicate limitation to 1 for non-diatoms
      DO jl = jpdia+1, jpdia+jppft-1
        xlim5(jl) = 1.
      END DO
!
      DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
!!
!! -------------------------------------------------------------------
!! Calculate growth rate for all PFTs. First calculate the nutrient 
!! limitation, then calculate the iron-light limitation
!! Note that Diatoms and N2 fixers have specific requirements
!! -------------------------------------------------------------------
!!
!       Michaelis-Menten Limitation term for nutrients
!       ----------------------------------------------
!
! Diatoms (DIA)
! -------------
!
            jl = jpdia
!
! michaelis menten P
            xlim4(jl) = trn(ji,jj,jk,jppo4)/(trn(ji,jj,jk,jppo4)+rn_kmpphy(jl))
!
! michaelis menten Si
            xlim5(jl) = trn(ji,jj,jk,jpsil)/(trn(ji,jj,jk,jpsil)+rn_sildia)
!
! michaelis menten N
            xlim6(jl) = trn(ji,jj,jk,jpdin)/(trn(ji,jj,jk,jpdin)+rn_kmnphy(jl))
!
            IF(ln_lop) THEN
              lim4po4(ji,jj,jk,jl-jpcru) = xlim4(jl)
              lim5si(ji,jj,jk,jl-jpcru)  = xlim5(jl)
              lim6din(ji,jj,jk,jl-jpcru) = xlim6(jl)
            END IF
!
! MIX, COC, PIC, PHA
! ------------------
!
            DO jl = jpdia+1, jpdia+jppft-2
!
! michaelis menten P
              xlim4(jl) = trn(ji,jj,jk,jppo4)/(trn(ji,jj,jk,jppo4)+rn_kmpphy(jl))
!
! michaelis menten N
              xlim6(jl) = trn(ji,jj,jk,jpdin)/(trn(ji,jj,jk,jpdin)+rn_kmnphy(jl))
!
              IF(ln_lop) THEN
                lim4po4(ji,jj,jk,jl-jpcru) = xlim4(jl)
                lim6din(ji,jj,jk,jl-jpcru) = xlim6(jl)
              END IF
            END DO
!
! Michaelis-Menten Limitation term for nutrients for bacteria -
! for use in degradation of DMS
            xlimb1 = trn(ji,jj,jk,jppo4)/(trn(ji,jj,jk,jppo4)+rn_kmpbac)
            xlimb2 = trn(ji,jj,jk,jpfer)/(trn(ji,jj,jk,jpfer)+rn_kmfbac)
            xlimb3 = trn(ji,jj,jk,jpdoc)/(trn(ji,jj,jk,jpdoc)+rn_kmobac)
            xlimbac(ji,jj,jk) = min(xlimb1,xlimb2,xlimb3)
!
! N2 fixers (FIX)
! ---------------
!
            jl = jpfix
!
! michaelis menten P
            xlim4(jl) = trn(ji,jj,jk,jppo4)/(trn(ji,jj,jk,jppo4)+rn_kmpphy(jl))
            dinlim    = trn(ji,jj,jk,jpdin)/(trn(ji,jj,jk,jpdin)+rn_kmnphy(jl))
!
! michaelis menten N
            xlim6(jl) = dinlim +rn_munfix*(1.-dinlim)
            dinpft(ji,jj,jk,jl) = dinlim/(xlim6(jl)+rtrn)*tmask(ji,jj,jk)
!
            IF(ln_lop) THEN
              lim4po4(ji,jj,jk,jl-jpcru) = xlim4(jl)
              lim6din(ji,jj,jk,jl-jpcru) = xlim6(jl)
            END IF 
!
! All phytoplankton PFTs
! ----------------------
!
            DO jl = jpdia, jpdia+jppft-1
!
!       quota model for Fe-light based on Buitenhuis and Geider, L&O 2010 
!       -----------------------------------------------------------------
!
              stofoo(ji,jj,jk,jl,2) = trn(ji,jj,jk,jl+jppft)          &
     &          /(trn(ji,jj,jk,jl)+rtrn)
              stofoo(ji,jj,jk,jl,3) = trn(ji,jj,jk,jl+2*jppft)        &
     &          /(trn(ji,jj,jk,jl)+rtrn)
              quopfe(jl) = max(min(stofoo(ji,jj,jk,jl,2),rn_qmaphy(jl)),rn_qmiphy(jl))
! 
              xlim1 = (rn_rhfphy(jl)*rn_qmaphy(jl)-rn_qmaphy(jl))*    &
     &          (rn_qmaphy(jl)-quopfe(jl))/                           &
     &          (rn_qmaphy(jl)-rn_qmiphy(jl))+rn_qmaphy(jl)
              xlim2(jl) = trn(ji,jj,jk,jpfer)/(trn(ji,jj,jk,jpfer)+       &
     &          rn_kmfphy(jl))
              xlim3 = min((quopfe(jl)-rn_qmiphy(jl))                   &
     &          /(rn_qopphy(jl)-rn_qmiphy(jl)),1.)
!
              IF(ln_lop) THEN
                lim3fe(ji,jj,jk,jl-jpcru) = xlim3
              END IF
!
              xlimpft(ji,jj,jk,jl) = min(xlim4(jl),xlim5(jl),xlim6(jl),xlim3)
! 
! Fe uptake rate 
              vcfer(jl) = rn_mumpft(jl)*(1.+rn_resphy(jl))            &
     &          *xlim1*min(xlim4(jl),xlim5(jl),xlim6(jl),xlim2(jl))
!
! 4.6 micromol photons/J at 550 nm
              perfrm(jl) = rn_alpphy(jl)*stofoo(ji,jj,jk,jl,3)          &
     &          *4.6*etot(ji,jj,jk)
              docpro = rn_docphy(jl)+(1.-xlimpft(ji,jj,jk,jl))*rn_domphy(jl)
              pctnut = pcmax(jl)*xlimpft(ji,jj,jk,jl)*tgfunc(ji,jj,jk,jl)

! light limitation
              xlim8      = (1.-exp(-perfrm(jl)/(pctnut+rtrn)))
              pcphot(jl) = pctnut*xlim8
              rhochl(jl) = rn_thmphy(jl)*pcphot(jl)/(perfrm(jl)+rtrn)
!
              IF(ln_lop) THEN
                lim8light(ji,jj,jk,jl-jpcru) = xlim8
              END IF
!
! synthesis rates
              prophy(ji,jj,jk,jl,3) = rhochl(jl)*pcphot(jl)           &
     &          *trn(ji,jj,jk,jl)*rfact
              prophy(ji,jj,jk,jl,2) = vcfer(jl)*xlim8*tgfunc(ji,jj,jk,jl) &
     &          *trn(ji,jj,jk,jl)*rfact/rjjss
              prophy(ji,jj,jk,jl,1) = pcphot(jl)*trn(ji,jj,jk,jl)*rfact
              docphy(ji,jj,jk,jl)   = prophy(ji,jj,jk,jl,1)*docpro
            END DO
!
            IF(ln_lop) THEN
              lim2mmfe(ji,jj,jk,1) = xlim2(jpdia)
              lim2mmfe(ji,jj,jk,2) = xlim2(jpmix)
              lim2mmfe(ji,jj,jk,3) = xlim2(jpcoc)
              lim2mmfe(ji,jj,jk,4) = xlim2(jppic)
              lim2mmfe(ji,jj,jk,5) = xlim2(jppha)
              lim2mmfe(ji,jj,jk,6) = xlim2(jpfix)
            END IF
!
!    FE/C and Si/C of diatoms
!    ------------------------
!    Si/C increases with iron stress and silicate availability
!    Si/C is increased for very high Si concentrations
!    to mimic the very high ratios observed in the Southern Ocean
!    (silpot2)
            silpot1 = 1.+rn_ferbsi*min(1.,trn(ji,jj,jk,jpsil)/rn_sildia)* &
     &        (1.-min(trn(ji,jj,jk,jpfer)/rn_kmfphy(jpdia),1.))
            silpot2 = rn_silbsi*trn(ji,jj,jk,jpsil)/(trn(ji,jj,jk,jpsil)+rn_kmsbsi)
            silfac  = max(silpot1,silpot2)
            ysopt   = rn_bsidia*silfac
            prorca3(ji,jj,jk) = prophy(ji,jj,jk,jpdia,1)*ysopt*tmask(ji,jj,jk)
          END DO
        END DO
      END DO

      IF(ln_lop) THEN
        CALL iom_put("lim2mmfe_dia",  lim2mmfe(:,:,:,1) )
        CALL iom_put("lim2mmfe_mix",  lim2mmfe(:,:,:,2) )
        CALL iom_put("lim2mmfe_coc",  lim2mmfe(:,:,:,3) )
        CALL iom_put("lim2mmfe_pic",  lim2mmfe(:,:,:,4) )
        CALL iom_put("lim2mmfe_pha",  lim2mmfe(:,:,:,5) )
        CALL iom_put("lim2mmfe_fix",  lim2mmfe(:,:,:,6) )
!
        CALL iom_put("lim3fe_dia",    lim3fe(:,:,:,1) )
        CALL iom_put("lim3fe_mix",    lim3fe(:,:,:,2) )
        CALL iom_put("lim3fe_coc",    lim3fe(:,:,:,3) )
        CALL iom_put("lim3fe_pic",    lim3fe(:,:,:,4) )
        CALL iom_put("lim3fe_pha",    lim3fe(:,:,:,5) )
        CALL iom_put("lim3fe_fix",    lim3fe(:,:,:,6) )
!
        CALL iom_put("lim4po4_dia",   lim4po4(:,:,:,1) )
        CALL iom_put("lim4po4_mix",   lim4po4(:,:,:,2) )
        CALL iom_put("lim4po4_coc",   lim4po4(:,:,:,3) )
        CALL iom_put("lim4po4_pic",   lim4po4(:,:,:,4) )
        CALL iom_put("lim4po4_pha",   lim4po4(:,:,:,5) )
        CALL iom_put("lim4po4_fix",   lim4po4(:,:,:,6) )
!
        CALL iom_put("lim5si_dia",    lim5si(:,:,:,1) )
!
        CALL iom_put("lim6din_dia",   lim6din(:,:,:,1) )
        CALL iom_put("lim6din_mix",   lim6din(:,:,:,2) )
        CALL iom_put("lim6din_coc",   lim6din(:,:,:,3) )
        CALL iom_put("lim6din_pic",   lim6din(:,:,:,4) )
        CALL iom_put("lim6din_pha",   lim6din(:,:,:,5) )
        CALL iom_put("lim6din_fix",   lim6din(:,:,:,6) )
!
        CALL iom_put("lim8light_dia", lim8light(:,:,:,1) )
        CALL iom_put("lim8light_mix", lim8light(:,:,:,2) )
        CALL iom_put("lim8light_coc", lim8light(:,:,:,3) )
        CALL iom_put("lim8light_pic", lim8light(:,:,:,4) )
        CALL iom_put("lim8light_pha", lim8light(:,:,:,5) )
        CALL iom_put("lim8light_fix", lim8light(:,:,:,6) )
      END IF

      RETURN
      END
#endif
