#if defined key_planktom && defined key_top
      SUBROUTINE bgcint(kt)
!!======================================================================
!!                         ***  ROUTINE bgcint  ***
!! TOP : PlankTOM ecosystem model 
!!======================================================================
!! History : see bgcbio
!!----------------------------------------------------------------------
!! ** Purpose : Initialise common variables at the start of PlankTOM 
!! ** Action  : 
!!              - Interpolate monthly dust fields
!!              - Calculate the temperature dependence of growth/grazing
!!              - Assign susceptibility to advection for crustaceans (code not used)
!!
!! References :  PlankTOM12 manual, Buitenhuis et al. 2023
!!                                 https://zenodo.org/records/8388158 
!!               Wiedenmann et al., Limnology and Oceanography 2009
!! ---------------------------------------------------------------------
!
! parameters and commons
      USE trp_trc
      USE sms_planktom
      USE trc
      USE oce_trc
      USE dom_oce , ONLY :   gphit     =>   gphit       !: latitude  of t-point (degre)
      USE dom_oce , ONLY :   nmonth    =>   nmonth      !: Current month
      USE sbc_oce , ONLY :   fr_i      =>   fr_i        !: ice fraction (between 0 to 1)
      IMPLICIT NONE
! 
! local variables
      INTEGER kt
      INTEGER ji, jj, jk, jl
      INTEGER nvit1t, nvit2t, icvit
      REAL    zt, zpdtmo, zdemi  
      REAL    zice1, zice2, zfacti, zfactl, zfactd, zmonth
!!
!! --------------------------- 
!! Initialisation of variables
!! ---------------------------
!! 
! assign number of readings (12 months in this case)
      icvit = 12
! 
! assign month of year and its middle point for interpolating monthly data
!
      zpdtmo = float(12*730*60*60/12)/rdttra(1)
      zdemi  = zpdtmo / 2.
!
! determine the monthly fraction
      zt     = ( float (kt) + zdemi ) / zpdtmo
!!      
!! -------------------------------
!! Interpolate monthly dust fields
!! -------------------------------
!!
! determine the relative weight of the two months 
      xtvit = zt - int (zt)
!
! determine the index of the two months needed for interpolating field 
      nvit1t = int (( float (kt) + zdemi ) / zpdtmo )
      nvit2t = nvit1t+1
      nvit1t = MOD (nvit1t-1, icvit)+1
      nvit2t = MOD (nvit2t-1, icvit)+1
!
      IF (nvit1t.eq.0) nvit1t=12 ! adjustment for lower month
!
! interpolate dust at specific time
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
          dust(ji,jj) = (1.-xtvit)*dustmo(ji,jj,nvit1t) + xtvit*dustmo(ji,jj,nvit2t)
        END DO
      END DO
!!
!! -------------------------------------------------
!! Initialisation of temperature diagnostic (biodts)
!! -------------------------------------------------
!!
      DO jk = 1, jpkm1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
            obstemp(ji,jj,jk) = (1.-xtvit)*obstempmo(ji,jj,jk,nvit1t) + xtvit*obstempmo(ji,jj,jk,nvit2t)
            obssal(ji,jj,jk)  = (1.-xtvit)*obssalmo(ji,jj,jk,nvit1t) + xtvit*obssalmo(ji,jj,jk,nvit2t)
          END DO
        END DO
      END DO
!
! Calculate the new value for tsn based on nn_biodts
!
      SELECT CASE ( nn_biodts )
        CASE ( 0 )  !  nn_biodts = 0 : using default modelled T and S'
          tsnbio = tsn 
          tsnche = tsn 

        CASE ( 1 )  !  nn_biodts = 1 : biological subroutines use obs T and S'
          tsnbio(:,:,:,1) = obstemp
          tsnbio(:,:,:,2) = obssal
          tsnche = tsn 

        CASE ( 2 )  !  nn_biodts = 2 : chemical subroutines use obs T and S'
          tsnbio = tsn 
          tsnche(:,:,:,1) = obstemp
          tsnche(:,:,:,2) = obssal

        CASE ( 3 )  !  nn_biodts = 3 : bio. and chem. subroutines use obs T and S'
          tsnbio(:,:,:,1) = obstemp
          tsnbio(:,:,:,2) = obssal
          tsnche(:,:,:,1) = obstemp
          tsnche(:,:,:,2) = obssal

      END SELECT
!!
!! -----------------------------------------------------------------------------------------
!! Calculate the temperature dependence of growth/grazing for all pft's - including bacteria
!! -----------------------------------------------------------------------------------------
!!
      DO jl = jpbac, jpdia+jppft-1
        DO jk = 1, jpkm1
          DO jj = 2, nlcj-1
            DO ji = 2, nlci-1
              tgfunc(ji,jj,jk,jl) = exp(-1.*(tsnbio(ji,jj,jk,1)-rn_mutpft(jl))**2./rn_mudpft(jl)**2.)
            END DO
          END DO
        END DO
      END DO
!
! Increase the growth rate of crustaceans when ice is between 0.1-0.3 to represent enhanced recruitment
! as in Wiedenmann et al., Limnology and Oceanography 2009
! ---------------------------------------------------------------------------------------------
!
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
!
! zfacti is 1 when ice is greater then 0.1 and less than 0.3 and 0 otherwise 
          zice1  = fr_i(ji,jj)-0.1
          zice2  = fr_i(ji,jj)-0.3
          zfacti = max(0.,zice1)/(zice1+rtrn) * min(0.,zice2)/(zice2+rtrn)
!
! zfactl is 1 during Aug-Nov in the SH and Feb-May in the NH and 0 otherwise
          zmonth = float(nmonth)
          zfactl =  min(gphit(ji,jj),0.)/(gphit(ji,jj)+rtrn) *          & 
     &                 (max(0.,(zmonth-7.))/(zmonth-7.+rtrn)) *         &
     &                 (min(0.,(zmonth-12.))/(zmonth-12.+rtrn)) +       &
     &              max(gphit(ji,jj),0.)/(gphit(ji,jj)+rtrn) *          &
     &                 (max(0.,(zmonth-1.))/(zmonth-1.+rtrn)) *         &
     &                 (min(0.,(zmonth-6.))/(zmonth-6.+rtrn))                    
!
! zfactd is 1 when the depth does not exceed 600 m and 0 otherwise
          zfactd = min(0.,mdept(ji,jj)-600.)/(mdept(ji,jj)-600.+rtrn)
!
! multiply growth by a factor rn_icecru (from Figure 4 in Wiedeman) when the above conditions are true
!
          DO jk = 1, jpkm1
            tgfunc(ji,jj,jk,jpcru) = max(1.,zfacti*zfactl*zfactd*rn_icecru)*    &
     &                                tgfunc(ji,jj,jk,jpcru)
!
          END DO
        END DO
      END DO
!
      RETURN
#endif
      END
