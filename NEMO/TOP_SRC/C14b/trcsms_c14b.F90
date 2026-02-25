MODULE trcsms_c14b
   !!======================================================================
   !!                      ***  MODULE trcsms_c14b  ***
   !! TOP : Bomb C14 main module
   !!======================================================================
   !! History     -   ! 1994-05 ( J. Orr ) original code
   !!            1.0  ! 2006-02 ( J.M. Molines )  Free form + modularity
   !!            2.0  ! 2008-12 ( C. Ethe ) reorganisation
   !!            4.0  ! 2011-02 ( A.R. Porter, STFC Daresbury ) Dynamic memory
   !!----------------------------------------------------------------------
#if defined key_c14b
   !!----------------------------------------------------------------------
   !!   'key_c14b'                                         Bomb C14 tracer
   !!----------------------------------------------------------------------
   !!   trc_sms_c14b :  compute and add C14 suface forcing to C14 trends
   !!----------------------------------------------------------------------
   USE oce_trc       ! Ocean variables
   USE par_trc       ! TOP parameters
   USE trc           ! TOP variables
   USE trd_oce
   USE trdtrc
   USE iom           ! I/O library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_c14b       ! called in trcsms.F90
   PUBLIC   trc_sms_c14b_alloc ! called in trcini_c14b.F90

   INTEGER , PUBLIC, PARAMETER ::   jpmaxrec  = 5412 ! TJ - has something to do with mbomb C14
   INTEGER , PUBLIC, PARAMETER ::   jpmaxrec2 = 5000 ! 2 * jpmaxrec ! TJ - sets no. of lines in co2fil should be less than this
   INTEGER , PUBLIC, PARAMETER ::   jpzon     = 3    ! number of zones

   INTEGER , PUBLIC    ::   nyear_start      ! initial year for C14 files
   INTEGER , PUBLIC    ::   nyear_end    , nco2ti
   REAL(wp), PUBLIC    ::   nobombatmco2, ryear, rn_bfrd14
   REAL(wp), PUBLIC,                    DIMENSION(jpmaxrec) :: bombyear,nobombyear
   REAL(wp), PUBLIC,                    DIMENSION(jpmaxrec,jpzon) ::   nobomb,bomb       !: C14 atm data (3 zones)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)          ::   fareaz     !: Spatial Interpolation Factors
   REAL(wp), PUBLIC,                    DIMENSION(jpmaxrec2)      ::   co2year,spco2      !: Atmospheric CO2
  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)            ::   qtr_c14,qtr_b14,qtr_d14    !: flux at surface
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)            ::   qint_c14   !: cumulative flux

   REAL(wp) ::   xlambda, xdecay, xaccum       ! C14 decay coef.  
   REAL(wp) ::   xconv1 = 1._wp                ! conversion from to 
   REAL(wp) ::   xconv2 = 0.01_wp / 3600._wp   ! conversion from cm/h to m/s: 
   REAL(wp) ::   xconv3 = 1.e+3_wp             ! conversion from mol/l/atm to mol/m3/atm

   !! * Substitutions
#  include "top_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_c14b.F90 5215 2015-04-15 16:11:56Z nicolasmartin $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_c14b( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_sms_c14b  ***
      !!
      !! ** Purpose :   Compute the surface boundary contition on C14bomb
      !!      passive tracer associated with air-mer fluxes and add it to 
      !!      the general trend of tracers equations.
      !!
      !! ** Original comments from J. Orr :
      !!
      !!      Calculates the input of Bomb C-14 to the surface layer of OPA
      !!
      !!      James Orr, LMCE, 28 October 1992
      !!
      !!      Initial approach mimics that of Toggweiler, Dixon, & Bryan (1989)
      !!      (hereafter referred to as TDB) with constant gas exchange,
      !!      although in this case, a perturbation approach is used for
      !!      bomb-C14 so that both the ocean and atmosphere begin at zero.
      !!      This saves tremendous amounts of computer time since no
      !!      equilibrum run is first required (i.e., for natural C-14).
      !!      Note: Many sensitivity tests can be run with this approach and
      !!            one never has to make a run for natural C-14; otherwise,
      !!            a run for natural C-14 must be run each time that one
      !!            changes a model parameter!
      !!
      !!
      !!      19 August 1993: Modified to consider Atmospheric C-14 fom IPCC.
      !!      That is, the IPCC has provided a C-14 atmospheric record (courtesy
      !!      of Martin Heimann) for model calibration.  This model spans from
      !!      preindustrial times to present, in a format different than that
      !!      given by TDB.  It must be converted to the ORR C-14 units used
      !!      here, although in this case, the perturbation includes not only
      !!      bomb C-14 but changes due to the Suess effect.
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !
      INTEGER :: ji, jj, jk, jz     ! dummy loop indices 
      INTEGER :: iyear_beg1, iyear_end1 
      INTEGER :: im1, in1 
      INTEGER :: imonth2, im2, in2 
      REAL(wp), DIMENSION(jpzon) :: zonbc14,zonbb14       !: time interp atm C14 
      REAL(wp)                   :: zpco2at       !: time interp atm C02 
      REAL(wp) :: zt, ztp, zsk      ! dummy variables
      REAL(wp) :: zsol              ! solubility
      REAL(wp) :: zsch              ! schmidt number
      REAL(wp) :: zv2               ! wind speed ( square)
      REAL(wp) :: zpv               ! piston velocity 
      REAL(wp) :: zdemi, ztra
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zatmbc14,zatmbb14
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zdecay,ydecay,wdecay
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_c14b')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj,      zatmbc14 )
      CALL wrk_alloc( jpi, jpj,      zatmbb14 )
      CALL wrk_alloc( jpi, jpj, jpk, zdecay   )
      CALL wrk_alloc( jpi, jpj, jpk, ydecay   )
      CALL wrk_alloc( jpi, jpj, jpk, wdecay   )
      IF( kt == nittrc000 )  THEN         ! Computation of decay coeffcient
         zdemi   = 5730._wp
         xlambda = LOG(2.) / zdemi / ( nyear_len(1) * rday )
         xdecay  = EXP( - xlambda * rdt )
         xaccum  = 1._wp -  xdecay
         !
         IF( ln_rsttr ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) ' Read specific variables from C14b model '
            CALL FLUSH(numout)
            CALL iom_get( numrtr, jpdom_autoglo, 'qint_c14', qint_c14 )
            IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~'
            CALL FLUSH(numout)
         ENDIF
      ENDIF
      ! Temporal interpolation
      ! ----------------------
      iyear_beg1 = nyear - nyear_start +1
      ! For Atmospheric C14 concentrations:
      IF ( nmonth <= 6 ) THEN
         iyear_beg1 = iyear_beg1 - 1
         im1       =  6 - nmonth + 1
         im2       =  6 + nmonth - 1
      ELSE
         im1       = 12 - nmonth + 7
         im2       =      nmonth - 7
      ENDIF
      iyear_end1 = iyear_beg1 + 1
      IF (iyear_beg1 .LE. 0 .OR. iyear_end1 .GT. jpmaxrec ) THEN
        WRITE(numout,*) 'trcsms_c14b',iyear_beg1,iyear_end1,jpmaxrec,nyear,nyear_start
        CALL FLUSH(numout)
        STOP "not enough years in atmc14.dat"
      ENDIF
      !  ----------------------------------------------------------------
      !  As explained by the TDB 89 papers, C-14/C-12 is the same
      !  as C-14 concentration in this special case (no fractionation
      !  effects in this model, which thus allow direct comparison
      !  to del-C-14, after unit change from above).
      ! -----------------------------------------------------------------------
      !  Calc C-14 in atmosphere based on interp of IPCC (M. Heimann) data
      !        - Compare input to calculated C-14 for each time in record
      !-------------------------------------------------------------------------
      !  Temporal and spatial interpolation at time k
      ! --------------------------------------------------
      ! Compute atmospheric C-14 for each zone (90-20S, 20S-20N, 20-90N)
      DO jz = 1, jpzon
         zonbc14(jz) = (  nobomb(iyear_beg1,jz) * FLOAT( im1 )  &
              &         + nobomb(iyear_end1,jz) * FLOAT( im2 ) ) / 12.
         ! C-14 exactly starts at zero :
         ! JMM +Zouhair : Slightly negative values are set to 0 (when perturbation approaches)
!         zonbc14(jz) = MAX( zonbc14(jz), 0. )
         zonbb14(jz) = (  bomb(iyear_beg1,jz) * FLOAT( im1 )  &
              &         + bomb(iyear_end1,jz) * FLOAT( im2 ) ) / 12.
         ! C-14 exactly starts at zero :
         ! JMM +Zouhair : Slightly negative values are set to 0 (when perturbation approaches)
!         zonbb14(jz) = MAX( zonbb14(jz), 0. )
      END DO
      !  For each (i,j)-box, with information from the fractional area
      !  (zonmean), computes area-weighted mean to give the atmospheric C-14
      !  ----------------------------------------------------------------
      zatmbc14(:,:) = zonbc14(1) * fareaz(:,:,1)   &
         &          + zonbc14(2) * fareaz(:,:,2)   &
         &          + zonbc14(3) * fareaz(:,:,3)
      zatmbb14(:,:) = zonbb14(1) * fareaz(:,:,1)   &
         &          + zonbb14(2) * fareaz(:,:,2)   &
         &          + zonbb14(3) * fareaz(:,:,3)
      ! time interpolation of CO2 concentrations
      ryear = float(nyear) + (float(nday_year)-0.5)/365.
      DO ji = nco2ti, jpmaxrec2-1
        IF(ryear.ge.co2year(ji) .AND. ryear.le.co2year(ji+1)) THEN
          zpco2at = spco2(ji)+(spco2(ji+1) - spco2(ji))*                &
     &      (ryear - co2year(ji)) /(co2year(ji+1) - co2year(ji))
          nco2ti = ji
          EXIT
        ENDIF
      END DO
!      zpco2at = (  spco2(iyear_beg2) * FLOAT( in1 )              &
!           &     + spco2(iyear_end2) * FLOAT( in2 ) ) / 6.

      IF(lwp .AND. mod(kt,15).EQ.1) THEN
          WRITE(numout, *) 'time : ', kt, ' year :',ryear,' pCO2 : ',zpco2at ,nco2ti
          WRITE(numout, *) ' C14 year begin/end :',iyear_beg1,'/',iyear_end1,nobombyear(iyear_beg1),bombyear(iyear_beg1),   &
          &                ' C14B concen (Z1/Z2/Z3) : ',zonbc14(1),'/',zonbc14(2),'/',zonbc14(3)
      ENDIF

      !  Determine seasonally variable gas exchange coefficient 
      !----------------------------------------------------------
      !  Computes the Schmidt number of CO2 in seawater using the formulation
      !  presented by Wanninkhof (1992, J. Geophys. Res., 97,7373-7382).  
      ! -------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi  
            !  Computation of solubility  
            IF (tmask(ji,jj,1) >  0.) THEN
               ztp  = ( tsn(ji,jj,1,jp_tem) + 273.16 ) * 0.01
               zsk  = 0.023517 + ztp * ( -0.023656 + 0.0047036 * ztp )
               zsol = EXP( -60.2409 + 93.4517 / ztp  + 23.3585 * LOG( ztp ) + zsk * tsn(ji,jj,1,jp_sal) )
               ! convert solubilities [mol/(l * atm)] -> [mol/(m^3 * ppm)]
               zsol = zsol * 1.e-03
            ELSE
               zsol = 0._wp
            ENDIF

            ! speed transfert : Formulation of Wanninkhof (1992, JGR, 97,7373-7382)
            ! JMM/Zouhair : coef of 0.25 rather than 0.3332 for CORe wind speed

            ! Computes the Schmidt number of CO2 in seawater
            zt   = tsn(ji,jj,1,jp_tem)
            zsch = 2073.1 + zt * ( -125.62 + zt * (3.6276 - 0.043219 * zt ) )

            ! Wanninkhof Piston velocity and convert from units [cm/hr] -> [m/s]
            zv2  = wndm(ji,jj) * wndm(ji,jj)
            zsch = zsch / 660.
            zpv  = ( 0.31  * zv2 / SQRT(zsch) ) * xconv2 * tmask(ji,jj,1)

            ! Flux of Bomb C-14 from air-sea : speed*(conc. at equil-conc. at surf)
            ! in C-14 (orr-model-units) / m**2 * s
            qtr_c14(ji,jj) = -zpv * zsol * nobombatmco2              &
     &                       * ( trb(ji,jj,1,jpc14) - zatmbc14(ji,jj) )  &
#if defined key_degrad
     &                       * facvol(ji,jj,1)                           &
#endif
     &        * tmask(ji,jj,1) * ( 1. - fr_i(ji,jj) ) / max(trn(ji,jj,1,jpdic)*1e3,rtrn)
            ! Add the surface flux to the trend
            tra(ji,jj,1,jpc14) = tra(ji,jj,1,jpc14) + qtr_c14(ji,jj) / fse3t(ji,jj,1) 
            
            ! cumulation of surface flux at each time step
            qint_c14(ji,jj) = qint_c14(ji,jj) + qtr_c14(ji,jj) * rdt
            qtr_b14(ji,jj) = -zpv * zsol * zpco2at                       &
     &                       * ( trb(ji,jj,1,jpb14) - zatmbb14(ji,jj) )  &
     &        * tmask(ji,jj,1) * ( 1. - fr_i(ji,jj) ) / max(trn(ji,jj,1,jpdic)*1e3,rtrn)
            ! Add the surface flux to the trend
            tra(ji,jj,1,jpb14) = tra(ji,jj,1,jpb14) + qtr_b14(ji,jj) / fse3t(ji,jj,1)
            qtr_d14(ji,jj) = -zpv * zsol * zpco2at                       &
     &                       * ( trb(ji,jj,1,jpd14) - zatmbb14(ji,jj) )  &
     &        * tmask(ji,jj,1) * ( 1. - fr_i(ji,jj) ) / max(trn(ji,jj,1,jpdic)*1e3,rtrn)
            ! Add the surface flux to the trend
            tra(ji,jj,1,jpd14) = tra(ji,jj,1,jpd14) + qtr_d14(ji,jj) / fse3t(ji,jj,1)
            !
         END DO
      END DO

      ! Computation of decay effects on tracer concentration
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
#if defined key_degrad
               zdecay(ji,jj,jk) = trn(ji,jj,jk,jpc14) * ( 1. - EXP( -xlambda * rdt * facvol(ji,jj,jk) ) )
#else
               zdecay(ji,jj,jk) = trn(ji,jj,jk,jpc14) * xaccum
#endif
               tra(ji,jj,jk,jpc14) = tra(ji,jj,jk,jpc14) - zdecay(ji,jj,jk) / rdt
               ydecay(ji,jj,jk) = trn(ji,jj,jk,jpb14) * xaccum
               tra(ji,jj,jk,jpb14) = tra(ji,jj,jk,jpb14) - ydecay(ji,jj,jk) / rdt
               wdecay(ji,jj,jk) = trn(ji,jj,jk,jpd14) * xaccum
               tra(ji,jj,jk,jpd14) = tra(ji,jj,jk,jpd14) - wdecay(ji,jj,jk) / rdt
               !
            END DO
         END DO
      END DO

      !
      IF( lrst_trc ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_sms_c14b : cumulated input function fields written in ocean restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         CALL iom_rstput( kt, nitrst, numrtw, 'qint_c14', qint_c14 )
      ENDIF
      !    
      IF( lk_iomput ) THEN
        CALL iom_put( "qtrC14b"  , qtr_c14  )
        CALL iom_put( "qintC14b" , qint_c14 )
        CALL iom_put( "fdecay"   , zdecay   )
        CALL iom_put( "qtrB14b"  , qtr_b14  )
        CALL iom_put( "bdecay"   , ydecay   )
        CALL iom_put( "qtrD14b"  , qtr_d14  )
        CALL iom_put( "ddecay"   , wdecay   )
      ELSE
         IF( ln_diatrc ) THEN
            trc2d(:,:  ,jp_c14b0_2d     ) = qtr_c14 (:,:)
            trc2d(:,:  ,jp_c14b0_2d + 1 ) = qint_c14(:,:)
            trc3d(:,:,:,jp_c14b0_3d     ) = zdecay  (:,:,:)
         ENDIF
      ENDIF

      IF( l_trdtrc )  CALL trd_trc( tra(:,:,:,jpc14), jpc14, jptra_sms, kt )   ! save trends
      CALL wrk_dealloc( jpi, jpj,      zatmbc14 )
      CALL wrk_dealloc( jpi, jpj,      zatmbb14 )
      CALL wrk_dealloc( jpi, jpj, jpk, zdecay   )
      CALL wrk_dealloc( jpi, jpj, jpk, ydecay   )
      CALL wrk_dealloc( jpi, jpj, jpk, wdecay   )
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_c14b')
      !
   END SUBROUTINE trc_sms_c14b


   INTEGER FUNCTION trc_sms_c14b_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_sms_c14b_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( fareaz  (jpi,jpj ,jpzon) ,     &
         &      qtr_b14 (jpi,jpj)        ,     &
         &      qtr_c14 (jpi,jpj)        ,     &
         &      qtr_d14 (jpi,jpj)        ,     &
         &      qint_c14(jpi,jpj)        , STAT=trc_sms_c14b_alloc )
         !
      IF( trc_sms_c14b_alloc /= 0 )   CALL ctl_warn('trc_sms_c14b_alloc: failed to allocate arrays')
      !
   END FUNCTION trc_sms_c14b_alloc

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Dummy module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_c14b( kt )       ! Empty routine
      WRITE(*,*) 'trc_freons: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_c14b
#endif

  !!======================================================================
END MODULE trcsms_c14b
