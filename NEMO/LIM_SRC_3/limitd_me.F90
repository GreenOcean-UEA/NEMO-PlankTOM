MODULE limitd_me
   !!======================================================================
   !!                       ***  MODULE limitd_me ***
   !! LIM-3 : Mechanical impact on ice thickness distribution      
   !!======================================================================
   !! History :  LIM  ! 2006-02  (M. Vancoppenolle) Original code 
   !!            3.2  ! 2009-07  (M. Vancoppenolle, Y. Aksenov, G. Madec) bug correction in smsw & sfx_dyn
   !!            4.0  ! 2011-02  (G. Madec) dynamical allocation
   !!----------------------------------------------------------------------
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3'                                      LIM-3 sea-ice model
   !!----------------------------------------------------------------------
   USE par_oce          ! ocean parameters
   USE dom_oce          ! ocean domain
   USE phycst           ! physical constants (ocean directory) 
   USE sbc_oce          ! surface boundary condition: ocean fields
   USE thd_ice          ! LIM thermodynamics
   USE ice              ! LIM variables
   USE dom_ice          ! LIM domain
   USE limvar           ! LIM
   USE lbclnk           ! lateral boundary condition - MPP exchanges
   USE lib_mpp          ! MPP library
   USE wrk_nemo         ! work arrays
   USE prtctl           ! Print control

   USE in_out_manager   ! I/O manager
   USE iom              ! I/O manager
   USE lib_fortran      ! glob_sum
   USE limdiahsb
   USE timing           ! Timing
   USE limcons          ! conservation tests

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_itd_me               ! called by ice_stp
   PUBLIC   lim_itd_me_icestrength
   PUBLIC   lim_itd_me_init
   PUBLIC   lim_itd_me_alloc        ! called by sbc_lim_init 

   !-----------------------------------------------------------------------
   ! Variables shared among ridging subroutines
   !-----------------------------------------------------------------------
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   asum     ! sum of total ice and open water area
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   aksum    ! ratio of area removed to area ridged
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   athorn   ! participation function; fraction of ridging/
   !                                                     ! closing associated w/ category n
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   hrmin    ! minimum ridge thickness
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   hrmax    ! maximum ridge thickness
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   hraft    ! thickness of rafted ice
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   krdg     ! mean ridge thickness/thickness of ridging ice 
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   aridge   ! participating ice ridging
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   araft    ! participating ice rafting

   REAL(wp), PARAMETER ::   krdgmin = 1.1_wp    ! min ridge thickness multiplier
   REAL(wp), PARAMETER ::   kraft   = 2.0_wp    ! rafting multipliyer
   REAL(wp), PARAMETER ::   kamax   = 1.0_wp    ! max of ice area authorized (clem: scheme is not stable if kamax <= 0.99)

   REAL(wp) ::   Cp                             ! 
   !
   !-----------------------------------------------------------------------
   ! Ridging diagnostic arrays for history files
   !-----------------------------------------------------------------------
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   dardg1dt   ! rate of fractional area loss by ridging ice (1/s)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   dardg2dt   ! rate of fractional area gain by new ridges (1/s)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   dvirdgdt   ! rate of ice volume ridged (m/s)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   opening    ! rate of opening due to divergence/shear (1/s)
   !
   !!----------------------------------------------------------------------
   !! NEMO/LIM3 3.3 , UCL - NEMO Consortium (2010)
   !! $Id: limitd_me.F90 5202 2015-04-07 15:40:16Z clem $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION lim_itd_me_alloc()
      !!---------------------------------------------------------------------!
      !!                ***  ROUTINE lim_itd_me_alloc ***
      !!---------------------------------------------------------------------!
      ALLOCATE(                                                                     &
         !* Variables shared among ridging subroutines
         &      asum (jpi,jpj)     , athorn(jpi,jpj,0:jpl)                    ,     &
         &      aksum(jpi,jpj)                                                ,     &
         !
         &      hrmin(jpi,jpj,jpl) , hraft(jpi,jpj,jpl) , aridge(jpi,jpj,jpl) ,     &
         &      hrmax(jpi,jpj,jpl) , krdg (jpi,jpj,jpl) , araft (jpi,jpj,jpl) ,     &
         !
         !* Ridging diagnostic arrays for history files
         &      dardg1dt(jpi,jpj)  , dardg2dt(jpi,jpj)                        ,     & 
         &      dvirdgdt(jpi,jpj)  , opening(jpi,jpj)                         , STAT=lim_itd_me_alloc )
         !
      IF( lim_itd_me_alloc /= 0 )   CALL ctl_warn( 'lim_itd_me_alloc: failed to allocate arrays' )
      !
   END FUNCTION lim_itd_me_alloc


   SUBROUTINE lim_itd_me
      !!---------------------------------------------------------------------!
      !!                ***  ROUTINE lim_itd_me ***
      !!
      !! ** Purpose :   computes the mechanical redistribution of ice thickness
      !!
      !! ** Method  :   Steps :
      !!       1) Thickness categories boundaries, ice / o.w. concentrations
      !!          Ridge preparation
      !!       2) Dynamical inputs (closing rate, divu_adv, opning)
      !!       3) Ridging iteration
      !!       4) Ridging diagnostics
      !!       5) Heat, salt and freshwater fluxes
      !!       6) Compute increments of tate variables and come back to old values
      !!
      !! References :   Flato, G. M., and W. D. Hibler III, 1995, JGR, 100, 18,611-18,626.
      !!                Hibler, W. D. III, 1980, MWR, 108, 1943-1973, 1980.
      !!                Rothrock, D. A., 1975: JGR, 80, 4514-4519.
      !!                Thorndike et al., 1975, JGR, 80, 4501-4513. 
      !!                Bitz et al., JGR, 2001
      !!                Amundrud and Melling, JGR 2005
      !!                Babko et al., JGR 2002 
      !!
      !!     This routine is based on CICE code and authors William H. Lipscomb,
      !!  and Elizabeth C. Hunke, LANL are gratefully acknowledged
      !!--------------------------------------------------------------------!
      INTEGER  ::   ji, jj, jk, jl        ! dummy loop index
      INTEGER  ::   niter                 ! local integer 
      INTEGER  ::   iterate_ridging       ! if true, repeat the ridging
      REAL(wp) ::   za, zfac              ! local scalar
      CHARACTER (len = 15) ::   fieldid
      REAL(wp), POINTER, DIMENSION(:,:)   ::   closing_net     ! net rate at which area is removed    (1/s)
                                                               ! (ridging ice area - area of new ridges) / dt
      REAL(wp), POINTER, DIMENSION(:,:)   ::   divu_adv        ! divu as implied by transport scheme  (1/s)
      REAL(wp), POINTER, DIMENSION(:,:)   ::   opning          ! rate of opening due to divergence/shear
      REAL(wp), POINTER, DIMENSION(:,:)   ::   closing_gross   ! rate at which area removed, not counting area of new ridges
      REAL(wp), POINTER, DIMENSION(:,:)   ::   msnow_mlt       ! mass of snow added to ocean (kg m-2)
      REAL(wp), POINTER, DIMENSION(:,:)   ::   esnow_mlt       ! energy needed to melt snow in ocean (J m-2)
      REAL(wp), POINTER, DIMENSION(:,:)   ::   vt_i_init, vt_i_final  !  ice volume summed over categories
      !
      INTEGER, PARAMETER ::   nitermax = 20    
      !
      REAL(wp) :: zvi_b, zsmv_b, zei_b, zfs_b, zfw_b, zft_b 
      !!-----------------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('limitd_me')

      CALL wrk_alloc( jpi,jpj, closing_net, divu_adv, opning, closing_gross, msnow_mlt, esnow_mlt, vt_i_init, vt_i_final )

      IF(ln_ctl) THEN
         CALL prt_ctl(tab2d_1=ato_i , clinfo1=' lim_itd_me: ato_i  : ', tab2d_2=at_i   , clinfo2=' at_i    : ')
         CALL prt_ctl(tab2d_1=divu_i, clinfo1=' lim_itd_me: divu_i : ', tab2d_2=delta_i, clinfo2=' delta_i : ')
      ENDIF

      IF( ln_limdyn ) THEN          !   Start ridging and rafting   !

      ! conservation test
      IF( ln_limdiahsb ) CALL lim_cons_hsm(0, 'limitd_me', zvi_b, zsmv_b, zei_b, zfw_b, zfs_b, zft_b)

      CALL lim_var_zapsmall
      CALL lim_var_glo2eqv            ! equivalent variables, requested for rafting

      !-----------------------------------------------------------------------------!
      ! 1) Thickness categories boundaries, ice / o.w. concentrations, init_ons
      !-----------------------------------------------------------------------------!
      Cp = 0.5 * grav * (rau0-rhoic) * rhoic * r1_rau0             ! proport const for PE
      !
      CALL lim_itd_me_ridgeprep                                    ! prepare ridging
      !
      IF( con_i)   CALL lim_column_sum( jpl, v_i, vt_i_init )      ! conservation check

      DO jj = 1, jpj                                     ! Initialize arrays.
         DO ji = 1, jpi
            msnow_mlt(ji,jj) = 0._wp
            esnow_mlt(ji,jj) = 0._wp
            dardg1dt (ji,jj) = 0._wp
            dardg2dt (ji,jj) = 0._wp
            dvirdgdt (ji,jj) = 0._wp
            opening  (ji,jj) = 0._wp

            !-----------------------------------------------------------------------------!
            ! 2) Dynamical inputs (closing rate, divu_adv, opning)
            !-----------------------------------------------------------------------------!
            !
            ! 2.1 closing_net
            !-----------------
            ! Compute the net rate of closing due to convergence 
            ! and shear, based on Flato and Hibler (1995).
            ! 
            ! The energy dissipation rate is equal to the net closing rate
            ! times the ice strength.
            !
            ! NOTE: The NET closing rate is equal to the rate that open water 
            !  area is removed, plus the rate at which ice area is removed by 
            !  ridging, minus the rate at which area is added in new ridges.
            !  The GROSS closing rate is equal to the first two terms (open
            !  water closing and thin ice ridging) without the third term
            !  (thick, newly ridged ice).

            closing_net(ji,jj) = rn_cs * 0.5 * ( delta_i(ji,jj) - ABS( divu_i(ji,jj) ) ) - MIN( divu_i(ji,jj), 0._wp )

            ! 2.2 divu_adv
            !--------------
            ! Compute divu_adv, the divergence rate given by the transport/
            ! advection scheme, which may not be equal to divu as computed 
            ! from the velocity field.
            !
            ! If divu_adv < 0, make sure the closing rate is large enough
            ! to give asum = 1.0 after ridging.

            divu_adv(ji,jj) = ( kamax - asum(ji,jj) ) * r1_rdtice  ! asum found in ridgeprep

            IF( divu_adv(ji,jj) < 0._wp )   closing_net(ji,jj) = MAX( closing_net(ji,jj), -divu_adv(ji,jj) )

            ! 2.3 opning
            !------------
            ! Compute the (non-negative) opening rate that will give asum = 1.0 after ridging.
            opning(ji,jj) = closing_net(ji,jj) + divu_adv(ji,jj)
         END DO
      END DO

      !-----------------------------------------------------------------------------!
      ! 3) Ridging iteration
      !-----------------------------------------------------------------------------!
      niter           = 1                 ! iteration counter
      iterate_ridging = 1

      DO WHILE ( iterate_ridging > 0 .AND. niter < nitermax )

         DO jj = 1, jpj
            DO ji = 1, jpi

               ! 3.2 closing_gross
               !-----------------------------------------------------------------------------!
               ! Based on the ITD of ridging and ridged ice, convert the net
               !  closing rate to a gross closing rate.  
               ! NOTE: 0 < aksum <= 1
               closing_gross(ji,jj) = closing_net(ji,jj) / aksum(ji,jj)

               ! correction to closing rate and opening if closing rate is excessive
               !---------------------------------------------------------------------
               ! Reduce the closing rate if more than 100% of the open water 
               ! would be removed.  Reduce the opening rate proportionately.
               za   = athorn(ji,jj,0) * closing_gross(ji,jj) * rdt_ice
               IF( za > epsi20 ) THEN
                  zfac = MIN( 1._wp, ato_i(ji,jj) / za )
                  closing_gross(ji,jj) = closing_gross(ji,jj) * zfac
                  opning       (ji,jj) = opning       (ji,jj) * zfac
               ENDIF

            END DO
         END DO

         ! correction to closing rate / opening if excessive ice removal
         !---------------------------------------------------------------
         ! Reduce the closing rate if more than 100% of any ice category 
         ! would be removed.  Reduce the opening rate proportionately.
         DO jl = 1, jpl
            DO jj = 1, jpj
               DO ji = 1, jpi
                  za = athorn(ji,jj,jl) * closing_gross(ji,jj) * rdt_ice
                  IF( za  >  epsi20 ) THEN
                     zfac = MIN( 1._wp, a_i(ji,jj,jl) / za )
                     closing_gross(ji,jj) = closing_gross(ji,jj) * zfac
                     opning       (ji,jj) = opning       (ji,jj) * zfac
                  ENDIF
               END DO
            END DO
         END DO

         ! 3.3 Redistribute area, volume, and energy.
         !-----------------------------------------------------------------------------!

         CALL lim_itd_me_ridgeshift( opning, closing_gross, msnow_mlt, esnow_mlt )

         ! 3.4 Compute total area of ice plus open water after ridging.
         !-----------------------------------------------------------------------------!
         ! This is in general not equal to one because of divergence during transport
         asum(:,:) = ato_i(:,:)
         DO jl = 1, jpl
            asum(:,:) = asum(:,:) + a_i(:,:,jl)
         END DO

         ! 3.5 Do we keep on iterating ???
         !-----------------------------------------------------------------------------!
         ! Check whether asum = 1.  If not (because the closing and opening
         ! rates were reduced above), ridge again with new rates.

         iterate_ridging = 0

         DO jj = 1, jpj
            DO ji = 1, jpi
               IF (ABS(asum(ji,jj) - kamax ) < epsi10) THEN
                  closing_net(ji,jj) = 0._wp
                  opning     (ji,jj) = 0._wp
               ELSE
                  iterate_ridging    = 1
                  divu_adv   (ji,jj) = ( kamax - asum(ji,jj) ) * r1_rdtice
                  closing_net(ji,jj) = MAX( 0._wp, -divu_adv(ji,jj) )
                  opning     (ji,jj) = MAX( 0._wp,  divu_adv(ji,jj) )
               ENDIF
            END DO
         END DO

         IF( lk_mpp )   CALL mpp_max( iterate_ridging )

         ! Repeat if necessary.
         ! NOTE: If strength smoothing is turned on, the ridging must be
         ! iterated globally because of the boundary update in the 
         ! smoothing.

         niter = niter + 1

         IF( iterate_ridging == 1 ) THEN
            IF( niter  >  nitermax ) THEN
               WRITE(numout,*) ' ALERTE : non-converging ridging scheme '
               WRITE(numout,*) ' niter, iterate_ridging ', niter, iterate_ridging
            ENDIF
            CALL lim_itd_me_ridgeprep
         ENDIF

      END DO !! on the do while over iter

      !-----------------------------------------------------------------------------!
      ! 4) Ridging diagnostics
      !-----------------------------------------------------------------------------!
      ! Convert ridging rate diagnostics to correct units.
      ! Update fresh water and heat fluxes due to snow melt.
      DO jj = 1, jpj
         DO ji = 1, jpi

            dardg1dt(ji,jj) = dardg1dt(ji,jj) * r1_rdtice
            dardg2dt(ji,jj) = dardg2dt(ji,jj) * r1_rdtice
            dvirdgdt(ji,jj) = dvirdgdt(ji,jj) * r1_rdtice
            opening (ji,jj) = opening (ji,jj) * r1_rdtice

            !-----------------------------------------------------------------------------!
            ! 5) Heat, salt and freshwater fluxes
            !-----------------------------------------------------------------------------!
            wfx_snw(ji,jj) = wfx_snw(ji,jj) + msnow_mlt(ji,jj) * r1_rdtice     ! fresh water source for ocean
            hfx_dyn(ji,jj) = hfx_dyn(ji,jj) + esnow_mlt(ji,jj) * r1_rdtice     ! heat sink for ocean (<0, W.m-2)

         END DO
      END DO

      ! Check if there is a ridging error
      IF( lwp ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( ABS( asum(ji,jj) - kamax)  >  epsi10 ) THEN   ! there is a bug
                  WRITE(numout,*) ' '
                  WRITE(numout,*) ' ALERTE : Ridging error: total area = ', asum(ji,jj)
                  WRITE(numout,*) ' limitd_me '
                  WRITE(numout,*) ' POINT : ', ji, jj
                  WRITE(numout,*) ' jpl, a_i, athorn '
                  WRITE(numout,*) 0, ato_i(ji,jj), athorn(ji,jj,0)
                  DO jl = 1, jpl
                     WRITE(numout,*) jl, a_i(ji,jj,jl), athorn(ji,jj,jl)
                  END DO
               ENDIF
            END DO
         END DO
      END IF

      ! Conservation check
      IF ( con_i ) THEN
         CALL lim_column_sum (jpl,   v_i, vt_i_final)
         fieldid = ' v_i : limitd_me '
         CALL lim_cons_check (vt_i_init, vt_i_final, 1.0e-6, fieldid) 
      ENDIF

      CALL lim_var_agg( 1 ) 

      !-----------------------------------------------------------------------------!
      ! control prints
      !-----------------------------------------------------------------------------!
      IF(ln_ctl) THEN 
         CALL lim_var_glo2eqv

         CALL prt_ctl_info(' ')
         CALL prt_ctl_info(' - Cell values : ')
         CALL prt_ctl_info('   ~~~~~~~~~~~~~ ')
         CALL prt_ctl(tab2d_1=e12t , clinfo1=' lim_itd_me  : cell area :')
         CALL prt_ctl(tab2d_1=at_i , clinfo1=' lim_itd_me  : at_i      :')
         CALL prt_ctl(tab2d_1=vt_i , clinfo1=' lim_itd_me  : vt_i      :')
         CALL prt_ctl(tab2d_1=vt_s , clinfo1=' lim_itd_me  : vt_s      :')
         DO jl = 1, jpl
            CALL prt_ctl_info(' ')
            CALL prt_ctl_info(' - Category : ', ivar1=jl)
            CALL prt_ctl_info('   ~~~~~~~~~~')
            CALL prt_ctl(tab2d_1=a_i   (:,:,jl)   , clinfo1= ' lim_itd_me  : a_i      : ')
            CALL prt_ctl(tab2d_1=ht_i  (:,:,jl)   , clinfo1= ' lim_itd_me  : ht_i     : ')
            CALL prt_ctl(tab2d_1=ht_s  (:,:,jl)   , clinfo1= ' lim_itd_me  : ht_s     : ')
            CALL prt_ctl(tab2d_1=v_i   (:,:,jl)   , clinfo1= ' lim_itd_me  : v_i      : ')
            CALL prt_ctl(tab2d_1=v_s   (:,:,jl)   , clinfo1= ' lim_itd_me  : v_s      : ')
            CALL prt_ctl(tab2d_1=e_s   (:,:,1,jl) , clinfo1= ' lim_itd_me  : e_s      : ')
            CALL prt_ctl(tab2d_1=t_su  (:,:,jl)   , clinfo1= ' lim_itd_me  : t_su     : ')
            CALL prt_ctl(tab2d_1=t_s   (:,:,1,jl) , clinfo1= ' lim_itd_me  : t_snow   : ')
            CALL prt_ctl(tab2d_1=sm_i  (:,:,jl)   , clinfo1= ' lim_itd_me  : sm_i     : ')
            CALL prt_ctl(tab2d_1=smv_i (:,:,jl)   , clinfo1= ' lim_itd_me  : smv_i    : ')
            DO jk = 1, nlay_i
               CALL prt_ctl_info(' ')
               CALL prt_ctl_info(' - Layer : ', ivar1=jk)
               CALL prt_ctl_info('   ~~~~~~~')
               CALL prt_ctl(tab2d_1=t_i(:,:,jk,jl) , clinfo1= ' lim_itd_me  : t_i      : ')
               CALL prt_ctl(tab2d_1=e_i(:,:,jk,jl) , clinfo1= ' lim_itd_me  : e_i      : ')
            END DO
         END DO
      ENDIF

      ! conservation test
      IF( ln_limdiahsb ) CALL lim_cons_hsm(1, 'limitd_me', zvi_b, zsmv_b, zei_b, zfw_b, zfs_b, zft_b)

      ENDIF  ! ln_limdyn=.true.
      !
      CALL wrk_dealloc( jpi, jpj, closing_net, divu_adv, opning, closing_gross, msnow_mlt, esnow_mlt, vt_i_init, vt_i_final )
      !
      IF( nn_timing == 1 )  CALL timing_stop('limitd_me')
   END SUBROUTINE lim_itd_me


   SUBROUTINE lim_itd_me_icestrength( kstrngth )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE lim_itd_me_icestrength ***
      !!
      !! ** Purpose :   computes ice strength used in dynamics routines of ice thickness
      !!
      !! ** Method  :   Compute the strength of the ice pack, defined as the energy (J m-2) 
      !!              dissipated per unit area removed from the ice pack under compression,
      !!              and assumed proportional to the change in potential energy caused
      !!              by ridging. Note that only Hibler's formulation is stable and that
      !!              ice strength has to be smoothed
      !!
      !! ** Inputs / Ouputs : kstrngth (what kind of ice strength we are using)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kstrngth    ! = 1 for Rothrock formulation, 0 for Hibler (1979)
      INTEGER             ::   ji,jj, jl   ! dummy loop indices
      INTEGER             ::   ksmooth     ! smoothing the resistance to deformation
      INTEGER             ::   numts_rm    ! number of time steps for the P smoothing
      REAL(wp)            ::   zhi, zp, z1_3  ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::   zworka   ! temporary array used here
      !!----------------------------------------------------------------------

      CALL wrk_alloc( jpi, jpj, zworka )

      !------------------------------------------------------------------------------!
      ! 1) Initialize
      !------------------------------------------------------------------------------!
      strength(:,:) = 0._wp

      !------------------------------------------------------------------------------!
      ! 2) Compute thickness distribution of ridging and ridged ice
      !------------------------------------------------------------------------------!
      CALL lim_itd_me_ridgeprep

      !------------------------------------------------------------------------------!
      ! 3) Rothrock(1975)'s method
      !------------------------------------------------------------------------------!
      IF( kstrngth == 1 ) THEN
         z1_3 = 1._wp / 3._wp
         DO jl = 1, jpl
            DO jj= 1, jpj
               DO ji = 1, jpi
                  !
                  IF( a_i(ji,jj,jl) > epsi10 .AND. athorn(ji,jj,jl) > 0._wp ) THEN
                     zhi = v_i(ji,jj,jl) / a_i(ji,jj,jl)
                     !----------------------------
                     ! PE loss from deforming ice
                     !----------------------------
                     strength(ji,jj) = strength(ji,jj) - athorn(ji,jj,jl) * zhi * zhi

                     !--------------------------
                     ! PE gain from rafting ice
                     !--------------------------
                     strength(ji,jj) = strength(ji,jj) + 2._wp * araft(ji,jj,jl) * zhi * zhi

                     !----------------------------
                     ! PE gain from ridging ice
                     !----------------------------
                     strength(ji,jj) = strength(ji,jj) + aridge(ji,jj,jl) / krdg(ji,jj,jl)     &
                        * z1_3 * ( hrmax(ji,jj,jl)**2 + hrmin(ji,jj,jl)**2 +  hrmax(ji,jj,jl) * hrmin(ji,jj,jl) )  
                        !!(a**3-b**3)/(a-b) = a*a+ab+b*b                      
                  ENDIF
                  !
               END DO
            END DO
         END DO
   
         strength(:,:) = rn_pe_rdg * Cp * strength(:,:) / aksum(:,:)
                         ! where Cp = (g/2)*(rhow-rhoi)*(rhoi/rhow) and rn_pe_rdg accounts for frictional dissipation
         ksmooth = 1

         !------------------------------------------------------------------------------!
         ! 4) Hibler (1979)' method
         !------------------------------------------------------------------------------!
      ELSE                      ! kstrngth ne 1:  Hibler (1979) form
         !
         strength(:,:) = rn_pstar * vt_i(:,:) * EXP( - rn_crhg * ( 1._wp - at_i(:,:) )  )
         !
         ksmooth = 1
         !
      ENDIF                     ! kstrngth

      !
      !------------------------------------------------------------------------------!
      ! 5) Impact of brine volume
      !------------------------------------------------------------------------------!
      ! CAN BE REMOVED
      !
      IF( ln_icestr_bvf ) THEN

         DO jj = 1, jpj
            DO ji = 1, jpi
               strength(ji,jj) = strength(ji,jj) * exp(-5.88*SQRT(MAX(bv_i(ji,jj),0.0)))
            END DO
         END DO

      ENDIF

      !
      !------------------------------------------------------------------------------!
      ! 6) Smoothing ice strength
      !------------------------------------------------------------------------------!
      !
      !-------------------
      ! Spatial smoothing
      !-------------------
      IF ( ksmooth == 1 ) THEN

         CALL lbc_lnk( strength, 'T', 1. )

         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               IF ( ( asum(ji,jj) - ato_i(ji,jj) ) > 0._wp) THEN 
                  zworka(ji,jj) = ( 4.0 * strength(ji,jj)              &
                     &                  + strength(ji-1,jj) * tmask(ji-1,jj,1) + strength(ji+1,jj) * tmask(ji+1,jj,1) &  
                     &                  + strength(ji,jj-1) * tmask(ji,jj-1,1) + strength(ji,jj+1) * tmask(ji,jj+1,1) &
                     &            ) / ( 4.0 + tmask(ji-1,jj,1) + tmask(ji+1,jj,1) + tmask(ji,jj-1,1) + tmask(ji,jj+1,1) )
               ELSE
                  zworka(ji,jj) = 0._wp
               ENDIF
            END DO
         END DO

         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               strength(ji,jj) = zworka(ji,jj)
            END DO
         END DO
         CALL lbc_lnk( strength, 'T', 1. )

      ENDIF

      !--------------------
      ! Temporal smoothing
      !--------------------
      IF ( numit == nit000 + nn_fsbc - 1 ) THEN
         strp1(:,:) = 0.0            
         strp2(:,:) = 0.0            
      ENDIF

      IF ( ksmooth == 2 ) THEN


         CALL lbc_lnk( strength, 'T', 1. )

         DO jj = 1, jpj - 1
            DO ji = 1, jpi - 1
               IF ( ( asum(ji,jj) - ato_i(ji,jj) ) > 0._wp) THEN 
                  numts_rm = 1 ! number of time steps for the running mean
                  IF ( strp1(ji,jj) > 0.0 ) numts_rm = numts_rm + 1
                  IF ( strp2(ji,jj) > 0.0 ) numts_rm = numts_rm + 1
                  zp = ( strength(ji,jj) + strp1(ji,jj) + strp2(ji,jj) ) / numts_rm
                  strp2(ji,jj) = strp1(ji,jj)
                  strp1(ji,jj) = strength(ji,jj)
                  strength(ji,jj) = zp

               ENDIF
            END DO
         END DO

      ENDIF ! ksmooth

      CALL lbc_lnk( strength, 'T', 1. )      ! Boundary conditions

      CALL wrk_dealloc( jpi, jpj, zworka )
      !
   END SUBROUTINE lim_itd_me_icestrength


   SUBROUTINE lim_itd_me_ridgeprep
      !!---------------------------------------------------------------------!
      !!                ***  ROUTINE lim_itd_me_ridgeprep ***
      !!
      !! ** Purpose :   preparation for ridging and strength calculations
      !!
      !! ** Method  :   Compute the thickness distribution of the ice and open water 
      !!              participating in ridging and of the resulting ridges.
      !!---------------------------------------------------------------------!
      INTEGER ::   ji,jj, jl    ! dummy loop indices
      REAL(wp) ::   Gstari, astari, zhi, hrmean, zdummy   ! local scalar
      REAL(wp), POINTER, DIMENSION(:,:)   ::   zworka    ! temporary array used here
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   Gsum      ! Gsum(n) = sum of areas in categories 0 to n
      !------------------------------------------------------------------------------!

      CALL wrk_alloc( jpi,jpj, zworka )
      CALL wrk_alloc( jpi,jpj,jpl+2, Gsum, kkstart = -1 )

      Gstari     = 1.0/rn_gstar    
      astari     = 1.0/rn_astar    
      aksum(:,:)    = 0.0
      athorn(:,:,:) = 0.0
      aridge(:,:,:) = 0.0
      araft (:,:,:) = 0.0
      hrmin(:,:,:)  = 0.0 
      hrmax(:,:,:)  = 0.0 
      hraft(:,:,:)  = 0.0 
      krdg (:,:,:)  = 1.0

      !     ! Zero out categories with very small areas
      CALL lim_var_zapsmall

      !------------------------------------------------------------------------------!
      ! 1) Participation function 
      !------------------------------------------------------------------------------!

      ! Compute total area of ice plus open water.
      ! This is in general not equal to one because of divergence during transport
      asum(:,:) = ato_i(:,:)
      DO jl = 1, jpl
         asum(:,:) = asum(:,:) + a_i(:,:,jl)
      END DO

      ! Compute cumulative thickness distribution function
      ! Compute the cumulative thickness distribution function Gsum,
      ! where Gsum(n) is the fractional area in categories 0 to n.
      ! initial value (in h = 0) equals open water area

      Gsum(:,:,-1) = 0._wp
      Gsum(:,:,0 ) = ato_i(:,:)

      ! for each value of h, you have to add ice concentration then
      DO jl = 1, jpl
         Gsum(:,:,jl) = Gsum(:,:,jl-1) + a_i(:,:,jl)
      END DO

      ! Normalize the cumulative distribution to 1
      zworka(:,:) = 1._wp / Gsum(:,:,jpl)
      DO jl = 0, jpl
         Gsum(:,:,jl) = Gsum(:,:,jl) * zworka(:,:)
      END DO

      ! 1.3 Compute participation function a(h) = b(h).g(h) (athorn)
      !--------------------------------------------------------------------------------------------------
      ! Compute the participation function athorn; this is analogous to
      ! a(h) = b(h)g(h) as defined in Thorndike et al. (1975).
      ! area lost from category n due to ridging/closing
      ! athorn(n)   = total area lost due to ridging/closing
      ! assume b(h) = (2/Gstar) * (1 - G(h)/Gstar). 
      !
      ! The expressions for athorn are found by integrating b(h)g(h) between
      ! the category boundaries.
      !-----------------------------------------------------------------

      IF( nn_partfun == 0 ) THEN       !--- Linear formulation (Thorndike et al., 1975)
         DO jl = 0, jpl    
            DO jj = 1, jpj 
               DO ji = 1, jpi
                  IF( Gsum(ji,jj,jl) < rn_gstar) THEN
                     athorn(ji,jj,jl) = Gstari * ( Gsum(ji,jj,jl) - Gsum(ji,jj,jl-1) ) * &
                        &                        ( 2.0 - (Gsum(ji,jj,jl-1) + Gsum(ji,jj,jl) ) * Gstari )
                  ELSEIF (Gsum(ji,jj,jl-1) < rn_gstar) THEN
                     athorn(ji,jj,jl) = Gstari * ( rn_gstar - Gsum(ji,jj,jl-1) ) *  &
                        &                        ( 2.0 - ( Gsum(ji,jj,jl-1) + rn_gstar ) * Gstari )
                  ELSE
                     athorn(ji,jj,jl) = 0.0
                  ENDIF
               END DO
            END DO
         END DO

      ELSE                             !--- Exponential, more stable formulation (Lipscomb et al, 2007)
         !                        
         zdummy = 1._wp / ( 1._wp - EXP(-astari) )        ! precompute exponential terms using Gsum as a work array
         DO jl = -1, jpl
            Gsum(:,:,jl) = EXP( -Gsum(:,:,jl) * astari ) * zdummy
         END DO
         DO jl = 0, jpl
             athorn(:,:,jl) = Gsum(:,:,jl-1) - Gsum(:,:,jl)
         END DO
         !
      ENDIF

      IF( ln_rafting ) THEN      ! Ridging and rafting ice participation functions
         !
         DO jl = 1, jpl
            DO jj = 1, jpj 
               DO ji = 1, jpi
                  IF ( athorn(ji,jj,jl) > 0._wp ) THEN
!!gm  TANH( -X ) = - TANH( X )  so can be computed only 1 time....
                     aridge(ji,jj,jl) = ( TANH (  rn_craft * ( ht_i(ji,jj,jl) - rn_hraft ) ) + 1.0 ) * 0.5 * athorn(ji,jj,jl)
                     araft (ji,jj,jl) = ( TANH ( -rn_craft * ( ht_i(ji,jj,jl) - rn_hraft ) ) + 1.0 ) * 0.5 * athorn(ji,jj,jl)
                     IF ( araft(ji,jj,jl) < epsi06 )   araft(ji,jj,jl)  = 0._wp
                     aridge(ji,jj,jl) = MAX( athorn(ji,jj,jl) - araft(ji,jj,jl), 0.0 )
                  ENDIF
               END DO
            END DO
         END DO

      ELSE
         !
         DO jl = 1, jpl
            aridge(:,:,jl) = athorn(:,:,jl)
         END DO
         !
      ENDIF

      IF( ln_rafting ) THEN

         IF( MAXVAL(aridge + araft - athorn(:,:,1:jpl)) > epsi10 .AND. lwp ) THEN
            DO jl = 1, jpl
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     IF ( aridge(ji,jj,jl) + araft(ji,jj,jl) - athorn(ji,jj,jl) > epsi10 ) THEN
                        WRITE(numout,*) ' ALERTE 96 : wrong participation function ... '
                        WRITE(numout,*) ' ji, jj, jl : ', ji, jj, jl
                        WRITE(numout,*) ' lat, lon   : ', gphit(ji,jj), glamt(ji,jj)
                        WRITE(numout,*) ' aridge     : ', aridge(ji,jj,1:jpl)
                        WRITE(numout,*) ' araft      : ', araft(ji,jj,1:jpl)
                        WRITE(numout,*) ' athorn     : ', athorn(ji,jj,1:jpl)
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF

      ENDIF

      !-----------------------------------------------------------------
      ! 2) Transfer function
      !-----------------------------------------------------------------
      ! Compute max and min ridged ice thickness for each ridging category.
      ! Assume ridged ice is uniformly distributed between hrmin and hrmax.
      ! 
      ! This parameterization is a modified version of Hibler (1980).
      ! The mean ridging thickness, hrmean, is proportional to hi^(0.5)
      !  and for very thick ridging ice must be >= krdgmin*hi
      !
      ! The minimum ridging thickness, hrmin, is equal to 2*hi 
      !  (i.e., rafting) and for very thick ridging ice is
      !  constrained by hrmin <= (hrmean + hi)/2.
      ! 
      ! The maximum ridging thickness, hrmax, is determined by
      !  hrmean and hrmin.
      !
      ! These modifications have the effect of reducing the ice strength
      ! (relative to the Hibler formulation) when very thick ice is
      ! ridging.
      !
      ! aksum = net area removed/ total area removed
      ! where total area removed = area of ice that ridges
      !         net area removed = total area removed - area of new ridges
      !-----------------------------------------------------------------

      ! Transfer function
      DO jl = 1, jpl !all categories have a specific transfer function
         DO jj = 1, jpj
            DO ji = 1, jpi

               IF (a_i(ji,jj,jl) > epsi10 .AND. athorn(ji,jj,jl) > 0.0 ) THEN
                  zhi = v_i(ji,jj,jl) / a_i(ji,jj,jl)
                  hrmean          = MAX(SQRT(rn_hstar*zhi), zhi*krdgmin)
                  hrmin(ji,jj,jl) = MIN(2.0*zhi, 0.5*(hrmean + zhi))
                  hrmax(ji,jj,jl) = 2.0*hrmean - hrmin(ji,jj,jl)
                  hraft(ji,jj,jl) = kraft*zhi
                  krdg(ji,jj,jl)  = hrmean / zhi
               ELSE
                  hraft(ji,jj,jl) = 0.0
                  hrmin(ji,jj,jl) = 0.0 
                  hrmax(ji,jj,jl) = 0.0 
                  krdg (ji,jj,jl) = 1.0
               ENDIF

            END DO
         END DO
      END DO

      ! Normalization factor : aksum, ensures mass conservation
      aksum(:,:) = athorn(:,:,0)
      DO jl = 1, jpl 
         aksum(:,:)    = aksum(:,:) + aridge(:,:,jl) * ( 1._wp - 1._wp / krdg(:,:,jl) )    &
            &                       + araft (:,:,jl) * ( 1._wp - 1._wp / kraft        )
      END DO
      !
      CALL wrk_dealloc( jpi,jpj, zworka )
      CALL wrk_dealloc( jpi,jpj,jpl+2, Gsum, kkstart = -1 )
      !
   END SUBROUTINE lim_itd_me_ridgeprep


   SUBROUTINE lim_itd_me_ridgeshift( opning, closing_gross, msnow_mlt, esnow_mlt )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE lim_itd_me_icestrength ***
      !!
      !! ** Purpose :   shift ridging ice among thickness categories of ice thickness
      !!
      !! ** Method  :   Remove area, volume, and energy from each ridging category
      !!              and add to thicker ice categories.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   opning         ! rate of opening due to divergence/shear
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   closing_gross  ! rate at which area removed, excluding area of new ridges
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   msnow_mlt      ! mass of snow added to ocean (kg m-2)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   esnow_mlt      ! energy needed to melt snow in ocean (J m-2)
      !
      CHARACTER (len=80) ::   fieldid   ! field identifier
      LOGICAL, PARAMETER ::   l_conservation_check = .true.  ! if true, check conservation (useful for debugging)
      !
      INTEGER ::   ji, jj, jl, jl1, jl2, jk   ! dummy loop indices
      INTEGER ::   ij                ! horizontal index, combines i and j loops
      INTEGER ::   icells            ! number of cells with aicen > puny
      REAL(wp) ::   hL, hR, farea, ztmelts    ! left and right limits of integration

      INTEGER , POINTER, DIMENSION(:) ::   indxi, indxj   ! compressed indices

      REAL(wp), POINTER, DIMENSION(:,:) ::   vice_init, vice_final   ! ice volume summed over categories
      REAL(wp), POINTER, DIMENSION(:,:) ::   eice_init, eice_final   ! ice energy summed over layers

      REAL(wp), POINTER, DIMENSION(:,:,:) ::   aicen_init, vicen_init   ! ice  area    & volume before ridging
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   vsnwn_init, esnwn_init   ! snow volume  & energy before ridging
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   smv_i_init, oa_i_init    ! ice salinity & age    before ridging

      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   eicen_init        ! ice energy before ridging

      REAL(wp), POINTER, DIMENSION(:,:) ::   afrac , fvol     ! fraction of category area ridged & new ridge volume going to n2
      REAL(wp), POINTER, DIMENSION(:,:) ::   ardg1 , ardg2    ! area of ice ridged & new ridges
      REAL(wp), POINTER, DIMENSION(:,:) ::   vsrdg , esrdg    ! snow volume & energy of ridging ice
      REAL(wp), POINTER, DIMENSION(:,:) ::   dhr   , dhr2     ! hrmax - hrmin  &  hrmax^2 - hrmin^2

      REAL(wp), POINTER, DIMENSION(:,:) ::   vrdg1   ! volume of ice ridged
      REAL(wp), POINTER, DIMENSION(:,:) ::   vrdg2   ! volume of new ridges
      REAL(wp), POINTER, DIMENSION(:,:) ::   vsw     ! volume of seawater trapped into ridges
      REAL(wp), POINTER, DIMENSION(:,:) ::   srdg1   ! sal*volume of ice ridged
      REAL(wp), POINTER, DIMENSION(:,:) ::   srdg2   ! sal*volume of new ridges
      REAL(wp), POINTER, DIMENSION(:,:) ::   smsw    ! sal*volume of water trapped into ridges
      REAL(wp), POINTER, DIMENSION(:,:) ::   oirdg1, oirdg2   ! ice age of ice ridged

      REAL(wp), POINTER, DIMENSION(:,:) ::   afrft            ! fraction of category area rafted
      REAL(wp), POINTER, DIMENSION(:,:) ::   arft1 , arft2    ! area of ice rafted and new rafted zone
      REAL(wp), POINTER, DIMENSION(:,:) ::   virft , vsrft    ! ice & snow volume of rafting ice
      REAL(wp), POINTER, DIMENSION(:,:) ::   esrft , smrft    ! snow energy & salinity of rafting ice
      REAL(wp), POINTER, DIMENSION(:,:) ::   oirft1, oirft2   ! ice age of ice rafted

      REAL(wp), POINTER, DIMENSION(:,:,:) ::   eirft      ! ice energy of rafting ice
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   erdg1      ! enth*volume of ice ridged
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   erdg2      ! enth*volume of new ridges
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   ersw       ! enth of water trapped into ridges
      !!----------------------------------------------------------------------

      CALL wrk_alloc( (jpi+1)*(jpj+1),       indxi, indxj )
      CALL wrk_alloc( jpi, jpj,              vice_init, vice_final, eice_init, eice_final )
      CALL wrk_alloc( jpi, jpj,              afrac, fvol , ardg1, ardg2, vsrdg, esrdg, dhr, dhr2 )
      CALL wrk_alloc( jpi, jpj,              vrdg1, vrdg2, vsw  , srdg1, srdg2, smsw, oirdg1, oirdg2 )
      CALL wrk_alloc( jpi, jpj,              afrft, arft1, arft2, virft, vsrft, esrft, smrft, oirft1, oirft2 )
      CALL wrk_alloc( jpi, jpj, jpl,         aicen_init, vicen_init, vsnwn_init, esnwn_init, smv_i_init, oa_i_init )
      CALL wrk_alloc( jpi, jpj, nlay_i,      eirft, erdg1, erdg2, ersw )
      CALL wrk_alloc( jpi, jpj, nlay_i, jpl, eicen_init )

      ! Conservation check
      eice_init(:,:) = 0._wp

      IF( con_i ) THEN
         CALL lim_column_sum        (jpl,    v_i,       vice_init )
         CALL lim_column_sum_energy (jpl, nlay_i,  e_i, eice_init )
         DO ji = mi0(iiceprt), mi1(iiceprt)
            DO jj = mj0(jiceprt), mj1(jiceprt)
               WRITE(numout,*) ' vice_init  : ', vice_init(ji,jj)
               WRITE(numout,*) ' eice_init  : ', eice_init(ji,jj)
            END DO
         END DO
      ENDIF

      !-------------------------------------------------------------------------------
      ! 1) Compute change in open water area due to closing and opening.
      !-------------------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            ato_i(ji,jj) = ato_i(ji,jj) - athorn(ji,jj,0) * closing_gross(ji,jj) * rdt_ice        &
               &                        + opning(ji,jj)                          * rdt_ice
            IF    ( ato_i(ji,jj) < -epsi10 ) THEN    ! there is a bug
               IF(lwp)   WRITE(numout,*) 'Ridging error: ato_i < 0 -- ato_i : ',ato_i(ji,jj)
            ELSEIF( ato_i(ji,jj) < 0._wp   ) THEN    ! roundoff error
               ato_i(ji,jj) = 0._wp
            ENDIF
         END DO
      END DO

      !-----------------------------------------------------------------
      ! 2) Save initial state variables
      !-----------------------------------------------------------------
      aicen_init(:,:,:)   = a_i  (:,:,:)
      vicen_init(:,:,:)   = v_i  (:,:,:)
      vsnwn_init(:,:,:)   = v_s  (:,:,:)
      smv_i_init(:,:,:)   = smv_i(:,:,:)
      esnwn_init(:,:,:)   = e_s  (:,:,1,:)
      eicen_init(:,:,:,:) = e_i  (:,:,:,:)
      oa_i_init (:,:,:)   = oa_i (:,:,:)

      !
      !-----------------------------------------------------------------
      ! 3) Pump everything from ice which is being ridged / rafted
      !-----------------------------------------------------------------
      ! Compute the area, volume, and energy of ice ridging in each
      ! category, along with the area of the resulting ridge.

      DO jl1 = 1, jpl !jl1 describes the ridging category

         !------------------------------------------------
         ! 3.1) Identify grid cells with nonzero ridging
         !------------------------------------------------

         icells = 0
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( aicen_init(ji,jj,jl1) > epsi10 .AND. athorn(ji,jj,jl1) > 0._wp  &
                  &   .AND. closing_gross(ji,jj) > 0._wp ) THEN
                  icells = icells + 1
                  indxi(icells) = ji
                  indxj(icells) = jj
               ENDIF
            END DO
         END DO

         DO ij = 1, icells
            ji = indxi(ij)
            jj = indxj(ij)

            !--------------------------------------------------------------------
            ! 3.2) Compute area of ridging ice (ardg1) and of new ridge (ardg2)
            !--------------------------------------------------------------------

            ardg1(ji,jj) = aridge(ji,jj,jl1)*closing_gross(ji,jj)*rdt_ice
            arft1(ji,jj) = araft (ji,jj,jl1)*closing_gross(ji,jj)*rdt_ice
            ardg2(ji,jj) = ardg1(ji,jj) / krdg(ji,jj,jl1)
            arft2(ji,jj) = arft1(ji,jj) / kraft

            !---------------------------------------------------------------
            ! 3.3) Compute ridging /rafting fractions, make sure afrac <=1 
            !---------------------------------------------------------------

            afrac(ji,jj) = ardg1(ji,jj) / aicen_init(ji,jj,jl1) !ridging
            afrft(ji,jj) = arft1(ji,jj) / aicen_init(ji,jj,jl1) !rafting

            IF( afrac(ji,jj) > kamax + epsi10 ) THEN  ! there is a bug
               IF(lwp)   WRITE(numout,*) ' ardg > a_i -- ardg, aicen_init : ', ardg1(ji,jj), aicen_init(ji,jj,jl1)
            ELSEIF( afrac(ji,jj) > kamax ) THEN       ! roundoff error
               afrac(ji,jj) = kamax
            ENDIF

            IF( afrft(ji,jj) > kamax + epsi10 ) THEN ! there is a bug
               IF(lwp)   WRITE(numout,*) ' arft > a_i -- arft, aicen_init : ', arft1(ji,jj), aicen_init(ji,jj,jl1) 
            ELSEIF( afrft(ji,jj) > kamax) THEN       ! roundoff error
               afrft(ji,jj) = kamax
            ENDIF

            !--------------------------------------------------------------------------
            ! 3.4) Subtract area, volume, and energy from ridging 
            !     / rafting category n1.
            !--------------------------------------------------------------------------
            vrdg1(ji,jj) = vicen_init(ji,jj,jl1) * afrac(ji,jj)
            vrdg2(ji,jj) = vrdg1(ji,jj) * ( 1. + rn_por_rdg )
            vsw  (ji,jj) = vrdg1(ji,jj) * rn_por_rdg

            vsrdg (ji,jj) = vsnwn_init(ji,jj,jl1) * afrac(ji,jj)
            esrdg (ji,jj) = esnwn_init(ji,jj,jl1) * afrac(ji,jj)
            srdg1 (ji,jj) = smv_i_init(ji,jj,jl1) * afrac(ji,jj)
            oirdg1(ji,jj) = oa_i_init (ji,jj,jl1) * afrac(ji,jj)
            oirdg2(ji,jj) = oa_i_init (ji,jj,jl1) * afrac(ji,jj) / krdg(ji,jj,jl1) 

            ! rafting volumes, heat contents ...
            virft (ji,jj) = vicen_init(ji,jj,jl1) * afrft(ji,jj)
            vsrft (ji,jj) = vsnwn_init(ji,jj,jl1) * afrft(ji,jj)
            esrft (ji,jj) = esnwn_init(ji,jj,jl1) * afrft(ji,jj)
            smrft (ji,jj) = smv_i_init(ji,jj,jl1) * afrft(ji,jj) 
            oirft1(ji,jj) = oa_i_init (ji,jj,jl1) * afrft(ji,jj) 
            oirft2(ji,jj) = oa_i_init (ji,jj,jl1) * afrft(ji,jj) / kraft 

            ! substract everything
            a_i(ji,jj,jl1)   = a_i(ji,jj,jl1)   - ardg1 (ji,jj) - arft1 (ji,jj)
            v_i(ji,jj,jl1)   = v_i(ji,jj,jl1)   - vrdg1 (ji,jj) - virft (ji,jj)
            v_s(ji,jj,jl1)   = v_s(ji,jj,jl1)   - vsrdg (ji,jj) - vsrft (ji,jj)
            e_s(ji,jj,1,jl1) = e_s(ji,jj,1,jl1) - esrdg (ji,jj) - esrft (ji,jj)
            smv_i(ji,jj,jl1) = smv_i(ji,jj,jl1) - srdg1 (ji,jj) - smrft (ji,jj)
            oa_i(ji,jj,jl1)  = oa_i(ji,jj,jl1)  - oirdg1(ji,jj) - oirft1(ji,jj)

            !-----------------------------------------------------------------
            ! 3.5) Compute properties of new ridges
            !-----------------------------------------------------------------
            !---------
            ! Salinity
            !---------
            smsw(ji,jj)  = vsw(ji,jj) * sss_m(ji,jj)                      ! salt content of seawater frozen in voids !! MV HC2014
            srdg2(ji,jj) = srdg1(ji,jj) + smsw(ji,jj)                     ! salt content of new ridge

            !srdg2(ji,jj) = MIN( rn_simax * vrdg2(ji,jj) , zsrdg2 )         ! impose a maximum salinity
            
            sfx_dyn(ji,jj) = sfx_dyn(ji,jj) - smsw(ji,jj) * rhoic * r1_rdtice
            wfx_dyn(ji,jj) = wfx_dyn(ji,jj) - vsw (ji,jj) * rhoic * r1_rdtice   ! increase in ice volume du to seawater frozen in voids             

            !------------------------------------            
            ! 3.6 Increment ridging diagnostics
            !------------------------------------            

            !        jl1 looping 1-jpl
            !           ij looping 1-icells

            dardg1dt   (ji,jj) = dardg1dt(ji,jj) + ardg1(ji,jj) + arft1(ji,jj)
            dardg2dt   (ji,jj) = dardg2dt(ji,jj) + ardg2(ji,jj) + arft2(ji,jj)
            opening    (ji,jj) = opening (ji,jj) + opning(ji,jj) * rdt_ice

            IF( con_i )   vice_init(ji,jj) = vice_init(ji,jj) + vrdg2(ji,jj) - vrdg1(ji,jj)

            !------------------------------------------            
            ! 3.7 Put the snow somewhere in the ocean
            !------------------------------------------            
            !  Place part of the snow lost by ridging into the ocean. 
            !  Note that esnow_mlt < 0; the ocean must cool to melt snow.
            !  If the ocean temp = Tf already, new ice must grow.
            !  During the next time step, thermo_rates will determine whether
            !  the ocean cools or new ice grows.
            !        jl1 looping 1-jpl
            !           ij looping 1-icells

            msnow_mlt(ji,jj) = msnow_mlt(ji,jj) + rhosn*vsrdg(ji,jj)*(1.0-rn_fsnowrdg)   &   ! rafting included
               &                                + rhosn*vsrft(ji,jj)*(1.0-rn_fsnowrft)

            ! in J/m2 (same as e_s)
            esnow_mlt(ji,jj) = esnow_mlt(ji,jj) - esrdg(ji,jj)*(1.0-rn_fsnowrdg)         &   !rafting included
               &                                - esrft(ji,jj)*(1.0-rn_fsnowrft)          

            !-----------------------------------------------------------------
            ! 3.8 Compute quantities used to apportion ice among categories
            ! in the n2 loop below
            !-----------------------------------------------------------------

            !        jl1 looping 1-jpl
            !           ij looping 1-icells

            dhr (ji,jj) = hrmax(ji,jj,jl1) - hrmin(ji,jj,jl1)
            dhr2(ji,jj) = hrmax(ji,jj,jl1) * hrmax(ji,jj,jl1) - hrmin(ji,jj,jl1) * hrmin(ji,jj,jl1)

         END DO

         !--------------------------------------------------------------------
         ! 3.9 Compute ridging ice enthalpy, remove it from ridging ice and
         !      compute ridged ice enthalpy 
         !--------------------------------------------------------------------
         DO jk = 1, nlay_i
            DO ij = 1, icells
               ji = indxi(ij)
               jj = indxj(ij)
               ! heat content of ridged ice
               erdg1(ji,jj,jk)      = eicen_init(ji,jj,jk,jl1) * afrac(ji,jj) 
               eirft(ji,jj,jk)      = eicen_init(ji,jj,jk,jl1) * afrft(ji,jj)
               e_i  (ji,jj,jk,jl1)  = e_i(ji,jj,jk,jl1) - erdg1(ji,jj,jk) - eirft(ji,jj,jk)
               
               
               ! enthalpy of the trapped seawater (J/m2, >0)
               ! clem: if sst>0, then ersw <0 (is that possible?)
               ersw(ji,jj,jk)   = - rhoic * vsw(ji,jj) * rcp * sst_m(ji,jj) * r1_nlay_i

               ! heat flux to the ocean
               hfx_dyn(ji,jj) = hfx_dyn(ji,jj) + ersw(ji,jj,jk) * r1_rdtice  ! > 0 [W.m-2] ocean->ice flux 

               ! it is added to sea ice because the sign convention is the opposite of the sign convention for the ocean
               erdg2(ji,jj,jk)  = erdg1(ji,jj,jk) + ersw(ji,jj,jk)

            END DO
         END DO


         IF( con_i ) THEN
            DO jk = 1, nlay_i
               DO ij = 1, icells
                  ji = indxi(ij)
                  jj = indxj(ij)
                  eice_init(ji,jj) = eice_init(ji,jj) + erdg2(ji,jj,jk) - erdg1(ji,jj,jk)
               END DO
            END DO
         ENDIF

         !-------------------------------------------------------------------------------
         ! 4) Add area, volume, and energy of new ridge to each category jl2
         !-------------------------------------------------------------------------------
         !        jl1 looping 1-jpl
         DO jl2  = 1, jpl 
            ! over categories to which ridged ice is transferred
            DO ij = 1, icells
               ji = indxi(ij)
               jj = indxj(ij)

               ! Compute the fraction of ridged ice area and volume going to 
               ! thickness category jl2.
               ! Transfer area, volume, and energy accordingly.

               IF( hrmin(ji,jj,jl1) >= hi_max(jl2) .OR. hrmax(ji,jj,jl1) <= hi_max(jl2-1) ) THEN
                  hL = 0._wp
                  hR = 0._wp
               ELSE
                  hL = MAX( hrmin(ji,jj,jl1), hi_max(jl2-1) )
                  hR = MIN( hrmax(ji,jj,jl1), hi_max(jl2)   )
               ENDIF

               ! fraction of ridged ice area and volume going to n2
               farea = ( hR - hL ) / dhr(ji,jj) 
               fvol(ji,jj) = ( hR*hR - hL*hL ) / dhr2(ji,jj)

               a_i  (ji,jj  ,jl2) = a_i  (ji,jj  ,jl2) + ardg2 (ji,jj) * farea
               v_i  (ji,jj  ,jl2) = v_i  (ji,jj  ,jl2) + vrdg2 (ji,jj) * fvol(ji,jj)
               v_s  (ji,jj  ,jl2) = v_s  (ji,jj  ,jl2) + vsrdg (ji,jj) * fvol(ji,jj) * rn_fsnowrdg
               e_s  (ji,jj,1,jl2) = e_s  (ji,jj,1,jl2) + esrdg (ji,jj) * fvol(ji,jj) * rn_fsnowrdg
               smv_i(ji,jj  ,jl2) = smv_i(ji,jj  ,jl2) + srdg2 (ji,jj) * fvol(ji,jj)
               oa_i (ji,jj  ,jl2) = oa_i (ji,jj  ,jl2) + oirdg2(ji,jj) * farea

            END DO

            ! Transfer ice energy to category jl2 by ridging
            DO jk = 1, nlay_i
               DO ij = 1, icells
                  ji = indxi(ij)
                  jj = indxj(ij)
                  e_i(ji,jj,jk,jl2) = e_i(ji,jj,jk,jl2) + fvol(ji,jj) * erdg2(ji,jj,jk)
               END DO
            END DO
            !
         END DO                 ! jl2 (new ridges)            

         DO jl2 = 1, jpl 

            DO ij = 1, icells
               ji = indxi(ij)
               jj = indxj(ij)
               ! Compute the fraction of rafted ice area and volume going to 
               ! thickness category jl2, transfer area, volume, and energy accordingly.
               !
               IF( hraft(ji,jj,jl1) <= hi_max(jl2) .AND. hraft(ji,jj,jl1) >  hi_max(jl2-1) ) THEN
                  a_i  (ji,jj  ,jl2) = a_i  (ji,jj  ,jl2) + arft2 (ji,jj)
                  v_i  (ji,jj  ,jl2) = v_i  (ji,jj  ,jl2) + virft (ji,jj)
                  v_s  (ji,jj  ,jl2) = v_s  (ji,jj  ,jl2) + vsrft (ji,jj) * rn_fsnowrft
                  e_s  (ji,jj,1,jl2) = e_s  (ji,jj,1,jl2) + esrft (ji,jj) * rn_fsnowrft
                  smv_i(ji,jj  ,jl2) = smv_i(ji,jj  ,jl2) + smrft (ji,jj)    
                  oa_i (ji,jj  ,jl2) = oa_i (ji,jj  ,jl2) + oirft2(ji,jj)
               ENDIF
               !
            END DO

            ! Transfer rafted ice energy to category jl2 
            DO jk = 1, nlay_i
               DO ij = 1, icells
                  ji = indxi(ij)
                  jj = indxj(ij)
                  IF( hraft(ji,jj,jl1) <= hi_max(jl2) .AND. hraft(ji,jj,jl1) > hi_max(jl2-1)  ) THEN
                     e_i(ji,jj,jk,jl2) = e_i(ji,jj,jk,jl2) + eirft(ji,jj,jk)
                  ENDIF
               END DO
            END DO

         END DO

      END DO ! jl1 (deforming categories)

      ! Conservation check
      IF ( con_i ) THEN
         CALL lim_column_sum (jpl,   v_i, vice_final)
         fieldid = ' v_i : limitd_me '
         CALL lim_cons_check (vice_init, vice_final, 1.0e-6, fieldid) 

         CALL lim_column_sum_energy (jpl, nlay_i,  e_i, eice_final )
         fieldid = ' e_i : limitd_me '
         CALL lim_cons_check (eice_init, eice_final, 1.0e-2, fieldid) 

         DO ji = mi0(iiceprt), mi1(iiceprt)
            DO jj = mj0(jiceprt), mj1(jiceprt)
               WRITE(numout,*) ' vice_init  : ', vice_init (ji,jj)
               WRITE(numout,*) ' vice_final : ', vice_final(ji,jj)
               WRITE(numout,*) ' eice_init  : ', eice_init (ji,jj)
               WRITE(numout,*) ' eice_final : ', eice_final(ji,jj)
            END DO
         END DO
      ENDIF
      !
      CALL wrk_dealloc( (jpi+1)*(jpj+1),        indxi, indxj )
      CALL wrk_dealloc( jpi, jpj,               vice_init, vice_final, eice_init, eice_final )
      CALL wrk_dealloc( jpi, jpj,               afrac, fvol , ardg1, ardg2, vsrdg, esrdg, dhr, dhr2 )
      CALL wrk_dealloc( jpi, jpj,               vrdg1, vrdg2, vsw  , srdg1, srdg2, smsw, oirdg1, oirdg2 )
      CALL wrk_dealloc( jpi, jpj,               afrft, arft1, arft2, virft, vsrft, esrft, smrft, oirft1, oirft2 )
      CALL wrk_dealloc( jpi, jpj, jpl,          aicen_init, vicen_init, vsnwn_init, esnwn_init, smv_i_init, oa_i_init )
      CALL wrk_dealloc( jpi, jpj, nlay_i,       eirft, erdg1, erdg2, ersw )
      CALL wrk_dealloc( jpi, jpj, nlay_i, jpl,  eicen_init )
      !
   END SUBROUTINE lim_itd_me_ridgeshift

   SUBROUTINE lim_itd_me_init
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE lim_itd_me_init ***
      !!
      !! ** Purpose :   Physical constants and parameters linked 
      !!                to the mechanical ice redistribution
      !!
      !! ** Method  :   Read the namiceitdme namelist 
      !!                and check the parameters values 
      !!                called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namiceitdme
      !!-------------------------------------------------------------------
      INTEGER :: ios                 ! Local integer output status for namelist read
      NAMELIST/namiceitdme/ rn_cs, rn_fsnowrdg, rn_fsnowrft,              & 
        &                   rn_gstar, rn_astar, rn_hstar, ln_rafting, rn_hraft, rn_craft, rn_por_rdg, &
        &                   nn_partfun
      !!-------------------------------------------------------------------
      !
      REWIND( numnam_ice_ref )              ! Namelist namicetdme in reference namelist : Ice mechanical ice redistribution
      READ  ( numnam_ice_ref, namiceitdme, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namiceitdme in reference namelist', lwp )

      REWIND( numnam_ice_cfg )              ! Namelist namiceitdme in configuration namelist : Ice mechanical ice redistribution
      READ  ( numnam_ice_cfg, namiceitdme, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namiceitdme in configuration namelist', lwp )
      IF(lwm) WRITE ( numoni, namiceitdme )
      !
      IF (lwp) THEN                          ! control print
         WRITE(numout,*)
         WRITE(numout,*)' lim_itd_me_init : ice parameters for mechanical ice redistribution '
         WRITE(numout,*)' ~~~~~~~~~~~~~~~'
         WRITE(numout,*)'   Fraction of shear energy contributing to ridging        rn_cs       = ', rn_cs 
         WRITE(numout,*)'   Fraction of snow volume conserved during ridging        rn_fsnowrdg = ', rn_fsnowrdg 
         WRITE(numout,*)'   Fraction of snow volume conserved during ridging        rn_fsnowrft = ', rn_fsnowrft 
         WRITE(numout,*)'   Fraction of total ice coverage contributing to ridging  rn_gstar    = ', rn_gstar
         WRITE(numout,*)'   Equivalent to G* for an exponential part function       rn_astar    = ', rn_astar
         WRITE(numout,*)'   Quantity playing a role in max ridged ice thickness     rn_hstar    = ', rn_hstar
         WRITE(numout,*)'   Rafting of ice sheets or not                            ln_rafting  = ', ln_rafting
         WRITE(numout,*)'   Parmeter thickness (threshold between ridge-raft)       rn_hraft    = ', rn_hraft
         WRITE(numout,*)'   Rafting hyperbolic tangent coefficient                  rn_craft    = ', rn_craft  
         WRITE(numout,*)'   Initial porosity of ridges                              rn_por_rdg  = ', rn_por_rdg
         WRITE(numout,*)'   Switch for part. function (0) linear (1) exponential    nn_partfun  = ', nn_partfun
      ENDIF
      !
   END SUBROUTINE lim_itd_me_init

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module          NO LIM-3 sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_itd_me           ! Empty routines
   END SUBROUTINE lim_itd_me
   SUBROUTINE lim_itd_me_icestrength
   END SUBROUTINE lim_itd_me_icestrength
   SUBROUTINE lim_itd_me_init
   END SUBROUTINE lim_itd_me_init
#endif
   !!======================================================================
END MODULE limitd_me
