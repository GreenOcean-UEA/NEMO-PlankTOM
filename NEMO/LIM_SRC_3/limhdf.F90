MODULE limhdf
   !!======================================================================
   !!                    ***  MODULE limhdf   ***
   !! LIM ice model : horizontal diffusion of sea-ice quantities
   !!======================================================================
   !! History :  LIM  !  2000-01 (LIM) Original code
   !!             -   !  2001-05 (G. Madec, R. Hordoir) opa norm
   !!            1.0  !  2002-08 (C. Ethe)  F90, free form
   !!----------------------------------------------------------------------
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3'                                      LIM3 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_hdf       : diffusion trend on sea-ice variable
   !!   lim_hdf_init  : initialisation of diffusion trend on sea-ice variable
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean domain
   USE ice            ! LIM-3: ice variables
   USE lbclnk         ! lateral boundary condition - MPP exchanges
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! work arrays
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_hdf         ! called by lim_trp
   PUBLIC   lim_hdf_init    ! called by sbc_lim_init

   LOGICAL  ::   linit = .TRUE.                             ! initialization flag (set to flase after the 1st call)
   INTEGER  ::   nn_convfrq                                 !:  convergence check frequency of the Crant-Nicholson scheme
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   efact   ! metric coefficient

   !! * Substitution 
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/LIM3 4.0 , UCL - NEMO Consortium (2010)
   !! $Id: limhdf.F90 5429 2015-06-16 09:57:07Z smasson $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_hdf( ptab )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE lim_hdf  ***
      !!
      !! ** purpose :   Compute and add the diffusive trend on sea-ice variables
      !!
      !! ** method  :   Second order diffusive operator evaluated using a
      !!              Cranck-Nicholson time Scheme.
      !!
      !! ** Action  :    update ptab with the diffusive contribution
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   ptab    ! Field on which the diffusion is applied
      !
      INTEGER                           ::  ji, jj                    ! dummy loop indices
      INTEGER                           ::  iter, ierr           ! local integers
      REAL(wp)                          ::  zrlxint, zconv     ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::  zrlx, zflu, zflv, zdiv0, zdiv, ztab0
      CHARACTER(lc)                     ::  charout                   ! local character
      REAL(wp), PARAMETER               ::  zrelax = 0.5_wp           ! relaxation constant for iterative procedure
      REAL(wp), PARAMETER               ::  zalfa  = 0.5_wp           ! =1.0/0.5/0.0 = implicit/Cranck-Nicholson/explicit
      INTEGER , PARAMETER               ::  its    = 100              ! Maximum number of iteration
      !!-------------------------------------------------------------------
      
      CALL wrk_alloc( jpi, jpj, zrlx, zflu, zflv, zdiv0, zdiv, ztab0 )

      !                       !==  Initialisation  ==!
      !
      IF( linit ) THEN              ! Metric coefficient (compute at the first call and saved in efact)
         ALLOCATE( efact(jpi,jpj) , STAT=ierr )
         IF( lk_mpp    )   CALL mpp_sum( ierr )
         IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'lim_hdf : unable to allocate arrays' )
         DO jj = 2, jpjm1  
            DO ji = fs_2 , fs_jpim1   ! vector opt.
               efact(ji,jj) = ( e2u(ji,jj) + e2u(ji-1,jj) + e1v(ji,jj) + e1v(ji,jj-1) ) * r1_e12t(ji,jj)
            END DO
         END DO
         linit = .FALSE.
      ENDIF
      !                             ! Time integration parameters
      !
      ztab0(:, : ) = ptab(:,:)      ! Arrays initialization
      zdiv0(:, 1 ) = 0._wp
      zdiv0(:,jpj) = 0._wp
      zflu (jpi,:) = 0._wp   
      zflv (jpi,:) = 0._wp
      zdiv0(1,  :) = 0._wp
      zdiv0(jpi,:) = 0._wp

      zconv = 1._wp           !==  horizontal diffusion using a Crant-Nicholson scheme  ==!
      iter  = 0
      !
      DO WHILE( zconv > ( 2._wp * 1.e-04 ) .AND. iter <= its )   ! Sub-time step loop
         !
         iter = iter + 1                                 ! incrementation of the sub-time step number
         !
         DO jj = 1, jpjm1                                ! diffusive fluxes in U- and V- direction
            DO ji = 1 , fs_jpim1   ! vector opt.
               zflu(ji,jj) = pahu(ji,jj) * e2u(ji,jj) * r1_e1u(ji,jj) * ( ptab(ji+1,jj) - ptab(ji,jj) )
               zflv(ji,jj) = pahv(ji,jj) * e1v(ji,jj) * r1_e2v(ji,jj) * ( ptab(ji,jj+1) - ptab(ji,jj) )
            END DO
         END DO
         !
         DO jj= 2, jpjm1                                 ! diffusive trend : divergence of the fluxes
            DO ji = fs_2 , fs_jpim1   ! vector opt. 
               zdiv(ji,jj) = ( zflu(ji,jj) - zflu(ji-1,jj) + zflv(ji,jj) - zflv(ji,jj-1) ) * r1_e12t(ji,jj)
            END DO
         END DO
         !
         IF( iter == 1 )   zdiv0(:,:) = zdiv(:,:)        ! save the 1st evaluation of the diffusive trend in zdiv0
         !
         DO jj = 2, jpjm1                                ! iterative evaluation
            DO ji = fs_2 , fs_jpim1   ! vector opt.
               zrlxint = (   ztab0(ji,jj)    &
                  &       +  rdt_ice * (           zalfa   * ( zdiv(ji,jj) + efact(ji,jj) * ptab(ji,jj) )   &
                  &                      + ( 1.0 - zalfa ) *   zdiv0(ji,jj) )                               & 
                  &      ) / ( 1.0 + zalfa * rdt_ice * efact(ji,jj) )
               zrlx(ji,jj) = ptab(ji,jj) + zrelax * ( zrlxint - ptab(ji,jj) )
            END DO
         END DO
         CALL lbc_lnk( zrlx, 'T', 1. )                   ! lateral boundary condition
         !
         IF ( MOD( iter, nn_convfrq ) == 0 )  THEN    ! convergence test every nn_convfrq iterations (perf. optimization)
            zconv = 0._wp
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  zconv = MAX( zconv, ABS( zrlx(ji,jj) - ptab(ji,jj) )  )
               END DO
            END DO
            IF( lk_mpp )   CALL mpp_max( zconv )      ! max over the global domain
         ENDIF
         !
         ptab(:,:) = zrlx(:,:)
         !
      END DO                                       ! end of sub-time step loop

      ! -----------------------
      !!! final step (clem) !!!
      DO jj = 1, jpjm1                                ! diffusive fluxes in U- and V- direction
         DO ji = 1 , fs_jpim1   ! vector opt.
            zflu(ji,jj) = pahu(ji,jj) * e2u(ji,jj) * r1_e1u(ji,jj) * ( ptab(ji+1,jj) - ptab(ji,jj) )
            zflv(ji,jj) = pahv(ji,jj) * e1v(ji,jj) * r1_e2v(ji,jj) * ( ptab(ji,jj+1) - ptab(ji,jj) )
         END DO
      END DO
      !
      DO jj= 2, jpjm1                                 ! diffusive trend : divergence of the fluxes
         DO ji = fs_2 , fs_jpim1   ! vector opt. 
            zdiv(ji,jj) = ( zflu(ji,jj) - zflu(ji-1,jj) + zflv(ji,jj) - zflv(ji,jj-1) ) * r1_e12t(ji,jj)
            ptab(ji,jj) = ztab0(ji,jj) + 0.5 * ( zdiv(ji,jj) + zdiv0(ji,jj) )
         END DO
      END DO
      CALL lbc_lnk( ptab, 'T', 1. )                   ! lateral boundary condition
      !!! final step (clem) !!!
      ! -----------------------

      IF(ln_ctl)   THEN
         zrlx(:,:) = ptab(:,:) - ztab0(:,:)
         WRITE(charout,FMT="(' lim_hdf  : zconv =',D23.16, ' iter =',I4,2X)") zconv, iter
         CALL prt_ctl( tab2d_1=zrlx, clinfo1=charout )
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, zrlx, zflu, zflv, zdiv0, zdiv, ztab0 )
      !
   END SUBROUTINE lim_hdf

   
   SUBROUTINE lim_hdf_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE lim_hdf_init  ***
      !!
      !! ** Purpose : Initialisation of horizontal diffusion of sea-ice 
      !!
      !! ** Method  : Read the namicehdf namelist
      !!
      !! ** input   : Namelist namicehdf
      !!-------------------------------------------------------------------
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      NAMELIST/namicehdf/ nn_convfrq
      !!-------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'lim_hdf : Ice horizontal diffusion'
         WRITE(numout,*) '~~~~~~~'
      ENDIF
      !
      REWIND( numnam_ice_ref )              ! Namelist namicehdf in reference namelist : Ice horizontal diffusion
      READ  ( numnam_ice_ref, namicehdf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namicehdf in reference namelist', lwp )

      REWIND( numnam_ice_cfg )              ! Namelist namicehdf in configuration namelist : Ice horizontal diffusion
      READ  ( numnam_ice_cfg, namicehdf, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namicehdf in configuration namelist', lwp )
      IF(lwm) WRITE ( numoni, namicehdf )
      !
      IF(lwp) THEN                          ! control print
         WRITE(numout,*)
         WRITE(numout,*)'   Namelist of ice parameters for ice horizontal diffusion computation '
         WRITE(numout,*)'      convergence check frequency of the Crant-Nicholson scheme   nn_convfrq   = ', nn_convfrq
      ENDIF
      !
   END SUBROUTINE lim_hdf_init
#else
   !!----------------------------------------------------------------------
   !!   Default option          Dummy module           NO LIM sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE limhdf
