MODULE limrhg
   !!======================================================================
   !!                     ***  MODULE  limrhg  ***
   !!   Ice rheology : sea ice rheology
   !!======================================================================
   !! History :   -   !  2007-03  (M.A. Morales Maqueda, S. Bouillon) Original code
   !!            3.0  !  2008-03  (M. Vancoppenolle) LIM3
   !!             -   !  2008-11  (M. Vancoppenolle, S. Bouillon, Y. Aksenov) add surface tilt in ice rheolohy 
   !!            3.3  !  2009-05  (G.Garric) addition of the lim2_evp cas
   !!            3.4  !  2011-01  (A. Porter)  dynamical allocation 
   !!            3.5  !  2012-08  (R. Benshila)  AGRIF 
   !!----------------------------------------------------------------------
#if defined key_lim3 || (  defined key_lim2 && ! defined key_lim2_vp )
   !!----------------------------------------------------------------------
   !!   'key_lim3'               OR                     LIM-3 sea-ice model
   !!   'key_lim2' AND NOT 'key_lim2_vp'            EVP LIM-2 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_rhg       : computes ice velocities
   !!----------------------------------------------------------------------
   USE phycst         ! Physical constant
   USE oce     , ONLY :  snwice_mass, snwice_mass_b
   USE par_oce        ! Ocean parameters
   USE dom_oce        ! Ocean domain
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbc_ice        ! Surface boundary condition: ice fields
#if defined key_lim3
   USE ice            ! LIM-3: ice variables
   USE dom_ice        ! LIM-3: ice domain
   USE limitd_me      ! LIM-3: 
#else
   USE ice_2          ! LIM-2: ice variables
   USE dom_ice_2      ! LIM-2: ice domain
#endif
   USE lbclnk         ! Lateral Boundary Condition / MPP link
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! work arrays
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  
#if defined key_agrif && defined key_lim2
   USE agrif_lim2_interp
#endif
#if defined key_bdy
   USE bdyice_lim
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_rhg        ! routine called by lim_dyn (or lim_dyn_2)

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/LIM3 4.0 , UCL - NEMO Consortium (2011)
   !! $Id: limrhg.F90 5429 2015-06-16 09:57:07Z smasson $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_rhg( k_j1, k_jpj )
      !!-------------------------------------------------------------------
      !!                 ***  SUBROUTINE lim_rhg  ***
      !!                          EVP-C-grid
      !!
      !! ** purpose : determines sea ice drift from wind stress, ice-ocean
      !!  stress and sea-surface slope. Ice-ice interaction is described by 
      !!  a non-linear elasto-viscous-plastic (EVP) law including shear 
      !!  strength and a bulk rheology (Hunke and Dukowicz, 2002).	
      !!
      !!  The points in the C-grid look like this, dear reader
      !!
      !!                              (ji,jj)
      !!                                 |
      !!                                 |
      !!                      (ji-1,jj)  |  (ji,jj)
      !!                             ---------   
      !!                            |         |
      !!                            | (ji,jj) |------(ji,jj)
      !!                            |         |
      !!                             ---------   
      !!                     (ji-1,jj-1)     (ji,jj-1)
      !!
      !! ** Inputs  : - wind forcing (stress), oceanic currents
      !!                ice total volume (vt_i) per unit area
      !!                snow total volume (vt_s) per unit area
      !!
      !! ** Action  : - compute u_ice, v_ice : the components of the 
      !!                sea-ice velocity vector
      !!              - compute delta_i, shear_i, divu_i, which are inputs
      !!                of the ice thickness distribution
      !!
      !! ** Steps   : 1) Compute ice snow mass, ice strength 
      !!              2) Compute wind, oceanic stresses, mass terms and
      !!                 coriolis terms of the momentum equation
      !!              3) Solve the momentum equation (iterative procedure)
      !!              4) Prevent high velocities if the ice is thin
      !!              5) Recompute invariants of the strain rate tensor
      !!                 which are inputs of the ITD, store stress
      !!                 for the next time step
      !!              6) Control prints of residual (convergence)
      !!                 and charge ellipse.
      !!                 The user should make sure that the parameters
      !!                 nn_nevp, elastic time scale and rn_creepl maintain stress state
      !!                 on the charge ellipse for plastic flow
      !!                 e.g. in the Canadian Archipelago
      !!
      !! References : Hunke and Dukowicz, JPO97
      !!              Bouillon et al., Ocean Modelling 2009
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   k_j1    ! southern j-index for ice computation
      INTEGER, INTENT(in) ::   k_jpj   ! northern j-index for ice computation
      !!
      INTEGER ::   ji, jj   ! dummy loop indices
      INTEGER ::   jter     ! local integers
      CHARACTER (len=50) ::   charout
      REAL(wp) ::   zt11, zt12, zt21, zt22, ztagnx, ztagny, delta                         !
      REAL(wp) ::   za, zstms          ! local scalars
      REAL(wp) ::   zc1, zc2, zc3      ! ice mass

      REAL(wp) ::   dtevp , z1_dtevp              ! time step for subcycling
      REAL(wp) ::   dtotel, z1_dtotel, ecc2, ecci ! square of yield ellipse eccenticity
      REAL(wp) ::   z0, zr, zcca, zccb            ! temporary scalars
      REAL(wp) ::   zu_ice2, zv_ice1              !
      REAL(wp) ::   zddc, zdtc                    ! delta on corners and on centre
      REAL(wp) ::   zdst                          ! shear at the center of the grid point
      REAL(wp) ::   zdsshx, zdsshy                ! term for the gradient of ocean surface
      REAL(wp) ::   sigma1, sigma2                ! internal ice stress

      REAL(wp) ::   zresm         ! Maximal error on ice velocity
      REAL(wp) ::   zintb, zintn  ! dummy argument

      REAL(wp), POINTER, DIMENSION(:,:) ::   zpresh           ! temporary array for ice strength
      REAL(wp), POINTER, DIMENSION(:,:) ::   zpreshc          ! Ice strength on grid cell corners (zpreshc)
      REAL(wp), POINTER, DIMENSION(:,:) ::   zfrld1, zfrld2   ! lead fraction on U/V points
      REAL(wp), POINTER, DIMENSION(:,:) ::   zmass1, zmass2   ! ice/snow mass on U/V points
      REAL(wp), POINTER, DIMENSION(:,:) ::   zcorl1, zcorl2   ! coriolis parameter on U/V points
      REAL(wp), POINTER, DIMENSION(:,:) ::   za1ct , za2ct    ! temporary arrays
      REAL(wp), POINTER, DIMENSION(:,:) ::   v_oce1           ! ocean u/v component on U points                           
      REAL(wp), POINTER, DIMENSION(:,:) ::   u_oce2           ! ocean u/v component on V points
      REAL(wp), POINTER, DIMENSION(:,:) ::   u_ice2, v_ice1   ! ice u/v component on V/U point
      REAL(wp), POINTER, DIMENSION(:,:) ::   zf1   , zf2      ! arrays for internal stresses
      REAL(wp), POINTER, DIMENSION(:,:) ::   zmask            ! mask ocean grid points
      
      REAL(wp), POINTER, DIMENSION(:,:) ::   zdt              ! tension at centre of grid cells
      REAL(wp), POINTER, DIMENSION(:,:) ::   zds              ! Shear on northeast corner of grid cells
      REAL(wp), POINTER, DIMENSION(:,:) ::   zs1   , zs2      ! Diagonal stress tensor components zs1 and zs2 
      REAL(wp), POINTER, DIMENSION(:,:) ::   zs12             ! Non-diagonal stress tensor component zs12
      REAL(wp), POINTER, DIMENSION(:,:) ::   zu_ice, zv_ice, zresr   ! Local error on velocity
      REAL(wp), POINTER, DIMENSION(:,:) ::   zpice            ! array used for the calculation of ice surface slope:
                                                              !   ocean surface (ssh_m) if ice is not embedded
                                                              !   ice top surface if ice is embedded   

      REAL(wp), PARAMETER               ::   zepsi = 1.0e-20_wp ! tolerance parameter
      REAL(wp), PARAMETER               ::   zvmin = 1.0e-03_wp ! ice volume below which ice velocity equals ocean velocity
      !!-------------------------------------------------------------------

      CALL wrk_alloc( jpi,jpj, zpresh, zfrld1, zmass1, zcorl1, za1ct , zpreshc, zfrld2, zmass2, zcorl2, za2ct )
      CALL wrk_alloc( jpi,jpj, u_oce2, u_ice2, v_oce1 , v_ice1 , zmask               )
      CALL wrk_alloc( jpi,jpj, zf1   , zu_ice, zf2   , zv_ice , zdt    , zds  )
      CALL wrk_alloc( jpi,jpj, zdt   , zds   , zs1   , zs2   , zs12   , zresr , zpice                 )

#if  defined key_lim2 && ! defined key_lim2_vp
# if defined key_agrif
      USE ice_2, vt_s => hsnm
      USE ice_2, vt_i => hicm
# else
      vt_s => hsnm
      vt_i => hicm
# endif
      at_i(:,:) = 1. - frld(:,:)
#endif
#if defined key_agrif && defined key_lim2 
      CALL agrif_rhg_lim2_load      ! First interpolation of coarse values
#endif
      !
      !------------------------------------------------------------------------------!
      ! 1) Ice strength (zpresh)                                !
      !------------------------------------------------------------------------------!
      !
      ! Put every vector to 0
      delta_i(:,:) = 0._wp   ;
      zpresh (:,:) = 0._wp   ;  
      zpreshc(:,:) = 0._wp
      u_ice2 (:,:) = 0._wp   ;   v_ice1(:,:) = 0._wp
      divu_i (:,:) = 0._wp   ;   zdt   (:,:) = 0._wp   ;   zds(:,:) = 0._wp
      shear_i(:,:) = 0._wp

#if defined key_lim3
      CALL lim_itd_me_icestrength( nn_icestr )      ! LIM-3: Ice strength on T-points
#endif

      DO jj = k_j1 , k_jpj       ! Ice mass and temp variables
         DO ji = 1 , jpi
#if defined key_lim3
            zpresh(ji,jj) = tmask(ji,jj,1) *  strength(ji,jj)
#endif
#if defined key_lim2
            zpresh(ji,jj) = tmask(ji,jj,1) *  pstar * vt_i(ji,jj) * EXP( -c_rhg * (1. - at_i(ji,jj) ) )
#endif
            ! zmask = 1 where there is ice or on land
            zmask(ji,jj)    = 1._wp - ( 1._wp - MAX( 0._wp , SIGN ( 1._wp , vt_i(ji,jj) - zepsi ) ) ) * tmask(ji,jj,1)
         END DO
      END DO

      ! Ice strength on grid cell corners (zpreshc)
      ! needed for calculation of shear stress 
      DO jj = k_j1+1, k_jpj-1
         DO ji = 2, jpim1 !RB caution no fs_ (ji+1,jj+1)
            zstms          =  tmask(ji+1,jj+1,1) * wght(ji+1,jj+1,2,2) + tmask(ji,jj+1,1) * wght(ji+1,jj+1,1,2) +   &
               &              tmask(ji+1,jj,1)   * wght(ji+1,jj+1,2,1) + tmask(ji,jj,1)   * wght(ji+1,jj+1,1,1)
            zpreshc(ji,jj) = ( zpresh(ji+1,jj+1) * wght(ji+1,jj+1,2,2) + zpresh(ji,jj+1) * wght(ji+1,jj+1,1,2) +   &
               &               zpresh(ji+1,jj)   * wght(ji+1,jj+1,2,1) + zpresh(ji,jj)   * wght(ji+1,jj+1,1,1)     &
               &             ) / MAX( zstms, zepsi )
         END DO
      END DO
      CALL lbc_lnk( zpreshc(:,:), 'F', 1. )
      !
      !------------------------------------------------------------------------------!
      ! 2) Wind / ocean stress, mass terms, coriolis terms
      !------------------------------------------------------------------------------!
      !
      !  Wind stress, coriolis and mass terms on the sides of the squares        
      !  zfrld1: lead fraction on U-points                                      
      !  zfrld2: lead fraction on V-points                                     
      !  zmass1: ice/snow mass on U-points                                    
      !  zmass2: ice/snow mass on V-points                                   
      !  zcorl1: Coriolis parameter on U-points                             
      !  zcorl2: Coriolis parameter on V-points                            
      !  (ztagnx,ztagny): wind stress on U/V points                       
      !  v_oce1: ocean v component on u points                          
      !  u_oce2: ocean u component on v points                         

      IF( nn_ice_embd == 2 ) THEN             !== embedded sea ice: compute representative ice top surface ==!
         !                                            
         ! average interpolation coeff as used in dynspg = (1/nn_fsbc) * {SUM[n/nn_fsbc], n=0,nn_fsbc-1}
         !                                               = (1/nn_fsbc)^2 * {SUM[n], n=0,nn_fsbc-1}
         zintn = REAL( nn_fsbc - 1 ) / REAL( nn_fsbc ) * 0.5_wp     
         !
         ! average interpolation coeff as used in dynspg = (1/nn_fsbc) * {SUM[1-n/nn_fsbc], n=0,nn_fsbc-1}
         !                                               = (1/nn_fsbc)^2 * (nn_fsbc^2 - {SUM[n], n=0,nn_fsbc-1})
         zintb = REAL( nn_fsbc + 1 ) / REAL( nn_fsbc ) * 0.5_wp
         !
         zpice(:,:) = ssh_m(:,:) + (  zintn * snwice_mass(:,:) +  zintb * snwice_mass_b(:,:)  ) * r1_rau0
         !
      ELSE                                    !== non-embedded sea ice: use ocean surface for slope calculation ==!
         zpice(:,:) = ssh_m(:,:)
      ENDIF

      DO jj = k_j1+1, k_jpj-1
         DO ji = fs_2, fs_jpim1

            zc1 = tmask(ji  ,jj  ,1) * ( rhosn * vt_s(ji  ,jj  ) + rhoic * vt_i(ji  ,jj  ) )
            zc2 = tmask(ji+1,jj  ,1) * ( rhosn * vt_s(ji+1,jj  ) + rhoic * vt_i(ji+1,jj  ) )
            zc3 = tmask(ji  ,jj+1,1) * ( rhosn * vt_s(ji  ,jj+1) + rhoic * vt_i(ji  ,jj+1) )

            zt11 = tmask(ji  ,jj,1) * e1t(ji  ,jj)
            zt12 = tmask(ji+1,jj,1) * e1t(ji+1,jj)
            zt21 = tmask(ji,jj  ,1) * e2t(ji,jj  )
            zt22 = tmask(ji,jj+1,1) * e2t(ji,jj+1)

            ! Leads area.
            zfrld1(ji,jj) = ( zt12 * ( 1.0 - at_i(ji,jj) ) + zt11 * ( 1.0 - at_i(ji+1,jj) ) ) / ( zt11 + zt12 + zepsi )
            zfrld2(ji,jj) = ( zt22 * ( 1.0 - at_i(ji,jj) ) + zt21 * ( 1.0 - at_i(ji,jj+1) ) ) / ( zt21 + zt22 + zepsi )

            ! Mass, coriolis coeff. and currents
            zmass1(ji,jj) = ( zt12 * zc1 + zt11 * zc2 ) / ( zt11 + zt12 + zepsi )
            zmass2(ji,jj) = ( zt22 * zc1 + zt21 * zc3 ) / ( zt21 + zt22 + zepsi )
            zcorl1(ji,jj) = zmass1(ji,jj) * ( e1t(ji+1,jj) * fcor(ji,jj) + e1t(ji,jj) * fcor(ji+1,jj) )   &
               &                          / ( e1t(ji,jj) + e1t(ji+1,jj) + zepsi )
            zcorl2(ji,jj) = zmass2(ji,jj) * ( e2t(ji,jj+1) * fcor(ji,jj) + e2t(ji,jj) * fcor(ji,jj+1) )   &
               &                          / ( e2t(ji,jj+1) + e2t(ji,jj) + zepsi )
            !
            ! Ocean has no slip boundary condition
            v_oce1(ji,jj)  = 0.5 * ( ( v_oce(ji  ,jj) + v_oce(ji  ,jj-1) ) * e1t(ji,jj)      &
               &                   + ( v_oce(ji+1,jj) + v_oce(ji+1,jj-1) ) * e1t(ji+1,jj) )  &
               &                   / ( e1t(ji+1,jj) + e1t(ji,jj) ) * umask(ji,jj,1)  

            u_oce2(ji,jj)  = 0.5 * ( ( u_oce(ji,jj  ) + u_oce(ji-1,jj  ) ) * e2t(ji,jj)      &
               &                   + ( u_oce(ji,jj+1) + u_oce(ji-1,jj+1) ) * e2t(ji,jj+1) )  &
               &                   / ( e2t(ji,jj+1) + e2t(ji,jj) ) * vmask(ji,jj,1)

            ! Wind stress at U,V-point
            ztagnx = ( 1. - zfrld1(ji,jj) ) * utau_ice(ji,jj)
            ztagny = ( 1. - zfrld2(ji,jj) ) * vtau_ice(ji,jj)

            ! Computation of the velocity field taking into account the ice internal interaction.
            ! Terms that are independent of the velocity field.

            ! SB On utilise maintenant le gradient de la pente de l'ocean
            ! include it later

            zdsshx =  ( zpice(ji+1,jj) - zpice(ji,jj) ) * r1_e1u(ji,jj)
            zdsshy =  ( zpice(ji,jj+1) - zpice(ji,jj) ) * r1_e2v(ji,jj)

            za1ct(ji,jj) = ztagnx - zmass1(ji,jj) * grav * zdsshx
            za2ct(ji,jj) = ztagny - zmass2(ji,jj) * grav * zdsshy

         END DO
      END DO

      !
      !------------------------------------------------------------------------------!
      ! 3) Solution of the momentum equation, iterative procedure
      !------------------------------------------------------------------------------!
      !
      ! Time step for subcycling
      dtevp  = rdt_ice / nn_nevp
#if defined key_lim3
      dtotel = dtevp / ( 2._wp * rn_relast * rdt_ice )
#else
      dtotel = dtevp / ( 2._wp * telast )
#endif
      z1_dtotel = 1._wp / ( 1._wp + dtotel )
      z1_dtevp  = 1._wp / dtevp
      !-ecc2: square of yield ellipse eccenticrity (reminder: must become a namelist parameter)
      ecc2 = rn_ecc * rn_ecc
      ecci = 1. / ecc2

      !-Initialise stress tensor 
      zs1 (:,:) = stress1_i (:,:) 
      zs2 (:,:) = stress2_i (:,:)
      zs12(:,:) = stress12_i(:,:)

      !                                               !----------------------!
      DO jter = 1 , nn_nevp                           !    loop over jter    !
         !                                            !----------------------!        
         DO jj = k_j1, k_jpj-1
            zu_ice(:,jj) = u_ice(:,jj)    ! velocity at previous time step
            zv_ice(:,jj) = v_ice(:,jj)
         END DO

         DO jj = k_j1+1, k_jpj-1
            DO ji = fs_2, fs_jpim1   !RB bug no vect opt due to zmask

               !  
               !- Divergence, tension and shear (Section a. Appendix B of Hunke & Dukowicz, 2002)
               !- divu_i(:,:), zdt(:,:): divergence and tension at centre of grid cells
               !- zds(:,:): shear on northeast corner of grid cells
               !
               !- IMPORTANT REMINDER: Dear Gurvan, note that, the way these terms are coded, 
               !                      there are many repeated calculations. 
               !                      Speed could be improved by regrouping terms. For
               !                      the moment, however, the stress is on clarity of coding to avoid
               !                      bugs (Martin, for Miguel).
               !
               !- ALSO: arrays zdt, zds and delta could 
               !  be removed in the future to minimise memory demand.
               !
               !- MORE NOTES: Note that we are calculating deformation rates and stresses on the corners of
               !              grid cells, exactly as in the B grid case. For simplicity, the indexation on
               !              the corners is the same as in the B grid.
               !
               !
               divu_i(ji,jj) = (  e2u(ji,jj) * u_ice(ji,jj) - e2u(ji-1,jj) * u_ice(ji-1,jj)   &
                  &             + e1v(ji,jj) * v_ice(ji,jj) - e1v(ji,jj-1) * v_ice(ji,jj-1)   &
                  &            ) * r1_e12t(ji,jj)

               zdt(ji,jj) = ( ( u_ice(ji,jj) * r1_e2u(ji,jj) - u_ice(ji-1,jj) * r1_e2u(ji-1,jj) ) * e2t(ji,jj) * e2t(ji,jj)   &
                  &         - ( v_ice(ji,jj) * r1_e1v(ji,jj) - v_ice(ji,jj-1) * r1_e1v(ji,jj-1) ) * e1t(ji,jj) * e1t(ji,jj)   &
                  &         ) * r1_e12t(ji,jj)

               !
               zds(ji,jj) = ( ( u_ice(ji,jj+1) * r1_e1u(ji,jj+1) - u_ice(ji,jj) * r1_e1u(ji,jj) ) * e1f(ji,jj) * e1f(ji,jj)   &
                  &         + ( v_ice(ji+1,jj) * r1_e2v(ji+1,jj) - v_ice(ji,jj) * r1_e2v(ji,jj) ) * e2f(ji,jj) * e2f(ji,jj)   &
                  &         ) * r1_e12f(ji,jj) * ( 2._wp - fmask(ji,jj,1) )   &
                  &         * zmask(ji,jj) * zmask(ji,jj+1) * zmask(ji+1,jj) * zmask(ji+1,jj+1)


               v_ice1(ji,jj)  = 0.5_wp * ( ( v_ice(ji  ,jj) + v_ice(ji  ,jj-1) ) * e1t(ji+1,jj)     &
                  &                      + ( v_ice(ji+1,jj) + v_ice(ji+1,jj-1) ) * e1t(ji  ,jj) )   &
                  &                      / ( e1t(ji+1,jj) + e1t(ji,jj) ) * umask(ji,jj,1) 

               u_ice2(ji,jj)  = 0.5_wp * ( ( u_ice(ji,jj  ) + u_ice(ji-1,jj  ) ) * e2t(ji,jj+1)     &
                  &                      + ( u_ice(ji,jj+1) + u_ice(ji-1,jj+1) ) * e2t(ji,jj  ) )   &
                  &                      / ( e2t(ji,jj+1) + e2t(ji,jj) ) * vmask(ji,jj,1)
            END DO
         END DO

         CALL lbc_lnk_multi( v_ice1, 'U', -1., u_ice2, 'V', -1. )      ! lateral boundary cond.
         
         DO jj = k_j1+1, k_jpj-1
            DO ji = fs_2, fs_jpim1

               !- Calculate Delta at centre of grid cells
               zdst          = ( e2u(ji,jj) * v_ice1(ji,jj) - e2u(ji-1,jj  ) * v_ice1(ji-1,jj  )   &
                  &            + e1v(ji,jj) * u_ice2(ji,jj) - e1v(ji  ,jj-1) * u_ice2(ji  ,jj-1)   &
                  &            ) * r1_e12t(ji,jj)

               delta          = SQRT( divu_i(ji,jj)**2 + ( zdt(ji,jj)**2 + zdst**2 ) * usecc2 )  
               delta_i(ji,jj) = delta + rn_creepl

               !- Calculate Delta on corners
               zddc  = (  ( v_ice1(ji,jj+1) * r1_e1u(ji,jj+1) - v_ice1(ji,jj) * r1_e1u(ji,jj) ) * e1f(ji,jj) * e1f(ji,jj)  &
                  &     + ( u_ice2(ji+1,jj) * r1_e2v(ji+1,jj) - u_ice2(ji,jj) * r1_e2v(ji,jj) ) * e2f(ji,jj) * e2f(ji,jj)  &
                  &    ) * r1_e12f(ji,jj)

               zdtc  = (- ( v_ice1(ji,jj+1) * r1_e1u(ji,jj+1) - v_ice1(ji,jj) * r1_e1u(ji,jj) ) * e1f(ji,jj) * e1f(ji,jj)  &
                  &     + ( u_ice2(ji+1,jj) * r1_e2v(ji+1,jj) - u_ice2(ji,jj) * r1_e2v(ji,jj) ) * e2f(ji,jj) * e2f(ji,jj)  &
                  &    ) * r1_e12f(ji,jj)

               zddc = SQRT( zddc**2 + ( zdtc**2 + zds(ji,jj)**2 ) * usecc2 ) + rn_creepl

               !-Calculate stress tensor components zs1 and zs2 at centre of grid cells (see section 3.5 of CICE user's guide).
               zs1(ji,jj)  = ( zs1 (ji,jj) + dtotel * ( divu_i(ji,jj) - delta ) / delta_i(ji,jj)   * zpresh(ji,jj)   &
                  &          ) * z1_dtotel
               zs2(ji,jj)  = ( zs2 (ji,jj) + dtotel *         ecci * zdt(ji,jj) / delta_i(ji,jj)   * zpresh(ji,jj)   &
                  &          ) * z1_dtotel
               !-Calculate stress tensor component zs12 at corners
               zs12(ji,jj) = ( zs12(ji,jj) + dtotel *         ecci * zds(ji,jj) / ( 2._wp * zddc ) * zpreshc(ji,jj)  &
                  &          ) * z1_dtotel 

            END DO
         END DO

         CALL lbc_lnk_multi( zs1 , 'T', 1., zs2, 'T', 1., zs12, 'F', 1. )
 
         ! Ice internal stresses (Appendix C of Hunke and Dukowicz, 2002)
         DO jj = k_j1+1, k_jpj-1
            DO ji = fs_2, fs_jpim1
               !- contribution of zs1, zs2 and zs12 to zf1
               zf1(ji,jj) = 0.5 * ( ( zs1(ji+1,jj) - zs1(ji,jj) ) * e2u(ji,jj)  &
                  &             + ( zs2(ji+1,jj) * e2t(ji+1,jj)**2 - zs2(ji,jj) * e2t(ji,jj)**2 ) * r1_e2u(ji,jj)          &
                  &             + 2.0 * ( zs12(ji,jj) * e1f(ji,jj)**2 - zs12(ji,jj-1) * e1f(ji,jj-1)**2 ) * r1_e1u(ji,jj)  &
                  &                ) * r1_e12u(ji,jj)
               ! contribution of zs1, zs2 and zs12 to zf2
               zf2(ji,jj) = 0.5 * ( ( zs1(ji,jj+1) - zs1(ji,jj) ) * e1v(ji,jj)  &
                  &             - ( zs2(ji,jj+1) * e1t(ji,jj+1)**2 - zs2(ji,jj) * e1t(ji,jj)**2 ) * r1_e1v(ji,jj)          &
                  &             + 2.0 * ( zs12(ji,jj) * e2f(ji,jj)**2 - zs12(ji-1,jj) * e2f(ji-1,jj)**2 ) * r1_e2v(ji,jj)  &
                  &               )  * r1_e12v(ji,jj)
            END DO
         END DO
         !
         ! Computation of ice velocity
         !
         ! Both the Coriolis term and the ice-ocean drag are solved semi-implicitly.
         !
         IF (MOD(jter,2).eq.0) THEN 

            DO jj = k_j1+1, k_jpj-1
               DO ji = fs_2, fs_jpim1
                  rswitch      = ( 1.0 - MAX( 0._wp, SIGN( 1._wp, -zmass1(ji,jj) ) ) ) * umask(ji,jj,1)
                  z0           = zmass1(ji,jj) * z1_dtevp

                  ! SB modif because ocean has no slip boundary condition
                  zv_ice1      = 0.5 * ( ( v_ice(ji  ,jj) + v_ice(ji  ,jj-1) ) * e1t(ji  ,jj)     &
                     &                 + ( v_ice(ji+1,jj) + v_ice(ji+1,jj-1) ) * e1t(ji+1,jj) )   &
                     &                 / ( e1t(ji+1,jj) + e1t(ji,jj) ) * umask(ji,jj,1)
                  za           = rhoco * SQRT( ( u_ice(ji,jj) - u_oce(ji,jj) )**2 +  &
                     &                         ( zv_ice1 - v_oce1(ji,jj) )**2 ) * ( 1.0 - zfrld1(ji,jj) )
                  zr           = z0 * u_ice(ji,jj) + zf1(ji,jj) + za1ct(ji,jj) + za * u_oce(ji,jj)
                  zcca         = z0 + za
                  zccb         = zcorl1(ji,jj)
                  u_ice(ji,jj) = ( zr + zccb * zv_ice1 ) / ( zcca + zepsi ) * rswitch 
               END DO
            END DO

            CALL lbc_lnk( u_ice(:,:), 'U', -1. )
#if defined key_agrif && defined key_lim2
            CALL agrif_rhg_lim2( jter, nn_nevp, 'U' )
#endif
#if defined key_bdy
         CALL bdy_ice_lim_dyn( 'U' )
#endif         

            DO jj = k_j1+1, k_jpj-1
               DO ji = fs_2, fs_jpim1

                  rswitch      = ( 1.0 - MAX( 0._wp, SIGN( 1._wp, -zmass2(ji,jj) ) ) ) * vmask(ji,jj,1)
                  z0           = zmass2(ji,jj) * z1_dtevp
                  ! SB modif because ocean has no slip boundary condition
                  zu_ice2      = 0.5 * ( ( u_ice(ji,jj  ) + u_ice(ji-1,jj  ) ) * e2t(ji,jj)       &
                     &                 + ( u_ice(ji,jj+1) + u_ice(ji-1,jj+1) ) * e2t(ji,jj+1) )   &
                     &                 / ( e2t(ji,jj+1) + e2t(ji,jj) ) * vmask(ji,jj,1)
                  za           = rhoco * SQRT( ( zu_ice2 - u_oce2(ji,jj) )**2 +  & 
                     &                         ( v_ice(ji,jj) - v_oce(ji,jj))**2 ) * ( 1.0 - zfrld2(ji,jj) )
                  zr           = z0 * v_ice(ji,jj) + zf2(ji,jj) + za2ct(ji,jj) + za * v_oce(ji,jj)
                  zcca         = z0 + za
                  zccb         = zcorl2(ji,jj)
                  v_ice(ji,jj) = ( zr - zccb * zu_ice2 ) / ( zcca + zepsi ) * rswitch
               END DO
            END DO

            CALL lbc_lnk( v_ice(:,:), 'V', -1. )
#if defined key_agrif && defined key_lim2
            CALL agrif_rhg_lim2( jter, nn_nevp, 'V' )
#endif
#if defined key_bdy
         CALL bdy_ice_lim_dyn( 'V' )
#endif         

         ELSE 
            DO jj = k_j1+1, k_jpj-1
               DO ji = fs_2, fs_jpim1
                  rswitch      = ( 1.0 - MAX( 0._wp, SIGN( 1._wp, -zmass2(ji,jj) ) ) ) * vmask(ji,jj,1)
                  z0           = zmass2(ji,jj) * z1_dtevp
                  ! SB modif because ocean has no slip boundary condition
                  zu_ice2      = 0.5 * ( ( u_ice(ji,jj  ) + u_ice(ji-1,jj  ) ) * e2t(ji,jj)       &
                     &                  +( u_ice(ji,jj+1) + u_ice(ji-1,jj+1) ) * e2t(ji,jj+1) )   &
                     &                 / ( e2t(ji,jj+1) + e2t(ji,jj) ) * vmask(ji,jj,1)   

                  za           = rhoco * SQRT( ( zu_ice2 - u_oce2(ji,jj) )**2 +  &
                     &                         ( v_ice(ji,jj) - v_oce(ji,jj) )**2 ) * ( 1.0 - zfrld2(ji,jj) )
                  zr           = z0 * v_ice(ji,jj) + zf2(ji,jj) + za2ct(ji,jj) + za * v_oce(ji,jj)
                  zcca         = z0 + za
                  zccb         = zcorl2(ji,jj)
                  v_ice(ji,jj) = ( zr - zccb * zu_ice2 ) / ( zcca + zepsi ) * rswitch
               END DO
            END DO

            CALL lbc_lnk( v_ice(:,:), 'V', -1. )
#if defined key_agrif && defined key_lim2
            CALL agrif_rhg_lim2( jter, nn_nevp, 'V' )
#endif
#if defined key_bdy
         CALL bdy_ice_lim_dyn( 'V' )
#endif         

            DO jj = k_j1+1, k_jpj-1
               DO ji = fs_2, fs_jpim1
                  rswitch      = ( 1.0 - MAX( 0._wp, SIGN( 1._wp, -zmass1(ji,jj) ) ) ) * umask(ji,jj,1)
                  z0           = zmass1(ji,jj) * z1_dtevp
                  zv_ice1      = 0.5 * ( ( v_ice(ji  ,jj) + v_ice(ji  ,jj-1) ) * e1t(ji,jj)       &
                     &                 + ( v_ice(ji+1,jj) + v_ice(ji+1,jj-1) ) * e1t(ji+1,jj) )   &
                     &                 / ( e1t(ji+1,jj) + e1t(ji,jj) ) * umask(ji,jj,1)

                  za           = rhoco * SQRT( ( u_ice(ji,jj) - u_oce(ji,jj) )**2 +  &
                     &                         ( zv_ice1 - v_oce1(ji,jj) )**2 ) * ( 1.0 - zfrld1(ji,jj) )
                  zr           = z0 * u_ice(ji,jj) + zf1(ji,jj) + za1ct(ji,jj) + za * u_oce(ji,jj)
                  zcca         = z0 + za
                  zccb         = zcorl1(ji,jj)
                  u_ice(ji,jj) = ( zr + zccb * zv_ice1 ) / ( zcca + zepsi ) * rswitch 
               END DO
            END DO

            CALL lbc_lnk( u_ice(:,:), 'U', -1. )
#if defined key_agrif && defined key_lim2
            CALL agrif_rhg_lim2( jter, nn_nevp, 'U' )
#endif
#if defined key_bdy
         CALL bdy_ice_lim_dyn( 'U' )
#endif         

         ENDIF
         
         IF(ln_ctl) THEN
            !---  Convergence test.
            DO jj = k_j1+1 , k_jpj-1
               zresr(:,jj) = MAX( ABS( u_ice(:,jj) - zu_ice(:,jj) ), ABS( v_ice(:,jj) - zv_ice(:,jj) ) )
            END DO
            zresm = MAXVAL( zresr( 1:jpi, k_j1+1:k_jpj-1 ) )
            IF( lk_mpp )   CALL mpp_max( zresm )   ! max over the global domain
         ENDIF

         !                                                ! ==================== !
      END DO                                              !  end loop over jter  !
      !                                                   ! ==================== !
      !
      !------------------------------------------------------------------------------!
      ! 4) Prevent ice velocities when the ice is thin
      !------------------------------------------------------------------------------!
      ! If the ice volume is below zvmin then ice velocity should equal the
      ! ocean velocity. This prevents high velocity when ice is thin
      DO jj = k_j1+1, k_jpj-1
         DO ji = fs_2, fs_jpim1
            IF ( vt_i(ji,jj) <= zvmin ) THEN
               u_ice(ji,jj) = u_oce(ji,jj)
               v_ice(ji,jj) = v_oce(ji,jj)
            ENDIF
         END DO
      END DO

      CALL lbc_lnk_multi( u_ice(:,:), 'U', -1., v_ice(:,:), 'V', -1. )

#if defined key_agrif && defined key_lim2
      CALL agrif_rhg_lim2( nn_nevp , nn_nevp, 'U' )
      CALL agrif_rhg_lim2( nn_nevp , nn_nevp, 'V' )
#endif
#if defined key_bdy
      CALL bdy_ice_lim_dyn( 'U' )
      CALL bdy_ice_lim_dyn( 'V' )
#endif         

      DO jj = k_j1+1, k_jpj-1 
         DO ji = fs_2, fs_jpim1
            IF ( vt_i(ji,jj) <= zvmin ) THEN
               v_ice1(ji,jj)  = 0.5_wp * ( ( v_ice(ji  ,jj) + v_ice(ji,  jj-1) ) * e1t(ji+1,jj)     &
                  &                      + ( v_ice(ji+1,jj) + v_ice(ji+1,jj-1) ) * e1t(ji  ,jj) )   &
                  &                      / ( e1t(ji+1,jj) + e1t(ji,jj) ) * umask(ji,jj,1)

               u_ice2(ji,jj)  = 0.5_wp * ( ( u_ice(ji,jj  ) + u_ice(ji-1,jj  ) ) * e2t(ji,jj+1)     &
                  &                      + ( u_ice(ji,jj+1) + u_ice(ji-1,jj+1) ) * e2t(ji,jj  ) )   &
                  &                      / ( e2t(ji,jj+1) + e2t(ji,jj) ) * vmask(ji,jj,1)
            ENDIF 
         END DO
      END DO

      CALL lbc_lnk_multi( u_ice2(:,:), 'V', -1., v_ice1(:,:), 'U', -1. )

      ! Recompute delta, shear and div, inputs for mechanical redistribution 
      DO jj = k_j1+1, k_jpj-1
         DO ji = fs_2, jpim1   !RB bug no vect opt due to zmask
            !- divu_i(:,:), zdt(:,:): divergence and tension at centre 
            !- zds(:,:): shear on northeast corner of grid cells
            IF ( vt_i(ji,jj) <= zvmin ) THEN

               divu_i(ji,jj) = (  e2u(ji,jj) * u_ice(ji,jj) - e2u(ji-1,jj  ) * u_ice(ji-1,jj  )   &
                  &             + e1v(ji,jj) * v_ice(ji,jj) - e1v(ji  ,jj-1) * v_ice(ji  ,jj-1)   &
                  &            ) * r1_e12t(ji,jj)

               zdt(ji,jj) = ( ( u_ice(ji,jj) * r1_e2u(ji,jj) - u_ice(ji-1,jj) * r1_e2u(ji-1,jj) ) * e2t(ji,jj) * e2t(ji,jj)  &
                  &          -( v_ice(ji,jj) * r1_e1v(ji,jj) - v_ice(ji,jj-1) * r1_e1v(ji,jj-1) ) * e1t(ji,jj) * e1t(ji,jj)  &
                  &         ) * r1_e12t(ji,jj)
               !
               ! SB modif because ocean has no slip boundary condition 
               zds(ji,jj) = ( ( u_ice(ji,jj+1) * r1_e1u(ji,jj+1) - u_ice(ji,jj) * r1_e1u(ji,jj) ) * e1f(ji,jj) * e1f(ji,jj)  &
                  &          +( v_ice(ji+1,jj) * r1_e2v(ji+1,jj) - v_ice(ji,jj) * r1_e2v(ji,jj) ) * e2f(ji,jj) * e2f(ji,jj)  &
                  &         ) * r1_e12f(ji,jj) * ( 2.0 - fmask(ji,jj,1) )                                     &
                  &         * zmask(ji,jj) * zmask(ji,jj+1) * zmask(ji+1,jj) * zmask(ji+1,jj+1)

               zdst = ( e2u(ji,jj) * v_ice1(ji,jj) - e2u(ji-1,jj  ) * v_ice1(ji-1,jj  )    &
                  &   + e1v(ji,jj) * u_ice2(ji,jj) - e1v(ji  ,jj-1) * u_ice2(ji  ,jj-1) ) * r1_e12t(ji,jj)

               delta = SQRT( divu_i(ji,jj)**2 + ( zdt(ji,jj)**2 + zdst**2 ) * usecc2 )  
               delta_i(ji,jj) = delta + rn_creepl
            
            ENDIF
         END DO
      END DO
      !
      !------------------------------------------------------------------------------!
      ! 5) Store stress tensor and its invariants
      !------------------------------------------------------------------------------!
      ! * Invariants of the stress tensor are required for limitd_me
      !   (accelerates convergence and improves stability)
      DO jj = k_j1+1, k_jpj-1
         DO ji = fs_2, fs_jpim1
            zdst           = (  e2u(ji,jj) * v_ice1(ji,jj) - e2u( ji-1, jj   ) * v_ice1(ji-1,jj)  &   
               &              + e1v(ji,jj) * u_ice2(ji,jj) - e1v( ji  , jj-1 ) * u_ice2(ji,jj-1) ) * r1_e12t(ji,jj) 
            shear_i(ji,jj) = SQRT( zdt(ji,jj) * zdt(ji,jj) + zdst * zdst )
         END DO
      END DO

      ! Lateral boundary condition
      CALL lbc_lnk_multi( divu_i (:,:), 'T', 1., delta_i(:,:), 'T', 1.,  shear_i(:,:), 'T', 1. )

      ! * Store the stress tensor for the next time step
      stress1_i (:,:) = zs1 (:,:)
      stress2_i (:,:) = zs2 (:,:)
      stress12_i(:,:) = zs12(:,:)

      !
      !------------------------------------------------------------------------------!
      ! 6) Control prints of residual and charge ellipse
      !------------------------------------------------------------------------------!
      !
      ! print the residual for convergence
      IF(ln_ctl) THEN
         WRITE(charout,FMT="('lim_rhg  : res =',D23.16, ' iter =',I4)") zresm, jter
         CALL prt_ctl_info(charout)
         CALL prt_ctl(tab2d_1=u_ice, clinfo1=' lim_rhg  : u_ice :', tab2d_2=v_ice, clinfo2=' v_ice :')
      ENDIF

      ! print charge ellipse
      ! This can be desactivated once the user is sure that the stress state
      ! lie on the charge ellipse. See Bouillon et al. 08 for more details
      IF(ln_ctl) THEN
         CALL prt_ctl_info('lim_rhg  : numit  :',ivar1=numit)
         CALL prt_ctl_info('lim_rhg  : nwrite :',ivar1=nwrite)
         CALL prt_ctl_info('lim_rhg  : MOD    :',ivar1=MOD(numit,nwrite))
         IF( MOD(numit,nwrite) .EQ. 0 ) THEN
            WRITE(charout,FMT="('lim_rhg  :', I4, I6, I1, I1, A10)") 1000, numit, 0, 0, ' ch. ell. '
            CALL prt_ctl_info(charout)
            DO jj = k_j1+1, k_jpj-1
               DO ji = 2, jpim1
                  IF (zpresh(ji,jj) > 1.0) THEN
                     sigma1 = ( zs1(ji,jj) + (zs2(ji,jj)**2 + 4*zs12(ji,jj)**2 )**0.5 ) / ( 2*zpresh(ji,jj) ) 
                     sigma2 = ( zs1(ji,jj) - (zs2(ji,jj)**2 + 4*zs12(ji,jj)**2 )**0.5 ) / ( 2*zpresh(ji,jj) )
                     WRITE(charout,FMT="('lim_rhg  :', I4, I4, D23.16, D23.16, D23.16, D23.16, A10)")
                     CALL prt_ctl_info(charout)
                  ENDIF
               END DO
            END DO
            WRITE(charout,FMT="('lim_rhg  :', I4, I6, I1, I1, A10)") 2000, numit, 0, 0, ' ch. ell. '
            CALL prt_ctl_info(charout)
         ENDIF
      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj, zpresh, zfrld1, zmass1, zcorl1, za1ct , zpreshc, zfrld2, zmass2, zcorl2, za2ct )
      CALL wrk_dealloc( jpi,jpj, u_oce2, u_ice2, v_oce1 , v_ice1 , zmask               )
      CALL wrk_dealloc( jpi,jpj, zf1   , zu_ice, zf2   , zv_ice , zdt    , zds  )
      CALL wrk_dealloc( jpi,jpj, zdt   , zds   , zs1   , zs2   , zs12   , zresr , zpice                 )

   END SUBROUTINE lim_rhg

#else
   !!----------------------------------------------------------------------
   !!   Default option          Dummy module           NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_rhg( k1 , k2 )         ! Dummy routine
      WRITE(*,*) 'lim_rhg: You should not have seen this print! error?', k1, k2
   END SUBROUTINE lim_rhg
#endif

   !!==============================================================================
END MODULE limrhg
