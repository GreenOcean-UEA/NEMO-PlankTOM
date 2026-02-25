MODULE limwri
   !!======================================================================
   !!                     ***  MODULE  limwri  ***
   !!         Ice diagnostics :  write ice output files
   !!======================================================================
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3'                                      LIM3 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_wri      : write of the diagnostics variables in ouput file 
   !!   lim_wri_state : write for initial state or/and abandon
   !!----------------------------------------------------------------------
   USE ioipsl
   USE dianam          ! build name of file (routine)
   USE phycst
   USE dom_oce
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE dom_ice
   USE ice
   USE limvar
   USE in_out_manager
   USE lbclnk
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! work arrays
   USE iom
   USE timing          ! Timing
   USE lib_fortran     ! Fortran utilities

   IMPLICIT NONE
   PRIVATE

   PUBLIC lim_wri        ! routine called by lim_step.F90
   PUBLIC lim_wri_state  ! called by dia_wri_state 

   !!----------------------------------------------------------------------
   !! NEMO/LIM3 4.0 , UCL - NEMO Consortium (2011)
   !! $Id: limwri.F90 5517 2015-06-30 13:09:58Z clem $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_dimgout
# include "limwri_dimg.h90"
#else

   SUBROUTINE lim_wri( kindic )
      !!-------------------------------------------------------------------
      !!  This routine computes the average of some variables and write it
      !!  on the ouput files.
      !!  ATTENTION cette routine n'est valable que si le pas de temps est
      !!  egale a une fraction entiere de 1 jours.
      !!  Diff 1-D 3-D : suppress common also included in etat
      !!                 suppress cmoymo 11-18
      !!  modif : 03/06/98
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kindic   ! if kindic < 0 there has been an error somewhere
      !
      INTEGER  ::  ji, jj, jk, jl  ! dummy loop indices
      REAL(wp) ::  z1_365
      REAL(wp) ::  ztmp
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zoi, zei, zt_i, zt_s
      REAL(wp), POINTER, DIMENSION(:,:)   ::  z2d, z2da, z2db, zswi    ! 2D workspace
      !!-------------------------------------------------------------------

      IF( nn_timing == 1 )  CALL timing_start('limwri')

      CALL wrk_alloc( jpi, jpj, jpl, zoi, zei, zt_i, zt_s )
      CALL wrk_alloc( jpi, jpj     , z2d, z2da, z2db, zswi )

      !-----------------------------
      ! Mean category values
      !-----------------------------
      z1_365 = 1._wp / 365._wp

      CALL lim_var_icetm      ! mean sea ice temperature

      CALL lim_var_bv         ! brine volume

      DO jj = 1, jpj          ! presence indicator of ice
         DO ji = 1, jpi
            zswi(ji,jj)  = MAX( 0._wp , SIGN( 1._wp , at_i(ji,jj) - epsi06 ) )
         END DO
      END DO
      !
      !
      !                                             
      IF ( iom_use( "icethic_cea" ) ) THEN                       ! mean ice thickness
         DO jj = 1, jpj 
            DO ji = 1, jpi
               z2d(ji,jj)  = vt_i(ji,jj) / MAX( at_i(ji,jj), epsi06 ) * zswi(ji,jj)
            END DO
         END DO
         CALL iom_put( "icethic_cea"  , z2d              )
      ENDIF

      IF ( iom_use( "snowthic_cea" ) ) THEN                      ! snow thickness = mean snow thickness over the cell 
         DO jj = 1, jpj                                            
            DO ji = 1, jpi
               z2d(ji,jj)  = vt_s(ji,jj) / MAX( at_i(ji,jj), epsi06 ) * zswi(ji,jj)
            END DO
         END DO
         CALL iom_put( "snowthic_cea" , z2d              )       
      ENDIF
      !
      IF ( iom_use( "uice_ipa" ) .OR. iom_use( "vice_ipa" ) .OR. iom_use( "icevel" ) ) THEN 
         DO jj = 2 , jpjm1
            DO ji = 2 , jpim1
               z2da(ji,jj)  = (  u_ice(ji,jj) * umask(ji,jj,1) + u_ice(ji-1,jj) * umask(ji-1,jj,1) ) * 0.5_wp
               z2db(ji,jj)  = (  v_ice(ji,jj) * vmask(ji,jj,1) + v_ice(ji,jj-1) * vmask(ji,jj-1,1) ) * 0.5_wp
           END DO
         END DO
         CALL lbc_lnk( z2da, 'T', -1. )
         CALL lbc_lnk( z2db, 'T', -1. )
         CALL iom_put( "uice_ipa"     , z2da             )       ! ice velocity u component
         CALL iom_put( "vice_ipa"     , z2db             )       ! ice velocity v component
         DO jj = 1, jpj                                 
            DO ji = 1, jpi
               z2d(ji,jj)  = SQRT( z2da(ji,jj) * z2da(ji,jj) + z2db(ji,jj) * z2db(ji,jj) ) 
            END DO
         END DO
         CALL iom_put( "icevel"       , z2d              )       ! ice velocity module
      ENDIF
      !
      IF ( iom_use( "miceage" ) ) THEN 
         z2d(:,:) = 0.e0
         DO jl = 1, jpl
            DO jj = 1, jpj
               DO ji = 1, jpi
                  rswitch    = MAX( 0._wp , SIGN( 1._wp , at_i(ji,jj) - 0.1 ) )
                  z2d(ji,jj) = z2d(ji,jj) + rswitch * oa_i(ji,jj,jl) / MAX( at_i(ji,jj), 0.1 )
               END DO
            END DO
         END DO
         CALL iom_put( "miceage"     , z2d * z1_365      )        ! mean ice age
      ENDIF

      IF ( iom_use( "micet" ) ) THEN 
         DO jj = 1, jpj
            DO ji = 1, jpi
               z2d(ji,jj) = ( tm_i(ji,jj) - rt0 ) * zswi(ji,jj)
            END DO
         END DO
         CALL iom_put( "micet"       , z2d               )        ! mean ice temperature
      ENDIF
      !
      IF ( iom_use( "icest" ) ) THEN 
         z2d(:,:) = 0.e0
         DO jl = 1, jpl
            DO jj = 1, jpj
               DO ji = 1, jpi
                  z2d(ji,jj) = z2d(ji,jj) + zswi(ji,jj) * ( t_su(ji,jj,jl) - rt0 ) * a_i(ji,jj,jl) / MAX( at_i(ji,jj) , epsi06 )
               END DO
            END DO
         END DO
         CALL iom_put( "icest"       , z2d              )        ! ice surface temperature
      ENDIF

      IF ( iom_use( "icecolf" ) ) THEN 
         DO jj = 1, jpj
            DO ji = 1, jpi
               rswitch  = MAX( 0._wp , SIGN( 1._wp , at_i(ji,jj) ) )
               z2d(ji,jj) = hicol(ji,jj) * rswitch
            END DO
         END DO
         CALL iom_put( "icecolf"     , z2d              )        ! frazil ice collection thickness
      ENDIF

      CALL iom_put( "isst"        , sst_m               )        ! sea surface temperature
      CALL iom_put( "isss"        , sss_m               )        ! sea surface salinity
      CALL iom_put( "iceconc"     , at_i                )        ! ice concentration
      CALL iom_put( "icevolu"     , vt_i                )        ! ice volume = mean ice thickness over the cell
      CALL iom_put( "icehc"       , et_i                )        ! ice total heat content
      CALL iom_put( "isnowhc"     , et_s                )        ! snow total heat content
      CALL iom_put( "ibrinv"      , bv_i * 100._wp      )        ! brine volume
      CALL iom_put( "utau_ice"    , utau_ice            )        ! wind stress over ice along i-axis at I-point
      CALL iom_put( "vtau_ice"    , vtau_ice            )        ! wind stress over ice along j-axis at I-point
      CALL iom_put( "snowpre"     , sprecip * 86400.    )        ! snow precipitation 
      CALL iom_put( "micesalt"    , smt_i               )        ! mean ice salinity

      CALL iom_put( "icestr"      , strength * 0.001    )        ! ice strength
      CALL iom_put( "idive"       , divu_i * 1.0e8      )        ! divergence
      CALL iom_put( "ishear"      , shear_i * 1.0e8     )        ! shear
      CALL iom_put( "snowvol"     , vt_s                )        ! snow volume
      
      CALL iom_put( "icetrp"      , diag_trp_vi * rday  )        ! ice volume transport
      CALL iom_put( "snwtrp"      , diag_trp_vs * rday  )        ! snw volume transport
      CALL iom_put( "saltrp"      , diag_trp_smv * rday * rhoic ) ! salt content transport
      CALL iom_put( "deitrp"      , diag_trp_ei         )        ! advected ice enthalpy (W/m2)
      CALL iom_put( "destrp"      , diag_trp_es         )        ! advected snw enthalpy (W/m2)

      CALL iom_put( "sfxbog"      , sfx_bog * rday      )        ! salt flux from brines
      CALL iom_put( "sfxbom"      , sfx_bom * rday      )        ! salt flux from brines
      CALL iom_put( "sfxsum"      , sfx_sum * rday      )        ! salt flux from brines
      CALL iom_put( "sfxsni"      , sfx_sni * rday      )        ! salt flux from brines
      CALL iom_put( "sfxopw"      , sfx_opw * rday      )        ! salt flux from brines
      CALL iom_put( "sfxdyn"      , sfx_dyn * rday      )        ! salt flux from ridging rafting
      CALL iom_put( "sfxres"      , sfx_res * rday      )        ! salt flux from limupdate (resultant)
      CALL iom_put( "sfxbri"      , sfx_bri * rday      )        ! salt flux from brines
      CALL iom_put( "sfx"         , sfx     * rday      )        ! total salt flux

      ztmp = rday / rhoic
      CALL iom_put( "vfxres"     , wfx_res * ztmp       )        ! daily prod./melting due to limupdate 
      CALL iom_put( "vfxopw"     , wfx_opw * ztmp       )        ! daily lateral thermodynamic ice production
      CALL iom_put( "vfxsni"     , wfx_sni * ztmp       )        ! daily snowice ice production
      CALL iom_put( "vfxbog"     , wfx_bog * ztmp       )        ! daily bottom thermodynamic ice production
      CALL iom_put( "vfxdyn"     , wfx_dyn * ztmp       )        ! daily dynamic ice production (rid/raft)
      CALL iom_put( "vfxsum"     , wfx_sum * ztmp       )        ! surface melt 
      CALL iom_put( "vfxbom"     , wfx_bom * ztmp       )        ! bottom melt 
      CALL iom_put( "vfxice"     , wfx_ice * ztmp       )        ! total ice growth/melt 
      CALL iom_put( "vfxsnw"     , wfx_snw * ztmp       )        ! total snw growth/melt 
      CALL iom_put( "vfxsub"     , wfx_sub * ztmp       )        ! sublimation (snow) 
      CALL iom_put( "vfxspr"     , wfx_spr * ztmp       )        ! precip (snow)
      
      CALL iom_put( "afxtot"     , afx_tot * rday       )        ! concentration tendency (total)
      CALL iom_put( "afxdyn"     , afx_dyn * rday       )        ! concentration tendency (dynamics)
      CALL iom_put( "afxthd"     , afx_thd * rday       )        ! concentration tendency (thermo)

      CALL iom_put ('hfxthd'     , hfx_thd(:,:)         )   !  
      CALL iom_put ('hfxdyn'     , hfx_dyn(:,:)         )   !  
      CALL iom_put ('hfxres'     , hfx_res(:,:)         )   !  
      CALL iom_put ('hfxout'     , hfx_out(:,:)         )   !  
      CALL iom_put ('hfxin'      , hfx_in(:,:)          )   !  
      CALL iom_put ('hfxsnw'     , hfx_snw(:,:)         )   !  
      CALL iom_put ('hfxsub'     , hfx_sub(:,:)         )   !  
      CALL iom_put ('hfxerr'     , hfx_err(:,:)         )   !  
      CALL iom_put ('hfxerr_rem' , hfx_err_rem(:,:)     )   !  
      
      CALL iom_put ('hfxsum'     , hfx_sum(:,:)         )   !  
      CALL iom_put ('hfxbom'     , hfx_bom(:,:)         )   !  
      CALL iom_put ('hfxbog'     , hfx_bog(:,:)         )   !  
      CALL iom_put ('hfxdif'     , hfx_dif(:,:)         )   !  
      CALL iom_put ('hfxopw'     , hfx_opw(:,:)         )   !  
      CALL iom_put ('hfxtur'     , fhtur(:,:) * SUM(a_i_b(:,:,:), dim=3) ) ! turbulent heat flux at ice base 
      CALL iom_put ('hfxdhc'     , diag_heat(:,:)       )   ! Heat content variation in snow and ice 
      CALL iom_put ('hfxspr'     , hfx_spr(:,:)         )   ! Heat content of snow precip 
      
      !--------------------------------
      ! Output values for each category
      !--------------------------------
      CALL iom_put( "iceconc_cat"      , a_i         )        ! area for categories
      CALL iom_put( "icethic_cat"      , ht_i        )        ! thickness for categories
      CALL iom_put( "snowthic_cat"     , ht_s        )        ! snow depth for categories
      CALL iom_put( "salinity_cat"     , sm_i        )        ! salinity for categories

      ! ice temperature
      IF ( iom_use( "icetemp_cat" ) ) THEN 
         zt_i(:,:,:) = SUM( t_i(:,:,:,:), dim=3 ) * r1_nlay_i
         CALL iom_put( "icetemp_cat"   , zt_i - rt0  )
      ENDIF
      
      ! snow temperature
      IF ( iom_use( "snwtemp_cat" ) ) THEN 
         zt_s(:,:,:) = SUM( t_s(:,:,:,:), dim=3 ) * r1_nlay_s
         CALL iom_put( "snwtemp_cat"   , zt_s - rt0  )
      ENDIF

      ! Compute ice age
      IF ( iom_use( "iceage_cat" ) ) THEN 
         DO jl = 1, jpl 
            DO jj = 1, jpj
               DO ji = 1, jpi
                  rswitch = MAX( 0._wp , SIGN( 1._wp , at_i(ji,jj) - 0.1 ) )
                  rswitch = rswitch * MAX( 0._wp , SIGN( 1._wp , a_i(ji,jj,jl) - 0.1 ) )
                  zoi(ji,jj,jl) = oa_i(ji,jj,jl)  / MAX( a_i(ji,jj,jl) , 0.1 ) * rswitch
               END DO
            END DO
         END DO
         CALL iom_put( "iceage_cat"   , zoi * z1_365 )        ! ice age for categories
      ENDIF

      ! Compute brine volume
      IF ( iom_use( "brinevol_cat" ) ) THEN 
         zei(:,:,:) = 0._wp
         DO jl = 1, jpl 
            DO jk = 1, nlay_i
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     rswitch = MAX( 0._wp , SIGN( 1._wp , a_i(ji,jj,jl) - epsi06 ) )
                     zei(ji,jj,jl) = zei(ji,jj,jl) + 100.0 *  &
                        ( - tmut * s_i(ji,jj,jk,jl) / MIN( ( t_i(ji,jj,jk,jl) - rt0 ), - epsi06 ) ) * &
                        rswitch * r1_nlay_i
                  END DO
               END DO
            END DO
         END DO
         CALL iom_put( "brinevol_cat"     , zei      )        ! brine volume for categories
      ENDIF

      !     !  Create an output files (output.lim.abort.nc) if S < 0 or u > 20 m/s
      !     IF( kindic < 0 )   CALL lim_wri_state( 'output.abort' )
      !     not yet implemented
      
      CALL wrk_dealloc( jpi, jpj, jpl, zoi, zei, zt_i, zt_s )
      CALL wrk_dealloc( jpi, jpj     , z2d, zswi, z2da, z2db )

      IF( nn_timing == 1 )  CALL timing_stop('limwri')
      
   END SUBROUTINE lim_wri
#endif

 
   SUBROUTINE lim_wri_state( kt, kid, kh_i )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE lim_wri_state  ***
      !!        
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains 
      !!      the instantaneous ice state and forcing fields for ice model
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! History :
      !!   4.1  !  2013-06  (C. Rousset)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt               ! ocean time-step index)
      INTEGER, INTENT( in ) ::   kid , kh_i       
      !!----------------------------------------------------------------------

      CALL histdef( kid, "iicethic", "Ice thickness"           , "m"      ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iiceconc", "Ice concentration"       , "%"      ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicetemp", "Ice temperature"         , "C"      ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicevelu", "i-Ice speed (I-point)"   , "m/s"    ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicevelv", "j-Ice speed (I-point)"   , "m/s"    ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "iicestru", "i-Wind stress over ice (I-pt)", "Pa",   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicestrv", "j-Wind stress over ice (I-pt)", "Pa",   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "iicesflx", "Solar flux over ocean"     , "w/m2" ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "iicenflx", "Non-solar flux over ocean" , "w/m2" ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "isnowpre", "Snow precipitation"      , "kg/m2/s",   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "iicesali", "Ice salinity"            , "PSU"    ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "iicevolu", "Ice volume"              , "m"      ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "iicedive", "Ice divergence"          , "10-8s-1",   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt ) 
      CALL histdef( kid, "iicebopr", "Ice bottom production"   , "m/s"    ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicedypr", "Ice dynamic production"  , "m/s"    ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicelapr", "Ice open water prod"     , "m/s"    ,   &
      &       jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicesipr", "Snow ice production "    , "m/s"    ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicerepr", "Ice prod from limupdate" , "m/s"    ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicebome", "Ice bottom melt"         , "m/s"    ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iicesume", "Ice surface melt"        , "m/s"    ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iisfxdyn", "Salt flux from dynmics"  , ""       ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )
      CALL histdef( kid, "iisfxres", "Salt flux from limupdate", ""       ,   &
      &      jpi, jpj, kh_i, 1, 1, 1, -99, 32, "inst(x)", rdt, rdt )

      CALL histend( kid, snc4set )   ! end of the file definition

      CALL histwrite( kid, "iicethic", kt, icethi        , jpi*jpj, (/1/) )    
      CALL histwrite( kid, "iiceconc", kt, at_i          , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicetemp", kt, tm_i - rt0    , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicevelu", kt, u_ice          , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicevelv", kt, v_ice          , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicestru", kt, utau_ice       , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicestrv", kt, vtau_ice       , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicesflx", kt, qsr , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicenflx", kt, qns , jpi*jpj, (/1/) )
      CALL histwrite( kid, "isnowpre", kt, sprecip        , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicesali", kt, smt_i          , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicevolu", kt, vt_i           , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicedive", kt, divu_i*1.0e8   , jpi*jpj, (/1/) )

      CALL histwrite( kid, "iicebopr", kt, wfx_bog        , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicedypr", kt, wfx_dyn        , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicelapr", kt, wfx_opw        , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicesipr", kt, wfx_sni        , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicerepr", kt, wfx_res        , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicebome", kt, wfx_bom        , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iicesume", kt, wfx_sum        , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iisfxdyn", kt, sfx_dyn        , jpi*jpj, (/1/) )
      CALL histwrite( kid, "iisfxres", kt, sfx_res        , jpi*jpj, (/1/) )

      ! Close the file
      ! -----------------
      !CALL histclo( kid )

    END SUBROUTINE lim_wri_state

#else
   !!----------------------------------------------------------------------
   !!   Default option :         Empty module          NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_wri          ! Empty routine
   END SUBROUTINE lim_wri
#endif

   !!======================================================================
END MODULE limwri
