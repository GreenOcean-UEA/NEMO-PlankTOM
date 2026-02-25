MODULE trcnam_planktom
#if defined key_planktom
!!======================================================================
!!                         ***  MODULE trcnam_planktom  ***
!! TOP : PlankTOM ecosystem model 
!!======================================================================
!! History : original  : 99-10 (M.A. Foujols, M. Levy) passive tracer
!!           addition  : 00-01 (L. Bopp) hamocc3,p3zd
!!           modification : 02 (E. Buitenhuis) dgom
!!           see bgcbio
!!----------------------------------------------------------------------
!!    'key_planktom'      :             use the PlankTOM ecosystem model
!!    'key_trc_n2o'       :                        calculates N2O fluxes
!!    'key_trc_ch4'       :                        calculates CH4 fluxes
!!    'key_trc_dms'       :            includes the DMS cycle and fluxes
!!    'key_trc_foram'     :      includes zooplankton calcifiers (foram) 
!!    'key_trc_diaadd'    :       save tracer diagnostic files 2D and 3D
!!----------------------------------------------------------------------
!! ** Purpose : READs options for the PlankTOM namelist
!!              Variables that affect phytoplankton have been put in arrays
!!
!! ** Action  : Read namelist input
!!              &natdet           : particulate aggregation parameters
!!              &natriv           : river elemental input parameters
!!              &natfer           : iron model parameters
!!              &natpre           : preference for food parameters (including for bacteria)
!!              &natzoo           : biological parameters related to heterotrophs
!!              &natphy           : biological parameters related to phyto
!!              &natpft           : biological parameters related to all PFTs
!!              &natquo           : iron quota model parameters
!!              &natlit           : light-iron model parameters
!!              &natzca           : calcite and aragonite parameters
!!              &natdms           : DMS parameters
!!              &natn2o           : N2O parameters
!!              &natch4           : CH4 parameters (with key_trc_ch4 only)
!!
!! References : PlankTOM12 manual, Buitenhuis et al. 2023
!!                                 https://zenodo.org/records/8388158 
!!======================================================================
!!
!! trc_nam_planktom       : PlankTOM model namelist read
!!----------------------------------------------------------------------
!
! parameters and commons
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE sms_planktom    ! sms trends
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_planktom   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_planktom.F90 2567 2011-01-25 09:36:27Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_planktom
!!----------------------------------------------------------------------
! 
! local variables
      INTEGER jl, jm, jn
!
! Initializations
! ---------------
!
       namelist/natdet/rn_ag1poc,rn_ag2poc,rn_ag3poc,rn_ag4poc,rn_ag5doc,rn_ag6doc,rn_ag7doc,          &
     &                 rn_readsi,rn_remdsi,rn_retdsi,                                                  &
     &                 rn_sisgoc,rn_siegoc,rn_snkhoc,rn_snkpoc,rn_snkmax
       namelist/natriv/rn_rivdic,rn_rivdoc,rn_rivfer,rn_rivnit,rn_rivpo4,rn_rivpoc,rn_rivsil
       namelist/natfer/rn_fersol,rn_silsol,rn_sedfer,rn_scofer,rn_scmfer,rn_ligdpn,rn_ligdps,rn_ligfrn,rn_liglat,rn_ligfrs
       namelist/natpre/rn_bmsdet,rn_gbadoc,rn_gbagoc,rn_gbagon,rn_gbahoc,rn_gbapoc,rn_bmspft,rn_prfzoo
       namelist/natzoo/nn_gonmin,nn_gonmax,nn_morsqu,nn_sizzoo,rn_ggebac,rn_ggtbac,rn_ggezoo,          &
     &                 rn_sigzoo,                                                                      &
     &                 rn_unazoo,                                                                      &
     &                 rn_grabac,rn_gramin,rn_grazoo,                                                  &
     &                 rn_kmobac,rn_grkzoo,                                                            &
     &                 rn_kmfbac,rn_kmpbac,rn_denno3,                                                  &
     &                 rn_morcru,rn_motcru,                                                            &
#if ! defined key_trc_foram
     &                 rn_morgel,rn_motgel,                                                            &
#endif
     &                 rn_rembac,rn_resbac,rn_reszoo,                                                  &
     &                 rn_retbac,rn_retzoo,                                                            &
     &                 rn_icecru
       namelist/natphy/rn_bsidia,rn_coccal,rn_discal,rn_ferbsi,rn_kmsbsi,rn_silbsi,                    &
     &                 rn_lyoco3,rn_lyscal,rn_sildia,rn_munfix,                                        &
     &                 rn_resphy,                                                                      &
     &                 rn_docphy,rn_domphy,                                                            &
     &                 rn_kmnphy,rn_kmpphy,rn_mumpft
       namelist/natpft/rn_mutpft,rn_mudpft
       namelist/natquo/rn_kmfphy,rn_rhfphy,rn_qmaphy,rn_qmiphy,rn_qopphy
       namelist/natlit/rn_ekwgrn,rn_ekwred,rn_alpphy,rn_kgrphy,rn_krdphy,rn_thmphy
       namelist/natzca/rn_disara,rn_disfor,rn_forcal,rn_pteara,rn_lysara
#    if defined key_trc_dms
       namelist/natdms/rn_dmsyld,rn_etomax,rn_rdddms,rn_xcldmd,rn_xpodms,                              &
     &                 rn_xkdms,rn_xkdmd,rn_assdms,rn_rphdmd,rn_xpldmd
#    endif
#    if defined key_trc_n2o
       namelist/natn2o/nn_deun2s,rn_aoun2s,rn_betn2s,rn_decn2s,rn_omxn2s
#    endif
#    if defined key_trc_ch4
       namelist/natch4/rn_proch1,rn_proch2,rn_proch3,rn_proch4,rn_proch5,                              &
     &                 rn_botch1,rn_botch2,rn_botch3,rn_botch4,rn_botch5,                              &
     &                 rn_conch1,rn_conch2,rn_conch3,rn_conch4,rn_conch5,                              &
     &                 rn_decch4
#    endif
#  if defined key_trc_diaadd
       namelist/natadd/ctrc3d,ctrc3l,ctrc2d,ctrc2l,ctrc3u,ctrc2u,nn_writedia
#  endif
       namelist/nat_biodts/nn_biodts
!
      CHARACTER (len=39) ::   clname
      INTEGER :: numnatpl
!
      IF(lwp) THEN
        WRITE(numout,*) ' '
        WRITE(numout,*) ' Read namelist for PlankTOM model'
        WRITE(numout,*) ' ********************************'
        WRITE(numout,*) ' '
        CALL FLUSH(numout)
      ENDIF
#  if defined key_trc_diaadd && ! key_iomput
      STOP 'key_trc_diaadd without key_iomput is no longer supported'
#  endif
!    
! open the namelist file 
      numnatpl=80
      clname='namelist.trc.sms'
      OPEN( numnatpl, FILE= clname, FORM='formatted', STATUS = 'old')
!
! Read the namelist 
! -----------------
!
      READ(numnatpl,natdet)
      READ(numnatpl,natriv)
      READ(numnatpl,natfer)
      READ(numnatpl,natpre)
      READ(numnatpl,natzoo)
      READ(numnatpl,natphy)
      READ(numnatpl,natpft)
      READ(numnatpl,natquo)
      READ(numnatpl,natlit)
      READ(numnatpl,natzca)
#    if defined key_trc_dms
      READ(numnatpl,natdms)
#    endif
#    if defined key_trc_n2o
      READ(numnatpl,natn2o)
#    endif
#    if defined key_trc_ch4
      READ(numnatpl,natch4)
#    endif
      READ(numnatpl,nat_biodts)
      REWIND(numnatpl) 
!
! Printout the namelist in the ocean output file 
! ----------------------------------------------
!
! &natdet           : particulate aggregation parameters
! ------------------------------------------------------
!
      IF(lwp) THEN
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natdet'
        WRITE(numout,*) ' '
        WRITE(numout,*) ' aggregation term 1                         (rn_ag1poc) =', rn_ag1poc
        WRITE(numout,*) ' aggregation term 2                         (rn_ag2poc) =', rn_ag2poc
        WRITE(numout,*) ' aggregation term 3                         (rn_ag3poc) =', rn_ag3poc
        WRITE(numout,*) ' aggregation term 4                         (rn_ag4poc) =', rn_ag4poc
        WRITE(numout,*) ' aggregation term 5                         (rn_ag5doc) =', rn_ag5doc
        WRITE(numout,*) ' aggregation term 6                         (rn_ag6doc) =', rn_ag6doc
        WRITE(numout,*) ' aggregation term 7                         (rn_ag7doc) =', rn_ag7doc
        WRITE(numout,*) ' maximum DSi degradation                    (rn_readsi) =', rn_readsi
        WRITE(numout,*) ' DSi degradation at infinite T              (rn_remdsi) =', rn_remdsi
        WRITE(numout,*) ' T dependence DSi degradation               (rn_retdsi) =', rn_retdsi
        WRITE(numout,*) ' POC sinking speed                          (rn_snkpoc) =', rn_snkpoc
        WRITE(numout,*) ' Big particles sinking multiplication       (rn_siegoc) =', rn_siegoc
        WRITE(numout,*) ' Big particles sinking division             (rn_sisgoc) =', rn_sisgoc
        WRITE(numout,*) ' Huge particles sinking speed               (rn_snkhoc) =', rn_snkhoc
        WRITE(numout,*) ' Maximum sink speed                         (rn_snkmax) =', rn_snkmax
!
! &natriv           : river elemental input parameters
! ----------------------------------------------------
!
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natriv'
        WRITE(numout,*) ' '
        WRITE(numout,*) ' river conversion of DIC (rn_rivdic) =', rn_rivdic
        WRITE(numout,*) ' river conversion of DOC (rn_rivdoc) =', rn_rivdoc
        WRITE(numout,*) ' river conversion of FER (rn_rivfer) =', rn_rivfer
        WRITE(numout,*) ' river conversion of NIT (rn_rivnit) =', rn_rivnit
        WRITE(numout,*) ' river conversion of PO4 (rn_rivpo4) =', rn_rivpo4
        WRITE(numout,*) ' river conversion of POC (rn_rivpoc) =', rn_rivpoc
        WRITE(numout,*) ' river conversion of SIL (rn_rivsil) =', rn_rivsil
!
! &natfer           : iron model parameters
! -----------------------------------------
!
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natfer'
        WRITE(numout,*) ' '
        WRITE(numout,*) ' solubility of iron in dust                             (rn_fersol) =', rn_fersol
        WRITE(numout,*) ' solubility of silica in dust                           (rn_silsol) =', rn_silsol
        WRITE(numout,*) ' coastal release of iron                                (rn_sedfer) =', rn_sedfer
        WRITE(numout,*) ' scavenging rate of iron                                (rn_scofer) =', rn_scofer
        WRITE(numout,*) ' minimum scavenging rate of iron                        (rn_scmfer) =', rn_scmfer
        WRITE(numout,*) ' ligand concentration of iron in the north              (rn_ligfrn) =', rn_ligfrn
        WRITE(numout,*) ' north boundary of SO to differentiate ligand conc.     (rn_liglat) =', rn_liglat
        WRITE(numout,*) ' depth iron conc. restored to ligand conc. in the north (rn_ligdpn) =', rn_ligdpn
        WRITE(numout,*) ' depth iron conc. restored to ligand conc. in the SO    (rn_ligdps) =', rn_ligdps
        WRITE(numout,*) ' ligand concentration of iron in the Southern Ocean     (rn_ligfrs) =', rn_ligfrs
!
! &natphy           : biological parameters related to phyto
! &natpft           : biological parameters related to all PFTs (phytoplankton ! listed only)
! &natzca           : calcite and aragonite parameters
! -------------------------------------------------------------------------------------------
!
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natphy, natpft, natzca'
        WRITE(numout,*) ' '
        WRITE(numout,*) ' maximum silification:photosynthesis diatoms (rn_bsidia) =', rn_bsidia
        WRITE(numout,*) ' calcification:photosynthesis coccolithoph.  (rn_coccal) =', rn_coccal
        WRITE(numout,*) ' calcite dissolution with coc loss           (rn_discal) =', rn_discal
        WRITE(numout,*) ' aragonite dissolution with pte loss         (rn_disara) =', rn_disara
#    if defined key_trc_foram
        WRITE(numout,*) ' calcite dissolution with for loss           (rn_disfor) =', rn_disfor
        WRITE(numout,*) ' calcification:grazing foraminifers          (rn_forcal) =', rn_forcal
#    endif
        WRITE(numout,*) ' calcification:grazing pteropods             (rn_pteara) =', rn_pteara
        WRITE(numout,*) ' subsaturated max aragonite dissolution rate (rn_lysara) =', rn_lysara
        WRITE(numout,*) ' subsaturated max calcite dissolution rate   (rn_lyscal) =', rn_lyscal
!
! loop over each phytoplankton
        DO jl = jpdia, jpdia+jppft-1
          WRITE(numout,*) ' '
          WRITE(numout,*) ctrcnm(jl)
          WRITE(numout,*) ' nitrogen half saturation concentrat.      (rn_kmnphy) =', rn_kmnphy(jl)
          WRITE(numout,*) ' phosphate half saturation concentration   (rn_kmpphy) =', rn_kmpphy(jl)
          WRITE(numout,*) ' phyto respiration as a fraction of growth (rn_resphy) =', rn_resphy(jl)
          WRITE(numout,*) ' excretion ratio                           (rn_docphy) =', rn_docphy(jl)
          WRITE(numout,*) ' maximum DOC excretion ratio               (rn_domphy) =', rn_domphy(jl)
          WRITE(numout,*) ' optimal growth rate                       (rn_mumpft) =', rn_mumpft(jl)
          WRITE(numout,*) ' temperature dependance growth rate        (rn_mutpft) =', rn_mutpft(jl)
          WRITE(numout,*) ' temperature width growth rate mudpft      (rn_mudpft) =', rn_mudpft(jl)
        END DO
!
        WRITE(numout,*) ' fraction of growth rate during N2fix      (rn_munfix) =', rn_munfix
        WRITE(numout,*) ' exponent for calcite dissolution rate     (rn_lyoco3) =', rn_lyoco3
!
! &natzoo           : biological parameters related to heterotrophs
! &natpft           : biological parameters related to all PFTs (zooplankton ! listed only)
! -----------------------------------------------------------------------------------------
!
! output global total amount of carbon in each food source
        WRITE(numout,*) ' '
        DO jl = jppoc,jpdia+jppft-1
          WRITE(numout,101) ' total amount of carbon in food source ',ctrcnm(jl),'  (rn_bmspft) =  ',rn_bmspft(jl)
          101 FORMAT(A40, A3, A17, F5.3)
        END DO
!
! loop over each zooplankton 
        DO jm = 1, jpzft
          WRITE(numout,*) ' '
          WRITE(numout,*) ctrcnm(jpbac+jm)
          WRITE(numout,*) ' maximal grazing rate                       (rn_grazoo) =', rn_grazoo(jm)
          WRITE(numout,*) ' half saturation grazing                    (rn_grkzoo) =', rn_grkzoo(jm)
          WRITE(numout,*) ' temperature dependance grazingrate         (rn_mutpft) =', rn_mutpft(jpbac+jm)
          WRITE(numout,*) ' temperature dependance grazingrate         (rn_mudpft) =', rn_mudpft(jpbac+jm)
          WRITE(numout,*) ' gross growth efficiency                    (rn_ggezoo) =', rn_ggezoo(jm)
          WRITE(numout,*) ' unassimilated fraction/fecal pell.         (rn_unazoo) =', rn_unazoo(jm)
          IF ( nn_sizzoo(jm).LT.0 .OR. nn_sizzoo(jm).GT.2)                      &
     &      CALL ctl_stop('E R R O R : nn_sizzoo element has inconsistent values')
          IF ( nn_sizzoo(jm).EQ.0 .AND. (nn_gonmin.LE.jm .OR. nn_gonmax.LE.jm)) &
     &      CALL ctl_stop('E R R O R : large zoo cannot produce sPOC')
          IF ( nn_sizzoo(jm).EQ.2 .AND. (nn_gonmin.GE.jm .OR. nn_gonmax.GE.jm)) &
     &      CALL ctl_stop('E R R O R : small zoo cannot produce HOC')
          WRITE(numout,*) ' losses go to POC=0, GOC=1 or HOC=2         (nn_sizzoo) =', nn_sizzoo(jm)
          WRITE(numout,*) ' fraction of respiration as PO4             (rn_sigzoo) =', rn_sigzoo(jm)
          WRITE(numout,*) ' respiration rate of zooplankton            (rn_reszoo) =', rn_reszoo(jm)
          WRITE(numout,*) ' temp. dep. zooplankton respiration         (rn_retzoo) =', rn_retzoo(jm)
!
! loop over each food source
          DO jl = jppoc,jpdia+jppft-1
            WRITE(numout,*) ' preference for ',ctrcnm(jl),'        (rn_prfzoo) = ',rn_prfzoo(jm,jl)
          END DO
        END DO
!
        WRITE(numout,*) ' '
        IF ( nn_gonmin.LT.1 .OR. nn_gonmin.GT.jpzft .OR.  nn_gonmax.LT.1 .OR. nn_gonmax.GT.jpzft ) &
     &      CALL ctl_stop('E R R O R : nn_gonmin or nn_gonmax has illegal value')
        WRITE(numout,*) ' zPFTs that feed into GON smallest, largest (nn_gonmin,nn_gonmax) =', nn_gonmin,nn_gonmax
        WRITE(numout,*) ' minimum food for grazing                             (rn_gramin) =', rn_gramin
!
        IF (nn_morsqu .EQ. 1) THEN
          WRITE(numout,*) ' square mortality GEL and CRU                         (nn_morsqu) =',nn_morsqu
        ELSE
          WRITE(numout,*) ' mortality of GEL and CRU uses        SUM(PFT biomass; nn_morsqu) =',nn_morsqu
        ENDIF
!
        WRITE(numout,*) ' crustacean macrozoo mortality rate                   (rn_morcru) =', rn_morcru
        WRITE(numout,*) ' temp. dep. crustacean macrozoo mortality rate        (rn_motcru) =', rn_motcru
#  if ! defined key_trc_foram
        WRITE(numout,*) ' gelatinouszooplankton mortality rate                 (rn_morgel) =', rn_morgel
        WRITE(numout,*) ' temp. dep. gelatinouszooplankton mortality rate      (rn_motgel) =', rn_motgel
#  endif
        WRITE(numout,*) ' crustacean macrozoo growth enhancement under ice     (rn_icecru) =', rn_icecru
!
! bacteria-related parameters
        WRITE(numout,*) ' '
        WRITE(numout,*) ' maximum growth rate bacteria                           (rn_grabac) =', rn_grabac
        WRITE(numout,*) ' temperature depend. growth bacteria                    (rn_mutpft) =', rn_mutpft(jpbac)
        WRITE(numout,*) ' temperature depend. growth bacteria                    (rn_mudpft) =', rn_mudpft(jpbac)
        WRITE(numout,*) ' Fe half-sat. of bacteria                               (rn_kmfbac) =', rn_kmfbac
        WRITE(numout,*) ' DOC half-sat. of bacteria                              (rn_kmobac) =', rn_kmobac
        WRITE(numout,*) ' PO4 half-sat. of bacteria                              (rn_kmpbac) =', rn_kmpbac
        WRITE(numout,*) ' global total amount of carbon in DOC, POC, GOC, HOC    (rn_bmsdet) =', rn_bmsdet(1),rn_bmsdet(2),rn_bmsdet(3),rn_bmsdet(4)
        WRITE(numout,*) ' bacteria relative preference for DOC                   (rn_gbadoc) =', rn_gbadoc
        WRITE(numout,*) ' bacteria relative preference for GOC                   (rn_gbagoc) =', rn_gbagoc
        WRITE(numout,*) ' bacteria relative preference for GON                   (rn_gbagon) =', rn_gbagon
        WRITE(numout,*) ' bacteria relative preference for POC                   (rn_gbapoc) =', rn_gbapoc
        WRITE(numout,*) ' bacteria relative preference for HOC                   (rn_gbahoc) =', rn_gbahoc
        WRITE(numout,*) ' bacteria growth efficiency and T depend.    (rn_ggebac, rn_ggtbac) =', rn_ggebac,rn_ggtbac
        WRITE(numout,*) ' bacteria respiration rate and T depend.     (rn_resbac, rn_retbac) =', rn_resbac,rn_retbac
        WRITE(numout,*) ' minimum concentration for bacteria in remineralisation (rn_rembac) =', rn_rembac
        WRITE(numout,*) ' oxygen dependence denitrification                      (rn_denno3) =', rn_denno3
!
        WRITE(numout,*) ' '
        WRITE(numout,*) ' silicate half saturation constant diatoms  (rn_sildia) =', rn_sildia
        WRITE(numout,*) ' iron limited multiplier of Si:C            (rn_ferbsi) =', rn_ferbsi
        WRITE(numout,*) ' silicate saturated multiplier of Si:C      (rn_silbsi) =', rn_silbsi
        WRITE(numout,*) ' silicate half saturation constant for Si:C (rn_kmsbsi) =', rn_kmsbsi
!
! &natquo           : iron quota model parameters
! -----------------------------------------------
!
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natquo'
        WRITE(numout,*) ' '
!
! loop over each phytoplankton
        DO jl = jpdia, jpdia+jppft-1
          WRITE(numout,*) ctrcnm(jl)
          WRITE(numout,*)' iron half saturation conc.              (rn_kmfphy) =', rn_kmfphy(jl)
          WRITE(numout,*)' threshold for nutrient uptake for Fe    (rn_rhfphy) =', rn_rhfphy(jl)
          WRITE(numout,*)' maximum quota for Fe                    (rn_qmaphy) =', rn_qmaphy(jl)
          WRITE(numout,*)' minimum quota for Fe                    (rn_qmiphy) =', rn_qmiphy(jl)
          WRITE(numout,*)' optimal quota for Fe                    (rn_qopphy) =', rn_qopphy(jl)
        END DO
!
! &natlit           : light-iron model parameters
! -----------------------------------------------
!
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natlit'
        WRITE(numout,*) ' '
        WRITE(numout,*) ' red light absorption coeff. of water   (rn_ekwred) =', rn_ekwred
        WRITE(numout,*) ' green light absorption coeff. of water (rn_ekwgrn) =', rn_ekwgrn
        DO jl = jpdia, jpdia+jppft-1
          WRITE(numout,*) ctrcnm(jl)
          WRITE(numout,*) ' initial slope PI curve                 (rn_alpphy) =', rn_alpphy(jl)
          WRITE(numout,*) ' maximum chl:C ratio                    (rn_thmphy) =', rn_thmphy(jl)
          WRITE(numout,*) ' light absorption in the red            (rn_krdphy) =', rn_krdphy(jl)
          WRITE(numout,*) ' light absorption in the blue-green     (rn_kgrphy) =', rn_kgrphy(jl)
        END DO
#    if defined key_trc_dms
!
! &natdms           : DMS parameters
! ----------------------------------
!
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natdms'
        WRITE(numout,*) ' '
        WRITE(numout,*) ' Ratio of DMS/DMSP released during grazing     (rn_rdddms) =', rn_rdddms
        WRITE(numout,*) ' DMSP-lyase DMSPd cleavage rate                (rn_xcldmd) =', rn_xcldmd
        WRITE(numout,*) ' microbial yield for DMS production            (rn_dmsyld) =', rn_dmsyld
        WRITE(numout,*) ' half saturation constant for bacteria on DMSPd (rn_xkdmd) =', rn_xkdmd
        WRITE(numout,*) ' half saturation constant for bacteria on DMS   (rn_xkdms) =', rn_xkdms
! 
! loop over each zooplankton
        DO jm = 1, jpzft
          WRITE(numout,*) ' % of DMSP in fecal pellets lost           (rn_assdms) =', rn_assdms(jm)
        END DO
!
        WRITE(numout,*) ' maximum surface insolation                    (rn_etomax) =', rn_etomax
!
! loop over each phytoplankton
        DO jl = jpdia, jpdia+jppft-1
          WRITE(numout,*) ctrcnm(jl)
          WRITE(numout,*) ' mean DMSP/C cell ratio                     (rn_rphdmd) =', rn_rphdmd(jl)
          WRITE(numout,*) ' DMSP leakage coefficient                   (rn_xpldmd) =', rn_rphdmd(jl)
        END DO
#    endif
!
! &natn2o           : N2O parameters
! ----------------------------------
!
#    if defined key_trc_n2o
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natn2o'
        WRITE(numout,*) ' '
        WRITE(numout,*) 'First depth with N2S production                     (nn_deun2s) =', nn_deun2s
        WRITE(numout,*) 'N2O/AOU production ratio                            (rn_aoun2s) =', rn_aoun2s
        WRITE(numout,*) 'hypoxic N2O/AOU production ratio                    (rn_betn2s) =', rn_betn2s
        WRITE(numout,*) 'Exponent of decreasing N2O production with O2       (rn_decn2s) =', rn_decn2s
        WRITE(numout,*) 'O2 concentration below which prod. stops increasing (rn_omxn2s) =', rn_omxn2s
#    endif
!
! &natch4           : CH4 parameters (with key_trc_ch4 only)
! ----------------------------------------------------------
!
#    if defined key_trc_ch4
        WRITE(numout,*) ' '
        WRITE(numout,*) 'natch4'
        WRITE(numout,*) ' '
        WRITE(numout,*) 'CH4 production ratios (rn_proch?) =', rn_proch1,rn_proch2,rn_proch3,rn_proch4,rn_proch5
#    endif
        CALL FLUSH(numout)
      ENDIF

   END SUBROUTINE trc_nam_planktom

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                   No PlankTOM bio-model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_planktom                      ! Empty routine
   END  SUBROUTINE  trc_nam_planktom
#endif  

   !!======================================================================
END MODULE trcnam_planktom
