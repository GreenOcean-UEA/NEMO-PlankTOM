MODULE par_planktom
!!======================================================================
!!                         ***  MODULE par_planktom  ***
!! TOP : PlankTOM ecosystem model 
!!======================================================================
!! History : PISCES       ! 2001 (O. Aumont) Original code 
!!           TOP 1.0      ! LOCEAN-IPSL (2005)
!!           see bgcbio
!!----------------------------------------------------------------------
!!    'key_planktom'      :             use the PlankTOM ecosystem model
!!    'key_trc_n2o'       :                        calculates N2O fluxes
!!    'key_trc_ch4'       :                        calculates CH4 fluxes
!!    'key_trc_dms'       :            includes the DMS cycle and fluxes
!!    'key_trc_cfc11'     :            includes the CFC11 fluxes
!!    'key_trc_foram'     :      includes zooplankton calcifiers (foram) 
!!    'key_trc_piic'      :    additional tracer with pre-industrial DIC 
!!----------------------------------------------------------------------
!! ** Purpose : set the PlankTOM ecosystem parameters
!! ** Action  : 
!!
!! References : PlankTOM12 manual, Buitenhuis et al. 2023
!!                                 https://zenodo.org/records/8388158 
!!
!! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
!!======================================================================
! 
#if defined key_planktom
   IMPLICIT NONE
   PUBLIC
!
!! ----------------------------------------------------------------------
!!    Assign an integer to name individual tracers.
!!    For planktom, the tracers should be in the order:
!!    chemical tracers, nutrients, detritus, zooplankton,
!!    phytoplankton (C,Fe,Chl,Si), other tracers and isotopes
!!    For documentation of variables see:
!!    https://doi.org/10.5281/zenodo.8388158
!! ----------------------------------------------------------------------
!
      LOGICAL, PUBLIC, PARAMETER :: lk_planktom     = .TRUE.  !: PlankTOM flag 
!
! fix the number of PFTs
      INTEGER, PUBLIC, PARAMETER :: jpzft=5            ! number of zooplankton PFTs
      INTEGER, PUBLIC, PARAMETER :: jppft=6            ! number of phytoplankton PFTs
      INTEGER, PUBLIC, PARAMETER :: jpfmx=jpzft*(3+jpzft+jppft) ! number of preference parameters expected
!
! set the carbonate tracers
      INTEGER, PUBLIC, PARAMETER :: jptal=1,jpoxy=2,jpdic=3  ! alkalinity, oxygen, DIC
!
! set the nutrient tracers
#  if defined key_trc_piic
!
! set one additional carbonate tracer, then the nutrient tracers
      INTEGER, PUBLIC, PARAMETER :: jppiic=4,jpdin=5
#  else
      INTEGER, PUBLIC, PARAMETER :: jpdin=4
#  endif
      INTEGER, PUBLIC, PARAMETER :: jpsil=jpdin+1,jppo4=jpdin+2
      INTEGER, PUBLIC, PARAMETER :: jpfer=jpdin+3
!
! set the dissolved and particulate organic tracers (detritus)
      INTEGER, PUBLIC, PARAMETER :: jpdoc=jpdin+4,jpdsi=jpdin+5,jpcal=jpdin+6
      INTEGER, PUBLIC, PARAMETER :: jpara=jpdin+7,jpgon=jpdin+8,jpsfe=jpdin+9
      INTEGER, PUBLIC, PARAMETER :: jpbfe=jpsfe+1,jpufe=jpsfe+2,jppoc=jpsfe+3,jpgoc=jpsfe+4,jphoc=jpsfe+5
!
! set the bacteria PFT
      INTEGER, PUBLIC, PARAMETER :: jpbac=jpsfe+6
!
! set the zooplankton PFTs
      INTEGER, PUBLIC, PARAMETER :: jppro=jpsfe+7
#    if defined key_trc_foram
      INTEGER, PUBLIC, PARAMETER :: jpfor=jpsfe+8,jppte=jpsfe+9,jpmes=jpsfe+10,jpcru=jpsfe+11
#    else
      INTEGER, PUBLIC, PARAMETER :: jppte=jpsfe+8,jpmes=jpsfe+9,jpgel=jpsfe+10,jpcru=jpsfe+11
#    endif
!
! set the phytoplankton PFTs
      INTEGER, PUBLIC, PARAMETER :: jpdia=jpcru+1 ,jpmix=jpcru+2 ,jpcoc=jpcru+3 ,jppic=jpcru+4 ,jppha=jpcru+5 ,jpfix=jpcru+6
!
! set the iron content of phytoplankton PFTs
      INTEGER, PUBLIC, PARAMETER :: jpdfe=jpcru+7 ,jpnfe=jpcru+8 ,jpcfe=jpcru+9 ,jppfe=jpcru+10,jphfe=jpcru+11,jpffe=jpcru+12
!
! set the chlorophyll content of phytoplankton PFTs
      INTEGER, PUBLIC, PARAMETER :: jpdch=jpcru+13,jpnch=jpcru+14,jpcch=jpcru+15,jppch=jpcru+16,jphch=jpcru+17,jpfch=jpcru+18
!
! set the silicon content of diatom PFT
      INTEGER, PUBLIC, PARAMETER :: jpbsi=jpcru+19
#  if defined key_trc_cfc11
!
! set the CFC11 tracer 
      INTEGER, PUBLIC, PARAMETER :: jpc11=jpcru+20
#  else
      INTEGER, PUBLIC, PARAMETER :: jpc11=jpcru+19
#  endif
#  if defined key_trc_dms
!
! set the DMS tracers
      INTEGER, PUBLIC, PARAMETER :: jpdms=jpc11+1,jpdmd=jpc11+2
#  else
      INTEGER, PUBLIC, PARAMETER :: jpdmd=jpc11
#  endif
#  if defined key_trc_n2o
!
! set the N2O tracers
      INTEGER, PUBLIC, PARAMETER :: jpn2s=jpdmd+1
#  else
      INTEGER, PUBLIC, PARAMETER :: jpn2s=jpdmd
#  endif
#  if defined key_trc_ch4
!
! set the CH4 tracers
      INTEGER, PUBLIC, PARAMETER :: jpch1 =jpn2s+1 ,jpch2 =jpn2s+2 ,jpch3 =jpn2s+3 ,jpch4 =jpn2s+4 ,jpch5 =jpn2s+5
      INTEGER, PUBLIC, PARAMETER :: jpch6 =jpn2s+6 ,jpch7 =jpn2s+7 ,jpch8 =jpn2s+8 ,jpch9 =jpn2s+9 ,jpch10=jpn2s+10
      INTEGER, PUBLIC, PARAMETER :: jpch11=jpn2s+11,jpch12=jpn2s+12,jpch13=jpn2s+13,jpch14=jpn2s+14,jpch15=jpn2s+15
      INTEGER, PUBLIC, PARAMETER :: jpch16=jpn2s+16,jpch17=jpn2s+17,jpch18=jpn2s+18,jpch19=jpn2s+19,jpch20=jpn2s+20
      INTEGER, PUBLIC, PARAMETER :: jpch21=jpn2s+21,jpch22=jpn2s+22,jpch23=jpn2s+23,jpch24=jpn2s+24,jpch25=jpn2s+25
      INTEGER, PUBLIC, PARAMETER :: jp_planktom=jpn2s+25
#  else
      INTEGER, PUBLIC, PARAMETER :: jp_planktom=jpn2s
#  endif
#else
      LOGICAL, PUBLIC, PARAMETER :: lk_planktom = .FALSE.
      INTEGER, PUBLIC, PARAMETER :: jp_planktom = 0
#endif
END MODULE par_planktom
