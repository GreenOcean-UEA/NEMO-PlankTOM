MODULE sms_planktom
!!======================================================================
!!                         ***  MODULE sms_planktom  ***
!! TOP : PlankTOM ecosystem model 
!!======================================================================
!! History : TOP 1.0      ! LOCEAN-IPSL (2005)
!!           see bgcbio
!!----------------------------------------------------------------------
!!    'key_planktom'      :             use the PlankTOM ecosystem model
!!    'key_trc_n2o'       :                        calculates N2O fluxes
!!    'key_trc_ch4'       :                        calculates CH4 fluxes
!!    'key_trc_dms'       :            includes the DMS cycle and fluxes
!!    'key_c14b'          :       calculate 14C fluxes and concentration
!!    'key_trc_piic'      :    additional tracer with pre-industrial DIC 
!!----------------------------------------------------------------------
!! ** Purpose : Allocate variables for the sources minus sinks (sms)
!! ** Action  : 
!!
!! References : PlankTOM12 manual, Buitenhuis et al. 2023
!!                                 https://zenodo.org/records/8388158 
!!
!! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
!!======================================================================
#if defined key_planktom && key_top
!
! parameters and commons
   USE par_oce
   USE par_trc
   USE par_planktom
   IMPLICIT NONE
   PUBLIC
!!
!!----------------------------------------------------------------------
!! Variable for chemistry of the CO2 cycle
!! ---------------------------------------------------------------------
!!
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: akb3, ak13, ak23
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: aksara, aksp, co3, hi
#  if defined key_trc_piic
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: pihi
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: pih2co3
#  endif
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: borat, hco3, discarb
!
      REAL, PUBLIC, SAVE :: atcox
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) ::  h2co3
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: akw3
!!
!!      VARIABLE           TYPE    PURPOSE.
!!      --------           ----    --------
!!
!!      *THIRD*            REAL    1./3. HALF PRECISION
!!      *THOUSI*           REAL    0.001 HALF PRECICION
!!      *SMICR*            REAL    1E-6 HALF PRECISION
!!      *SALCHL*           REAL    CONVERSION FACTOR TO CALCULATE CHLORID
!!                                 CONCENTRATION
!!                                 S(O/OO)=1.80655*CL(O/OO)
!!                                 KALLE/DIETRICH , P. 60, OR
!!                                 WOOSTER ET AL., 1969.
!!      *TEMZER*           REAL    ZERO DEG C EXPRESSED AS ABS. TEMPERATURE
!!                                 (DEG KELVIN)
!!      *PERC*             REAL    0.01 HALF PRECISION
!!      *AKCC1*            REAL    COEFFICIENT FOR EMPIRICALLY DETERMINED
!!                                 APPARENT SOLUBILITY PRODUCT K'SP OF
!!                                 CALCITE (INGLE, 1800, EQ. 6)
!!                                 (CF. BROECKER ET AL., 1982)
!!      *AKCC2*            REAL    COEFFICIENT FOR EMPIRICALLY DETERMINED
!!                                 APPARENT SOLUBILITY PRODUCT K'SP OF
!!                                 CALCITE (INGLE, 1800, EQ. 6)
!!                                 (CF. BROECKER ET AL., 1982)
!!      *AKCC3*            REAL    COEFFICIENT FOR EMPIRICALLY DETERMINED
!!                                 APPARENT SOLUBILITY PRODUCT K'SP OF
!!                                 CALCITE (INGLE, 1800, EQ. 6)
!!                                 (CF. BROECKER ET AL., 1982)
!!      *AKCC4*            REAL    COEFFICIENT FOR EMPIRICALLY DETERMINED
!!                                 APPARENT SOLUBILITY PRODUCT K'SP OF
!!                                 CALCITE (INGLE, 1800, EQ. 6)
!!                                 (CF. BROECKER ET AL., 1982)
!!      *ARAFRA*           REAL    FRACTION OF ARAGONITE IN BIOGENIC CACO3
!!                                 PARTICLES (E.G. 0.3 FOR 30 PERCENT)
!!      *CALFRA*           REAL    FRACTION OF CALCITE IN BIOGENIC CACO3
!!                                 PARTICLES (E.G. 0.3 FOR 30 PERCENT)
!!      *ARACAL*           REAL    FACTOR TO CONVERT APP. CALCITE SOLUBILITY
!!                                 PRODUCT (0 DBAR) TO THE APP. SOLUBILITY
!!                                 PRODUCT OF ARAGONITE (BERNER, 1976;
!!                                 CF. BROECKER ET AL., 1982)
!!      *DEVK1*            REAL    COEFFICIENT FOR SEAWATER PRESSURE CORRECTION
!!                                 OF 1. DISSOCIATION CONSTANT OF CARBONIC
!!                                 ACID AFTER CULBERSON AND PYTKOWICZ, 1968
!!                                 (CF. BROECKER ET AL., 1982)
!!      *DEVK2*            REAL    COEFFICIENT FOR SEAWATER PRESSURE CORRECTION
!!                                 OF 2. DISSOCIATION CONSTANT OF CARBONIC
!!                                 ACID AFTER CULBERSON AND PYTKOWICZ, 1968
!!                                 (CF. BROECKER ET AL., 1982)
!!      *DEVKB*            REAL    COEFFICIENT FOR SEAWATER PRESSURE CORRECTION
!!                                 OF 1. DISSOCIATION CONSTANT OF BORIC
!!                                 ACID AFTER CULBERSON AND PYTKOWICZ, 1968
!!                                 (CF. BROECKER ET AL., 1982)
!!      *DEVK1T*           REAL    COEFFICIENT FOR SEAWATER PRESSURE CORRECTION
!!                                 OF FIRST DISSOCIATION CONSTANT OF CARBONIC
!!                                 ACID AFTER CULBERSON AND PYTKOWICZ, 1968
!!                                 (CF. BROECKER ET AL., 1982)
!!      *DEVK2T*           REAL    COEFFICIENT FOR SEAWATER PRESSURE CORRECTION
!!                                 OF SECOND DISSOCIATION CONSTANT OF CARBONIC
!!                                 ACID AFTER CULBERSON AND PYTKOWICZ, 1968
!!                                 (CF. BROECKER ET AL., 1982)
!!      *DEVKBT*           REAL    COEFFICIENT FOR SEAWATER PRESSURE CORRECTION
!!                                 OF DISSOCIATION CONSTANT OF BORIC
!!                                 ACID AFTER CULBERSON AND PYTKOWICZ, 1968
!!                                 (CF. BROECKER ET AL., 1982)
!!      *DEVKS*            REAL    COEFFICIENT FOR PRESSURE CORRECTION OF
!!                                 SOLUBILITY PRODUCT OF CALCITE OR ARAGONITE
!!                                 AFTER EDMOND AND GIESKES (1970), P. 1285
!!                                 (REFERENCE TO CULBERSON AND PYTKOWICZ, 1968,
!!                                 AS DONE IN BROECKER ET AL., 1982, IS
!!                                 NOT CORRECT)
!!      *DEVKST*           REAL    COEFFICIENT FOR PRESSURE CORRECTION OF
!!                                 SOLUBILITY PRODUCT OF CALCITE OR ARAGONITE
!!                                 AFTER EDMOND AND GIESKES (1970), P. 1285,
!!                                 IN TERM WITH TEMPERATURE
!!                                 (REFERENCE TO CULBERSON AND PYTKOWICZ, 1968,
!!                                 AS DONE IN BROECKER ET AL., 1982, IS
!!                                 NOT CORRECT)
!!      *RGAS*             REAL    UNIVERSAL GAS CONSTANT (BOLTZMANN'S CONSTANT
!!                                 TIMES AVOGADRO'S CONSTANT =
!!                                 1.3804E-16*6.023*10E+23=8.3143E+7 ERG/GRD*MOL=
!!                                 83.143E+6 ERG/GRD*MOL=8.3143 J/K*MOL)
!!                                 MULTIPLIED WITH 10 (TO ACCOUNT FOR
!!                                 CHANGE FROM BAR TO DBAR)
!!                                 (CF. EDMOND AND GIESKES, 1970. P. 1285,
!!                                 BROECKER ET AL., 1982, P. 79)
!!      *BOR1*             REAL    TOTAL BORON CONTENT IN G/KG AT CL=19 O/OO
!!                                 (S=35)
!!                                 (CF. RILEY AND SKIRROW, VOL. 1, P. 648)
!!      *BOR2*             REAL    INVERSE OF ATOMIC WEIGHT OF BORON FOR
!!                                 CONVERTING SPECIFIC TOTAL BORATE IN
!!                                 CONCENTRATIONS
!!      *OXYCO*            REAL    INVERS OF NORMAL MOLAL VOLUME OF AN
!!                                 IDEAL GAS [CM**-3]
!!      *C00*              REAL    VOLUMETRIC SOLUBILITY CONSTANT A1 FOR
!!                                 THE SOLUBILITY OF CO2 IN ML/L FROM AIR
!!                                 AT ONE ATMOSPHERE (WEISS, 1974)
!!      *C01*              REAL    VOLUMETRIC SOLUBILITY CONSTANT A2 FOR
!!                                 THE SOLUBILITY OF CO2 IN ML/L FROM AIR
!!                                 AT ONE ATMOSPHERE (WEISS, 1974)
!!      *C02*              REAL    VOLUMETRIC SOLUBILITY CONSTANT A3 FOR
!!                                 THE SOLUBILITY OF CO2 IN ML/L FROM AIR
!!                                 AT ONE ATMOSPHERE (WEISS, 1974)
!!      *C03*              REAL    VOLUMETRIC SOLUBILITY CONSTANT B1 FOR
!!                                 THE SOLUBILITY OF CO2 IN ML/L FROM AIR
!!                                 AT ONE ATMOSPHERE (WEISS, 1974)
!!      *C04*              REAL    VOLUMETRIC SOLUBILITY CONSTANT B2 FOR
!!                                 THE SOLUBILITY OF CO2 IN ML/L FROM AIR
!!                                 AT ONE ATMOSPHERE (WEISS, 1974)
!!      *C05*              REAL    VOLUMETRIC SOLUBILITY CONSTANT B3 FOR
!!                                 THE SOLUBILITY OF CO2 IN ML/L FROM AIR
!!                                 AT ONE ATMOSPHERE (WEISS, 1974)
!!      *C10*              REAL    COEFF. FOR 1. H2CO3 DISSOC. CONST. AFTER
!!                                 EDMOND AND GIESKES (1970)
!!      *C11*              REAL    COEFF. FOR F1. H2CO3 DISSOC. CONST. AFTER
!!                                 EDMOND AND GIESKES (1970)
!!      *C12*              REAL    COEFF. FOR 1. H2CO3 DISSOC. CONST. AFTER
!!                                 EDMOND AND GIESKES (1970)
!!      *C13*              REAL    COEFF. FOR 1. H2CO3 DISSOC. CONST. AFTER
!!                                 EDMOND AND GIESKES (1970)
!!      *C20*              REAL    COEFF. FOR 2. H2CO3 DISSOC. CONST. AFTER
!!                                 EDMOND AND GIESKES (1970)
!!      *C21*              REAL    COEFF. FOR 2. H2CO3 DISSOC. CONST. AFTER
!!                                 EDMOND AND GIESKES (1970)
!!      *C22*              REAL    COEFF. FOR 2. H2CO3 DISSOC. CONST. AFTER
!!                                 EDMOND AND GIESKES (1970)
!!      *C23*              REAL    COEFF. FOR 2. H2CO3 DISSOC. CONST. AFTER
!!                                 EDMOND AND GIESKES (1970)
!!      *CB0*              REAL    COEFF. FOR 1. H3BO3 DISSOC. CONST. AFTER
!!                                 EDMOND AND GIESKES (1970)
!!      *CB1*              REAL    COEFF. FOR 1. H3BO3 DISSOC. CONST. AFTER
!!                                 EDMOND AND GIESKES (1970)
!!      *CB2*              REAL    COEFF. FOR 1. H3BO3 DISSOC. CONST. AFTER
!!                                 EDMOND AND GIESKES (1970)
!!      *CB3*              REAL    COEFF. FOR 1. H3BO3 DISSOC. CONST. AFTER
!!                                 EDMOND AND GIESKES (1970)
!!      *CW0*              REAL    COEFF. FOR KW (DICKSON AND RILEY, 1979)
!!      *CW1*              REAL    COEFF. FOR KW (DICKSON AND RILEY, 1979)
!!      *CW2*              REAL    COEFF. FOR KW (DICKSON AND RILEY, 1979)
!!                                 (CORRECTED ACCORDING TO B. BACASTOW,
!!                                 PERS. COMMUN., 1988)
!!      *OX0*              REAL    VOLUMETRIC SOLUBILITY CONSTANT A1 FOR
!!                                 THE SOLUBILITY OF OXYGEN IN ML/L FROM
!!                                 MOIST AIR AT ONE ATMOSPHERE (WEISS, 1970)
!!      *OX1*              REAL    VOLUMETRIC SOLUBILITY CONSTANT A2 FOR
!!                                 THE SOLUBILITY OF OXYGEN IN ML/L FROM
!!                                 MOIST AIR AT ONE ATMOSPHERE (WEISS, 1970)
!!      *OX2*              REAL    VOLUMETRIC SOLUBILITY CONSTANT A3 FOR
!!                                 THE SOLUBILITY OF OXYGEN IN ML/L FROM
!!                                 MOIST AIR AT ONE ATMOSPHERE (WEISS, 1970)
!!      *OX3*              REAL    VOLUMETRIC SOLUBILITY CONSTANT A4 FOR
!!                                 THE SOLUBILITY OF OXYGEN IN ML/L FROM
!!                                 MOIST AIR AT ONE ATMOSPHERE (WEISS, 1970)
!!      *OX4*              REAL    VOLUMETRIC SOLUBILITY CONSTANT B1 FOR
!!                                 THE SOLUBILITY OF OXYGEN IN ML/L FROM
!!                                 MOIST AIR AT ONE ATMOSPHERE (WEISS, 1970)
!!      *OX5*              REAL    VOLUMETRIC SOLUBILITY CONSTANT B2 FOR
!!                                 THE SOLUBILITY OF OXYGEN IN ML/L FROM
!!                                 MOIST AIR AT ONE ATMOSPHERE (WEISS, 1970)
!!      *OX6*              REAL    VOLUMETRIC SOLUBILITY CONSTANT B3 FOR
!!                                 THE SOLUBILITY OF OXYGEN IN ML/L FROM
!!                                 MOIST AIR AT ONE ATMOSPHERE (WEISS, 1970)
!!      *T*                REAL    DUMMY VARIABLE, ABSOLUTE SEAWATER TEMP.
!!      *qtt*                REAL    DUMMY VARIABLE, ABSOLUTE SEAWATER TEMP.,
!!                                 DIVIDED BY 100.
!!      *S*                REAL    DUMMY VARIABLE, SALINITY
!!      *CL*               REAL    CHLORINITY (CL(O/OO)=S(O/OO)/1.80655)
!!                                 AFTER WOOSTER ET AL., 1969
!!                                 (C.F. KALLE/DIETRICH , P. 60)
!!      *CEK0*             REAL    LN(K0), LOGARITHM OF CO2 SOLUBILITY IN
!!                                 SEAWATER IN VOLUMETRIC UNITS (EQ. 12 IN
!!                                 WEISS, 1974)
!!      *CK1*              REAL    PK1-VALUE (K1= 1. H2CO3 DISSOC. CONST..),
!!                                 AFTER EDMOND AND GIESKES (1970)
!!      *CK2*              REAL    PK2-VALUE (K2= 2. H2CO3 DISSOC. CONST..),
!!                                 AFTER EDMOND AND GIESKES (1970)
!!      *CKB*              REAL    PKB-VALUE (KB= 1. H3BO3 DISSOC. CONST..)
!!                                 AFTER EDMOND AND GIESKES (1970)
!!      *CKW*              REAL    PKW-VALUE (KW=H2O DISSOC. CONST.) AFTER
!!                                 DICKSON AND RILEY (1979)
!!      *OXY*              REAL    LN(C*), LOGARITHM OF O2 SOLUBILITY IN
!!                                 SEAWATER IN VOLUMETRIC UNITS (EQ. 4 IN
!!                                 WEISS, 1970)
!!      *AK1*              REAL    K1, 1. H2CO3 DISSOC. CONSTANT
!!                                 (EDMOND AND GIESKES, 1970)
!!      *AK2*              REAL    K2, 2. H2CO3 DISSOC. CONSTANT
!!                                 (EDMOND AND GIESKES, 1970)
!!      *AKB*              REAL    KB, 1. H3BO3 DISSOC. CONSTANT
!!                                 (EDMOND AND GIESKES, 1970)
!!      *AKW*              REAL    KW, H2O DISSOC. CONSTANT, LIT ?
!!      *AK0*              REAL    EXP(LN(K0))=K0 CO2 SOLUBILITY IN SEAWATER
!!                                 IN VOLUMETRIC UNITS (ML/L)(WEISS, 1974,
!!                                 CF. EQ. 12)
!!      *RRR*              REAL    SIGMA-T IN OCEAN MODEL, DUMMY VARIABLE
!!                                 (USED FOR CALCULATION OF TOTAL BORAT
!!                                 CONCENTRATION) (SIGMA-T=RHO(S,T,0)/1000.)
!!      *BOR*              REAL    TOTAL BORAT CONCENTRATION , DUMMY VAR.
!!      *P*                REAL    APPROXIMATE PRESSURE AT DEPTH OF U-POINTS
!!                                 IN BAR, DUMMY VARIABLE
!!      *AKSP0*            REAL    CACO3 SOLUBILITY PRODUCT AT P=0 DBAR
!!                                 ACCORDING TO INGLE (1800), EQ. 6; THE
!!                                 CITATION OF CULBERSON AND PYTKOWICZ, 1968,
!!                                 IN BROECKER ET AL., 1982, IS PRESUMABLY
!!                                 NOT CORRECT)
!!      *AKSP(jpi,jpj,jpk)*   REAL    CACO3 SOLUBILITY PRODUCT AT IN SITU PRESSURE
!!                                 FOLLOWING THE PROCEDURE DESCRIBED IN EDMOND
!!                                 AND GIESKES (1970), P. 1285
!!      *CP*               REAL    TERM IN EXPONENT OF EQUATIONS FOR PRESSURE
!!                                 CORRECTION OF DISSOC. CONSTANTS (CARB.,
!!                                 BOR. ACID) AND CALCITE SOLUB. PRODUCT
!!                                 (CF. BROECKER ET AL., 1982, EDMOND AND
!!                                 GIESKES, 1970)
!!      *TC*               REAL    TEMPERATURE AT OCEAN GRID POINTS (DEG C),
!!                                 DUMMY VARIABLE
!!      *KI*               INTEGER COUNTS ITERATIONS FOR NEWTON-RAPHSON METHOD
!!                                 FOR INITIATION OF [CO3--] AND [H+]
!!      *H*                REAL    [H+], DUMMY VARIABLE
!!      *R*                REAL    [CO3--] [MOLE/L], DUMMY VARIABLE
!!      *ALKA*             REAL    GIVEN ALKALINITY [EQV/L], DUMMY VARIABLE
!!      *C*                REAL    GIVEN [SUM(12C)O2] [MOLE/L], DUMMY VARIABLE
!!      *A*                REAL    ALKALINITY [EQV/L] AS FUNCTION OF [CO3--]
!!                                 AND [H+], DUMMY VARIABLE
!!      *DCDS*             REAL    LOCAL DERIVATIVE
!!                                 [SUM(CO2)]([CO3--],H+]) -> [CO3--]
!!      *DADS*             REAL    LOCAL DERIVATIVE
!!                                 [ALK]([CO3--],H+]) -> [CO3--]
!!      *DCDH*             REAL    LOCAL DERIVATIVE
!!                                 [SUM(CO2)]([CO3--],H+]) -> [H+]
!!      *DADH*             REAL    LOCAL DERIVATIVE
!!                                 [ALK]([CO3--],H+]) -> [H+]
!!      *F*                REAL    FUNCTION [ALK] GIVEN MINUS [ALK] IN TERMS
!!                                 OF [CO3--] AND [H+]
!!      *CC*               REAL    [SUM(CO2)] [MOLE/L] AS FUNCTION OF
!!                                 [CO3--] AND [H+], DUMMY VARIABLE
!!      *GG*               REAL    FUNCTION [SUM(CO2)] GIVEN MINUS [SUM(CO2)]
!!                                 IN TERMS OF [CO3--] AND [H+]
!!      *DETI*             REAL    DETERMINANT WITH LOCAL DERIVATIVES FOR
!!                                 NEWTON-RAPHSON ITERATION
!!
      REAL, PUBLIC, SAVE :: akcc1, akcc2, akcc3, akcc4
      REAL, PUBLIC, SAVE :: arafra, calfra, aracal, devk1, devk2, devkb
      REAL, PUBLIC, SAVE :: devk1t, devk2t, devkbt, devkst, devks
      REAL, PUBLIC, SAVE :: bor1, bor2, c00, c01, c02, c03, c04, c05, c10, c11
      REAL, PUBLIC, SAVE :: c12, c13, c20, c21, c22, c23, cb0, cb1, cb2, cb3
      REAL, PUBLIC, SAVE :: c14, c15, c16, c17, c24, c25, c26, c27
      REAL, PUBLIC, SAVE :: cb4, cb5, cb6, cb7, cb8, cb9, cb10, cb11
      REAL, PUBLIC, SAVE :: cw3, cw4, cw5, cw6
      REAL, PUBLIC, SAVE :: cw0, cw1, cw2, ox0, ox1, ox2, ox3, ox4
      REAL, PUBLIC, SAVE :: cek0, ckb, ck1, ck2, ckw, ak1, ak2, ak0
      REAL, PUBLIC, SAVE :: bor, aksp0, ah
      REAL, PUBLIC, SAVE :: smicr, thousi, perc, third
      REAL, PUBLIC, SAVE :: salchl, temzer, rgas, oxyco, ox5
      REAL, PUBLIC, SAVE :: oxy, rrr
      REAL, PUBLIC, SAVE :: qtt
      REAL, PUBLIC, SAVE :: cpexp
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: chemc
!!
!!----------------------------------------------------------------------
!! Time variables
!!
!!    rjjss = seconds in day (=rjjhh*rhhmm*rmmss defined in OPA_SRC/DOM)
!!    raajj = number of days in year
!!    raass =  number of seconds in year
!!    rmoss = number of seconds in month (= raass/12)
!!----------------------------------------------------------------------
!!
      REAL, PUBLIC, SAVE :: rfact, rfactr, xtvit, rjjss, raajj, raass, rmoss
!!
!!----------------------------------------------------------------------
!! Gas exchange variables
!! ---------------------------------------------------------------------
!!
      REAL, PUBLIC, SAVE :: qcumul(jptra)
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: kgwanin
#if defined key_trc_atmco2
!!
!!----------------------------------------------------------------------
!! Atmospheric CO2 data variables
!! -----------------------------------------------------------------
!!
      INTEGER, PUBLIC, PARAMETER :: nmaxrec=5000
      REAL, PUBLIC, SAVE :: yrco2(nmaxrec),sipco2(nmaxrec),pco2at,pn2oat,pch4at
      INTEGER, PUBLIC, SAVE :: nutatm
      REAL, PUBLIC :: yrcfc(336),sipcfc(336,3,2)
#endif
!!----------------------------------------------------------------------
!! Biogeochemical parameters in namelist.trc.sms
!!
!! with rn_km?phy = MICHAELIS MENTEN CONSTANTS OF PRODUCTIVITY (HALF 
!!           SATURATION CONSTANT; NUTRIENT CONCENTRATION, WHERE NUTRIENT
!!           UPTAKE VELOCITY IS HALF ITS POTENTIAL MAXIMUM VALUE;
!!           MICHAELIS & MENTEN(1913)
!!----------------------------------------------------------------------
!
! aggregation parameters
      REAL, PUBLIC, SAVE :: rn_ag1poc,rn_ag2poc,rn_ag3poc,rn_ag4poc,rn_ag5doc,rn_ag6doc,rn_ag7doc
      REAL, PUBLIC, SAVE :: rn_sisgoc,rn_siegoc,rn_snkhoc,rn_snkpoc,rn_snkmax
!
! river fluxes
      REAL, PUBLIC, SAVE :: rn_rivdic,rn_rivdoc,rn_rivfer,rn_rivnit,rn_rivpo4,rn_rivpoc,rn_rivsil
!
! iron parameters
      REAL, PUBLIC, SAVE :: rn_fersol,rn_silsol,rn_sedfer,rn_scofer,rn_scmfer
      REAL, PUBLIC, SAVE :: rn_ligdpn,rn_ligdps,rn_ligfrn,rn_liglat,rn_ligfrs
!
! biological parameters for preferences 
      INTEGER, PUBLIC, SAVE :: jpfoo
      REAL, PUBLIC, SAVE :: rn_bmsdet(4),rn_gbadoc,rn_gbapoc,rn_gbagoc,rn_gbagon,rn_gbahoc
      REAL, PUBLIC, SAVE :: gbadoc,gbapoc,gbagoc,gbagon,gbahoc
      REAL, PUBLIC, SAVE :: rn_readsi,rn_remdsi,rn_retdsi
      REAL, PUBLIC, SAVE :: rn_bmspft(jppoc:jpdia+jppft-1),rn_prfzoo(jpzft,jppoc:jpdia+jppft-1)
      REAL, PUBLIC, SAVE :: prfzoo(jpzft,jppoc:jpdia+jppft-1)
!
! biological parameters for heterotrophs 
      INTEGER, PUBLIC, SAVE :: nn_sizzoo(jpzft),nn_gonmin,nn_gonmax,nn_morsqu,ngochoc
      REAL, PUBLIC, SAVE :: rn_ggebac,rn_ggtbac,rn_ggezoo(jpzft)
      REAL, PUBLIC, SAVE :: rn_sigzoo(jpzft)
      REAL, PUBLIC, SAVE :: rn_unazoo(jpzft)
      REAL, PUBLIC, SAVE :: rn_gramin,rn_grazoo(jpzft),rn_grkzoo(jpzft)
      REAL, PUBLIC, SAVE :: rn_grabac
      REAL, PUBLIC, SAVE :: rn_kmobac,rn_kmfbac,rn_kmpbac,rn_denno3
      REAL, PUBLIC, SAVE :: rn_morcru,rn_motcru
      REAL, PUBLIC, SAVE :: rn_morgel,rn_mosgel,rn_motgel
      REAL, PUBLIC, SAVE :: rn_resbac,rn_reszoo(jpzft)
      REAL, PUBLIC, SAVE :: rn_retbac,rn_retzoo(jpzft)
      REAL, PUBLIC, SAVE :: rn_icecru
!
! parameters for zooplankton calcifiers
      REAL, PUBLIC, SAVE :: rn_forcal,rn_lysara,rn_pteara
!
! biological parameters for phytoplankton 
      REAL, PUBLIC, SAVE :: rn_bsidia,rn_coccal,rn_disara,rn_discal,rn_disfor
      REAL, PUBLIC, SAVE :: rn_ferbsi,rn_kmsbsi,rn_lyoco3,rn_lyscal,rn_munfix
      REAL, PUBLIC, SAVE :: rn_rembac,rn_silbsi,rn_sildia
      REAL, PUBLIC, SAVE :: rn_resphy(jpdia:jpdia+jppft-1)

      REAL, PUBLIC, SAVE :: rn_docphy(jpdia:jpdia+jppft-1),rn_domphy(jpdia:jpdia+jppft-1)
      REAL, PUBLIC, SAVE :: rn_kmnphy(jpdia:jpdia+jppft-1),rn_kmpphy(jpdia:jpdia+jppft-1)
      REAL, PUBLIC, SAVE :: rn_mumpft(jpdia:jpdia+jppft-1)
      REAL, PUBLIC, SAVE :: pcmax(jpdia:jpdia+jppft-1)
!
! temperature-dependence for all PFTs
      REAL, PUBLIC, SAVE :: rn_mutpft(jpbac:jpdia+jppft-1),rn_mudpft(jpbac:jpdia+jppft-1)
!
! parameters for the Fe-light model 
      REAL, PUBLIC, SAVE :: rn_kmfphy(jpdia:jpdia+jppft-1), rn_rhfphy(jpdia:jpdia+jppft-1)      
      REAL, PUBLIC, SAVE :: rn_qmaphy(jpdia:jpdia+jppft-1), rn_qmiphy(jpdia:jpdia+jppft-1)    
      REAL, PUBLIC, SAVE :: rn_qopphy(jpdia:jpdia+jppft-1)
!
! other light parameters 
      REAL, PUBLIC, SAVE :: rn_ekwgrn,rn_ekwred
      REAL, PUBLIC, SAVE :: rn_alpphy(jpdia:jpdia+jppft-1)
      REAL, PUBLIC, SAVE :: rn_kgrphy(jpdia:jpdia+jppft-1),rn_krdphy(jpdia:jpdia+jppft-1)
      REAL, PUBLIC, SAVE :: rn_thmphy(jpdia:jpdia+jppft-1)
#    if defined key_trc_dms
!
! parameters for the DMS model
      REAL, PUBLIC, SAVE :: rn_dmsyld, rn_etomax
      REAL, PUBLIC, SAVE :: rn_rdddms 
      REAL, PUBLIC, SAVE :: rn_xcldmd, rn_xpodms
      REAL, PUBLIC, SAVE :: rn_xkdms, rn_xkdmd
      REAL, PUBLIC, SAVE :: rn_assdms(jpzft)
      REAL, PUBLIC, SAVE :: rn_rphdmd(jpdia:jpdia+jppft-1)
      REAL, PUBLIC, SAVE :: rn_xpldmd(jpdia:jpdia+jppft-1)
#    endif
#    if defined key_trc_n2o
!
! parameters for the N2O model
      INTEGER, PUBLIC, SAVE :: nn_deun2s
      REAL, PUBLIC, SAVE :: rn_aoun2s,rn_betn2s,rn_decn2s,rn_omxn2s
      REAL, PUBLIC, SAVE :: sipn2o(nmaxrec)
#    endif
#    if defined key_trc_ch4
      REAL, PUBLIC, SAVE :: sipch4(nmaxrec)
      REAL, PUBLIC, SAVE :: rn_proch1,rn_proch2,rn_proch3,rn_proch4,rn_proch5
      REAL, PUBLIC, SAVE :: rn_botch1,rn_botch2,rn_botch3,rn_botch4,rn_botch5
      REAL, PUBLIC, SAVE :: rn_conch1,rn_conch2,rn_conch3,rn_conch4,rn_conch5
      REAL, PUBLIC, SAVE :: rn_decch4
#      if ! defined key_trc_n2o
      REAL, PUBLIC, SAVE :: sipn2o(nmaxrec)
#      endif
#    endif
!!
!!----------------------------------------------------------------------
!! Biological parameters inherited from HAMOCC3
!!
!! calcon = MEAN TOTAL [CA++] IN SEAWATER (MOLES/KG) 
!!             (SEE BROECKER A. PENG, 1982, P. 26)
!!             ([CA++](MOLES/KG)=1.026E-2*(S/35.) AFTER
!!             CULKIN(1965), CF. BROECKER ET AL. 1982)
!! alknut = REDFIELD RATIO MOLES N + TPO4 / MOLES C
!!             (FOR CHANGE IN ALKALINITY DUE TO PRODUCTION/REMINE-
!!             RALIZATION OF ORGANIC MATTER)
!!             N:C=16:122, SEE Wolf-Gladrow et al. 2007
!!
!!             COMPOSITION OF PLANKTONIC MATERIAL AND RESPIRATION:
!!             C106H263O110N16P + 138 O2 ->
!!            -> 106 CO2 + 16 NO3- + HPO4(2-) + 122 H2O + 18 H+  ,
!!             SEE DEGENS ET AL. (1984), P. 152
!! *PO4R*      REAL  RATIO (MOLES P)/(MOLES C) (REDFIELD RATIO P:C)
!! ---------------------------------------------------------------------
!! 
      REAL, PUBLIC, SAVE :: alknut,ratn2c, ratc2n, rato2c
      REAL, PUBLIC, SAVE :: calcon, plafr13, pdb, minfer
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: fld,flu,flu16,zkgco2
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: zkgo2,dpco2,pco2,cflx
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: prcaca,proara,prorca3,denitr
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: losara,loscal,losfor
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:  ) :: dust,mdept
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: dustmo
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: tgfunc
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: depdic,depdoc,deppoc
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: depnit,deppo4,depsil
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: depfer
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: atmdin
      REAL, PUBLIC, SAVE :: extinp(7),sedcor(7)
!!
!!----------------------------------------------------------------------
!! Biological parameters inherited from PISCES
!!----------------------------------------------------------------------
!!
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: bactge,delo2,etot
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: snkara,snkbfe,snkcal
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: snkdsi,snkgoc,snkgon,snkhoc
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: snkpoc,snksfe,snkufe
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: xagg,xaggfe
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: xaggdoc
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: xscave,remdoc
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: rempoc,remgoc,remgon,remhoc
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: remsfe,rembfe,remufe
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: remdsi,xaggdoc2
!!
!!----------------------------------------------------------------------
!! Biological parameters specific to PlankTOM
!!----------------------------------------------------------------------
!!
      REAL, PUBLIC, SAVE :: dnsmin,ferat3,grizoo(2,jpfmx),ncf_fill
      REAL, PUBLIC, SAVE :: moworg,mowco3,mwsio2,rhoorg,rhoco3,rhsio2
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:):: snkmax
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: grazoo
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: grazoc,grazof
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: grapoc,grarem,grafer
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: losbsi
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: irondep,sidep
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: mgezoo
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: docphy
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: prophy
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: resphy
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: resbac
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: stofoo
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: reszoo
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: torcru,torgel
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: rbafer,ubafer
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: dinpft
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: xlimbac
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: volumt
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: icecru,trncru,rcoast
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: prodt,ppt,pptdoc
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: tchl,out3d
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: dyphy,eliblu,elired
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: xvsink,zstre1
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: grazing,trophic
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: out2d,ppint
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: xlimpft
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: ligfer
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: trnsed,remsed
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: lim2mmfe,lim3fe
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: lim4po4,lim5si
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: lim6din,lim8light
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: obstempmo,obssalmo
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: tsnbio,tsnche
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:)   :: obstemp,obssal
      INTEGER, PUBLIC, SAVE :: nn_biodts
!
#    if defined key_trc_dms
!!
!!-----------------------------------------------------------------
!! Variables for the prognostic DMS model of Vogt et al. JGR (2010)
!!-----------------------------------------------------------------
!!
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: fludms,zkgdms
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: prodms, prodmd
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: degdms, degdmd
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: dmddms
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: dmspp
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: rphdmd
      REAL, PUBLIC, SAVE :: dms_snkcount,dmd_snkcount
#    endif
#    if defined key_trc_n2o
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: zkgn2o,flun2s,dpn2s
#    endif
#    if defined key_trc_ch4
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: zkgch4,fluch1,fluch2,fluch3,fluch4,fluch5,dpch4
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: fluch6,fluch7,fluch8,fluch9,fluch10
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: fluch11,fluch12,fluch13,fluch14,fluch15
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: fluch16,fluch17,fluch18,fluch19,fluch20
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: fluch21,fluch22,fluch23,fluch24,fluch25
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: proch1,proch2,proch3,proch4,proch5,proch6,proch7,proch8
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: proch9,proch10,proch11,proch12,proch13,proch14
#    endif
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:) :: zkgcfc11,fluc11
#    if defined key_c14b
      REAL, PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: d14pro,d14res
#    endif

CONTAINS

   INTEGER FUNCTION sms_planktom_alloc()
      !!----------------------------------------------------------------------
      !!        *** ROUTINE sms_planktom_alloc ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_warn
      INTEGER ::   ierr(10)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      !* biological fluxes for light
      ALLOCATE( etot(jpi,jpj,jpk), STAT=ierr(1) )
      !
      !* CO2 cycle
      ALLOCATE( akb3(jpi,jpj,jpk)  , ak13(jpi,jpj,jpk), ak23(jpi,jpj,jpk)   ,  &
             &  aksara(jpi,jpj,jpk), aksp(jpi,jpj,jpk), co3(jpi,jpj,jpk)    ,  &
             &  hi(jpi,jpj,jpk)    ,                                           &
#  if defined key_trc_piic
             &  pihi(jpi,jpj,jpk)  , pih2co3(jpi,jpj) ,                        &
#  endif
             &  borat(jpi,jpj,jpk) , hco3(jpi,jpj,jpk), h2co3(jpi,jpj)      ,  &
             &  akw3(jpi,jpj,jpk)  , chemc(jpi,jpj,6) , discarb(jpi,jpj,jpk), STAT=ierr(2) )
      !
      !* gas exchange
      ALLOCATE( kgwanin(jpi,jpj),  STAT=ierr(3) )
      !
      !* biological parameters inherited from HAMOCC3
      ALLOCATE( fld(jpi,jpj)       , flu(jpi,jpj)     , flu16(jpi,jpj)     ,   &
             &  zkgco2(jpi,jpj)    , zkgo2(jpi,jpj)   ,                        &
             &  dpco2(jpi,jpj)     , pco2(jpi,jpj)    , cflx(jpi,jpj)      ,   &
             &  prcaca(jpi,jpj,jpk), proara(jpi,jpj,jpk),                      &
             &  prorca3(jpi,jpj,jpk),                                          &
             &  losara(jpi,jpj,jpk), loscal(jpi,jpj,jpk),losfor(jpi,jpj,jpk) , &
             &  denitr(jpi,jpj,jpk), dust(jpi,jpj),dustmo(jpi,jpj,12),         &
             &  mdept(jpi,jpj),      tgfunc(jpi,jpj,jpk,jpbac:jpdia+jppft-1) , &
             &  depdic(jpi,jpj,jpk), depdoc(jpi,jpj,jpk), deppoc(jpi,jpj,jpk), &
             &  depnit(jpi,jpj,jpk), deppo4(jpi,jpj,jpk), depsil(jpi,jpj,jpk), &
             &  depfer(jpi,jpj,jpk), atmdin(jpi,jpj,jpk), STAT=ierr(4) )
      !
      !* biological parameters inherited from pisces
      ALLOCATE( snkgoc(jpi,jpj,jpk), snkcal(jpi,jpj,jpk), delo2(jpi,jpj,jpk),  &
             &  snkara(jpi,jpj,jpk), snkpoc(jpi,jpj,jpk),                      &
             &  snkgon(jpi,jpj,jpk), snkhoc(jpi,jpj,jpk), snksfe(jpi,jpj,jpk), &
             &  snkbfe(jpi,jpj,jpk), snkufe(jpi,jpj,jpk), snkdsi(jpi,jpj,jpk), &
             &  xagg(jpi,jpj,jpk),   xaggfe(jpi,jpj,jpk),                      &
             &  xaggdoc(jpi,jpj,jpk), bactge(jpi,jpj,jpk),                     &
             &  xscave(jpi,jpj,jpk), remdoc(jpi,jpj,jpk),                      &
             &  rempoc(jpi,jpj,jpk), remgoc(jpi,jpj,jpk), remgon(jpi,jpj,jpk), &
             &  remhoc(jpi,jpj,jpk), remufe(jpi,jpj,jpk),                      &
             &  remsfe(jpi,jpj,jpk), rembfe(jpi,jpj,jpk),                      &
             &  remdsi(jpi,jpj,jpk), xaggdoc2(jpi,jpj,jpk), STAT=ierr(5) )
     !
     !* biological parameters specific to PlankTOM
      ALLOCATE( snkmax(jpk)        , grafer(jpi,jpj,jpk,jpzft),                &
             &  grazoo(jpi,jpj,jpk,jpzft,jppoc:jpdia+jppft-1),                 &
             &  grazoc(jpi,jpj,jpk,jpzft), grazof(jpi,jpj,jpk,jpzft),          &
             &  grapoc(jpi,jpj,jpk,jpzft),                                     &
             &  grarem(jpi,jpj,jpk,jpzft),                                     &
             &  losbsi(jpi,jpj,jpk),                                           &
             &  irondep(jpi,jpj,jpk), sidep(jpi,jpj,jpk),                      &
             &  mgezoo(jpi,jpj,jpk,jpzft),                                     &
             &  docphy(jpi,jpj,jpk,jpdia:jpdia+jppft-1),                       &
             &  prophy(jpi,jpj,jpk,jpdia:jpdia+jppft-1,3),                     &
             &  resphy(jpi,jpj,jpk,jpdia:jpdia+jppft-1,3),                     &
             &  resbac(jpi,jpj,jpk),                                           &
             &  stofoo(jpi,jpj,jpk,jppoc:jpdia+jppft-1,3),                     &
             &  reszoo(jpi,jpj,jpk,jpzft),                                     &
             &  torcru(jpi,jpj,jpk),                                           &
             &  torgel(jpi,jpj,jpk),                                           &
             &  rbafer(jpi,jpj,jpk), ubafer(jpi,jpj,jpk),                      &
             &  dinpft(jpi,jpj,jpk,jpdia:jpdia+jppft-1),                       &
             &  lim2mmfe(jpi,jpj,jpk,jppft),                                   &
             &  lim3fe(jpi,jpj,jpk,jppft),                                     &
             &  lim4po4(jpi,jpj,jpk,jppft),                                    &
             &  lim5si(jpi,jpj,jpk,jppft),                                     &
             &  lim6din(jpi,jpj,jpk,jppft),                                    &
             &  lim8light(jpi,jpj,jpk,jppft),                                  &
             &  xlimbac(jpi,jpj,jpk), volumt(jpi,jpj,jpk),                     &
             &  grazing(jpi,jpj,jpk,5),                                        &
             &  prodt(jpi,jpj,jpk) , ppint(jpi,jpj)     , ppt(jpi,jpj,jpk),    &
             &  pptdoc(jpi,jpj,jpk), trophic(jpi,jpj,jpk,3),                   &
             &  tchl(jpi,jpj,jpk)  , out2d(jpi,jpj)     , out3d(jpi,jpj,jpk),  &
             &  dyphy(jpi,jpj,jpk) , eliblu(jpi,jpj,jpk), elired(jpi,jpj,jpk), &
             &  xlimpft(jpi,jpj,jpk,jpdia:jpdia+jppft-1), ligfer(jpj,jpk),     &
             &  zstre1(jpi,jpj,jpk), xvsink(jpi,jpj,jpk),                      &
             &  trnsed(jpi,jpj,jpdsi:jphoc),                                   &
             &  remsed(jpi,jpj,jpdsi:jphoc),                                   &
             &  obstempmo(jpi,jpj,jpk,12), obssalmo(jpi,jpj,jpk,12),           &
             &  tsnbio(jpi,jpj,jpk,2), tsnche(jpi,jpj,jpk,2),                  &
             &  obstemp(jpi,jpj,jpk), obssal(jpi,jpj,jpk),                     &
             &  icecru(jpi,jpj),trncru(jpi,jpj),rcoast(jpi,jpj),               &
             &  zkgcfc11(jpi,jpj),fluc11(jpi,jpj),  STAT=ierr(6) )
     !
#    if defined key_trc_dms
     !* variables for DMSD model
      ALLOCATE( fludms(jpi,jpj)    , zkgdms(jpi,jpj)    ,                      &
             &  prodms(jpi,jpj,jpk), prodmd(jpi,jpj,jpk),                      &
             &  degdms(jpi,jpj,jpk), degdmd(jpi,jpj,jpk),                      &
             &  dmddms(jpi,jpj,jpk), dmspp(jpi,jpj,jpk) ,                      &
             &  rphdmd(jpi,jpj,jpk,jpdia:jpdia+jppft-1) ,    STAT=ierr(7) )
#    endif
     !
#    if defined key_trc_n2o
     !* variables for N2O model
      ALLOCATE( zkgn2o(jpi,jpj)    , flun2s(jpi,jpj)    , dpn2s(jpi,jpj), STAT=ierr(8) )
#    endif
     !
#    if defined key_trc_ch4
     !* variables for Ch4 model
      ALLOCATE( zkgch4(jpi,jpj)    , fluch1(jpi,jpj)    , fluch2(jpi,jpj)  ,   &
      &         fluch3(jpi,jpj)    , fluch4(jpi,jpj)    , fluch5(jpi,jpj)  ,   &
      &         fluch6(jpi,jpj)    , fluch7(jpi,jpj)    , fluch8(jpi,jpj)  ,   &
      &         fluch9(jpi,jpj)    , fluch10(jpi,jpj)   ,                      &
      &         fluch11(jpi,jpj)   , fluch12(jpi,jpj)   , fluch13(jpi,jpj) ,   &
      &         fluch14(jpi,jpj)   , fluch15(jpi,jpj)   ,                      &
      &         fluch16(jpi,jpj)   , fluch17(jpi,jpj)   , fluch18(jpi,jpj) ,   &
      &         fluch19(jpi,jpj)   , fluch20(jpi,jpj)   ,                      &
      &         fluch21(jpi,jpj)   , fluch22(jpi,jpj)   , fluch23(jpi,jpj) ,   &
      &         fluch24(jpi,jpj)   , fluch25(jpi,jpj)   ,                      &
      &         dpch4(jpi,jpj)     , proch1(jpi,jpj,jpk),proch2(jpi,jpj,jpk),  &
      &         proch3(jpi,jpj,jpk), proch4(jpi,jpj,jpk),                      &
      &         proch5(jpi,jpj,jpk), proch6(jpi,jpj,jpk),proch7(jpi,jpj,jpk),  &
      &         proch8(jpi,jpj,jpk), proch9(jpi,jpj,jpk),                      &
      &         proch10(jpi,jpj,jpk),proch11(jpi,jpj,jpk),proch12(jpi,jpj,jpk),&
      &         proch13(jpi,jpj,jpk),proch14(jpi,jpj,jpk),   STAT=ierr(9) )
#    endif
!
#    if defined key_c14b
     !* variables for 14C isotopes
      ALLOCATE( d14pro(jpi,jpj,jpk), d14res(jpi,jpj,jpk),  STAT=ierr(10) )
#    endif
      sms_planktom_alloc = MAXVAL( ierr )
      !
      IF( sms_planktom_alloc /= 0 )   CALL ctl_warn('sms_planktom_alloc: failed to allocate arrays')
 END FUNCTION sms_planktom_alloc
#endif
   !!======================================================================  
END MODULE sms_planktom
