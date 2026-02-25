      SUBROUTINE bgcche
#if defined key_planktom && defined key_top
!!======================================================================
!!                         ***  ROUTINE bgcche  ***
!! TOP : PlankTOM ecosystem model 
!!======================================================================
!! History : original h3cche :      1988 E. Maier-Reimer
!!           additions :     1998 O. Aumont
!!           modifications : 1999 C. Le Quere
!!----------------------------------------------------------------------
!!    'key_trc_n2o'       :                        calculates N2O fluxes
!!    'key_trc_ch4'       :                        calculates CH4 fluxes
!!----------------------------------------------------------------------
!! ** Purpose : sets chemical constants
!! ** Action  : 
!!
!!    Method: 
!!          1) set constants for carbonate chemistry as described in
!!             in Broecker et al. (1982, geosecs) and Edmond A. Gieskes
!!             (1970)
!!          2) initiate [co3--] and ph-value by iteration
!!             (Newton-Raphson method for solving nonlinear simultaneous
!!              equations, see E.G. Scarborough, J. (1958))
!!
!!    Externals: 
!!          *RHO*     - half precision function, eq. of state oF
!!                      seawater
!!
!!    References: 
!!
!!          Berner, R. A. (1976)
!!          The solubility of calcite and aragonite in sea water
!!          at atmospheric pressure and 34.5 o/oo salinity.
!!          American Journal of Science, vol. 276, 713-730.
!!          (k'sp(aragonite)=1.45 k'sp(calcite))
!!
!!          Broecker, W.S., D.W. Spencer, AND H. Craig (1982)
!!          Geosecs Pacific Expedition. VOL. 3.. Hydrographic Data
!!          1973-1974, Superintendant of documents, U.S. Government
!!          printing office, Washington, D.C., 137 PP..
!!
!!          Culberson, C.H., AND R.M. Pytkowicz (1968)
!!          effect on pressure on carbonic acid, boric acid and the ph
!!          in sea water.
!!          Limnology and Oceanography, Vol. 13, 403-417.
!!
!!          Dickson, A.G., AND J.P. Riley (1979)
!!          The estimation of acid dissociation constants in seawater
!!          media from potentiometric titrations with strong base.
!!          I. The ionic product of water - KW.
!!          Marine Chemistry, vol. 7, 89-99.
!!
!!          Edmond, J.M., AND J.M.T.M. Gieskes (1970)
!!          on the calculation of the degree of saturation of sea water
!!          with respect to calcium carbonate under in situ conditions.
!!          Geochim. et Cosmochim. ACTA, 34, 1261-1291.
!!
!!          Ingle, S.E. (1800)
!!          Solubility of calcite in the ocean.
!!          Marine Chemistry, vol. 3, 301-319.
!!
!!          Ingle, S.E., C.H. Culberson, J.E. Hawley, AND R.M. Pytkowicz
!!          (1973) The solubility of calcite in seawater at atmospheric
!!          pressure and 35 O/OO salinity.
!!          Marine Chemistry, vol. 1, 295-307.
!!
!!          Riley, J. P., AND G. Skirrow, eds. (1965)
!!          Chemical Oceanography. vol. 1, 712 PP., Academic Press,
!!          London A. New York.
!!
!!          Scarborough, J. (1958) Numerical mathematical analysis.
!!          Oxford University Press, London, 4TH ED., 576 PP..
!!
!!          Weiss, R. F. (1970) The solubility of nitrogen
!!          oxygen and argon in water and seawater.
!!          Deep-Sea Research, vol. 17, 721-735.
!!
!!          Weiss, R. F. (1974)
!!          Carbon dioxide in water and seawater: the solubility of a
!!          non ideal gas. Marine Chemistry, vol. 2, 203-215.
!!
!!          Wooster, W.S., A.J. Lee, and G. Dietrich (1969)
!!          redefinition of salinity. Z. Geophys., vol.35, 611-613.
!!
!!          Broecker, W.S., D.W. Spencer, AND H. Craig (1982)
!!          GEOSECS Pacific expedition. vol. 3.. Hydrographic data
!!          1973-1974, Superintendant of documents, U.S. government
!!          printing office, Washington, D.C., 137 PP..
!!
!!
!!      VARIABLE           TYPE    PURPOSE.
!!      --------           ----    --------
!!
!!       pres              REAL    approximate pressure at depth of u-points
!!                                 in bar, dummy variablE
!!       TC                REAL    temperature at ocean grid points (deg c),
!!                                 dummy variable
!!       cl                REAL    chlorinity (cl(o/oo)=s(o/oo)/1.80655)
!!                                 after wooster et al., 1969
!!                                 (c.f. kalle/dietrich , p. 60)
!!       akw               REAL    kw, h2o dissoc. constant, lit ?
!!       h                 REAL    [h+], dummy variable
!!       rrr               REAL    [co3--] [mole/l], dummy variable
!!       c                 REAL    given [sum(12c)o2] [mole/l], dummy variable
!!       a                 REAL    alkalinity [eqv/l] as function of [co3--]
!!                                 and [h+], dummy variable
!!----------------------------------------------------------------------
! parameters and commons
      USE trp_trc
      USE sms_planktom
      USE oce_trc
      IMPLICIT NONE
! 
! local variables
      INTEGER ji, jj, jk
      REAL tkel, sal
      REAL pres, tc, r, cl
      REAL akw, a, c, h, akb
      REAL zsqrt, ztr, zlogt
      REAL zqtt,qtt2,sal2,sal15,tscale
#    if defined key_trc_n2o
      REAL n2osol
#    endif
#    if defined key_trc_ch4
      REAL ch4sol
#    endif
      REAL ctcsol,cfc11sol,cfc12sol
!
! Chemical constants - surface layer
! ----------------------------------
!
! vertical slab
! -------------
!
      DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
!
! set absolute temperature
          tkel = tsnche(ji,jj,1,1)+temzer
          qtt  = tkel*perc
          qtt2 = qtt**2
          sal  = tsnche(ji,jj,1,2) + (1.-tmask(ji,jj,1))*35.
          zqtt = log(qtt)
!
! chlorinity (Wooster et al., 1969)
          cl = sal*salchl
!
! ln(k0) of solubility of co2 (eq. 12, Weiss, 1974)
          cek0 = c00+c01/qtt+c02*zqtt+sal*(c03+c04*qtt+c05*qtt2)
          ak0 = exp(cek0)*smicr
          chemc(ji,jj,1) = ak0
!
! ln(k0) of solubility of o2 and n2 (eq. 4, Weiss, 1970)
          oxy = ox0+ox1/qtt+ox2*zqtt+sal*(ox3+ox4*qtt+ox5*qtt2)
          chemc(ji,jj,3) = exp(oxy)*oxyco
#    if defined key_trc_n2o
!
! solubility of N2O (Weiss & Price 1980)
          n2osol = -62.7062+97.3066/qtt+24.1406*zqtt+sal*(-0.058420     &
     &      +0.0331983*qtt-0.0051313*qtt2)
          chemc(ji,jj,4) = exp(n2osol)
#    endif
#    if defined key_trc_ch4
!
! solubility of CH4 (Wanninkhof et al. 2014)
          ch4sol = -415.2807+596.8104/qtt+379.2599*zqtt-62.0757*qtt     &
     &      +sal*(-0.059160+0.032174*qtt-0.0048198*qtt2)
          chemc(ji,jj,5) = exp(ch4sol)*1e-9
#    endif
          cfc11sol = -134.1536+203.2156/qtt+56.232*zqtt+ &
       &     sal*(-0.144449+0.092952*qtt-0.0159977*qtt2)
          chemc(ji,jj,6) = exp(cfc11sol)
!
! density of seawater and total borate in mole/
          rrr = rhop(ji,jj,1) *thousi
          bor = bor1*rrr*cl*bor2
          chemc(ji,jj,2) = bor
        END DO
      END DO
!
! chemical constants - deep ocean
      DO jk = 1,jpk
!
! approx. seawater pressure at u-point depth (bar)
        pres = 1.025e-1*gdept_0(ji,jj,jk)
!
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
!
! set limits for seawater temp. and salinity
! (this is done to avoid computational crash at dry
! points during calculation of chemical constants)
!
! set absolute temperature
            tkel   = tsnche(ji,jj,jk,1)+temzer
            qtt    = tkel*perc
            sal    = tsnche(ji,jj,jk,2) + (1.-tmask(ji,jj,jk))*35.
            zsqrt  = sqrt(sal)
            sal15  = sal**1.5
            zlogt  = log(tkel)
            ztr    = 1./tkel
!
! chlorinity (Wooster et al., 1969)
            cl = sal*salchl
!
! ln(k0) of solubility of co2 (eq. 12, Weiss, 1974)
            cek0 = c00+c01/qtt+c02*alog(qtt)+ &
     &          sal*(c03+c04*qtt+c05*qtt**2)
!
!  coefficient ocmip 
            ckb = (cb0+cb1*zsqrt+cb2*sal+cb3*sal15+cb4*sal**2)*ztr &
     &          +(cb5+cb6*zsqrt+cb7*sal)+ &
     &          (cb8+cb9*zsqrt+cb10*sal)*zlogt+cb11*zsqrt*tkel
            ck1 = c10*ztr+c11+c12*zlogt+(c13*ztr+c14)*zsqrt+ &
     &          c15*sal+c16*sal15+log(1.+c17*sal)
            ck2 = c20*ztr+c21+c22*zlogt+(c23*ztr+c24)*zsqrt+c25*sal &
     &          +c26*sal15+log(1.+c27*sal)
!
! pkw (h2o) (Dickson and Riley, 1979)
            ckw = cw0*ztr+cw1+cw2*zlogt+(cw3*ztr+cw4+cw5*zlogt)* &
     &          zsqrt+cw6*sal
!
! ln(k0) of solubility of o2 (eq. 4, Weiss, 1970)
            oxy = ox0+ox1/qtt+ox2*alog(qtt)+sal*(ox3+ox4*qtt+ox5*qtt**2)
!
! k1, k2 of carbonic acid, kb of boric acid, kw (h2o) (lit.?)
            ak1 = exp(ck1)
            ak2 = exp(ck2)
            akb = exp(ckb)
            akw3(ji,jj,jk) = exp(ckw)
!
! apparent solubility product k'sp of calcite in seawater
!       (S=27-43, T=2-25 DEG C) at pres=0 (atmosph. pressure)
!       (Ingle, 1800, eq. 6)
            aksp0 = 1.E-7*(akcc1+akcc2*sal**(third)+akcc3*alog10(sal) &
     &          +akcc4*tkel**2)
!
! formula for cpexp after Edmond and Gieskes (1970)
!        (reference to culberson and pytkoqicz (1968) as made
!        in broecker et al. (1982) is incorrect; here rgas is
!        taken tenfold to correct for the notation of pres  in
!        dbar instead of bar and the expression for cpexp is
!        multiplied by ln(10.) to allow use of exp-function
!        with basis e in the formula for akspp (cf. edmond
!        and gieskes (1970), p. 1285 and p. 1286 (the small
!        formula on p. 1286 is right and consistent with the
!        sign in partial molar volume change as shown on
!        p. 1285))
            cpexp = pres /(rgas*tkel)
!
! kb of boric acid, k1,k2 of carbonic acid pressure
!        correction after culberson and pytkowicz (1968)
!        (cf. broecker et al., 1982)
            tc = tsnche(ji,jj,jk,1) + (1.-tmask(ji,jj,jk))*20.
            akb3(ji,jj,jk) = akb*exp(cpexp*(devkb-devkbt*tc))
            ak13(ji,jj,jk) = ak1*exp(cpexp*(devk1-devk1t*tc))
            ak23(ji,jj,jk) = ak2*exp(cpexp*(devk2-devk2t*tc))
!
! apparent solubility product k'sp of calcite (or aragonite)
!        as function of pressure follwing edmond and gieskes (1970)
!        (p. 1285) and berner (1976)
            aksp(ji,jj,jk) = aracal*aksp0*exp(cpexp*(devks-devkst*tc))
!
! solubility Aragonite, based on arafra=1. see trcini.dgom.F90
            aksara(ji,jj,jk) = 1.45*aksp0*exp(cpexp*(32.8-devkst*tc))
!
! density of seawater and total borate concentr. [moles/l]
            rrr = rhop(ji,jj,jk)*thousi
            bor = bor1*rrr*cl*bor2
            borat(ji,jj,jk) = bor
!
          END DO
        END DO
      END DO
#endif
      RETURN
      END
