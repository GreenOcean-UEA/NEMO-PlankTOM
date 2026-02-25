MODULE trp_trc

#if defined key_planktom
   !!======================================================================
   !! Module trp_trc
   !!======================================================================
   !!  TOP 1.0,  LOCEAN-IPSL (2005)
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! passive tracers number
   USE par_trc , ONLY : &
      jptra    =>   jptra       !!: number of passive tracers

#if defined key_trc_diaadd
   USE par_trc , ONLY : &
      jpdia2d  =>  jpdia2d , &  !!: number of 2d outputs
      jpdia3d  =>  jpdia3d
#endif

   !! passive tracers fields 
   USE trc , ONLY :  &
      trai   =>   trai , &  !!: initial total tracer
      trb    =>   trb  , &  !!: tracer field (before)
      tra    =>   tra  , &  !!: tracer field (now)
      trn    =>   trn       !!: tracer field (after)

#if defined key_trc_diaadd && ! key_iomput
   USE trc , ONLY :  &
      trc2d   =>   trc2d , &  !!: additional 2D variable for outputs
      trc3d   =>   trc3d      !!: additional 3D variable for outputs
#endif
   !! time step - not used in PlankTOM so not required
   !!   USE trc , ONLY :  &
   !!   nn_dttrc =>   nn_dttrc    !!: frequency of step on passive tracers (NAMELIST)

   !! non-centered advection scheme (smolarkiewicz)
   USE trc , ONLY : &
      rtrn   =>   rtrn      !!: value for truncation (NAMELIST)

   USE trc , ONLY : &
      ctrcnm   =>   ctrcnm      !!: value for truncation (NAMELIST)
#else
   !!======================================================================
   !!  Empty module : No passive tracer 
   !!======================================================================
#endif

END MODULE trp_trc
