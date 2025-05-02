MODULE trcisf
   !!==============================================================================
   !!                       ***  MODULE  traisf  ***
   !! Ocean active tracers:  ice shelf boundary condition
   !!==============================================================================
   !! History :    4.0  !  2019-09  (P. Mathiot) original file
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   trc_isf       : update the tracer trend at ocean surface
   !!----------------------------------------------------------------------
   USE isf_oce                                       ! Ice shelf variables
   USE isftrc_oce, ONLY : risfcpl_trc, risfcpl_cons_trc ! Ice shelf variables
   USE par_oce   , ONLY : nijtile, ntile, ntsi, ntei, ntsj, ntej
   USE dom_oce                                       ! ocean space domain variables
   USE isfutils  , ONLY : debug                      ! debug option
   USE timing    , ONLY : timing_start, timing_stop  ! Timing
   USE in_out_manager                                ! I/O manager
   USE trc       , ONLY : ctrcnm
   USE par_trc   , ONLY : jptra

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_isf   ! routine called by step.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: trasbc.F90 10499 2019-01-10 15:12:24Z deazer $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_isf ( kt, Kmm, ptr, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_isf  ***
      !!
      !! ** Purpose :  Compute the temperature trend due to the ice shelf melting (qhoce + qhc)
      !!
      !! ** Action  : - update pts(:,:,:,:,Krhs) for cav, par and cpl case
      !!----------------------------------------------------------------------
      INTEGER                                  , INTENT(in   ) :: kt        ! ocean time step
      INTEGER                                  , INTENT(in   ) :: Kmm, Krhs ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk,jptra,jpt), INTENT(inout) :: ptr       ! passive tracers and RHS of tracer equation
      !!----------------------------------------------------------------------
      INTEGER       :: jn      ! passive tracer indices
      CHARACTER(35) :: charout ! text for debug call
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_isf')
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'trc_isf : Ice shelf -- passive tracer update '
            IF(lwp) WRITE(numout,*) '~~~~~~~ '
         ENDIF
      ENDIF
      !
      ! ice sheet coupling case by default
      !
      ! Dynamical stability at start up after change in under ice shelf cavity geometry is achieve by correcting the divergence.
      ! This is achieved by applying a volume flux in order to keep the horizontal divergence after remapping
      ! the same as at the end of the latest time step. So correction need to be apply at nit000 (euler time step) and
      ! half of it at nit000+1 (leap frog time step).
      ! in accordance to this, the heat content flux due to injected water need to be added in the temperature and salt trend
      ! at time step nit000 and nit000+1
      IF ( kt == nit000  ) CALL trc_isf_cpl(Kmm, risfcpl_trc       , ptr(:,:,:,:,Krhs))
      IF ( kt == nit000+1) CALL trc_isf_cpl(Kmm, risfcpl_trc*0.5_wp, ptr(:,:,:,:,Krhs))
      !!
      !!----------------------
      ! ensure 0 trend due to unconservation of the ice shelf coupling
      IF ( ln_isfcpl_cons ) CALL trc_isf_cpl(Kmm, risfcpl_cons_trc, ptr(:,:,:,:,Krhs))
      !!
      IF ( ln_isfdebug ) THEN
         IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN                       ! Do only for the full domain
            DO jn = 1, jptra
               WRITE(charout,*) 'trc_isf: tr(:,:,:,:,Krhs )', trim(ctrcnm(jn))
               CALL debug( trim(charout), ptr(:,:,:,jn,Krhs))
            END DO
         ENDIF
      END IF
      !
      IF( ln_timing )   CALL timing_stop('trc_isf')
      !
   END SUBROUTINE trc_isf
   !
   SUBROUTINE trc_isf_cpl( Kmm, ptsc, ptsa )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_isf_cpl  ***
      !!
      !! *** Action :: Update pts(:,:,:,:,Krhs) with the ice shelf coupling trend
      !!
      !!----------------------------------------------------------------------
      REAL(dp), DIMENSION(jpi,jpj,jpk,jptra), INTENT(inout) :: ptsa
      !!----------------------------------------------------------------------
      INTEGER                               , INTENT(in   ) :: Kmm   ! ocean time level index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra), INTENT(in   ) :: ptsc
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj, jk, jn
      !!----------------------------------------------------------------------
      !
      DO_3D( 0, 0, 0, 0, 1, jpk )
         DO jn = 1, jptra
          ptsa(ji,jj,jk,jn) = ptsa(ji,jj,jk,jn) + ptsc(ji,jj,jk,jp_tem) * r1_e1e2t(ji,jj) / e3t(ji,jj,jk,Kmm)
         END DO
      END_3D
      !
   END SUBROUTINE trc_isf_cpl
   !
END MODULE trcisf
