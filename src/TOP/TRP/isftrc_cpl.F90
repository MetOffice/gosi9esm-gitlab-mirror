MODULE isftrc_cpl
   !!======================================================================
   !!                       ***  MODULE  isfcpl  ***
   !!
   !! iceshelf coupling module : module managing the coupling between NEMO and an ice sheet model
   !!
   !!======================================================================
   !! History :  4.1  !  2019-07  (P. Mathiot) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   isfrst        : read/write iceshelf variables in/from restart
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   !
   USE domutl  , ONLY : dom_ngb          ! find the closest grid point from a given lon/lat position
   USE isf_oce        ! ice shelf variable
   USE isftrc_oce     ! trc shelf variable 
   USE isfutils, ONLY : debug
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lib_mpp , ONLY : mpp_sum, mpp_max ! mpp routine
   !
   USE par_trc , ONLY : jptra, jpmaxtrc
   !

   IMPLICIT NONE

   PRIVATE
   !
   PUBLIC update_isfptr,get_correction_pt  ! iceshelf restart read and write
   !
   !!---------------------------------------------------------------------
   !
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "single_precision_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcisf.F90 10536 2019-01-16 19:21:09Z mathiot $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE update_isfptr(sisfpts, kpts, ki, kj, kk, pdtr, pratio, kfind)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE update_isfpts  ***
      !!
      !! ** Purpose : if a cell become dry, we need to put the corrective increment elsewhere
      !!
      !! ** Action  : update the list of point
      !!
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      TYPE(isfconspt), DIMENSION(:), INTENT(inout) :: sisfpts
      INTEGER                      , INTENT(inout) :: kpts
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   )           :: ki, kj, kk !    target location (kfind=0)
      !                                                                 ! or source location (kfind=1)
      INTEGER                   , INTENT(in   ), OPTIONAL :: kfind      ! 0  target cell already found
      !                                                                 ! 1  target to be determined
      REAL(wp), DIMENSION(jpmaxtrc), INTENT(in   )        :: pdtr       ! passive tracer increment
      REAL(wp)                  , INTENT(in   )           :: pratio     ! and ratio in case increment span over multiple cells.
      !!----------------------------------------------------------------------
      INTEGER :: ifind, jn
      !!----------------------------------------------------------------------
      !
      ! increment position
      kpts = kpts + 1
      !
      ! define if we need to look for closest valid wet cell (no neighbours or neighbourg on halo)
      IF ( PRESENT(kfind) ) THEN
         ifind = kfind
      ELSE
         ifind = ( 1 - tmask_i(ki,kj) ) * tmask(ki,kj,kk)
      END IF
      !
      ! update isfpts structure
      sisfpts(kpts) = isfconspt(mig(ki), mjg(kj), kk, pratio * pdtr, glamt(ki,kj), gphit(ki,kj), ifind )
      !
   END SUBROUTINE update_isfptr
   !
   SUBROUTINE get_correction_pt( ki, kj, kk, plon, plat, ptrlinc, kfind)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE get_correction  ***
      !!
      !! ** Action : - Find the closest valid cell if needed (wet and not on the halo)
      !!             - Scale the correction depending of pratio (case where multiple wet neighbourgs)
      !!             - Fill the correction array
      !!
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in) :: ki, kj, kk, kfind        ! target point indices
      REAL(wp)                  , INTENT(in) :: plon, plat               ! target point lon/lat
      REAL(wp), DIMENSION(jpmaxtrc), INTENT(in) :: ptrlinc                  ! correction increment for vol/temp/salt
      !!----------------------------------------------------------------------
      INTEGER :: jj, ji, jn, iig, ijg
      !!----------------------------------------------------------------------
      !
      ! define global indice of correction location
      iig = ki ; ijg = kj
      IF ( kfind == 1 ) CALL dom_ngb( plon, plat, iig, ijg,'T', kk)
      !
      ! fill the correction array
      DO jj = mj0(ijg),mj1(ijg)
         DO ji = mi0(iig),mi1(iig)
            ! correct the vol_flx and corresponding heat/salt flx in the closest cell
            DO jn = 1,jptra
               risfcpl_cons_trc(ji,jj,kk,jn) =  risfcpl_cons_trc(ji,jj,kk,jn) + ptrlinc(jn)
            ENDDO
         END DO
      END DO

   END SUBROUTINE get_correction_pt
   !
END MODULE isftrc_cpl
