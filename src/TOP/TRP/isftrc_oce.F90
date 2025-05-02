MODULE isftrc_oce
   !!======================================================================
   !!                       ***  MODULE  isf_oce  ***
   !! Ice shelves :  ice shelves variables defined in memory
   !!======================================================================
   !! History :  3.2  !  2011-02  (C.Harris  ) Original code isf cav
   !!            X.X  !  2006-02  (C. Wang   ) Original code bg03
   !!            3.4  !  2013-03  (P. Mathiot) Merging + parametrization
   !!            4.1  !  2019-09  (P. Mathiot) Split param/explicit ice shelf and re-organisation
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   isf          : define and allocate ice shelf variables
   !!----------------------------------------------------------------------
   USE par_kind
   USE par_oce       , ONLY: jpi, jpj, jpk
   USE in_out_manager
   USE lib_mpp       , ONLY: ctl_stop, mpp_sum      ! MPP library
   USE fldread        ! read input fields
   USE isf_oce       , ONLY: ln_isfcpl_cons
   USE par_trc       , ONLY: jptra, jpmaxtrc

   IMPLICIT NONE

   PRIVATE

   PUBLIC   isftrc_alloc_cpl, isftrc_dealloc_cpl
   !
   !! Jpalm -- 12-03-2025 -- need a real parameter constant for
   !!          the array dimension in the structure bellow
   !!          the allocation does not work otherwise
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   risfcpl_trc, risfcpl_cons_trc  !:
   !
   TYPE, PUBLIC :: isfconspt                          !! pt for Passive Tracers
      INTEGER                             ::   ii     ! i global
      INTEGER                             ::   jj     ! j global
      INTEGER                             ::   kk     ! k level
      REAL(wp), DIMENSION(jpmaxtrc)       ::   dpt    ! number of passive tracers
      REAL(wp)                            ::   lon    ! lon
      REAL(wp)                            ::   lat    ! lat
      INTEGER                             ::   ngb    ! 0/1 (valid location or not (ie on halo or no neighbourg))
   END TYPE

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcisf.F90 10536 2019-01-16 19:21:09Z mathiot $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE isftrc_alloc_cpl()
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_alloc_cpl  ***
      !!
      !! ** Purpose : allocate array use for the ice sheet coupling
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ierr, ialloc
      !!----------------------------------------------------------------------
      ierr = 0
      !
      ALLOCATE( risfcpl_trc(jpi,jpj,jpk,jptra) ,STAT=ialloc )
      ierr = ierr + ialloc
      !
      risfcpl_trc(:,:,:,:) = 0._wp
      !
      IF(lwp) WRITE(numout,*) ' isftrc_oce :  risfcpl_trc allocated'
      IF(lwp) WRITE(numout,*) '            just before allocating ln_isfcpl_cons'
      IF(lwp) WRITE(numout,*) ' '
      CALL flush(numout)
      !!
      IF ( ln_isfcpl_cons ) THEN
         ALLOCATE( risfcpl_cons_trc(jpi,jpj,jpk,jptra), STAT=ialloc )
         ierr = ierr + ialloc
         !
         risfcpl_cons_trc(:,:,:,:) = 0._wp
      ENDIF
      !
      IF(lwp) WRITE(numout,*) ' isftrc_oce :  risfcpl_cons_trc allocated'
      IF(lwp) WRITE(numout,*) ' '
      CALL flush(numout)
      !
      CALL mpp_sum ( 'isf', ierr )
      IF( ierr /= 0 )   CALL ctl_stop('STOP','isfcpl: failed to allocate arrays.')
      !
   END SUBROUTINE isftrc_alloc_cpl

   
   SUBROUTINE isftrc_dealloc_cpl()
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE isf_dealloc_cpl  ***
      !!
      !! ** Purpose : de-allocate useless public 3d array used for ice sheet coupling
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ierr, ialloc
      !!----------------------------------------------------------------------
      ierr = 0
      !
      DEALLOCATE( risfcpl_trc , STAT=ialloc )
      ierr = ierr + ialloc
      !
      CALL mpp_sum ( 'isf', ierr )
      IF( ierr /= 0 )   CALL ctl_stop('STOP','isfcpl: failed to deallocate arrays.')
      !
   END SUBROUTINE isftrc_dealloc_cpl

   !!======================================================================
END MODULE isftrc_oce
