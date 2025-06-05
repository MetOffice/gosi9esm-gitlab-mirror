MODULE isfcpl
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
#if defined key_qco
   USE domqco  , ONLY : dom_qco_zgr      ! vertical scale factor interpolation
#elif defined key_linssh
   !                                     ! fix in time coordinate
#else
   USE domvvl  , ONLY : dom_vvl_zgr      ! vertical scale factor interpolation
#endif
   USE domutl  , ONLY : dom_ngb          ! find the closest grid point from a given lon/lat position
   USE isf_oce        ! ice shelf variable
   USE isfutils, ONLY : debug
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lib_mpp , ONLY : mpp_sum, mpp_max ! mpp routine
   !
#if defined key_top
   USE par_trc,    ONLY : jptra, jpmaxtrc
   USE isftrc_oce
#endif
   !
   IMPLICIT NONE

   PRIVATE

   PUBLIC isfcpl_rst_write, isfcpl_init                    ! iceshelf restart read and write
   PUBLIC isfcpl_ssh, isfcpl_tr, isfcpl_vol, isfcpl_cons  ! iceshelf correction for ssh, tra, dyn and conservation

   TYPE isfcons
      INTEGER                   ::   ii     ! i global
      INTEGER                   ::   jj     ! j global
      INTEGER                   ::   kk     ! k level
      REAL(wp)                  ::   dvol   ! volume increment
      REAL(wp), DIMENSION(jpts) ::   dts    ! 1-heat increment; 2-salt increment
      REAL(wp)                  ::   lon    ! lon
      REAL(wp)                  ::   lat    ! lat
      INTEGER                   ::   ngb    ! 0/1 (valid location or not (ie on halo or no neighbourg))
   END TYPE
   !!---------------------------------------------------------------------
   INTEGER, PUBLIC :: id
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
   SUBROUTINE isfcpl_init(Kbb, Kmm, Kaa)
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE iscpl_init  ***
      !!
      !! ** Purpose : correct ocean state for new wet cell and horizontal divergence
      !!              correction for the dynamical adjustement
      !!
      !! ** Action : - compute ssh on new wet cell
      !!             - compute T/S on new wet cell
      !!             - compute horizontal divergence correction as a volume flux
      !!             - compute the T/S/vol correction increment to keep trend to 0
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: Kbb, Kmm, Kaa      ! ocean time level indices
      !!----------------------------------------------------------------------
      !
      ! start on an euler time step
      l_1st_euler = .TRUE.
      !
      ! allocation and initialisation to 0
      CALL isf_alloc_cpl()
      !
      ! check presence of variable needed for coupling
      ! iom_varid return 0 if not found
      id = 1
      id = id * iom_varid(numror, 'ssmask', ldstop = .false.)
      id = id * iom_varid(numror, 'tmask' , ldstop = .false.)
      id = id * iom_varid(numror, 'e3t_n' , ldstop = .false.)
      id = id * iom_varid(numror, 'e3u_n' , ldstop = .false.)
      id = id * iom_varid(numror, 'e3v_n' , ldstop = .false.)
      IF(lwp) WRITE(numout,*) ' isfcpl_init:', id
      IF (id == 0) THEN
         IF(lwp) WRITE(numout,*) ' isfcpl_init: restart variables for ice sheet coupling are missing, skip coupling for this leg '
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~'
         IF(lwp) WRITE(numout,*) ''
      ELSE
         ! extrapolation ssh
         CALL isfcpl_ssh(Kbb, Kmm, Kaa)
         !
         ! extrapolation tracer properties
         CALL isfcpl_tr(Kmm,'TRA',ts,jpts)
         !
         ! correction of the horizontal divergence and associated temp. and salt content flux
         ! Need to : - include in the cpl cons the risfcpl_vol/tsc contribution
         !           - decide how to manage thickness level change in conservation
         CALL isfcpl_vol(Kmm)
         !
         ! apply the 'conservation' method
         IF ( ln_isfcpl_cons ) CALL isfcpl_cons(Kmm,'TRA',ts,jpts)
         !
      END IF
      !
      ! mask velocity properly (mask used in restart not compatible with new mask)
      uu(:,:,:,Kmm) = uu(:,:,:,Kmm) * umask(:,:,:)
      vv(:,:,:,Kmm) = vv(:,:,:,Kmm) * vmask(:,:,:)
      !
      ! all before fields set to now values
      ts  (:,:,:,:,Kbb) = ts  (:,:,:,:,Kmm)
      uu   (:,:,:,Kbb)   = uu   (:,:,:,Kmm)
      vv   (:,:,:,Kbb)   = vv   (:,:,:,Kmm)
      ssh (:,:,Kbb)     = ssh (:,:,Kmm)
#if ! defined key_qco   &&   ! defined key_linssh
      e3t(:,:,:,Kbb)   = e3t(:,:,:,Kmm)
#endif
   END SUBROUTINE isfcpl_init


   SUBROUTINE isfcpl_rst_write( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE iscpl_rst_write  ***
      !!
      !! ** Purpose : write icesheet coupling variables in restart
      !!
      !!-------------------------- IN  --------------------------------------
      INTEGER, INTENT(in) :: kt
      INTEGER, INTENT(in) :: Kmm    ! ocean time level index
      !!----------------------------------------------------------------------
      INTEGER :: jk                               ! loop index
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ze3t, ze3u, ze3v, zgdepw  ! for qco substitution
      !!----------------------------------------------------------------------
      !
      DO jk = 1, jpk
         ze3t(:,:,jk) = e3t(:,:,jk,Kmm)
         ze3u(:,:,jk) = e3u(:,:,jk,Kmm)
         ze3v(:,:,jk) = e3v(:,:,jk,Kmm)
         !
         zgdepw(:,:,jk) = gdepw(:,:,jk,Kmm)
      END DO
      !
      CALL iom_rstput( kt, nitrst, numrow, 'tmask'  , tmask  )
      CALL iom_rstput( kt, nitrst, numrow, 'ssmask' , ssmask )
      CALL iom_rstput( kt, nitrst, numrow, 'e3t_n'  , ze3t   )
      CALL iom_rstput( kt, nitrst, numrow, 'e3u_n'  , ze3u   )
      CALL iom_rstput( kt, nitrst, numrow, 'e3v_n'  , ze3v   )
      CALL iom_rstput( kt, nitrst, numrow, 'gdepw_n', zgdepw )
      !
   END SUBROUTINE isfcpl_rst_write


   SUBROUTINE isfcpl_ssh(Kbb, Kmm, Kaa)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE iscpl_ssh  ***
      !!
      !! ** Purpose :   basic guess of ssh in new wet cell
      !!
      !! ** Method  :   basic extrapolation from neighbouring cells
      !!
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) :: Kbb, Kmm, Kaa    ! ocean time level indices
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj, jd, jk      !! loop index
      INTEGER :: jip1, jim1, jjp1, jjm1
      !!
      REAL(wp):: zsummsk
      REAL(wp), DIMENSION(jpi,jpj) :: zdssmask, zssmask0, zssmask_b, zssh
      !!----------------------------------------------------------------------
      !
      CALL iom_get( numror, jpdom_auto, 'ssmask'  , zssmask_b   ) ! need to extrapolate T/S

      ! compute new ssh if we open a full water column
      ! crude average of the closest neighbourgs (e1e2t not taking into account)
      !
      zssh(:,:)     = ssh(:,:,Kmm)
      zssmask0(:,:) = zssmask_b(:,:)
      !
      DO jd = 1, nn_drown
         !
         zdssmask(:,:) = ssmask(:,:) - zssmask0(:,:)
         DO_2D( 0, 0, 0, 0 )
            jip1=ji+1   ;   jim1=ji-1
            jjp1=jj+1   ;   jjm1=jj-1
            !
            zsummsk = zssmask0(jip1,jj) + zssmask0(jim1,jj) + zssmask0(ji,jjp1) + zssmask0(ji,jjm1)
            !
            IF (zdssmask(ji,jj) == 1._wp .AND. zsummsk /= 0._wp) THEN
               ssh(ji,jj,Kmm)=( zssh(jip1,jj)*zssmask0(jip1,jj)     &
                  &           + zssh(jim1,jj)*zssmask0(jim1,jj)     &
                  &           + zssh(ji,jjp1)*zssmask0(ji,jjp1)     &
                  &           + zssh(ji,jjm1)*zssmask0(ji,jjm1)) / zsummsk
               zssmask_b(ji,jj) = 1._wp
            ENDIF
         END_2D
         CALL lbc_lnk( 'isfcpl', ssh(:,:,Kmm), 'T', 1.0_dp)
         CALL lbc_lnk( 'isfcpl', zssmask_b(:,:), 'T', 1.0_wp )
         !
         zssh(:,:) = ssh(:,:,Kmm)
         zssmask0(:,:) = zssmask_b(:,:)
         !
         !
      END DO
      !
      ! update ssh(:,:,Kmm)
      ssh(:,:,Kmm) = zssh(:,:) * ssmask(:,:)
      !
      ssh(:,:,Kbb) = ssh(:,:,Kmm)
      !
      IF ( ln_isfdebug ) CALL debug('isfcpl_ssh: sshn',CASTSP(ssh(:,:,Kmm)))
      !
      ! recompute the vertical scale factor, depth and water thickness
      IF(lwp) write(numout,*) 'isfcpl_ssh : recompute scale factor from ssh (new wet cell,Kmm)'
      IF(lwp) write(numout,*) '~~~~~~~~~~~'
#if defined key_qco
      CALL dom_qco_zgr(Kbb, Kmm)
#elif defined key_linssh
      ! linear ssh : fix in time coord.
#else
      DO jk = 1, jpk
         e3t(:,:,jk,Kmm) = e3t_0(:,:,jk) * ( 1._wp + (ht_0(:,:) + ssh(:,:,Kmm)) * r1_ht_0(:,:) )
      END DO
      e3t(:,:,:,Kbb) = e3t(:,:,:,Kmm)
      CALL dom_vvl_zgr(Kbb, Kmm, Kaa)
#endif
      !
   END SUBROUTINE isfcpl_ssh


   SUBROUTINE isfcpl_tr(Kmm,cdtype,pt,kjpt)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE iscpl_tr  ***
      !!
      !! ** Purpose :   compute new tn, sn in case of evolving geometry of ice shelves
      !!
      !! ** Method  :   tn, sn : basic extrapolation from neighbouring cells
      !!
      !!----------------------------------------------------------------------
      INTEGER                                  , INTENT(in   ) ::   Kmm    ! ocean time level index
      CHARACTER(len=3)                         , INTENT(in   ) ::   cdtype ! =TRA or TRC (tracer indicator)
      INTEGER                                  , INTENT(in   ) ::   kjpt   ! number of tracers = jpts (=2 - TRA) or jptra (TRC)
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt,jpt), INTENT(inout) ::   pt     ! tracers and RHS of tracer equation
                                                                           ! pt is ts (TRA) or the tr array (TRC)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ztmask_b
      !REAL(wp), DIMENSION(:,:,:  ), INTENT(in ) :: pdepw_b                         !! depth w before
      !!
      INTEGER :: ji, jj, jk, jn, jd          !! loop index
      INTEGER :: jip1, jim1, jjp1, jjm1, jkp1, jkm1
      !!
      REAL(wp):: zsummsk
      REAL(wp):: zdz, zdzm1, zdzp1
      !!
      REAL(wp), DIMENSION(jpi,jpj)          :: zdmask
      REAL(wp), DIMENSION(jpi,jpj,jpk)      :: ztmask0, zwmaskn
      REAL(wp), DIMENSION(jpi,jpj,jpk)      :: ztmask1, zwmaskb, ztmp3d
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt) :: zpt0
      !!----------------------------------------------------------------------
      !
      ! Julpal --Check  
      IF(lwp) write(numout,*)'isfcpl_tr -- ', cdtype,' running on ', kjpt,' variables'
      IF(lwp) write(numout,*)''
      !
      CALL iom_get( numror, jpdom_auto, 'tmask'  , ztmask_b   ) ! need to extrapolate T/S
      !CALL iom_get( numror, jpdom_auto, 'wmask'  , zwmask_b  ) ! need to extrapolate T/S
      !CALL iom_get( numror, jpdom_auto, 'gdepw_n', zdepw_b(:,:,:) ) ! need to interpol vertical profile (vvl)
      !
      !
      ! compute new T/S (interpolation) if vvl only for common wet cell in before and after wmask
      !PM: Is this IF needed since change to VVL by default
      !bugged : to be corrected (PM)
      ! back up original t/s/mask
      !tsb (:,:,:,:) = ts(:,:,:,:,Kmm)
      !
     ! compute new T/S (interpolation) if vvl only for common wet cell in before and after wmask

!      IF (.NOT.ln_linssh) THEN
!         DO jk = 2,jpk-1
!            DO jj = 1,jpj
!               DO ji = 1,jpi
!                  IF (wmask(ji,jj,jk) * zwmaskb(ji,jj,jk) == 1._wp .AND. (tmask(ji,jj,1)==0._wp .OR. ztmask_b(ji,jj,1)==0._wp) ) THEN
!
!                     !compute weight
!                     zdzp1 = MAX(0._wp,pdepw_b(ji,jj,jk+1) - gdepw(ji,jj,jk+1,Kmm))
!                     zdzm1 = MAX(0._wp,gdepw(ji,jj,jk  ,Kmm) - pdepw_b(ji,jj,jk  ))
!                     zdz   = e3t(ji,jj,jk,Kmm) - zdzp1 - zdzm1 ! if isf : e3t = gdepw(ji,jj,jk+1,Kmm)- gdepw(ji,jj,jk,Kmm)
!
!                     IF (zdz .LT. 0._wp) THEN
!                        CALL ctl_stop( 'STOP', 'rst_iscpl : unable to compute the interpolation' )
!                     END IF
!
!                     ts(ji,jj,jk,jp_tem,Kmm) = ( zdzp1*ts(ji,jj,jk+1,jp_tem,Kbb) &
!                        &                   + zdz  *ts(ji,jj,jk  ,jp_tem,Kbb) &
!                        &                   + zdzm1*ts(ji,jj,jk-1,jp_tem,Kbb) )/e3t(ji,jj,jk,Kmm)
!
!                     ts(ji,jj,jk,jp_sal,Kmm) = ( zdzp1*ts(ji,jj,jk+1,jp_sal,Kbb) &
!                        &                   + zdz  *ts(ji,jj,jk  ,jp_sal,Kbb) &
!                        &                   + zdzm1*ts(ji,jj,jk-1,jp_sal,Kbb) )/e3t(ji,jj,jk,Kmm)
!
!                  END IF
!               END DO
!            END DO
!         END DO
!      END IF

      zpt0(:,:,:,:)  = pt(:,:,:,:,Kmm)
      ztmask0(:,:,:) = ztmask_b(:,:,:)
      ztmask1(:,:,:) = ztmask_b(:,:,:)
      !
      ! iterate the extrapolation processes nn_drown times
      DO jd = 1,nn_drown ! resolution dependent (OK for ISOMIP+ case)
         DO jk = 1,jpk-1
            !
            ! define new wet cell
            zdmask(:,:) = tmask(:,:,jk) - ztmask0(:,:,jk);
            !
            DO_2D( 0, 0, 0, 0 )
               jip1=ji+1; jim1=ji-1;
               jjp1=jj+1; jjm1=jj-1;
               !
               ! check if a wet neighbouring cell is present
               zsummsk = ztmask0(jip1,jj  ,jk) + ztmask0(jim1,jj  ,jk) &
                       + ztmask0(ji  ,jjp1,jk) + ztmask0(ji  ,jjm1,jk)
               !
               ! if neighbouring wet cell available at the same level
               IF ( zdmask(ji,jj) == 1._wp  .AND. zsummsk /= 0._wp ) THEN
                  !
                  ! horizontal basic extrapolation
                  DO jn = 1,kjpt
                    pt(ji,jj,jk,jn,Kmm)=( zpt0(jip1,jj  ,jk,jn) * ztmask0(jip1,jj  ,jk) &
                       &                + zpt0(jim1,jj  ,jk,jn) * ztmask0(jim1,jj  ,jk) &
                       &                + zpt0(ji  ,jjp1,jk,jn) * ztmask0(ji  ,jjp1,jk) &
                       &                + zpt0(ji  ,jjm1,jk,jn) * ztmask0(ji  ,jjm1,jk) ) / zsummsk
                  ENDDO
                  !
                  ! update mask for next pass
                  ztmask1(ji,jj,jk)=1
                  !
               ! in case no neighbouring wet cell available at the same level
               ! check if a wet cell is available below
               ELSEIF (zdmask(ji,jj) == 1._wp .AND. zsummsk == 0._wp) THEN
                  !
                  ! vertical extrapolation if horizontal extrapolation failed
                  jkm1=max(1,jk-1) ; jkp1=min(jpk,jk+1)
                  !
                  ! check if a wet neighbouring cell is present
                  zsummsk = ztmask0(ji,jj,jkm1) + ztmask0(ji,jj,jkp1)
                  IF (zdmask(ji,jj) == 1._wp .AND. zsummsk /= 0._wp ) THEN
                     DO jn = 1,kjpt
                        pt(ji,jj,jk,jn,Kmm)=( zpt0(ji,jj,jkp1,jn)*ztmask0(ji,jj,jkp1)     &
                            &               + zpt0(ji,jj,jkm1,jn)*ztmask0(ji,jj,jkm1)) / zsummsk
                     ENDDO
                     !
                     ! update mask for next pass
                     ztmask1(ji,jj,jk)=1._wp
                  END IF
               END IF
            END_2D
         END DO
         !
         CALL lbc_lnk( 'isfcpl', pt(:,:,:,:,Kmm), 'T', 1.0_dp)
         CALL lbc_lnk( 'isfcpl', ztmask1, 'T', 1.0_wp)
         !
         ! update temperature and salinity and mask
         zpt0(:,:,:,:)  = pt(:,:,:,:,Kmm)
         ztmask0(:,:,:) = ztmask1(:,:,:)
         !
         !
      END DO  ! nn_drown
      !
      ! mask new ts(:,:,:,:,Kmm) field
      DO jn = 1,kjpt
         pt(:,:,:,jn,Kmm) = zpt0(:,:,:,jn) * tmask(:,:,:)
      ENDDO
      !
      ! sanity check
      ! -----------------------------------------------------------------------------------------
      ! case we open a cell but no neigbour cells available to get an estimate of T and S
      IF( cdtype == 'TRA' ) THEN
         DO_3D( 0, 0, 0, 0, 1,jpk-1 )
            IF (tmask(ji,jj,jk) == 1._wp .AND. pt(ji,jj,jk,jp_sal,Kmm) == 0._wp) &
               &   CALL ctl_stop('STOP', 'failing to fill all new weet cell,     &
               &                          try increase nn_drown or activate XXXX &
               &                         in your domain cfg computation'         )
         END_3D
      ENDIF
      !
   END SUBROUTINE isfcpl_tr


   SUBROUTINE isfcpl_vol(Kmm)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE iscpl_vol  ***
      !!
      !! ** Purpose : compute the correction of the local divergence to apply
      !!              during the first time step after the coupling.
      !!
      !! ** Method  : - compute horizontal vol div. before/after coupling
      !!              - compute vertical input
      !!              - compute correction
      !!
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) :: Kmm    ! ocean time level index
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj, jk
      INTEGER :: ikb, ikt
      !!
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zqvolb, zqvoln  ! vol flux div.         before/after coupling
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ze3u_b, ze3v_b  ! vertical scale factor before/after coupling
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ztmask_b        ! mask                  before       coupling
      !!----------------------------------------------------------------------
      !
      CALL iom_get( numror, jpdom_auto, 'tmask'  , ztmask_b )
      CALL iom_get( numror, jpdom_auto, 'e3u_n'  , ze3u_b   )
      CALL iom_get( numror, jpdom_auto, 'e3v_n'  , ze3v_b   )
      !
      ! 1.0: compute horizontal volume flux divergence difference before-after coupling
      !
      DO jk = 1, jpk                                 ! Horizontal slab
         ! 1.1: get volume flux before coupling (>0 out)
         DO_2D( 0, 0, 0, 0 )
            zqvolb(ji,jj,jk) =    &
               &  (   e2u(ji  ,jj  ) * ze3u_b(ji  ,jj  ,jk) * uu(ji  ,jj  ,jk,Kmm)      &
               &    - e2u(ji-1,jj  ) * ze3u_b(ji-1,jj  ,jk) * uu(ji-1,jj  ,jk,Kmm)      &
               &    + e1v(ji  ,jj  ) * ze3v_b(ji  ,jj  ,jk) * vv(ji  ,jj  ,jk,Kmm)      &
               &    - e1v(ji  ,jj-1) * ze3v_b(ji  ,jj-1,jk) * vv(ji  ,jj-1,jk,Kmm)  )   &
               &                * ztmask_b(ji,jj,jk)
         END_2D
         !
         ! 1.2: get volume flux after coupling (>0 out)
         ! properly mask velocity
         ! (velocity are still mask with old mask at this stage)
         uu(:,:,jk,Kmm) = uu(:,:,jk,Kmm) * umask(:,:,jk)
         vv(:,:,jk,Kmm) = vv(:,:,jk,Kmm) * vmask(:,:,jk)
         ! compute volume flux divergence after coupling
         DO_2D( 0, 0, 0, 0 )
            zqvoln(ji,jj,jk) =   &
               &  (   e2u(ji  ,jj  ) * e3u(ji  ,jj  ,jk,Kmm) * uu(ji  ,jj  ,jk,Kmm)    &
               &    - e2u(ji-1,jj  ) * e3u(ji-1,jj  ,jk,Kmm) * uu(ji-1,jj  ,jk,Kmm)    &
               &    + e1v(ji  ,jj  ) * e3v(ji  ,jj  ,jk,Kmm) * vv(ji  ,jj  ,jk,Kmm)    &
               &    - e1v(ji  ,jj-1) * e3v(ji  ,jj-1,jk,Kmm) * vv(ji  ,jj-1,jk,Kmm)  ) &
               &               * tmask(ji,jj,jk)
            !
            ! 1.3: get 3d volume flux difference (before - after cpl) (>0 out)
            !      correction to add is _b - _n
            risfcpl_vol(ji,jj,jk) = zqvolb(ji,jj,jk) - zqvoln(ji,jj,jk)
         END_2D
      END DO
      !
      ! 2.0: include the contribution of the vertical velocity in the volume flux correction
      !
      DO_2D( 0, 0, 0, 0 )
         !
         ikt = mikt(ji,jj)
         IF ( ikt > 1 .AND. ssmask(ji,jj) == 1 ) THEN
            risfcpl_vol(ji,jj,ikt) = risfcpl_vol(ji,jj,ikt) + SUM(zqvolb(ji,jj,1:ikt-1))  ! test sign
         ENDIF
         !
      END_2D
      !
      CALL lbc_lnk( 'isfcpl', risfcpl_vol, 'T', 1.0_wp )
      !
      ! 3.0: set total correction (div, tr(:,:,:,:,Krhs), ssh)
      !
      ! 3.1: mask volume flux divergence correction
      risfcpl_vol(:,:,:) = risfcpl_vol(:,:,:) * tmask(:,:,:)
      !
      ! 3.2: get 3d tr(:,:,:,:,Krhs) increment to apply at the first time step
      ! temperature and salt content flux computed using local ts(:,:,:,:,Kmm)
      ! (very simple advection scheme)
      ! (>0 out)
      risfcpl_tsc(:,:,:,jp_tem) = -risfcpl_vol(:,:,:) * ts(:,:,:,jp_tem,Kmm)
      risfcpl_tsc(:,:,:,jp_sal) = -risfcpl_vol(:,:,:) * ts(:,:,:,jp_sal,Kmm)
      !
      ! 3.3: ssh correction (for dynspg_ts)
      risfcpl_ssh(:,:) = 0.0
      DO jk = 1,jpk
         risfcpl_ssh(:,:) = risfcpl_ssh(:,:) + risfcpl_vol(:,:,jk) * r1_e1e2t(:,:)
      END DO
      !
   END SUBROUTINE isfcpl_vol


   SUBROUTINE isfcpl_cons(Kmm,cdtype,pt,kjpt)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE iscpl_cons  ***
      !!
      !! ** Purpose :   compute the corrective increment in volume/salt/heat to put back the vol/heat/salt
      !!                removed or added during the coupling processes (wet or dry new cell)
      !!
      !! ** Method  :   - compare volume/heat/salt before and after
      !!                - look for the closest wet cells (share amoung neighbourgs if there are)
      !!                - build the correction increment to applied at each time step
      !!
      !!----------------------------------------------------------------------
      !
#if defined key_top
      USE trc,        ONLY : ctrcnm
      USE isftrc_cpl, ONLY : update_isfptr, get_correction_pt
#endif
      !!----------------------------------------------------------------------
      TYPE(isfcons)  , DIMENSION(:),ALLOCATABLE :: zisfpts ! list of point receiving a correction
#if defined key_top
      TYPE(isfconspt), DIMENSION(:),ALLOCATABLE :: zisfptr ! list of point receiving a correction
#endif
      !
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)             :: Kmm    ! ocean time level index
      CHARACTER(len=3), INTENT(in   ) :: cdtype ! =TRA or TRC (tracer indicator)
      INTEGER         , INTENT(in   ) :: kjpt   ! number of tracers = jpts (=2 - TRA) or jptra (TRC)
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt,jpt), INTENT(inout) ::   pt
      !!----------------------------------------------------------------------
      INTEGER  ::   ji   , jj  , jk  , jn, jproc      ! loop index
      INTEGER  ::   jip1 , jim1, jjp1, jjm1           ! dummy indices
      INTEGER  ::   iig  , ijg, ik                    ! dummy indices
      INTEGER  ::   jisf                              ! start, end and current position in the increment array
      INTEGER  ::   ingb, ifind                       ! 0/1 target found or need to be found
      INTEGER  ::   nisfl_area                        ! global number of cell concerned by the wet->dry case
      INTEGER, DIMENSION(jpnij) :: nisfl              ! local  number of cell concerned by the wet->dry case
      !
      REAL(wp) ::   z1_sum, z1_rdtiscpl
      REAL(wp) ::   zdvol, zratio                     ! vol increment
      REAL(wp), ALLOCATABLE, DIMENSION(:) :: zdpt     ! 1-heat increment; 2-salt (TRA)
                                                      ! Or pass trc (TRC) increment
      REAL(wp) ::   zlon , zlat                       ! target location
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ztmask_b    ! mask before
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ze3t_b      ! scale factor before
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ztmp        ! scale factor before
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: ztmpt  ! local array to read restart var
      !!----------------------------------------------------------------------

      !==============================================================================
      ! 1.0: initialisation
      !==============================================================================

      !
      ALLOCATE( ztmpt(jpi,jpj,jpk,kjpt))
      IF( cdtype == 'TRA' ) THEN
         ALLOCATE( zdpt(kjpt))
#if defined key_top
      ELSEIF( cdtype == 'TRC' ) THEN
         ALLOCATE( zdpt(jpmaxtrc)) !! needs jpmaxtrc to match the structure array dimension
      ENDIF
#endif
      ! init zdpt
      zdpt(:) = 0.0
      IF(lwp) WRITE(numout,*) ' isfcpl_cons -- zdpt allocated and initialized'
      CALL flush(numout)
      !
      ! get restart variable
      CALL iom_get( numror, jpdom_auto, 'tmask'  , ztmask_b(:,:,:) ) ! need to extrapolate T/S
      CALL iom_get( numror, jpdom_auto, 'e3t_n'  , ze3t_b(:,:,:)   )
      !

      ! compute run length
      nstp_iscpl  = nitend - nit000 + 1
      rdt_iscpl   = nstp_iscpl * rn_Dt
      z1_rdtiscpl = 1._wp / rdt_iscpl

      IF (lwp) WRITE(numout,*) '            nb of stp for cons  = ', nstp_iscpl
      IF (lwp) WRITE(numout,*) '            coupling time step  = ', rdt_iscpl

      !
      !! get to TRA or TRC specific calculations :
      IF( cdtype == 'TRA' ) THEN
         CALL iom_get( numror, jpdom_auto, 'tn'     , ztmpt(:,:,:,1)   ) !! T for jn=1
         CALL iom_get( numror, jpdom_auto, 'sn'     , ztmpt(:,:,:,2)   )
         !CALL iom_get( numror, jpdom_auto, 'tn'     , zt_b(:,:,:)     )
         !CALL iom_get( numror, jpdom_auto, 'sn'     , zs_b(:,:,:)     )

         ! initialisation correction
         risfcpl_cons_vol = 0.0
         risfcpl_cons_ssh = 0.0
         risfcpl_cons_tsc = 0.0

#if defined key_top
      ELSEIF( cdtype == 'TRC' ) THEN
         DO jn = 1,kjpt
            CALL iom_get( numrtr, jpdom_auto, 'TRN'//ctrcnm(jn), ztmpt(:,:,:,jn) )
         ENDDO
         risfcpl_cons_trc = 0.0
#endif
      ENDIF
      !
      IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' risfcpl_cons_tr initialized'
      CALL flush(numout)

      !==============================================================================
      ! 2.0: diagnose the heat, salt and volume input and compute the correction variable
      !      for case where we wet a cell or cell still wet (no change in cell status)
      !==============================================================================

      IF( cdtype == 'TRA' ) THEN
         DO jk = 1,jpk-1
            DO jj = Njs0,Nje0
               DO ji = Nis0,Nie0
                  ! volume diff
                  zdvol =   e3t(ji,jj,jk,Kmm) *  tmask  (ji,jj,jk)   &
                     &   - ze3t_b(ji,jj,jk    ) * ztmask_b(ji,jj,jk)
                  DO jn = 1,kjpt
                     ! 1: heat diff - 2: Sal diff
                     zdpt(jn) = pt (ji,jj,jk,jn,Kmm) *  e3t(ji,jj,jk,Kmm) *  tmask  (ji,jj,jk)   &
                           - ztmpt(ji,jj,jk,jn)        * ze3t_b(ji,jj,jk) * ztmask_b(ji,jj,jk)
                  ENDDO
                  ! volume, heat and salt differences in each cell (>0 means correction is an outward flux)
                  ! in addition to the geometry change unconservation, need to add the divergence correction as it is flux across the boundary
                  risfcpl_cons_vol(ji,jj,jk)        = (   zdvol * e1e2t(ji,jj) + risfcpl_vol(ji,jj,jk)        ) * z1_rdtiscpl
                  DO jn = 1,kjpt
                     risfcpl_cons_tsc(ji,jj,jk,jn) = ( - zdpt(jn) * e1e2t(ji,jj) + risfcpl_tsc(ji,jj,jk,jn) ) * z1_rdtiscpl
                  ENDDO
               END DO

            END DO
         END DO
#if defined key_top
     ELSEIF( cdtype == 'TRC' ) THEN
        !! need the vol update as well as done for T-S in isfcpl_vol() :
        !! now moved out into isfcpl_vol_pt
        ! DO jn = 1,kjpt
        !   risfcpl_trc(:,:,:,jn) = -risfcpl_vol(:,:,:) * pt(:,:,:,jn,Kmm)
        ! END DO
        !!
        DO jk = 1,jpk-1
            DO jj = Njs0,Nje0
               DO ji = Nis0,Nie0
                  DO jn = 1,kjpt
                     zdpt(jn) = pt (ji,jj,jk,jn,Kmm) *  e3t(ji,jj,jk,Kmm) *  tmask  (ji,jj,jk)   &
                           - ztmpt(ji,jj,jk,jn)        * ze3t_b(ji,jj,jk) * ztmask_b(ji,jj,jk)
                     risfcpl_cons_trc(ji,jj,jk,jn) = ( - zdpt(jn) * e1e2t(ji,jj) + risfcpl_trc(ji,jj,jk,jn) ) * z1_rdtiscpl
                  ENDDO
               END DO
            END DO
         END DO
#endif
      ENDIF
      !
      IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' risfcpl_cons_tr filled'
      CALL flush(numout)
      !
      !==============================================================================
      ! 3.0: diagnose the heat, salt and volume input and compute the correction variable
      !      for case where we close a cell
      !==============================================================================
      !
      ! compute the total number of point receiving a correction increment for each processor
      ! local
      nisfl(:)=0
      DO jk = 1,jpk-1
         DO jj = Njs0,Nje0
            DO ji = Nis0,Nie0
               jip1=MIN(ji+1,jpi) ; jim1=MAX(ji-1,1) ; jjp1=MIN(jj+1,jpj) ; jjm1=MAX(jj-1,1) ;
               IF ( tmask(ji,jj,jk) == 0._wp .AND. ztmask_b(ji,jj,jk) == 1._wp ) THEN
                  nisfl(narea) = nisfl(narea) + MAX(SUM(tmask(jim1:jip1,jjm1:jjp1,jk)),1._wp)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      !
      ! global
      CALL mpp_sum('isfcpl',nisfl  )
      !
      IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' isfcpl nisfl = ',nisfl
      CALL flush(numout)
      !
      IF( cdtype == 'TRA' ) THEN
         ! allocate list of point receiving correction
         ALLOCATE(zisfpts(nisfl(narea)))
         !
         IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' zisfpts allocated'
         CALL flush(numout)
         !!TYPE isfcons   (ii,jj,kk,dvol,dts(T,S),  lon,lat,ngb)
         !! re-use zdpt. was use for risfcpl_cons_tsc calc - but is not needed anymore and has the correct dim
         zdpt(:) = -HUGE(1.0)
         !! now init zisfpts
         zisfpts(:) = isfcons(0,0,0,-HUGE(1.0), zdpt, -HUGE(1.0), -HUGE(1.0), 0)
#if defined key_top
      ELSEIF( cdtype == 'TRC' ) THEN
         !
         IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' isfcpl nisfl = ',nisfl
         CALL flush(numout)
         !
         ALLOCATE(zisfptr(nisfl(narea)))
         !
         IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' zisfptr allocated'
         CALL flush(numout)
         !!TYPE isfconspt (ii,jj,kk,dts(jpmaxtrc),lon,lat,ngb)
         !! re-use zdpt. was use for risfcpl_cons_tsc calc - but is not needed anymore and has the correct dim
         zdpt(:) = -HUGE(1.0)
         !! now init zisfptr
         zisfptr(:) = isfconspt(0,0,0, zdpt, -HUGE(1.0), -HUGE(1.0), 0)
#endif
      ENDIF
      !
      IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' zisfpts/r initialized'
      CALL flush(numout)
      !
      !
      !
      ! start computing the correction and fill zisfpts
      ! local
      jisf = 0
      DO jk = 1,jpk-1
         DO jj = Njs0,Nje0
            DO ji = Nis0,Nie0
               IF ( tmask(ji,jj,jk) == 0._wp .AND. ztmask_b(ji,jj,jk) == 1._wp ) THEN

                  jip1=MIN(ji+1,jpi) ; jim1=MAX(ji-1,1) ; jjp1=MIN(jj+1,jpj) ; jjm1=MAX(jj-1,1) ;

                  zdvol = risfcpl_cons_vol(ji,jj,jk       )
                  IF( cdtype == 'TRA' ) THEN
                     DO jn = 1,kjpt
                        zdpt(jn) = risfcpl_cons_tsc(ji,jj,jk,jn)
                     ENDDO
#if defined key_top
                  ELSEIF( cdtype == 'TRC' ) THEN
                     DO jn = 1,kjpt
                        zdpt(jn) = risfcpl_cons_trc(ji,jj,jk,jn)
                     ENDDO
#endif
                  ENDIF
                  !
                  IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' correction calculations -- zdpt filled from risfcpl_cons_tr'
                  IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' now will call update isfpt...'
                  CALL flush(numout)
                  !!
                  IF( cdtype == 'TRA' ) THEN
                     IF ( SUM( tmask(jim1:jip1,jjm1:jjp1,jk) ) > 0._wp ) THEN
                        ! spread correction amoung neighbouring wet cells (horizontal direction first)
                        ! as it is a crude correction corner and lateral cell have the same weight
                        !
                        z1_sum =  1._wp / SUM( tmask(jim1:jip1,jjm1:jjp1,jk) )
                        !
                        ! lateral cells
                        IF (tmask(jip1,jj  ,jk) == 1) CALL update_isfpts(zisfpts, jisf, jip1, jj  , jk, zdvol, zdpt, z1_sum)
                        IF (tmask(jim1,jj  ,jk) == 1) CALL update_isfpts(zisfpts, jisf, jim1, jj  , jk, zdvol, zdpt, z1_sum)
                        IF (tmask(ji  ,jjp1,jk) == 1) CALL update_isfpts(zisfpts, jisf, ji  , jjp1, jk, zdvol, zdpt, z1_sum)
                        IF (tmask(ji  ,jjm1,jk) == 1) CALL update_isfpts(zisfpts, jisf, ji  , jjm1, jk, zdvol, zdpt, z1_sum)
                        !
                        ! corner  cells
                        IF (tmask(jip1,jjm1,jk) == 1) CALL update_isfpts(zisfpts, jisf, jip1, jjm1, jk, zdvol, zdpt, z1_sum)
                        IF (tmask(jim1,jjm1,jk) == 1) CALL update_isfpts(zisfpts, jisf, jim1, jjm1, jk, zdvol, zdpt, z1_sum)
                        IF (tmask(jim1,jjp1,jk) == 1) CALL update_isfpts(zisfpts, jisf, jim1, jjp1, jk, zdvol, zdpt, z1_sum)
                        IF (tmask(jip1,jjp1,jk) == 1) CALL update_isfpts(zisfpts, jisf, jip1, jjp1, jk, zdvol, zdpt, z1_sum)
                        !
                     ELSE IF ( tmask(ji,jj,jk+1) == 1._wp ) THEN
                        ! spread correction amoung neighbouring wet cells (vertical direction)
                        CALL update_isfpts(zisfpts, jisf, ji  , jj  , jk+1, zdvol, zdpt, 1.0_wp, 0)
                     ELSE
                        ! need to find where to put correction in later on
                        CALL update_isfpts(zisfpts, jisf, ji  , jj  , jk  , zdvol, zdpt, 1.0_wp, 1)
                     END IF
#if defined key_top
                  ELSEIF( cdtype == 'TRC' ) THEN
                      IF ( SUM( tmask(jim1:jip1,jjm1:jjp1,jk) ) > 0._wp ) THEN
                        ! spread correction amoung neighbouring wet cells (horizontal direction first)
                        ! as it is a crude correction corner and lateral cell have the same weight
                        !
                        z1_sum =  1._wp / SUM( tmask(jim1:jip1,jjm1:jjp1,jk) )
                        !
                        ! lateral cells
                        IF (tmask(jip1,jj  ,jk) == 1) CALL update_isfptr(zisfptr, jisf, jip1, jj  , jk, zdpt, z1_sum)
                        IF (tmask(jim1,jj  ,jk) == 1) CALL update_isfptr(zisfptr, jisf, jim1, jj  , jk, zdpt, z1_sum)
                        IF (tmask(ji  ,jjp1,jk) == 1) CALL update_isfptr(zisfptr, jisf, ji  , jjp1, jk, zdpt, z1_sum)
                        IF (tmask(ji  ,jjm1,jk) == 1) CALL update_isfptr(zisfptr, jisf, ji  , jjm1, jk, zdpt, z1_sum)
                        !
                        ! corner  cells
                        IF (tmask(jip1,jjm1,jk) == 1) CALL update_isfptr(zisfptr, jisf, jip1, jjm1, jk, zdpt, z1_sum)
                        IF (tmask(jim1,jjm1,jk) == 1) CALL update_isfptr(zisfptr, jisf, jim1, jjm1, jk, zdpt, z1_sum)
                        IF (tmask(jim1,jjp1,jk) == 1) CALL update_isfptr(zisfptr, jisf, jim1, jjp1, jk, zdpt, z1_sum)
                        IF (tmask(jip1,jjp1,jk) == 1) CALL update_isfptr(zisfptr, jisf, jip1, jjp1, jk, zdpt, z1_sum)
                        !
                     ELSE IF ( tmask(ji,jj,jk+1) == 1._wp ) THEN
                        ! spread correction amoung neighbouring wet cells (vertical direction)
                        CALL update_isfptr(zisfptr, jisf, ji  , jj  , jk+1, zdpt, 1.0_wp, 0)
                     ELSE
                        ! need to find where to put correction in later on
                        CALL update_isfptr(zisfptr, jisf, ji  , jj  , jk  , zdpt, 1.0_wp, 1)
                     END IF
#endif
                  ENDIF
               END IF
            END DO
         END DO
      END DO
      !
      IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' update isfpt DONE'
      IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' now will get_corrections'
      CALL flush(numout)
      !
      ! share data among all processes because for some point we need to find the closest wet point (could be on other process)
      DO jproc=1,jpnij
         !
         ! share total number of isf point treated for proc jproc
         IF (jproc==narea) THEN
            nisfl_area=nisfl(jproc)
         ELSE
            nisfl_area=0
         END IF
         CALL mpp_max('isfcpl',nisfl_area)
         !
         DO jisf = 1,nisfl_area
            IF( cdtype == 'TRA' ) THEN
            !
               IF (jproc==narea) THEN
                  ! indices (conversion to global indices and sharing)
                  iig = zisfpts(jisf)%ii       ; ijg = zisfpts(jisf)%jj       ; ik = zisfpts(jisf)%kk
                  !
                  ! data
                  zdvol = zisfpts(jisf)%dvol
                  DO jn = 1,kjpt
                     zdpt(jn) = zisfpts(jisf)%dts(jn)
                  ENDDO
                  !
                  ! location
                  zlat = zisfpts(jisf)%lat     ; zlon = zisfpts(jisf)%lon
                  !
                  ! find flag
                  ingb = zisfpts(jisf)%ngb
               ELSE
                  iig  =0   ; ijg  =0   ; ik   =0
                  zdvol=-HUGE(1.0) ; zdpt(:)=-HUGE(1.0)
                  zlat =-HUGE(1.0) ; zlon =-HUGE(1.0)
                  ingb = 0
               END IF
               !
               ! share data (need synchronisation of data as get_correction call a global com)
               CALL mpp_max('isfcpl',iig)   ; CALL mpp_max('isfcpl',ijg)   ; CALL mpp_max('isfcpl',ik)
               CALL mpp_max('isfcpl',zdvol) ;
               DO jn = 1,kjpt
                  CALL mpp_max('isfcpl',zdpt(jn))
               ENDDO
               CALL mpp_max('isfcpl',zlat)  ; CALL mpp_max('isfcpl',zlon)
               CALL mpp_max('isfcpl',ingb)
               !
               ! fill the 3d correction array
               CALL get_correction(iig, ijg, ik, zlon, zlat, zdvol, zdpt, ingb)
#if defined key_top
            ELSEIF( cdtype == 'TRC' ) THEN
            !
               IF (jproc==narea) THEN
                  ! indices (conversion to global indices and sharing)
                  iig = zisfptr(jisf)%ii       ; ijg = zisfptr(jisf)%jj       ; ik = zisfptr(jisf)%kk
                  !
                  ! data
                  DO jn = 1,kjpt
                     zdpt(jn) = zisfptr(jisf)%dpt(jn)
                  ENDDO
                  !
                  ! location
                  zlat = zisfptr(jisf)%lat     ; zlon = zisfptr(jisf)%lon
                  !
                  ! find flag
                  ingb = zisfptr(jisf)%ngb
               ELSE
                  iig  =0   ; ijg  =0   ; ik   =0
                  zdpt(:)=-HUGE(1.0)
                  zlat =-HUGE(1.0) ; zlon =-HUGE(1.0)
                  ingb = 0
               END IF
               !
               ! share data (need synchronisation of data as get_correction call a global com)
               CALL mpp_max('isfcpl',iig)   ; CALL mpp_max('isfcpl',ijg)   ; CALL mpp_max('isfcpl',ik)
               DO jn = 1,kjpt
                  CALL mpp_max('isfcpl',zdpt(jn))
               ENDDO
               CALL mpp_max('isfcpl',zlat)  ; CALL mpp_max('isfcpl',zlon)
               CALL mpp_max('isfcpl',ingb)
               !
               ! fill the 3d correction arraydptdptdpt
               CALL get_correction_pt(iig, ijg, ik, zlon, zlat, zdpt, ingb)
#endif
            ENDIF
         END DO
      END DO
      !
      IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' get_corrections DONE'
      IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' finalize risfcpl_cons_trc '
      CALL flush(numout)
      !
      !==============================================================================
      ! 4.0: finalisation and compute ssh equivalent of the volume correction
      !==============================================================================
      !
      IF( cdtype == 'TRA' ) THEN
         ! mask
         risfcpl_cons_vol(:,:,:       ) = risfcpl_cons_vol(:,:,:       ) * tmask(:,:,:)
         DO jn = 1,kjpt
            risfcpl_cons_tsc(:,:,:,jn) = risfcpl_cons_tsc(:,:,:,jn) * tmask(:,:,:)
         ENDDO
         !
         ! add lbclnk
         CALL lbc_lnk( 'isfcpl', risfcpl_cons_tsc(:,:,:,jp_tem), 'T', 1.0_wp, risfcpl_cons_tsc(:,:,:,jp_sal), 'T', 1.0_wp, &
            &                    risfcpl_cons_vol(:,:,:)       , 'T', 1.0_wp)
         !
         ! ssh correction (for dynspg_ts)
         DO jk = 1,jpk
            risfcpl_cons_ssh(:,:) = risfcpl_cons_ssh(:,:) + risfcpl_cons_vol(:,:,jk)
         END DO
         risfcpl_cons_ssh(:,:) = risfcpl_cons_ssh(:,:) * r1_e1e2t(:,:)
         !
#if defined key_top
      ELSEIF( cdtype == 'TRC' ) THEN
         ! mask
         DO jn = 1,kjpt
            risfcpl_cons_trc(:,:,:,jn) = risfcpl_cons_trc(:,:,:,jn) * tmask(:,:,:)
         ENDDO
         !
         ! add lbclnk
         CALL lbc_lnk( 'isfcpl', risfcpl_cons_trc(:,:,:,:), 'T', 1.0_wp)
         !
#endif
      ENDIF
      !
      IF(lwp) WRITE(numout,*) ' isfcpl_cons -- ',cdtype,' CALL finished -- all done'
      CALL flush(numout)
      !
   END SUBROUTINE isfcpl_cons
   !
   SUBROUTINE update_isfpts(sisfpts, kpts, ki, kj, kk, pdvol, pdts, pratio, kfind)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE update_isfpts  ***
      !!
      !! ** Purpose : if a cell become dry, we need to put the corrective increment elsewhere
      !!
      !! ** Action  : update the list of point
      !!
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      TYPE(isfcons), DIMENSION(:), INTENT(inout) :: sisfpts
      INTEGER                    , INTENT(inout) :: kpts
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in   )           :: ki, kj, kk      !    target location (kfind=0)
      !                                                                   ! or source location (kfind=1)
      INTEGER                , INTENT(in   ), OPTIONAL :: kfind           ! 0  target cell already found
      !                                                                   ! 1  target to be determined
      REAL(wp), DIMENSION(jpts) , INTENT(in   )        :: pdts            ! sal/tem increment
      REAL(wp)                  , INTENT(in   )        :: pdvol, pratio   ! vol increment
      !                                                                   ! and ratio in case increment span over multiple cells.
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
      sisfpts(kpts) = isfcons(mig(ki), mjg(kj), kk, pratio * pdvol, pratio * pdts, glamt(ki,kj), gphit(ki,kj), ifind )
      !
   END SUBROUTINE update_isfpts
   !
   SUBROUTINE get_correction( ki, kj, kk, plon, plat, pvolinc, ptslinc, kfind)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE get_correction  ***
      !!
      !! ** Action : - Find the closest valid cell if needed (wet and not on the halo)
      !!             - Scale the correction depending of pratio (case where multiple wet neighbourgs)
      !!             - Fill the correction array
      !!
      !!----------------------------------------------------------------------
      INTEGER                 , INTENT(in) :: ki, kj, kk, kfind ! target point indices
      REAL(wp)                , INTENT(in) :: plon, plat        ! target point lon/lat
      REAL(wp)                , INTENT(in) :: pvolinc           ! correction increment for vol
      REAL(wp),DIMENSION(jpts), INTENT(in) :: ptslinc           ! correction increment for temp/salt
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
            risfcpl_cons_vol(ji,jj,kk)        =  risfcpl_cons_vol(ji,jj,kk       ) + pvolinc
            DO jn = 1,2
               risfcpl_cons_tsc(ji,jj,kk,jn) =  risfcpl_cons_tsc(ji,jj,kk,jn) + ptslinc(jn)
            END DO
         END DO
      END DO

   END SUBROUTINE get_correction

END MODULE isfcpl
