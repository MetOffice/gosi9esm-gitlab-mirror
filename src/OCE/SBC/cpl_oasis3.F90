MODULE cpl_oasis3
   !!======================================================================
   !!                    ***  MODULE cpl_oasis  ***
   !! Coupled O/A : coupled ocean-atmosphere case using OASIS3-MCT
   !!=====================================================================
   !! History :  1.0  !  2004-06  (R. Redler, NEC Laboratories Europe, Germany) Original code
   !!             -   !  2004-11  (R. Redler, NEC Laboratories Europe; N. Keenlyside, W. Park, IFM-GEOMAR, Germany) revision
   !!             -   !  2004-11  (V. Gayler, MPI M&D) Grid writing
   !!            2.0  !  2005-08  (R. Redler, W. Park) frld initialization, paral(2) revision
   !!             -   !  2005-09  (R. Redler) extended to allow for communication over root only
   !!             -   !  2006-01  (W. Park) modification of physical part
   !!             -   !  2006-02  (R. Redler, W. Park) buffer array fix for root exchange
   !!            3.4  !  2011-11  (C. Harris) Changes to allow mutiple category fields
   !!            3.6  !  2014-11  (S. Masson) OASIS3-MCT
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_oasis3'                    coupled Ocean/Atmosphere via OASIS3-MCT
   !!----------------------------------------------------------------------
   !!   cpl_init     : initialization of coupled mode communication
   !!   cpl_define   : definition of grid and fields
   !!   cpl_snd      : snd out fields in coupled mode
   !!   cpl_rcv      : receive fields in coupled mode
   !!   cpl_finalize : finalize the coupled mode communication
   !!----------------------------------------------------------------------
#if defined key_oasis3
   USE mod_oasis                    ! OASIS3-MCT module
#endif
#if defined key_xios
   USE xios                         ! I/O server
#endif
   USE par_oce                      ! ocean parameters

   USE cpl_rnf_1d, ONLY: nn_cpl_river   ! Variables used in 1D river outflow 

   USE dom_oce                      ! ocean space and time domain
   USE in_out_manager               ! I/O manager
   USE lbclnk                       ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp
#if defined key_agrif || ! defined key_mpi_off
   USE MPI
#endif

   IMPLICIT NONE
   PRIVATE
#if ! defined key_oasis3
   ! Dummy interface to oasis_get if not using oasis
   INTERFACE oasis_get
      MODULE PROCEDURE oasis_get_1d, oasis_get_2d
   END INTERFACE
#endif 

   PUBLIC   cpl_init
   PUBLIC   cpl_define
   PUBLIC   cpl_snd
   PUBLIC   cpl_rcv
   PUBLIC   cpl_rcv_1d
   PUBLIC   cpl_freq
   PUBLIC   cpl_finalize

   INTEGER, PARAMETER         ::   localRoot  = 0
   INTEGER, PUBLIC            ::   OASIS_Rcv  = 1    !: return code if received field
   INTEGER, PUBLIC            ::   OASIS_idle = 0    !: return code if nothing done by oasis
   INTEGER                    ::   ncomp_id          ! id returned by oasis_init_comp
   INTEGER                    ::   nerror            ! return error code
#if ! defined key_oasis3
   ! OASIS Variables not used. defined only for compilation purpose
   INTEGER                    ::   OASIS_Out         = -1
   INTEGER                    ::   OASIS_REAL        = -1
   INTEGER                    ::   OASIS_Ok          = -1
   INTEGER                    ::   OASIS_In          = -1
   INTEGER                    ::   OASIS_Sent        = -1
   INTEGER                    ::   OASIS_SentOut     = -1
   INTEGER                    ::   OASIS_ToRest      = -1
   INTEGER                    ::   OASIS_ToRestOut   = -1
   INTEGER                    ::   OASIS_Recvd       = -1
   INTEGER                    ::   OASIS_RecvOut     = -1
   INTEGER                    ::   OASIS_FromRest    = -1
   INTEGER                    ::   OASIS_FromRestOut = -1
#endif

   INTEGER                    ::   nrcv         ! total number of fields received
   INTEGER                    ::   nsnd         ! total number of fields sent
   INTEGER                    ::   ncplmodel    ! Maximum number of models to/from which NEMO is potentialy sending/receiving data

  INTEGER, PUBLIC, PARAMETER ::   nmaxfld=100  ! Maximum number of coupling fields

   INTEGER, PUBLIC, PARAMETER ::   nmaxcat=5    ! Maximum number of coupling fields
   INTEGER, PUBLIC, PARAMETER ::   nmaxcpl=5    ! Maximum number of coupling fields
   LOGICAL, PARAMETER         ::   ltmp_wapatch = .FALSE.   ! patch to restore wraparound rows in cpl_send, cpl_rcv, cpl_define 
   LOGICAL, PARAMETER         ::   ltmp_landproc = .TRUE.   ! patch to restrict coupled area to non halo cells

   TYPE, PUBLIC ::   FLD_CPL               !: Type for coupling field information
      LOGICAL               ::   laction   ! To be coupled or not
      CHARACTER(len = 8)    ::   clname    ! Name of the coupling field
      CHARACTER(len = 1)    ::   clgrid    ! Grid type
      REAL(wp)              ::   nsgn      ! Control of the sign change
      INTEGER, DIMENSION(nmaxcat,nmaxcpl) ::   nid   ! Id of the field (no more than 9 categories and 9 extrena models)
      INTEGER               ::   nct       ! Number of categories in field
      INTEGER               ::   ncplmodel ! Maximum number of models to/from which this variable may be sent/received
      INTEGER               ::   dimensions ! Number of dimensions of coupling field 
   END TYPE FLD_CPL

   TYPE(FLD_CPL), DIMENSION(nmaxfld), PUBLIC ::   srcv, ssnd   !: Coupling fields

   REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   exfld   ! Temporary buffer for receiving
   REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   exfld_ext  ! Temporary buffer for receiving with wrap points

   INTEGER :: ishape_ext(4)  ! shape of 2D arrays passed to PSMILe extended for wrap points in weights data
   INTEGER :: ishape(4)      ! shape of 2D arrays passed to PSMILe 

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: cpl_oasis3.F90 14434 2021-02-11 08:20:52Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE cpl_init( cd_modname, kl_comm )
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_init  ***
      !!
      !! ** Purpose :   Initialize coupled mode communication for ocean
      !!    exchange between AGCM, OGCM and COUPLER. (OASIS3 software)
      !!
      !! ** Method  :   OASIS3 MPI communication
      !!--------------------------------------------------------------------
      CHARACTER(len = *), INTENT(in   ) ::   cd_modname   ! model name as set in namcouple file
      INTEGER           , INTENT(  out) ::   kl_comm      ! local communicator of the model
      INTEGER                           ::   error
      !!--------------------------------------------------------------------

      ! WARNING: No write in numout in this routine
      !============================================

      !------------------------------------------------------------------
      ! 1st Initialize the OASIS system for the application
      !------------------------------------------------------------------
      CALL oasis_init_comp ( ncomp_id, TRIM(cd_modname), nerror )
      IF( nerror /= OASIS_Ok ) &
         CALL oasis_abort (ncomp_id, 'cpl_init', 'Failure in oasis_init_comp')

      !------------------------------------------------------------------
      ! 3rd Get an MPI communicator for OCE local communication
      !------------------------------------------------------------------

      CALL oasis_get_localcomm ( kl_comm, nerror )
      IF( nerror /= OASIS_Ok ) &
         CALL oasis_abort (ncomp_id, 'cpl_init','Failure in oasis_get_localcomm' )
      !

   END SUBROUTINE cpl_init


   SUBROUTINE cpl_define( krcv, ksnd, kcplmodel )
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_define  ***
      !!
      !! ** Purpose :   Define grid and field information for ocean
      !!    exchange between AGCM, OGCM and COUPLER. (OASIS3 software)
      !!
      !! ** Method  :   OASIS3 MPI communication
      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) ::   krcv, ksnd     ! Number of received and sent coupling fields
      INTEGER, INTENT(in) ::   kcplmodel      ! Maximum number of models to/from which NEMO is potentialy sending/receiving data
      !
      INTEGER :: id_part_0d     ! Partition for 0d fields
      INTEGER :: id_part_rnf_1d ! Partition for 1d river outflow fields
      INTEGER :: id_part_2d     ! Partition for 2d fields
      INTEGER :: id_part_2d_ext ! Partition for 2d fields extended for old (pre vn4.2) style remapping weights!
      INTEGER :: id_part_temp   ! Temperary partition used to choose either 0d or 1d partitions
      INTEGER :: paral(5)       ! OASIS3 box partition

      INTEGER :: paral_ext(5)       ! OASIS3 box partition extended
      INTEGER :: ishape0d1d(2)  ! Shape of 0D or 1D arrays passed to PSMILe.
      INTEGER :: var_nodims(2)  ! Number of coupling field dimensions.
                                ! var_nodims(1) is redundant from OASIS3-MCT vn4.0 onwards
                                ! but retained for backward compatibility.
                                ! var_nodims(2) is the number of fields in a bundle
                                ! or 1 for unbundled fields (bundles are not yet catered for
                                ! in NEMO hence we default to 1).        
      INTEGER :: ji,jc,jm       ! local loop indicees
      LOGICAL :: llenddef       ! should we call xios_oasis_enddef and oasis_enddef?
      CHARACTER(LEN=64) :: zclname
      CHARACTER(LEN=2) :: cli2
      INTEGER :: i_offset       ! Used in calculating offset for extended partition.
      !!--------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'cpl_define : initialization in coupled ocean/atmosphere case'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*)

      ncplmodel = kcplmodel
      IF( kcplmodel > nmaxcpl ) THEN
         CALL oasis_abort ( ncomp_id, 'cpl_define', 'ncplmodel is larger than nmaxcpl, increase nmaxcpl')   ;   RETURN
      ENDIF

      nrcv = krcv
      IF( nrcv > nmaxfld ) THEN
         CALL oasis_abort ( ncomp_id, 'cpl_define', 'nrcv is larger than nmaxfld, increase nmaxfld')   ;   RETURN
      ENDIF

      nsnd = ksnd
      IF( nsnd > nmaxfld ) THEN
         CALL oasis_abort ( ncomp_id, 'cpl_define', 'nsnd is larger than nmaxfld, increase nmaxfld')   ;   RETURN
      ENDIF
      !
      ! ... Define the shape for the area that excludes the halo as we don't want them to be "seen" by oasis
      !
      ishape(1) = 1
      ishape(2) = Ni_0
      ishape(3) = 1
      ishape(4) = Nj_0

      ishape0d1d(1) = 0
      ishape0d1d(2) = 0       !

      ! ... Allocate memory for data exchange
      !
      ALLOCATE(exfld(Ni_0, Nj_0), stat = nerror)                 ! allocate full domain (without halos)
      IF( nerror > 0 ) THEN
         CALL oasis_abort ( ncomp_id, 'cpl_define', 'Failure in allocating exfld')   ;   RETURN
      ENDIF
      !
      ! -----------------------------------------------------------------
      ! ... Define the partition, excluding halos as we don't want them to be "seen" by oasis
      ! -----------------------------------------------------------------

      paral(1) = 2                                      ! box partitioning
      paral(2) = Ni0glo * mjg0(nn_hls) + mig0(nn_hls)   ! NEMO lower left corner global offset, without halos
      paral(3) = Ni_0                                   ! local extent in i, excluding halos
      paral(4) = Nj_0                                   ! local extent in j, excluding halos
      paral(5) = Ni0glo                                 ! global extent in x, excluding halos

      IF( sn_cfctl%l_oasout ) THEN
         WRITE(numout,*) ' multiexchg: paral (1:5)', paral
         WRITE(numout,*) ' multiexchg: Ni_0, Nj_0 =', Ni_0, Nj_0
         WRITE(numout,*) ' multiexchg: Nis0, Nie0, nimpp =', Nis0, Nie0, nimpp
         WRITE(numout,*) ' multiexchg: Njs0, Nje0, njmpp =', Njs0, Nje0, njmpp
      ENDIF


      ! We still set up the new vn4.2 style box partition for reference, though it doesn't actually get used,
      ! we can easily swap back to it if we ever manage to successfully generate vn4.2 compatible weights, or introduce 
      ! RTL controls to distinguish between onl and new style weights.  

      CALL oasis_def_partition ( id_part_2d, paral, nerror, Ni0glo*Nj0glo )   ! global number of points, excluding halos

      ! RSRH Set up 2D box partition for compatibility with existing weights files
      ! so we don't have to generate and manage multiple sets of weights purely because of 
      ! the changes to nemo 4.2+ code!

      ! This is just a hack for global cyclic models for the time being
      Ni0glo_ext = jpiglo
      Nj0glo_ext = Nj0glo +1 ! We can't use jpjglo here because for some reason at 4.2 this is bigger
                             ! than at 4.0.... e.g. for ORCA1 it is 333 when it should only be 332!

      ! RSRH extended shapes for old style dimensioning. Allows backwards compatibility with existing weights files, 
      ! which the new code DOES NOT, causing headaches not only for users but also for management of weights files. 
      ishape_ext(:) = ishape(:)
      IF (mig(1) == 1 .OR. mig(jpi)==jpiglo) THEN
         ! Extra columns in PEs dealing with wrap points
         ishape_ext(2) = ishape_ext(2) + 1
      ENDIF

      ! Workout any extra offset in the i dimension
      IF (mig(1) == 1 ) THEN
         i_offset = mig0(nn_hls)
      ELSE
         i_offset = mig(nn_hls)
      ENDIF
       
      ALLOCATE(exfld_ext(ishape_ext(2), ishape_ext(4)), stat = nerror)        ! allocate full domain (with wrap pts)
      IF( nerror > 0 ) THEN
         CALL oasis_abort ( ncomp_id, 'cpl_define', 'Failure in allocating exfld_ext')   ;   RETURN
      ENDIF


      ! Now we have the appropriate dimensions, we can set up the partition array for the old-style extended grid
      paral_ext(1) = 2                                      ! box partitioning
      paral_ext(2) = (Ni0glo_ext * mjg0(nn_hls)) + i_offset ! NEMO lower left corner global offset, with wrap pts
      paral_ext(3) = Ni_0_ext                               ! local extent in i, including halos
      paral_ext(4) = Nj_0_ext                               ! local extent in j, including halos
      paral_ext(5) = Ni0glo_ext                             ! global extent in x, including halos

      IF( sn_cfctl%l_oasout ) THEN
         WRITE(numout,*) ' multiexchg: paral_ext (1:5)', paral_ext, jpiglo, jpjglo, Ni0glo_ext, Nj0glo_ext
         WRITE(numout,*) ' multiexchg: Ni_0_ext, Nj_0_ext i_offset =', Ni_0_ext, Nj_0_ext, i_offset
         WRITE(numout,*) ' multiexchg: Nis0_ext, Nie0_ext =', Nis0_ext, Nie0_ext
         WRITE(numout,*) ' multiexchg: Njs0_ext, Nje0_ext =', Njs0_ext, Nje0_ext
      ENDIF

      ! Define our extended grid
      CALL oasis_def_partition ( id_part_2d_ext, paral_ext, nerror, Ni0glo_ext*Nj0glo_ext ) 
  
      ! OK so now we should have a usable 2d partition for fields defined WITH redundant points. 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! A special partition is needed for 0D fields
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      paral(1) = 0                                       ! serial partitioning
      paral(2) = 0   
      IF ( mpprank == 0) THEN
         paral(3) = 1                   ! Size of array to couple (scalar)
      ELSE
         paral(3) = 0                   ! Dummy size for PE's not involved
      END IF
      paral(4) = 0
      paral(5) = 0
             
      CALL oasis_def_partition ( id_part_0d, paral, nerror, 1 )


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Another special partition is needed for 1D river routing fields
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      paral(1) = 0                                       ! serial partitioning
      paral(2) = 0   
      IF ( mpprank == 0) THEN
         paral(3) = nn_cpl_river                   ! Size of array to couple (vector)
      ELSE
         paral(3) = 0                   ! Dummy size for PE's not involved
      END IF
      paral(4) = 0
      paral(5) = 0


      CALL oasis_def_partition ( id_part_rnf_1d, paral, nerror, nn_cpl_river )

      !
      ! ... Announce send variables.
      !
      ssnd(:)%ncplmodel = kcplmodel
      !
      DO ji = 1, ksnd
         IF( ssnd(ji)%laction ) THEN

            IF( ssnd(ji)%nct > nmaxcat ) THEN
               CALL oasis_abort ( ncomp_id, 'cpl_define', 'Number of categories of '//   &
                  &              TRIM(ssnd(ji)%clname)//' is larger than nmaxcat, increase nmaxcat' )
               RETURN
            ENDIF

            DO jc = 1, ssnd(ji)%nct
               DO jm = 1, kcplmodel

                  IF( ssnd(ji)%nct .GT. 1 ) THEN
                     WRITE(cli2,'(i2.2)') jc
                     zclname = TRIM(ssnd(ji)%clname)//'_cat'//cli2
                  ELSE
                     zclname = ssnd(ji)%clname
                  ENDIF
                  IF( kcplmodel  > 1 ) THEN
                     WRITE(cli2,'(i2.2)') jm
                     zclname = 'model'//cli2//'_'//TRIM(zclname)
                  ENDIF
#if defined key_agrif
                  IF( agrif_fixed() /= 0 ) THEN
                     zclname=TRIM(Agrif_CFixed())//'_'//TRIM(zclname)
                  ENDIF
#endif
                  IF( sn_cfctl%l_oasout ) WRITE(numout,*) "Define", ji, jc, jm, " "//TRIM(zclname), " for ", OASIS_Out

                  flush(numout)


                  IF( sn_cfctl%l_oasout ) WRITE(numout,*) "Define", ji, jc, jm, " "//TRIM(zclname), " for ", OASIS_Out
                     CALL oasis_def_var (ssnd(ji)%nid(jc,jm), zclname, id_part_2d   , (/ 2, 1 /),   &
                     &                OASIS_Out           , ishape , OASIS_REAL, nerror )

                  IF( nerror /= OASIS_Ok ) THEN
                     WRITE(numout,*) 'Failed to define transient ', ji, jc, jm, " "//TRIM(zclname)
                     CALL oasis_abort ( ssnd(ji)%nid(jc,jm), 'cpl_define', 'Failure in oasis_def_var' )
                  ENDIF
                  IF( sn_cfctl%l_oasout .AND. ssnd(ji)%nid(jc,jm) /= -1 ) WRITE(numout,*) "put variable defined in the namcouple"
                  IF( sn_cfctl%l_oasout .AND. ssnd(ji)%nid(jc,jm) == -1 ) WRITE(numout,*) "put variable NOT defined in the namcouple"




               END DO
            END DO
         ENDIF
      END DO
      !
      ! ... Announce received variables.
      !
      srcv(:)%ncplmodel = kcplmodel
      !
      DO ji = 1, krcv
         IF( srcv(ji)%laction ) THEN

            IF( srcv(ji)%nct > nmaxcat ) THEN
               CALL oasis_abort ( ncomp_id, 'cpl_define', 'Number of categories of '//   &
                  &              TRIM(srcv(ji)%clname)//' is larger than nmaxcat, increase nmaxcat' )
               RETURN
            ENDIF

            DO jc = 1, srcv(ji)%nct
               DO jm = 1, kcplmodel

                  IF( srcv(ji)%nct .GT. 1 ) THEN
                     WRITE(cli2,'(i2.2)') jc
                     zclname = TRIM(srcv(ji)%clname)//'_cat'//cli2
                  ELSE
                     zclname = srcv(ji)%clname
                  ENDIF
                  IF( kcplmodel  > 1 ) THEN
                     WRITE(cli2,'(i2.2)') jm
                     zclname = 'model'//cli2//'_'//TRIM(zclname)
                  ENDIF
#if defined key_agrif
                  IF( agrif_fixed() /= 0 ) THEN
                     zclname=TRIM(Agrif_CFixed())//'_'//TRIM(zclname)
                  ENDIF
#endif

                  IF( sn_cfctl%l_oasout ) WRITE(numout,*) "Define", ji, jc, jm, " "//TRIM(zclname), " for ", OASIS_In

                  flush(numout)

                  ! Define 0D (Greenland or Antarctic ice mass) or 1D (river outflow) coupling fields
                  IF (srcv(ji)%dimensions <= 1) THEN
                     var_nodims(1) = 1
                     var_nodims(2) = 1 ! Modify this value to cater for bundled fields.
                     IF (mpprank == 0) THEN
	
                       IF (srcv(ji)%dimensions == 0) THEN

                          ! If 0D then set temporary variables to 0D components
                          id_part_temp = id_part_0d
                          ishape0d1d(2) = 1
                       ELSE

                          ! If 1D then set temporary variables to river outflow components
                          id_part_temp = id_part_rnf_1d
                          ishape0d1d(2)= nn_cpl_river

                       END IF

                       CALL oasis_def_var (srcv(ji)%nid(jc,jm), zclname, id_part_temp   , var_nodims,   &
                                        OASIS_In           , ishape0d1d(1:2) , OASIS_REAL, nerror )
                    ELSE
                       ! Dummy call to keep OASIS3-MCT happy.
                       CALL oasis_def_var (srcv(ji)%nid(jc,jm), zclname, id_part_0d   , var_nodims,   &
                                        OASIS_In           , ishape0d1d(1:2) , OASIS_REAL, nerror )
                    END IF
                  ELSE
                    ! It's a "normal" 2D (or pseudo 3D) coupling field.
                    ! ... Set the field dimension and bundle count
                    var_nodims(1) = 2
                    var_nodims(2) = 1 ! Modify this value to cater for bundled fields.

                    CALL oasis_def_var (srcv(ji)%nid(jc,jm), zclname, id_part_2d   , var_nodims,   &
                                        OASIS_In           , ishape , OASIS_REAL, nerror )

                  ENDIF 



                  IF( nerror /= OASIS_Ok ) THEN
                     WRITE(numout,*) 'Failed to define transient ', ji, jc, jm, " "//TRIM(zclname)
                     CALL oasis_abort ( srcv(ji)%nid(jc,jm), 'cpl_define', 'Failure in oasis_def_var' )
                  ENDIF
                  IF( sn_cfctl%l_oasout .AND. srcv(ji)%nid(jc,jm) /= -1 ) WRITE(numout,*) "get variable defined in the namcouple"
                  IF( sn_cfctl%l_oasout .AND. srcv(ji)%nid(jc,jm) == -1 ) WRITE(numout,*) "get variable NOT defined in the namcouple"


               END DO
            END DO
         ENDIF
      END DO

      !------------------------------------------------------------------
      ! End of definition phase
      !------------------------------------------------------------------
      !
#if defined key_agrif
      IF( Agrif_Root() ) THEN   ! Warning: Agrif_Nb_Fine_Grids not yet defined -> must use Agrif_Root_Only()
         llenddef = Agrif_Root_Only()   ! true of no nested grid
      ELSE                      ! Is it the last nested grid?
         llenddef = agrif_fixed() == Agrif_Nb_Fine_Grids()
      ENDIF
#else
      llenddef = .TRUE.
#endif
      IF( llenddef ) THEN
#if defined key_xios
         CALL xios_oasis_enddef()   ! see "Joint_usage_OASIS3-MCT_XIOS.pdf" on XIOS wiki page
#endif
         CALL oasis_enddef(nerror)
         IF( nerror /= OASIS_Ok )   CALL oasis_abort ( ncomp_id, 'cpl_define', 'Failure in oasis_enddef')
      ENDIF
      !
   END SUBROUTINE cpl_define


   SUBROUTINE cpl_snd( kid, kstep, pdata, kinfo )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_snd  ***
      !!
      !! ** Purpose : - At each coupling time-step,this routine sends fields
      !!      like sst or ice cover to the coupler or remote application.
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   kid       ! variable index in the array
      INTEGER                   , INTENT(  out) ::   kinfo     ! OASIS3 info argument
      INTEGER                   , INTENT(in   ) ::   kstep     ! ocean time-step in seconds
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdata
      !!
      INTEGER                                   ::   jc,jm     ! local loop index
      !!--------------------------------------------------------------------
      !

      ! snd data to OASIS3
      !
      DO jc = 1, ssnd(kid)%nct
         DO jm = 1, ssnd(kid)%ncplmodel

            IF( ssnd(kid)%nid(jc,jm) /= -1 ) THEN   ! exclude halos from data sent to oasis

               ! The field is "put" directly, using appropriate start/end indexing - i.e. we don't
               ! copy it to an intermediate buffer. 
               CALL oasis_put ( ssnd(kid)%nid(jc,jm), kstep, pdata(Nis0:Nie0,Njs0:Nje0,jc), kinfo )

               IF ( sn_cfctl%l_oasout ) THEN
                  IF ( kinfo == OASIS_Sent     .OR. kinfo == OASIS_ToRest .OR.   &
                     & kinfo == OASIS_SentOut  .OR. kinfo == OASIS_ToRestOut ) THEN
                     WRITE(numout,*) '****************'
                     WRITE(numout,*) 'oasis_put: Outgoing ', ssnd(kid)%clname
                     WRITE(numout,*) 'oasis_put: ivarid ', ssnd(kid)%nid(jc,jm)
                     WRITE(numout,*) 'oasis_put:  kstep ', kstep
                     WRITE(numout,*) 'oasis_put:   info ', kinfo
                     WRITE(numout,*) '     - Minimum value is ', MINVAL(pdata(Nis0:Nie0,Njs0:Nje0,jc))
                     WRITE(numout,*) '     - Maximum value is ', MAXVAL(pdata(Nis0:Nie0,Njs0:Nje0,jc))
                     WRITE(numout,*) '     -     Sum value is ',    SUM(pdata(Nis0:Nie0,Njs0:Nje0,jc))
                     WRITE(numout,*) '****************'
                     CALL FLUSH(numout)
                  ENDIF
               ENDIF

            ENDIF

         ENDDO
      ENDDO

      !
    END SUBROUTINE cpl_snd


   SUBROUTINE cpl_rcv( kid, kstep, pdata, pmask, kinfo )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_rcv  ***
      !!
      !! ** Purpose : - At each coupling time-step,this routine receives fields
      !!      like stresses and fluxes from the coupler or remote application.
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   kid       ! variable index in the array
      INTEGER                   , INTENT(in   ) ::   kstep     ! ocean time-step in seconds
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pdata     ! IN to keep the value if nothing is done
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pmask     ! coupling mask
      INTEGER                   , INTENT(  out) ::   kinfo     ! OASIS3 info argument
      !!
      INTEGER                                   ::   jc,jm     ! local loop index
      LOGICAL                                   ::   llaction, ll_1st, lrcv

      !!--------------------------------------------------------------------
      !
      ! receive local data from OASIS3 on every process
      !
      kinfo = OASIS_idle
      lrcv=.FALSE.
      !

      DO jc = 1, srcv(kid)%nct
         ll_1st = .TRUE.

         DO jm = 1, srcv(kid)%ncplmodel

            IF( srcv(kid)%nid(jc,jm) /= -1 ) THEN

               CALL oasis_get ( srcv(kid)%nid(jc,jm), kstep, exfld, kinfo )

               llaction =  kinfo == OASIS_Recvd   .OR. kinfo == OASIS_FromRest .OR.   &
                  &        kinfo == OASIS_RecvOut .OR. kinfo == OASIS_FromRestOut

               IF ( sn_cfctl%l_oasout )   &
                  &  WRITE(numout,*) "llaction, kinfo, kstep, ivarid: " , llaction, kinfo, kstep, srcv(kid)%nid(jc,jm)

               IF( llaction ) THEN   ! data received from oasis do not include halos
                                     ! but DO still cater for wrap columns when using pre vn4.2 compatible remapping weights. 

                  kinfo = OASIS_Rcv
                  lrcv=.TRUE. 
                  IF( ll_1st ) THEN
                     pdata(Nis0:Nie0,Njs0:Nje0,jc) =   exfld(:,:) * pmask(Nis0:Nie0,Njs0:Nje0,jm)

                     ll_1st = .FALSE.
                  ELSE

                     pdata(Nis0:Nie0,Njs0:Nje0,jc) = pdata(Nis0:Nie0,Njs0:Nje0,jc)   &
                        &                                + exfld(:,:) * pmask(Nis0:Nie0,Njs0:Nje0,jm)
                  ENDIF

                  IF ( sn_cfctl%l_oasout ) THEN
                     WRITE(numout,*) '****************'
                     WRITE(numout,*) 'oasis_get: Incoming ', srcv(kid)%clname
                     WRITE(numout,*) 'oasis_get: ivarid '  , srcv(kid)%nid(jc,jm)
                     WRITE(numout,*) 'oasis_get:   kstep', kstep
                     WRITE(numout,*) 'oasis_get:   info ', kinfo
                     WRITE(numout,*) '     - Minimum value is ', MINVAL(pdata(Nis0:Nie0,Njs0:Nje0,jc))
                     WRITE(numout,*) '     - Maximum value is ', MAXVAL(pdata(Nis0:Nie0,Njs0:Nje0,jc))
                     WRITE(numout,*) '     -     Sum value is ',    SUM(pdata(Nis0:Nie0,Njs0:Nje0,jc))
                     WRITE(numout,*) '****************'
                     CALL FLUSH(numout)
                  ENDIF

               ENDIF

            ENDIF

         ENDDO

      ENDDO

      ! RSRH I've changed this since:
      ! 1) it seems multi cat fields may not be properly updated in halos when called on a per 
      !    category basis(?) 
      ! 2) it's more efficient to have a single call (and simpler to understand) to update ALL 
      !    categories at the same time!
      !--- Call lbc_lnk to populate halos of received fields.
      IF (lrcv) then
         CALL lbc_lnk( 'cpl_oasis3', pdata(:,:,:), srcv(kid)%clgrid, srcv(kid)%nsgn )
      endif
      !
   END SUBROUTINE cpl_rcv

  SUBROUTINE cpl_rcv_1d( kid, kstep, pdata, nitems, kinfo )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_rcv_1d  ***
      !!
      !! ** Purpose : - A special version of cpl_rcv to deal exclusively with
      !! receipt of 0D or 1D fields.
      !! The fields are recieved into a 1D array buffer which is simply a
      !! dynamically sized sized array (which may be of size 1)
      !! of 0 dimensional fields. This allows us to pass miltiple 0D
      !! fields via a single put/get operation. 
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in   ) ::   nitems      ! Number of 0D items to recieve
                                               ! during this get operation. i.e.
                                               ! The size of the 1D array in which
                                               ! 0D items are passed.   
      INTEGER , INTENT(in   ) ::   kid         ! ID index of the incoming
                                               ! data. 
      INTEGER , INTENT(in   ) ::   kstep       ! ocean time-step in seconds
      REAL(wp), INTENT(inout) ::   pdata(1:nitems) ! The original value(s), 
                                                   ! unchanged if nothing is
                                                   ! received
      INTEGER , INTENT(  out) ::   kinfo       ! OASIS3 info argument
      !!
      REAL(wp) ::   recvfld(1:nitems)          ! Local receive field buffer
      INTEGER  ::   jc,jm     ! local loop index
      INTEGER  ::   ierr
      LOGICAL  ::   llaction
      INTEGER  ::   MPI_WORKING_PRECISION
      INTEGER  ::   number_to_print
      !!--------------------------------------------------------------------
      !
      ! receive local data from OASIS3 on every process
      !
      kinfo = OASIS_idle
      !
      ! 0D and 1D fields won't have categories or any other form of "pseudo level"
      ! so we only cater for a single set of values and thus don't bother
      ! with a loop over the jc index
      jc = 1

      DO jm = 1, srcv(kid)%ncplmodel

         IF( srcv(kid)%nid(jc,jm) /= -1 ) THEN

            IF ( ( srcv(kid)%dimensions <= 1) .AND. (mpprank == 0) ) THEN
                ! Since there is no concept of data decomposition for zero
                ! dimension fields, they must only be exchanged through the master PE,
                ! unlike "normal" 2D field cases where every PE is involved.
 
                CALL oasis_get ( srcv(kid)%nid(jc,jm), kstep, recvfld, kinfo )   
                
                llaction =  kinfo == OASIS_Recvd   .OR. kinfo == OASIS_FromRest .OR.   &
                            kinfo == OASIS_RecvOut .OR. kinfo == OASIS_FromRestOut
                
                IF ( sn_cfctl%l_oasout ) WRITE(numout,*) "llaction, kinfo, kstep, ivarid: " , &
                                      llaction, kinfo, kstep, srcv(kid)%nid(jc,jm)
                
                IF ( llaction ) THEN
	                 
                   kinfo = OASIS_Rcv
	           pdata(1:nitems) = recvfld(1:nitems)
            
                   IF( sn_cfctl%l_oasout ) THEN     
	                     number_to_print = 10
	                     IF ( nitems < number_to_print ) number_to_print = nitems
	                     WRITE(numout,*) '****************'
	                     WRITE(numout,*) 'oasis_get: Incoming ', srcv(kid)%clname
	                     WRITE(numout,*) 'oasis_get: ivarid '  , srcv(kid)%nid(jc,jm)
	                     WRITE(numout,*) 'oasis_get:   kstep', kstep
	                     WRITE(numout,*) 'oasis_get:   info ', kinfo
	                     WRITE(numout,*) '     - Minimum Value is ', MINVAL(pdata(:))
	                     WRITE(numout,*) '     - Maximum value is ', MAXVAL(pdata(:))
	                     WRITE(numout,*) '     - Start of data is ', pdata(1:number_to_print)
	                     WRITE(numout,*) '****************'
                             CALL FLUSH(numout)
                  ENDIF
	                 
               ENDIF
            ENDIF   
          ENDIF
	           
       ENDDO

#if defined key_mpi_off

       CALL oasis_abort( ncomp_id, "cpl_rcv_1d", "Unable to use mpi_bcast with key_mpi_off in your list of NEMO keys." )

#else
       ! Set the precision that we want to broadcast using MPI_BCAST
       SELECT CASE( wp )
       CASE( sp )
         MPI_WORKING_PRECISION = MPI_REAL                ! Single precision
       CASE( dp )
         MPI_WORKING_PRECISION = MPI_DOUBLE_PRECISION    ! Double precision
       CASE default
         CALL oasis_abort( ncomp_id, "cpl_rcv_1d", "Could not find precision for coupling 0d or 1d field" )
       END SELECT

       ! We have to broadcast (potentially) received values from PE 0 to all
       ! the others. If no new data has been received we're just
       ! broadcasting the existing values but there's no more efficient way
       ! to deal with that w/o NEMO adopting a UM-style test mechanism
       ! to determine active put/get timesteps.
       CALL mpi_bcast( pdata, nitems, MPI_WORKING_PRECISION, localRoot, mpi_comm_oce, ierr )

#endif
	
        !
END SUBROUTINE cpl_rcv_1d


   INTEGER FUNCTION cpl_freq( cdfieldname )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_freq  ***
      !!
      !! ** Purpose : - send back the coupling frequency for a particular field
      !!----------------------------------------------------------------------
      CHARACTER(len = *), INTENT(in) ::   cdfieldname    ! field name as set in namcouple file
      !!
      INTEGER               :: id
      INTEGER               :: info
      INTEGER, DIMENSION(1) :: itmp
      INTEGER               :: ji,jm     ! local loop index
      INTEGER               :: mop
      !!----------------------------------------------------------------------
      cpl_freq = 0   ! defaut definition
      id = -1        ! defaut definition
      !
      DO ji = 1, nsnd
         IF(ssnd(ji)%laction ) THEN
            DO jm = 1, ncplmodel
               IF( ssnd(ji)%nid(1,jm) /= -1 ) THEN
                  IF( TRIM(cdfieldname) == TRIM(ssnd(ji)%clname) ) THEN
                     id = ssnd(ji)%nid(1,jm)
                     mop = OASIS_Out
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      DO ji = 1, nrcv
         IF(srcv(ji)%laction ) THEN
            DO jm = 1, ncplmodel
               IF( srcv(ji)%nid(1,jm) /= -1 ) THEN
                  IF( TRIM(cdfieldname) == TRIM(srcv(ji)%clname) ) THEN
                     id = srcv(ji)%nid(1,jm)
                     mop = OASIS_In
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      !
      IF( id /= -1 ) THEN
#if ! defined key_oa3mct_v1v2
         CALL oasis_get_freqs(id, mop, 1, itmp, info)
#else
         CALL oasis_get_freqs(id,      1, itmp, info)
#endif
         cpl_freq = itmp(1)
      ENDIF
      !
   END FUNCTION cpl_freq


   SUBROUTINE cpl_finalize
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_finalize  ***
      !!
      !! ** Purpose : - Finalizes the coupling. If MPI_init has not been
      !!      called explicitly before cpl_init it will also close
      !!      MPI communication.
      !!----------------------------------------------------------------------
      !
      DEALLOCATE( exfld )
      DEALLOCATE( exfld_ext )
      IF(nstop == 0) THEN
         CALL oasis_terminate( nerror )
      ELSE
         CALL oasis_abort( ncomp_id, "cpl_finalize", "NEMO ABORT STOP" )
      ENDIF
      !
   END SUBROUTINE cpl_finalize

#if ! defined key_oasis3

   !!----------------------------------------------------------------------
   !!   No OASIS Library          OASIS3 Dummy module...
   !!----------------------------------------------------------------------

   SUBROUTINE oasis_init_comp(k1,cd1,k2)
      CHARACTER(*), INTENT(in   ) ::  cd1
      INTEGER     , INTENT(  out) ::  k1,k2
      k1 = -1 ; k2 = -1
      WRITE(numout,*) 'oasis_init_comp: Error you sould not be there...', cd1
   END SUBROUTINE oasis_init_comp

   SUBROUTINE oasis_abort(k1,cd1,cd2)
      INTEGER     , INTENT(in   ) ::  k1
      CHARACTER(*), INTENT(in   ) ::  cd1,cd2
      WRITE(numout,*) 'oasis_abort: Error you sould not be there...', cd1, cd2
   END SUBROUTINE oasis_abort

   SUBROUTINE oasis_get_localcomm(k1,k2)
      INTEGER     , INTENT(  out) ::  k1,k2
      k1 = -1 ; k2 = -1
      WRITE(numout,*) 'oasis_get_localcomm: Error you sould not be there...'
   END SUBROUTINE oasis_get_localcomm

   SUBROUTINE oasis_def_partition(k1,k2,k3,k4)
      INTEGER     , INTENT(  out) ::  k1,k3
      INTEGER     , INTENT(in   ) ::  k2(5)
      INTEGER     , INTENT(in   ) ::  k4
      k1 = k2(1) ; k3 = k2(5)+k4
      WRITE(numout,*) 'oasis_def_partition: Error you sould not be there...'
   END SUBROUTINE oasis_def_partition

   SUBROUTINE oasis_def_var(k1,cd1,k2,k3,k4,k5,k6,k7)
      CHARACTER(*), INTENT(in   ) ::  cd1
      INTEGER     , INTENT(in   ) ::  k2,k3(2),k4,k5(*),k6
      INTEGER     , INTENT(  out) ::  k1,k7
      k1 = -1 ; k7 = -1
      WRITE(numout,*) 'oasis_def_var: Error you sould not be there...', cd1
   END SUBROUTINE oasis_def_var

   SUBROUTINE oasis_enddef(k1)
      INTEGER     , INTENT(  out) ::  k1
      k1 = -1
      WRITE(numout,*) 'oasis_enddef: Error you sould not be there...'
   END SUBROUTINE oasis_enddef

   SUBROUTINE oasis_put(k1,k2,p1,k3)
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::  p1
      INTEGER                 , INTENT(in   ) ::  k1,k2
      INTEGER                 , INTENT(  out) ::  k3
      k3 = -1
      WRITE(numout,*) 'oasis_put: Error you sould not be there...'
   END SUBROUTINE oasis_put

   SUBROUTINE oasis_get_2d(k1,k2,p1,k3)
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::  p1
      INTEGER                 , INTENT(in   ) ::  k1,k2
      INTEGER                 , INTENT(  out) ::  k3
      p1(1,1) = -1. ; k3 = -1
      WRITE(numout,*) 'oasis_get_2d: Error you sould not be there...'
   END SUBROUTINE oasis_get_2d

   SUBROUTINE oasis_get_1d(k1,k2,p1,k3)
      REAL(wp), DIMENSION(:)  , INTENT(  out) ::  p1
      INTEGER                 , INTENT(in   ) ::  k1,k2
      INTEGER                 , INTENT(  out) ::  k3
      p1(1) = -1. ; k3 = -1
      WRITE(numout,*) 'oasis_get_1d: Error you sould not be there...'
   END SUBROUTINE oasis_get_1d

   SUBROUTINE oasis_get_freqs(k1,k5,k2,k3,k4)
      INTEGER              , INTENT(in   ) ::  k1,k2
      INTEGER, DIMENSION(1), INTENT(  out) ::  k3
      INTEGER              , INTENT(  out) ::  k4,k5
      k3(1) = k1 ; k4 = k2 ; k5 = k2
      WRITE(numout,*) 'oasis_get_freqs: Error you sould not be there...'
   END SUBROUTINE oasis_get_freqs

   SUBROUTINE oasis_terminate(k1)
      INTEGER     , INTENT(  out) ::  k1
      k1 = -1
      WRITE(numout,*) 'oasis_terminate: Error you sould not be there...'
   END SUBROUTINE oasis_terminate

#endif

   !!=====================================================================
END MODULE cpl_oasis3
