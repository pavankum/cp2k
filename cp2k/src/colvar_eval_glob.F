!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief defines collective variables s({R}) and the derivative of this variable wrt R
!>      these can then be used in constraints, restraints and metadynamics ...
!> \par History
!>      04.2004 created
!>      01.2006 Refactored [Joost VandeVondele]
!>      11.2015 Moved force_env dependent stuff into separate module [Ole Schuett]
!> \author Alessandro Laio,Fawzi Mohamed
! *****************************************************************************
MODULE colvar_eval_glob
  USE cell_types,                      ONLY: cell_type
  USE colvar_methods,                  ONLY: check_fixed_atom_cns_colv,&
                                             colvar_eval_mol_f,&
                                             get_coordinates,&
                                             put_derivative
  USE colvar_types,                    ONLY: HBP_colvar_id,&
                                             Wc_colvar_id,&
                                             colvar_type,&
                                             u_colvar_id
  USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                             cp_to_string
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_p_type,&
                                             cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_mixed_force
  USE force_fields_util,               ONLY: get_generic_info
  USE fparser,                         ONLY: EvalErrType,&
                                             evalf,&
                                             evalfd,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE mixed_energy_types,              ONLY: mixed_force_type
  USE mixed_environment_utils,         ONLY: get_subsys_map_index
  USE particle_list_types,             ONLY: particle_list_p_type,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE string_utilities,                ONLY: compress
  USE wannier_states_types,            ONLY: wannier_centres_type
#include "./base/base_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'colvar_eval_glob'
  REAL(KIND=dp), PRIVATE, PARAMETER    :: tolerance_acos = 1.0E-5_dp

  PUBLIC :: colvar_eval_glob_f

CONTAINS

! *****************************************************************************
!> \brief evaluates the derivatives (dsdr) given and due to the given colvar
!> \param icolvar the collective variable to evaluate
!> \param force_env ...
!> \author Alessandro Laio and fawzi
!> \note
!>      The torsion that generally is defined without the continuity problem
!>      here (for free energy calculations) is defined only for (-pi,pi]
! *****************************************************************************
  SUBROUTINE colvar_eval_glob_f(icolvar,force_env)
    INTEGER                                  :: icolvar
    TYPE(force_env_type), POINTER            :: force_env

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_eval_glob_f', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: colvar_ok
    TYPE(cell_type), POINTER                 :: cell
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(qs_environment_type), POINTER       :: qs_env

    NULLIFY(subsys,cell,colvar,qs_env)
    CALL force_env_get(force_env,subsys=subsys,cell=cell,qs_env=qs_env)
    colvar_ok=ASSOCIATED(subsys%colvar_p)
    CPASSERT(colvar_ok)

    CALL cp_subsys_get(subsys,particle_set=particles)

    colvar => subsys%colvar_p(icolvar)%colvar
    ! Initialize the content of the derivative
    colvar%dsdr=0.0_dp

    SELECT CASE(colvar%type_id)
    CASE (u_colvar_id)
       CALL u_colvar(colvar,force_env=force_env)
    CASE (Wc_colvar_id)
       CALL Wc_colvar(colvar,cell,subsys=subsys, qs_env=qs_env)
    CASE (HBP_colvar_id)
       CALL HBP_colvar(colvar,cell,subsys=subsys, qs_env=qs_env)
    CASE DEFAULT
       CALL colvar_eval_mol_f(colvar, cell, particles=particles)
    END SELECT

    ! Check for fixed atom constraints
    CALL check_fixed_atom_cns_colv(subsys%gci%fixd_list,colvar)
  END SUBROUTINE colvar_eval_glob_f

! *****************************************************************************
!> \brief evaluates the force due (and on) the energy as collective variable
!> \param colvar ...
!> \param force_env ...
!> \par History Modified to allow functions of energy in a mixed_env environment
!>              Teodoro Laino [tlaino] - 02.2011
!> \author Sebastiano Caravati
! *****************************************************************************
  SUBROUTINE u_colvar(colvar,force_env)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env

    CHARACTER(len=*), PARAMETER :: routineN = 'u_colvar', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: coupling_function
    CHARACTER(LEN=default_string_length)     :: def_error, this_error
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: parameters
    INTEGER                                  :: iatom, iforce_eval, &
                                                iparticle, jparticle, natom, &
                                                natom_iforce, nforce_eval
    INTEGER, DIMENSION(:), POINTER           :: glob_natoms, map_index
    REAL(dp)                                 :: dedf, dx, err, fi(3), lerr, &
                                                potential_energy
    REAL(KIND=dp), DIMENSION(:), POINTER     :: values
    TYPE(cp_subsys_p_type), DIMENSION(:), &
      POINTER                                :: subsystems
    TYPE(cp_subsys_type), POINTER            :: subsys_main
    TYPE(mixed_force_type), DIMENSION(:), &
      POINTER                                :: global_forces
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles
    TYPE(particle_list_type), POINTER        :: particles_main
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                mapping_section, wrk_section

    IF (PRESENT(force_env)) THEN
       NULLIFY(particles_main, subsys_main)
       CALL force_env_get(force_env=force_env,subsys=subsys_main)
       CALL cp_subsys_get(subsys=subsys_main, particles=particles_main)
       natom = SIZE(particles_main%els)
       colvar%n_atom_s      = natom
       colvar%u_param%natom = natom
       CALL reallocate(colvar%i_atom,  1,natom)
       CALL reallocate(colvar%dsdr,1,3,1,natom)
       DO iatom=1,natom
          colvar%i_atom(iatom)=iatom
       ENDDO

       IF (.NOT.ASSOCIATED(colvar%u_param%mixed_energy_section)) THEN
          CALL force_env_get(force_env,potential_energy=potential_energy)
          colvar%ss            = potential_energy

          DO iatom=1,natom
             ! store derivative
             fi(:)=-particles_main%els(iatom)%f
             CALL put_derivative(colvar, iatom, fi)
          ENDDO
       ELSE
          IF(force_env%in_use/=use_mixed_force)&
             CALL cp_abort(__LOCATION__,&
                  'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
                  ' A combination of mixed force_eval energies has been requested as '//&
                  ' collective variable, but the MIXED env is not in use! Aborting.')
          CALL force_env_get(force_env, force_env_section=force_env_section)
          mapping_section => section_vals_get_subs_vals(force_env_section,"MIXED%MAPPING")
          NULLIFY(values, parameters, subsystems, particles, global_forces, map_index, glob_natoms)
          nforce_eval = SIZE(force_env%sub_force_env)
          ALLOCATE(glob_natoms(nforce_eval))
          ALLOCATE(subsystems(nforce_eval))
          ALLOCATE(particles(nforce_eval))
          ! Local Info to sync
          ALLOCATE(global_forces(nforce_eval))

          glob_natoms=0
          DO iforce_eval = 1, nforce_eval
             NULLIFY(subsystems(iforce_eval)%subsys, particles(iforce_eval)%list)
             IF (.NOT.ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
             ! Get all available subsys
             CALL force_env_get(force_env=force_env%sub_force_env(iforce_eval)%force_env,&
                                subsys=subsystems(iforce_eval)%subsys)
             ! Get available particles
             CALL cp_subsys_get(subsys=subsystems(iforce_eval)%subsys,&
                                particles=particles(iforce_eval)%list)

             ! Get Mapping index array
             natom_iforce = SIZE(particles(iforce_eval)%list%els)

             ! Only the rank 0 process collect info for each computation
             IF ( force_env%sub_force_env(iforce_eval)%force_env%para_env%mepos==&
                  force_env%sub_force_env(iforce_eval)%force_env%para_env%source) THEN
                glob_natoms(iforce_eval) = natom_iforce
             END IF
          END DO

          ! Handling Parallel execution
          CALL mp_sync(force_env%para_env%group)
          CALL mp_sum(glob_natoms, force_env%para_env%group)

          ! Transfer forces
          DO iforce_eval = 1, nforce_eval
             ALLOCATE(global_forces(iforce_eval)%forces(3,glob_natoms(iforce_eval)))
             global_forces(iforce_eval)%forces = 0.0_dp
             IF (ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) THEN
                IF ( force_env%sub_force_env(iforce_eval)%force_env%para_env%mepos==&
                     force_env%sub_force_env(iforce_eval)%force_env%para_env%source) THEN
                   ! Forces
                   DO iparticle = 1, glob_natoms(iforce_eval)
                      global_forces(iforce_eval)%forces(:,iparticle) = &
                           particles(iforce_eval)%list%els(iparticle)%f
                   END DO
                END IF
             END IF
             CALL mp_sum(global_forces(iforce_eval)%forces, force_env%para_env%group)
          END DO

          wrk_section => colvar%u_param%mixed_energy_section
          ! Support any number of force_eval sections
          CALL get_generic_info(wrk_section, "ENERGY_FUNCTION", coupling_function, parameters,&
               values, force_env%mixed_env%energies)
          CALL initf(1)
          CALL parsef(1,TRIM(coupling_function),parameters)
          ! Store the value of the COLVAR
          colvar%ss = evalf(1,values)
          CPASSERT(EvalErrType<=0)

          DO iforce_eval = 1, nforce_eval
             CALL section_vals_val_get(wrk_section,"DX",r_val=dx)
             CALL section_vals_val_get(wrk_section,"ERROR_LIMIT",r_val=lerr)
             dedf = evalfd(1,iforce_eval,values,dx,err)
             IF (ABS(err)>lerr) THEN
                WRITE(this_error,"(A,G12.6,A)")"(",err,")"
                WRITE(def_error,"(A,G12.6,A)")"(",lerr,")"
                CALL compress(this_error,.TRUE.)
                CALL compress(def_error,.TRUE.)
                CALL cp_warn(__LOCATION__,&
                     'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
                     ' Error '//TRIM(this_error)//' in computing numerical derivatives larger then'//&
                     TRIM(def_error)//' .')
             END IF
             ! General Mapping of forces...
             ! First: Get Mapping index array
             CALL get_subsys_map_index(mapping_section, glob_natoms(iforce_eval), iforce_eval,&
                                       nforce_eval, map_index)

             ! Second: store derivatives
             DO iparticle = 1, glob_natoms(iforce_eval)
                jparticle = map_index(iparticle)
                fi= - dedf * global_forces(iforce_eval)%forces(:,iparticle)
                CALL put_derivative(colvar, jparticle, fi)
             END DO
             ! Deallocate map_index array
             IF (ASSOCIATED(map_index)) THEN
                DEALLOCATE(map_index)
             END IF
          END DO
          CALL finalizef()
          DO iforce_eval = 1, nforce_eval
             DEALLOCATE(global_forces(iforce_eval)%forces)
          END DO
          DEALLOCATE(glob_natoms)
          DEALLOCATE(values)
          DEALLOCATE(parameters)
          DEALLOCATE(global_forces)
          DEALLOCATE(subsystems)
          DEALLOCATE(particles)
       END IF
    ELSE
        CPABORT("need force_env!")
    ENDIF
  END SUBROUTINE u_colvar

  
! *****************************************************************************
!> \brief evaluates the collective variable associated with a hydrogen bond
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param qs_env ...
!> \author alin m elena
! *****************************************************************************
  SUBROUTINE Wc_colvar(colvar,cell,subsys,particles,qs_env)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(qs_environment_type), POINTER, OPTIONAL       :: qs_env ! optional just because I am lazy... but I should get rid of it...

    CHARACTER(len=*), PARAMETER :: routineN = 'Wc_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: Od, H, Oa
    REAL(dp)                                 :: rOd(3), rOa(3), rH(3), &
                                                x,y,s(3),xv(3),dmin,amin
    INTEGER                                  :: idmin, iamin,i,j
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles
    TYPE(wannier_centres_type), DIMENSION(:),POINTER :: wc
    INTEGER, ALLOCATABLE                     :: wcai(:),wcdi(:) ! contains the indeces of the wannier centres closed to the donor and acceptor
    INTEGER                                  :: nwca,nwcd
    REAL(dp)                                 :: rcut

    NULLIFY(particles_i,wc)

    CPASSERT(colvar%type_id==Wc_colvar_id)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPASSERT(PRESENT(subsys))
       CALL cp_subsys_get(subsys,particles=particles_i)
       my_particles => particles_i%els
    END IF
    CALL get_qs_env(qs_env,WannierCentres=wc)
    rcut =colvar%Wc%rcut ! distances are in bohr as far as I remember
    Od=colvar%Wc%ids(1)
    H=colvar%Wc%ids(2)
    Oa=colvar%Wc%ids(3)
    CALL get_coordinates(colvar, Od, rOd, my_particles)
    CALL get_coordinates(colvar, H, rH, my_particles)
    CALL get_coordinates(colvar, Oa, rOa, my_particles)
    ALLOCATE(wcai(SIZE(wc(1)%WannierHamDiag)))
    ALLOCATE(wcdi(SIZE(wc(1)%WannierHamDiag)))
    nwca=0
    nwcd=0
    DO j=1,SIZE(wc(1)%WannierHamDiag)
      x=distance(rOd-wc(1)%centres(:,j))
      y=distance(rOa-wc(1)%centres(:,j))
      IF (x<rcut) THEN
        nwcd=nwcd+1
        wcdi(nwcd)=j
        CYCLE
      ENDIF
      IF (y<rcut) THEN
        nwca=nwca+1
        wcai(nwca)=j
      ENDIF
    ENDDO

    dmin=distance(rH-wc(1)%centres(:,wcdi(1)))
    amin=distance(rH-wc(1)%centres(:,wcai(1)))
    idmin=wcdi(1)
    iamin=wcai(1)
    !dmin constains the smallest numer, amin the next smallest
    DO i=2,nwcd
      x=distance(rH-wc(1)%centres(:,wcdi(i)))
      IF (x<dmin) THEN
        dmin=x
        idmin=wcdi(i)
      ENDIF
    ENDDO
    DO i=2,nwca
      x=distance(rH-wc(1)%centres(:,wcai(i)))
      IF (x<amin) THEN
        amin=x
        iamin=wcai(i)
      ENDIF
    ENDDO
!     zero=0.0_dp
!     CALL put_derivative(colvar, 1, zero)
!     CALL put_derivative(colvar, 2,zero)
!     CALL put_derivative(colvar, 3, zero)

!     write(*,'(2(i0,1x),4(f16.8,1x))')idmin,iamin,wc(1)%WannierHamDiag(idmin),wc(1)%WannierHamDiag(iamin),dmin,amin
    colvar%ss = wc(1)%WannierHamDiag(idmin)-wc(1)%WannierHamDiag(iamin)
    DEALLOCATE(wcai)
    DEALLOCATE(wcdi)

    CONTAINS
! *****************************************************************************
!> \brief ...
!> \param rij ...
!> \retval distance ...
! *****************************************************************************
    REAL(dp) FUNCTION distance(rij)
    REAL(dp), INTENT(in) :: rij(3)

      s=MATMUL(cell%h_inv,rij)
      s=s-NINT(s)
      xv=MATMUL(cell%hmat,s)
      distance=SQRT(DOT_PRODUCT(xv,xv))
    END FUNCTION distance

  END SUBROUTINE Wc_colvar

! *****************************************************************************
!> \brief evaluates the collective variable associated with a hydrogen bond wire
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param qs_env ...
!> \author alin m elena
! *****************************************************************************
  SUBROUTINE HBP_colvar(colvar,cell,subsys,particles,qs_env)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(qs_environment_type), POINTER, OPTIONAL       :: qs_env ! optional just because I am lazy... but I should get rid of it...

    CHARACTER(len=*), PARAMETER :: routineN = 'HBP_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: Od, H, Oa
    REAL(dp)                                 :: rOd(3), rOa(3), rH(3), &
                                                x,y,s(3),xv(3),dmin,amin
    INTEGER                                  :: idmin, iamin,i,j,il,output_unit
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles
    TYPE(wannier_centres_type),&
                        DIMENSION(:),POINTER :: wc
    INTEGER, ALLOCATABLE                     :: wcai(:),wcdi(:) ! contains the indeces of the wannier centres closed to the donor and acceptor
    INTEGER                                  :: nwca,nwcd
    REAL(dp)                                 :: rcut



    NULLIFY(particles_i,wc)
    output_unit= cp_logger_get_default_io_unit()

    CPASSERT(colvar%type_id==HBP_colvar_id)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPASSERT(PRESENT(subsys))
       CALL cp_subsys_get(subsys,particles=particles_i)
       my_particles => particles_i%els
    END IF
    CALL get_qs_env(qs_env,WannierCentres=wc)
    rcut = colvar%HBP%rcut ! distances are in bohr as far as I remember
    ALLOCATE(wcai(SIZE(wc(1)%WannierHamDiag)))
    ALLOCATE(wcdi(SIZE(wc(1)%WannierHamDiag)))
    colvar%ss=0.0_dp
    DO il=1,colvar%HBP%nPoints
      Od=colvar%HBP%ids(il,1)
      H=colvar%HBP%ids(il,2)
      Oa=colvar%HBP%ids(il,3)
      CALL get_coordinates(colvar, Od, rOd, my_particles)
      CALL get_coordinates(colvar, H, rH, my_particles)
      CALL get_coordinates(colvar, Oa, rOa, my_particles)
      nwca=0
      nwcd=0
      DO j=1,SIZE(wc(1)%WannierHamDiag)
        x=distance(rOd-wc(1)%centres(:,j))
        y=distance(rOa-wc(1)%centres(:,j))
        IF (x<rcut) THEN
          nwcd=nwcd+1
          wcdi(nwcd)=j
          CYCLE
        ENDIF
        IF (y<rcut) THEN
          nwca=nwca+1
          wcai(nwca)=j
        ENDIF
      ENDDO

      dmin=distance(rH-wc(1)%centres(:,wcdi(1)))
      amin=distance(rH-wc(1)%centres(:,wcai(1)))
      idmin=wcdi(1)
      iamin=wcai(1)
      !dmin constains the smallest numer, amin the next smallest
      DO i=2,nwcd
        x=distance(rH-wc(1)%centres(:,wcdi(i)))
        IF (x<dmin) THEN
          dmin=x
          idmin=wcdi(i)
        ENDIF
      ENDDO
      DO i=2,nwca
        x=distance(rH-wc(1)%centres(:,wcai(i)))
        IF (x<amin) THEN
          amin=x
          iamin=wcai(i)
        ENDIF
      ENDDO
!       if (output_unit>0) then
!          write(output_unit,'(a,2(i0,1x),5(f16.8,1x))')"HBP|",idmin,iamin,wc(1)%WannierHamDiag(idmin),wc(1)%WannierHamDiag(iamin),dmin,amin, wc(1)%WannierHamDiag(idmin)-wc(1)%WannierHamDiag(iamin)
!       endif
      colvar%HBP%ewc(il) = colvar%HBP%shift + wc(1)%WannierHamDiag(idmin)-wc(1)%WannierHamDiag(iamin)
      colvar%ss = colvar%ss + colvar%HBP%shift + wc(1)%WannierHamDiag(idmin)-wc(1)%WannierHamDiag(iamin)
    ENDDO
     IF (output_unit>0) THEN
       DO il=1,colvar%HBP%nPoints
         WRITE(output_unit,'(a,1(f16.8,1x))')"HBP| = ",colvar%HBP%ewc(il)
       ENDDO
       WRITE(output_unit,'(a,1(f16.8,1x))')"HBP|\theta(x) = ",colvar%ss
    ENDIF
    DEALLOCATE(wcai)
    DEALLOCATE(wcdi)

  CONTAINS
! *****************************************************************************
!> \brief ...
!> \param rij ...
!> \retval distance ...
! *****************************************************************************
    REAL(dp) FUNCTION distance(rij)
    REAL(dp), INTENT(in) :: rij(3)

      s=MATMUL(cell%h_inv,rij)
      s=s-NINT(s)
      xv=MATMUL(cell%hmat,s)
      distance=SQRT(DOT_PRODUCT(xv,xv))
    END FUNCTION distance

  END SUBROUTINE HBP_colvar

END MODULE colvar_eval_glob
