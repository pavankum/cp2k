!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Reading of input parameters for the pw_poisson-modules.
!> \par History
!>      01.2014 Code moved into separate module to make pw_poisson-modules
!>              independet from input_section_types and input_constants.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pw_poisson_read_input
   USE cell_types,                      ONLY: use_perd_none,&
                                              use_perd_x,&
                                              use_perd_xy,&
                                              use_perd_xyz,&
                                              use_perd_xz,&
                                              use_perd_y,&
                                              use_perd_yz,&
                                              use_perd_z
   USE dirichlet_bc_types,              ONLY: xy_plane,&
                                              xz_plane,&
                                              yz_plane
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE ps_wavelet_types,                ONLY: WAVELET0D,&
                                              WAVELET2D,&
                                              WAVELET3D
   USE pw_poisson_types,                ONLY: &
        do_ewald_none, pw_poisson_analytic, pw_poisson_implicit, pw_poisson_mt, &
        pw_poisson_multipole, pw_poisson_none, pw_poisson_parameter_type, pw_poisson_periodic, &
        pw_poisson_wavelet
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_poisson_read_input'

   PUBLIC :: pw_poisson_read_parameters

CONTAINS

! **************************************************************************************************
!> \brief Reads the POISSON input-section and into pw_poisson_parameter_type.
!> \param poisson_section ...
!> \param params ...
!> \par History
!>      01.2014 Code moved into separate module from pw_poisson_types,
!>              pw_poisson_methods and ps_wavelet_types.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE pw_poisson_read_parameters(poisson_section, params)
      TYPE(section_vals_type), POINTER                   :: poisson_section
      TYPE(pw_poisson_parameter_type), INTENT(INOUT)     :: params

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_poisson_read_parameters', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: periodic
      TYPE(section_vals_type), POINTER                   :: ewald_section, mt_section, &
                                                            wavelet_section

      NULLIFY (ewald_section, mt_section, wavelet_section)

      CALL section_vals_val_get(poisson_section, "POISSON_SOLVER", i_val=params%solver)

      ! Decoding PERIODIC depending on chosen solver,
      ! because not all solvers support every possible periodicity
      CALL section_vals_val_get(poisson_section, "PERIODIC", i_val=periodic)
      SELECT CASE (params%solver)
      CASE (pw_poisson_periodic, pw_poisson_analytic, pw_poisson_mt, pw_poisson_multipole, &
            pw_poisson_implicit)
         CALL decode_periodic_green(periodic, params)
      CASE (pw_poisson_wavelet)
         CALL decode_periodic_wavelet(periodic, params)
      CASE (pw_poisson_none)
      CASE default
         CPABORT("")
      END SELECT

      ! parsing EWALD subsection
      params%ewald_type = do_ewald_none
      ewald_section => section_vals_get_subs_vals(poisson_section, "EWALD", can_return_null=.TRUE.)
      IF (ASSOCIATED(ewald_section)) THEN
         CALL section_vals_val_get(ewald_section, "EWALD_TYPE", i_val=params%ewald_type)
         CALL section_vals_val_get(ewald_section, "o_spline", i_val=params%ewald_o_spline)
         CALL section_vals_val_get(ewald_section, "alpha", r_val=params%ewald_alpha)
      ENDIF

      ! parsing MT subsection
      mt_section => section_vals_get_subs_vals(poisson_section, "MT")
      CALL section_vals_val_get(mt_section, "REL_CUTOFF", r_val=params%mt_rel_cutoff)
      CALL section_vals_val_get(mt_section, "ALPHA", r_val=params%mt_alpha)

      ! parsing WAVELET subsection
      wavelet_section => section_vals_get_subs_vals(poisson_section, "WAVELET")
      CALL section_vals_val_get(wavelet_section, "SCF_TYPE", i_val=params%wavelet_scf_type)

      ! parsing IMPLICIT subsection
      CALL ps_implicit_read_parameters(poisson_section, params)

   END SUBROUTINE pw_poisson_read_parameters

! **************************************************************************************************
!> \brief Helper routien for pw_poisson_read_parameters
!> \param periodic ...
!> \param params ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE decode_periodic_green(periodic, params)
      INTEGER, INTENT(IN)                                :: periodic
      TYPE(pw_poisson_parameter_type), INTENT(INOUT)     :: params

      CHARACTER(len=*), PARAMETER :: routineN = 'decode_periodic_green', &
         routineP = moduleN//':'//routineN

      SELECT CASE (periodic)
      CASE (use_perd_x)
         params%periodic = (/1, 0, 0/)
      CASE (use_perd_y)
         params%periodic = (/0, 1, 0/)
      CASE (use_perd_z)
         params%periodic = (/0, 0, 1/)
      CASE (use_perd_xy)
         params%periodic = (/1, 1, 0/)
      CASE (use_perd_xz)
         params%periodic = (/1, 0, 1/)
      CASE (use_perd_yz)
         params%periodic = (/0, 1, 1/)
      CASE (use_perd_xyz)
         params%periodic = (/1, 1, 1/)
      CASE (use_perd_none)
         params%periodic = (/0, 0, 0/)
      CASE DEFAULT
         CPABORT("")
      END SELECT
      ! check for consistent use of periodicity (cell <-> Poisson solver)
      !CPPostcondition(ALL(perd == cell%perd),cp_fatal_level,routineP,failure)

   END SUBROUTINE decode_periodic_green

! **************************************************************************************************
!> \brief Helper routien for pw_poisson_read_parameters
!> \param periodic ...
!> \param params ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE decode_periodic_wavelet(periodic, params)
      INTEGER, INTENT(IN)                                :: periodic
      TYPE(pw_poisson_parameter_type), INTENT(INOUT)     :: params

      CHARACTER(len=*), PARAMETER :: routineN = 'decode_periodic_wavelet', &
         routineP = moduleN//':'//routineN

      params%wavelet_special_dimension = 0

      SELECT CASE (periodic)
      CASE (use_perd_none)
         params%wavelet_geocode = "F"
         params%wavelet_method = WAVELET0D
      CASE (use_perd_xz)
         params%wavelet_geocode = "S"
         params%wavelet_method = WAVELET2D
         params%wavelet_special_dimension = 2
      CASE (use_perd_xyz)
         params%wavelet_geocode = "P"
         params%wavelet_method = WAVELET3D
      CASE (use_perd_x, use_perd_y, use_perd_z, use_perd_xy, use_perd_yz)
         CPABORT("Poisson solver for this periodicity not yet implemented")
      CASE DEFAULT
         CPABORT("")
      END SELECT

   END SUBROUTINE decode_periodic_wavelet

! **************************************************************************************************
!> \brief Reads the subsection IMPLICIT and initializes corresponding parameters in
!>        pw_poisson_parameter_type
!> \param poisson_section poisson section to be read from input
!> \param params poisson_env parameters
!> \par History
!>      08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE ps_implicit_read_parameters(poisson_section, params)
      TYPE(section_vals_type), POINTER                   :: poisson_section
      TYPE(pw_poisson_parameter_type), INTENT(INOUT)     :: params

      CHARACTER(len=*), PARAMETER :: routineN = 'ps_implicit_read_parameters', &
         routineP = moduleN//':'//routineN

      LOGICAL                                            :: has_dielectric
      TYPE(section_vals_type), POINTER                   :: dbc_section, dielectric_section, &
                                                            ps_implicit_section

      NULLIFY (ps_implicit_section, dielectric_section, dbc_section)

      ! parsing IMPLICIT subsection
      ps_implicit_section => section_vals_get_subs_vals(poisson_section, "IMPLICIT")
      CALL section_vals_val_get(ps_implicit_section, "BOUNDARY_CONDITIONS", &
                                i_val=params%ps_implicit_params%boundary_condition)
      CALL section_vals_val_get(ps_implicit_section, "ZERO_INITIAL_GUESS", &
                                l_val=params%ps_implicit_params%zero_initial_guess)
      CALL section_vals_val_get(ps_implicit_section, "max_iter", i_val=params%ps_implicit_params%max_iter)
      CALL section_vals_val_get(ps_implicit_section, "tol", r_val=params%ps_implicit_params%tol)
      CALL section_vals_val_get(ps_implicit_section, "omega", r_val=params%ps_implicit_params%omega)
      CALL section_vals_val_get(ps_implicit_section, "neumann_directions", &
                                i_val=params%ps_implicit_params%neumann_directions)

      ! parsing DIELECTRIC subsection
      dielectric_section => section_vals_get_subs_vals(ps_implicit_section, "DIELECTRIC")
      CALL section_vals_get(dielectric_section, explicit=has_dielectric)
      params%has_dielectric = has_dielectric
      CALL dielectric_read_parameters(dielectric_section, params)

      ! parsing DIRICHLET_BC subsection
      dbc_section => section_vals_get_subs_vals(ps_implicit_section, "DIRICHLET_BC")
      CALL dirichlet_bc_read_parameters(dbc_section, params)

   END SUBROUTINE ps_implicit_read_parameters

! **************************************************************************************************
!> \brief Reads the subsection DIELECTRIC and initializes corresponding parameters in
!>        pw_poisson_parameter_type
!> \param dielectric_section dielectric section to be read from input
!> \param params poisson_env parameters
!> \par History
!>      07.2015 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE dielectric_read_parameters(dielectric_section, params)
      TYPE(section_vals_type), POINTER                   :: dielectric_section
      TYPE(pw_poisson_parameter_type), INTENT(INOUT)     :: params

      CHARACTER(len=*), PARAMETER :: routineN = 'dielectric_read_parameters', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, n_aac_rep, n_xaaa_rep
      LOGICAL                                            :: aa_cuboidal_explicit, &
                                                            xaa_annular_explicit
      REAL(dp)                                           :: eps, zeta
      REAL(dp), DIMENSION(:), POINTER :: aa_cuboidal_xxtnt, aa_cuboidal_yxtnt, aa_cuboidal_zxtnt, &
         xaa_annular_bctr, xaa_annular_brad, xaa_annular_xxtnt
      TYPE(section_vals_type), POINTER                   :: aa_cuboidal_section, xaa_annular_section

      CALL section_vals_val_get(dielectric_section, "DIELECTRIC_CORE_CORRECTION", &
                                l_val=params%dielectric_params%dielec_core_correction)
      CALL section_vals_val_get(dielectric_section, "DIELECTRIC_FUNCTION_TYPE", &
                                i_val=params%dielectric_params%dielec_functiontype)
      CALL section_vals_val_get(dielectric_section, "epsilon", r_val=params%dielectric_params%eps0)
      CALL section_vals_val_get(dielectric_section, "rho_min", r_val=params%dielectric_params%rho_min)
      CALL section_vals_val_get(dielectric_section, "rho_max", r_val=params%dielectric_params%rho_max)
      CALL section_vals_val_get(dielectric_section, "DERIVATIVE_METHOD", &
                                i_val=params%dielectric_params%derivative_method)

      aa_cuboidal_section => section_vals_get_subs_vals(dielectric_section, "DIELEC_AA_CUBOIDAL")
      xaa_annular_section => section_vals_get_subs_vals(dielectric_section, "DIELEC_XAA_ANNULAR")
      CALL section_vals_get(aa_cuboidal_section, explicit=aa_cuboidal_explicit, n_repetition=n_aac_rep)
      CALL section_vals_get(xaa_annular_section, explicit=xaa_annular_explicit, n_repetition=n_xaaa_rep)

      IF (params%solver .EQ. pw_poisson_implicit) THEN

         IF (aa_cuboidal_explicit) THEN
            params%dielectric_params%n_aa_cuboidal = n_aac_rep
            ALLOCATE (params%dielectric_params%aa_cuboidal_xxtnt(2, n_aac_rep), &
                      params%dielectric_params%aa_cuboidal_yxtnt(2, n_aac_rep), &
                      params%dielectric_params%aa_cuboidal_zxtnt(2, n_aac_rep), &
                      params%dielectric_params%aa_cuboidal_eps(n_aac_rep), &
                      params%dielectric_params%aa_cuboidal_zeta(n_aac_rep))
            NULLIFY (aa_cuboidal_xxtnt, aa_cuboidal_yxtnt, aa_cuboidal_zxtnt)
            DO i = 1, n_aac_rep
               CALL section_vals_val_get(aa_cuboidal_section, "epsilon", i_rep_section=i, r_val=eps)
               CALL section_vals_val_get(aa_cuboidal_section, "zeta", i_rep_section=i, r_val=zeta)
               CALL section_vals_val_get(aa_cuboidal_section, "X_xtnt", i_rep_section=i, r_vals=aa_cuboidal_xxtnt)
               CALL section_vals_val_get(aa_cuboidal_section, "Y_xtnt", i_rep_section=i, r_vals=aa_cuboidal_yxtnt)
               CALL section_vals_val_get(aa_cuboidal_section, "Z_xtnt", i_rep_section=i, r_vals=aa_cuboidal_zxtnt)
               params%dielectric_params%aa_cuboidal_eps(i) = eps
               params%dielectric_params%aa_cuboidal_zeta(i) = zeta
               params%dielectric_params%aa_cuboidal_xxtnt(:, i) = aa_cuboidal_xxtnt
               params%dielectric_params%aa_cuboidal_yxtnt(:, i) = aa_cuboidal_yxtnt
               params%dielectric_params%aa_cuboidal_zxtnt(:, i) = aa_cuboidal_zxtnt
            END DO
         ELSE
            params%dielectric_params%n_aa_cuboidal = 0
         END IF

         IF (xaa_annular_explicit) THEN
            params%dielectric_params%n_xaa_annular = n_xaaa_rep
            ALLOCATE (params%dielectric_params%xaa_annular_xxtnt(2, n_xaaa_rep), &
                      params%dielectric_params%xaa_annular_bctr(2, n_xaaa_rep), &
                      params%dielectric_params%xaa_annular_brad(2, n_xaaa_rep), &
                      params%dielectric_params%xaa_annular_eps(n_xaaa_rep), &
                      params%dielectric_params%xaa_annular_zeta(n_xaaa_rep))
            NULLIFY (xaa_annular_xxtnt, xaa_annular_bctr, xaa_annular_brad)
            DO i = 1, n_xaaa_rep
               CALL section_vals_val_get(xaa_annular_section, "epsilon", i_rep_section=i, r_val=eps)
               CALL section_vals_val_get(xaa_annular_section, "zeta", i_rep_section=i, r_val=zeta)
               CALL section_vals_val_get(xaa_annular_section, "X_xtnt", i_rep_section=i, r_vals=xaa_annular_xxtnt)
               CALL section_vals_val_get(xaa_annular_section, "BASE_CENTER", i_rep_section=i, r_vals=xaa_annular_bctr)
               CALL section_vals_val_get(xaa_annular_section, "BASE_RADII", i_rep_section=i, r_vals=xaa_annular_brad)
               params%dielectric_params%xaa_annular_eps(i) = eps
               params%dielectric_params%xaa_annular_zeta(i) = zeta
               params%dielectric_params%xaa_annular_xxtnt(:, i) = xaa_annular_xxtnt
               params%dielectric_params%xaa_annular_bctr(:, i) = xaa_annular_bctr
               params%dielectric_params%xaa_annular_brad(:, i) = xaa_annular_brad
            END DO
         ELSE
            params%dielectric_params%n_xaa_annular = 0
         END IF

      END IF

   END SUBROUTINE dielectric_read_parameters

! **************************************************************************************************
!> \brief Reads the subsection DIRICHLET_BC and initializes corresponding parameters in
!>        pw_poisson_parameter_type
!> \param dbc_section dirichlet_bc section to be read from input
!> \param params poisson_env parameters
!> \par History
!>      08.2014 created [Hossein Bani-Hashemian]
!>      07.2015 refactored [Hossein Bani-Hashemian]
!>      10.2015 revised [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE dirichlet_bc_read_parameters(dbc_section, params)
      TYPE(section_vals_type), POINTER                   :: dbc_section
      TYPE(pw_poisson_parameter_type), INTENT(INOUT)     :: params

      CHARACTER(len=*), PARAMETER :: routineN = 'dirichlet_bc_read_parameters', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, n_aac_rep, n_aap_rep, n_p_rep, &
                                                            n_xaac_rep, parallel_plane, &
                                                            xaa_cylindrical_apxtyp, &
                                                            xaa_cylindrical_nsides
      INTEGER, DIMENSION(:), POINTER                     :: aa_cuboidal_nprtn, aa_planar_nprtn, &
                                                            planar_nprtn, xaa_cylindrical_nprtn
      LOGICAL                                            :: aa_cuboidal_explicit, &
                                                            aa_planar_explicit, planar_explicit, &
                                                            xaa_cylindrical_explicit
      REAL(dp)                                           :: delta_alpha, intercept, sigma, &
                                                            thickness, v_D, xaa_cylindrical_brad
      REAL(dp), DIMENSION(:), POINTER :: aa_cuboidal_xxtnt, aa_cuboidal_yxtnt, aa_cuboidal_zxtnt, &
         aa_planar_xxtnt, aa_planar_yxtnt, aa_planar_zxtnt, planar_Avtx, planar_Bvtx, planar_Cvtx, &
         xaa_cylindrical_bctr, xaa_cylindrical_xxtnt
      TYPE(section_vals_type), POINTER                   :: aa_cuboidal_section, aa_planar_section, &
                                                            planar_section, xaa_cylindrical_section

      CALL section_vals_val_get(dbc_section, "VERBOSE_OUTPUT", l_val=params%dbc_params%verbose_output)
      aa_planar_section => section_vals_get_subs_vals(dbc_section, "AA_PLANAR")
      planar_section => section_vals_get_subs_vals(dbc_section, "PLANAR")
      xaa_cylindrical_section => section_vals_get_subs_vals(dbc_section, "XAA_CYLINDRICAL")
      aa_cuboidal_section => section_vals_get_subs_vals(dbc_section, "AA_CUBOIDAL")
      CALL section_vals_get(aa_planar_section, explicit=aa_planar_explicit, n_repetition=n_aap_rep)
      CALL section_vals_get(planar_section, explicit=planar_explicit, n_repetition=n_p_rep)
      CALL section_vals_get(xaa_cylindrical_section, explicit=xaa_cylindrical_explicit, n_repetition=n_xaac_rep)
      CALL section_vals_get(aa_cuboidal_section, explicit=aa_cuboidal_explicit, n_repetition=n_aac_rep)

      IF (params%solver .EQ. pw_poisson_implicit) THEN

         IF (aa_planar_explicit) THEN
            params%dbc_params%n_aa_planar = n_aap_rep
            ALLOCATE (params%dbc_params%aa_planar_nprtn(3, n_aap_rep), &
                      params%dbc_params%aa_planar_pplane(n_aap_rep), &
                      params%dbc_params%aa_planar_xxtnt(2, n_aap_rep), &
                      params%dbc_params%aa_planar_yxtnt(2, n_aap_rep), &
                      params%dbc_params%aa_planar_zxtnt(2, n_aap_rep), &
                      params%dbc_params%aa_planar_vD(n_aap_rep), &
                      params%dbc_params%aa_planar_sigma(n_aap_rep), &
                      params%dbc_params%aa_planar_thickness(n_aap_rep))
            NULLIFY (aa_planar_nprtn, aa_planar_xxtnt, aa_planar_yxtnt, aa_planar_zxtnt)
            DO i = 1, n_aap_rep
               CALL section_vals_val_get(aa_planar_section, "v_D", i_rep_section=i, r_val=v_D)
               CALL section_vals_val_get(aa_planar_section, "SIGMA", i_rep_section=i, r_val=sigma)
               CALL section_vals_val_get(aa_planar_section, "THICKNESS", i_rep_section=i, r_val=thickness)

               params%dbc_params%aa_planar_vD(i) = v_D
               params%dbc_params%aa_planar_sigma(i) = sigma
               params%dbc_params%aa_planar_thickness(i) = thickness

               CALL section_vals_val_get(aa_planar_section, "PARALLEL_PLANE", i_rep_section=i, i_val=parallel_plane)
               CALL section_vals_val_get(aa_planar_section, "INTERCEPT", i_rep_section=i, r_val=intercept)
               SELECT CASE (parallel_plane)
               CASE (xy_plane)
                  params%dbc_params%aa_planar_pplane(i) = xy_plane
                  CALL section_vals_val_get(aa_planar_section, "X_xtnt", i_rep_section=i, r_vals=aa_planar_xxtnt)
                  CALL section_vals_val_get(aa_planar_section, "Y_xtnt", i_rep_section=i, r_vals=aa_planar_yxtnt)
                  params%dbc_params%aa_planar_xxtnt(:, i) = aa_planar_xxtnt
                  params%dbc_params%aa_planar_yxtnt(:, i) = aa_planar_yxtnt
                  params%dbc_params%aa_planar_zxtnt(:, i) = intercept

                  CALL section_vals_val_get(aa_planar_section, "N_PRTN", i_rep_section=i, i_vals=aa_planar_nprtn)
                  params%dbc_params%aa_planar_nprtn(1, i) = aa_planar_nprtn(1)
                  params%dbc_params%aa_planar_nprtn(2, i) = aa_planar_nprtn(2)
                  params%dbc_params%aa_planar_nprtn(3, i) = 1
               CASE (yz_plane)
                  params%dbc_params%aa_planar_pplane(i) = yz_plane
                  CALL section_vals_val_get(aa_planar_section, "Y_xtnt", i_rep_section=i, r_vals=aa_planar_yxtnt)
                  CALL section_vals_val_get(aa_planar_section, "Z_xtnt", i_rep_section=i, r_vals=aa_planar_zxtnt)
                  params%dbc_params%aa_planar_xxtnt(:, i) = intercept
                  params%dbc_params%aa_planar_yxtnt(:, i) = aa_planar_yxtnt
                  params%dbc_params%aa_planar_zxtnt(:, i) = aa_planar_zxtnt

                  CALL section_vals_val_get(aa_planar_section, "N_PRTN", i_rep_section=i, i_vals=aa_planar_nprtn)
                  params%dbc_params%aa_planar_nprtn(1, i) = 1
                  params%dbc_params%aa_planar_nprtn(2, i) = aa_planar_nprtn(1)
                  params%dbc_params%aa_planar_nprtn(3, i) = aa_planar_nprtn(2)
               CASE (xz_plane)
                  params%dbc_params%aa_planar_pplane(i) = xz_plane
                  CALL section_vals_val_get(aa_planar_section, "X_xtnt", i_rep_section=i, r_vals=aa_planar_xxtnt)
                  CALL section_vals_val_get(aa_planar_section, "Z_xtnt", i_rep_section=i, r_vals=aa_planar_zxtnt)
                  params%dbc_params%aa_planar_xxtnt(:, i) = aa_planar_xxtnt
                  params%dbc_params%aa_planar_yxtnt(:, i) = intercept
                  params%dbc_params%aa_planar_zxtnt(:, i) = aa_planar_zxtnt

                  CALL section_vals_val_get(aa_planar_section, "N_PRTN", i_rep_section=i, i_vals=aa_planar_nprtn)
                  params%dbc_params%aa_planar_nprtn(1, i) = aa_planar_nprtn(1)
                  params%dbc_params%aa_planar_nprtn(2, i) = 1
                  params%dbc_params%aa_planar_nprtn(3, i) = aa_planar_nprtn(2)
               END SELECT

            END DO
         ELSE
            params%dbc_params%n_aa_planar = 0
         END IF

         IF (planar_explicit) THEN
            params%dbc_params%n_planar = n_p_rep
            ALLOCATE (params%dbc_params%planar_nprtn(2, n_p_rep), &
                      params%dbc_params%planar_Avtx(3, n_p_rep), &
                      params%dbc_params%planar_Bvtx(3, n_p_rep), &
                      params%dbc_params%planar_Cvtx(3, n_p_rep), &
                      params%dbc_params%planar_vD(n_p_rep), &
                      params%dbc_params%planar_sigma(n_p_rep), &
                      params%dbc_params%planar_thickness(n_p_rep))
            NULLIFY (planar_nprtn, planar_Avtx, planar_Bvtx, planar_Cvtx)
            DO i = 1, n_p_rep
               CALL section_vals_val_get(planar_section, "N_PRTN", i_rep_section=i, i_vals=planar_nprtn)
               CALL section_vals_val_get(planar_section, "A", i_rep_section=i, r_vals=planar_Avtx)
               CALL section_vals_val_get(planar_section, "B", i_rep_section=i, r_vals=planar_Bvtx)
               CALL section_vals_val_get(planar_section, "C", i_rep_section=i, r_vals=planar_Cvtx)
               CALL section_vals_val_get(planar_section, "v_D", i_rep_section=i, r_val=v_D)
               CALL section_vals_val_get(planar_section, "SIGMA", i_rep_section=i, r_val=sigma)
               CALL section_vals_val_get(planar_section, "THICKNESS", i_rep_section=i, r_val=thickness)
               params%dbc_params%planar_nprtn(:, i) = planar_nprtn
               params%dbc_params%planar_Avtx(:, i) = planar_Avtx
               params%dbc_params%planar_Bvtx(:, i) = planar_Bvtx
               params%dbc_params%planar_Cvtx(:, i) = planar_Cvtx
               params%dbc_params%planar_vD(i) = v_D
               params%dbc_params%planar_sigma(i) = sigma
               params%dbc_params%planar_thickness(i) = thickness
            END DO
         ELSE
            params%dbc_params%n_planar = 0
         END IF

         IF (xaa_cylindrical_explicit) THEN
            params%dbc_params%n_xaa_cylindrical = n_xaac_rep
            ALLOCATE (params%dbc_params%xaa_cylindrical_nprtn(2, n_xaac_rep), &
                      params%dbc_params%xaa_cylindrical_nsides(n_xaac_rep), &
                      params%dbc_params%xaa_cylindrical_apxtyp(n_xaac_rep), &
                      params%dbc_params%xaa_cylindrical_xxtnt(2, n_xaac_rep), &
                      params%dbc_params%xaa_cylindrical_bctr(2, n_xaac_rep), &
                      params%dbc_params%xaa_cylindrical_brad(n_xaac_rep), &
                      params%dbc_params%xaa_cylindrical_vD(n_xaac_rep), &
                      params%dbc_params%xaa_cylindrical_sigma(n_xaac_rep), &
                      params%dbc_params%xaa_cylindrical_thickness(n_xaac_rep), &
                      params%dbc_params%xaa_cylindrical_sgap(n_xaac_rep))
            NULLIFY (xaa_cylindrical_nprtn, xaa_cylindrical_xxtnt, xaa_cylindrical_bctr)
            DO i = 1, n_xaac_rep
               CALL section_vals_val_get(xaa_cylindrical_section, "N_PRTN", i_rep_section=i, i_vals=xaa_cylindrical_nprtn)
               CALL section_vals_val_get(xaa_cylindrical_section, "N_SIDES", i_rep_section=i, i_val=xaa_cylindrical_nsides)
               CALL section_vals_val_get(xaa_cylindrical_section, "APX_TYPE", i_rep_section=i, i_val=xaa_cylindrical_apxtyp)
               CALL section_vals_val_get(xaa_cylindrical_section, "X_xtnt", i_rep_section=i, r_vals=xaa_cylindrical_xxtnt)
               CALL section_vals_val_get(xaa_cylindrical_section, "BASE_CENTER", i_rep_section=i, r_vals=xaa_cylindrical_bctr)
               CALL section_vals_val_get(xaa_cylindrical_section, "BASE_RADIUS", i_rep_section=i, r_val=xaa_cylindrical_brad)
               CALL section_vals_val_get(xaa_cylindrical_section, "v_D", i_rep_section=i, r_val=v_D)
               CALL section_vals_val_get(xaa_cylindrical_section, "SIGMA", i_rep_section=i, r_val=sigma)
               CALL section_vals_val_get(xaa_cylindrical_section, "THICKNESS", i_rep_section=i, r_val=thickness)
               CALL section_vals_val_get(xaa_cylindrical_section, "delta_alpha", i_rep_section=i, r_val=delta_alpha)
               params%dbc_params%xaa_cylindrical_nprtn(:, i) = xaa_cylindrical_nprtn
               params%dbc_params%xaa_cylindrical_nsides(i) = xaa_cylindrical_nsides
               params%dbc_params%xaa_cylindrical_apxtyp(i) = xaa_cylindrical_apxtyp
               params%dbc_params%xaa_cylindrical_xxtnt(:, i) = xaa_cylindrical_xxtnt
               params%dbc_params%xaa_cylindrical_bctr(:, i) = xaa_cylindrical_bctr
               params%dbc_params%xaa_cylindrical_brad(i) = xaa_cylindrical_brad
               params%dbc_params%xaa_cylindrical_vD(i) = v_D
               params%dbc_params%xaa_cylindrical_sigma(i) = sigma
               params%dbc_params%xaa_cylindrical_thickness(i) = thickness
               params%dbc_params%xaa_cylindrical_sgap(i) = delta_alpha
            END DO
         ELSE
            params%dbc_params%n_xaa_cylindrical = 0
            ALLOCATE (params%dbc_params%xaa_cylindrical_nsides(n_xaac_rep))
         END IF

         IF (aa_cuboidal_explicit) THEN
            params%dbc_params%n_aa_cuboidal = n_aac_rep
            ALLOCATE (params%dbc_params%aa_cuboidal_nprtn(3, n_aac_rep), &
                      params%dbc_params%aa_cuboidal_xxtnt(2, n_aac_rep), &
                      params%dbc_params%aa_cuboidal_yxtnt(2, n_aac_rep), &
                      params%dbc_params%aa_cuboidal_zxtnt(2, n_aac_rep), &
                      params%dbc_params%aa_cuboidal_vD(n_aac_rep), &
                      params%dbc_params%aa_cuboidal_sigma(n_aac_rep))
            NULLIFY (aa_cuboidal_nprtn, aa_cuboidal_xxtnt, aa_cuboidal_yxtnt, aa_cuboidal_zxtnt)
            DO i = 1, n_aac_rep
               CALL section_vals_val_get(aa_cuboidal_section, "N_PRTN", i_rep_section=i, i_vals=aa_cuboidal_nprtn)
               CALL section_vals_val_get(aa_cuboidal_section, "X_xtnt", i_rep_section=i, r_vals=aa_cuboidal_xxtnt)
               CALL section_vals_val_get(aa_cuboidal_section, "Y_xtnt", i_rep_section=i, r_vals=aa_cuboidal_yxtnt)
               CALL section_vals_val_get(aa_cuboidal_section, "Z_xtnt", i_rep_section=i, r_vals=aa_cuboidal_zxtnt)
               CALL section_vals_val_get(aa_cuboidal_section, "v_D", i_rep_section=i, r_val=v_D)
               CALL section_vals_val_get(aa_cuboidal_section, "SIGMA", i_rep_section=i, r_val=sigma)
               params%dbc_params%aa_cuboidal_nprtn(:, i) = aa_cuboidal_nprtn
               params%dbc_params%aa_cuboidal_xxtnt(:, i) = aa_cuboidal_xxtnt
               params%dbc_params%aa_cuboidal_yxtnt(:, i) = aa_cuboidal_yxtnt
               params%dbc_params%aa_cuboidal_zxtnt(:, i) = aa_cuboidal_zxtnt
               params%dbc_params%aa_cuboidal_vD(i) = v_D
               params%dbc_params%aa_cuboidal_sigma(i) = sigma
            END DO
         ELSE
            params%dbc_params%n_aa_cuboidal = 0
         END IF

      END IF

   END SUBROUTINE dirichlet_bc_read_parameters

END MODULE pw_poisson_read_input
