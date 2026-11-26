module user
  use neko
  implicit none

  ! The mesh (turb_pipe.nmsh) has its pipe axis along z. It is rotated so
  ! that x is the axial direction and z the vertical direction, using the
  ! cyclic (right-handed) permutation of the mesh coordinates
  !   x = z_mesh (axial), y = x_mesh (horizontal), z = y_mesh (vertical).

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%mesh_setup => user_mesh_setup
  end subroutine user_setup

  ! Rotate the mesh right after it is read, before the dofmap is built: the
  ! vertices and the midpoints of the curved (second-order) edges are
  ! permuted, so the GLL points and everything that uses the mesh vertices
  ! (e.g. the periodic shift of the Lagrangian particles) see the rotated
  ! pipe. Connectivity and periodicity are defined by ids, not coordinates.
  subroutine user_mesh_setup(msh, time)
    type(mesh_t), intent(inout) :: msh
    type(time_state_t), intent(in) :: time
    integer :: i, j

    do i = 1, msh%mpts
       call rotate_pipe_to_x(msh%points(i)%x)
    end do

    do i = 1, msh%curve%size
       do j = 1, 12
          select case (msh%curve%curve_el(i)%curve_type(j))
          case (0)
          case (4) ! edge midpoint
             call rotate_pipe_to_x(msh%curve%curve_el(i)%curve_data(1:3, j))
          case default
             call neko_error('turb_pipe: unsupported curve type in the mesh')
          end select
       end do
    end do

    call neko_log%message('Pipe rotated: x axial, z vertical')
  end subroutine user_mesh_setup

  ! (x, y, z)_mesh -> (z, x, y)_mesh: x = axial, y = horizontal, z = vertical
  subroutine rotate_pipe_to_x(p)
    real(kind=dp), intent(inout) :: p(3)
    p = [p(3), p(1), p(2)]
  end subroutine rotate_pipe_to_x

  ! User defined initial condition
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    integer :: i
    real(kind=rp) :: uvw(3), x, y, z
    type (field_t), pointer :: u, v, w, p

    u => fields%get("u")
    v => fields%get("v")
    w => fields%get("w")
    p => fields%get("p")

    ! Initialize pressure field to zero
    call field_rzero(p)

    do i = 1, u%size()
       x = u%dof%x%x(i,1,1,1)
       y = u%dof%y%x(i,1,1,1)
       z = u%dof%z%x(i,1,1,1)
       ! pipe_ic works in the original mesh frame (axis along z):
       ! pass (y, z, x) and rotate the result back (axial -> u)
       uvw = pipe_ic(y, z, x)

       u%x(i,1,1,1) = uvw(3)
       v%x(i,1,1,1) = uvw(1)
       w%x(i,1,1,1) = uvw(2)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call u%copy_from(HOST_TO_DEVICE, sync=.false.)
       call v%copy_from(HOST_TO_DEVICE, sync=.false.)
       call w%copy_from(HOST_TO_DEVICE, sync=.true.)
    end if
  end subroutine initial_conditions

  function pipe_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: rand, ux, uy, uz, xr, yr, rr, zo
    real(kind=rp) :: amp_z, freq_z, freq_t, amp_tht, amp_clip, blt
    real(kind=rp) :: phase_z, arg_tht, amp_sin, pi, th

    pi = 4_rp * atan(1.0_rp)
    xr = x
    yr = y
    rr = xr*xr + yr*yr
    if (rr .gt. 0.0_rp) rr = sqrt(rr)
    th = atan2(y,x)
    zo = 2.0_rp * pi * z / 25.0_rp

    uz = 6.0_rp * (1.0_rp - rr**6.0_rp) / 5.0_rp

    ! Assign a wiggly shear layer near the wall
    amp_z = 35e-2_rp ! Fraction of 2pi for z-based phase modification
    freq_z = 5.0_rp ! Number of wiggles in axial- (z-) direction
    freq_t = 6.0_rp ! Frequency of wiggles in azimuthal-direction

    amp_tht = 10.0_rp ! Amplification factor for clipped sine function
    amp_clip = 4e-1_rp ! Clipped amplitude

    blt = 3.5e-1_rp ! Fraction of boundary layer with momentum deficit

    phase_z = amp_z * (2.0_rp * pi) * sin(freq_z * zo)

    arg_tht = freq_t * th + phase_z
    amp_sin = 5.0_rp * sin(arg_tht)
    if (amp_sin .gt. amp_clip) amp_sin = amp_clip
    if (amp_sin .lt. -amp_clip) amp_sin = -amp_clip

    if (rr .gt. (1-blt)) uz = uz + amp_sin
    call random_number(rand)

    ux = 5e-2_rp * rand * rand
    uy = 1e-1_rp * rand * rand * rand
    uz = uz + 1e-2_rp * rand

    uvw(1) = ux
    uvw(2) = uy
    uvw(3) = uz

  end function pipe_ic
end module user
