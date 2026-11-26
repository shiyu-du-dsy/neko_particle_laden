"""Helper functions for mapping a straight pipe onto a helical pipe."""
import numpy as np


def helix_frame(s, curvature, torsion):
    """Centreline point c and Frenet frame (t, n, b) at arc length s.

    Same helix as turb_pipe_curved_z.f90: the helix axis is vertical (+z)
    and passes through (r_c, 0, 0), the coils lie in x-y planes, the pitch
    advances along +z and the centreline starts at the origin (s = 0).
    Returns arrays of shape (3, *s.shape).
    """
    k2 = curvature**2 + torsion**2
    r_c = curvature / k2          # helix radius
    p_h = torsion / k2            # pitch per radian
    l_h = 1.0 / np.sqrt(k2)       # arc length per radian
    phi = np.asarray(s) / l_h     # helix angle
    cos, sin = np.cos(phi), np.sin(phi)
    zero, one = np.zeros_like(phi), np.ones_like(phi)

    c = np.array([r_c * (1 - cos), -r_c * sin, p_h * phi])
    t = np.array([r_c * sin, -r_c * cos, p_h * one]) / l_h
    n = np.array([cos, sin, zero])                    # towards the helix axis
    b = np.array([-p_h * sin, p_h * cos, r_c * one]) / l_h  # mostly vertical
    return c, t, n, b


def straight_to_helix(x_axial, y_cs, z_cs, curvature, torsion):
    """Map straight-pipe coordinates onto the helical pipe.

    The axial coordinate is the arc length along the centreline, the
    horizontal offset y_cs goes along the normal (towards the helix axis)
    and the vertical offset z_cs along the binormal, so (x, y, z) of the
    straight pipe maps to (t, n, b) and stays right-handed.
    """
    c, t, n, b = helix_frame(x_axial, curvature, torsion)
    return c + y_cs * n + z_cs * b


def pipe_wall(pipe_length, curvature=0.0, torsion=0.0, n_s=100, n_theta=40):
    """Points on the pipe wall (radius 1) for plotting, shape (n_s, n_theta)
    each. Straight pipe along x if curvature == 0, otherwise the helical
    pipe of straight_to_helix."""
    s, theta = np.meshgrid(np.linspace(0, pipe_length, n_s),
                           np.linspace(0, 2 * np.pi, n_theta), indexing='ij')
    y_cs, z_cs = np.cos(theta), np.sin(theta)
    if curvature == 0.0:
        return s, y_cs, z_cs
    return straight_to_helix(s, y_cs, z_cs, curvature, torsion)


def plot_particles(x, y, z, filename=None, color=None, color_label=None,
                   pipe_length=25.0, curvature=0.0, torsion=0.0):
    """Plot the particle positions as three separate figures: a 3D view, a
    top view (x-y plane) and the cross-section (y-z plane, all particles
    projected along x).

    filename: base name of the output, e.g. 'tracers.png' is saved as
    tracers_3d.png, tracers_top.png and tracers_cross_section.png.
    If None, the figures are shown instead.
    color: optional value per particle (e.g. a velocity component) used to
    colour the markers, with color_label as the colour-bar label.
    pipe_length, curvature, torsion: pipe geometry drawn in the 3D view
    (straight pipe along x if curvature == 0, no pipe if pipe_length is None).
    Returns the three figures.
    """
    import os
    import matplotlib.pyplot as plt

    def finish(fig, ax, sc, suffix, pad=0.05):
        if color is not None:
            fig.colorbar(sc, ax=ax, shrink=0.8, pad=pad, label=color_label)
        if filename is not None:
            base, ext = os.path.splitext(filename)
            fig.savefig(f'{base}_{suffix}{ext or ".png"}', dpi=150,
                        bbox_inches='tight')

    # 3D view with equal scaling
    fig3d = plt.figure(figsize=(12, 5))
    ax3d = fig3d.add_subplot(projection='3d')
    sc = ax3d.scatter(x, y, z, s=1, c=color)
    ax3d.set_xlabel('x')
    ax3d.set_ylabel('y')
    ax3d.set_zlabel('z (vertical)')
    if pipe_length is not None:
        xw, yw, zw = pipe_wall(pipe_length, curvature, torsion)
        ax3d.plot_surface(xw, yw, zw, color='grey', alpha=0.15,
                          linewidth=0, shade=False)
        for i in (0, -1):  # inlet and outlet rings
            ax3d.plot(xw[i], yw[i], zw[i], color='k', lw=0.8)
        xa = np.concatenate([np.ravel(x), xw.ravel()])
        ya = np.concatenate([np.ravel(y), yw.ravel()])
        za = np.concatenate([np.ravel(z), zw.ravel()])
    else:
        xa, ya, za = x, y, z
    ax3d.set_box_aspect((np.ptp(xa), np.ptp(ya), np.ptp(za)))
    ax3d.grid(False)
    if curvature == 0.0:  # straight pipe: y and z span only the pipe radius
        ax3d.set_yticks([-1, 0, 1])
        ax3d.set_zticks([-1, 0, 1])
    finish(fig3d, ax3d, sc, '3d', pad=0.15)  # room for the z label

    # top view: shows the curvature of the centreline
    figtop, axtop = plt.subplots(figsize=(12, 4))
    sc = axtop.scatter(x, y, s=2, c=color)
    axtop.set_xlabel('x')
    axtop.set_ylabel('y')
    axtop.set_title('top view')
    axtop.set_aspect('equal')
    finish(figtop, axtop, sc, 'top')

    # cross-section: y-z plane, viewed along the pipe axis
    figcs, axcs = plt.subplots(figsize=(6, 5))
    sc = axcs.scatter(y, z, s=4, c=color)
    axcs.add_patch(plt.Circle((0, 0), 1, fill=False, color='k'))  # pipe wall
    axcs.set_xlim(-1.1, 1.1)
    axcs.set_ylim(-1.1, 1.1)
    axcs.set_xticks([-1, 0, 1])
    axcs.set_yticks([-1, 0, 1])
    axcs.set_xlabel('y')
    axcs.set_ylabel('z (vertical)')
    axcs.set_title('cross-section')
    axcs.set_aspect('equal')
    finish(figcs, axcs, sc, 'cross_section')

    if filename is None:
        plt.show()
    return fig3d, figtop, figcs
