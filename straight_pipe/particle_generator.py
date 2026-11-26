import numpy as np

from helix_tool import straight_to_helix, plot_particles

#### parameter setup for the helix
curvature = 0.0
torsion = 0.0

#### parameter setup for the cross section
radius_cs = 1.0
l_axial = 25.0

#### number of particles in each direction
n_theta = 16
n_r = 5
n_axial = 35

#### setup inertial particles
d_p_value = 0.00131750*2.0 # particle diameter
rho_p_value = 1.0e+3 # particle density

print("total number of particles: ", \
      n_theta * n_r * n_axial - n_axial*(n_theta - 1)) # subtract the duplicates on the centreline

###################################################################
#### first setup the particles in a straight pipe
theta_1d = np.linspace(0, 2 * np.pi, n_theta, endpoint=False) # azimuthal coordinate
h_nearwall = 1e-2 # the first point close to the wall
r_1d = np.linspace(0, radius_cs *(1 - h_nearwall), n_r) # radial coordinate
x_axial_1d = np.linspace(0, l_axial, n_axial, endpoint=False) # axial coordinate

theta, r, x_axial = np.meshgrid(theta_1d, r_1d, x_axial_1d, indexing='ij')

# at r = 0 every theta maps to the same point on the centreline:
# keep only the first theta there to remove the duplicates.
# (the arrays become 1D lists of particles after this masking)
keep = (r > 0) | (theta == theta_1d[0])
theta, r, x_axial = theta[keep], r[keep], x_axial[keep]

# local Cartesian coordinates of the cross-section (straight pipe):
# y_cs horizontal, z_cs vertical, theta measured from +y_cs towards +z_cs
y_cs = r * np.cos(theta)
z_cs = r * np.sin(theta)

####################################################################
#### given the curvature and torsion, map the particles to a helical pipe
# x, y, z = straight_to_helix(x_axial, y_cs, z_cs, curvature, torsion)

##### test with straight pipe
x, y, z = x_axial, y_cs, z_cs

plot_particles(x, y, z, filename='particles_ic.png')

u = np.zeros_like(x)
v = np.zeros_like(y)
w = np.zeros_like(z)
d_p = np.ones_like(x) # particle diameter
d_p *= d_p_value
rho_p = np.ones_like(x) # particle density
rho_p *= rho_p_value

data = np.stack([x, y, z, u, v, w, d_p, rho_p], axis=-1)
np.savetxt("particles_ic.csv", data, delimiter=",")
