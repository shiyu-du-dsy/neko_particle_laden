import numpy as np

from helix_tool import plot_particles

# particle output written by neko: tstep,time,particle_id,x,y,z,u,v,w
filename = 'data/test/tracers_5.csv'
data = np.genfromtxt(filename, delimiter=',', names=True)

tsteps = np.unique(data['tstep'])
print(f'{filename}: {len(data)} rows, {len(tsteps)} snapshots, '
      f'time steps {tsteps}')


def snapshot(tstep):
    return data[data['tstep'] == tstep]


# 3D view, top view and cross-section of the last written snapshot
tstep_choosen = -1
snap = snapshot(tsteps[tstep_choosen])
print(f'3D plot: tstep {tsteps[tstep_choosen]:g} (t = {snap["time"][0]:g}), '
      f'{len(snap)} particles')
plot_particles(snap['x'], snap['y'], snap['z'],
               filename='tracers.png', color=snap['u'], color_label='u (axial)')
