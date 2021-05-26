#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import warnings
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
warnings.filterwarnings("ignore")


# # <center>Lecture 2: Handling Trajectory Data</center>

# <center>MDAnalysis workshop - PRACE/SURF 2021</center>

# ## Aims and objectives <a class="tocSkip">

# - Give a general overview of how MDAnalysis reads coordinate data
# - How to programmatically access / use trajectory data with the MDAnalysis API
# - How to write out coordinates

# ## Overview of File IO in MDAnalysis <a class="tocSkip">

# - Main file IO in MDAnalysis
#   - Topologies
#   - Coordinates (single frame & trajectories)
#   - Auxilliary (session 6)
#   
# - Topology & coordinates form the basis of Universe creation

# <center><img src="imgs/MDA-overview.png" alt="MDA-overview" style="width: 700px;"/></center>

# ## Topology files <a class="tocSkip">

# - 27 formats supported
# - Provides base information about the system
#   - e.g. atom names, elements, charges, bonds
# - MDAnalysis **does not** parse energy function related terms
# - Read only
#   - Except for topologies that contain coordinate data (PDB, GRO, etc...)
# - Only read once at Universe creation

# In[ ]:


import MDAnalysis as mda
from MDAnalysis.tests.datafiles import TPR

u = mda.Universe(TPR)

print(u.atoms.names)

print(u.atoms.charges)


# ## Coordinates files <a class="tocSkip">

# - 28 coordinate formats supported
#   - Extra file formats via [chemfiles](https://github.com/chemfiles/chemfiles.py) converter interface
# - Reads both single frame (e.g. MMTF, PQR) and trajectory files (e.g. XTC, DCD, NetCDF)
# - MDAnalysis can both read & write
#   - 15 / 28 coordinate file types have writers
# - Iterative IO model
#   - Data read from disk frame at a time as requested

# <center><img src="imgs/mda-coordinates.png" alt="MDA-orview-of-coordinate-readers" style="width: 700px;"/></center>

# ## General overview of coordinate IO <a class="tocSkip">

# <center><img src="imgs/MDA-CoordinateReading.png" alt="MDA-coordinate-workflow" style="width: 1000px;"/></center>

# ## Trajectory class: `Universe.trajectory` <a class="tocSkip">

# - Main class for trajectory access
#   - Makes the Reader fetch the relevant data from disk
#   - Works on the basis of frames
# - Works as an iterator
# - Can be randomly seeked

# In[ ]:


from MDAnalysis.tests.datafiles import (PRM_NCBOX,
                                        TRJ_NCBOX)

u = mda.Universe(PRM_NCBOX, TRJ_NCBOX)

u.trajectory


# ## Trajectory class: `Universe.trajectory` <a class="tocSkip">

# - Main class for trajectory access
#   - Makes the Reader fetch the relevant data from disk
#   - Works on the basis of frames
#   - Time is in units of **picoseconds**
# - Works as an iterator
# - Can be randomly seeked

# In[ ]:


print("number of frames: ", u.trajectory.n_frames)

print("current frame: ", u.trajectory.frame)

print("current time: ", u.trajectory.time)


# ## Trajectory class: `Universe.trajectory` <a class="tocSkip">

# - Main class for trajectory access
#   - Makes the Reader fetch the relevant data from disk
#   - Works on the basis of frames
#   - Time is in units of **picoseconds**
# - Works as an iterator
# - Can be randomly seeked

# In[ ]:


# trajectories can be seeked
u.trajectory[5]

print("current frame: ", u.trajectory.frame)
print("current time: ", u.trajectory.time)


# ## Trajectory class: `Universe.trajectory` <a class="tocSkip">

# - Main class for trajectory access
#   - Makes the Reader fetch the relevant data from disk
#   - Works on the basis of frames
#   - Time is in units of **picoseconds**
# - Works as an iterator
# - Can be randomly seeked

# In[ ]:


# trajectories can be iterated
for ts in u.trajectory:
    print('current frame / time: ', u.trajectory.frame,
          u.trajectory.time)


# ## Timestep class `Universe.trajectory.ts` <a class="tocSkip">

# - Class that holds the current frame's data
#   - Trajectory tells Reader to populate Timestep with a given frame
#   - Changes every time the trajectory seeks
# - Contains time-dependent information
#   - `positions`, `dimensions`, `time`, `frame number`, etc..
# - Data is then fed to other parts of the Universe
#   - i.e. `AtomGroup.positions` takes a copy of `ts.positions` given atoms

# In[ ]:


u.trajectory.ts


# ## Timestep class `Universe.trajectory.ts` <a class="tocSkip">

# - Class that holds the current frame's data
#   - Trajectory tells Reader to populate Timestep with a given frame
#   - Changes every time the trajectory seeks
# - Contains time-dependent information
#   - `positions`, `dimensions`, `time`, `frame number`, etc..
# - Data is then fed to other parts of the Universe
#   - i.e. `AtomGroup.positions` takes a copy of `ts.positions` given atoms

# In[ ]:


# ts holds time-dependent data like positions and forces

print('positions:\n', u.trajectory.ts.positions)

print('forces:\n', u.trajectory.ts.forces)


# ## Timestep class `Universe.trajectory.ts` <a class="tocSkip">

# - Class that holds the current frame's data
#   - Trajectory tells Reader to populate Timestep with a given frame
#   - Changes every time the trajectory seeks
# - Contains time-dependent information
#   - `positions`, `dimensions`, `time`, `frame number`, etc..
# - Data is then fed to other parts of the Universe
#   - i.e. `AtomGroup.positions` takes a copy of `ts.positions` given atoms

# In[ ]:


# positions are passed to things like AtomGroups
print('position:\n', u.trajectory.ts.positions[0])

print('atomgroup position:\n', u.atoms[0].position)


# ## AtomGroups and UpdatingAtomGroups <a class="tocSkip">

# - Main entry point for data is AtomGroup
#   - API preferred over `ts`
# 
# - AtomGroups are statically defined
#   - Atoms don't change
#   - Only time-dependent variables change
#   
# - UpdatingAtomGroups
#   - Time dependent selection of atoms
#   - Expensive, but can be useful

# In[ ]:


u.trajectory[0]

# atoms with y coordinates > 1 Angstrom
ag = u.select_atoms(
    'prop y > 1',
    updating=True)

print('atoms: ', ag.n_atoms)

u.trajectory[5]

print('atoms: ', ag.n_atoms)


# ## Example: TC5b Trp Cage <a class="tocSkip">

# - Implicit solvent simulation of TC5b
#   - Makes different folds over time
# 
# 
# - Can use NGLView to visualise the trajectory

# In[ ]:


import nglview as nv

u = mda.Universe('data/TC5b.ff14igb8.HMR.parm7',
                 'data/TC5b.nc')

nv.show_mdanalysis(u)


# ## Example: TC5b Trp Cage <a class="tocSkip">

# - Let's calculate a simple end-to-end alpha-carbon distance over the trajectory
#   - Very coarse evaluation of folding

# In[ ]:


from MDAnalysis.lib.distances import calc_bonds

nterm = u.select_atoms('name CA')[0]
cterm = u.select_atoms('name CA')[-1]

times = np.empty(u.trajectory.n_frames)
distances = np.empty(u.trajectory.n_frames)

for i, ts in enumerate(u.trajectory):
    times[i] = u.trajectory.time / 1000
    distances[i] = calc_bonds(nterm.position,
                              cterm.position)
    
plt.plot(times, distances)
plt.xlabel('time (ns)')
plt.ylabel('distance (Angstroms)')


# ## I/O Performance Considerations <a class="tocSkip">

# - I/O is expensive - minimize wherever possible
# - Random seeking can be much more expensive than iterating
# - Binary are faster to read than ASCII
# - Lots of concurrent I/O calls can cause issues
#   - See parallelism session

# ## MemoryReader <a class="tocSkip">

# - Transferring data to memory speeds up access to trajectory data
#   - ~ 1-2 orders of magnitude speedups
#   - Requires being able to store data in RAM
# 
# - Changes to trajectory data are "permanent"
#   - Directly changes the memory object
#   - Can be useful
#   - example: trajectory aligning (`mda.align.AlignTraj`)
#   
# - Converters are usually MemoryReader instances
#   - See session 6

# ## Writing coordinates <a class="tocSkip">

# - MDAnalysis implements writers for most coordinate types
# - Can directly write from `AtomGroup` or use `MDAnalysis.Writer` class
# - Writers will automatically guess file type based on extension
#   - Can also pass a `format` keyword to enfore the file type

# In[ ]:


# AtomGroup of all atoms
ag = u.atoms

# write out a single frame
ag.write('test.pdb')

# write out all frames
ag.write('test.xtc', frame=all)


# ## Writing coordinates <a class="tocSkip">

# - MDAnalysis implements writers for most coordinate types
# - Can directly write from `AtomGroup` or use `MDAnalysis.Writer` class
# - Writers will automatically guess file type based on extension
#   - Can also pass a `format` keyword to enfore the file type

# In[ ]:


# Common to use mda.Writer as a context manager
with mda.Writer('test.dcd', n_atoms=u.atoms.n_atoms,
                format='DCD') as W:
    # Iterate through trajectory to feed the writer
    # in this case - downsample by skipping by 10
    for ts in u.trajectory[::10]:
        W.write(u)


# ## Review <a class="tocSkip">

# - Universes are associated with a Reader when loading coordinates
# - `Universe.trajectory` allows for iterating through a trajectory
# - `Universe.trajectory.ts` is assigned time-dependent data which is passed to other parts of the Universe
# - Writing is achieved through either `AtomGroup.write()` or `MDAnalysis.Write`
# - I/O is expensive, transferring to memory can be helpful, but generally try to limit I/O calls

# ## Session 2 - Practical <a class="tocSkip">

# - Notebook under Day1-Session2-Practical/session2.ipynb
#   - Covers concepts described here
#   - Also gives case study examples of how to iterate through trajectories to carry out analyses
