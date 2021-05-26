#!/usr/bin/env python
# coding: utf-8

# # <center>Day 1, Session 1: Molecules</center>

# <center>MDAnalysis Workshop - May 2021 - PRACE/SURF</center>

# ## Getting started with MDAnalysis

# **Installing MDAnalysis**
# 
# We've done the hard work for you of installing MDAnalysis into the JupyterHub environment.
# 
# If you want to install on your own computer, this is normally done through **pip** or **conda** (see https://www.mdanalysis.org/pages/installation_quick_start/). Note though that in this workshop, we're using the new MDAnalysis 2.0.0 beta version, which is installed using pip by:
# 
# `pip install MDAnalysis==2.0.0b0`
# 
# If you want to use the example data used here, you'll also need MDAnalysisTests:
# 
# `pip isntall MDAnalysis==2.0.0b0`

# ### The general structure of MDAnalysis
# 
# The two fundamental classes of MDAnalysis are the `Universe` and the `AtomGroup`.
# 

# 
# <center><img src="imgs/mdaclasses.png" alt="mda" style="width: 1000px;"/></center>

# - **The `Universe`** contains everything about a molecular dynamics system
#   - Static information: atoms and their connectivities
#   - Dynamic information: The trajectory
#   
# - The atoms in a `Universe` can be accessed through a hierarchy of containers:
#  - *Atoms* can be grouped together into **an `AtomGroup`** 
#     - *Residues* are made up of *atoms*. They can be grouped into `ResidueGroups`
#       - *Segments* are made up of *residues*. They can be grouped into `SegmentGroups`.

# **A (very) basic workflow for an analysis in MDAnalysis:**

# 1. import MDAnalysis
# 2. load a Universe
# 3. define an atomgroup
# 4. collect position data
# 5. analyse!

# # `Universe`

# **The basic command for loading a universe is:**
# 
#  `u = mda.Universe(topology, trajectory)`

# - The *topology* file must contain the atom information 
# - The (optional) *trajectory* file(s) contains the positions of atoms with time - more on this next Session. 
# 
# Note that some files can double as both a *topology* and a *trajectory*.  
# 
# MDanalysis supports [over 40 input file types](https://userguide.mdanalysis.org/2.0.0-dev0/formats/index.html#formats)

# In[1]:


# First we import MDAnalysis
import MDAnalysis as mda

# Let's get some example data
from MDAnalysis.tests.datafiles import PSF, DCD

# and now load our universe!
u = mda.Universe(PSF, DCD)
print(u)


# **Key properties of a `Universe`:**

# - `atoms`: an `AtomGroup` containing all of the system's atoms
#     - similarly, `segments` and `residues`; a `SegmentGroup` and a `ResidueGroup`, respectively
#     
# - Various bond and angle information, as `TopologyGroups`: `bonds`, `angles`, `dihedrals`, `impropers` (if found in the topology file)
# 
# - `trajectory`: next Session!

# In[2]:


u.bonds


# # AtomGroups

# **An `AtomGroup` is an "array" of atoms.**
# 
# We can get various properties of each atom contained in an `AtomGroup` through attribues, e.g.:

#    - `names`
#    - `resnames`
#    - `resids` 
#    - `charges`
#    - `masses`
# 
# Exactly which properties you can get depend on what is read from the topology (see the [documentation](https://userguide.mdanalysis.org/2.0.0-dev0/formats/index.html#formats))

# In[3]:


ag = u.atoms
ag.names


# **`ResidueGroup`s and `SegmentGroup`s work similarly.**
# 

# We can get various properties of each atom contained in it through attribues, using `atoms`, `residues` and `segments`

# In[4]:


rg = u.residues
rg.resids


# We don't ususally want to work with the whole set of atoms in a trajectory. We need a way to create `AtomGroups` containing selected atoms...

# # Atom selections

# ### But first... visualising selections

# [nglview](https://github.com/nglviewer/nglview#usage) will allow us to view MDAnalysis Universes and AtomGroups inside Jupyter notebooks.
# 
# We've intalled if for you on the workshop Jupyter hub, but you can also install it on your own computer (*Beware: it's not always straightforward getting it to actually work. Try making sure to install nglview and launch the notebook from the same conda environment*)

# In[5]:


# first, import nglview
import nglview as nv
  
# add a universe (or atomgroup)
view_u = nv.show_mdanalysis(u)

# launch the viewer
view_u


# ## Selecting atoms to create AtomGroups

#  - indexing

# In[6]:


ag = u.atoms[0:2]


# In[7]:


view = nv.show_mdanalysis(ag)
view


# **Selection strings and `select_atoms`**
# 
# We can use the `select_atoms()` method of an `AtomGroup` or `Universe` to return an `AtomGroup` based on a selection string.

# There's a lot of options for selection strings (see the  [UserGuide]( https://userguide.mdanalysis.org/2.0.0-dev0/selections.html)); including:
# 
#  - selection by attribute (e.g. residue name (`resname`)), including presets like `protein`
#  - wildcard matching (`*`)
#  - boolean operators (`and`, `or`, `not`)
#  - geometric (e.g. `around`, `sphzone`, ...)
#  - and more!
#  
#  

# In[8]:


ag = u.select_atoms('protein')
view_ag = nv.show_mdanalysis(ag)
view_ag


# In[9]:


view_ag.add_licorice()


# # Working with coordinates

# **The most useful attribute of out atoms are their coordinates, available in the `positions` attribute of an `AtomGroup`**

# The positions are returned as a numpy array, which we can then readily manipulate.
# 
# There are some built-in functions based on position data, e.g. `center_of_mass()`, `center_of_geometry()`

# In[10]:


pos = u.atoms.positions
print(pos)


# This is just data from one frame - in the next session, you'll learn how to work with trajectories to get data across a whole simulation.

# ## A summary of Lecture 1
# 
# Most simulation analysis will involve extracting position data from certain atoms.
# 
# - A `Universe` contains all information about a simulation system
# 
# - An `AtomGroup` contains information about a group of atoms
# 
# - We can use `Universe.select_atoms()` to create an `AtomGroup` containing specific atoms from a `Universe`
# 
# - Positions of atoms in an AtomGroup are accessed through `AtomGroup.positions`

# ### Now - on to the first tutorial!
# 
# Find the tutorial notebook at https://jupyter.lisa.surfsara.nl/jhlsrf005/hub/ - look for Day1-Session1-Practical/session1_practical.ipynb
# 
# **Remember:**
# - Go at your own pace!
# - Ask questions!
# - Be respectful!
# - Take breaks!

# In[ ]:




