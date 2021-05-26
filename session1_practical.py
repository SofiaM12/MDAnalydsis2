#!/usr/bin/env python
# coding: utf-8

# # Day 1, Practical 1
# ## Introducting MDAnalysis and working with molecules

# **This notebook introduces the fundamental objects of `MDAnalysis`:**
# 
#  - the `Universe`,
#  - `AtomGroup`, and `Atom` objects.
#  - selecting and manipulating `AtomGroup`s
#  - `Residue`s and `Segment`s
#  - `Bond` `Angle` and `Dihedral`
# 
# 
# **Additional resources**
#  - During the workshop, feel free to ask questions at any time
#  - For more on how to use MDAnalysis, see the [User Guide](https://userguide.mdanalysis.org/2.0.0-dev0/) and [documentation](https://docs.mdanalysis.org/2.0.0-dev0/)
#  - Ask questions on the [user mailing list](https://groups.google.com/group/mdnalysis-discussion) or on [Discord](https://discord.gg/fXTSfDJyxE)
#  - Report bugs on [GitHub](https://github.com/MDAnalysis/mdanalysis/issues?)

# # 1. Fundamental MDAnalysis objects
# 
# ## Universe loading 101
# 
# > "If you wish to make an apple pie from scratch, you must first invent the Universe." 
# 
# > ~ Carl Sagan
# 
# First, we need to import `MDAnalysis`, giving us access to all the components in its namespace:

# In[2]:


import MDAnalysis as mda


# One of the most fundamental objects in the `MDAnalysis` data model is the `Universe` object.
# A `Universe` can be thought of as an interface to all the data of a simulation;
# it contains all of a simulations' topology information (names, charges, masses etc) at the least,
# but usually also includes trajectory information (positions, velocities etc) as well.

# In order to do anything, we do need some actual molecular dynamics data to work with. Let's load an example (the protein adenylate kinase, or AdK) from the MDAnalysis tests data:

# In[3]:


from MDAnalysis.tests.datafiles import PSF, DCD


# To make a `Universe`, we need at the very least a topology file - see the [topology readers](https://docs.mdanalysis.org/2.0.0-dev0/documentation_pages/topology/init.html) documentation for a list of supported file formats. 
# 
# Since the type of topology file we're using in this example (a PSF file) doesn't contain coordinates, we'll also need to load a trajectory file (in this case a DCD file) so we have some position data to work with later. You'll learn more about working with trajectories in the next session. 
# 

# In[4]:


u = mda.Universe(PSF, DCD)

print(u)


# ## AtomGroups

# We now have a `Universe` object. Since the topology (PSF) file we loaded contained both atom identities and bond information, the `Universe` is able to access these details.
# 
# 

# We can access all atoms in the `Universe` through the `Universe.atoms` attribute.
# This returns an `AtomGroup`, which is probably the most important class we will learn about.

# In[5]:


ag = u.atoms
type(ag)


# An `AtomGroup` is like an array of atoms, and offers access to the data of these atoms through various attributes:

# In[6]:


ag.indices


# In[7]:


ag.names


# In[8]:


ag.resnames


# In[9]:


ag.resids


# In[10]:


ag.charges


# In[11]:


ag.masses


# In[12]:


ag.types


# All of these attributes of an `AtomGroup` return numpy arrays of the same length as the `AtomGroup` itself;
# that is, each element corresponds to each atom in the `AtomGroup`, in order.

# In[13]:


print(ag.n_atoms)
print(len(ag.names))


# In general, `MDAnalysis` will try and extract as much information as possible from the files given to `Universe`. The [topology readers](https://docs.mdanalysis.org/2.0.0-dev0/documentation_pages/topology/init.html) documentation lists the attributes that are extracted from each filetype.

# ## Working with individual atoms
# 
# By slicing an `AtomGroup` we can access individual `Atom` objects.
# These `Atom` objects will have singular versions of the various attributes of `AtomGroup`s.
# 
# In general working with individual `Atom` objects is discouraged as it is inefficient and will lead to poor performance.

# In[14]:


u = mda.Universe(PSF, DCD)
a = u.atoms[0]
print(a)


# In[15]:


print("name:", a.name)
print("resid:", a.resid)
print("resname:", a.resname)


# ## ResidueGroups and SegmentGroups

# The `Universe` also gives higher-order topology objects, including `ResidueGroups` and `SegmentGroups`. We can access all residues in the `Universe` with:

# In[16]:


u.residues


# And all segments with:

# In[17]:


u.segments


# `ResidueGroups` and `SegmentGroups` also behave similarly to `AtomGroups`, with many of their methods returning `numpy` arrays with each element corresponding to a single residue or segment, respectively.

# In[18]:


u.residues.resnames


# In[19]:


u.segments.segids


# You can also create a `ResidueGroup` from an `AtomGroup`:

# In[20]:


ag.residues


# ## Exercise 1
# 
# **1a. Load the GRO topology file from `MDAnalysis.tests.datafiles` and count how many atoms, residues and segments it contains.**

# In[21]:


from MDAnalysis.tests.datafiles import GRO
u_gro = mda.Universe(GRO)
print('Atoms: ', u_gro.atoms.n_atoms)
print(len(u_gro.atoms))
print('Residues: ', u_gro.residues.n_residues)
print('Segments: ', u_gro.segments.n_segments)


# In[30]:


print('Segment 1: ', u_gro.segments[0].segid)
print('last at: ', u_gro.atoms[-1].name)
print('10th residue: ', u_gro.residues[10].resname)


# **1b. From the above universe, find the name of the:**
#   - first segment
#   - last atom
#   - 10th residue.

# In[22]:


print('First segment:', u_gro.segments[0].segid)
print('Last atom: ', u_gro.atoms[-1].name)
print('10th residue: ', u_gro.residues[10].resname)


# In[23]:


# Exercise 1b


# # 2. Selecting atoms
# 
# It is also rare that we want to operate on all atoms in the system!
# 
# `MDAnalysis` offers a few different ways to select atoms,
# in this section we will go over the most useful methods.

# ## Numpy style selections
# 
# As previously mentioned, an `AtomGroup` is like an array of atoms,
# and therefore we can slice it exactly like we would slice a `numpy` array.
# 
# ### Fancy indexing
# 
# The simplest option to select specific atom is to use fancy indexing. You can specify the atoms in a list

# In[24]:


u = mda.Universe(PSF, DCD)
u.atoms[[1, 4, 5 , 0]]


# or as a range

# In[25]:


u.atoms[1:10]


# ### Boolean indexing
# 
# You can also create an array with `True`/`False` values of the same length as the `AtomGroup`. Every atom for which the array is set to `True` will be selected.
# 
# We can create such a boolean array by doing a comparison of a `numpy` array:

# In[26]:


selection_ar = u.atoms.resnames == 'GLY'
print("selection array = ", selection_ar)
u.atoms[selection_ar]


# or shorter:

# In[27]:


u.atoms[u.atoms.resnames == 'GLY']


# We can also do this with `ResidueGroup`s and `SegmentGroup`s, e.g.:

# In[28]:


u.residues[u.residues.resnames == 'GLY']


# ## Selection Strings and `select_atoms`
# 
# We've already seen that complex selections can be performed on `AtomGroups` using numpy style indexing.
# However, `MDAnalysis` also features a CHARMM-style atom selection mechanism that is often more convenient.
# We can consult the docstring for `select_atoms` to see the available selection keywords:

# In[29]:


get_ipython().run_line_magic('pinfo', 'u.atoms.select_atoms')


# You can also see them in the docs [here](https://docs.mdanalysis.org/2.0.0-dev0/documentation_pages/selections.html).
# 
# Although boolean selections work well enough for selecting out atoms from AtomGroups, the selection language makes more complex selections possible with probably less effort.
# 
# For example, we can select all glycines by:

# In[30]:


glycines = u.select_atoms("resname GLY")
glycines


# If selecting by residue number, we can specify a range:

# In[31]:


first10resids = u.select_atoms("resid 1-10")
first10resids


# We can use `and`/`or`/`not` with [logical conjunctions](https://en.wikipedia.org/wiki/Logical_conjunction):

# In[32]:


acidic = u.select_atoms("resname GLU or resname ASP")
acidic


# For name like selections, we can also Unix shell-style wildcards such as `*`. Here for example `name OD*` would select atoms named `OD1, OD2, OD3` etc:

# In[33]:


acidic_o = acidic.select_atoms('name OD* or name OE*')
acidic_o


# As a shortcut, multiple values can be given and these will be implicitly OR'd together.
# To select all atoms with name NZ or NH* in residues named LYS or ARG:

# In[34]:


basic_n = u.select_atoms("(resname LYS ARG) and (name NZ NH*)")
basic_n


# There are also several preset keywords for useful selections such as `backbone`, which selects all CA, C, O and N atoms:

# In[35]:


backbone_1 = u.select_atoms('backbone')
backbone_2 = u.select_atoms('name CA C O N')
backbone_1 == backbone_2


# ### Geometric selections
# 
# The `select_atoms` method also has various geometric keywords that make selecting atoms based on geometric criteria much easier.
# 
# For example, we can look for salt bridges by using the `around` selection operator to specify only atoms within 4 angstroms of a particular selection. Note we can simplify selection strings by referring to previous `AtomGroups`.

# In[36]:


acidic = u.select_atoms("group acidic and around 4 group basic", acidic=acidic_o, basic=basic_n)
acidic


# You can also select atoms based on absolute position using `prop`, e.g.:

# In[37]:


upper_z = u.select_atoms('prop z > 10')
upper_z


# ## Set operations
# 
# `AtomGroup`s can also be combined using `&` for "and" and `|` for "or":

# In[38]:


heavy = u.atoms[u.atoms.masses > 10]
arg = u.atoms[u.atoms.resnames == 'ARG']
print(len(heavy), len(arg))
print(len(heavy & arg))


# These two selections are identical:

# In[39]:


ag1 = u.atoms.select_atoms('resname LYS ARG')
ag2 = u.atoms.select_atoms('resname LYS') | u.atoms.select_atoms('resname ARG')

ag1 == ag2


# We can also concatenate atomgroups with `+` and subtract with `-`; the following two selections are also identical: 

# In[40]:


no_H1 = u.atoms.select_atoms('resname LYS ARG and not name H*')
no_H2 = u.atoms.select_atoms('resname LYS ARG') - u.atoms.select_atoms('name H*')

no_H1 == no_H2


# By design, an `AtomGroup` can have repeats of the same atom, for example through this selection:

# In[41]:


ag = u.atoms[[0, 0, 1, 2, 4, 4, 5]]
ag


# The `unique` property will return a version of the `AtomGroup` with only one of each Atom:

# In[42]:


ag.unique


# ## Exercise 2
# 
# **2a. Select residues 100 to 200, first using indexing and then using a selection string, and confirm you get the same selection.**

# In[43]:


ag1 = u.residues[99:200]
ag2 = u.select_atoms("resid 100-200").residues
ag1 == ag2


# In[44]:


# Exercise 2a


# What might be some reasons to prefer one method over the other?

# **2b. Count the number of arginine residues**

# In[45]:


ag = u.select_atoms('resname ARG')
print(len(ag.residues))


# In[46]:


# Exercise 2b


# **2c. Select all nitrogen atoms within 5.0A of a carbon atom**

# In[47]:


ag = u.select_atoms('name N* and around 5 name CA')
print(ag)


# In[48]:


# Exercise 2c


# **2d. Select all hydrogens that are bonded to an alpha carbon**
# 
# *Hint: Look through the `select_atoms` docstring above or [here](https://docs.mdanalysis.org/2.0.0-dev0/documentation_pages/selections.html) for keywords that might help!*

# In[49]:


ag = u.select_atoms('name H* and bonded name CA')
print(ag)


# In[50]:


# Exercise 2d


# # 3. Visualization with nglview

# [nglview](https://github.com/nglviewer/nglview#usage) is a package that allows to visualize a MDAnalysis `Universe` or `AtomGroup` directly on the jupyter notebook. First we load the structure:

# In[51]:


import nglview as nv

view_u = nv.show_mdanalysis(u)
view_u.camera = 'orthographic'


# And now we can launch the viewer:

# In[52]:


view_u


# Have a play around - you can rotate, zoom (and play through the trajectory) to view the structure. You can also select atoms to display their information.

# nglview can be useful to check that we are indeed selecting the correct `AtomGroup`:

# In[53]:


helix_1 = u.select_atoms('resid 12-25')
view_ag = nv.show_mdanalysis(helix_1)

view_ag


# To add another selection to the same view:

# In[54]:


resid_20 = u.select_atoms('resid 20')
view_ag.add_trajectory(resid_20)


# ### Exercise 3

# Show all atoms that are below the plane x=4.0

# In[55]:


ag = u.select_atoms('prop x < 4')
view_ex3 = nv.show_mdanalysis(ag)
view_ex3


# In[56]:


# Exercise 3


# # 4. Accessing coordinates
# 
# The most important attribute of your atoms is undoubtedly their positions! The position information is made available via an `AtomGroup` in the `positions` attribute:

# In[57]:


Ca = u.select_atoms('type CA')

Ca.positions


# This returns a `numpy` array, which can be easily manipulated. For example, we can could calulate the center of geometry by:

# In[58]:


import numpy as np
ca = u.select_atoms('name CA')
cog = np.mean(ca.positions, axis=0)
print(cog)


# Alternatively, we can use the built-in method:

# In[59]:


ca.center_of_geometry()


# Other convenient methods for common calculations based on positions include `center_of_mass()`, `radius_of_gyration()` and `principal_axes()`.

# ### Exercise 4

# 
# Let's look in more detail at the AdK protein. AdK has three domains:
# 
#  - CORE (residues 1-29, 60-121, 160-214)
#  - NMP (residues 30-59)
#  - LID (residues 122-159)
# 

# <center><img src="imgs/adk.png" alt="mda" style="width: 300px;"/></center>

# **4a. Calculate the center of mass of each of the three domains**

# In[60]:


domains = {
           'CORE': u.select_atoms("resid 1-29 60-121 160-214"),
           'NMP': u.select_atoms("resid 30-59"),
           'LID': u.select_atoms("resid 122-159")
          }
cogs = {name: ag.center_of_mass() for name,ag in domains.items()} 
print(cogs)


# In[61]:


# Exercise 4a


# How do these differ from the centers of geometry?

# Angles between these domains can be used to distinguish open and closed states of the protein (see the figure above) [1]. These angles are defined between the center of geometry of the backbone and C$_\beta$ atoms of the following groups of atoms:
# 
#  - $\theta_{NMP}$ is defined between residues:
#    - A: 115-125 
#    - B: 90-100
#    - C: 35-55
#  - $\theta_{LIC}$ is defined between residues:
#    - A: 179-185
#    - B: 112-125
#    - C: 125-153 
#  
# **4b. Calculate the current values of these angles**
# 
# Hints:
#  - The angle between two vectors is given by:
# $$
# \theta = arccos\left( \frac{\vec{BA}\cdot\vec{BC}}{|\vec{BA}||\vec{BC}|} \right)
# $$
# 
# - You can use numpy `numpy.linalg.norm()` to calculate the norm of a vector. Numpy also has functions `numpy.arccos()` and `numpy.dot()`.

# In[62]:


import numpy as np
from numpy.linalg import norm

A_NMP = u.select_atoms('resid 115-125 and (backbone or name CB)').center_of_geometry()
B_NMP = u.select_atoms('resid 90-100 and (backbone or name CB)').center_of_geometry()
C_NMP = u.select_atoms('resid 35-55 and (backbone or name CB)').center_of_geometry()
BA_NMP = A_NMP - B_NMP
BC_NMP = C_NMP - B_NMP
theta_NMP = np.arccos(np.dot(BA_NMP, BC_NMP)/(norm(BA_NMP)*norm(BC_NMP)))
print('theta_NMP: ', np.rad2deg(theta_NMP))


A_LID = u.select_atoms('resid 179-185 and (backbone or name CB)').center_of_geometry()
B_LID = u.select_atoms('resid 112-125 and (backbone or name CB)').center_of_geometry()
C_LID = u.select_atoms('resid 125-153 and (backbone or name CB)').center_of_geometry()
BA_LID = A_LID - B_LID
BC_LID = C_LID - B_LID
theta_LID = np.arccos(np.dot(BA_LID, BC_LID)/(norm(BA_LID)*norm(BC_LID)))
print('theta_LID: ', np.rad2deg(theta_LID))


# In[63]:


# Exercise 4b


# Currently, we're just looking at a static, single frame, so we only get one set of angle values. In the next Session you'll learn how we can get a timeseries of data from a trajectory.

# [1] O Beckstein. EJ Denning, JR Perilla, and TB Woolf. Zipping and Unzipping of Adenylate Kinase: Atomistic Insights into the Ensemble of Open/Closed Transitions. J Mol Biol 394 (2009), 160–176. doi:10.1016/j.jmb.2009.09.009

# # 5. Bonds, angles, and dihedrals
# 
# We can also get at connectivity information between atoms, such as bonds, angles, and dihedrals 

# In[64]:


u.atoms.bonds


# In[65]:


u.atoms.angles


# In[66]:


u.atoms.dihedrals


# In[67]:


u.bonds[3]


# To get the actual value:

# In[68]:


u.bonds[3].value()


# We can plot the distribution of bond lengths using matplotlib:

# In[69]:


import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

plt.hist(u.bonds.values(), bins=30)


# These work the same way as `AtomGroup`s. They're sliceable, and indexing them works too to give individual bonds, angles, dihedrals.

# # 6. Optional Tasks

# - If you haven't already, read through the [selection documentation](https://docs.mdanalysis.org/2.0.0-dev0/documentation_pages/selections.html) and play around with all the selection options not covered above
# 
# 
# - [nglview](https://github.com/nglviewer/nglview#usage) has a lot of options for visualising - you can add multiple selections to one view, change their colour and representation style and more - look through their documentaiton and see what you can create!

# In[ ]:




