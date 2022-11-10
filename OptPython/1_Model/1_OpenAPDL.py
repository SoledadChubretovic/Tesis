
#%% Open APDL license

from ansys.mapdl.core import Mapdl, launch_mapdl
mapdl = launch_mapdl()
Mapdl(start_instance=True)
print(mapdl)

# %%