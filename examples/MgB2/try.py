# this script is taken from the github page of ifermi
# to run as test for the MgB2 result as shown here

from pymatgen.io.vasp.outputs import Vasprun
from ifermi.surface import FermiSurface
from ifermi.interpolate import FourierInterpolator
from ifermi.plot import FermiSlicePlotter, FermiSurfacePlotter, save_plot, show_plot
from ifermi.kpoints import kpoints_from_bandstructure

# load VASP calculation outputs
vr = Vasprun("vasprun.xml")
bs = vr.get_band_structure()

# interpolate the energies onto a dense k-point mesh
interpolator = FourierInterpolator(bs)
dense_bs, velocities = interpolator.interpolate_bands(return_velocities=True)

# generate the Fermi surface and calculate the dimensionality
fs = FermiSurface.from_band_structure(
  dense_bs, mu=0.0, wigner_seitz=True, calculate_dimensionality=True
)

# generate the Fermi surface and calculate the group velocity at the
# center of each triangular face
dense_kpoints = kpoints_from_bandstructure(dense_bs)
fs = FermiSurface.from_band_structure(
  dense_bs, mu=0.0, wigner_seitz=True, calculate_dimensionality=True,
  property_data=velocities, property_kpoints=dense_kpoints
)

# number of isosurfaces in the Fermi surface
#fs.n_surfaces

# number of isosurfaces for each Spin channel
#fs.n_surfaces_per_spin

# the total area of the Fermi surface
#fs.area

# the area of each isosurface
#fs.area_surfaces

# loop over all isosurfaces and check their properties
# the isosurfaces are given as a list for each spin channel
#for spin, isosurfaces in fs.isosurfaces.items():
#    for isosurface in isosurfaces:
        # the dimensionality (does the surface cross periodic boundaries)
#        isosurface.dimensionality

        # what is the orientation
#        isosurface.orientation

        # does the surface have face properties
#        isosurface.has_properties

        # calculate the norms of the properties
#        isosurface.properties_norms

        # calculate scalar projection of properties on to [0 0 1] vector
#        isosurface.scalar_projection((0, 0, 1))

        # uniformly sample the surface faces to a consistent density
#        isosurface.sample_uniform(0.1)

# plot the Fermi surface
fs_plotter = FermiSurfacePlotter(fs)
plot = fs_plotter.get_plot()

save_plot(plot, "fermi-surface.png")  # saves the plot to a file
show_plot(plot)  # displays an interactive plot

# generate Fermi slice along the (0 0 1) plane going through the Γ-point.
fermi_slice = fs.get_fermi_slice((0, 0, 1))

# number of isolines in the slice
fermi_slice.n_lines

# do the lines have segment properties
fermi_slice.has_properties

# plot slice
slice_plotter = FermiSlicePlotter(fermi_slice)
plot = slice_plotter.get_plot()

save_plot(plot, "fermi-slice.png")  # saves the plot to a file
show_plot(plot)  # displays an interactive plot