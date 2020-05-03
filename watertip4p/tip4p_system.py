import numpy as np
from copy import deepcopy
from pele.angleaxis import RBTopology, RBSystem
 # For rigid fragment
import numpy as np
from pele.potentials import LJ
from pele.utils import xyz
from pele.angleaxis import rigidbody
from pele.optimize import lbfgs_py
 
# read in coordinates from xyz file
ref = xyz.read_xyz(open("TIP4P-16.xyz"))
xyz.write_xyz(open("test.xyz", "w"), coords=ref.coords)
# lookup table for atom masses
mass_lookup = {'O': 16., 'H': 1.}


def water(xyzfile): # TIP4P-16.xyz"
    ref = xyz.read_xyz(open(xyzfile))
    xyz.write_xyz(open("test.xyz", "w"), coords=ref.coords)
# lookup table for atom masses
    mass_lookup = {'O': 16., 'H': 1.}
    rb_sites = []
    for atomtype, x, i in zip(ref.atomtypes, ref.coords, xrange(len(ref.atomtypes))):
        # every 3rd atom, define a new reigid molecule
        if i % 3 == 0:
            rb = rigidbody.RigidFragment()
            rb_sites.append(rb)
        rb.add_atom(atomtype, x, mass_lookup[atomtype])
# finalize the rigid body setup
    for rb in rb_sites:
        rb.finalize_setup()
    return rb_sites, ref 

class TIP4PSystem(RBSystem):
    def __init__(self, xyzfile):
        self.xyzfile = xyzfile
        RBSystem.__init__(self)
    def setup_aatopology(self):
        water_sites, ref = water(self.xyzfile)
        rbsystem = RBTopology()
        rbsystem.add_sites(water_sites)
        #print len(rbsystem.sites), len(rbsystem.indices)
        print "I have %d water molecules in the system" % len(rbsystem.sites)
        rbcoords = rbsystem.coords_adapter(np.zeros(len(rbsystem.sites)*6))
        for site, com in zip(rbsystem.sites, rbcoords.posRigid):
            com[:] = ref.coords[site.atom_indices[0]] - site.atom_positions[0]
        pot = LJ(eps=0.1550, sig=3.1536)
        # get the flattened coordinate array
        print "The initial energy is", pot.getEnergy(ref.coords.flatten())
        rbpot = rigidbody.RBPotentialWrapper(rbsystem, pot) 
        print "rbpot.getEnergy(rbcoords.coords)", rbpot.getEnergy(rbcoords.coords)
        e, g = rbpot.getEnergyGradient(rbcoords.coords)
        g_n = rbpot.NumericalDerivative(rbcoords.coords, eps=1e-4)
        cg = rbsystem.coords_adapter(g-g_n)

        # coords = rbpot.getCoords()
        nrigid = rbcoords.size / 6
        print "nrigid", nrigid
        self.potential = rbpot
        self.nrigid = nrigid
        self.render_scale = 0.3
        self.atom_types = rbsystem.get_atomtypes()

        self.draw_bonds = []
        for i in xrange(nrigid):
            self.draw_bonds.append((3*i, 3*i+1))
            self.draw_bonds.append((3*i, 3*i+2))
        return system

    def get_potential(self):
        return self.potential
