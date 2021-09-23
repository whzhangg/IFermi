import re
import numpy as np
from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin

# Wenhao Zhang, 2021/09/23 wenhao997@outlook.com
# examples/MgB2_QE contains test output for spin non-polarized and spin polarized calculation, which reproduce the 
# result as in MgB2

Bohr2A = 0.529177210 # A

class nscfOut:
    """
    a class written to praser nscf output of Quantum espresso pw.x

    requirement:
        verbosity = 'high'   (output energies for each k point)

    parsing assumes a fixed output format of PWSCF v. 6.4

    """
    def __init__(self, filename: str):
        self.output = filename
        
        self.spin_polarized = False
        self.efermi = None

        # crystal structure
        self.positions = None
        self.symbols = None
        self.lattice = None

        # reciprocal lattice
        self.reciprocal = None

        # kpoints
        self.kpoints = None
        self.eig = None

        self._parse_nscf()

    def _parse_nscf(self) -> None:
        """
        extract crystal structure, reciprocal cell,
        kpoints fraction coordinate and eigenvalue
        """
        alat = 0
        lattice = np.zeros((3,3))
        recip = np.zeros((3,3))
        nbnd = 0
        natom = 0
        positions = []
        nk = 0
        symbols = []
        k_frac = []
        efermi = 0

        energy = {"spinup" : [],
                  "spindown" : []
                  }

        which = "spinup"  # remember if we are reading spin up or spin down
        
        with open(self.output,'r') as f:
            aline=f.readline()

            while aline:
                # read information by checking the flags
                if "lattice parameter (alat)  =" in aline:
                    data = aline.split('=')[1]
                    data = data.split()
                    alat = float(data[0])  # in Bohr

                if "number of Kohn-Sham states" in aline:
                    data = aline.split()[-1]
                    nbnd = int(data)

                if "number of atoms/cell" in aline:
                    data = aline.split()[-1]
                    natom = int(data)

                if "crystal axes: (cart. coord. in units of alat)" in aline:
                    for i in range(3):
                        data = f.readline().split()[3:6]
                        lattice[i] = np.array(data, dtype = float) 
                    lattice *= alat * Bohr2A

                if "reciprocal axes: (cart. coord. in units 2 pi/alat)" in aline:
                    for i in range(3):
                        data = f.readline().split()[3:6]
                        recip[i] = np.array(data, dtype = float)
                    recip *= 2 * np.pi / (alat * Bohr2A)

                if "site n.     atom                  positions (cryst. coord.)" in aline:
                    for i in range(natom):
                        data = f.readline()
                        symbols.append(re.findall(r'[A-Z][a-z]*', data)[0])
                        positions.append(np.array(re.findall('-?\d+\.\d+', data), dtype = float))
                
                if "number of k points= " in aline:
                    nk = int( re.findall(r'\d+', aline)[0] )
                    k_frac = np.zeros((nk,3))

                if re.search(r'k\(.+\)\s+=\s+\(.+\)', aline) != None:
                    parts = aline.split('=')
                    ik = int( re.findall(r'\d+', parts[0])[0] )
                    pos = np.array(re.findall(r'-?\d+\.\d+', parts[1]), dtype = float)
                    k_frac[ik-1] = pos

                if "the Fermi energy is" in aline:
                    efermi = float(re.findall(r'-?\d+\.\d+', aline)[0])

                if "------ SPIN UP ------------" in aline:
                    which = "spinup"

                if "------ SPIN DOWN ----------" in aline:
                    which = "spindown"

                if re.search('k\s+=\s*-?\d+\.\d+\s*-?\d+\.\d+\s*-?\d+\.\d+\s',aline) != None:
                    kstr=re.findall(r'-?\d+\.\d+',aline)

                    f.readline()

                    lenergy = [] # local energy for each k point
                    while len(lenergy) < nbnd:
                        aline = f.readline()
                        data = np.array(aline.split(), dtype = float)
                        for d in data:
                            lenergy.append(d)

                    if len(lenergy) > nbnd:
                        raise "length of energy > nbnd"

                    energy[which].append(lenergy)
                
                aline = f.readline()

        self.efermi = efermi
        self.lattice = lattice
        self.symbols = symbols 
        self.positions = np.array(positions)
        self.reciprocal = recip
        self.kpoints = k_frac

        self.eig = {}
        self.eig[Spin.up] = np.array(energy["spinup"]).T

        if energy["spindown"]:
            self.spin_polarized = True
            self.eig[Spin.down] = np.array(energy["spindown"]).T


def QE_get_band_structure(nscfout: str) -> BandStructure:
    """
    get pymatgen.electronic_structure.bandstructure object from a nscf output
    
    args:
        file name of the nscf calculation
    """
    result = nscfOut(nscfout)

    lattice_new = Lattice(result.reciprocal)
    structure = Structure(result.lattice, result.symbols, result.positions)

    return BandStructure(result.kpoints, 
                         result.eig, 
                         lattice_new, 
                         result.efermi, 
                         structure = structure)

if __name__ == "__main__":
    bs = QE_get_band_structure("../examples/MgB2_QE/MgB2/nscf.nspin1.out")
    print(bs.bands)