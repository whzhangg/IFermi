# this is supposed to be an adaptor that read an nscf.out calculation and return an band structure class


import re
import numpy as np
from ase.data import atomic_numbers
from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin

Bohr2A = 0.529177210 # A

class nscfOut:

    def __init__(self, filename: str):
        self.output = filename
        
        self.spin_polarized = False
        self.efermi = None

        # crystal structure
        self.positions = None
        self.types = None
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
        energy_up = []
        energy_down = []
        symbols = []
        k_frac = []
        efermi = 0
        
        with open(self.output,'r') as f:
            aline=f.readline()

            while aline:

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

                elif re.search('k\s+=\s*-?\d+\.\d+\s*-?\d+\.\d+\s*-?\d+\.\d+\s',aline) != None:
                    kstr=re.findall(r'-?\d+\.\d+',aline)

                    f.readline()

                    energy = []
                    while len(energy) < nbnd:
                        aline = f.readline()
                        data = np.array(aline.split(), dtype = float)
                        for d in data:
                            energy.append(d)
                    if len(energy) > nbnd:
                        raise "length of energy > nbnd"
                    #self.energy.append(energy)
                    energy_up.append(energy)
                
                aline = f.readline()

        self.efermi = efermi
        self.lattice = lattice
        self.types = [ atomic_numbers[s] for s in symbols ]
        self.symbols = symbols 
        self.positions = np.array(positions)

        self.reciprocal = recip

        self.kpoints = k_frac

        self.eig = {}
        energy_spin = [np.array(energy_up).T, np.array(energy_down).T]
        for i, s in enumerate(Spin):
            self.eig[s] = energy_spin[i]

    def __distance(self,k1,k2):
        dk = np.array(k1) - np.array(k2)
        return np.sqrt(dk.dot(dk))

from collections import defaultdict

def QE_get_band_structure(nscfout: str) -> BandStructure:
    #
    result = nscfOut(nscfout)
    lattice_new = Lattice(result.reciprocal)

    #eigenvals: DefaultDict[Spin, list] = defaultdict(list)
    structure = Structure(result.lattice, result.symbols, result.positions)

    band = BandStructure(result.kpoints, 
                         result.eig, 
                         lattice_new, 
                         result.efermi, 
                         structure = structure)
    return band

QE_get_band_structure("MgB2/wannier/nscf.out")