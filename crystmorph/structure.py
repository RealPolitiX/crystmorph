#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: R. Patrick Xian
"""

from . import polyhedron as ph, utils as u
import numpy as np
# Geometry packages
import vg

import json
from copy import deepcopy
from funcy import project
try:
    import pymatgen
    from pymatgen.core.structure import Structure
    from pymatgen.core.sites import PeriodicSite
    from pymatgen.symmetry import analyzer as sym
except:
    pass

try:
    from CifFile import ReadCif
except:
    pass

import warnings as w
w.filterwarnings("ignore")


class CIFParser(object):
    """ Parser class for crystallographic information files (CIFs).
    """

    def __init__(self, filepath, backend='pymatgen'):

        if backend == 'pymatgen':
            self.struct = Structure.from_file(filepath)
        elif backend == 'pycifrw':
            self.struct = ReadCif(filepath)


class JSONParser(object):
    """ Parser class for crystallographic information in json format.
    """

    def __init__(self, filepath, **kwargs):
        
        with open(filepath) as f:
            self.struct = json.load(f, **kwargs)


class CIFExporter(object):
    """ Exporter class for crystallographic information files (CIFs).
    """

    def __init__(self, obj, filepath, backend='pymatgen'):

        if backend == 'pymatgen':
            obj.to(filepath)
        elif backend == 'pycifrw':
            with open(filepath) as f:
                f.write(obj.WriteOut())


class JSONExporter(object):
    """ Exporter class for crystallographic information in json format.
    """

    def __init__(self, obj, filepath, **kwargs):

        with open(filepath, 'w') as f:
            json.dump(obj, f, **kwargs)


class StructureParser(object):
    """ Parser class for crystal structures. The class funnels the atomistic information
    of a crystal structure into a searchable format using a combination of new methods and
    those present in ``pymatgen.core.structure.Structure`` class.
    """

    def __init__(self, filepath, form='cif', parser_kwargs={}, **kwargs):

        if form == 'cif':
            parser = CIFParser(filepath=filepath, **parser_kwargs)
        elif form == 'json':
            parser = JSONParser(filepath=filepath, **parser_kwargs)
        self.struct = parser.struct

        # Retrieve essential chemical information (default is the atom name, number, Cartesian and fractional coordinates)
        self.atoms_info_list = kwargs.pop('atoms_info_list', ['specie.name', 'specie.number', 'coords', 'frac_coords'])
        self.atoms_dict = u.multiretrieve(self.atoms_info_list, self.struct.sites)
        indices = np.array(list(range(self.n_atoms)))
        self.atoms_dict['atom_id'] = indices
        
        self.unique_elements = [s.name for s in self.struct.types_of_specie]
        self.composition = self.struct.composition.as_dict()

    @property
    def n_atoms(self):
        """ Number of atoms in the structure.
        """

        try:
            return self.struct.num_sites
        except:
            return 0

    def filter_elements(self, keep=[], remove=[], update_structure=False):
        """ Remove atoms according to the element type.
        """

        if len(keep) > 0:
            pass
        elif len(remove) > 0:
            pass

        if update_structure:
            self.struct

        return 

    def symmetry_analyze(self, **kwargs):
        """ Produce symmetry analysis.
        """

        self.symmetry = sym.SpacegroupAnalyzer(self.struct, **kwargs)
        self.crystal_system = self.symmetry.get_crystal_system()

    def find_atoms(self, names=None):
        """ Locate the atomic species according to specification.
        """

        atom_names = self.atoms_dict['specie.name']
        if type(names) == str:
            names = [names]
        for name in names:
            loc = np.where(atom_names == name)
            selected_atoms_dict = {ai_name: ai_value[loc] for ai_name, ai_value in self.atoms_dict.items()}

        return selected_atoms_dict

    def find_neighbors(self, radius, idx=None, include_center_atom=False):
        """ Locate the neighboring atom species according to specification.
        """

        # Retrieve the atomic site instance
        site = self.struct.sites[idx]
        neighbors_sites = self.struct.get_neighbors(site, radius)
        neighbors_dict = u.multiretrieve(self.atoms_info_list, neighbors_sites)
        if include_center_atom:
            return site.as_dict(), neighbors_dict
        else:
            return neighbors_dict

    def find_local_environment(self, r_cutoff, a_cutoff=0.3, excluded_atoms=[]):
        """ Determination of local coordination environment.
        """

        from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
        from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
        from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments

        lgf = LocalGeometryFinder()
        lgf.setup_parameters(centering_type='centroid', include_central_site_in_centroid=False)
        lgf.setup_structure(structure=self.struct)
        self.struct_env = lgf.compute_structure_environments(excluded_atoms=excluded_atoms, maximum_distance_factor=1.41)
        strategy = SimplestChemenvStrategy(distance_cutoff=r_cutoff, angle_cutoff=a_cutoff)
        self.lse = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=self.struct_env)

    def perturb(self, distance, min_distance=None, keep=False, ret=True):
        """ Apply random perturbation to the atomic coordinates.
        """

        struct = deepcopy(self.struct)
        struct_perturbed = self.struct.perturb(distance, min_distance)
        if not keep:
            self.struct = struct

        if ret:
            return struct_perturbed

    def calculate_distances(atom1=[], atom2=[], atom_pairs=None):
        """ Calculate distances between specified coordinates of atom pairs.
        """

        if atom_pairs is not None:
            atom1, atom2 = atom_pairs

        d = vg.euclidean_distance(atom1, atom2)

        return d        

    def calculate_angles(atom1=[], atom2=[], atom3=[], atom_trios=None):
        """ Calculate angles between specified coordinates of atom trios.
        """

        if atom_trios is not None:
            atom1, atom2, atom3 = atom_trios
        
        ag = vg.angle(atom2-atom1, atom2-atom3)

        return ag

class PerovskiteParser(StructureParser):
    """ Parser class for perovskite structures. The single perovskite chemical formula
    follows ABX3 (A, B are cations, X is the anion). The double perovskite chemical formula
    follows AA'B2X6, A2BB'X6, or AA'BB'X6 (A, A', B, B' are cations).
    """

    def __init__(self, filepath, form='cif', category='single', parser_kwargs={}):

        super().__init__(filepath=filepath, form=form, parser_kwargs=parser_kwargs)
        self.category = category
        self.A_cation = None
        self.B_cation = None

    def count(self, propname):
        """ Count the number of things.
        """

        propval = getattr(self, propname)
        if propval is None:
            return 0
        else:
            return len(propval)

    @property
    def nA(self):
        """ Number of A cations.
        """

        return self.count('A_cation')

    @property
    def nB(self):
        """ Number of B cations.
        """

        return self.count('B_cation')

    def find_cation(self, label, names=None):
        """ Find the cation and their basic information from the known structure.
        """

        if label == 'A':
            self.A_cation = self.find_atoms(names=names)
        elif label == 'B':
            self.B_cation = self.find_atoms(names=names)

    def find_octahedron(self, radius=3.5, B_idx=0):
        """ Determine the octahedral coordination of a B cation.

        **Parameters**:
        radius: numeric
            Radial cutoff for finding the atoms within the octahedra.
        """

        if self.B_cation is None:
            raise ValueError('B-site cation information is missing.')
        else:
            octahedron_dict = self.find_neighbors(radius=radius, idx=B_idx)
            octahedron = ph.Octahedron(n_vertex=6, n_edge=12)

            return octahedron

    def get_octahedral_components(self, rcutoff=3.5):
        """ Retrieve all symmetry inequivalent octahedral components in the crystal structure.
        """

        pass

    def get_organic_components(self, method='filter_elements'):
        """ Retrive all symmetry inequivalent organic components in the crystal structure.
        The current algorithm includes filtering out the elemenets or filtering specific structural
        components (e.g. octahedra).
        """

        if self.category == 'hybrid':
            self.struct.remove_species

        pass

    def adjust_A_cation(self, idx=None, shifts=None):
        """ Adjust individual or all A-site cations (either type of atom or the position).
        """

        site_subsitute = PeriodicSite
        if idx is None: # Replace all A cations
            self.struct.specie[0] = []
        elif type(idx) in (list, tuple):
            pass
            
        else:
            raise ValueError('Index should either be None or list or tuple.')

    def adjust_B_cation(self, idx=None, shifts=None):
        """ Adjust individual or all B-site cations (either type of atom or the position).
        """

        pass

    def adjust_X_anion(self, idx=None, shifts=None):
        """ Adjust individual or all X anions.
        """

        pass

    def estimate_tilting(self):
        """ Perform tilting angle estimation.
        """

        pass

    def adjust_octahedral_tilting(self):
        """ Adjust the octahedral tilting angles.
        """

        pass

    def save_structure(self, struct, filepath, form='json', **kwargs):
        """ Save the adjusted crystal structure.
        """

        if form == 'json':
            if (type(struct) == Structure) or (type(struct) == PeriodicSite):
                dic = struct.as_dict()
            JSONExporter(obj=dic, filepath=filepath)
        elif form == 'cif':
            CIFExporter(obj=struct, filepath=filepath, **kwargs)

