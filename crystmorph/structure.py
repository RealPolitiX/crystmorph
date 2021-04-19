#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: R. Patrick Xian
"""

from . import polyhedron as ph, utils as u
import numpy as np
import json
from copy import deepcopy
from funcy import project
try:
    import pymatgen
    from pymatgen.core.structure import Structure
    from pymatgen.core.sites import PeriodicSite
except:
    pass

try:
    from CifFile import ReadCif
except:
    pass


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
        
        self.unique_elements = [s.name for s in self.struct.types_of_specie]
        self.composition = self.struct.composition.as_dict()

    def find_atoms(self, names=None):
        """ Locate the atomic species according to specification.
        """

        atom_names = self.atoms_dict['specie.name']
        if type(names) == str:
            names = [names]
        for name in names:
            loc = np.where(atom_names == name)
        selected_atoms_dict = {ai_name: ai_value[loc] for ai_name, ai_value in self.atoms_dict}

        return selected_atoms_dict

    def find_neighbors(self, radius, name=None, idx=None, include_center_atom=False):
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

    def perturb(self, distance, min_distance=None, keep=False, ret=True):
        """ Apply random perturbation to the atomic coordinates.
        """

        struct = deepcopy(self.struct)
        struct_perturbed = self.struct.perturb(distance, min_distance)
        if not keep:
            self.struct = struct

        if ret:
            return struct_perturbed



        
