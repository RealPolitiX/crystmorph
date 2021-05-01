#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: R. Patrick Xian
"""

import numpy as np
from operator import attrgetter


def multiretrieve(info_list, iterable):
    """ Retrieve multiple pieces of information at different levels of iterable object.
    """

    info = list(map(lambda s: attrgetter(*info_list)(s), iterable))
    info = np.asarray(info).T
    info_dict = {name: value for name, value in zip(info_list, info)}

    return info_dict


def multidict_merge(composites):
    """ Merge multiple dictionaries with the same keys.
    """
    
    if len(composites) in [0, 1]:
        return composites
    else:
        # Initialize an empty dictionary with only keys
        merged = dict.fromkeys(list(composites[0].keys()), [])
        for k in merged.keys():
            merged[k] = list(d[k] for d in composites)
    
    return merged


def csm(coords, model):
    """ Continuous symmetry measure for ordered polyhedral vertices.
    """
    
    numerator = np.sum(np.linalg.norm(coords - model, axis=1))
    center = coords.mean(axis=0)
    denominator = np.sum(np.linalg.norm(coords - center, axis=1))
    
    metric = numerator / denominator
    
    return metric