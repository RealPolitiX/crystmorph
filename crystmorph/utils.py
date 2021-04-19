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

    info = list(map(lambda s: attrgetter(info_list)(s), iterable))
    info = np.array(info).T
    info_dict = {name: value for name, value in zip(info_list, info)}

    return info_dict