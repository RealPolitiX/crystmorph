#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: R. Patrick Xian
"""

import numpy as np

sin = np.sin
cos = np.cos


###############################################
# 3D transformation (homogeneous coordinates) #
###############################################

def translation3D(tx=0, ty=0, tz=0):
    """ Translation matrix in 3D.
    """

    mat = np.array([[1, 0, 0, tx],
                    [0, 1, 0, ty],
                    [0, 0, 1, tz],
                    [0, 0, 0,  1]])

    return mat


def scaling3D(sx, sy, sz):
    """ Scaling matrix in 3D.
    """

    mat = np.array([[sz, 0, 0, 0],
                    [0, sy, 0, 0],
                    [0, 0, sz, 0],
                    [0, 0, 0,  1]])

    return mat


def rotation_x(rx, radian=False):
    """ Rotation matrix about x axis.
    """

    if not radian:
        rx = np.radians(rx)
    
    mat = np.array([[1,       0,        0, 0],
                    [0, cos(rx), -sin(rx), 0],
                    [0, sin(rx),  cos(rx), 0],
                    [0,       0,        0, 1]])
    
    return mat


def rotation_y(ry, radian=False):
    """ Rotation matrix about y axis.
    """

    if not radian:
        ry = np.radians(ry)

    mat = np.array([[cos(ry), 0, -sin(ry), 0],
                    [0,       1,        0, 0],
                    [sin(ry), 0,  cos(ry), 0],
                    [0,       0,        0, 1]])

    return mat


def rotation_z(rz, radian=False):
    """ Rotation matrix about z axis.
    """

    if not radian:
        rz = np.radians(rz)

    mat = np.array([[cos(rz), -sin(rz), 0, 0],
                    [sin(rz),  cos(rz), 0, 0],
                    [0,              0, 1, 0],
                    [0,              0, 0, 1]])

    return mat


def Glazer_deform(a, b, c):
    """ Deformation matrix formulation of Glazer tilt system (single/double perovskites).
    """

    mat = np.array([[0, -c, b],
                    [-c, 0, a],
                    [b, -a, 0]])

    return mat


def cart2homo(coords_cart):
    """ Transform from Cartesian to homogeneous coordinates.
    """

    ones = np.ones((len(coords_cart), 1))
    coords_homo = np.squeeze(np.concatenate((coords_cart, ones), axis=1))

    return coords_homo


def homo2cart(coords_homo):
    """ Transform from homogeneous to Cartesian coordinates.
    """

    try:
        coords_cart = np.squeeze(coords_homo)[:,:2]
    except:
        coords_cart = np.squeeze(coords_homo)[:2]

    return coords_cart


###################################
# Crystallographic transformation #
###################################

