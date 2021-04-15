#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: R. Patrick Xian
"""

import numpy as np


class ConvexPolyhedron(object):
    """ Root class for convex polyhedra.
    """
    
    def __init__(self, center, n_vertex=None, n_edge=None):
        
        self.center = center
        self.n_vertex = n_vertex
        self.n_edge = n_edge
    
    @property
    def n_face(self):
        """ Number of faces (calculated with Euler's formula).
        """
        
        try:
            n_face = self.n_vertex + self.n_edge - 2
            return n_face
        except:
            return None
    
    @property
    def is_proper_shape(self):
        """ Check if the shape is proper.
        """
        
        try:
            if (self.n_vertex > 0) and (self.n_edge > 0) and (self.n_face > 0):
                return True
        except:
            return False


class Face3D(object):
    """ Polygon data structure featuring an ordered vertex list.
    """
    
    def __init__(self, vertices, order='cw'):
        
        self.vertices = vertices
        self.n_vertex = len(vertices)
        if self.n_vertex < 3:
            raise ValueError('Number of vertices should be more than 3!')
        else:
            # Calculate basic properties of the polygon
            self.vertex_order(order=order)
            self.calculate_edges()
            self.calculate_normal()
        
    def vertex_order(self, order='cw'):
        """ Vertex list ordering ('cw' as clockwise or 'ccw' as counterclockwise).
        """
        
        keys = list(range(1, self.n_vertex+1))
        self.vertdict = dict(zip(keys, self.vertices))
        
    def calculate_normal(self):
        """ Calculate the unit normal vector of the face in 3D.
        """
        
        normal = np.cross(self.edges[1], self.edges[2])
        self.unit_normal = normal / np.linalg.norm(normal)
    
    def calculate_edges(self):
        """ Calculate the edge properties of vertices.
        """
        
        self.edges = np.diff(self.vertices, axis=0, append=self.vertices[:1])
        self.edgelengths = np.linalg.norm(self.edges, axis=1)
        
    def calculate_area(self):
        """ Area of the face in 3D (polygon area).
        """
        
        total = np.array([0, 0, 0])
        for i in range(self.n_vertex):
            vi1 = self.vertices[i]
            if i == self.n_vertex-1:
                vi2 = self.vertices[0]
            else:
                vi2 = self.vertices[i+1]
            prod = np.cross(vi1, vi2)
            total += prod
        
        self.area = abs(np.dot(total, self.unit_normal) / 2)