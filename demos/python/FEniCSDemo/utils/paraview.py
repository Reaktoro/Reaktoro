'''
Created on 13 Feb 2015

@author: Allan Leal
'''

import sys
import xml.etree.ElementTree as et

def xdmf(filename, **kwargs):
    ''' Convert a xdmf file from dolfin to a xdmf that is correctly parsed by Paraview'''

    output = kwargs.get('output', filename)
    suffix = kwargs.get('suffix', None)

    tree = et.parse(filename)
    root = tree.getroot()
    domain = root.find('Domain')
    grids = domain.findall('Grid')

    timeseries = grids[0]

    other_timeseries = grids[1:]

    attributes = []

    for grid in other_timeseries:
        innergrids = grid.findall('Grid')
        attributes_aux = []
        for x in innergrids:
            attribute = x.find('Attribute')
            attributes_aux.append(attribute)
        attributes.append(attributes_aux)

    timeseries_grids = timeseries.findall('Grid')

    for grid_attributes in attributes:
        for (grid, attribute) in zip(timeseries_grids, grid_attributes):
            grid.append(attribute)

    for x in other_timeseries:
        domain.remove(x)

    if suffix is not None:
        for attribute in root.iter('Attribute'):
            attribute.attrib['Name'] += suffix

    tree.write(output, encoding="utf-8", xml_declaration=True)
