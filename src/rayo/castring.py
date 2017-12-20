'''
Created on 18/12/2017

@author: 

XXX: https://code.google.com/codejam/contest/1460488/dashboard#s=p3&a=3
'''
#XXX: https://stackoverflow.com/questions/30844482/what-is-most-efficient-way-to-find-the-intersection-of-a-line-and-a-circle-in-py

import logging
import sys
from matplotlib import pyplot
from descartes import PolygonPatch
import math
import json
from shapely.geometry import Point, Polygon, MultiLineString
import math
import geojson
from shapely.geometry.linestring import LineString
import numpy as np
from shapely.ops import linemerge

nivel_log = logging.ERROR
nivel_log = logging.DEBUG
logger_cagada = None


house_of_pain_movimientos_esquinas = [(0, 0), (0, 1), (1, 1), (1, 0)]
house_of_pain_movimientos_matrix = [(-1, 0), (0, 1), (1, 0), (0, -1)]


def house_of_pain_genera_lineas_de_celda(celda):
    putos = map(lambda x:Point(celda[0] + x[0], celda[1], x[1]), house_of_pain_movimientos_esquinas)
    lineas = []
    puto_ant = putos[0]
    for puto in putos[1:] + [putos[0]]:
        lineas.append(LineString(sorted([puto_ant, puto])))
        puto_ant = puto
    return lineas

# XXX: https://gis.stackexchange.com/questions/223447/weld-individual-line-segments-into-one-linestring-using-shapely
def house_of_pain_crea_poligono_de_lineas(lineas):
    multicaca = MultiLineString(lineas)
    contorno = linemerge(multicaca)
    poligono = Polygon(contorno)
    assert isinstance(poligono, Polygon)
    return poligono
    
def house_of_pain_genera_poligono_y_putos_de_celda_dfs(matrix, celda_inicial, mapa_puto_a_forma, mapa_linea_a_forma, mapa_celda_a_poligono, celdas_ya_visitadas):
    pila = [celda_inicial]
    duenio = celda_inicial
    lineas_poligono = set()
    while pila:
        celda_actual = pila.pop()
        lineas = house_of_pain_genera_lineas_de_celda(celda_actual)
        assert len(lineas) == 4
        celdas_ya_visitadas.add(celda_inicial)
        for linea in lineas:
# XXX: https://stackoverflow.com/questions/20474549/extract-points-coordinates-from-python-shapely-polygon
            putos = np.array(linea)
            assert len(putos) == 2
            
            if linea in mapa_linea_a_forma:
                mapa_linea_a_forma.pop(linea)
                lineas_poligono.remove(linea)
            else:
                mapa_linea_a_forma[linea] = duenio
                lineas_poligono.add(linea)
    
    for linea in list(lineas_poligono):
        for puto in list(linea.coords):
            mapa_puto_a_forma[puto].append(duenio)
    
    poligono = house_of_pain_crea_poligono_de_lineas(lineas_poligono)
    
# XXX: https://gis.stackexchange.com/questions/216745/get-polygon-shapefile-in-python-shapely-by-clipping-linearring-with-linestring
    logger_cagada.debug("de celda {} se generaro poligono {}".format(celda_inicial, list(poligono.exterior.coords)))
                    
def house_of_pain_genera_poligonos_y_putos(matrix, mapa_puto_a_forma, mapa_linea_a_forma):
    celdas_ya_visitadas = set()
    for i, fila in enumerate(matrix[1:-1], 1):
        for j, cerda in enumerate(fila[1, -1], 1):
            if(cerda == "#"):
                putos = [Point(i, j)]
                for mov in house_of_pain_movimientos_esquinas:
                    putos.append(Point(i + mov[0], j + mov[1]))
                


if __name__ == '__main__':
        FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
        logging.basicConfig(level=nivel_log, format=FORMAT)
        logger_cagada = logging.getLogger("asa")
        logger_cagada.setLevel(nivel_log)
