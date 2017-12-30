'''
Created on 18/12/2017

@author: 

XXX: https://code.google.com/codejam/contest/1460488/dashboard#s=p3&a=3
'''
# XXX: https://stackoverflow.com/questions/30844482/what-is-most-efficient-way-to-find-the-intersection-of-a-line-and-a-circle-in-py

import logging
from sys import stdin
from matplotlib import pyplot
from descartes import PolygonPatch
import math
import json
from shapely.geometry import Point, Polygon, MultiLineString, mapping
import math
import geojson
from shapely.geometry.linestring import LineString
import numpy as np
from shapely.ops import linemerge
from collections import defaultdict
from functools import total_ordering

nivel_log = logging.ERROR
nivel_log = logging.DEBUG
logger_cagada = None

immutable_types = set((int, str))


# XXX: http://code.activestate.com/recipes/576527-freeze-make-any-object-immutable/
@total_ordering
class Frozen(object):

    def __init__(self, value):
        self._value = value

    def __getattribute__(self, name):
        if name == '_value': return super(Frozen, self).__getattribute__(name)
        v = getattr(self._value, name)
#        return v if v.__class__ in immutable_types else freeze(v)
        return v
  
    def __setattr__(self, name, value):
        if name == '_value': super(Frozen, self).__setattr__(name, value)
        else: raise Exception("Can't modify frozen object {0}".format(self._value))
    
    def __repr__(self):
        return self._value.__str__()
    
    def __hash__(self):
        putos = None
        mierda = self._value
# XXX: https://stackoverflow.com/questions/20474549/extract-points-coordinates-from-python-shapely-polygon
        putos = mapping(mierda)["coordinates"]
#        logger_cagada.debug("de {} el hash es {}".format(mierda, putos.__hash__()))
        return putos.__hash__()

    def __lt__(self, orto):
        mierda = self._value
        putos = mapping(mierda)["coordinates"]
        ass = orto._value
        fuck = mapping(ass)["coordinates"]
        return putos < fuck
    
    def __eq__(self, orto):
        mierda = self._value
        putos = mapping(mierda)["coordinates"]
        ass = orto._value
        fuck = mapping(ass)["coordinates"]
        return putos == fuck

    
def freeze(value):
  return Frozen(value)


house_of_pain_movimientos_esquinas = [(0, 0), (0, 1), (1, 1), (1, 0)]
house_of_pain_movimientos_matrix = [(-1, 0), (0, 1), (1, 0), (0, -1)]


def house_of_pain_genera_lineas_de_celda(celda):
    putos = list(map(lambda x:Point(celda[0] + x[0], celda[1] + x[1]), house_of_pain_movimientos_esquinas))
    lineas = []
    puto_ant = putos[0]
    for puto in putos[1:] + [putos[0]]:
        linea_nueva = freeze(LineString(sorted([puto_ant, puto], key=lambda puto:(puto.x, puto.y))))
#        logger_cagada.debug("linea creada {} de {} y {}".format(linea_nueva, puto_ant, puto))
        lineas.append(linea_nueva)
        puto_ant = puto
    return lineas


def house_of_pain_descongela_lista(cacas):
    return list(map(lambda x:x._value, cacas))


# XXX: https://gis.stackexchange.com/questions/223447/weld-individual-line-segments-into-one-linestring-using-shapely
def house_of_pain_crea_poligono_de_lineas(lineas):
    multicaca = MultiLineString(house_of_pain_descongela_lista(lineas))
    contorno = linemerge(multicaca)
    logger_cagada.debug("el contorno es {}".format(contorno))
    poligono = Polygon(contorno)
    assert isinstance(poligono, Polygon)
    return poligono

    
def house_of_pain_genera_poligono_y_putos_de_celda_dfs(matrix, celda_inicial, mapa_puto_a_forma, mapa_linea_a_forma, mapa_celda_a_poligono, celdas_ya_visitadas):
    pila = [celda_inicial]
    duenio = celda_inicial
    lineas_poligono = set()
    celdas_ya_visitadas_int = set()
    celdas_ya_visitadas_int.add(celda_inicial)
    while pila:
        celda_actual = pila.pop()
        logger_cagada.debug("celda act {}".format(celda_actual))
        lineas = house_of_pain_genera_lineas_de_celda(celda_actual)
        assert len(lineas) == 4
#        logger_cagada.debug("las celdas ia visitadas {}".format(celdas_ya_visitadas_int))
        assert lineas[0] != lineas[1]
        assert lineas[0] < lineas[1]
        for linea in lineas:
            logger_cagada.debug("linea act {}".format(linea))
# XXX: https://stackoverflow.com/questions/20474549/extract-points-coordinates-from-python-shapely-polygon
            putos = list(linea.coords)
            assert len(putos) == 2
            
            if linea in lineas_poligono:
                lineas_poligono.remove(linea)
                logger_cagada.debug("linea removida")
            else:
                lineas_poligono.add(linea)
                logger_cagada.debug("linea anadida")
                
        logger_cagada.debug("las lineas aora son {}".format(sorted(list(lineas_poligono))))
        
        for celda_aledana in map(lambda mov: (celda_actual[0] + mov[0], celda_actual[1] + mov[1]), house_of_pain_movimientos_matrix):
 #           logger_cagada.debug("proc celda ale {}".format(celda_aledana))
            if celda_aledana not in celdas_ya_visitadas_int and matrix[celda_aledana[0]][celda_aledana[1]] == '#':
 #               logger_cagada.debug("anadida")
                pila.append(celda_aledana)
                celdas_ya_visitadas_int.add(celda_aledana)
    
    logger_cagada.debug("las lineas del pol {}".format(lineas_poligono))
    for linea in list(lineas_poligono):
        mapa_linea_a_forma[linea] = duenio
        for puto in list(linea.coords):
            mapa_puto_a_forma[puto].append(duenio)
    
    poligono = house_of_pain_crea_poligono_de_lineas(lineas_poligono)
    
    for celda_visitada in list(celdas_ya_visitadas_int):
        mapa_celda_a_poligono[celda_visitada] = poligono
        celdas_ya_visitadas.add(celda_visitada)
    
# XXX: https://gis.stackexchange.com/questions/216745/get-polygon-shapefile-in-python-shapely-by-clipping-linearring-with-linestring
    logger_cagada.debug("de celda {} se generaro poligono {}".format(celda_inicial, list(poligono.exterior.coords)))

                    
def house_of_pain_genera_poligonos_y_putos(matrix, mapa_puto_a_forma, mapa_linea_a_forma, mapa_celda_a_poligono):
    celdas_ya_visitadas = set()
    for i, fila in enumerate(matrix[1:-1], 1):
        logger_cagada.debug("fila es {}".format(fila))
        for j, cerda in enumerate(fila[1:-1], 1):
            celda_act = (i, j)
            if(cerda == "#" and celda_act not in celdas_ya_visitadas):
                house_of_pain_genera_poligono_y_putos_de_celda_dfs(matrix, celda_act, mapa_puto_a_forma, mapa_linea_a_forma, mapa_celda_a_poligono, celdas_ya_visitadas)

                
def house_of_pain_core(matrix):
    mapa_puto_a_forma = defaultdict(lambda:[])
    mapa_linea_a_forma = {}
    mapa_celda_a_poligono = {}
    
    house_of_pain_genera_poligonos_y_putos(matrix, mapa_puto_a_forma, mapa_linea_a_forma, mapa_celda_a_poligono)


def caca_comun_lee_linea_como_num():
    return int(stdin.readline().strip())


def caca_comun_lee_linea_como_monton_de_numeros():
    return list(map(int, stdin.readline().strip().split(" ")))


def caca_comun_imprime_matrix(matrix):
    return "\n".join(matrix)


def house_of_pain_main():
    num_cacasos = caca_comun_lee_linea_como_num()
    for num_cacaso in range(num_cacasos):
        alto, ancho, limite_rayo_vallecano = caca_comun_lee_linea_como_monton_de_numeros()
        matrix = []
        matrix.append("%" * (ancho + 2))
        for _ in range(alto):
            matrix.append("%" + stdin.readline().strip() + "%")
        matrix.append("%" * (ancho + 2))
        
        logger_cagada.debug("la matrix leida\n{}".format(caca_comun_imprime_matrix(matrix)))
        house_of_pain_core(matrix)


if __name__ == '__main__':
        FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
        logging.basicConfig(level=nivel_log, format=FORMAT)
        logger_cagada = logging.getLogger("asa")
        logger_cagada.setLevel(nivel_log)
        house_of_pain_main()
