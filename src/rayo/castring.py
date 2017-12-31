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
from math import sqrt
from numpy.ma.core import arctan2, sin, cos
from shapely.geometry.collection import GeometryCollection

nivel_log = logging.ERROR
nivel_log = logging.DEBUG
logger_cagada = None

immutable_types = set((int, str))

BLUE = '#6699cc'
GRAY = '#999999'

cont_figs = 1
fig = None
ax = None


# XXX: http://code.activestate.com/recipes/576527-freeze-make-any-object-immutable/
@total_ordering
class Frozen(object):

    def __init__(self, value):
        assert isinstance(value, LineString), "{} no es la instancia esperada".format(type(value))
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

# XXX: https://stackoverflow.com/questions/3942303/how-does-a-python-set-check-if-two-objects-are-equal-what-methods-does-an-o
    def __hash__(self):
        putos = None
        mierda = self._value
# XXX: https://stackoverflow.com/questions/20474549/extract-points-coordinates-from-python-shapely-polygon
        putos = mapping(mierda)["coordinates"]
#        logger_cagada.debug("de {} el hash es {}".format(mierda, putos.__hash__()))
        return putos.__hash__()

    def __lt__(self, orto):
        mierda = self._value
# XXX: https://stackoverflow.com/questions/20474549/extract-points-coordinates-from-python-shapely-polygon
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


def posicion_suma(pos_1, pos_2):
    return (pos_1[0] + pos_2[0], pos_1[1] + pos_2[1])


def posicion_resta(pos_1, pos_2):
    return posicion_suma(pos_1, (-pos_2[0], -pos_2[1]))

    
def posicion_a_posicion_polar(centro, posi):
    distancias = posicion_resta(posi, centro)
    angulo = arctan2(distancias[0], distancias[1])
    radio = sqrt(pow(distancias[0], 2) + pow(distancias[1], 2))
    return (angulo, radio)


def posicion_polar_a_posicion(centro, pos_pol):
    dist_x = sin(pos_pol[0]) * pos_pol[1]
    dist_y = cos(pos_pol[0]) * pos_pol[1]
    poli = posicion_suma(centro, (dist_x, dist_y))
    return poli


class sektor_cir_culo():

    def __init__(self, centro, radio, pos_lim_1, pos_lim_2):
        self.centro = centro
        self.radio = radio
        self.pos_lim_1 = min(pos_lim_1, pos_lim_2)
        self.pos_lim_2 = max(pos_lim_1, pos_lim_2)
        self.angulo = None
        self.segmento_sektor_1 = None
        self.segmento_sektor_2 = None
        self.circulo = None
        self.pos_polar_seg_1 = None
        self.pos_polar_seg_2 = None
    
    @classmethod
    def calcula_angulo_sektor(cls, centro, pos_lim_1, pos_lim_2):
        pos_lim_1_int = posicion_resta(pos_lim_1, centro)
        pos_lim_2_int = posicion_resta(pos_lim_2, centro)
        centro_int = posicion_resta(centro, centro)
        assert centro_int == (0, 0)
        
        angulo_1 = arctan2(pos_lim_1_int[0], pos_lim_1_int[1]) * 180 / np.pi
        angulo_2 = arctan2(pos_lim_2_int[0], pos_lim_2_int[1]) * 180 / np.pi
        
        if(angulo_1 < 0):
            angulo_1 = 360 + angulo_1
            
        if(angulo_2 < 0):
            angulo_2 = 360 + angulo_2
        
        angulo = max(angulo_1, angulo_2) - min(angulo_1, angulo_2)
        
        if angulo > 180:
            angulo = 360 - angulo
        
        return angulo
    
    @classmethod
    def genera_segmento_linea(cls, centro, pos, tam):
        polar = posicion_a_posicion_polar(centro, pos)
        polar[1] = tam
        final_segmento = posicion_polar_a_posicion(centro, polar)
        return final_segmento
    
    def _inicializa(self):
        self.angulo = sektor_cir_culo.calcula_angulo_sektor(self.centro, self.pos_lim_1, self.pos_lim_2)
        final_segmento_1 = sektor_cir_culo.genera_segmento_linea(self.centro, self.pos_lim_1, self.radio)
        final_segmento_2 = sektor_cir_culo.genera_segmento_linea(self.centro, self.pos_lim_2, self.radio)
        
        self.segmento_sektor_1 = LineString([Point(self.centro), Point(final_segmento_1)])
        self.segmento_sektor_2 = LineString(map(lambda pos:Point(pos[0], pos[1]), [self.centro, final_segmento_2]))
        
        self.circulo = Point(self.centro).buffer(self.radio)
        
        pos_polar_seg_1 = posicion_a_posicion_polar(self.centro, final_segmento_1)
        pos_polar_seg_2 = posicion_a_posicion_polar(self.centro, final_segmento_2)
        self.pos_polar_seg_1 = min(pos_polar_seg_1, pos_polar_seg_2)
        self.pos_polar_seg_2 = max(pos_polar_seg_1, pos_polar_seg_2)
    
    def calcula_interseccion_segmento_de_linea_con_semicirculo_sektor(self, segmento):
        intersex = self.circulo.intersection(segmento)
        assert not isinstance(intersex, GeometryCollection)
        if not isinstance(intersex, Point):
            assert isinstance(intersex, LineString)
            for posicion in list(intersex.coords):
                pos_pol = posicion_a_posicion_polar(self.centro, posicion)
                
    
    def calcula_interseccion_con_segmento_de_linea(self, segmento, con_segmento_sektor_1):
        if con_segmento_sektor_1:
            segmento_sektor = self.segmento_sektor_1
        else:
            segmento_sektor = self.segmento_sektor_2
        
        intersex = segmento.intersection(segmento_sektor)
        
        assert (not isinstance(intersex, LineString))
        
        if not isinstance(intersex, Point):
            pass
    
    def __repr__(self):
        return "centro {} segmentos sektor {} y {}".format(self.centro, self.segmento_sektor_1, self.segmento_sektor_2)

    
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
    if isinstance(contorno, LineString):
        contorno_externo = list(contorno.coords)
        contornos_internos = []
    else:
        if isinstance(contorno, MultiLineString):
            contorno_externo = list(contorno[0].coords)
            contornos_internos = list(map(lambda multlinea:list(multlinea.coords), contorno[1:]))
            pass
        else:
            assert False, "no es linea ni multiplea {}".format(contorno)
    logger_cagada.debug("el contorno ext es {}, los int {}".format(contorno_externo, contornos_internos))
# XXX: https://gis.stackexchange.com/questions/72306/does-shapely-within-function-identify-inner-holes
    poligono = Polygon(contorno_externo, contornos_internos)
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
    logger_cagada.debug("de celda {} se generaro poligono {} de cirds {}".format(celda_inicial, poligono, list(poligono.exterior.coords)))
    logger_cagada.debug("cekdas invol {}".format(sorted(list(celdas_ya_visitadas_int))))
    
    global cont_figs
    global fig
    global ax
    poly = mapping(poligono)
    patch = PolygonPatch(poly, fc=BLUE, ec=GRAY, alpha=0.5, zorder=2)
    ax.add_patch(patch)
    ax.set_xlim(xmin=-10, xmax=10)
    ax.set_ylim(ymin=-10, ymax=10)
    cont_figs += 1

                    
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
    
    global fig
    global ax
    fig = pyplot.figure(1, dpi=180)
    ax = fig.add_subplot(121)
    house_of_pain_genera_poligonos_y_putos(matrix, mapa_puto_a_forma, mapa_linea_a_forma, mapa_celda_a_poligono)
    pyplot.show()


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
