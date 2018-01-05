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
from math import sqrt, ceil, floor
from numpy.ma.core import arctan2, sin, cos
from shapely.geometry.collection import GeometryCollection
from functional import compose, partial
from functools import reduce

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

def puto_a_posicion(puto):
    return puto.x, puto.y

def puto_a_posicion_polar(centro, puto):
    return posicion_a_posicion_polar(centro, puto_a_posicion(puto))

def posicion_valida_dimension_igual(posiciones, para_dimension_x):
    idx_dimension = 0
    
    assert len(posiciones)
    
    if not para_dimension_x:
        idx_dimension = 1
    
    referencia = posiciones[0][idx_dimension]
        
    return all(map(lambda posi:posi[idx_dimension] == referencia, posiciones))

def posicion_valida_colinearidad(posiciones):
    return posicion_valida_dimension_igual(posiciones, True) and posicion_valida_dimension_igual(posiciones, False)

# XXX: https://stackoverflow.com/questions/21291725/determine-if-shapely-point-is-within-a-linestring-multilinestring
def puto_valida_colinearidad(segmento, putos):
    return all(map(lambda puto:segmento.distance(puto) < 1e-8, putos))

# XXX: https://stackoverflow.com/questions/16739290/composing-functions-in-python
def caca_comun_composicion(functiones):
    return partial(reduce, compose)(functiones)

class sektor_cir_culo():

    angulos_extremos_verticales = [np.pi / 2, -np.pi / 2]
    angulos_extremos_horizontales = [0, np.pi]
    
    def __init__(self, centro, radio, pos_lim_1, pos_lim_2, lim_x, lim_y):
        self.centro = centro
        self.radio = radio
        self.lim_x = lim_x
        self.lim_y = lim_y
        self.pos_lim_1 = min(pos_lim_1, pos_lim_2)
        self.pos_lim_2 = max(pos_lim_1, pos_lim_2)
        self.angulo = None
        self.segmento_sektor_1 = None
        self.segmento_sektor_2 = None
        self.circulo = None
        self.pos_polar_seg_1 = None
        self.pos_polar_seg_2 = None
        self.normalizar_angulos = False
        self.pos_pol_final_segmento_1 = None
        self.pos_pol_final_segmento_2 = None
        self.pos_final_segmento_1 = None
        self.pos_final_segmento_2 = None
        self.valores_x_abarcados = None
        self.valores_y_abarcados = None
        self.extremo_caja_vertical_min = None
        self.extremo_caja_vertical_max = None
        self.extremo_caja_horizontal_min = None
        self.extremo_caja_horizontal_max = None
    
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
        self.pos_final_segmento_1 = pos_final_segmento_1 = sektor_cir_culo.genera_segmento_linea(self.centro, self.pos_lim_1, self.radio)
        self.pos_final_segmento_2 = pos_final_segmento_2 = sektor_cir_culo.genera_segmento_linea(self.centro, self.pos_lim_2, self.radio)
        
        self.pos_pol_final_segmento_1 = self.posicion_a_posicion_polar_de_sektor(self.pos_final_segmento_1)
        self.pos_pol_final_segmento_2 = self.posicion_a_posicion_polar_de_sektor(self.pos_final_segmento_2)
        
        self.segmento_sektor_1 = LineString([Point(self.centro), Point(pos_final_segmento_1)])
        self.segmento_sektor_2 = LineString(map(Point, [self.centro, pos_final_segmento_2]))
        
        self.circulo = Point(self.centro).buffer(self.radio)
        
        self._inicializa_posiciones_polares()
        
        logger_cagada.debug("las posiciones polares de {} {} son {} {}".format(self.segmento_sektor_1, self.segmento_sektor_2, self.pos_polar_seg_1, self.pos_polar_seg_2))
    
    def _inicializa_posiciones_polares(self):
        pos_polar_seg_1 = self.pos_pol_final_segmento_1
        pos_polar_seg_2 = self.pos_pol_final_segmento_2
        
        angulo_1 = pos_polar_seg_1[0]
        angulo_2 = pos_polar_seg_2[0]
        
        angulo_1, angulo_2 = sorted([angulo_1, angulo_2])
        
        es_ang_negativo_1 = angulo_1 < 0
        es_ang_negativo_2 = angulo_2 < 0
        
        if es_ang_negativo_1 != es_ang_negativo_2:
            dif_angulos = angulo_2 - angulo_1
            if dif_angulos > np.pi:
                angulo_1, angulo_2 = angulo_2, np.pi * 2 + angulo_1
                self.normalizar_angulos = True
        
        pos_polar_seg_1[0] = angulo_1
        pos_polar_seg_2[0] = angulo_2
        self.pos_polar_seg_1 = pos_polar_seg_1
        self.pos_polar_seg_2 = pos_polar_seg_2
    
    @property
    def posiciones_finales_segmentos(self):
        return [self.pos_final_segmento_1, self.pos_final_segmento_2]
    
    @property
    def posiciones_polares_finales_segmentos(self):
        return [self.pos_pol_final_segmento_1, self.pos_pol_final_segmento_2]
    
    def _calcula_extremos_caja(self, calcula_extremos_verticales):
        angulos_extremos = []
        pos_pol_extremas = []
        pos_pol_extremas_en_sektor = []
        limite_max = None
        limite_min = None
        idx_dimension = 0
        posiciones_potenciales_extremas = self.posiciones_finales_segmentos + [self.centro]
        
        if calcula_extremos_verticales:
            angulos_extremos = sektor_cir_culo.angulos_extremos_verticales
            idx_dimension = 0
        else:
            angulos_extremos = sektor_cir_culo.angulos_extremos_horizontales
            idx_dimension = 1
        
        pos_pol_extremas = map(lambda angulo:(angulo, self.radio), angulos_extremos)
        
        pos_pol_extremas_en_sektor = filter(self.posicion_polar_en_cuerda_sektor, pos_pol_extremas)
        
        assert len(pos_pol_extremas_en_sektor) < 2
        
        if len(pos_pol_extremas_en_sektor):
            posiciones_potenciales_extremas.append(self.posicion_a_posicion_polar_de_sektor(pos_pol_extremas_en_sektor[0]))
        
        limite_max = max(posiciones_potenciales_extremas, key=lambda posi:posi[idx_dimension])[idx_dimension]
        limite_min = min(posiciones_potenciales_extremas, key=lambda posi:posi[idx_dimension])[idx_dimension]
        
        return limite_min, limite_max

    def _calcula_caja_sektor(self):
        self.extremo_caja_vertical_min, self.extremo_caja_vertical_max = self._calcula_extremos_caja(True)
        self.extremo_caja_horizontal_min, self.extremo_caja_horizontal_max = self._calcula_extremos_caja(False)
        
        self.valores_x_abarcados = range(max(ceil(self.extremo_caja_vertical_min), 0), min(floor(self.extremo_caja_vertical_max), self.lim_x))
        self.valores_y_abarcados = range(max(ceil(self.extremo_caja_horizontal_min), 0), min(floor(self.extremo_caja_horizontal_max), self.lim_y))
    
    def normaliza_posicion_polar(self, pos_pol):
        if self.normalizar_angulos:
            angulo_pol = pos_pol[0]
            if angulo_pol < 0:
                angulo_pol = np.pi * 2 + angulo_pol
            pos_pol[0] = angulo_pol
    
    def posicion_polar_dentro_de_sektor(self, pos_pol):
        if self.pos_polar_seg_1[0] <= pos_pol[0] <= self.pos_polar_seg_2[0]:
            return True
        else:
            return False
    
    def posicion_polar_en_cuerda_sektor(self, pos_pol):
        return self.posicion_polar_dentro_de_sektor(pos_pol) and pos_pol[1] == self.radio
    
    def posicion_en_cuerda_sektor(self, posi):
        return self.posicion_polar_en_cuerda_sektor(self.posicion_a_posicion_polar_de_sektor(posi))
    
    def puto_en_cuerda_sektor(self, puto):
        return self.posicion_en_cuerda_sektor(puto_a_posicion(puto))
    
    def posicion_a_posicion_polar_de_sektor(self, posi):
        return posicion_a_posicion_polar(self.centro, posi)
    
    def posicion_polar_a_posicion_de_sektor(self, pos_pol):
        return posicion_polar_a_posicion(self.centro, pos_pol)
    
    @property
    def segmentos_sektor(self):
        return [self.segmento_sektor_1, self.segmento_sektor_2]
    
    def valida_colinearidad(self, putos):
        return posicion_valida_colinearidad(map(lambda puto:(puto.x, puto.y), self.segmentos_sektor + putos)) or any(map(lambda segmento:puto_valida_colinearidad(segmento, putos), self.segmentos_sektor))
    
    def calcula_interseccion_con_arco_sektor(self, segmento):
        putos_intersex = []
        puto_tangencial = False
        intersex = self.circulo.intersection(segmento)
        assert not isinstance(intersex, GeometryCollection)
        if not isinstance(intersex, Point):
            assert isinstance(intersex, LineString)
            
#            putos_intersex = map(lambda pos_pol:Point(self.posicion_polar_a_posicion_de_sektor(pos_pol)), filter(lambda pos_pol:self.posicion_polar_en_cuerda_sektor(pos_pol), map(lambda posi:self.normaliza_posicion_polar(self.posicion_a_posicion_polar_de_sektor(posi)), list(intersex.coords))))
            putos_intersex = map(Point, filter(caca_comun_composicion(self.posicion_a_posicion_polar_de_sektor, self.normaliza_posicion_polar, self.posicion_polar_en_cuerda_sektor) , list(intersex.coords)))
            
#            for posicion in list(intersex.coords):
#                pos_pol = posicion_a_posicion_polar(self.centro, posicion)
#                self.normaliza_posicion_polar(pos_pol)
#                esta_en_sector = self.posicion_polar_en_cuerda_sektor(pos_pol)
#                if esta_en_sector:
#                    putos_intersex.append(Point(posicion))
        else:
            putos_intersex.append(intersex)
            assert self.puto_en_cuerda_sektor(intersex)
            puto_tangencial = True
        
        assert self.valida_colinearidad(putos_intersex)
        
        return putos_intersex, puto_tangencial
    
    def calcula_interseccion_con_segmento_de_linea(self, segmento, con_segmento_sektor_1):
        if con_segmento_sektor_1:
            segmento_sektor = self.segmento_sektor_1
        else:
            segmento_sektor = self.segmento_sektor_2
        
        intersex = segmento.intersection(segmento_sektor)
        
        assert (not isinstance(intersex, LineString))
        
        if isinstance(intersex, Point):
            return intersex
        else:
            return None
    
    def calcula_intersexiones_de_segmento_con_sektor(self, segmento):
        putos_intersexion = []

        intersex_seg_sektor_1 = self.calcula_interseccion_con_segmento_de_linea(segmento, True)
        intersex_seg_sektor_2 = self.calcula_interseccion_con_segmento_de_linea(segmento, False)
        intersexs_arco, puto_tangencial = self.calcula_interseccion_con_arco_sektor(segmento)
        
        intersex_arco_cnt = len(intersexs_arco)
        if(intersex_arco_cnt == 2 or puto_tangencial):
            assert (puto_tangencial and intersex_arco_cnt == 1) or (not puto_tangencial and intersex_arco_cnt == 2)
            
            putos_intersexion = intersexs_arco
            
            if(puto_tangencial):
                putos_intersexion += intersexs_arco
        else:
            if(intersex_arco_cnt == 1):
                assert (intersex_seg_sektor_1 is None) != (intersex_seg_sektor_2 is None)
                putos_intersexion = intersexs_arco + [intersex_seg_sektor_1 if intersex_seg_sektor_1 else intersex_seg_sektor_2]
            else:
                assert intersex_seg_sektor_1 is not None and intersex_seg_sektor_2 is not None
                
                putos_intersexion = [intersex_seg_sektor_1, intersex_seg_sektor_2]
        
        assert self.valida_colinearidad(putos_intersexion)
        
        assert (puto_tangencial and len(putos_intersexion) == 1) or (not puto_tangencial and len(putos_intersexion) == 2)
        
        return map(puto_a_posicion, putos_intersexion)
    
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
