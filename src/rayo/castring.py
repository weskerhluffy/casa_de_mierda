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
from collections import defaultdict, namedtuple
from functools import total_ordering
from math import sqrt, ceil, floor
from numpy.ma.core import arctan2, sin, cos
from shapely.geometry.collection import GeometryCollection
from functional import compose, partial
from functools import reduce
from bisect import bisect_left, bisect_right
from operator import add, sub
from asyncio.log import logger
from numpy import degrees
from shapely.affinity import rotate
from itertools import chain

nivel_log = logging.ERROR
nivel_log = logging.DEBUG
logger_cagada = None

immutable_types = set((int, str))

BLUE = '#6699cc'
GRAY = '#999999'
ROJO = "#FC2C00"
VERDE = "#00ACFC"
AMARILLO = "#F7DC6F"
NARANJA = "#FFA500"

cont_figs = 1
fig = None
ax = None


# XXX: http://code.activestate.com/recipes/576527-freeze-make-any-object-immutable/
@total_ordering
class Frozen(object):

    def __init__(self, value):
        assert isinstance(value, LineString) or isinstance(value, Polygon), "{} no es la instancia esperada".format(type(value))
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
        return "i si acemos un acaca " + self._value.__str__()

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

# XXX: http://code.activestate.com/recipes/577197-sortedcollection/
class SortedCollection(object):
    '''Sequence sorted by a key function.

    SortedCollection() is much easier to work with than using bisect() directly.
    It supports key functions like those use in sorted(), min(), and max().
    The result of the key function call is saved so that keys can be searched
    efficiently.

    Instead of returning an insertion-point which can be hard to interpret, the
    five find-methods return a specific item in the sequence. They can scan for
    exact matches, the last item less-than-or-equal to a key, or the first item
    greater-than-or-equal to a key.

    Once found, an item's ordinal position can be located with the index() method.
    New items can be added with the insert() and insert_right() methods.
    Old items can be deleted with the remove() method.

    The usual sequence methods are provided to support indexing, slicing,
    length lookup, clearing, copying, forward and reverse iteration, contains
    checking, item counts, item removal, and a nice looking repr.

    Finding and indexing are O(log n) operations while iteration and insertion
    are O(n).  The initial sort is O(n log n).

    The key function is stored in the 'key' attibute for easy introspection or
    so that you can assign a new key function (triggering an automatic re-sort).

    In short, the class was designed to handle all of the common use cases for
    bisect but with a simpler API and support for key functions.

    >>> from pprint import pprint
    >>> from operator import itemgetter

    >>> s = SortedCollection(key=itemgetter(2))
    >>> for record in [
    ...         ('roger', 'young', 30),
    ...         ('angela', 'jones', 28),
    ...         ('bill', 'smith', 22),
    ...         ('david', 'thomas', 32)]:
    ...     s.insert(record)

    >>> pprint(list(s))         # show records sorted by age
    [('bill', 'smith', 22),
     ('angela', 'jones', 28),
     ('roger', 'young', 30),
     ('david', 'thomas', 32)]

    >>> s.find_le(29)           # find oldest person aged 29 or younger
    ('angela', 'jones', 28)
    >>> s.find_lt(28)           # find oldest person under 28
    ('bill', 'smith', 22)
    >>> s.find_gt(28)           # find youngest person over 28
    ('roger', 'young', 30)

    >>> r = s.find_ge(32)       # find youngest person aged 32 or older
    >>> s.index(r)              # get the index of their record
    3
    >>> s[3]                    # fetch the record at that index
    ('david', 'thomas', 32)

    >>> s.key = itemgetter(0)   # now sort by first name
    >>> pprint(list(s))
    [('angela', 'jones', 28),
     ('bill', 'smith', 22),
     ('david', 'thomas', 32),
     ('roger', 'young', 30)]

    '''

    def __init__(self, iterable=(), key=None):
        self._given_key = key
        key = (lambda x: x) if key is None else key
        decorated = sorted((key(item), item) for item in iterable)
        self._keys = [k for k, item in decorated]
        self._items = [item for k, item in decorated]
        self._key = key

    def _getkey(self):
        return self._key

    def _setkey(self, key):
        if key is not self._key:
            self.__init__(self._items, key=key)

    def _delkey(self):
        self._setkey(None)

    key = property(_getkey, _setkey, _delkey, 'key function')

    def clear(self):
        self.__init__([], self._key)

    def copy(self):
        return self.__class__(self, self._key)

    def __len__(self):
        return len(self._items)

    def __getitem__(self, i):
        return self._items[i]

    def __iter__(self):
        return iter(self._items)

    def __reversed__(self):
        return reversed(self._items)

    def __repr__(self):
        return '%s(%r, key=%s)' % (
            self.__class__.__name__,
            self._items,
            getattr(self._given_key, '__name__', repr(self._given_key))
        )

    def __reduce__(self):
        return self.__class__, (self._items, self._given_key)

    def __contains__(self, item):
        k = self._key(item)
        i = bisect_left(self._keys, k)
        j = bisect_right(self._keys, k)
        return item in self._items[i:j]

    def index(self, item):
        'Find the position of an item.  Raise ValueError if not found.'
        k = self._key(item)
        i = bisect_left(self._keys, k)
        j = bisect_right(self._keys, k)
        return self._items[i:j].index(item) + i

    def count(self, item):
        'Return number of occurrences of item'
        k = self._key(item)
        i = bisect_left(self._keys, k)
        j = bisect_right(self._keys, k)
        return self._items[i:j].count(item)

    def insert(self, item):
        'Insert a new item.  If equal keys are found, add to the left'
        k = self._key(item)
        i = bisect_left(self._keys, k)
        self._keys.insert(i, k)
        self._items.insert(i, item)

    def insert_right(self, item):
        'Insert a new item.  If equal keys are found, add to the right'
        k = self._key(item)
        i = bisect_right(self._keys, k)
        self._keys.insert(i, k)
        self._items.insert(i, item)

    def remove(self, item):
        'Remove first occurence of item.  Raise ValueError if not found'
        i = self.index(item)
        del self._keys[i]
        del self._items[i]

    def find(self, k):
        'Return first item with a key == k.  Raise ValueError if not found.'
        i = bisect_left(self._keys, k)
        if i != len(self) and self._keys[i] == k:
            return self._items[i]
        raise ValueError('No item found with key equal to: %r' % (k,))

    def find_le(self, k):
        'Return last item with a key <= k.  Raise ValueError if not found.'
        i = bisect_right(self._keys, k)
        if i:
            return self._items[i - 1]
        raise ValueError('No item found with key at or below: %r' % (k,))

    def find_lt(self, k):
        'Return last item with a key < k.  Raise ValueError if not found.'
        i = bisect_left(self._keys, k)
        if i:
            return self._items[i - 1]
        raise ValueError('No item found with key below: %r' % (k,))

    def find_ge(self, k):
        'Return first item with a key >= equal to k.  Raise ValueError if not found'
        i = bisect_left(self._keys, k)
        if i != len(self):
            return self._items[i]
        raise ValueError('No item found with key at or above: %r' % (k,))

    def find_gt(self, k):
        'Return first item with a key > k.  Raise ValueError if not found'
        i = bisect_right(self._keys, k)
        if i != len(self):
            return self._items[i]
        raise ValueError('No item found with key above: %r' % (k,))
    
    def encuentra_interfalo(self, llave_inicial, llave_final):
        k_ini = self._key(llave_inicial)
        k_fin = self._key(llave_final)
        assert k_ini <= k_fin
        i = bisect_left(self._keys, k_ini)
        logger_cagada.debug("todo el arreglo {}".format(self._keys))
        logger_cagada.debug("todo el arreglo {}".format(self._items))
        logger_cagada.debug("llave ini {} llave fin {}".format(k_ini, k_fin))
        logger_cagada.debug("pos ini {}".format(i))
        valores = []
        if i != len(self):
            j = bisect_right(self._keys, k_fin)
            logger_cagada.debug("pos fin {} lo q ai {}".format(j, self._items[j - 1]))
            assert i <= j
            valores = self._items[i:j]
        return valores

# XXX: https://www.accelebrate.com/blog/named-tuples-in-python/
conjunto_ordenado_dimension = namedtuple("conjunto_ordenado_dimension", "idx_dimension llave")
class conjunto_ordenado_en_dimensiones():
    def __init__(self, dimensiones):
        self.dimensiones = dimensiones
        self.conjunto = None
        self.conjuntos_ordenados = {}
        self._inicializa_conjuntos_ordenados()
    
    def _inicializa_conjuntos_ordenados(self):
        if not self.dimensiones:
            self.conjuntos_ordenados[0] = SortedCollection([])
        else:
            for idx_dime, llave in self.dimensiones:
                self.conjuntos_ordenados[idx_dime] = SortedCollection([], key=llave)
    
    def insertar(self, valor):
        for colleccion in self.conjuntos_ordenados.values():
            colleccion.insert(valor)
    
    def encuentra_interfalo(self, valor_inicial, valor_final, idx_dimension):
        conjunto = self.conjuntos_ordenados[idx_dimension]
        valores = conjunto.encuentra_interfalo(valor_inicial, valor_final)
        logger_cagada.debug("de llave {} se consigio ladim {} i los res {}".format(idx_dimension, conjunto.key, valores))
        return valores

def posicion_suma(pos_1, pos_2):
    return [pos_1[0] + pos_2[0], pos_1[1] + pos_2[1]]


def posicion_resta(pos_1, pos_2):
    return posicion_suma(pos_1, (-pos_2[0], -pos_2[1]))

    
def posicion_a_posicion_polar(centro, posi):
    logger_cagada.debug("calculando pos pol de {}".format(posi))
    distancias = posicion_resta(posi, centro)
    logger_cagada.debug("las distancias {}".format(distancias))
    angulo = arctan2(distancias[0], distancias[1])
    radio = sqrt(pow(distancias[0], 2) + pow(distancias[1], 2))
    logger_cagada.debug("pos pol es {}".format([degrees(angulo), radio]))
    return [angulo, radio]


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
        
    return all(map(lambda posi:posi[idx_dimension] == referencia or abs(posi[idx_dimension] - referencia) < 1e-8, posiciones))

def posicion_valida_colinearidad(posiciones):
    logger_cagada.debug("validando pos {}".format(posiciones))
    return posicion_valida_dimension_igual(posiciones, True) or posicion_valida_dimension_igual(posiciones, False)

# XXX: https://stackoverflow.com/questions/21291725/determine-if-shapely-point-is-within-a-linestring-multilinestring
def puto_valida_colinearidad(segmento, putos):
    return all(map(lambda puto:segmento.distance(puto) < 1e-8, putos))

# XXX: https://stackoverflow.com/questions/16739290/composing-functions-in-python
def caca_comun_composicion(functiones):
#    logger_cagada.debug("compocaca")
    caca = partial(reduce, compose)(reversed(functiones))
#    logger_cagada.debug("se regresa {}".format(caca))
    return caca

def posicion_chicharronera_intersexion(dx, dy, dr, D, raiz_discriminante, calcular_x):
  if calcular_x:
    term_var_1 = D * dy
    term_var_2 = dx * (-1 if dy < 0 else 1)
  else:
    term_var_1 = -D * dx
    term_var_2 = abs(dy)
  
  formula = lambda opera:(opera(term_var_1, term_var_2 * raiz_discriminante) / pow(dr, 2))
  res1, res2 = formula(add), formula(sub)
  return res1, res2

# XXX: http://mathworld.wolfram.com/Circle-LineIntersection.html
def posicion_calcula_intersexion_segmento_con_cir_culo(extremo_segmento_1, extremo_segmento_2, centro, radio):
    logger_cagada.debug("sacando intersex de {}-{} con circ {} r {}".format(extremo_segmento_1, extremo_segmento_2, centro, radio))
    ext_seg_1 = posicion_resta(extremo_segmento_1, centro)
    ext_seg_2 = posicion_resta(extremo_segmento_2, centro)
    logger_cagada.debug("segmento normalizado {}-{}".format(ext_seg_1, ext_seg_2))
    x1 = ext_seg_1[0]
    y1 = ext_seg_1[1]
    x2 = ext_seg_2[0]
    y2 = ext_seg_2[1]
    dx = x2 - x1
    dy = y2 - y1
    logger_cagada.debug("distancias x {} y {}".format(dx, dy))
    dr = sqrt(pow(dx, 2) + pow(dy, 2))
    logger_cagada.debug("dr {}".format(dr))
    D = x1 * y2 - x2 * y1
    discriminante = pow(radio, 2) * pow(dr, 2) - pow(D, 2)
    raiz_discriminante = sqrt(discriminante)
    intersexiones = None
    if discriminante >= 0:
        interx = posicion_chicharronera_intersexion(dx, dy, dr, D, raiz_discriminante, True)
        intery = posicion_chicharronera_intersexion(dx, dy, dr, D, raiz_discriminante, False)
        interx = list(map(lambda x:x + centro[0], interx))
        intery = list(map(lambda y:y + centro[1], intery))
        if discriminante == 0:
            assert interx[0] == interx[1]
            assert intery[0] == intery[1]
            interx = interx[0]
            intery = intery[0]
        intersexiones = list(zip(interx, intery))

    logger_cagada.debug("las intersex son {}".format(intersexiones))
    return intersexiones

def puto_intersexion_cir_culo_segmento(segmento, pos_centro, radio):
    intersex = None
    pos_seg = list(segmento.coords)
    pos_intersex = posicion_calcula_intersexion_segmento_con_cir_culo(pos_seg[0], pos_seg[1], pos_centro, radio)
    if pos_intersex is not None:
        if len(pos_intersex) == 1:
            intersex = Point(pos_intersex[0])
        else:
            assert len(pos_intersex) == 2
            intersex = LineString(map(Point, pos_intersex))

    return intersex

def caca_comun_redondea_intervalo_a_enteros(inter_1, inter_2, hacia_dentro):
    limite_menor, limite_mayor = sorted([inter_1, inter_2])
    
    limites = [limite_menor, limite_mayor]
    # XXX: https://stackoverflow.com/questions/44116557/how-to-extend-concatenate-two-iterators-in-python
    limites_posibles = list(sorted(chain(map(ceil, limites) , map(floor, limites))))
    
    if hacia_dentro:
        idx_lim_men = 1
        idx_lim_may = 2
    else:
        idx_lim_men = 0
        idx_lim_may = 3
    
#    logger_cagada.debug("intervalo {} se redondeo a {}".format([inter_1, inter_2], [limites_posibles[idx_lim_men], limites_posibles[idx_lim_may]]))
    
    return limites_posibles[idx_lim_men], limites_posibles[idx_lim_may]
    
def caca_comun_topa_intervalo_a_limites(inter, limite_1, limite_2):
        limite_1, limite_2 = sorted([limite_1, limite_2])
        filtro_pasa_banda = caca_comun_composicion([partial(max, limite_1), partial(min, limite_2)])
        intervalo = [inter[0], inter[1]]
        
        intervalo_topadp = list(map(filtro_pasa_banda, intervalo))
        logger_cagada.debug("intervalo {} se topo a {}".format(inter, intervalo_topadp))
        return intervalo_topadp

def caca_comun_redondea_intervalo_a_enteros_topado(inter_1, inter_2, limite_1, limite_2, hacia_dentro):
    caca = caca_comun_composicion([partial(caca_comun_redondea_intervalo_a_enteros, hacia_dentro=hacia_dentro), partial(caca_comun_topa_intervalo_a_limites, limite_1=limite_1, limite_2=limite_2)])
#    caca = partial(caca_comun_redondea_intervalo_a_enteros, hacia_dentro=hacia_dentro)
#    caca = partial(caca_comun_topa_intervalo_a_limites, limite_1=limite_1, limite_2=limite_2)
#    logger_cagada.debug("pero q mierda {}".format(caca))
    inter_1_ent, inter_2_ent = caca(inter_1=inter_1, inter_2=inter_2)
    return inter_1_ent, inter_2_ent
    

class sektor_cir_culo():

    angulos_extremos_verticales = [np.pi / 2, -np.pi / 2]
    angulos_extremos_horizontales = [0, np.pi]
    
    def __init__(self, centro, radio, pos_lim_1, pos_lim_2, lim_x, lim_y, putos_ordenados, mapa_posicion_a_forma, mapa_celda_a_poligono):
        self.centro = centro
        self.radio = radio
        self.pos_lim_1 = min(pos_lim_1, pos_lim_2)
        self.pos_lim_2 = max(pos_lim_1, pos_lim_2)
        self.lim_x = lim_x
        self.lim_y = lim_y
        self.putos_ordenados = putos_ordenados
        self.mapa_posicion_a_forma = mapa_posicion_a_forma
        self.mapa_celda_a_poligono = mapa_celda_a_poligono
        self.angulo = None
        self.segmento_sektor_1 = None
        self.segmento_sektor_2 = None
        self.circulo = None
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
        self._inicializa()
    
    @classmethod
    def calcula_angulo_sektor(cls, centro, pos_lim_1, pos_lim_2):
        pos_lim_1_int = posicion_resta(pos_lim_1, centro)
        pos_lim_2_int = posicion_resta(pos_lim_2, centro)
        centro_int = posicion_resta(centro, centro)
        assert centro_int == [0, 0]
        
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
        assert not posicion_valida_colinearidad([self.centro, self.pos_lim_1, self.pos_lim_2])
        self.angulo = sektor_cir_culo.calcula_angulo_sektor(self.centro, self.pos_lim_1, self.pos_lim_2)
        self.pos_final_segmento_1 = pos_final_segmento_1 = sektor_cir_culo.genera_segmento_linea(self.centro, self.pos_lim_1, self.radio)
        self.pos_final_segmento_2 = pos_final_segmento_2 = sektor_cir_culo.genera_segmento_linea(self.centro, self.pos_lim_2, self.radio)
        
        self.pos_pol_final_segmento_1 = self.posicion_a_posicion_polar_de_sektor(self.pos_final_segmento_1)
        self.pos_pol_final_segmento_2 = self.posicion_a_posicion_polar_de_sektor(self.pos_final_segmento_2)
        
        self.segmento_sektor_1 = LineString([Point(self.centro), Point(pos_final_segmento_1)])
        self.segmento_sektor_2 = LineString(map(Point, [self.centro, pos_final_segmento_2]))
        
        self.circulo = Point(self.centro).buffer(self.radio, 3600)
        
        self._inicializa_posiciones_polares()
        
        logger_cagada.debug("las posiciones polares de {} {} son {} {}".format(self.segmento_sektor_1, self.segmento_sektor_2, self.pos_pol_final_segmento_1, self.pos_pol_final_segmento_2))
        house_of_pain_pinta_figura(self.circulo)
        house_of_pain_pinta_figura(self.segmento_sektor_1)
        house_of_pain_pinta_figura(self.segmento_sektor_2)
        
        self._calcula_caja_sektor()
        self._determina_formas_tokadas()
    
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
            assert dif_angulos < np.pi
        
        pos_polar_seg_1[0] = angulo_1
        pos_polar_seg_2[0] = angulo_2
        self.pos_pol_final_segmento_1 = pos_polar_seg_1
        self.pos_pol_final_segmento_2 = pos_polar_seg_2
    
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
        
        logger_cagada.debug("las pos extr pot {}".format(posiciones_potenciales_extremas))
        
        if calcula_extremos_verticales:
            angulos_extremos = sektor_cir_culo.angulos_extremos_verticales
            idx_dimension = 0
        else:
            angulos_extremos = sektor_cir_culo.angulos_extremos_horizontales
            idx_dimension = 1
        
        pos_pol_extremas = map(lambda angulo:(angulo, self.radio), angulos_extremos)
#        logger_cagada.debug("las pos pol extr {}".format(list(pos_pol_extremas)))
        
        pos_pol_extremas_en_sektor = list(filter(self.posicion_polar_en_cuerda_sektor, pos_pol_extremas))
        
        assert posicion_valida_colinearidad(list(map(self.posicion_polar_a_posicion_de_sektor, pos_pol_extremas_en_sektor)) + posiciones_potenciales_extremas) or len(pos_pol_extremas_en_sektor) < 2, "los pos pol extre en sektor son {}".format(list(map(self.posicion_polar_a_posicion_de_sektor, pos_pol_extremas_en_sektor)))
        
        if len(pos_pol_extremas_en_sektor):
            posiciones_potenciales_extremas.append(self.posicion_a_posicion_polar_de_sektor(pos_pol_extremas_en_sektor[0]))
        
        limite_max = max(posiciones_potenciales_extremas, key=lambda posi:posi[idx_dimension])[idx_dimension]
        limite_min = min(posiciones_potenciales_extremas, key=lambda posi:posi[idx_dimension])[idx_dimension]
        
        return limite_min, limite_max
    
    def _calcula_limites_enteros_topados_por_caja_sektor(self, limite_1, limite_2, para_x, hacia_dentro):
        if hacia_dentro:
            offset = 2
        else:
            offset = 1
        if para_x:
            limite_superior_dim_var = self.lim_x - offset
        else:
            limite_superior_dim_var = self.lim_y - offset
        
        return caca_comun_redondea_intervalo_a_enteros_topado(limite_1, limite_2, offset, limite_superior_dim_var, hacia_dentro)
    
    def _calcula_caja_sektor(self):
        extremos_caja_vertical = self._calcula_extremos_caja(True)
        extremos_caja_horizontal = self._calcula_extremos_caja(False)
        self.extremo_caja_vertical_min, self.extremo_caja_vertical_max = extremos_caja_vertical
        self.extremo_caja_horizontal_min, self.extremo_caja_horizontal_max = extremos_caja_horizontal
        logger_cagada.debug("los extremos v {},{} h {},{}".format(self.extremo_caja_vertical_min, self.extremo_caja_vertical_max, self.extremo_caja_horizontal_min, self.extremo_caja_horizontal_max))
        logger_cagada.debug("lim x {} y {}".format(self.lim_x, self.lim_y))
        
        limites_x = list(sorted(self._calcula_limites_enteros_topados_por_caja_sektor(extremos_caja_vertical[0], extremos_caja_vertical[1], True, True)))
        limites_y = list(sorted(self._calcula_limites_enteros_topados_por_caja_sektor(extremos_caja_horizontal[0], extremos_caja_horizontal[1], False, True)))
        limites_x[1] += 1
        limites_y[1] += 1
        self.valores_x_abarcados = range(*limites_x)
        self.valores_y_abarcados = range(*limites_y)
        
        logger_cagada.debug("valors x {} y {}".format(self.valores_x_abarcados, self.valores_y_abarcados))
    
    def _determina_celdas_tokadas_por_caja(self):
        celdas_tokadas = []
        limites_valores_x = self._calcula_limites_enteros_topados_por_caja_sektor(self.extremo_caja_vertical_min, self.extremo_caja_vertical_max, True , False)
        limites_valores_y = self._calcula_limites_enteros_topados_por_caja_sektor(self.extremo_caja_horizontal_min, self.extremo_caja_horizontal_max, False, False)
        logger_cagada.debug("los limites de la caja son {} {}".format(limites_valores_x, limites_valores_y))
        
        for valor_abarcado_x in range(limites_valores_x[0], limites_valores_x[1]):
            for valor_abarcado_y in range(limites_valores_y[0], limites_valores_y[1]):
                celda = [float(valor_abarcado_x), float(valor_abarcado_y)]
                celdas_tokadas.append(tuple(celda))
                
        return celdas_tokadas
    
    def _determina_putos_celdas_tokados(self, para_x):
        putos_intersextados = []
        if para_x:
            valores_abarcados = self.valores_x_abarcados
            idx_dim_fija = 0
            idx_dim_var = 1
            dim_a_buscar = "x"
        else:
            valores_abarcados = self.valores_y_abarcados
            idx_dim_fija = 1
            idx_dim_var = 0
            dim_a_buscar = "y"
        
        for valor in valores_abarcados:
            logger_cagada.debug("en valor {}".format(valor))
            posi_inicio_raya = [0, 0]
            posi_fin_raya = [0, 0]
            posi_inicio_raya[idx_dim_fija] = posi_fin_raya[idx_dim_fija] = valor
            posi_inicio_raya[idx_dim_var] = self.centro[1] - self.radio
            posi_fin_raya[idx_dim_var] = self.centro[1] + self.radio
            raya = LineString((posi_inicio_raya, posi_fin_raya))
            intersexs = self.calcula_intersexiones_de_segmento_con_sektor(raya)
            logger_cagada.debug("las intersex {}".format(intersexs))
            putos_en_intersex = self.putos_ordenados.encuentra_interfalo(intersexs[0], intersexs[1], dim_a_buscar)
            logger_cagada.debug("los putos enc {}".format(putos_en_intersex))
            putos_intersextados += putos_en_intersex
        
        return putos_intersextados
    
    def _determina_formas_tokadas(self):
        putos_intersextados = []
        poligonos_tocados = set()
        celdas_tokadas = set()
        
        putos_intersextados_x = self._determina_putos_celdas_tokados(True)
        putos_intersextados_y = self._determina_putos_celdas_tokados(False)
        celdas_tokadas |= set(self._determina_celdas_tokadas_por_caja())
        logger_cagada.debug("las celdas tokadiscos son {}".format(celdas_tokadas))
        
        putos_intersextados.extend(putos_intersextados_x + putos_intersextados_y)
        
        for puto in putos_intersextados:
            formas = self.mapa_posicion_a_forma[puto]
            for forma in formas:
                logger_cagada.debug("de puto {} sale caca {}".format(puto, forma))
                poligonos_tocados.add(forma)
            house_of_pain_pinta_figura(Point(puto), ROJO)
        
        for celda in list(celdas_tokadas):
            if celda in self.mapa_celda_a_poligono:
                forma = self.mapa_celda_a_poligono[celda]
                poligonos_tocados.add(forma)
                cuadro = house_of_pain_crea_poligono_de_celda(celda)
                house_of_pain_pinta_figura(cuadro, NARANJA)
        
        for poli in poligonos_tocados:
            house_of_pain_pinta_figura(poli, AMARILLO)
        
    def normaliza_posicion_polar(self, pos_pol):
        logger_cagada.debug("normalizando {}".format(pos_pol))
        if self.normalizar_angulos:
            angulo_pol = pos_pol[0]
            if angulo_pol < 0:
                angulo_pol = np.pi * 2 + angulo_pol
            pos_pol[0] = angulo_pol
        return pos_pol
    
    def posicion_polar_dentro_de_sektor(self, pos_pol):
        logger_cagada.debug("compradano {} <= {} <= {}".format(self.pos_pol_final_segmento_1[0], pos_pol[0], self.pos_pol_final_segmento_2[0]))
        if self.pos_pol_final_segmento_1[0] <= pos_pol[0] <= self.pos_pol_final_segmento_2[0]:
            logger_cagada.debug("si sta dentro")
            return True
        else:
            return False
    
    def posicion_polar_en_cuerda_sektor(self, pos_pol):
        return self.posicion_polar_dentro_de_sektor(pos_pol) and abs(pos_pol[1] - self.radio) < 1e-8
    
    def posicion_en_cuerda_sektor(self, posi):
        return self.posicion_polar_en_cuerda_sektor(self.posicion_a_posicion_polar_de_sektor(posi))
    
    def puto_en_cuerda_sektor(self, puto):
        return self.posicion_en_cuerda_sektor(puto_a_posicion(puto))
    
    def posicion_a_posicion_polar_de_sektor(self, posi):
        logger_cagada.debug("sacando pos pol de {} con centro {}".format(posi, self.centro))
        return posicion_a_posicion_polar(self.centro, posi)
    
    def posicion_polar_a_posicion_de_sektor(self, pos_pol):
        return posicion_polar_a_posicion(self.centro, pos_pol)
    
    @property
    def segmentos_sektor(self):
        return [self.segmento_sektor_1, self.segmento_sektor_2]
    
    @property
    def posiciones_de_segmentos_sektor(self):
        return reduce(add, (map(lambda segmento:list(segmento.coords), self.segmentos_sektor)), [])
    
    def valida_colinearidad(self, putos):
        logger_cagada.debug("colineadirdad d  {}".format(putos))
        return posicion_valida_colinearidad(list(map(lambda puto:(puto.x, puto.y), putos))) or any(map(lambda segmento:puto_valida_colinearidad(segmento, putos), self.segmentos_sektor))
    
    def calcula_interseccion_con_arco_sektor(self, segmento):
        putos_intersex = []
        puto_tangencial = False
        logger_cagada.debug("calculando intersex de circ centro {} radio {} con seg {}".format(self.centro, self.radio, segmento))
 #       intersex = self.circulo.intersection(segmento)
        intersex = puto_intersexion_cir_culo_segmento(segmento, self.centro, self.radio)
        logger_cagada.debug("intersex es {}".format(intersex))
        assert not isinstance(intersex, GeometryCollection)
        if not isinstance(intersex, Point):
            assert isinstance(intersex, LineString)
            
#            putos_intersex = map(lambda pos_pol:Point(self.posicion_polar_a_posicion_de_sektor(pos_pol)), filter(lambda pos_pol:self.posicion_polar_en_cuerda_sektor(pos_pol), map(lambda posi:self.normaliza_posicion_polar(self.posicion_a_posicion_polar_de_sektor(posi)), list(intersex.coords))))
            putos_intersex = list(map(Point, filter(caca_comun_composicion((self.posicion_a_posicion_polar_de_sektor, self.normaliza_posicion_polar, self.posicion_polar_en_cuerda_sektor)) , list(intersex.coords))))
            
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
        
        logger_cagada.debug("los putos d intersex con el arco {}".format(putos_intersex))
        
        assert not putos_intersex or self.valida_colinearidad(putos_intersex)
        
        return putos_intersex, puto_tangencial
    
    def calcula_interseccion_con_segmento_de_linea(self, segmento, con_segmento_sektor_1):
        if con_segmento_sektor_1:
            segmento_sektor = self.segmento_sektor_1
        else:
            segmento_sektor = self.segmento_sektor_2
        
        logger_cagada.debug("calculando intersex de seg {} con seg {}".format(segmento, segmento_sektor))
        intersex = segmento.intersection(segmento_sektor)
        logger_cagada.debug("la intersex es {}".format(intersex))
        
        assert (not isinstance(intersex, LineString))
        
        if isinstance(intersex, Point):
            return intersex
        else:
            return None
    
    def calcula_intersexiones_de_segmento_con_sektor(self, segmento):
        putos_intersexion = []

        logger_cagada.debug("calculando intersex de {}".format(segmento))
        intersex_seg_sektor_1 = self.calcula_interseccion_con_segmento_de_linea(segmento, True)
        intersex_seg_sektor_2 = self.calcula_interseccion_con_segmento_de_linea(segmento, False)
        intersexs_arco, puto_tangencial = self.calcula_interseccion_con_arco_sektor(segmento)
        logger_cagada.debug("intersex con sektor 1 {}".format(intersex_seg_sektor_1))
        logger_cagada.debug("intersex con sektor 2 {}".format(intersex_seg_sektor_2))
        logger_cagada.debug("intersex con arco {},{}".format(intersexs_arco, puto_tangencial))
        
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
        
        return list(map(puto_a_posicion, putos_intersexion))
    
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
    return Frozen(poligono)

def house_of_pain_crea_poligono_de_celda(celda):
    lineas_cerda = house_of_pain_genera_lineas_de_celda(celda)
    policaca = house_of_pain_crea_poligono_de_lineas(lineas_cerda)
    return policaca
    
def house_of_pain_genera_poligono_y_putos_de_celda_dfs(matrix, celda_inicial, mapa_posicion_a_forma, mapa_linea_a_forma, mapa_celda_a_poligono, celdas_ya_visitadas, putos_ordenados):
    pila = [celda_inicial]
    duenio = celda_inicial
    lineas_poligono = set()
    celdas_ya_visitadas_int = set()
    celdas_ya_visitadas_int.add(celda_inicial)
    while pila:
        celda_actual = pila.pop()
#        logger_cagada.debug("celda act {}".format(celda_actual))
        lineas = house_of_pain_genera_lineas_de_celda(celda_actual)
        assert len(lineas) == 4
#        logger_cagada.debug("las celdas ia visitadas {}".format(celdas_ya_visitadas_int))
        assert lineas[0] != lineas[1]
        assert lineas[0] < lineas[1]
        for linea in lineas:
 #           logger_cagada.debug("linea act {}".format(linea))
# XXX: https://stackoverflow.com/questions/20474549/extract-points-coordinates-from-python-shapely-polygon
            putos = list(linea.coords)
            assert len(putos) == 2
            
            if linea in lineas_poligono:
                lineas_poligono.remove(linea)
#                logger_cagada.debug("linea removida")
            else:
                lineas_poligono.add(linea)
 #               logger_cagada.debug("linea anadida")
                
#        logger_cagada.debug("las lineas aora son {}".format(sorted(list(lineas_poligono))))
        
        for celda_aledana in map(lambda mov: (celda_actual[0] + mov[0], celda_actual[1] + mov[1]), house_of_pain_movimientos_matrix):
#            logger_cagada.debug("proc celda ale {}".format(celda_aledana))
            if celda_aledana not in celdas_ya_visitadas_int and matrix[celda_aledana[0]][celda_aledana[1]] == '#':
#                logger_cagada.debug("anadida")
                pila.append(celda_aledana)
                celdas_ya_visitadas_int.add(celda_aledana)
    
#    logger_cagada.debug("las lineas del pol {}".format(lineas_poligono))
    
    poligono = house_of_pain_crea_poligono_de_lineas(lineas_poligono)
    
    for linea in list(lineas_poligono):
        mapa_linea_a_forma[linea] = poligono
        for puto in list(linea.coords):
            mapa_posicion_a_forma[puto].append(poligono)
            putos_ordenados.insertar(puto)
    
    
    for celda_visitada in list(celdas_ya_visitadas_int):
        mapa_celda_a_poligono[celda_visitada] = poligono
        celdas_ya_visitadas.add(celda_visitada)
    
# XXX: https://gis.stackexchange.com/questions/216745/get-polygon-shapefile-in-python-shapely-by-clipping-linearring-with-linestring
    logger_cagada.debug("de celda {} se generaro poligono {} de cirds {}".format(celda_inicial, poligono, list(poligono.exterior.coords)))
    logger_cagada.debug("cekdas invol {}".format(sorted(list(celdas_ya_visitadas_int))))
    
    house_of_pain_pinta_figura(poligono)
    
def house_of_pain_pinta_figura(figura, color=GRAY):
    global cont_figs
    global fig
    global ax
    if isinstance(figura, Frozen):
        figura = figura._value
    logger_cagada.debug("intentando pintar {}".format(figura))
    if isinstance(figura, Polygon):
        figura = rotate(figura, -90, origin=(0, 0))
        poly = mapping(figura)
#    poly = {"type": "Polygon", "coordinates": mapping(figura)["coordinates"]}
        patch = PolygonPatch(poly, fc=color, ec=BLUE, alpha=0.5, zorder=2)
        ax.add_patch(patch)
    else:
        if isinstance(figura, LineString):
            figura = rotate(figura, -90, origin=(0, 0))
            x, y = figura.xy
            ax.plot(x, y, color=color, linewidth=1, solid_capstyle='round', zorder=1)
        else:
            if isinstance(figura, Point):
# XXX: https://stackoverflow.com/questions/27779845/how-to-plot-one-single-data-point
                figura = rotate(figura, -90, origin=(0, 0))
                x, y = figura.xy
                ax.plot(x, y, color=color, marker='o', markersize=3)
    ax.set_xlim(xmin=-10, xmax=10)
    ax.set_ylim(ymin=-10, ymax=10)
    cont_figs += 1

                    
def house_of_pain_genera_poligonos_y_putos(matrix, mapa_posicion_a_forma, mapa_linea_a_forma, mapa_celda_a_poligono, putos_ordenados):
    celdas_ya_visitadas = set()
    for i, fila in enumerate(matrix[1:-1], 1):
        logger_cagada.debug("fila es {}".format(fila))
        for j, cerda in enumerate(fila[1:-1], 1):
            celda_act = (i, j)
            if(cerda == "#" and celda_act not in celdas_ya_visitadas):
                house_of_pain_genera_poligono_y_putos_de_celda_dfs(matrix, celda_act, mapa_posicion_a_forma, mapa_linea_a_forma, mapa_celda_a_poligono, celdas_ya_visitadas, putos_ordenados)

                
def house_of_pain_core(matrix, pos_inicio):
    mapa_posicion_a_forma = defaultdict(lambda:[])
    mapa_linea_a_forma = {}
    mapa_celda_a_poligono = {}
    putos_ordenados = conjunto_ordenado_en_dimensiones([conjunto_ordenado_dimension("x", None), conjunto_ordenado_dimension("y", lambda posi:(posi[1], posi[0]))])
    
    lim_x = len(matrix)
    lim_y = len(matrix[0])
   
    global fig
    global ax
    fig = pyplot.figure(1, dpi=180)
    ax = fig.add_subplot(121)
    house_of_pain_genera_poligonos_y_putos(matrix, mapa_posicion_a_forma, mapa_linea_a_forma, mapa_celda_a_poligono, putos_ordenados)
    
    sektor = sektor_cir_culo(pos_inicio, 8, (2, 2), (2, 3), lim_x, lim_y, putos_ordenados, mapa_posicion_a_forma, mapa_celda_a_poligono)
    
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
        pos_inicio = None
        for i in range(alto):
            linea = stdin.readline().strip() 
            if "X" in linea:
                j = linea.index("X")
                pos_inicio = (i + 1.5, j + 1.5)
            matrix.append("%" + linea + "%")
        matrix.append("%" * (ancho + 2))
        
        assert pos_inicio
        logger_cagada.debug("la matrix leida\n{}".format(caca_comun_imprime_matrix(matrix)))
        logger_cagada.debug("la pos ini {}".format(pos_inicio))
        house_of_pain_core(matrix, pos_inicio)


if __name__ == '__main__':
        FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
        logging.basicConfig(level=nivel_log, format=FORMAT)
        logger_cagada = logging.getLogger("asa")
        logger_cagada.setLevel(nivel_log)
        house_of_pain_main()
