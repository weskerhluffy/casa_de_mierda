import pytest
from rayo.castring import Poligono, nivel_log
import logging
from rayo import castring

# XXX: https://docs.pytest.org/en/latest/fixture.html#fixture
@pytest.fixture(scope="module")
def mierda():
    if not castring.logger_cagada:
        FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
        logging.basicConfig(level=nivel_log, format=FORMAT)
        castring.logger_cagada = logging.getLogger("asa")
        castring.logger_cagada.setLevel(nivel_log)

# TODO: https://stackoverflow.com/questions/21824157/how-to-extract-interior-polygon-coordinates-using-shapely
# XXX: https://stackoverflow.com/questions/26405380/how-do-i-correctly-setup-and-teardown-my-pytest-class-with-tests
# XXX: http://pythontesting.net/framework/pytest/pytest-fixtures-easy-example/
@pytest.fixture
def caca():
    putos = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
    caca = Poligono(putos)
    return caca

@pytest.fixture
def caca_hueca():
    putos_e = [(0, 0), (0, 8), (4, 8), (4, 0), (0, 0)]
    putos_i = [[(1, 1), (1, 3), (3, 3), (3, 1)], [(1, 5), (1, 7), (3, 7), (3, 5)]]
    caca = Poligono(putos_e, putos_i)
    return caca

putos_chuecos = [
        (0, 0),
        (100, 0),
        (100, 100),
        (0, 100),
        (0, 20),
        (20, 20),
        (20, 30),
        (10, 30),
        (10, 70),
        (70, 70),
        (70, 50),
        (50, 50),
        (50, 55),
        (40, 55),
        (40, 25),
        (55, 25),
        (55, 30),
        (69, 30),
        (69, 10),
        (59, 10),
        (59, 15),
        (56, 15),
        (56, 17),
        (48, 17),
        (48, 15),
        (41, 15),
        (41, 10),
        (0, 10),
        (0, 0)
        ]
@pytest.fixture
def caca_chueca():
    putos_e = putos_chuecos
    putos_i = []
    caca = Poligono(putos_e, putos_i)
    return caca
    

def test_caca(mierda, caca):
    assert caca
    
def test_caca_hueca(mierda, caca_hueca):
    exte, inte = caca_hueca.extract_poly_coords()
    assert exte and len(inte) == 2      
    
def test_caca_chueca(mierda, caca_chueca):
    exte, inte = caca_chueca.extract_poly_coords()
#    castring.logger_cagada.debug("las mierdas {}".format(exte))
    assert exte and not len(inte) and exte == putos_chuecos
