'''
Created on 12/01/2018

@author: ernesto
'''
import unittest
from rayo.castring import conjunto_ordenado_en_dimensiones, \
    conjunto_ordenado_dimension, nivel_log
import logging
from rayo import castring

class Test(unittest.TestCase):

    def setUp(self):
        FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
        logging.basicConfig(level=nivel_log, format=FORMAT)
        castring.logger_cagada = logging.getLogger("asa")
        castring.logger_cagada.setLevel(nivel_log)
        self.mierda = conjunto_ordenado_en_dimensiones([conjunto_ordenado_dimension("x", None), conjunto_ordenado_dimension("y", lambda posi:(posi[1], posi[0]))])
        for i in range(0, 10):
            for j in range(0, 10):
                self.mierda.insertar((i, j))

    def testCacaX(self):
        ass = self.mierda.encuentra_interfalo((3, 4), (3, 9), "x")
        self.assertTrue(all(map(lambda posi:posi[0] == 3, ass)), "puta veerga lo q regreso {}".format(ass))

    def testCacaY(self):
        ass = self.mierda.encuentra_interfalo((4, 3), (9, 3), "y")
        self.assertTrue(all(map(lambda posi:posi[1] == 3, ass)), "puta veerga lo q regreso {}".format(ass))

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
