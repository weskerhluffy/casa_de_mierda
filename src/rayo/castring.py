'''
Created on 18/12/2017

@author: 

XXX: https://code.google.com/codejam/contest/1460488/dashboard#s=p3&a=3
'''

import logging
import sys

nivel_log = logging.ERROR
nivel_log = logging.DEBUG
logger_cagada = None

if __name__ == '__main__':
        FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
        logging.basicConfig(level=nivel_log, format=FORMAT)
        logger_cagada = logging.getLogger("asa")
        logger_cagada.setLevel(nivel_log)