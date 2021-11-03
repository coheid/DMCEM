import argparse
import os 

from library.model import *

base   = os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser(description='Market simulation.')
parser.add_argument('path', metavar='path', type=str,
                    help='path to a valid configuration')
parser.add_argument('-d'  , metavar='base', type=str, default=base,
                    help='path to the base directory')
args = parser.parse_args()

m = Model("Model", args.path, args.d)
m.procedure()
