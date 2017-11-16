from __future__ import absolute_import, division, print_function

import os
import sys

path_already_set = False
pyPlanetPath = os.getenv('PYPLANETPATH')
if pyPlanetPath is None:
    print('No PYPLANETPATH environment variable, trying .paths.json')
    import code_path
    pyPlanetPath = code_path.set_and_get('pyPlanet')
    if pyPlanetPath is None:
        print("No .paths.json file either.  Setting to ./")
        pyPlanetPath = './'
    else:
        path_already_set = True
if not path_already_set:
    sys.path.append(pyPlanetPath)
