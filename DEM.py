#! /usr/bin/env python3
#

import os.path
import os
import bz2
import numpy as np
from scipy.interpolate import RectBivariateSpline

class DEM:
    def __init__(self, fn:str) -> None:
        if not os.path.isfile(fn): raise Exception('File {} does not exist'.format(fn))
        (root, ext) = os.path.splitext(fn) # Split off the extension
        binfn = root + ".npz" # Numpy binary file type
        if os.path.isfile(binfn) and (os.stat(fn).st_mtime < os.stat(binfn).st_mtime):
            (x, y, z) = self.loadBinary(binfn)
        else:
            (x, y, z) = self.loadESRI(fn)
            self.saveBinary(binfn, x, y, z)
        self.__func = RectBivariateSpline(y, x, z)

    def loadBinary(self, fn:str) -> None: # Load the binary 
        info = np.load(fn)
        return (info['x'], info['y'], info['z'])

    def saveBinary(self, fn:str, x:np.array, y:np.array, z:np.array) -> None: # ->Binary file
        np.savez_compressed(fn, x=x, y=y, z=z)


    def loadESRI(self, fn): # Load a GEBCO ESRI file
        z = None
        with bz2.open(fn) as fp:
            intKeys = {'ncols', 'nrows'}
            info = {}
            cnt = 0
            rowIndex = 0
            for line in fp:
                cnt += 1
                fields = line.split() # Split on whitespace
                fields = [x.decode('utf-8') for x in fields]
                if len(fields) == 2: # Key/Value pair
                    (key, val) = fields
                    if key in intKeys:
                        val = int(val)
                    else:
                        val = float(val)
                    info[key] = val
                elif ('ncols' in info) and (len(fields) == info['ncols']):
                    if z is None:
                        z = np.zeros((info['nrows'], info['ncols']))
                    z[rowIndex] = np.asarray(fields)
                    rowIndex += 1
                else:
                    raise Exception('Unrecognized line({}) in {}'.format(cnt, fn))

        for item in ['xllcorner', 'yllcorner', 'nrows', 'ncols', 'cellsize']:
            if item not in info: raise Exception(item + ' not in ' + fn)
        if z is None: raise Exception('Data is empty in ' + fn)

        delta = info['cellsize']
        nx = info['ncols']
        ny = info['nrows']
        xll = info['xllcorner']
        yll = info['yllcorner']

        x = np.arange(xll, xll + delta * (nx - 0.5), delta)
        y = np.arange(yll, yll + delta * (ny - 0.5), delta)
        return (x, y, z)

    def depth(self, lat:float, lon:float) -> float:
        return None if self.__func is None else -self.__func(lat, lon)[0][0]
