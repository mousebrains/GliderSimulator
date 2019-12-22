# Load JSON parameter files

from logging import Logger
from argparse import ArgumentParser
import json
import math
import os.path

class BaseParam:
    def __init__(self, fn:str, logger:Logger, outdir:str) -> None:
        self.logger = logger
        if fn is None: 
            self.info = None
            return
        with open(fn, 'rb') as fp:
            content = fp.read() # Load as a string
            self.info = json.loads(content)
            if outdir is not None:
                ofn = os.path.join(outdir, os.path.basename(fn))
                with open(ofn, 'wb') as ofp:
                    ofp.write(content)

    def __repr__(self):
        return str(self.info)

class Glider(BaseParam):
    def __init__(self, fn:str, logger:Logger, outdir:str) -> None:
        BaseParam.__init__(self, fn, logger, outdir)
        self.intervalIndex = 0 # Which index to use in the call_interval list
        if self.info is not None:
            self.surface_affine = self.info['surface_inflate_time'] \
                    + self.info['surface_pre_gps_time'] \
                    + self.info['surface_call_setup_time'] \
                    + self.info['surface_call_drop_time'] \
                    + self.info['surface_post_gps_time'] \
                    + self.info['surface_deflate_time']

    def currentUseCorrection(self) -> bool:
        return False if self.info is None else self.info['current_use_correction']

    def currentMaxFraction(self) -> float:
        return 1 if self.info is None else self.info['current_max_perpendicular_fraction']

    def currentToFastFraction(self) -> float:
        return 0.5 if self.info is None else self.info['current_to_fast_fraction']

    def waypointApproachRadius(self) -> bool:
        return 10 if self.info is None else self.info['waypoint_approach_radius']

    def waypointSurfaceOnHitting(self) -> bool:
        return False if self.info is None else self.info['waypoint_surface_on']

    def yosBetweenSurfacings(self) -> float:
        return 1e20 if self.info is None else self.info['call_interval_yos']

    def average_horizontal_speed(self) -> float:
        """ Expected average horizontal speed """
        # From l_vert = t_dive * W_dive = t_climb * W_climb
        # l_horiz = t_dive * U_dive + t_climb * U_climb
        # <U> = l_horiz / (t_dive + t_climb)
        if self.info is None: return self.__calcU(0.2, 26 + 2)
        W_dive  = self.info['dive_Wg']
        W_climb = self.info['climb_Wg']
        U_dive = self.__calcU(W_dive, self.info['dive_pitch'] + self.info['dive_aoa'])
        U_climb= self.__calcU(W_climb, self.info['climb_pitch'] + self.info['climb_aoa'])
        num   = U_climb * W_dive + U_dive * W_climb
        denom = W_climb + W_dive
        return num / denom


        

    def nextSurfacingTime(self, tNow:float, tNext:float) -> float:
        """ When should the next surfacing time be """
        if self.info is None:
            dt = 7200 # Every 2 hours
        else:
            index = self.intervalIndex
            intervals = self.info['call_interval']
            n = len(intervals) - 1
            if index >= len(intervals): index = n
            self.intervalIndex = index + 1
            dt = intervals[index]
        if self.info['call_interval_type'] == 0:
            tNext = tNow
        while tNext <= tNow: tNext += dt
        return tNext

    def surfaceVolume(self) -> float:
        return 430 if self.info is None else self.info['surface_volume']

    def surfaceTime(self, dt:float) -> float:
        """ How long was the total surface time """
        if self.info is None: return 0
        callTime = (dt / 60) * self.info['surface_call_time_per_minute']
        return (self.surface_affine + callTime, callTime)

    def diveTime(self, sDepth:float, bathDepth:float, dtMax:float) -> float:
        """ Dive to bathDepth-altimeter or crush depth """
        if self.info is None: return self.__makeTimeDistance(sDepth, bathDepth-10, 0.2, 26, 2, 260)
        eDepth = min(self.info['max_depth'], bathDepth - self.info['min_altimeter'])
        self.logger.debug('sDepth=%s eDepth=%s', sDepth, eDepth)
        return self.__makeTimeDistance(sDepth, eDepth,
                self.info['dive_Wg'],
                self.info['dive_pitch'],
                self.info['dive_aoa'],
                self.info['dive_turn_time'],
                self.info['dive_turn_factor'],
                self.info['dive_volume'],
                dtMax)

    def climbTime(self, sDepth:float, dtMax:float) -> float:
        """ Climb to topDepth """
        if self.info is None: return self.__makeTimeDistance(sDepth, 10, 0.2, 26, 2, 260)
        return self.__makeTimeDistance(sDepth, self.info['top_depth'],
                self.info['climb_Wg'],
                self.info['climb_pitch'],
                self.info['climb_aoa'],
                self.info['climb_turn_time'],
                self.info['climb_turn_factor'],
                self.info['climb_volume'],
                dtMax)

    def climbSurfaceTime(self, sDepth:float, dtMax:float) -> float:
        """ Climb to the surface """
        if self.info is None: return self.__makeTimeDistance(sDepth, 0, 0.2, 26, 2, 260)
        return self.__makeTimeDistance(sDepth, 0, 
                self.info['surface_Wg'], 
                self.info['surface_pitch'], 
                self.info['surface_aoa'], 
                0,
                0,
                self.info['surface_volume'],
                dtMax)

    def __calcU(self, W:float, theta:float) -> float:
        """ Horizontal speed """
        return abs(W / math.tan(math.radians(theta)))

    def __makeTimeDistance(self, sDepth:float, eDepth:float, 
            W:float, pitch:float, aoa:float, 
            tTurn:float,
            turnFactor:float,
            volume:float, dtMax:float) -> tuple:
        """ Return a duration in seconds and distance in meters """
        dt = min(dtMax, abs(eDepth - sDepth) / abs(W)) # t to change the depth by dDepth at speed W
        dDepth = W * dt
        finalDepth = sDepth + dDepth * (1 if eDepth >= sDepth else -1)
        U = self.__calcU(W, pitch + aoa) # Horizontal speed
        return (dt + tTurn, (dt + tTurn*turnFactor) * U, volume, finalDepth)

class Energy(BaseParam):
    def __init__(self, fn:str, logger:Logger, outdir:str) -> None:
        BaseParam.__init__(self, fn, logger, outdir)
        # Get the available charge in Joules
        self._charge = 1e12 if self.info is None else \
                (3600 * self.info['battery_full_charge'] * self.info['available_fraction'])

    def __repr__(self) -> str:
        return str(round(self._charge, 0))

    def __bool__(self) -> bool:
        """ Is there any energy left or not? """
        return bool(self._charge > 0)

    def __isub__(self, delta:float):
        """ Subtract some energy usage """
        if self._charge is not None: self._charge -= delta
        return self

    def surface(self, dtSurface:float, dtCall:float) -> float:
        """ dtSurface is the time at the surface, and dtCall is the iridium call duration """
        return 0 if self.info is None else \
                self.info['air_pump'] \
                + 2 * self.info['battery_surface'] \
                + self.info['surface'] * dtSurface \
                + self.info['iridium'] * dtCall

    def dive(self, dt:float, dVol:float) -> float:
        """ dt is the dive duration, dVol is the buoyancy throw """
        return 0 if self.info is None else \
                self.info['battery_shift'] \
                + self.info['deep_buoyancy_per_CC'] * dVol \
                + self.info['diving'] * dt

    def climb(self, dt:float, dVol:float) -> float:
        """ dt is the climb duration, dVol is the buoyancy throw """
        return 0 if self.info is None else \
                self.info['battery_shift'] \
                + self.info['shallow_buoyancy_per_CC'] * dVol \
                + self.info['climbing'] * dt

class Waypoints(BaseParam):
    def __init__(self, fn:str, logger:Logger, outdir:str) -> None:
        BaseParam.__init__(self, fn, logger, outdir)
        self.index = 0
        self.n = None if self.info is None else len(self.info)
        
    def initial(self) -> dict:
        if self.info is None: return {'name':'GotMe', 'lat': 0, 'lon': 0}
        wpt = self.info[0]
        self.logger.info('Initial location %s', wpt)
        return wpt

    def next(self) -> dict:
        self.index += 1
        if self.index >= self.n: self.index = 0
        wpt = self.info[self.index]
        self.logger.info('Next waypoint%s', wpt)
        return wpt


class Mission(BaseParam):
    def __init__(self, fn:str, logger:Logger, outdir:str) -> None:
        BaseParam.__init__(self, fn, logger, outdir)

class Drift(BaseParam): # For both wind and current
    def __init__(self, fn:str, logger:Logger, outdir:str) -> None:
        BaseParam.__init__(self, fn, logger, outdir)
        if self.info is None: return
        for item in self.info['speed']:
            item['latMin'] = min(item['latLimits'])
            item['latMax'] = max(item['latLimits'])
            item['lonMin'] = min(item['lonLimits'])
            item['lonMax'] = max(item['lonLimits'])

    def drift(self, lat:float, lon:float, dt:float) -> tuple:
        """ Drift distance in meters for x and y direction """
        if (dt <= 0) or (self.info is None): return (0, 0)

        for item in self.info['speed']:
            if (item['latMin'] <= lat) and \
                    (lat <= item['latMax']) and \
                    (item['lonMin'] <= lon) and \
                    (lon <= item['lonMax']):
                return self.__calcDrift(item, lat, lon, dt)
        return (0, 0)

    def __calcDrift(self, item:dict, lat:float, lon:float, dt:float) -> tuple:
        dLat = lat - item["latMin"]
        dLon = lon - item["lonMin"]
        norm = self.info['normalization']
        spd = item["speedLL"] + dLat * item["latSlope"] + dLon * item["lonSlope"]
        distance = norm * spd * dt
        hdg = math.radians(item["direction"])
        self.logger.debug('Drift spd=%s distance=%s hdg=%s',
                round(spd,2), round(distance,1), math.degrees(hdg))
        return (distance * math.sin(hdg), distance * math.cos(hdg))
