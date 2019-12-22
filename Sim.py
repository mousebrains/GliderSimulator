#
# Main part of the simulator
#
# Dec-2019, Pat Welch, pat@mousebrains.com

import DEM
import MyGeo
import math
from logging import Logger
from argparse import ArgumentParser
import SimParams
import os
import datetime

def addArgs(parser:ArgumentParser):
    """ Add my options to an argparse object """

    parser.add_argument('--dem', type=str, required=True, help='DEM filename in ESRI format')
    parser.add_argument('--outdir', type=str, help='Directory to write results to')

    grp = parser.add_argument_group('JSON related options')
    grp.add_argument('--glider', type=str, required=True, help='Glider config JSON')
    grp.add_argument('--energy', type=str, help='Glider energy config JSON')
    grp.add_argument('--waypoints', type=str, required=True, help='Waypoints to fly JSON')
    grp.add_argument('--mission', type=str, required=True, help='Define mission to fly JSON')
    grp.add_argument('--current', type=str, help='Define current to fly through JSON')

    grp.add_argument('--wind', type=str, help='Define surface wind to operate in JSON')
    grp = parser.add_argument_group('Simulation related options')
    grp.add_argument('--maxDuration', type=float, default=1000,
            help='Maximum deployment duration in days')
    grp.add_argument('--maxYos', type=int, default=10000, help='Maximum number of yos')
    grp.add_argument('--nYos', type=int, help='Number of yos between surfacings')

def runSimulation(args:ArgumentParser, logger:Logger) -> None:
    sim = Simulate(args, logger)

    while sim: # Keep on running
        sim.step() # Step to the next state

    sim.close() # Clean up

class Simulate: # A simulation
    def __init__(self, args:ArgumentParser, logger:Logger):
        self.args = args
        self.logger = logger
        self.dem = DEM.DEM(args.dem)
        self.__initialize(args)
        self.__dataFP = None

    def __bool__(self) -> bool:
        """ Can the simulation continue or not, i.e. time and energy left """
        return bool(self.t < self.tMax) \
                and (self.state != 'done') \
                and bool(self.energy) \
                and (self.totalYos < self.maxYos)

    def __initialize(self, args:ArgumentParser) -> None: 
        """ Set up the initial data structures """
        if (args.outdir is not None) and not os.path.isdir(args.outdir): 
            self.logger.info('Creating %s for the output directory', args.outdir)
            os.makedirs(args.outdir)
        self.glider = SimParams.Glider(args.glider, self.logger, args.outdir)
        self.energy = SimParams.Energy(args.energy, self.logger, args.outdir)
        self.waypoints = SimParams.Waypoints(args.waypoints, self.logger, args.outdir)
        self.mission = SimParams.Mission(args.mission, self.logger, args.outdir)
        self.current = SimParams.Drift(args.current, self.logger, args.outdir)
        self.wind = SimParams.Drift(args.wind, self.logger, args.outdir)
        self.tMax = args.maxDuration * 86400 # Maximum deployment duration length
        self.state = 'surface' # Next state to be processed
        self.t = 0 # Initial mission time in seconds
        self.tLastSurface = 0 # Time of the previous surfacing
        self.tNextSurfacing = self.t + 900 # This is not actually used
        self.nYosBetweenSurfacings = self.glider.yosBetweenSurfacings() \
                if args.nYos is None else args.nYos
        self.volume = self.glider.surfaceVolume()
        self.waypointRadius = self.glider.waypointApproachRadius()
        self.qSurfaceOnWaypoint = self.glider.waypointSurfaceOnHitting()
        self.depth = 0 # Initially on the surface
        self.totalYos = 0 # How many yos have been done in total
        self.nYos = 0 # How many yos since the last surfacing
        self.maxYos = args.maxYos
        wpt = self.waypoints.initial()
        self.wpt = self.waypoints.next()
        self.lat = wpt['lat'] if wpt is not None else 0
        self.lon = wpt['lon'] if wpt is not None else 0
        self.qUseCurrentCorrection = self.glider.currentUseCorrection();
        self.currentMaxPerpFrac = self.glider.currentMaxFraction()
        self.currentToFastFrac = self.glider.currentToFastFraction()
        self.avgHorizontalSpeed = self.glider.average_horizontal_speed()
        self.latGlider = self.lat # Where the glider thinks it is
        self.lonGlider = self.lon
        self.water_vx = 0 # Depth averaged water velocity in x direction (m/s)
        self.water_vy = 0 # Depth averaged water velocity in y direction (m/s)
        self.heading = 0 # Heading the glider is on to the next waypoint (deg)

    def close(self) -> None:
        """ Close the open file handle """
        if self.__dataFP is not None: self.__dataFP.close()

    def __dumpState(self) -> None:
        """ Dump out the current state to a file """
        if self.__dataFP is None: # Not opened yet
            if self.args.outdir is None: return # Nowhere to write to
            ofp = os.path.join(self.args.outdir, 'sim.csv')
            fp = open(ofp, 'w')
            self.__dataFP = fp
            fp.write('state' \
                    + ',time' \
                    + ',lat,lon' \
                    + ',m_lat,m_lon' \
                    + ',latWpt,lonWpt' \
                    + ',water_vx,water_vy' \
                    + ',depth' \
                    + ',volume' \
                    + ',heading' \
                    + ',energy' \
                    + '\n')

        fp = self.__dataFP
        fp.write(self.state \
                + ',' + str(self.t) \
                + ',' + str(self.lat) + ',' + str(self.lon) \
                + ',' + str(self.latGlider) + ',' + str(self.lonGlider) \
                + ',' + str(self.wpt['lat']) + ',' + str(self.wpt['lon']) \
                + ',' + str(self.water_vx) + ',' + str(self.water_vy) \
                + ',' + str(self.depth) \
                + ',' + str(self.volume) \
                + ',' + str(self.heading) \
                + ',' + str(self.energy) \
                + "\n")

    def __time2string(self) -> str:
        t = self.t
        days = math.floor(t / 86400)
        hours = math.floor((t % 86400) / 3600)
        minutes = math.floor((t % 3600) / 60)
        seconds = math.floor(t % 60)
        if days > 0: return '{} {:02d}:{:02d}:{:02d}'.format(days, hours, minutes, seconds)
        return '{:02d}:{:02d}:{:02d}'.format(hours, minutes, seconds)

    def step(self) -> None: # Walk to the next step
        self.logger.info('STEP state=%s t=%s energy=%s volume=%s depth=%s lat=%s lon=%s', 
                self.state,
                self.__time2string(),
                self.energy, 
                self.volume, round(self.depth, 1),
                round(self.lat,6), round(self.lon,6))
        if self.state == 'surface':
            self.__stepSurface()
        elif self.state == 'dive':
            self.__stepDive()
        elif self.state == 'climb':
            self.__stepClimb()
        elif self.state == 'climbSurface':
            self.__stepClimbSurface()
        else:
            raise Exception('Unrecognized state, "' + self.state + '"')

    def __stepSurface(self) -> None:
        """ I'm at the surface """
        self.__calcWaterVelocities()

        (dtSurf, dtCall) = self.glider.surfaceTime(self.tLastSurface)
        self.t += dtSurf
        self.tLastSurface = self.t # For an estimate of data and water velocities
        self.tNextSurfacing = self.glider.nextSurfacingTime(self.t, self.tNextSurfacing)
        self.nYos = 0 # Reset number of yos since the last surfacing
        dE = self.energy.surface(dtSurf, dtCall)
        self.energy -= dE

        # Drift due to the wind
        (xDrift, yDrift) = self.wind.drift(self.lat, self.lon, dtSurf)
        (self.lat, self.lon) = MyGeo.drift(self.lat, self.lon, xDrift, yDrift) 

        self.latGlider = self.lat # Rest where the glider thinks it is
        self.lonGlider = self.lon # Rest where the glider thinks it is

        self.__nextWaypoint()

        self.__dumpState() # Write out the state to a CSV file

        self.state = 'dive'
        self.logger.info('SURFACE t=%s tNxt=%s dt=%s dtCall=%s dE=%s drift %s,%s', 
                round(self.t,1), round(self.tNextSurfacing,1), 
                round(dtSurf,1), round(dtCall,1), round(dE,1),
                round(xDrift, 1), round(yDrift))


    def __stepDive(self) -> None:
        """ Dive to the target depth """
        nxtState = 'climb'
        self.totalYos += 1
        self.nYos += 1
        depth = self.dem.depth(self.lat, self.lon)
        dtMax = max(0, self.tNextSurfacing - self.t) # How long until next surfacing triggered
        (dt, dist, volume, eDepth) = self.glider.diveTime(self.depth, depth, dtMax)

        # Check if I am going to hit a waypoint
        distCheck = self.__checkForWaypoint(dt, dist)
        if distCheck is not None:
            dtMax = dt * (distCheck / dist) # How much time for this action
            (dt, dist, volume, eDepth) = self.glider.diveTime(self.depth, dtMax)
            nxtState = self.state

        self.__adjustPosition(dt, dist) # Change the lat/lon by the distance traveled

        self.logger.info('DIVE t=%s dt=%s dist=%s eDepth=%s',
                round(self.t,1), round(dt,1), round(dist,1), round(eDepth,1))

        self.__adjustState(dt, dist, volume, eDepth, 
                self.energy.dive(dt, abs(self.volume - volume)))

        if (self.t + dt) >= self.tNextSurfacing:
            self.logger.debug('Aborting dive depth=%s eDepth=%s', self.depth, eDepth)
            nxtState = 'climbSurface'

        if distCheck is not None: # Advance the waypoint
            self.wpt = self.waypoints.next()
            self.__nextWaypoint()

        self.state = nxtState

    def __stepClimb(self) -> None:
        """ Climb to the target depth """
        nxtState = 'dive'

        dtMax = max(0, self.tNextSurfacing - self.t) # How long until next surfacing triggered
        (dt, dist, volume, eDepth) = self.glider.climbTime(self.depth, dtMax)

        # Check if I am going to hit a waypoint
        distCheck = self.__checkForWaypoint(dt, dist)
        if distCheck is not None:
            dtMax = dt * (distCheck / dist) # How much time for this action
            (dt, dist, volume, eDepth) = self.glider.climbTime(self.depth, dtMax)
            nxtState = self.state

        self.__adjustPosition(dt, dist) # Change the lat/lon by the distance traveled

        self.logger.info('CLIMB t=%s dt=%s dDepth=%s dist=%s n=%s',
                round(self.t,1), round(dt,1), 
                round(self.depth - eDepth, 1),
                round(dist,1), self.nYos)

        self.__adjustState(dt, dist, volume, eDepth,
                self.energy.climb(dt, abs(self.volume - volume)))

        if ((self.t + dt) >= self.tNextSurfacing) or (self.nYos >= self.nYosBetweenSurfacings):
            self.logger.debug('Aborting climb depth=%s eDepth=%s', self.depth, eDepth)

        if distCheck is not None: # Advance the waypoint
            self.wpt = self.waypoints.next()
            self.__nextWaypoint()

        self.state = nxtState

    def __stepClimbSurface(self) -> None:
        """ Climb to the surface """
        nxtState = 'surface'

        (dt, dist, volume, eDepth) = self.glider.climbSurfaceTime(self.depth, 1e90)

        # Check if I am going to hit a waypoint
        distCheck = self.__checkForWaypoint(dt, dist)
        if distCheck is not None:
            dtMax = dt * (distCheck / dist) # How much time for this action
            (dt, dist, volume, eDepth) = self.glider.climbSurfaceTime(self.depth, dtMax)
            nxtState = self.state

        self.__adjustPosition(dt, dist) # Change the lat/lon by the distance traveled

        self.logger.info('CLIMB2SURFACE dt=%s dDepth=%s dist=%s',
                round(dt,1), round(eDepth,1), round(dist,1))

        self.__adjustState(dt, dist, volume, eDepth, 
                self.energy.climb(dt, abs(self.volume - volume)))

        if distCheck is not None: # Advance the waypoint
            self.wpt = self.waypoints.next()
            self.__nextWaypoint()
            nxtState = self.state

        self.state = nxtState

    def __adjustState(self, dt:float, dist:float, volume:float, eDepth:float, dE:float) -> None:
        self.logger.debug('Adj %s dt=%s dist=%s vol=%s eDepth=%s dE=%s', 
                self.state, round(dt, 1), round(dist, 1), volume, round(eDepth, 1), round(dE))
        self.t += dt
        self.volume = volume
        self.depth = eDepth
        self.energy -= dE
        self.__dumpState() # Write out the state to a CSV file

    def __nextWaypoint(self) -> None:
        """ If current correction is on, then adjust the aiming point """
        self.latAim = self.wpt['lat']
        self.lonAim = self.wpt['lon']
        hdg = math.radians(MyGeo.bearing(self.lat, self.lon, self.latAim, self.lonAim))
        q = self.qUseCurrentCorrection # No current correction
        wx = self.water_vx if q else 0
        wy = self.water_vy if q else 0
        U = self.avgHorizontalSpeed
        Ux = U * math.sin(hdg)
        Uy = U * math.cos(hdg)
        # Take the cross product of (wx,wy) and (Ux,Uy)/U to get the perpendicular
        # component of current to our desired direction of travel
        wPerp = abs(wx * Uy - wy * Ux) / U
        if wPerp >= (U * self.currentMaxPerpFrac): # invalid solution
            alpha = self.currentToFastFrac * U / wPerp
            wx *= alpha
            wy *= alpha
        hdg = math.atan2(Ux - wx, Uy - wy)
        self.heading = math.degrees(hdg)
        self.cheading = math.cos(hdg)
        self.sheading = math.sin(hdg)
        dist = MyGeo.distance(self.lat, self.lon, self.latAim, self.lonAim)
        (self.latAim, self.lonAim) = MyGeo.dreckon(self.lat, self.lon, dist, self.heading)
        self.logger.info('Next wpt adj %s %s hdg %s', 
                self.latAim, self.lonAim, round(self.heading,3))

    def __checkForWaypoint(self, dt:float, dist:float) -> float:
        """ Check if we would pass a waypoint close enough to trigger advancing to the next """
        hdg = math.radians(self.heading) # What direction I'm heading towards
        dx = dist * math.sin(hdg)
        dy = dist * math.cos(hdg)
        if self.qUseCurrentCorrection:
            dx += dt * self.water_vx
            dy += dt * self.water_vy
        lat0 = self.latGlider # Where the glider thought it was at the beginning
        lon0 = self.lonGlider
        (lat1, lon1) = MyGeo.drift(lat0, lon0, dx, dy)

        latWpt = self.wpt['lat']
        lonWpt = self.wpt['lon']

        (latC, lonC, distC) = self.__closestApproach(lat0, lon0, lat1, lon1, latWpt, lonWpt)

        # self.logger.info('CLS wpt %s %s C %s %s %s',
        #         latWpt, lonWpt, round(latC, 6), round(lonC, 6), round(distC, 1))

        if distC >= self.waypointRadius: return None # Still too far

        d = MyGeo.distance(lat0, lon0, latC, lonC) # Distance to point of closest approach
        self.logger.info('CHK distance to point of closest approach %s from %s', d, dist)

        return d


    def __distanceToEnds(self,
            lat0:float, lon0:float, 
            lat1:float, lon1:float, 
            latWpt:float, lonWpt:float) -> tuple:
        """ Return the endpoint that is closest to waypoint """
        dist0 = MyGeo.distance(lat0, lon0, latWpt, lonWpt)
        dist1 = MyGeo.distance(lat1, lon1, latWpt, lonWpt)
        return (lat0, lon0, dist0) if dist0 <= dist1 else (lat1, lon1, dist1) 

    def __closestApproach(self, 
            lat0:float, lon0:float, 
            lat1:float, lon1:float, 
            latWpt:float, lonWpt:float) -> tuple:
        """ Calculate a trajectory's closest approach to a waypoint 
            and return the distance and point. 

            At the point of closest approach the dot product will be zero.

        Parameterize lat/lon0 -> lat/lon1 as:
            let dx = x1-x0 and dy=y1-y0
            y(x) = y0 + (x-x0) dy/dx
                 = (y0 - x0 dy/dx) + x dy/dx
                 = affine + x slope

            lon is x and lat is y
        """

        dy = lat1 - lat0 
        dx = lon1 - lon0

        if dx == 0: # Vertical trajectory
            if (latWpt < min(lat0, lat1)) or (latWpt > max(lat0, lat1)):
                return self.__distanceToEnds(lat0, lon0, lat1, lon0, latWpt, lonWpt)
            dist = MyGeo.distance(latWpt, lon0, latWpt, lonWpt)
            return (latWpt, lon0, dist)

        b = dy / dx # Slope
        a = lat0 - lon0 * b # Affine

        # The trajectory vector is then
        # (x0, y0) -> (x1, y1)
        # 
        # v0 = (x1-x0, y1-y0)
        #
        # The perpendicular vector from (x,y(x)) to the waypoint is
        # v1 = (xwpt-x, ywpt-y(x))
        #
        # The dot product is 
        # 0 = v0 dot v1
        # 0 = dx (xwpt-x) + dy (ywpt-y(x))
        #   = dx xwpt - x dx + dy ywpt - dy (a + b x)
        #   = xwpt dx + ywpt dy - a dy - x (dx + b dy)
        # x (dx + b dy) = xwpt dx + ywpt dy - a dy
        # x = (xwpt dx + ywpt dy - a dy) / (dx + b dy)
        #
        # N.B. dx!=0 and b = dy/dx => dx + b dy = dx + dy^2 / dx = (dx^2 + dy^2) / dx != 0

        xWpt = lonWpt
        yWpt = latWpt

        x = (xWpt * dx + yWpt * dy - a * dy) / (dx + b * dy)

        if (x < lon0) or (x > lon1): # outside trajectory
            return self.__distanceToEnds(lat0, lon0, lat1, lon1, latWpt, lonWpt)
        
        y = a + b * x
        dist = MyGeo.distance(y, x, latWpt, lonWpt)
        self.logger.info('Inside x %s y %s dist %s', y, x, dist)
        return (y, x, dist)

    def __adjustPosition(self, dt:float, dist:float) -> None:
        """ Calculate a new lat/lon in reality and for what the glider believes """
        x = dist * self.sheading
        y = dist * self.cheading
        
        (xDrift, yDrift) = self.current.drift(self.lat, self.lon, dt)
        (self.lat, self.lon) = MyGeo.drift(self.lat, self.lon, x + xDrift, y + yDrift) 

        # Now where does the glider think it is?

        q = self.qUseCurrentCorrection
        dx = (self.water_vx * dt) if q else 0
        dy = (self.water_vy * dt) if q else 0

        (latGlider, lonGlider) = MyGeo.drift(self.latGlider, self.lonGlider, x + dx, y + dy)
        self.__closestApproach(self.latGlider, self.lonGlider, latGlider, lonGlider, 
                self.wpt["lat"], self.wpt["lon"])
        self.latGlider = latGlider
        self.lonGlider = lonGlider

    def __calcWaterVelocities(self) -> None:
        """ Look at the difference in lat/lon in reality and for what the glider thinks """
        lat1 = self.lat
        lon1 = self.lon
        lat0 = self.latGlider
        lon0 = self.lonGlider

        dist = MyGeo.distance(lat0, lon0, lat1, lon1)
        hdg = math.radians(MyGeo.bearing(lat0, lon0, lat1, lon1))
        dx = dist * math.sin(hdg)
        dy = dist * math.cos(hdg)
        dt = self.t - self.tLastSurface
        self.water_vx = (dx/dt) if dt > 0 else 0
        self.water_vy = (dy/dt) if dt > 0 else 0
