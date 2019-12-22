#
# Calculate distance, azimuthal angle, and dead reckoning
# on the earth. 
# All lat/lon are in decimal degrees
# All angles/headings are in degrees true.
# All distances are in meters
#

from math import sin, cos, degrees, radians, atan2, asin, sqrt

radiusEarth = 6371.0088e3 # in meters

def distance(lat0:float, lon0:float, lat1:float, lon1:float) -> float:
    """ Calculate the distance in meters between (lat0,lon0) and (lat1,lon1) 
        using the haversine formula
    """
    lat0 = radians(lat0)
    lat1 = radians(lat1)
    lon0 = radians(lon0)
    lon1 = radians(lon1)
    dLat = lat1 - lat0
    dLon = lon1 - lon0
    d = sin(dLat/2)**2  + cos(lat0) * cos(lat1) * sin(dLon/2)**2
    return 2 * radiusEarth * asin(sqrt(d))

def bearing(lat0:float, lon0:float, lat1:float, lon1:float) -> float:
    """ Initial bearing to go from pt0 to pt1 in degrees true """
    lat0 = radians(lat0)
    lat1 = radians(lat1)
    lon0 = radians(lon0)
    lon1 = radians(lon1)
    dLat = lat1 - lat0
    dLon = lon1 - lon0
    x = sin(dLon) * cos(lat1)
    y = cos(lat0) * sin(lat1) - sin(lat0) * cos(lat1) * cos(dLon)
    return degrees(atan2(x, y))

def dreckon(lat:float, lon:float, dist:float, bearing:float) -> tuple:
    """ Shift lat/lon by a distance and bearing """
    hdg = radians(bearing)
    return drift(lat, lon, dist * sin(hdg), dist * cos(hdg))

def drift(lat:float, lon:float, dx:float, dy:float) -> tuple:
    """ Shift lat/lon by distances dx/dy """
    mPerLat = distance(lat-0.5, lon, lat+0.5, lon)
    mPerLon = distance(lat, lon-0.5, lat, lon+0.5)
    return (lat + dy / mPerLat, lon + dx / mPerLon) 
