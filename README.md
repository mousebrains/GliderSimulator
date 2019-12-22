# GliderSimulator

A detailed simulator for an AUV, including energy consumption.

Specifically this is designed with a TWR Slocum glider in mind.

The simulator is controlled by a set of JSON files, whose names are specified on the command line:
 1) glider.json is the main parameter file for how the glider will operate
 2) waypoints.json is a list of waypoints to fly, once the end of the list is hit, it wraps to the front.
 3) energy.json defines the energy usage for various phased of the flight
 4) current.json defines the depth averaged water currents the glider encounters
 5) wind.json defines the wind the glider encounters on the surface

The simulation is pretty detailed:
 1) Surface with the following included:
    1) Inflate the air bladder
    1) Get an initial GPS fix
    1) Iridium call setup
    1) Iridium data transfer (time dependent on previous dive duration)
    1) Iridium call take down
    1) Get final GPS fix
    1) Deflate the air bladder
    1) Adjust the position due to the winds
 2) Dive, normal climb, and climb to surface
    1) Dive speed specified in JSON file
    1) Horizontal speed is determined by Wg / tan(pitch + angleOfAttack)
    1) The bottom depth comes from bathymetry information. I use GEBCO.
    1) Dive time and distance are calculated.
    1) The dive time includes the inflection upwards at the bottom.
    1) A check is run to see if a waypoint was passed. If so then the glider advances to the next waypoint.
    1) A new location is estimated based on distance, heading, and currents.
