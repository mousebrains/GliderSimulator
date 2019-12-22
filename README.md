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
 #) Surface with the following included:
   .) Inflate the air bladder
   .) Get an initial GPS fix
   .) Iridium call setup
   .) Iridium data transfer (time dependent on previous dive duration)
   .) Iridium call take down
   .) Get final GPS fix
   .) Deflate the air bladder
   .) Adjust the position due to the winds
 #) Dive, normal climb, and climb to surface
   .) Dive speed specified in JSON file
   .) Horizontal speed is determined by Wg / tan(pitch + angleOfAttack)
   .) The bottom depth comes from bathymetry information. I use GEBCO.
   .) Dive time and distance are calculated.
   .) The dive time includes the inflection upwards at the bottom.
   .) A check is run to see if a waypoint was passed. If so then the glider advances to the next waypoint.
   .) A new location is estimated based on distance, heading, and currents.
