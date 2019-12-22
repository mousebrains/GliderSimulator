# GliderSimulator

A detailed simulator for an AUV, including energy consumption.

Specifically this is designed with a TWR Slocum glider in mind.

The simulator is controlled by a set of JSON files, whose names are specified on the command line:
 1) glider.json is the main parameter file for how the glider will operate
 2) waypoints.json is a list of waypoints to fly, once the end of the list is hit, it wraps to the front.
 3) energy.json defines the energy usage for various phased of the flight
 4) current.json defines the depth averaged water currents the glider encounters
 5) wind.json defines the wind the glider encounters on the surface
