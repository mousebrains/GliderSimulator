#! /usr/bin/env python3
#
# Simulate a glider's flight path and energy usage.
#
# The configuration is done through a set of JSON parameter files:
#
# --glider is a JSON file with parameters for the Glider
# --current is a JSON file with piecewise linear equations for the current
# --energy is a JSON file defining how much energy is used by various actions
# --waypoints is a JSON file with a sequence of waypoints
# --mission is a JSON file with initial parameters for the mission
#
# Dec-2019, Pat Welch, pat@mousebrains.com

import argparse
import MyLogger
import Sim

parser = argparse.ArgumentParser(description='Glider Simulator')
parser.add_argument('--log', type=str, help='Filename to write log records to')
MyLogger.addArgs(parser)
Sim.addArgs(parser)
args = parser.parse_args()

logger = MyLogger.mkLogger(args, 'Sim')

logger.info('args=%s', args)

Sim.runSimulation(args, logger)

logger.info('Finished with simulation')
