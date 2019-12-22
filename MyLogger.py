#
# Construct a logger given a standard set of command line arguments
#
# Oct-2019, Pat Welch, pat@mousebrains.com

import logging
import logging.handlers
import argparse

def addArgs(parser: argparse.ArgumentParser):
    """ Add my options to an argparse object """
    grp = parser.add_argument_group('Log related options')
    grp.add_argument('--logfile', type=str, help='Log filename')
    grp.add_argument('--verbose', action='store_true', help='Enable debug messages')
    grp.add_argument('--maxlogsize', type=int, default=10000000, help='Maximum logfile size')
    grp.add_argument('--backupcount', type=int, default=7, help='Number of logfiles to keep')

def mkLogger(args:argparse.ArgumentParser, name:str, 
        fmt:str=None) -> logging.Logger:
    """ Make a logging object based on the options in args """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    if args.logfile is None:
        ch = logging.StreamHandler()
    else:
        ch = logging.handlers.RotatingFileHandler(args.logfile,
                                maxBytes=args.maxlogsize, 
                                backupCount=args.backupcount)

    ch.setLevel(logging.DEBUG if args.verbose else logging.INFO)

    if fmt is None:
        fmt = '%(asctime)s: %(threadName)s:%(levelname)s - %(message)s'
    formatter = logging.Formatter(fmt)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger
