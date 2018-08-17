#!/usr/bin/env python
import json
import argparse
import numpy
import sys
import copy
from astropy.coordinates import SkyCoord
from astropy import units
import logging

# These are the exposure times set in PH2 on CFHT (or must be) so that we get the correct ones.
# First exposure time is I1, then I2 etc.
# Besure these exposure times match whats in the phase 2.
IC_exptimes = [40, 80, 120, 160, 200, 240, 300, 340, 380, 420, 440, 480, 500]

def exposure_time_index(mag):
    """
    Compute the exposure time required, in seconds, scaling off a 300s exposure needed for a 24.5 mag source.

    Mimimum exposure time is 40s (CFHT Overhead)
    :return: float
    """
    cuts = numpy.array(IC_exptimes)
    exact_exptime = min(499, max(40, 300.0 / ((10 ** ((24.5 - mag) / 2.5)) ** 2)))
    return len(cuts) - (exact_exptime < cuts).sum()


def exposure_time(mag):
    return IC_exptimes[exposure_time_index(mag)]


def instrument_configuration_identifier(mag):
    return "I{}".format(exposure_time_index(mag)+1)


class Program(object):
    def __init__(self, runid="17BC08", pi_login="kavelaars"):
        self.config = {"runid": runid,
                       "pi_login": pi_login,
                       "program_configuration": {"targets": [],
                                                 "observing_blocks": [],
                                                 "observing_groups": []
                                                 }}

    def add_target(self, target):
        self.config["program_configuration"]["targets"].append(target)

    def add_observing_block(self, observing_block):
        self.config["program_configuration"]["observing_blocks"].append(observing_block)

    def add_observing_group(self, observing_group):
        self.config["program_configuration"]["observing_groups"].append(observing_group)


class Target(object):
    def __init__(self, filename=None):
        self.config = json.load(open(filename))

    @property
    def name(self):
        return self.config["name"]

    @property
    def token(self):
        return self.config["identifier"]["client_token"]

    @property
    def mag(self):
        return self.config["moving_target"]["ephemeris_points"][0]["mag"]

    @property
    def coordinate(self):
        return SkyCoord(self.config["moving_target"]["ephemeris_points"][0]["coordinate"]["ra"],
                        self.config["moving_target"]["ephemeris_points"][0]["coordinate"]["dec"],
                        unit='degree')


class ObservingBlock(object):
    def __init__(self, client_token, target_token):
        self.config = {"identifier": {"client_token": client_token},
                       "target_identifier": {"client_token": target_token},
                       "constraint_identifiers": [{"server_token": "C1"}],
                       "instrument_config_identifiers": [{"server_token": "I1"}]}

    @property
    def token(self):
        return self.config["identifier"]["client_token"]


class ObservingGroup(object):
    def __init__(self, client_token):
        self.config = {"identifier": {"client_token": client_token},
                       "observing_block_identifiers": []}

    def add_ob(self, client_token):
        self.config["observing_block_identifiers"].append({"client_token": client_token})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('runid')
    parser.add_argument('qrunid')
    parser.add_argument('targets', nargs='+')
    parser.add_argument('--verbose', help="Verbose message reporting.", action="store_true", default=False)
    parser.add_argument('--debug', help="Provide debuging information.", action="store_true", default=False)
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    logging.basicConfig(level=logging.ERROR)

    program = Program(args.runid, pi_login="kavelaars")
    ob_tokens = []
    mags = {}
    ob_coordinate = {}
    for filename in args.targets:
        target = Target(filename)
        logging.info("Programming: {} V:{}".format(target.name, target.mag))
        try:
            mag = target.mag
        except:
            mag = 25.0
        program.add_target(target.config)
        ob_token = "OB-{}-{}".format(args.qrunid, target.token)
        ob = ObservingBlock(ob_token, target.token)
        idx = instrument_configuration_identifier(mag)
        ob.config["instrument_config_identifiers"] = [{"server_token": "I{}".format(idx)}]
        program.add_observing_block(ob.config)
        ob_tokens.append(ob_token)
        mags[ob_token] = mag
        ob_coordinate[ob_token] = target.coordinate

    # Order the tokens by RA of the target.
    sf = lambda x, y: cmp(x.ra, y.ra)
    order_tokens = sorted(ob_coordinate, cmp=sf, key=ob_coordinate.get)

    # Keep track of total observing time, which OGs have been built and which OBs are in an OG.
    total_itime = 0
    ogs = {}
    scheduled = {}

    # start the OG index at 1.
    og_idx = 0
    while len(scheduled) < len(ob_tokens):
        og_idx += 1
        og_coord = None
        og_itime = 0
        repeat = 0
        og_token = "OG_{}_{}_{}_{}".format(args.runid, args.qrunid, og_idx, repeat)
        og = ObservingGroup(og_token)

        for ob_token in order_tokens:
            if ob_token not in scheduled:
                if og_coord is None:
                    og_coord = ob_coordinate[ob_token]
                if ob_coordinate[ob_token].separation(og_coord) > 10 * units.degree:
                    # Skip over targets that are more that 10 degrees away.
                    # They get their own OG.
                    continue
                og.add_ob(ob_token)
                scheduled[ob_token] = True
                og_itime += exposure_time(mags[ob_token]) + 40
                if og_itime > 3000.0:
                    break
                break

        total_itime += og_itime
        sys.stdout.write("OG {} is {}s in duration.\n".format(og_token, og_itime))
        program.add_observing_group(og.config)
        nrepeats = 2  # do each OG twice, for tracking.
        for repeat in range(nrepeats):
            total_itime += og_itime
            og_token = "OG_{}_{}_{}_{}".format(args.runid, args.qrunid, og_idx, repeat+1)
            og = copy.deepcopy(og)
            og.config["identifier"]["client_token"] = og_token
            program.add_observing_group(og.config)

    json.dump(program.config, open('PH2_{}_{}.json'.format(args.runid, args.qrunid), 'w'), indent=4, sort_keys=True)
