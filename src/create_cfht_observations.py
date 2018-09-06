#!/usr/bin/env python

import argparse
import json
import logging
import sys

import numpy
from astropy import units
from astropy.coordinates import SkyCoord
from cfht.sdk.http_client import HttpClient
from cfht.sdk.repository_strategy import ApiStrategy, PrintStrategy
from cfhtapi.observing.api_pb2 import PiObservingGroup as ObservingGroup
from cfhtapi.observing.persistence_pb2 import TargetData as Target
from google.protobuf import json_format

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
    return "I{}".format(exposure_time_index(mag) + 1)


def sky_coord(data_coord):
    return SkyCoord(data_coord.ra, data_coord.dec)

class CfhtUpload():
    def __init__(self, runid, qrunid, targets, bearer_token, host, mode):

        if mode == 'print':
            self.cfht_strategy = PrintStrategy(sys.stdout)
        elif mode == 'api':
            self.cfht_strategy = ApiStrategy(HttpClient(host, bearer_token))
        else:
            raise Exception("Unknown mode %s, must be one of (console, api)" % args.mode)

        self.observing_group_repository = self.cfht_strategy.observing_group_repository(runid)
        self.observing_template_repository = self.cfht_strategy.observing_template_repository(runid)
        self.target_repository = self.cfht_strategy.observing_template_repository(runid)
        self.generate_target_token = self.cfht_strategy.target_repository(runid).generate_token
        self.generate_observing_group_token = self.observing_group_repository.generate_token
        self.targets = targets
        self.qrunid = qrunid

    def execute(self):

        # Collect available observing templates by exposure time
        observing_template_by_exptime = {}
        for observing_template in self.observing_template_repository.list():
            etime_sec = int(observing_template.instrument_configuration.etime_or_snr.exposure_time_ms / 1000.0)
            observing_template_by_exptime[etime_sec] = observing_template

        # Create a list of targets to persist, and the target -> observing template mappings to persist into ogs
        targets_to_create_or_update = []
        for filename in self.targets:
            target = json_format.Parse(json.load(open(filename)), Target())
            target.token = self.generate_target_token(target.name)
            etime_sec = exposure_time(target.magnitude.magnitude)
            if etime_sec not in observing_template_by_exptime:
                raise Exception(
                        "Could not find observing template for etime %s, you must configure an IC in ph2" % etime_sec)

            targets_to_create_or_update.append(target)

        for target in targets_to_create_or_update:
            logging.info("Programming: {} V:{}".format(target.name, target.magnitude.magnitude))
            self.target_repository.persist(target.token, target)

        # Order the tokens by RA of the target.
        sf = lambda x, y: cmp(x[0].moving_target.ephemeris_point[0].coordinate.ra,
                              y[0].moving_target.ephemeris_point[0].coordinate.ra)
        targets_by_ra = sorted(targets_to_create_or_update, cmp=sf)

        # Keep track of total observing time, which OGs have been built and which OBs are in an OG.
        total_itime = 0
        ogs = {}
        scheduled = {}

        # start the OG index at 1.
        og_idx = 0
        while len(scheduled) < len(targets_by_ra):
            og_idx += 1
            og_coord = None
            og_itime = 0

            observing_group = ObservingGroup()
            observing_block = observing_group.read_observing_block.add()

            for target in targets_by_ra:
                if target.token not in scheduled:
                    target_coord = sky_coord(target.moving_target.ephemeris_point[0].coordinate)
                    # First time through loop
                    if og_coord is None:
                        og_coord = target_coord
                    if target_coord.separation(og_coord) > 10 * units.degree:
                        # Skip over targets that are more that 10 degrees away.
                        # They get their own OG.
                        continue

                    scheduled[target.token] = True
                    etime_sec = exposure_time(target.magnitude.magnitude)#  + 40 (Why + 40... probably not needed now)

                    og_itime += etime_sec
                    observing_component = observing_block.observing_component.add()
                    observing_component.target_token = target.token
                    observing_component.observing_teplate_token = observing_template_by_exptime[etime_sec].token
                    if og_itime > 3000.0:
                        break
                    break

            for repeat in range(3): # do each OG twice, for tracking. ... this is now 3 to mimic previous behavior...
                og_token = self.generate_observing_group_token("%s-%s-%s" % (self.qrunid, og_idx, repeat))
                total_itime += og_itime
                sys.stdout.write("OG {} is {}s in duration.\n".format(og_token, og_itime))
                self.observing_group_repository.persist(og_token, observing_group)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('runid')
    parser.add_argument('qrunid',
                        help="guarantees a new OG is created for the upcoming run, even if it's been observed")
    parser.add_argument('targets', nargs='+')
    parser.add_argument('--bearer-token', help="Authorization Token")
    parser.add_argument('--verbose', help="Verbose message reporting.", action="store_true", default=False)
    parser.add_argument('--debug', help="Provide debuging information.", action="store_true", default=False)
    parser.add_argument('--host', help="URL of the CFHT api", default="https://api.cfht.hawaii.edu")
    parser.add_argument('--mode', help="(print, api) indicates if the output should go to the console, or the api",
                        default="api")

    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    logging.basicConfig(level=logging.WARN)

    CfhtUpload(args.runid, args.qrunid, args.targets, args.bearer_token, args.host, args.mode).execute()
