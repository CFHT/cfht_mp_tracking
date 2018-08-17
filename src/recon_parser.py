"""
Parse the TNORecon tables to determine which sources should be scheduled at CFHT.
"""
import argparse
import sys
from bs4 import BeautifulSoup
import requests
import pandas as pd
from datetime import date
from astropy.time import Time
from astropy import units
from astropy.table import Table
from datetime import datetime
from cStringIO import StringIO
import numpy
import ephem
import math
import logging
from mp_ephem import horizons

DESCRIPTION = """Connects to the web server at SWRI to retrieve various lists of occultation and apulse predictions.
Parses through the table on those pages to deliver a list of targets that would be suitable for tracking with CFHT
during the given time period.
"""

MINIMUM_UP_TIME = 1.0 * units.hour
DES_CLASSES = ['CENTAURR',
               'ERR2LARGE',
               'RESONANT',
               'CLASSICAL',
               'SCATNEAR',
               ]
SERVICE_URL = "http://www.boulder.swri.edu/~buie/recon/"
DEFAULT_ORB_CLASS = 'RESONANT'
MINIMUM_ELEVATION = 40*units.degree
# Column Names in RECON tables.
OBJ_ID = 'Desig'
ORB_CLASS = 'DES Classification'
EPHEM_UNCERTAINTY = 'TNO pos err'
EVENT_TIME = "ET"


class HTMLTableParser(object):
    """A Parser for an HTML Table, based on example from the BeautifulSoup cookbook."""

    def parse_url(self, url):
        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'html5lib')
        return [(table.get('id', ''), self.parse_html_table(table)) \
                for table in soup.find_all('table')]

    def parse_html_table(self, table):
        """
        Parse the given HTML Table DOM.

        :param table:
        :return: pandas dataframe
        """

        n_columns = 0
        n_rows = 0
        column_names = []
        skip_first_tr = False
        # Find number of rows and columns
        # we also find the column titles if we can
        for row in table.find_all('tr'):
            # Determine the number of rows in the table
            td_tags = row.find_all('td')
            if len(td_tags) > 0:
                n_rows += 1
            if n_columns == 0:
                # Set the number of columns for our table
                n_columns = len(td_tags)

                # Handle column names if we find them
            th_tags = row.find_all('th')
            if len(th_tags) > 0 and len(column_names) == 0:
                for th in th_tags:
                    column_names.append(th.get_text())
                skip_first_tr = False
            # Assume first row is names as <td> objects if not set
            if len(column_names) == 0:
                td_tags = row.find_all('td')
                for td in td_tags:
                    column_names.append(td.get_text())
                n_rows -= 1
                skip_first_tr = True

        # Safeguard on Column Titles
        if len(column_names) > 0 and len(column_names) != n_columns:
            raise Exception("Column titles do not match the number of columns")

        columns = column_names if len(column_names) > 0 else range(0, n_columns)
        df = pd.DataFrame(columns=columns,
                          index=range(0, n_rows))
        row_marker = 0
        for row in table.find_all('tr'):
            if skip_first_tr:
                skip_first_tr = False
                continue
            column_marker = 0
            columns = row.find_all('td')
            for column in columns:
                df.iat[row_marker, column_marker] = column.get_text()
                column_marker += 1
            if len(columns) > 0:
                row_marker += 1

        # Convert to float if possible
        for col in df:
            try:
                df[col] = df[col].astype(float)
            except ValueError:
                pass
        return df


def main():

    # Configure the source of the recon tno lists.
    event_list_url = {}

    # allevents.html lists the best candidates for TNO occultation for the next two years
    # without regard for observing location.
    # Adding some observations might raise the probability of an event for one of these.
    event_list_url['all'] = "{}/allevents.html".format(SERVICE_URL)

    # reconlist contains the best candidates for TNO occultation
    # for the next two years that are visible from the entire network.
    #  These events all have a minimum success probability of 30%.
    #  Near term observation will help confirm the probability
    event_list_url['best'] = "{}/reconlist.html".format(SERVICE_URL)

    # reconwatch.html lists the candidates for TNO occultation for the next
    # two years that are visible from the entire network.
    # observations in the near term will help firm up the occultation probabilities.
    event_list_url['watch'] = "{}/reconwatch.html".format(SERVICE_URL)

    # longlist.html This is a list of TNOs with current positional errors between 2 and 240 arcseconds.
    # These objects are easy to find but have errors too large to permit useful
    # predictions of occultation opportunities.
    event_list_url['long'] = "{}/longlist.html".format(SERVICE_URL)

    # there is also a CSV list, but its everything bunched together.
    # url = "http://www.boulder.swri.edu/~buie/recon/reconlist.csv"

    # These are the classes of TNOs that are in the Buie classification.

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--verbose', help="Verbose message reporting.", action="store_true", default=False)
    parser.add_argument('--debug', help="Provide debuging information.", action="store_true", default=False)
    parser.add_argument('--list', help="Which recon candidate list do you want to check?",
                        choices=['all', 'best', 'watch', 'long'],
                        default='all')
    parser.add_argument('--classes', nargs='*',
                        help='List of classes of objects to select',
                        default=['ERR2LARGE',
                                 'RESONANT',
                                 'CLASSICAL',
                                 'SCATNEAR'],
                        choices=DES_CLASSES)
    parser.add_argument('--min_uncertainty',
                        help="Minimum uncertainty in orbit required to trigger tracking (in arsec)",
                        default=0.1,
                        type=float)
    parser.add_argument('start_time', help="Start of period to look for events.", type=Time)
    parser.add_argument('stop_time', help="End of period to check for events.", type=Time)
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    logging.basicConfig(level=logging.ERROR)

    url = event_list_url[args.list]

    logging.info("Working on events in list: {}".format(url))
    for target in parse_recon_table(url, start_time=args.start_time, end_time=args.stop_time,
                                    orbit_classes=args.classes, min_uncertainty=args.min_uncertainty):
        sys.stdout.write("{}\t{}\t{}\t{}\n".format(target.name, target.mag, target.ra, target.dec))


def parse_recon_table(url, start_time, end_time, orbit_classes, min_uncertainty):
    """Parse the HTML tables distributed by the RECON project.

    """

    # there are some differences in the column used by reconlist.csv and the other lists.
    col_name_mapping = {'Object ID': OBJ_ID,
                        'Type': ORB_CLASS,
                        'PosErr': EPHEM_UNCERTAINTY}

    good_targets = []

    if url.endswith('html'):
        hp = HTMLTableParser()
        ptable = hp.parse_url(url)[0][1]
        table = Table.from_pandas(ptable)
    else:
        if url.startswith('http'):
            fobj = StringIO(requests.get(url).content)
        else:
            fobj = open(url)
        fobj.seek(0)
        table = Table.read(fobj, format='csv')

    # Some of the files have different column names
    for old_name in col_name_mapping:
        if old_name in table.colnames:
            table.rename_column(old_name, col_name_mapping[old_name])

    for row in table:
        if row[ORB_CLASS] not in DES_CLASSES:
            row[ORB_CLASS] = DEFAULT_ORB_CLASS

    # These are the sources that I'm interested in tracking: Kuiper belt objects with larger position uncertainty.

    conds = []
    for orbit_class in orbit_classes:
        conds.append(table[ORB_CLASS] == orbit_class)

    cond = numpy.any(conds, axis=0)
    cond = numpy.all((cond, table[EPHEM_UNCERTAINTY] > min_uncertainty), axis=0)

    # If these are possible Events (ie, not just ephemeris improvement) then only do nearby events
    if EVENT_TIME in table.colnames:
        event_time = [datetime.strptime(date_str, '%Y %b %d %H:%M:%S') for date_str in table[EVENT_TIME]]
        event_time = Time(event_time, scale='utc')
        cond = numpy.all((event_time > start_time,
                          event_time < end_time,
                          cond), axis=0)
    else:
        table[EVENT_TIME] = Time(str(date.today())).iso

    table = table[cond]

    cfht = ephem.Observer()
    cfht.lat = 0.344
    cfht.lon = -2.707
    cfht.elevation = 4100
    cfht.date = '2017/10/13 20:00:00'
    cfht.date = start_time.iso.replace("-", "/")
    cfht.horizon = math.radians(-7)

    sun_set_time = Time(str(ephem.Sun(cfht).set_time).replace("/", "-"), scale='utc')
    sun_rise_time = Time(str(ephem.Sun(cfht).rise_time).replace("/", "-"), scale='utc')

    cfht.horizon = MINIMUM_ELEVATION.to('radian').value

    logging.info("Table at {} contains {} matching entries.".format(url, len(table)))

    count = 0
    for row in table:
        target = ephem.FixedBody()
        count += 1
        try:
            name = "{:.0f}".format(float(row[OBJ_ID].split()[0]))
        except ValueError as ve:
            logging.debug(str(ve))
            if float(row[OBJ_ID][0:2]) < 19:
                century = 20
            else:
                century = 19
            name = '{}{} {}'.format(century, row[OBJ_ID][0:2], row[OBJ_ID][2:])
        logging.info("Doing object {} from row: {}".format(name, count))

        o = horizons.Body(name, start_time=start_time, stop_time=end_time, center='@568')
        o.predict(sun_set_time)
        logging.debug("Getting coordinates from Horizons.")
        target._ra = o.coordinate.ra.radian
        target._dec = o.coordinate.dec.radian
        target.name = name
        target.mag = o.mag
        target.compute(cfht)
        target_rise_time = Time(str(target.rise_time).replace("/", "-"), scale='utc')
        target_set_time = Time(str(target.set_time).replace("/", "-"), scale='utc')

        start = end = None
        if target_rise_time > sun_rise_time:
            if target_set_time > sun_set_time:
                # target rises after sun_rise but target_set after sun_set
                start = sun_set_time
                if target_set_time < sun_rise_time:
                    end = target_set_time
                else:
                    end = sun_rise_time
        elif target_rise_time < sun_set_time:
            start = sun_set_time
            if target_set_time < sun_rise_time:
                end = target_set_time
            else:
                end = sun_rise_time

        elif sun_set_time < target_rise_time < sun_rise_time:
            start = target_rise_time
            if target_set_time < sun_rise_time:
                end = target_set_time
            else:
                end = sun_rise_time

        duration = (end - start).to(units.hour)
        if duration < MINIMUM_UP_TIME:
            logging.info("Skipping traget {}:  only up for {} hours ".format(name, duration))
            logging.debug("Rise time: {}, Set time: {}".format(start, end))
            continue
        good_targets.append(target)
        logging.debug("{:12s} {:12s} {:10s} {:12s} {:5.2f} {:12s} {:12s} {:5.1f}".format(str(target.ra),
                                                                                         str(target.dec),
                                                                                         row[OBJ_ID],
                                                                                         row[EVENT_TIME],
                                                                                         row[EPHEM_UNCERTAINTY],
                                                                                         str(target_rise_time),
                                                                                         str(target_set_time),
                                                                                         duration))
    return good_targets


main.__doc__ = DESCRIPTION

if __name__ == '__main__':
    sys.exit(main())
