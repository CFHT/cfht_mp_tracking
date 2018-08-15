import argparse
import minor_planet_ephemeris
import sys
from bs4 import BeautifulSoup
import requests
import pandas as pd
from datetime import date
import time
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

max_per_iter = 7

class HTMLTableParser(object):

    def parse_url(self, url):
        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'html5lib')
        return [(table.get('id', ''), self.parse_html_table(table)) \
                for table in soup.find_all('table')]

    def parse_html_table(self, table):
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

### Configure the source of the recon tno lists.
eventListUrl = {}
# allevents lists the best candidates for TNO occultation for the next two years
# without regard for observing location.
# Adding some observations might riase the probabilitiy abouve one.
eventListUrl['all'] = "http://www.boulder.swri.edu/~buie/recon/allevents.html"

# reconlist contains the best candidates for TNO occultation
# for the next two years that are visible from the entire network.
#  These events all have a minimum success probability of 30%.
#  Near term observation will help confirm the probabilty
eventListUrl['best'] = "http://www.boulder.swri.edu/~buie/recon/reconlist.html"

# reconwatch lists the candidates for TNO occultation for the next two years that are visible from the entire network.
# observations in the near term will help firm up the occultation probabilities.
eventListUrl['watch'] = "http://www.boulder.swri.edu/~buie/recon/reconwatch.html"

# This is a list of TNOs with current positional errors between 2 and 240 arcseconds.
# These objects are easy to find but have errors too large to permit useful predictions of occultation opportunities.
eventListUrl['long'] = "http://www.boulder.swri.edu/~buie/recon/longlist.html"
# url = "http://www.boulder.swri.edu/~buie/recon/reconlist.csv"

# there are some differences in the column used by reconlist and the other lists.
col_name_mapping = {'Object ID': 'Desig',
                    'Type': 'DES Classification',
                    'PosErr': 'TNO pos err'}


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('start_time', help="Start of current observing period.")
    parser.add_argument('stop_time', help="End of current observing period.")
    parser.add_argument('--target', help="Which recon candidate list do you want to check?")
    args = parser.parse_args()
    import mpcread
    target = args.target
    minor_planet_ephemeris.build_ephem_files(target, args.start_time, 
                                             args.stop_time, 
                                             ephem_format="CFHT API", 
                                             runid='17BC08')



def parse_recon_table(url, iter=0):
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

    # These are the sources that I'm interested in tracking
    cond = numpy.all((table['DES Classification'] != 'CENTAURR',
                      table['TNO pos err'] > 0.10), axis=0)

    # If these are possible Events (ie, not just ephemeris improvement) then only do nearby events
    if 'ET' in table.colnames:
        event_time = [datetime.strptime(date_str, '%Y %b %d %H:%M:%S') for date_str in table['ET']]
        event_time = Time(event_time, scale='utc')
        cond = numpy.all((event_time > Time('2018-01-15 00:00:00', scale='utc'),
                          event_time < Time('2018-07-15 00:00:00', scale='utc'),
                          cond), axis=0)
    else:
        table['ET'] = Time(str(date.today())).iso

    table = table[cond]

    cfht = ephem.Observer()
    cfht.lat = 0.344
    cfht.lon = -2.707
    cfht.elevation = 4100
    cfht.date = '2017/10/13 20:00:00'
    cfht.horizon = math.radians(-7)

    sun_set_time = Time(str(ephem.Sun(cfht).set_time).replace("/", "-"), scale='utc')
    sun_rise_time = Time(str(ephem.Sun(cfht).rise_time).replace("/", "-"), scale='utc')

    target = ephem.FixedBody()
    cfht.horizon = math.radians(40)

    print("Table at {} contains {} entries.".format(url, len(table)))

    count = 0
    for row in table[iter*max_per_iter:]:
        count += 1
        try:
            name = "{:.0f}".format(float(row['Desig'].split()[0]))
        except Exception as e:
            if float(row['Desig'][0:    2]) < 17:
                century = 20
            else:
                century = 19
            name = '{}{} {}'.format(century, row['Desig'][0:2], row['Desig'][2:])
        logging.info("Doing {} from row: {}".format(name, count+iter*max_per_iter-1))

        o = horizons.Body(name, start_time=sun_set_time, stop_time=sun_rise_time)
        o.predict(sun_set_time)
        logging.info("Getting coordinates from Horizons.")
        target._ra = o.coordinate.ra.radian
        target._dec = o.coordinate.dec.radian
        target.compute(cfht)
        target_rise_time = Time(str(target.rise_time).replace("/", "-"), scale='utc')
        target_set_time = Time(str(target.set_time).replace("/", "-"), scale='utc')

        if target_rise_time > sun_rise_time:
            if target_set_time > sun_set_time:
                # target rises after sun_rise but target_set after sun_set
                start = sun_set_time
                if target_set_time < sun_rise_time:
                    end = target_set_time
                else:
                    end = sun_rise_time
            else:
                start = end = target_set_time
        elif target_rise_time < sun_set_time:
            start = sun_set_time
            if target_set_time < sun_rise_time:
                end = target_set_time
            else:
                end = sun_rise_time

        elif target_rise_time > sun_set_time and target_rise_time < sun_rise_time:
            start = target_rise_time
            if target_set_time < sun_rise_time:
                end = target_set_time
            else:
                end = sun_rise_time

        duration = (end - start).to(units.hour)
        if duration < 1.0 * units.hour:
            print("Skipping traget {}:  only up for {} hours ".format(name, duration))
            print("Rise time: {}, Set time: {}".format(start, end))
            if count == max_per_iter:
                break
            else:
                continue
        good_targets.append(name)
        sys.stdout.write("{:12s} {:12s} {:10s} {:12s} {:5.2f} {:12s} {:12s} {:5.1f}\n".format(str(target.ra),
                                                                                            str(target.dec),
                                                                                            row['Desig'],
                                                                                            row['ET'],
                                                                                            row['TNO pos err'],
                                                                                            str(target_rise_time),
                                                                                            str(target_set_time),
                                                                                            duration))
        if count == max_per_iter:
            break
    return good_targets

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
