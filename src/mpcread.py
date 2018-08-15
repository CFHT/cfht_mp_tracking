#!/usr/bin/env python

import sys,re
import math, numpy
from astropy.io import ascii
from astropy.time import Time


uncertainty =  numpy.array([1.0, 4.4, 19.6, 86.5, 382, 1692, 7488, 33121, 146502, 146502*2])/10.0


current_time = Time("2017-09-01").jd

Number_Mil={'B': 110000, 'C': 120000, 'D': 130000, 'E': 140000, 'F': 150000}
Number_Cent={'J': 1900, 'K': 2000}
Ncode='0123456789ABCDEFGHIJKLMNOPQRSTUV'
Ncode='0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
Kilo='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
YY={'I': 1800, 'J': 1900, 'K': 2000}


def date_unpack(pdate):
    yyyy=YY[pdate[0]]+int(pdate[1:3])
    mm=Ncode.rindex(pdate[3])
    dd=float(Ncode.rindex(pdate[4]))
    return (yyyy, mm, dd)

def desig_pack(desig):
    try:
        f = int(desig)
        if f < 100000:
            return desig
        else:
            digits = f % 10000
            pow = Ncode[(f - digits)/10000]
            return "{}{:04d}".format(pow, digits)
    except:
        c = Ncode[int(desig[0:2])]
        y = desig[2:4]
        month = desig[5]
        l2 = desig[6]
        cycle = int(desig[7:])
        digit = cycle % 10
        return "{}{}{}{}{}{}".format(c, y, month, Ncode[(cycle - digit)/10], digit, l2)


def desig_unpack(desig):

    import re
    if re.match('\d+',desig):
        return str(int(desig))
    j = Kilo.rfind(desig[0])
    if j != -1 :
      if re.match('^\d+$',desig[1:7]):
        return str(100000+j*10000+int(desig[1:7]))
    try:
       yyyy=str(YY[desig[0]]+int(desig[1:3]))
    except KeyError as e:
       return desig
    Mcode=desig[3]+desig[6]
    cycle=Ncode.rindex(desig[4])*10+int(desig[5])
    if (cycle>0) :
	cycle=str(cycle)
    else:
	cycle=''
    return yyyy+' '+Mcode+cycle
    

def main():
    f=open('/Users/kavelaarsj/MPCORB.DAT')
    #f=open('schwamb_orbits.dat')
    lines=f.readlines()
    f.close()
    
    import ephem,sys
    kbo=ephem.EllipticalBody()
    line=lines.pop(0)
    while (line[0:3]!="---") : 
       line=lines.pop(0)
    
    lines.append(line)
    nobj=0
    lineCount=0
    for line in lines:
        lineCount=lineCount+1
        if lineCount %1000 == 0 : 
    	    sys.stderr.write("# Line: %d \n" % ( lineCount)) 
        if line[0]=='#' or len(line) < 103 or line[0:3]=='---':
            #sys.stderr.write(line)
            continue
        
        if len(line[8:13].strip()):
            kbo._H=float(line[8:13])
            kbo._G=float(line[14:19])
        else:
            kbo._H=20
            kbo._G=0
        arc = line[127:136]
        try:
            if 'days' in arc:
    	        arc = int(arc.split()[0])/365.25
            else:  
                dt = -eval(arc)
        except:
            sys.stderr.write("Error parsing the arc length value: {}".format(arc))
            continue
        kbo._epoch_M=date_unpack(line[20:25].strip())
        kbo._M=float(line[26:35].strip())
        kbo._om=float(line[37:46].strip())
        kbo._Om=float(line[48:57].strip())
        kbo._inc=float(line[59:68].strip())
        kbo._e=float(line[70:79].strip())
        kbo._epoch='2017/09/04'
        kbo._a=float(line[92:103].strip())
        try:
           U = int(line[105])
           nobs = int(line[117:122])
           last_obs = Time("{}-{}-{}".format(line[194:198],line[198:200],line[200:202]))
        except:
           print line
           U = 9
           nobs = 3
           last_obs = Time("2017-01-01")
        pU = uncertainty[U]
        pU =  (current_time - last_obs.jd)*U/365.25
        #print last_obs, U, uncertainty[U], pU
        a= kbo._a
        e= kbo._e
        H= kbo._H
        i= kbo._inc
        kbo.compute(ephem.date('2011/12/22'))
        kbo.name=desig_unpack(line[0:7].strip()) 
        T_J = (5.2/a) + 2.0 * math.sqrt((1-e**2)*(a/5.2)) * math.cos(i)
        if eval(cond): 
           if line[0]=='P' or line[0]=='T':
              # Ignore the PLS and T (?) astroid surveys.
              continue
           #try : 
           kbo.name=desig_unpack(line[0:7].strip()) 
           #except : 
           #   sys.stderr.write(line[0:7].strip())
           #nobj = nobj+1
           print "%20s %5.1f %5.1f %5.1f %f %f %f" % ( kbo.name.replace(" ","_"), a, e, math.degrees(i), H, math.degrees(kbo.ra), math.degrees(kbo.dec) )
           for column in columns:
               out_data[column].append(eval(column))
           print "%20s %5.1f %5.1f %5.1f %f %f %f" % ( kbo.name.replace(" ","_"), a, e, math.degrees(i), H, math.degrees(kbo.ra), math.degrees(kbo.dec) )
           # print "%20s %5.1f %5.1f %5.3f %5.1f %5.1f %5.1f %5.1f %s %s" % (kbo.name,H ,a,e,a*(1-e),i*57.3,kbo.sun_distance,kbo.mag,kbo.ra,kbo.dec)
           # print "grep \"%s\" schwamb_T2.dat" % ( kbo.name)
           # print math.degrees(kbo._Om)
           # print "# ",kbo.sun_distance, kbo._a, kbo._e, kbo._inc, kbo._a, kbo.sun_distance
           # print "%10.5f,%10.5f,%-20s,%10.2f,%10.2f" % ( math.degrees(kbo.ra),math.degrees(kbo.dec),"'"+kbo.name+"'", kbo._a, kbo.sun_distance)
           # print H,kbo._a, kbo._e, math.degrees(kbo._inc), kbo.sun_distance,a*(1-e),a*(1+e), kbo.name.replace(" ","_")
           # print kbo.writedb()  
    
    ascii.write(out_data, 'mpcread.dat', names=columns)
        
    #print nobj

if __name__ == '__main__':
    cond=sys.argv[1]
    columns = sys.argv[1:]
    out_data = {}
    for column in columns:
       out_data[column] = []
    print "# "+cond
    main()
