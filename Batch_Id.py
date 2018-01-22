import os
import sys

Gauges = []
for i in range(2,16):
    if i != 3:
        Gauges.append("WOR_domain_%g" % (i))  
print Gauges

for gauge in Gauges:
    print 'launching the marsh detector for %s' % (gauge)
    os.system('python MarshPlatformAnalysis.py -dir /home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshAnalysis/Input/ -sites ' + gauge + ' -MID True -MIDP False &')