""" BlueSky plugin template. The text you put here will be visible
    in BlueSky as the description of your plugin. """
# Import the global bluesky objects. Uncomment the ones you need
from bluesky import stack, sim, traf  #, settings, navdb, traf, sim, scr, tools
import time
import random
from bluesky.tools import areafilter, datalog, plotter, geo, TrafficArrays
from bluesky.tools.aero import nm, ft
import matplotlib.pyplot as plt
import bluesky as bs
import bluesky.traffic.performance.openap as openap
import numpy as np
import pandas as pd
from pathlib import Path
import bluesky as bs
import random as rnd
import matplotlib as matplotlib
from bluesky.traffic.asas import ConflictDetection
from bluesky.traffic.asas import ConflictResolution
from bluesky.traffic.asas import StateBased
import os


### Initialization function of your plugin. Do not change the name of this
### function, as it is the way BlueSky recognises this file as a plugin.
def init_plugin():
    global recording, prevconfpairs, first_confframe, resos, reso_index, first_delframe, ntrafframe, first_ntrafframe, locframe, first_locframe, first_logger, headings, densities, density_index, heading_index, reso, superfirst, i
    recording = False
    prevconfpairs = []
    first_confframe = True
    first_delframe = True
    first_logger = True
    first_ntrafframe = True
    first_locframe = True
    headings = [15,20,25]
    #densities = [5,10,20,30,50]
    densities = [5]
    density_index = 0
    heading_index = 0
    i = 1
    superfirst = True
    resos = ['MVPSPD','MVPSPDFTR','MVPVREL','MVPVRELFTR']
    reso_index = 0
    # Configuration parameters
    config = {
        # The name of your plugin
        'plugin_name':     'LS5Logger2',

        # The type of this plugin. For now, only simulation plugins are possible.
        'plugin_type':     'sim',

        # Update interval in seconds. By default, your plugin's update function(s)
        # are called every timestep of the simulation. If your plugin needs less
        # frequent updates provide an update interval.
        'update_interval': 1,

        # The update function is called after traffic is updated. Use this if you
        # want to do things as a result of what happens in traffic. If you need to
        # something before traffic is updated please use preupdate.
        'update':          update,

        # The preupdate function is called before traffic is updated. Use this
        # function to provide settings that need to be used by traffic in the current
        # timestep. Examples are ASAS, which can give autopilot commands to resolve
        # a conflict.
        'preupdate':       preupdate,

        # If your plugin has a state, you will probably need a reset function to
        # clear the state in between simulations.
        'reset':         reset
        }

    stackfunctions = {
        # The command name for your function
        }


    # init_plugin() should always return these two dicts.
    return config, stackfunctions


### Periodic update functions that are called by the simulation. You can replace
### this by anything, so long as you communicate this in init_plugin

def update():
    global recording, prevconfpairs, confframe, delframe, resos, reso_index, ntrafframe, first_ntrafframe, first_confframe, locframe, first_locframe, first_delframe, first_logger, headings, densities, heading_index, density_index, reso, superfirst, i
    stack.stack('FF')
    if superfirst:
        superfirst = False
        stack.stack('PLUGIN MVPSPD')
        stack.stack('PLUGIN MVPSPDFTR')
        stack.stack('PLUGIN MVPVREL')
        stack.stack('PLUGIN MVPVRELFTR')
        stack.stack('ASAS ON')
        stack.stack('CDMETHOD CStateBased')
        stack.stack('RESO '+resos[reso_index])
        stack.stack('PCALL '+'Wouter_MTM/'+str(headings[heading_index])+'_'+str(densities[density_index])+'_'+str(i)+'.SCN')
        stack.stack('FF')
    if traf.ntraf > 0 and not recording:                                                                          # Voor low density case.
        recording = True
    else:
        pass

    tlist_id = [x for x in traf.id if x[0] == 'T']
    tlist_idx = traf.id2idx(tlist_id)

    if first_ntrafframe:
        ntrafframe = np.array([sim.simt, traf.ntraf])
        first_ntrafframe = False
    else:
        ntrafframe = np.vstack((ntrafframe, np.array([sim.simt, traf.ntraf])))

    if first_locframe:
        locframe = np.array([sim.simt, traf.lat, traf.lon])
        first_locframe = False
    else:
        if sim.simt % 5*60 < 1:
            locframe = np.vstack((locframe, np.array([sim.simt, traf.lat, traf.lon])))
        else:
            pass

    if recording:                                                                                               # Omgekeerd Statement, de not moet eigenlijk weg.
        confpairs, lospairs, inconf, tcpamax, qdr, dist, dcpa, tcpa, tLOS = \
        StateBased.detect(StateBased, traf, traf, bs.traf.cd.rpz, bs.traf.cd.hpz, bs.traf.cd.dtlookahead)
        #newgscapped = ConflictResolution.resolve(traf.asas,traf)
        if len([x for x in confpairs if x not in prevconfpairs]) > 0:
            newcomers = [confpairs.index(i) for i in confpairs if i not in prevconfpairs]
            for n in newcomers:
                confname = confpairs[n]
                ac1 = confname[0]
                ac2 = confname[1]
                ac1idx = traf.id2idx(ac1)
                ac2idx = traf.id2idx(ac2)
                bearing = qdr[n]
                distance = dist[n]
                distancecpa = dcpa[n]
                timecpa = tcpa[n]
                timeLOS = tLOS[n]
                initialtas = traf.tas[ac1idx]
                initialhdg = traf.hdg[ac1idx]
                initialtasi = traf.tas[ac2idx]
                initialhdgi = traf.hdg[ac2idx]
                latitude = traf.lat[ac1idx]
                longitude = traf.lon[ac1idx]
                latitudei = traf.lat[ac2idx]
                longitudei = traf.lon[ac2idx]

                if first_confframe:
                    confframe = np.array([sim.simt, confname, ac1idx, ac1, bearing, distance, distancecpa, timecpa, timeLOS, initialtas, initialhdg, initialtasi, initialhdgi, latitude, longitude, latitudei, longitudei])
                    first_confframe = False
                else:
                    confframe = np.vstack((confframe, np.array([sim.simt, confname, ac1idx, ac1, bearing, distance, distancecpa, timecpa, timeLOS, initialtas, initialhdg, initialtasi, initialhdgi, latitude, longitude, latitudei, longitudei])))
        prevconfpairs = confpairs
        if sim.simt > 3.5*3600 and not first_confframe and not first_delframe: # Moet 4 ofzo zijn
            if not os.path.exists('/BSData2/Dens'+str(densities[density_index])+'/Head'+str(headings[heading_index])):
                os.makedirs('/BSData2/Dens'+str(densities[density_index])+'/Head'+str(headings[heading_index]), exist_ok=True)

            #pd.DataFrame(confframe).to_csv(str(headings[heading_index])+'_'+str(densities[density_index])+'_'+str(i)+'_'+reso+'_ResultsDF.csv', sep=';')
            pd.DataFrame(confframe).to_csv('/BSData2/Dens'+str(densities[density_index])+'/Head'+str(headings[heading_index])+'/'+str(i)+'-'+resos[reso_index]+'-confframe.csv' ,sep=';')

            #pd.DataFrame(delframe).to_csv(str(headings[heading_index])+'_'+str(densities[density_index])+'_'+str(i)+'_'+reso+'_Results-DelsDF.csv', sep=';')
            pd.DataFrame(delframe).to_csv('/BSData2/Dens'+str(densities[density_index])+'/Head'+str(headings[heading_index])+'/'+str(i)+'-'+resos[reso_index]+'-delframe.csv' ,sep=';')

            #pd.DataFrame(ntrafframe).to_csv(str(headings[heading_index])+'_'+str(densities[density_index])+'_'+str(i)+'_'+reso+'_Results-ntrafDF.csv', sep=';')
            pd.DataFrame(ntrafframe).to_csv('/BSData2/Dens'+str(densities[density_index])+'/Head'+str(headings[heading_index])+'/'+str(i)+'-'+resos[reso_index]+'-ntrafframe.csv' ,sep=';')

            #pd.DataFrame(locframe).to_csv(str(headings[heading_index])+'_'+str(densities[density_index])+'_'+str(i)+'_'+reso+'_Results-locsDF.csv', sep=';')
            pd.DataFrame(locframe).to_csv('/BSData2/Dens'+str(densities[density_index])+'/Head'+str(headings[heading_index])+'/'+str(i)+'-'+resos[reso_index]+'-locframe.csv' ,sep=';')

            #stack.stack('HOLD')
            stack.stack('RESET')

            reso_index += 1

            c = False
            if reso_index > (len(resos) - 1):
                reso_index = 0
                i += 1
                if i > 10: #10
                    i = 1
                    heading_index += 1
                    if heading_index > 9: #9
                        heading_index = 0
                        density_index += 1
                        if density_index > 0: #9
                            stack.stack('HOLD')
                            print('complete')
                            c = True

            if not c:
                stack.stack('CDMETHOD CStateBased')
                stack.stack('RESO '+resos[reso_index])
                stack.stack('PCALL '+'Wouter_MTM/'+str(headings[heading_index])+'_'+str(densities[density_index])+'_'+str(i)+'.SCN')
                stack.stack('FF')
            else:
                stack.stack('HOLD')
            first_confframe = True
            first_delframe = True
            first_ntrafframe = True
            recording = False
        
    for aircraft in traf.id:
        if traf.actwp.lon[traf.id2idx(aircraft)] > 10:
            if first_delframe:
                delframe = np.array([sim.simt, aircraft, traf.lat[traf.id2idx(aircraft)], traf.lon[traf.id2idx(aircraft)], traf.distflown[traf.id2idx(aircraft)]])
                first_delframe = False       
            else:
                delframe = np.vstack((delframe, np.array([sim.simt, aircraft, traf.lat[traf.id2idx(aircraft)], traf.lon[traf.id2idx(aircraft)], traf.distflown[traf.id2idx(aircraft)]])))                                                      
            
            traf.delete(traf.id2idx(aircraft))
        else:
            pass

def dvsaver(ids, dv):
    global df_dv
     
def preupdate():
    pass

def reset():
    pass

