''' Eby ConflictResolution implementation plugin. '''
import numpy as np
from bluesky.tools.aero import vtas2eas
from bluesky.tools import geo
from bluesky.traffic.asas import MVP
from bluesky.traffic.asas import StateBased
import bluesky as bs


def init_plugin():

    # Addtional initilisation code

    # Configuration parameters
    config = {
        # The name of your plugin
        'plugin_name':     'MVPFTR',

        # The type of this plugin. For now, only simulation plugins are possible.
        'plugin_type':     'sim'
    }

    # init_plugin() should always return these two dicts.
    return config, {}


class MVPFTR(MVP):
    def resumenav(self, conf, ownship, intruder):
        '''
            Decide for each aircraft in the conflict list whether the ASAS
            should be followed or not, based on if the aircraft pairs passed
            their CPA.
        '''
        pass
        # Add new conflicts to resopairs and confpairs_all and new losses to lospairs_all
        self.resopairs.update(conf.confpairs)

        # Conflict pairs to be deleted
        delpairs = set()
        changeactive = dict()

        # Look at all conflicts, also the ones that are solved but CPA is yet to come
        for conflict in self.resopairs:
            idx1, idx2 = bs.traf.id2idx(conflict)
            # If the ownship aircraft is deleted remove its conflict from the list
            if idx1 < 0:
                delpairs.add(conflict)
                continue

            if idx2 >= 0:
                # Distance vector using flat earth approximation
                re = 6371000.
                dist = re * np.array([np.radians(intruder.lon[idx2] - ownship.lon[idx1]) *
                                      np.cos(0.5 * np.radians(intruder.lat[idx2] +
                                                              ownship.lat[idx1])),
                                      np.radians(intruder.lat[idx2] - ownship.lat[idx1])])

                # Relative velocity vector
                vrel = np.array([intruder.gseast[idx2] - ownship.gseast[idx1],
                                 intruder.gsnorth[idx2] - ownship.gsnorth[idx1]])

                # Check if conflict is past CPA
                past_cpa = np.dot(dist, vrel) > 0.0

                # Free to Resume Method
                des_aphdg_1 = bs.traf.ap.trk[idx1]
                des_aphdg_2 = bs.traf.ap.trk[idx2]

                des_apspd_1 = bs.traf.ap.tas[idx1]
                des_apspd_2 = bs.traf.ap.tas[idx2]

                des_gsnorth_1 = des_apspd_1 * np.cos(np.radians(des_aphdg_1))
                des_gseast_1 = des_apspd_1 * np.sin(np.radians(des_aphdg_1))

                des_gsnorth_2 = des_apspd_2 * np.cos(np.radians(des_aphdg_2))
                des_gseast_2 = des_apspd_2 * np.sin(np.radians(des_aphdg_2))
                  
                des_vrel = np.array([des_gseast_2 - des_gseast_1,
                                     des_gsnorth_2 - des_gsnorth_1])

                des_tcpa_1 = np.maximum(-(des_vrel[0] * dist[0] + des_vrel[1] * dist[1]) / (des_vrel[0]*des_vrel[0] + des_vrel[1]*des_vrel[1]),0.0)
                des_dcpa_1 = dist + des_vrel*des_tcpa_1
                free = np.linalg.norm(des_dcpa_1) > bs.traf.cd.rpz
                # End of free to resume method

                # hor_los:
                # Aircraft should continue to resolve until there is no horizontal
                # LOS. This is particularly relevant when vertical resolutions
                # are used.
                hdist = np.linalg.norm(dist)
                hor_los = hdist < conf.rpz

                # Bouncing conflicts:
                # If two aircraft are getting in and out of conflict continously,
                # then they it is a bouncing conflict. ASAS should stay active until
                # the bouncing stops.
                is_bouncing = abs(
                    ownship.trk[idx1] - intruder.trk[idx2]) < 30.0 and hdist < conf.rpz * self.resofach

            # Start recovery for ownship if intruder is deleted, or if past CPA
            # and not in horizontal LOS or a bouncing conflict
            if idx2 >= 0 and (not free or hor_los or is_bouncing):
                # Enable ASAS for this aircraft
                changeactive[idx1] = True
            else:
                # Switch ASAS off for ownship if there are no other conflicts
                # that this aircraft is involved in.
                changeactive[idx1] = changeactive.get(idx1, False)
                # If conflict is solved, remove it from the resopairs list
                delpairs.add(conflict)

        for idx, active in changeactive.items():
            # Loop a second time: this is to avoid that ASAS resolution is
            # turned off for an aircraft that is involved simultaneously in
            # multiple conflicts, where the first, but not all conflicts are
            # resolved.
            self.active[idx] = active
            if not active:
                # Waypoint recovery after conflict: Find the next active waypoint
                # and send the aircraft to that waypoint.
                iwpid = bs.traf.ap.route[idx].findact(idx)
                if iwpid != -1:  # To avoid problems if there are no waypoints
                    bs.traf.ap.route[idx].direct(
                        idx, bs.traf.ap.route[idx].wpname[iwpid])

        # Remove pairs from the list that are past CPA or have deleted aircraft
        self.resopairs -= delpairs