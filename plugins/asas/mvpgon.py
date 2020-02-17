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
        'plugin_name':     'MVPGON',

        # The type of this plugin. For now, only simulation plugins are possible.
        'plugin_type':     'sim'
    }

    # init_plugin() should always return these two dicts.
    return config, {}


class MVPGON(MVP):
    def MVP(self, ownship, intruder, conf, qdr, dist, tcpa, tLOS, idx1, idx2):
        """Modified Voltage Potential (MVP) resolution method"""
        # Preliminary calculations-------------------------------------------------

        # Convert qdr from degrees to radians
        qdr = np.radians(qdr)

        # Relative position vector between id1 and id2
        drel = np.array([np.sin(qdr)*dist, \
                        np.cos(qdr)*dist, \
                        intruder.alt[idx2]-ownship.alt[idx1]])

        # Write velocities as vectors and find relative velocity vector
        v1 = np.array([ownship.gseast[idx1], ownship.gsnorth[idx1], ownship.vs[idx1]])
        v2 = np.array([intruder.gseast[idx2], intruder.gsnorth[idx2], intruder.vs[idx2]])
        vrel = np.array(v2-v1)


        # Horizontal resolution----------------------------------------------------

        # Find horizontal distance at the tcpa (min horizontal distance)
        dcpa  = drel + vrel*tcpa
        dabsH = np.sqrt(dcpa[0]*dcpa[0]+dcpa[1]*dcpa[1])

        # Compute horizontal intrusion
        iH = (conf.rpz * self.resofach) - dabsH

        # Exception handlers for head-on conflicts
        # This is done to prevent division by zero in the next step
        if dabsH <= 10.:
            dabsH = 10.
            dcpa[0] = drel[1] / dist * dabsH
            dcpa[1] = -drel[0] / dist * dabsH

        R = np.linalg.norm(drel)/2
        d = dabsH
        i = (conf.rpz * self.resofach) - d
        h = np.sqrt((2*R)**2-d**2)
        alpha = np.arcsin(d/(2*R))

        x1 = (h-R*np.sin(np.arccos(np.sin(alpha)+i/R))-R*np.cos(alpha)) * np.tan((1/4)*np.pi - alpha/2 - (1/2)*np.arccos(np.sin(alpha)+i/R))
        if x1 == x1:
            x = x1
        else:
            x = 0

        # If intruder is outside the ownship PZ, then apply extra factor
        # to make sure that resolution does not graze IPZ
        if (conf.rpz * self.resofach) < dist and dabsH < dist:
            # Compute the resolution velocity vector in horizontal direction.
            # abs(tcpa) because it bcomes negative during intrusion
            if x != 0:
                dv1 = ((i + x)*dcpa[0])/(abs(tcpa)*dabsH)
                dv2 = ((i + x)*dcpa[1])/(abs(tcpa)*dabsH)
            else:
                erratum=np.cos(np.arcsin((conf.rpz * self.resofach)/dist)-np.arcsin(dabsH/dist))
                dv1 = (((conf.rpz * self.resofach)/erratum - dabsH)*dcpa[0])/(abs(tcpa)*dabsH)
                dv2 = (((conf.rpz * self.resofach )/erratum - dabsH)*dcpa[1])/(abs(tcpa)*dabsH)
        else:
            dv1 = (iH * dcpa[0]) / (abs(tcpa) * dabsH)
            dv2 = (iH * dcpa[1]) / (abs(tcpa) * dabsH)

        # Vertical resolution------------------------------------------------------

        # Compute the  vertical intrusion
        # Amount of vertical intrusion dependent on vertical relative velocity
        iV = (conf.hpz * self.resofacv) if abs(vrel[2])>0.0 else (conf.hpz * self.resofacv)-abs(drel[2])

        # Get the time to solve the conflict vertically - tsolveV
        tsolV = abs(drel[2]/vrel[2]) if abs(vrel[2])>0.0 else tLOS

        # If the time to solve the conflict vertically is longer than the look-ahead time,
        # because the the relative vertical speed is very small, then solve the intrusion
        # within tinconf
        if tsolV>conf.dtlookahead:
            tsolV = tLOS
            iV    = (conf.hpz * self.resofacv)

        # Compute the resolution velocity vector in the vertical direction
        # The direction of the vertical resolution is such that the aircraft with
        # higher climb/decent rate reduces their climb/decent rate
        dv3 = np.where(abs(vrel[2])>0.0,  (iV/tsolV)*(-vrel[2]/abs(vrel[2])), (iV/tsolV))

        # It is necessary to cap dv3 to prevent that a vertical conflict
        # is solved in 1 timestep, leading to a vertical separation that is too
        # high (high vs assumed in traf). If vertical dynamics are included to
        # aircraft  model in traffic.py, the below three lines should be deleted.
    #    mindv3 = -400*fpm# ~ 2.016 [m/s]
    #    maxdv3 = 400*fpm
    #    dv3 = np.maximum(mindv3,np.minimum(maxdv3,dv3))


        # Combine resolutions------------------------------------------------------

        # combine the dv components
        dv = np.array([dv1,dv2,dv3])

        return dv, tsolV