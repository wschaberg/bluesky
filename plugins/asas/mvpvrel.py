''' Eby ConflictResolution implementation plugin. '''
import numpy as np
from bluesky.tools.aero import vtas2eas
from bluesky.tools import geo
from bluesky.traffic.asas import MVP
from bluesky.traffic.asas import StateBased
import bluesky as bs
from bluesky.tools import aero
import bluesky.traffic.performance.openap as ap


def init_plugin():

    # Addtional initilisation code

    # Configuration parameters
    config = {
        # The name of your plugin
        'plugin_name':     'MVPVREL',

        # The type of this plugin. For now, only simulation plugins are possible.
        'plugin_type':     'sim'
    }

    # init_plugin() should always return these two dicts.
    return config, {}


class MVPVREL(MVP):
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

        def length(v):
            return np.sqrt(v[0]**2+v[1]**2)
        def dot_product(v,w):
           return v[0]*w[0]+v[1]*w[1]
        def determinant(v,w):
           return v[0]*w[1]-v[1]*w[0]
        def inner_angle(v,w):
           cosx = dot_product(v,w)/(length(v)*length(w))
           rad=np.arccos(cosx) # in radians
           return rad*180/np.pi # returns degrees
        def angle_clockwise(A, B):
            inner = inner_angle(A,B)
            det = determinant(A,B)
            if det<0: #this is a property of the det. If the det < 0 then B is clockwise of A
                return inner
            else: # if the det > 0 then A is immediately clockwise of B
                return -inner

        vavg = np.array([np.mean(np.delete(ownship.gseast, [idx1, idx2], 0)), np.mean(np.delete(ownship.gsnorth, [idx1, idx2], 0)), np.mean(np.delete(ownship.vs, [idx1, idx2], 0))])
        
        # if angle_clockwise(vrel, drel) > 180:
        #     phi = np.radians(angle_clockwise(dcpa[:-1], vavg[:-1]))
        #     if abs(phi) > 0.5*np.pi:
        #         phi = np.radians(angle_clockwise(dcpa[:-1], vavg[:-1]))
        #     else:
        #         pass
        # else:
        #     phi = -np.radians(angle_clockwise(dcpa[:-1], vavg[:-1]))
        #     if abs(phi) > 0.5*np.pi:
        #         phi = -np.radians(angle_clockwise(dcpa[:-1], vavg[:-1]))
        #     else:
        #         pass

        # absvrel = np.linalg.norm(vrel)
        # tcpa_s = tcpa + (dabsH * np.tan(phi))/absvrel
        # dcpa_s  = drel + vrel*tcpa_s
        # dabsH_s = np.sqrt(dcpa_s[0]*dcpa_s[0]+dcpa_s[1]*dcpa_s[1])
        # iH_s = (conf.rpz * self.resofach) - dabsH_s

        # If intruder is outside the ownship PZ, then apply extra factor
        # to make sure that resolution does not graze IPZ
        if (conf.rpz * self.resofach) < dist and dabsH < dist:
            # Compute the resolution velocity vector in horizontal direction.
            # abs(tcpa) because it bcomes negative during intrusion.
            erratum=np.cos(np.arcsin((conf.rpz * self.resofach)/dist)-np.arcsin(dabsH/dist))
            # erratum_s=np.cos(np.arcsin((conf.rpz * self.resofach)/dist)-np.arcsin(dabsH_s/dist) + phi)

            dv1 = (((conf.rpz * self.resofach)/erratum - dabsH)*dcpa[0])/(abs(tcpa)*dabsH)
            dv2 = (((conf.rpz * self.resofach)/erratum - dabsH)*dcpa[1])/(abs(tcpa)*dabsH)

            # dv1_s = (((conf.rpz * self.resofach)/erratum_s - dabsH_s)*dcpa_s[0])/(abs(tcpa_s)*dabsH_s)
            # dv2_s = (((conf.rpz * self.resofach)/erratum_s - dabsH_s)*dcpa_s[1])/(abs(tcpa_s)*dabsH_s)
        else:
            dv1 = (iH * dcpa[0]) / (abs(tcpa) * dabsH)
            dv2 = (iH * dcpa[1]) / (abs(tcpa) * dabsH)

            # dv1_s = (iH_s * dcpa_s[0]) / (abs(tcpa_s) * dabsH_s)
            # dv2_s = (iH_s * dcpa_s[1]) / (abs(tcpa_s) * dabsH_s)

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

        ## Check if choosing vavg as solution leads to conflict ##

        
        vrel_vavg = np.array(vavg - v2)
        vrel = np.array(v1-v2)

        tcpa_vavg = -(vrel_vavg[0]*drel[0] + vrel_vavg[1]*drel[1])/(vrel_vavg[0]**2 + vrel_vavg[1]**2)
        dcpa_vavg = drel + vrel_vavg*tcpa_vavg
        if np.linalg.norm(-dcpa_vavg) > (conf.rpz * self.resofach):
            if np.sign(angle_clockwise(drel, vrel)) == np.sign(angle_clockwise(drel, vrel_vavg)):
                dv1 = vavg[0] - v1[0]
                dv2 = vavg[1] - v1[1]
            else:
                dv1 = -dv1
                dv2 = -dv2
        else:
            dv1 = -dv1
            dv2 = -dv2
        ##

        # Combine resolutions------------------------------------------------------

        # combine the dv components
        #dv = np.array([dv1,dv2,dv3])
        # if idx1 == 0:
        #     print(np.linalg.norm(v1[:-1] + np.array([dv1_s/2, dv2_s/2])), bs.traf.perf.show_vminmax()['vmax'][idx1]*0.9, bs.traf.perf.show_vminmax()['vmin'][idx1]*1.1)
        # if (np.linalg.norm(v1[:-1] + np.array([dv1_s/2, dv2_s/2])) < bs.traf.perf.show_vminmax()['vmax'][idx1]*0.9) and (np.linalg.norm(v1[:-1] + np.array([dv1_s/2, dv2_s/2])) > bs.traf.perf.show_vminmax()['vmin'][idx1]*1.1):
        #     dv = np.array([dv1_s,dv2_s,dv3])
        #     if idx1 == 0:
        #         print(idx1, 'spd', np.linalg.norm(v1[:-1] + np.array([dv1_s/2, dv2_s/2])), 'vmin', bs.traf.perf.show_vminmax()['vmin'][idx1], 'vmax', bs.traf.perf.show_vminmax()['vmax'][idx1])
        #         print(idx1, 'adapted')
        # else:
        #     if idx1 == 0:
        #         print(idx1, 'spd', np.linalg.norm(v1[:-1] + np.array([dv1_s/2, dv2_s/2])), 'vmin',bs.traf.perf.show_vminmax()['vmin'][idx1], 'vmax', bs.traf.perf.show_vminmax()['vmax'][idx1])
        #         print(idx1, 'original')

        dv = np.array([dv1,dv2,dv3])

        return dv, tsolV

    def resolve(self, conf, ownship, intruder):
        ''' Resolve all current conflicts '''
        # Initialize an array to store the resolution velocity vector for all A/C
        dv = np.zeros((ownship.ntraf, 3))

        # Initialize an array to store time needed to resolve vertically
        timesolveV = np.ones(ownship.ntraf) * 1e9

        # Call MVP function to resolve conflicts-----------------------------------
        for ((ac1, ac2), qdr, dist, tcpa, tLOS) in zip(conf.confpairs, conf.qdr, conf.dist, conf.tcpa, conf.tLOS):
            idx1 = ownship.id.index(ac1)
            idx2 = intruder.id.index(ac2)

            # If A/C indexes are found, then apply MVP on this conflict pair
            # Because ADSB is ON, this is done for each aircraft separately
            if idx1 >-1 and idx2 > -1:
                dv_mvp, tsolV = self.MVP(ownship, intruder, conf, qdr, dist, tcpa, tLOS, idx1, idx2)
                if idx1 == 5:
                    print('dv_mvp:',dv_mvp)
                if tsolV < timesolveV[idx1]:
                    timesolveV[idx1] = tsolV

                # Use priority rules if activated
                if self.swprio:
                    dv[idx1], _ = self.applyprio(dv_mvp, dv[idx1], dv[idx2], ownship.vs[idx1], intruder.vs[idx2])
                else:
                    # since cooperative, the vertical resolution component can be halved, and then dv_mvp can be added
                    dv_mvp[2] = 0.5 * dv_mvp[2]
                    dv[idx1] = dv[idx1] + dv_mvp
                    if idx1 == 5:
                        print('dv:',dv)

                # Check the noreso aircraft. Nobody avoids noreso aircraft.
                # But noreso aircraft will avoid other aircraft
                if self.noresoac[idx2]:
                    dv[idx1] = dv[idx1] + dv_mvp

                # Check the resooff aircraft. These aircraft will not do resolutions.
                if self.resooffac[idx1]:
                    dv[idx1] = 0.0


        # Determine new speed and limit resolution direction for all aicraft-------

        # Resolution vector for all aircraft, cartesian coordinates
        dv = np.transpose(dv)

        # The old speed vector, cartesian coordinates
        v = np.array([ownship.gseast, ownship.gsnorth, ownship.vs])

        # The new speed vector, cartesian coordinates
        newv = v + dv

        # Limit resolution direction if required-----------------------------------

        # Compute new speed vector in polar coordinates based on desired resolution
        if self.swresohoriz: # horizontal resolutions
            if self.swresospd and not self.swresohdg: # SPD only
                newtrack = ownship.trk
                newgs    = np.sqrt(newv[0,:]**2 + newv[1,:]**2)
                newvs    = ownship.vs
            elif self.swresohdg and not self.swresospd: # HDG only
                newtrack = (np.arctan2(newv[0,:],newv[1,:])*180/np.pi) % 360
                newgs    = ownship.gs
                newvs    = ownship.vs
            else: # SPD + HDG
                newtrack = (np.arctan2(newv[0,:],newv[1,:])*180/np.pi) %360
                newgs    = np.sqrt(newv[0,:]**2 + newv[1,:]**2)
                newvs    = ownship.vs
        elif self.swresovert: # vertical resolutions
            newtrack = ownship.trk
            newgs    = ownship.gs
            newvs    = newv[2,:]
        else: # horizontal + vertical
            newtrack = (np.arctan2(newv[0,:],newv[1,:])*180/np.pi) %360
            newgs    = np.sqrt(newv[0,:]**2 + newv[1,:]**2)
            newvs    = newv[2,:]

        # Determine ASAS module commands for all aircraft--------------------------

        # Cap the velocity
        newgscapped = np.maximum(ownship.perf.vmin,np.minimum(ownship.perf.vmax,newgs))

        # Cap the vertical speed
        vscapped = np.maximum(ownship.perf.vsmin,np.minimum(ownship.perf.vsmax,newvs))

        # Calculate if Autopilot selected altitude should be followed. This avoids ASAS from
        # climbing or descending longer than it needs to if the autopilot leveloff
        # altitude also resolves the conflict. Because asasalttemp is calculated using
        # the time to resolve, it may result in climbing or descending more than the selected
        # altitude.
        asasalttemp = vscapped * timesolveV + ownship.alt
        signdvs = np.sign(vscapped - ownship.ap.vs * np.sign(ownship.selalt - ownship.alt))
        signalt = np.sign(asasalttemp - ownship.selalt)
        alt = np.where(np.logical_or(signdvs == 0, signdvs == signalt), asasalttemp, ownship.selalt)

        # To compute asas alt, timesolveV is used. timesolveV is a really big value (1e9)
        # when there is no conflict. Therefore asas alt is only updated when its
        # value is less than the look-ahead time, because for those aircraft are in conflict
        altCondition = np.logical_and(timesolveV<conf.dtlookahead, np.abs(dv[2,:])>0.0)
        alt[altCondition] = asasalttemp[altCondition]

        # If resolutions are limited in the horizontal direction, then asasalt should
        # be equal to auto pilot alt (aalt). This is to prevent a new asasalt being computed
        # using the auto pilot vertical speed (ownship.avs) using the code in line 106 (asasalttemp) when only
        # horizontal resolutions are allowed.
        alt = alt * (1 - self.swresohoriz) + ownship.selalt * self.swresohoriz
        return newtrack, newgscapped, vscapped, alt