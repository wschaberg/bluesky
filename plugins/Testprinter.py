""" BlueSky plugin template. The text you put here will be visible
    in BlueSky as the description of your plugin. """
# Import the global bluesky objects. Uncomment the ones you need
import bluesky as bs
from bluesky import stack, traffic, traf  #, settings, navdb, traf, sim, scr, tools
from bluesky.traffic.asas import ConflictDetection
from bluesky.traffic.asas import ConflictResolution
from bluesky.traffic.asas import MVP

### Initialization function of your plugin. Do not change the name of this
### function, as it is the way BlueSky recognises this file as a plugin.
def init_plugin():

    # Addtional initilisation code
    stack.stack('ASAS ON')
    stack.stack('PLUGIN MVPGON')
    stack.stack('RESO MVPGON')
    stack.stack('CRE 1 B737 52 2 0 FL100 280')
    stack.stack('CRECONFS 2 B737 1 30 4.5 302')
    # Configuration parameters
    config = {
        # The name of your plugin
        'plugin_name':     'Testprinter',

        # The type of this plugin. For now, only simulation plugins are possible.
        'plugin_type':     'sim',

        # Update interval in seconds. By default, your plugin's update function(s)
        # are called every timestep of the simulation. If your plugin needs less
        # frequent updates provide an update interval.
        'update_interval': 2.5,

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
        'MYFUN': [
            # A short usage string. This will be printed if you type HELP <name> in the BlueSky console
            'MYFUN ON/OFF',

            # A list of the argument types your function accepts. For a description of this, see ...
            '[onoff]',

            # The name of your function in this plugin
            myfun,

            # a longer help text of your function.
            'Print something to the bluesky console based on the flag passed to MYFUN.']
    }

    # init_plugin() should always return these two dicts.
    return config, stackfunctions

### Periodic update functions that are called by the simulation. You can replace
### this by anything, so long as you communicate this in init_plugin

def update():
    pass

def preupdate():
    pass

def reset():
    pass

### Other functions of your plugin
def myfun(flag=True):
    return True, 'My plugin received an o%s flag.' % ('n' if flag else 'ff')
