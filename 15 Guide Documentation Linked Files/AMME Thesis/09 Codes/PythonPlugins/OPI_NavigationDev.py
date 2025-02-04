"""
Navigation.py

Ported to Python by Sandy Barbour - 28/04/2005
Ported to XPPython3 by Peter Buckner - 2-Aug-2020

This example demonstrates how to use the FMC and the navigation databases in
X-Plane.
"""
import xp
import numpy as np
from Autopilot.FlightDirector import FlightDirector, NavigationFunctions
nearestAirport = 1
programFMC = 2


class PythonInterface:
    def __init__(self):
        self.TargetAirport = "YSSY"
        self.FlightDirector = FlightDirector()
        self.Navigate = NavigationFunctions()

    def XPluginStart(self):
        self.Name = "Test Navigation1 v1.0"
        self.Sig = "navigation1.demos.xppython3"
        self.Desc = "A plugin that controls the FMC."

        
        mySubMenuItem = xp.appendMenuItem(xp.findPluginsMenu(), "Python - Navigation 1", 0)
        self.myMenu = xp.createMenu("Navigation1", xp.findPluginsMenu(), mySubMenuItem, self.MyMenuHandlerCallback, 0)
        xp.appendMenuItem(self.myMenu, "Program FMC", programFMC)
        xp.appendMenuItem(self.myMenu, "Say next waypoint", nearestAirport)


        return self.Name, self.Sig, self.Desc

    def XPluginStop(self):
        xp.destroyMenu(self.myMenu)

    def XPluginEnable(self):
        return 1

    def XPluginDisable(self):
        pass

    def XPluginReceiveMessage(self, inFromWho, inMessage, inParam):
        pass

    def MyMenuHandlerCallback(self, inMenuRef, inItemRef):
       
        if (inItemRef == programFMC):
            self.FlightDirector.read_waypoints_from_file()

        if (inItemRef == nearestAirport):
            target_lla = self.FlightDirector.current_waypoint()
            
            if target_lla is not None:
                lat = np.deg2rad(xp.getDataf(xp.findDataRef("sim/flightmodel/position/latitude")))  
                lon = np.deg2rad(xp.getDataf(xp.findDataRef("sim/flightmodel/position/longitude")))
                alt = xp.getDataf(xp.findDataRef("sim/flightmodel/position/y_agl"))


                hdgc, altc = self.Navigate.navigation_commanded_input([lat,lon,alt], target_lla)

                buffer = "The next waypoint is located at, Lat %.2f degrees, Lon %.2f degrees.  Commanded heading: %.2f degrees, Commanded altitude %.2f feet" % (np.rad2deg(target_lla[0]), np.rad2deg(target_lla[1]), np.rad2deg(hdgc), target_lla[2])
                xp.speakString(buffer)

                target_achieved = self.Navigate.target_achieved([lat,lon,alt], target_lla)

                if target_achieved:
                    self.FlightDirector.delete_waypoint()



        # if (inItemRef == nearestAirport):
        #     # First find the plane's position.
        #     lat = xp.getDataf(xp.findDataRef("sim/flightmodel/position/latitude"))
        #     lon = xp.getDataf(xp.findDataRef("sim/flightmodel/position/longitude"))
        #     # Find the nearest airport to us.
        #     ref = xp.findNavAid(None, self.TargetAirport , lat, lon, None, xp.Nav_Airport)
        #     if (ref != xp.NAV_NOT_FOUND):
        #         navAidInfo = xp.getNavAidInfo(ref)
        #         buf = "The nearest airport is %s, %s, Lat %f, Lon %f" % (navAidInfo.navAidID, navAidInfo.name, navAidInfo.latitude, navAidInfo.longitude)
        #         xp.speakString(buf)
        #         print(buf)
        #     else:
        #         xp.speakString("No airports were found!")
        #         print("No airports were found!")

        # if (inItemRef == programFMC):
        #     # This code programs the flight management computer.  We simply set each entry to a navaid
        #     # that we find by searching by name or ID.
        #     xp.setFMSEntryInfo(0, xp.findNavAid(None, "KBOS", None, None, None, xp.Nav_Airport), 3000)
        #     xp.setFMSEntryInfo(1, xp.findNavAid(None, "LUCOS", None, None, None, xp.Nav_Fix), 20000)
        #     xp.setFMSEntryInfo(3, xp.findNavAid(None, "PARCH", None, None, None, xp.Nav_Fix), 20000)
        #     xp.setFMSEntryInfo(5, xp.findNavAid(None, "ROBER", None, None, None, xp.Nav_Fix), 9000)
        #     xp.setFMSEntryInfo(6, xp.findNavAid(None, "KJFK", None, None, None, xp.Nav_Airport), 3000)
        #     xp.clearFMSEntry(7)
        #     xp.clearFMSEntry(8)



















# class PythonInterface:
#     def __init__(self):
#         self.TargetAirport = "YSSY"

#     def XPluginStart(self):
#         self.Name = "Navigation1 v1.0"
#         self.Sig = "navigation1.demos.xppython3"
#         self.Desc = "A plugin that controls the FMC."
#         mySubMenuItem = xp.appendMenuItem(xp.findPluginsMenu(), "Python - Navigation 1", 0)
#         self.myMenu = xp.createMenu("Navigation1", xp.findPluginsMenu(), mySubMenuItem, self.MyMenuHandlerCallback, 0)
#         xp.appendMenuItem(self.myMenu, "Say nearest airport", nearestAirport)
#         xp.appendMenuItem(self.myMenu, "Program FMC", programFMC)
#         return self.Name, self.Sig, self.Desc

#     def XPluginStop(self):
#         xp.destroyMenu(self.myMenu)

#     def XPluginEnable(self):
#         return 1

#     def XPluginDisable(self):
#         pass

#     def XPluginReceiveMessage(self, inFromWho, inMessage, inParam):
#         pass

#     def MyMenuHandlerCallback(self, inMenuRef, inItemRef):
#         if (inItemRef == nearestAirport):
#             # First find the plane's position.
#             lat = xp.getDataf(xp.findDataRef("sim/flightmodel/position/latitude"))
#             lon = xp.getDataf(xp.findDataRef("sim/flightmodel/position/longitude"))
#             # Find the nearest airport to us.
#             ref = xp.findNavAid(None, self.TargetAirport , lat, lon, None, xp.Nav_Airport)
#             if (ref != xp.NAV_NOT_FOUND):
#                 navAidInfo = xp.getNavAidInfo(ref)
#                 buf = "The nearest airport is %s, %s, Lat %f, Lon %f" % (navAidInfo.navAidID, navAidInfo.name, navAidInfo.latitude, navAidInfo.longitude)
#                 xp.speakString(buf)
#                 print(buf)
#             else:
#                 xp.speakString("No airports were found!")
#                 print("No airports were found!")

#         if (inItemRef == programFMC):
#             # This code programs the flight management computer.  We simply set each entry to a navaid
#             # that we find by searching by name or ID.
#             xp.setFMSEntryInfo(0, xp.findNavAid(None, "KBOS", None, None, None, xp.Nav_Airport), 3000)
#             xp.setFMSEntryInfo(1, xp.findNavAid(None, "LUCOS", None, None, None, xp.Nav_Fix), 20000)
#             xp.setFMSEntryInfo(3, xp.findNavAid(None, "PARCH", None, None, None, xp.Nav_Fix), 20000)
#             xp.setFMSEntryInfo(5, xp.findNavAid(None, "ROBER", None, None, None, xp.Nav_Fix), 9000)
#             xp.setFMSEntryInfo(6, xp.findNavAid(None, "KJFK", None, None, None, xp.Nav_Airport), 3000)
#             xp.clearFMSEntry(7)
#             xp.clearFMSEntry(8)

#             # general form of xp,findNavAid:
#             # XPLM_API XPLMNavRef XPLMFindNavAid(
#             #              const char *         inNameFragment,    /* Can be NULL */
#             #              const char *         inIDFragment,    /* Can be NULL */
#             #              float *              inLat,    /* Can be NULL */
#             #              float *              inLon,    /* Can be NULL */
#             #              int *                inFrequency,    /* Can be NULL */
#             #              XPLMNavType          inType);
