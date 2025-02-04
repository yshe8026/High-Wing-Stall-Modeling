from XPPython3 import xp


class UIElements:
    def __init__(self):
        # General Parameters:
        self.message = None
        # Message Parameters:
        self.message_start      = None
        self.message_duration   = 5.0
        self.windowID     = None
        self.message_loop = None
        self.menuEntry = {}

    def ShowMessage(self, message):
        if self.message_loop is not None:
            self.DestroyWindow()
            xp.destroyFlightLoop(self.message_loop)
            self.message_loop = None

        self.message = message
        xp.speakString(self.message)
        # This is a bit of a hack to get the message length in pixels
        message_length = len(self.message) * 8
        if self.windowID is None: 
            self.windowID = xp.createWindowEx(self.GetWindow(25, 75, message_length, 25))
            self.message_start = xp.getElapsedTime()
            self.message_loop = xp.createFlightLoop(self.message_distroy_loop, 
                                                xp.FlightLoop_Phase_AfterFlightModel)
            xp.scheduleFlightLoop(self.message_loop, 0.1, 1)

    def message_expired(self, elapsedTime):
        return elapsedTime - self.message_start > self.message_duration
    
    def createMenu(self, menuName, handler, parentMenu=None):
        if parentMenu is None:
            parentMenu = xp.findPluginsMenu()
        self.menuEntry[menuName] = xp.appendMenuItem(parentMenu, menuName, 0)
        self.menuHandler = xp.createMenu(menuName, parentMenu, self.menuEntry[menuName], handler, 0)
        return self.menuHandler
    
    def attachWindowToDestroy(self, windowID, duration=5.0):
        if self.message_loop is not None:
            self.DestroyWindow()
            xp.destroyFlightLoop(self.message_loop)
            self.message_loop = None
        self.windowID = windowID
        self.message_duration = duration
        self.message_start = xp.getElapsedTime()
        self.message_loop = xp.createFlightLoop(self.message_distroy_loop, 
                                            xp.FlightLoop_Phase_AfterFlightModel)
        xp.scheduleFlightLoop(self.message_loop, 0.1, 1)

    def DestroyWindow(self):
        '''
        Distroy current window
        '''
        if self.windowID is not None:
            xp.destroyWindow(self.windowID)
            self.windowID = None
            self.message = ""
 
    def GetWindow(self, xmin, ymax, xmax, ymin):
        return (xmin, ymax, xmax, ymin,
                1,
                self.DrawWindowCallback,
                self.MouseClickCallback,
                self.KeyCallback,
                self.CursorCallback,
                self.MouseWheelCallback,
                0,
                xp.WindowDecorationRoundRectangle,
                xp.WindowLayerFloatingWindows,
                None)
     
    """
    This callback draws the pop-up window with custom text once per sim cycle 
    until cancelled. 
    """
    def DrawWindowCallback(self, inWindowID, inRefcon):
        # First we get the location of the window passed in to us.
        (left, top, right, bottom) = xp.getWindowGeometry(inWindowID)
        
        xp.drawTranslucentDarkBox(left, top, right, bottom)
        color = 1.0, 1.0, 1.0

        xp.drawString(color, left + 5, top - 20, self.message, 0, xp.Font_Basic)

    """
    Keyboard input callback
    """
    def KeyCallback(self, inWindowID, inKey, inFlags, inVirtualKey, inRefcon, losingFocus):
        pass

    """
    Mouse click callback
    """
    def MouseClickCallback(self, inWindowID, x, y, inMouse, inRefcon):
        pass

    """
    Cursor status callback.
    """
    def CursorCallback(self, inWindowID, x, y, inRefcon):
        return xp.CursorDefault

    """
    Mouse wheel callback; return 1 stop mouse wheel movement, 0 to pass them to lower window.
    """
    def MouseWheelCallback(self, inWindowID, x, y, wheel, clicks, inRefcon):
        return 1
    
    """
    ======================================================================================================================
        \brief  Destroys pop-up message for custom physics on/off after x-seconds.
    ======================================================================================================================
    """
    def message_distroy_loop(self, sinceLast, elapsedTime, counter, refcon):
        # Store the time at first exection
        if self.message_expired(elapsedTime):
            # Distroy window
            self.DestroyWindow()
            # set duration back to default
            self.message_duration = 5.0
            return 0
        return 1