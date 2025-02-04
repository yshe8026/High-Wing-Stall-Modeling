# XPlanePlugins

This repository contains a collection of plugins for the X-Plane flight simulator for custom Physics models and autopilot controllers. The plugins provide a easy interface to implement custom physics models and autopilot controllers for aircrafts in X-Plane. Currently, the physics models are based on linearised flight derivatives and the autopilot controllers are based on PID controllers. The plugins are written in Python and are compatible with X-Plane 11 (tested) and X-Plane 12 (untested).

## Installation

### X-Plane 11 
1. Install X-Plane 11 from [here](https://www.x-plane.com/desktop/try-it/older/)

2. Install Python (3.10.x) from [here](https://www.python.org/downloads/release/python-31013/).  It is important to install the correct version of python for the plugin to work correctly.
    + During the installation step:
        - Ensure you tick the box saying 'Add python.exe to path'
        - Ensure you tick the box saying 'Install with admin privelages'

    + During post installation, click on 'disable path length restriction'.

3. Download and install the XPPython plugin from [here](https://xppython3.readthedocs.io/en/3.1.5/usage/installation_plugin.html). Please follow the instructions on the website to install the plugin correctly. Make sure the plugin loads and works as per documentation, before proceeding to the next step.

4. Download and Copy the plugin repository (request access permission by contacting authors) to 'PythonPlugins' folder in the X-Plane 11 directory. The folder structure should look like this:

```
X-Plane 11
├── Resources
│   ├── plugins
│   │   ├── PythonPlugins
│   │   │   ├── PI_custom_physics.py
            |── ...
```

5. Install the required Python packages. Execute the following command in the terminal to install the required packages:

```
python.exe -m pip install --upgrade pip
pip install -r requirements.txt
```

Restart the PC after installing the requirements.

If there is an error installing Numba, one might need to update the Microsoft C and C++ (MSVC) runtime libraries.  To do this, download and execute the architecture specific executable from [here](https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170)

6. Install X-PlaneConnect. Download the latest version of X-PlaneConnect from [here](
    https://github.com/nasa/XPlaneConnect/releases). Copy the X-PlaneConnect plugin to the 'Resources/plugins' folder in the X-Plane 11 directory. The folder structure should look like this:

```
X-Plane 11
├── Resources
│   ├── plugins
│   │   ├── XPlaneConnect
│   │   │   ├── 64
...
```

7. Copy any custom aircraft files to the 'Aircraft' folder in the X-Plane 11 directory.
```
X-Plane 11
├── Aircraft
│   ├── Custom Aircraft Directory
```

8. (Optional) Install Visual Studio Code from [here](https://code.visualstudio.com/download). This is not a critical step, but is a highly recomended code editor.

9. Start X-Plane 11 and load the aircraft. The plugin should load automatically. Once the plane is loaded, the custom physics model can be enabled by pressing the Key 9 on the keyboard. This should open up a dialog box on the bottom right corner saying "Custom Physics Activated". You should also hear a voice saying "Custom Physics Activated". If you do not hear the voice, please check the log file (X-Plane 11/XPPython3.log) for any errors. If you do not see a log file, please check the installation of the XPPython plugin.

10. Additional Hardware.  The plugin currently supports the 'Logitech G Flight Simulator Autopilot Multipanel'.  To set up the autopilot multipanel:
    1. Install Xsaitekpanels. Download the latest version of Xsaitekpanels from [here](https://forums.x-plane.org/index.php?/files/file/14646-xsaitekpanels-linwinmac/). Copy the Xsaitekpanels plugin to the 'Resources/plugins' folder in the X-Plane 11 directory. The folder structure should look like this:

    ```
    X-Plane 11
    ├── Resources
    │   ├── plugins
    │   │   ├── Xsaitekpanels
    ...
    ```

    2. Download and install the multipanel driver from [here](https://support.logi.com/hc/en-au/articles/360024692354--Downloads-Flight-Multi-Panel)

    The panel should now work assuming correct installation of xsaitekpanels in the X-Plane -> Resources -> plugins folder. For Windows users: if the panel powers on momentarily, before powering off, the power saving mode is interfering with the USB port. To over come this, open ‘Device Manager’. For each ‘Generic USB Hub’ and/or ‘USB Root HUB (USB 3.0)’ under the drop down menu for ‘Universal Serial Bus controllers’:
        i.   Right click on ‘Generic USB Hub’ or ‘USB Root HUB (USB 3.0)’
        ii.  Navigate: Properties -> Power Management
        iii. Untick: ‘Allow the computer to turn off this device to save power’

    3.  The autopilot has a yaw damper control loop.  The 'Logitech G Flight Simulator Autopilot Multipanel' does not have a yaw damper switch, thus you have to manually add one by assigning an X-Plane Command to an unused button (i.e. not on the Logitech G Flight Simulator Autopilot Multipanel).  To do this:
        i.   Open X-Plane Settings.  
        ii.  Locate Hardware Calibration.
        iii. Set a Yaw Damper On/Off Switch.

### X-Plane 12

For X-Plane 12, the installation process is similar to X-Plane 11. The only difference is that the a different version of XPPython plugin is required. The XPPython plugin for X-Plane 12 can be downloaded from [here](https://xppython3.readthedocs.io/en/latest/usage/installation_plugin.html). It also requires a different Python version. Please follow the instructions on the website to install the plugin correctly. As of today (5/10/2023), Numba is not supported for Python 3.12. This might change in the future. 

NOTE: The plugin is currently designed for X-Plane 11.  There are some known differences with the X-Plane 12 SDK, which might result in undefined plugin behaviour.  In future, the plugin will have a branch specific to X-Plane 12 to overcome this issue. 


## Authors
- [Brendan Waters](brendan.waters@sydney.edu.au)
- [Ankith Anil Das](ankith.anildas@sydney.edu.au)


