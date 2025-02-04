import numpy as np

"""
======================================================================================================================
    \brief   Computes the transformation matrix from the local OpenGL to standard Aircraft Body Coordernates
                +x-axis points to the aircrafts nose
                +y-axis points to the right wing
                +z-axis points down

            DO NOT CHANGE THIS FUNCTION UNLESS YOU KNOW WHAT YOU ARE DOING!
            CHANGE THE FUNCTION CAN RESULT IN UNEXPECTED BEHAVIOUR OF THE AIRCRAFT!
======================================================================================================================
"""
def LocalToBody(xyz_acft : np.ndarray, angles : np.ndarray):
    '''
    Converts from local OpenGL frame coordinates to standard aircraft body coordinates

    Parameters:
    ------------
    xyz_acft: 3x1 numpy array of floats 
        Vector in local OpenGL frame coordinates
    angles: 3x1 numpy array of floats
        Euler angles in radians
    
    Returns:
    ------------
    xyz_body: 3x1 numpy array of floats
        Vector in standard aircraft body coordinates

    Here the first 3 transformations are standard rotation matrices and the last one changes the axes from OpenGL to
    standard aircraft body coordinates. 
    '''

    phi = angles[0]
    theta = angles[1]
    psi = angles[2]

    R_z = np.array([[np.cos(phi), np.sin(phi), 0],
                    [-np.sin(phi), np.cos(phi), 0],
                    [0, 0, 1]])
    R_x = np.array([[1, 0,0],
                    [0, np.cos(theta), -np.sin(theta)],
                    [0, np.sin(theta), np.cos(theta)]])
    R_y = np.array([[np.cos(psi), 0, -np.sin(psi)],
                    [0, 1, 0],
                    [np.sin(psi), 0, np.cos(psi)]])
    R_t = np.array([[0,0,-1], 
                [1, 0, 0],
                [0, -1, 0]])

    R = R_y @ R_x @ R_z 
    R = R_t @ R.T
    return R @ xyz_acft

"""
======================================================================================================================
    \brief   Computes the transformation matrix standard Aircraft Body Coordernates to the local OpenGL frame 
             of reference.

             Inverse of LocalToBody()

             DO NOT CHANGE THIS FUNCTION UNLESS YOU KNOW WHAT YOU ARE DOING!
             CHANGE THE FUNCTION CAN RESULT IN UNEXPECTED BEHAVIOUR OF THE AIRCRAFT!
======================================================================================================================
"""
def BodyToLocal(xyz_body : np.ndarray, angles : np.ndarray):
    '''
    Converts from standard aircraft body coordinates to local OpenGL frame coordinates. 
    It is the inverse of LocalToBody()
    
    Parameters:
    ------------
    xyz_body: 3x1 numpy array of floats
        Vector in standard aircraft body coordinates
    angles: 3x1 numpy array of floats
        Euler angles in radians
    
    Returns:
    ------------
    xyz_acft: 3x1 numpy array of floats
        Vector in local OpenGL frame coordinates
    
    '''

    phi = angles[0]
    theta = angles[1]
    psi = angles[2]

    R_z = np.array([[np.cos(phi), np.sin(phi), 0],
                    [-np.sin(phi), np.cos(phi), 0],
                    [0, 0, 1]])
    R_x = np.array([[1, 0,0],
                    [0, np.cos(theta), -np.sin(theta)],
                    [0, np.sin(theta), np.cos(theta)]])
    R_y = np.array([[np.cos(psi), 0, -np.sin(psi)],
                    [0, 1, 0],
                    [np.sin(psi), 0, np.cos(psi)]])
    R_t = np.array([[0,0,-1], 
                [1, 0, 0],
                [0, -1, 0]])

    R = R_y @ R_x @ R_z 
    R = R @ R_t.T
    return R @ xyz_body


"""
======================================================================================================================
    \brief   The PythonInterface class is the plug-in interface to X-Plane 11.

             This script is executed within the X-Plane 11 kernel itself, thus enabling direct access to X-Plane's
             data structures (DRef's). These Data References can be read-only, write-only, or read-and-write.

             Writing to a dataref directly changes the aircraft's state or function.

             The available data references can be found in X-Planes documentation:
                - https://developer.x-plane.com/datarefs/    
======================================================================================================================
"""