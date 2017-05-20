# -----------------------------------------------------------------------------
# Class to handle pdb transformations
#

# -----------------------------------------------------------------------------
# Apply a rotation and translation to atoms.
#
def transform_atom_coordinates_adp(atoms, xform):

    for a in atoms:
        a.setCoord(xform.apply(a.coord()))

# -----------------------------------------------------------------------------
# Apply a rotation and translation to model coordinate axes.
#
def transform_coordinate_axes_adp(model, xform):

    model.openState.localXform(xform)

# -----------------------------------------------------------------------------
# Rotation applied first, then translation.
#
def euler_xform_adp(euler_angles, translation):

    xf = euler_rotation_adp(*euler_angles)
    from chimera import Xform
    xf.premultiply(Xform.translation(*translation))
    return xf

# -----------------------------------------------------------------------------
# Convert Euler angles to an equivalent Chimera transformation matrix.
#
# Angles must be in degrees, not radians.
#
# Uses the most common Euler angle convention z-x-z (the chi-convention)
# described at
#
#   http://mathworld.wolfram.com/EulerAngles.html
#
def euler_rotation_adp(phi, theta, psi):

    from chimera import Xform, Vector
    xf1 = Xform.rotation(Vector(0,0,1), phi)    # Rotate about z-axis
    xp = xf1.apply(Vector(1,0,0))               # New x-axis
    xf2 = Xform.rotation(xp, theta)             # Rotate about new x-axis
    zp = xf2.apply(Vector(0,0,1))               # New z-axis
    xf3 = Xform.rotation(zp, psi)               # Rotate about new z-axis

    xf = Xform()
    xf.premultiply(xf1)
    xf.premultiply(xf2)
    xf.premultiply(xf3)

    return xf
