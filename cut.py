import numpy as np
from math import pi, tan, cos, sin, sqrt
import sys
import argparse
render = True
try:
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
except:
    render = False

parser = argparse.ArgumentParser(description=\
    "Calculates cutting path around cylinder for certain angled cut.")
parser.add_argument("-r", "--radius", type=int, default=100,
    help="Radius of the cylinder in mm")
parser.add_argument("-a", "--angle", type=int, default=45,
    help="Angle of the cut in degrees from the cylinder axis")
parser.add_argument("-i", "--interval", type=float, default=.5,
    help="Cylinder intersection interval in proportion of circumference (0.0-1.0)")
parser.add_argument('--display', dest='display', action='store_true',
    help="Render cut")
parser.add_argument('--no-display', dest='display', action='store_false',
    help="Do not render cut")
parser.set_defaults(display=True)
parser.add_argument("-f", "--file", type=str, default='cut.csv',
    help="CSV file to write into cut mark positions (around cylinder and along)")
args = parser.parse_args()

radius = args.radius
assert radius > 15, "Radius must be positive and in mm."
angle = args.angle
assert 90 >= angle > 0, "Angle must be between 0 and 90 degrees."
angle = (angle * pi) / 180
interval = args.interval
assert 0.25 >= interval >= 0.005, "Interval must be <= 0.25 and >= 0.005"
render = render and args.display
filename = args.file
assert len(filename) > 0, "Filename must be at least one character long"

circumference = (int)(radius * 2 * pi)
interval = circumference * interval
cyl_length = 2 * radius / tan(angle)
cut_length = 2 * radius / sin(angle)

print("Calculating {0} degree cut of {1}mm radius cylinder. "
      "Approximating at {2} mm arc intervals".format(args.angle, radius, interval))

def rotation_matrix(axis,theta):
    '''Create a rotation matrix for a theta radian
       rotation around the axis given.'''
    axis = axis/sqrt(np.dot(axis,axis))
    a = cos(theta/2)
    b,c,d = -axis*sin(theta/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def vertical_plane_normal(p1,p2):
    '''Compute a normal to the cutting plane'''
    p3 = p1 + [0,0,1]
    return np.cross(p2-p1,p3-p1)		

# Approximate cylinder with parallel lines
lines_cylinder = []
for i in range(0,int(circumference/interval)):
    ''' Builds a cylinder intersection approximated as a set of parallel lines
        each separated by an arc of length 'interval' '''
    theta = (2 * pi) * (i / (circumference/interval))
    rotmat = rotation_matrix(np.array([0, 1, 0]), -theta)
    lines_cylinder.append(np.dot(rotmat, np.array([0, -cyl_length/2, radius])))
    lines_cylinder.append(np.dot(rotmat, np.array([0, cyl_length/2, radius])))

# Create cutting plane (a line will do for now)
rotmat = rotation_matrix(np.array([0,0,1]),angle)
cutting_line_st = np.dot(rotmat, np.array([0, -cut_length/2, 0]))
cutting_line_end = np.dot(rotmat, np.array([0, cut_length/2, 0]))

# Calculate cutting plane/cylinder intersection points.
# Only computes the first 180 degrees as the other 180
# is just a mirror of it.
ixs = []
for i in range(0, len(lines_cylinder), 2):
    N = np.array(vertical_plane_normal(lines_cylinder[i], lines_cylinder[i+1]))
    ix = cutting_line_st + (np.dot(N, lines_cylinder[i] - cutting_line_st) /
             np.dot(N, cutting_line_end - cutting_line_st)) * (cutting_line_end - cutting_line_st)
    ix = [lines_cylinder[i][0], ix[1], lines_cylinder[i][2]];
    ixs.append(ix)

# Flatten cylinder intersections to give cuts on a 2D plane.
# These can be applied to the real cylinder by wrapping
# this 2D plane around the cylinder. The best way to do this
# is either by printing (to correct scale) the markers and
# wrapping the 2D paper around the cylinder or drawing these
# marks on graph paper and wrapping this around the cylinder. 
ixs_flat = []
for i in range(int(len(ixs)/2)):
    point = [i * interval, ixs[i][1]]
    ixs_flat.append(point)
for i in range(int(len(ixs)/2)):
    point = [circumference/2 + (i * interval), - ixs[i][1]]
    ixs_flat.append(point)
	
f4 = np.poly1d(np.polyfit([ix[0] for ix in ixs_flat] , [ix[1] for ix in ixs_flat], 8))
xp = np.linspace(0, circumference, 100)

if render:
    # Render 3D cut
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.axis('equal')

    for i in range(0, len(lines_cylinder), 2):
	l_st = lines_cylinder[i]
	l_en = lines_cylinder[i+1]
	ax.plot([l_st[0],l_en[0]], [l_st[1],l_en[1]],zs=[l_st[2],l_en[2]])

    ax.plot([cutting_line_st[0], cutting_line_end[0]], [cutting_line_st[1], cutting_line_end[1]], [cutting_line_st[2],cutting_line_end[2]])
    ax.scatter([ix[0] for ix in ixs], [ix[1] for ix in ixs], zs=[ix[2] for ix in ixs])
    ax.set_ylabel('Cylinder Axis (mm)')
    ax.set_xlabel('mm')
    ax.set_zlabel('mm')
    plt.show()

    # Render cut marker positions
    fig = plt.plot([ix[0] for ix in ixs_flat], [ix[1] for ix in ixs_flat], '.', xp, f4(xp), '-')
    plt.ylim(min([ix[1] for ix in ixs_flat]), max([ix[1] for ix in ixs_flat]))
    plt.xlabel('Around the Cylinder (mm)')
    plt.ylabel('Along the Cylinder (mm)')
    plt.title('Unwrapped cylinder cut marker positions (printed and wrapped around cylinder).')
    plt.axis('equal')
    plt.show()

# Write cut markers to file
print("Writing cut marker positions to {}".format(filename))
file = open(filename, 'w')
file.write("arc pos (mm), length pos (mm)")
for ix in ixs_flat:
    file.write("{0[0]:.3f}, {0[1]:.3f}\n".format(ix))
file.close()
print("Finished writing to file")
