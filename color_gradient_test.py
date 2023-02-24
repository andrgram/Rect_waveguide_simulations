import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.widgets import Slider, RadioButtons, CheckButtons



# TODO Runtime optimalizations....
# TODO Examine quiver plots and streamplots


# Number of sample points
points = 10
x_points = points
y_points = points
z_points = int(points * 2 * np.pi)
t_points = int(points * 2 * np.pi)

x_range = np.linspace(0, 1, x_points)
y_range = np.linspace(0, 1, y_points)
z_range = np.linspace(0, 2*np.pi, z_points)
t_range = np.linspace(0, 2*np.pi, t_points)

XYplane = np.meshgrid(x_range, y_range)
XZplane = np.meshgrid(z_range, x_range)
YZplane = np.meshgrid(z_range, y_range)


# ********************************************
# ****** Setup figure and axes ***************
# ********************************************
fig, ((axXY, axYZ),(axEmpty, axXZ)) = plt.subplots(2,2, width_ratios=[1,3], figsize=(16,8))
fig.subplots_adjust(bottom=0.25)
fig.delaxes(axEmpty)



# ********************************************
# ****** Setup for the sliders and buttons ***
# ********************************************
axXcoord = fig.add_axes([0.15, 0.15, 0.7, 0.05]) # x,y offset for the center of the slider, width and height.
xCoord_slider = Slider(
    ax=axXcoord,
    label = 'x-coordinate',
    valmin=0,
    valmax=1,
    valinit=0.5,
    valstep=1/x_points
)
axYcoord = fig.add_axes([0.15, 0.1, 0.7, 0.05])
yCoord_slider = Slider(axYcoord, 'y-coordinate', 0, 1, valinit=0.5, valstep=1/y_points)
axZcoord = fig.add_axes([0.15, 0.05, 0.7, 0.05])
zCoord_slider = Slider(axZcoord, 'z-coordinate', 0, 2*np.pi, valinit=0, valstep=2*np.pi/z_points)
axTcoord = fig.add_axes([0.15, 0, 0.7, 0.05])
tCoord_slider = Slider(axTcoord, 'time', 0, 2*np.pi, valinit=0, valstep=2*np.pi/t_points)
axMmode = fig.add_axes([axXY.get_position().x0, axXY.get_position().y0 - 0.1, axXY.get_position().width, 0.05])
mMode_slider = Slider(axMmode, 'modes in x-direction', 0, 4, valinit=1, valstep=1)
axNmode = fig.add_axes([axXY.get_position().x0, axXY.get_position().y0 - 0.15, axXY.get_position().width, 0.05])
nMode_slider = Slider(axNmode, 'modes in y-direction', 0, 4, valinit=0, valstep=1)
axFieldButton = fig.add_axes([axXY.get_position().x0 - 0.05, axXY.get_position().y0 - 0.3, 0.05, 0.1])
fieldButton = RadioButtons(axFieldButton, ('E-field', 'H-field'))
axPlotButton = fig.add_axes([axXY.get_position().x0 + 0.01, axXY.get_position().y0 - 0.3, 0.07, 0.1])
plotButton = RadioButtons(axPlotButton, ('Arrow plot', 'Streamplot', 'No Arrows'))
axBackgButton = fig.add_axes([axXY.get_position().x0 + 0.09, axXY.get_position().y0 - 0.3, 0.11, 0.1])
backgButton = CheckButtons(axBackgButton, ['Background fields', 'Crosshairs'], actives=[True, True])



# ********************************************
# ****** Functions for the fields ************
# ********************************************
# Field and wave Electromagnetics p.553, TE propagation in rectangular waveguides
def H(x,y,m=1,n=0):
    H_x = 1j * m * np.sin(m * np.pi * x) * np.cos(n * np.pi * y)
    H_y = 1j * n * np.cos(m * np.pi * x) * np.sin(n * np.pi * y)
    H_z = np.cos(m * np.pi * x) * np.cos(n * np.pi * y)
    return [H_x, H_y, H_z]

def H_ins(x,y,z,t,m=1,n=0):
    H_plane = H(x,y,m,n)
    H_x = np.real(H_plane[0] * np.e ** ((t-z)*1j))
    H_y = np.real(H_plane[1] * np.e ** ((t-z)*1j))
    H_z = np.real(H_plane[2] * np.e ** ((t-z)*1j))
    return [H_x, H_y, H_z]

def E(x,y,m=1,n=0):
    E_x =  1j * n * np.cos(m * np.pi * x) * np.sin(n * np.pi * y)
    E_y = -1j * m * np.sin(m * np.pi * x) * np.cos(n * np.pi * y)
    E_z =  0
    return [E_x, E_y, E_z]

def E_ins(x,y,z,t,m=1,n=0):
    E_plane = E(x,y,m,n)
    E_x = np.real(E_plane[0] * np.e ** ((t-z)*1j))
    E_y = np.real(E_plane[1] * np.e ** ((t-z)*1j))
    E_z = np.real(E_plane[2] * np.e ** ((t-z)*1j))
    return [E_x, E_y, E_z]




# ***************************************************************
# ****** Calculate the fields in the planes examined ************
# ***************************************************************
def calculateFields():
    # Calculating the fields in either E or H
    if fieldButton.value_selected == 'E-field':
        XYplaneField = np.array([[[(E_mag := E_ins(x,y,zCoord_slider.val,tCoord_slider.val,mMode_slider.val,nMode_slider.val))[0], E_mag[1], E_mag[2]] for x in x_range] for y in y_range])
        XZplaneField = np.array([[[(E_mag := E_ins(x,yCoord_slider.val,z,tCoord_slider.val,mMode_slider.val,nMode_slider.val))[0], E_mag[1], E_mag[2]] for z in z_range] for x in x_range])
        YZplaneField = np.array([[[(E_mag := E_ins(xCoord_slider.val,y,z,tCoord_slider.val,mMode_slider.val,nMode_slider.val))[0], E_mag[1], E_mag[2]] for z in z_range] for y in y_range])
        
    else:
        XYplaneField = np.array([[[(H_mag := H_ins(x,y,zCoord_slider.val,tCoord_slider.val,mMode_slider.val,nMode_slider.val))[0], H_mag[1], H_mag[2]] for x in x_range] for y in y_range])
        XZplaneField = np.array([[[(H_mag := H_ins(x,yCoord_slider.val,z,tCoord_slider.val,mMode_slider.val,nMode_slider.val))[0], H_mag[1], H_mag[2]] for z in z_range] for x in x_range])
        YZplaneField = np.array([[[(H_mag := H_ins(xCoord_slider.val,y,z,tCoord_slider.val,mMode_slider.val,nMode_slider.val))[0], H_mag[1], H_mag[2]] for z in z_range] for y in y_range])

    # Generate the field magnitude plot from the field vector plots 
    XYplaneMag = np.array([[np.sqrt(x**2 + y**2 + z**2) for x,y,z in row] for row in XYplaneField])
    XZplaneMag = np.array([[np.sqrt(x**2 + y**2 + z**2) for x,y,z in row] for row in XZplaneField])
    YZplaneMag = np.array([[np.sqrt(x**2 + y**2 + z**2) for x,y,z in row] for row in YZplaneField])
    
    # Transposing the array of field vectors to make it easier to split out the x, y and z components as needed
    XYplaneField = np.transpose(XYplaneField, (2,0,1))
    XZplaneField = np.transpose(XZplaneField, (2,0,1))
    YZplaneField = np.transpose(YZplaneField, (2,0,1))

    return  XYplaneField, XZplaneField, YZplaneField, XYplaneMag, XZplaneMag, YZplaneMag



# ********************************************
# ****** Drawing the windows *****************
# ********************************************
def drawContour(XYplaneMag, XZplaneMag, YZplaneMag, Norm):
    axXY.contourf(XYplane[0], XYplane[1], XYplaneMag, 100, norm=Norm)
    axXZ.contourf(XZplane[0], XZplane[1], XZplaneMag, 100, norm=Norm)
    axYZ.contourf(YZplane[0], YZplane[1], YZplaneMag, 100, norm=Norm)

def drawQuiver(XYplaneField, XZplaneField, YZplaneField, XYplaneMag, XZplaneMag, YZplaneMag, Norm):
    axXY.quiver(XYplane[0], XYplane[1], XYplaneField[0], XYplaneField[1], XYplaneMag, cmap='turbo', norm=Norm)
    axXZ.quiver(XZplane[0], XZplane[1], XZplaneField[2], XZplaneField[0], XZplaneMag, cmap='turbo', norm=Norm)
    axYZ.quiver(YZplane[0], YZplane[1], YZplaneField[2], YZplaneField[1], YZplaneMag, cmap='turbo', norm=Norm)

def drawStreamlines(XYplaneField, XZplaneField, YZplaneField, XYplaneMag, XZplaneMag, YZplaneMag, Norm):
    axXY.streamplot(XYplane[0], XYplane[1], XYplaneField[0], XYplaneField[1], color=XYplaneMag, cmap='turbo', norm=Norm)
    axXZ.streamplot(XZplane[0], XZplane[1], XZplaneField[2], XZplaneField[0], color=XZplaneMag, cmap='turbo', norm=Norm)
    axYZ.streamplot(YZplane[0], YZplane[1], YZplaneField[2], YZplaneField[1], color=YZplaneMag, cmap='turbo', norm=Norm)

def drawCrosshairs(x,y,z):
    axXY.axhline(y, color='red')
    axXY.axvline(x, color='red')
    axXZ.axhline(x, color='red')
    axXZ.axvline(z, color='red')
    axYZ.axhline(y, color='red')
    axYZ.axvline(z, color='red')


def update(val):

    axXY.clear()
    axXZ.clear()
    axYZ.clear()

    XYplaneField, XZplaneField, YZplaneField, XYplaneMag, XZplaneMag,YZplaneMag = calculateFields()

    # Color normalization
    vMax = max(nMode_slider.val,mMode_slider.val)
    Norm = colors.Normalize(vmin=0, vmax=vMax)

    # Contour plot for the background
    if backgButton.get_status()[0]:
        drawContour(XYplaneMag, XZplaneMag, YZplaneMag, Norm)

    # Quiver plot
    if plotButton.value_selected == 'Arrow plot':
        drawQuiver(XYplaneField, XZplaneField, YZplaneField, XYplaneMag, XZplaneMag,YZplaneMag, Norm)
    
    # Stream plot (This one takes a lot longer to draw. The plot becomes very unresponsive...)
    elif plotButton.value_selected == 'Streamplot':
        drawStreamlines(XYplaneField, XZplaneField, YZplaneField, XYplaneMag, XZplaneMag,YZplaneMag, Norm)

    if backgButton.get_status()[1]:
        drawCrosshairs(xCoord_slider.val, yCoord_slider.val, zCoord_slider.val)
    
    axXY.set_title('x-y plane')
    axYZ.set_title('y-z plane')
    axXZ.set_title('x-z plane')

    plt.draw()


update(0)


# ********************************************
# ****** Change handlers *********************
# ********************************************
xCoord_slider.on_changed(update)
yCoord_slider.on_changed(update)
zCoord_slider.on_changed(update)
tCoord_slider.on_changed(update)
mMode_slider.on_changed(update)
nMode_slider.on_changed(update)
fieldButton.on_clicked(update)
plotButton.on_clicked(update)
backgButton.on_clicked(update)

plt.show()
