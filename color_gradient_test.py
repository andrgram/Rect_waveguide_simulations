import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
import matplotlib.streamplot as spl

import timeit
startTime = timeit.default_timer()
endTime = timeit.default_timer()

# Prints status every second
def printProgress(var, minVal, maxVal):
    global endTime
    currTime = timeit.default_timer()
    if currTime - endTime > 1:
        length = maxVal - minVal
        progress = int((var-minVal)/ length*100)
        timeLeft = int((currTime - startTime) * ((maxVal-(var+1))/((var+1)-minVal)))
        yearString = str(timeLeft // 31536000)+"y, " if timeLeft > 31536000 else ""
        dayString = str(timeLeft % 31536000 // 86400)+"d, " if timeLeft > 86400 else ""
        hourString = str(timeLeft % 86400 // 3600)+"h, " if timeLeft > 3600 else ""
        minuteString = str(timeLeft % 3600 // 60)+"m, " if timeLeft > 60 else ""
        secondString = str(timeLeft % 60)+"s"
        timeLeftString = yearString + dayString + hourString + minuteString + secondString
        print("runTime: "+str(int((currTime - startTime)))+"s. Progress: "+str(progress)+"%, "+str((var)-minVal)+"/"+str(length)+" calculations complete. Estimated time remaining: "+timeLeftString+"                          ", end="\r")
        endTime = currTime


# TODO Find out about the slider color gradient...
# TODO Runtime optimalizations....
# TODO Boundary cases.
# TODO Axis text
# TODO Examine quiver plots and streamplots
# TODO Button for crosshairs and for the background contour
# TODO clim for the streamplot




quiverPoints = 10
points = 10
x_points = points
y_points = points
z_points = int(points * 2 * np.pi)
t_points = int(points * 2 * np.pi)
# x_points = 10
# y_points = 10
# z_points = 50
# t_points = 50

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
# ****** Setup for the sliders ***************
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
yCoord_slider = Slider(axYcoord, 'y-coordinate', 0, 1, 0.5, valstep=1/y_points)
axZcoord = fig.add_axes([0.15, 0.05, 0.7, 0.05])
zCoord_slider = Slider(axZcoord, 'z-coordinate', 0, 2*np.pi, 0, valstep=2*np.pi/z_points)
axTcoord = fig.add_axes([0.15, 0, 0.7, 0.05])
tCoord_slider = Slider(axTcoord, 'time', 0, 2*np.pi, 0, valstep=2*np.pi/t_points)
axMmode = fig.add_axes([axXY.get_position().x0, axXY.get_position().y0 - 0.1, axXY.get_position().width, 0.05])
mMode_slider = Slider(axMmode, 'modes in x-direction', 0, 4, 1, valstep=1)
axNmode = fig.add_axes([axXY.get_position().x0, axXY.get_position().y0 - 0.15, axXY.get_position().width, 0.05])
nMode_slider = Slider(axNmode, 'modes in y-direction', 0, 4, 0, valstep=1)
axFieldButton = fig.add_axes([axXY.get_position().x0, axXY.get_position().y0 - 0.3, 0.05, 0.1])
fieldButton = RadioButtons(axFieldButton, ('E-field', 'H-field'))
axPlotButton = fig.add_axes([axXY.get_position().x0 + 0.1, axXY.get_position().y0 - 0.3, 0.05, 0.1])
plotButton = RadioButtons(axPlotButton, ('Arrow plot', 'Streamplot'))



# ********************************************
# ****** Functions for the fields ************
# ********************************************
# Field and wave Electromagnetics p.553
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




# There should be a button for activating or deactivating this!
def drawCrosshairs(x,y,z):
    axXY.axhline(y, color='red')
    axXY.axvline(x, color='red')
    axXZ.axhline(x, color='red')
    axXZ.axvline(z, color='red')
    axYZ.axhline(y, color='red')
    axYZ.axvline(z, color='red')



# Are the clear operations needed?
# They are needed when using the crosshairs!
def update(val):

    startTime = timeit.default_timer()

    axXY.clear()
    axXZ.clear()
    axYZ.clear()

    clearTime = timeit.default_timer()

    # Deciding whether to plot the E-field or the H-field
    if fieldButton.value_selected == 'E-field':
        field_XY = np.array([[[(E_mag := E_ins(x,y,zCoord_slider.val,tCoord_slider.val,mMode_slider.val,nMode_slider.val))[0], E_mag[1], E_mag[2]] for x in x_range] for y in y_range])
        field_XZ = np.array([[[(E_mag := E_ins(x,yCoord_slider.val,z,tCoord_slider.val,mMode_slider.val,nMode_slider.val))[0], E_mag[1], E_mag[2]] for z in z_range] for x in x_range])
        field_YZ = np.array([[[(E_mag := E_ins(xCoord_slider.val,y,z,tCoord_slider.val,mMode_slider.val,nMode_slider.val))[0], E_mag[1], E_mag[2]] for z in z_range] for y in y_range])
        
        XYmagPlane = np.array([[np.sqrt(coord[0]**2 + coord[1]**2 + coord[2]**2) for coord in row] for row in field_XY])
        XZmagPlane = np.array([[np.sqrt(coord[0]**2 + coord[1]**2 + coord[2]**2) for coord in row] for row in field_XZ])
        YZmagPlane = np.array([[np.sqrt(coord[0]**2 + coord[1]**2 + coord[2]**2) for coord in row] for row in field_YZ])

        field_XY = np.transpose(field_XY, (2,0,1))
        field_XZ = np.transpose(field_XZ, (2,0,1))
        field_YZ = np.transpose(field_YZ, (2,0,1))

    else:
        field_XY = np.array([[[(H_mag := H_ins(x,y,zCoord_slider.val,tCoord_slider.val,mMode_slider.val,nMode_slider.val))[0], H_mag[1], H_mag[2]] for x in x_range] for y in y_range])
        field_XZ = np.array([[[(H_mag := H_ins(x,yCoord_slider.val,z,tCoord_slider.val,mMode_slider.val,nMode_slider.val))[0], H_mag[1], H_mag[2]] for z in z_range] for x in x_range])
        field_YZ = np.array([[[(H_mag := H_ins(xCoord_slider.val,y,z,tCoord_slider.val,mMode_slider.val,nMode_slider.val))[0], H_mag[1], H_mag[2]] for z in z_range] for y in y_range])
        
        XYmagPlane = np.array([[np.sqrt(coord[0]**2 + coord[1]**2 + coord[2]**2) for coord in row] for row in field_XY])
        XZmagPlane = np.array([[np.sqrt(coord[0]**2 + coord[1]**2 + coord[2]**2) for coord in row] for row in field_XZ])
        YZmagPlane = np.array([[np.sqrt(coord[0]**2 + coord[1]**2 + coord[2]**2) for coord in row] for row in field_YZ])

        field_XY = np.transpose(field_XY, (2,0,1))
        field_XZ = np.transpose(field_XZ, (2,0,1))
        field_YZ = np.transpose(field_YZ, (2,0,1))


    calculateTime = timeit.default_timer()


    # Contour plot for the background
    axXY.contourf(XYplane[0], XYplane[1], XYmagPlane, 100, vmin=0, vmax=1)
    axXZ.contourf(XZplane[0], XZplane[1], XZmagPlane, 100, vmin=0, vmax=1)
    axYZ.contourf(YZplane[0], YZplane[1], YZmagPlane, 100, vmin=0, vmax=1)

    if plotButton.value_selected == 'Arrow plot':
        # Quiver plot
        axXY.quiver(XYplane[0], XYplane[1], field_XY[0], field_XY[1], XYmagPlane, cmap='turbo', clim=[0,1])
        axYZ.quiver(YZplane[0], YZplane[1], field_YZ[2], field_YZ[1], YZmagPlane, cmap='turbo', clim=[0,1])
        axXZ.quiver(XZplane[0], XZplane[1], field_XZ[2], field_XZ[0], XZmagPlane, cmap='turbo', clim=[0,1])

    else:
        # Stream plot
        axXY.streamplot(XYplane[0], XYplane[1], field_XY[0], field_XY[1], color=XYmagPlane, cmap='turbo')#, density=0.5)
        axYZ.streamplot(YZplane[0], YZplane[1], field_YZ[2], field_YZ[1], color=YZmagPlane, cmap='turbo')#, density=0.5)
        axXZ.streamplot(XZplane[0], XZplane[1], field_XZ[2], field_XZ[0], color=XZmagPlane, cmap='turbo')#, density=0.5)

    # TODO Check if there is an error on the x/y density in streamplot


    drawTime = timeit.default_timer()

    drawCrosshairs(xCoord_slider.val, yCoord_slider.val, zCoord_slider.val)

    crosshairTime = timeit.default_timer()
    
    axXY.set_title('x-y plane')
    axYZ.set_title('y-z plane')
    axXZ.set_title('x-z plane')

    titleTime = timeit.default_timer()

    plt.draw()

    totalTime = timeit.default_timer()
    print('Clearing took', round((clearTime - startTime)*1000), 'ms')
    print('Calculations took', round((calculateTime - clearTime)*1000), 'ms')
    print('Drawing took', round((drawTime - calculateTime)*1000), 'ms')
    print('Drawing crosshairs took', round((crosshairTime - drawTime)*1000), 'ms')
    print('Drawing titles took', round((titleTime - crosshairTime)*1000), 'ms')
    print('The draw command took', round((totalTime - titleTime)*1000), 'ms')
    print('Total function time:', round((totalTime - startTime)*1000), 'ms\n')

E_field = 1

update(0)


# ********************************************
# ****** Change handler **********************
# ********************************************
xCoord_slider.on_changed(update)
yCoord_slider.on_changed(update)
zCoord_slider.on_changed(update)
tCoord_slider.on_changed(update)
mMode_slider.on_changed(update)
nMode_slider.on_changed(update)
fieldButton.on_clicked(update)
plotButton.on_clicked(update)

plt.show()
