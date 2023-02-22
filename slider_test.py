import numpy as np
import matplotlib.pyplot as plt
import matplotlib.text as txt
from matplotlib.widgets import Slider, Button
from matplotlib.patches import Rectangle



# Examining all modes up to the 2,2 mode
TE_modes_examined = [[1,0], [0,1], [1,1], [2,0], [2,1], [0,2], [1,2], [2,2]]
TM_modes_examined = [[1,0], [0,1], [1,1], [2,0], [2,1], [0,2], [1,2], [2,2]]

start_width = 10
start_height = 10

c = 3 * 10**8

# figure setup
fig, ax = plt.subplots()
ax.add_patch(Rectangle((0,0), start_width, start_width, fill=False))
# fig.subplots_adjust(left=0.25, bottom=0.25)
# t1 = txt.Text(
#     x=0,
#     y=0,
#     text='test1'
# )
# t2 = txt.Text(
#     x=0,
#     y=2,
#     text='test2'
# )


# Cut-off frequency for specific modes in a rectangular waveguide sized m*n
def fmn_c(m, n, a, b):
    return c / 2 * np.sqrt((m / a)**2 + (n / b)**2)

# Setting up the sliders
axWidth = fig.add_axes([0.15, 0, 0.7, 0.05])
width_slider = Slider(
    ax=axWidth,
    label='Waveguide width [mm]',
    valmin=1,
    valmax=30,
    valinit=start_width,
    valstep=1
)
axHeight = fig.add_axes([0, 0.15, 0.05, 0.7])
height_slider = Slider(
    ax=axHeight,
    label='Waveguide height [mm]',
    valmin=1,
    valmax=30,
    valinit=start_height,
    valstep=1,
    orientation='vertical'
)
axXModes = fig.add_axes([0.15, 0.05, 0.7, 0.05])
XModes_slider = Slider(axXModes, 'x-modes', 0, 4, 1, valstep=1)
axYModes = fig.add_axes([0.05, 0.15, 0.05, 0.7])
YModes_slider = Slider(axYModes, 'y-modes', 0, 4, 0, valstep=1, orientation='vertical')


texts = [ax.text(0, y/10, 'init') for y in range(len(TE_modes_examined))]

def update(val):
    frequencies = [fmn_c(mode[0], mode[1], width_slider.val * 10**-3, height_slider.val * 10**-3) for mode in TE_modes_examined]
    for i, f in enumerate(frequencies):
        texts[i].set_text(('TE' + str(TE_modes_examined[i][0]) + str(TE_modes_examined[i][1]) + ' ' + str(round(f*10**-9,1)) + 'GHz'))


width_slider.on_changed(update)
height_slider.on_changed(update)

ax.margins(0)
ax.set_axis_off()
plt.show()

testA = 1.1 * 10**-2
testB = 1.0 * 10**-2

TE10 = fmn_c(1, 0, testA, testB)
TE01 = fmn_c(0, 1, testA, testB)
TE11 = fmn_c(1, 1, testA, testB)
TE20 = fmn_c(2, 0, testA, testB)

print('Mode \t Cut-off')
print('TE10 \t', round(TE10 * 10**-9, 1), 'GHz')
print('TE01 \t', round(TE01 * 10**-9, 1), 'GHz')
print('TE11 \t', round(TE11 * 10**-9, 1), 'GHz')
print('TE20 \t', round(TE20 * 10**-9, 1), 'GHz')
