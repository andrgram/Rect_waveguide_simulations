import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button


# The parametrised function to be plotted
def f(t, amplitude, frequency):
    return amplitude * np.sin(2 * np.pi * frequency * t)

t = np.linspace(0, 1, 1000)

# Define the initial parameters
initAmplitude = 5
initFrequency = 3

# Create the figure and the line that we will manipulate
fig, ax = plt.subplots()
line, = ax.plot(t, f(t, initAmplitude, initFrequency), lw=2)
ax.set_xlabel('Time [s]')

# Adjust the main plot to make room for the sliders
fig.subplots_adjust(left=0.25, bottom=0.25)

# Make a horizontal slider to control the frequency
axfreq = fig.add_axes([0.25, 0.1, 0.65, 0.03])
freqSlider = Slider(
    ax=axfreq,
    label='Frequency [Hz]',
    valmin=0.1,
    valmax=30,
    valinit=initFrequency,
)

# Make a vertical slider to control the amplitude
axamp = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
ampSlider = Slider(
    ax=axamp,
    label='Amplitude',
    valmin=0,
    valmax=10,
    valinit=initAmplitude,
    orientation='vertical'
)

# The function to be called anytime a slider's value changes
def update(val):
    line.set_ydata(f(t, ampSlider.val, freqSlider.val))
    fig.canvas.draw_idle()

# Register the update function with each slider
freqSlider.on_changed(update)
ampSlider.on_changed(update)

# Create a 'matplotlib.widgets.Button' to reset the sliders to initial values
resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')

def reset(event):
    freqSlider.reset()
    ampSlider.reset()

button.on_clicked(reset)

plt.show()