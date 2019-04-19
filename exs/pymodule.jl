using PyCall # This enables call many python module
pygui(:qt5) # Choose backend
using PyPlot # This is alternative to get matplotlib by pyplot

# This is a way to call python module
# Example to call matplotlib

plt = pyimport("matplotlib.pyplot")
# Then we can use same syntax as matplotlib
x = range(0;stop=2*pi,length=1000); y = sin.(3*x + 4*cos.(2*x));
plt.plot(x, y, color="red", linewidth=2.0, linestyle="--")
plt.show()

# Alternatively, we can use pyplot since we call it:
# x = range(0; stop=2*pi, length=1000); y = sin.(3 * x + 4 * cos.(2 * x));
#plot(x, y, color="red", linewidth=2.0, linestyle="--")

# Example call numpy
np = pyimport("numpy")

# Some parameters
timesteps = 100
N = 100

phi = np.zeros((timesteps,N)) # This makes 100 x 100 array in Float64
