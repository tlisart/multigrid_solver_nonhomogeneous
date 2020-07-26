from numpy import *
import matplotlib.pyplot as plt

data = loadtxt("data/timeData.txt")
time = data[:]
it = []

for i in range(len(data)) :
    it.append(((i*8 + 1) - 1)*2 + 1)

plt.plot(time, it)
plt.xlabel('time (s)')
plt.ylabel('fine matrix discretization m')
plt.title('Graph of the time per resolution in function of the amount of discretization')

plt.grid();
plt.show()
