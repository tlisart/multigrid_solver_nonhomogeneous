from numpy import *
import matplotlib.pyplot as plt

data = loadtxt("data/plotRes.txt")

resNorm = data[:]
it = []

for i in range(len(data)) :
    it.append(i)

plt.plot(it, log(resNorm))

plt.xlabel('amount of iterations')
plt.ylabel('residual vector norm (log)')
plt.title('Graph of the residual norm in function of the amount of iterations (semilog)')

plt.grid()
plt.show()
