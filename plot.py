import matplotlib.pyplot as plt

coordinates=[[113.69, 31.5261], [112.813, 31.0048], [113.632, 30.9871], [113.172, 31.0079]]
coordinates.append(coordinates[0])

x, y = zip(*coordinates)

plt.figure()
plt.plot(x, y, linestyle='dashed', marker='o', markerfacecolor='red', markersize=10)
plt.show()