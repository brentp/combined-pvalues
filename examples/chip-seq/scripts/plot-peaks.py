from matplotlib import pyplot
import sys

vals = [map(int, l.split()) for l in open(sys.argv[1]) if l[0] != 's']
xs, ys = zip(*vals)
pyplot.plot(xs, ys)
pyplot.show()
