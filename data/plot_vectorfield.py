from matplotlib import pyplot as pt

def read_from_file (name):
	f = open(name, "r")
	data = f.read()
	f.close()
	return eval(data)

u = read_from_file("horizontal0.dat")
v = read_from_file("vertical0.dat")
p = read_from_file("pressure0.dat")

pt.quiver(v, u, p)
pt.show()
