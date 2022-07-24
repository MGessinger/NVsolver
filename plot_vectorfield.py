from matplotlib import pyplot as pt

def read_from_file (name):
	with open(name, "r") as f:
		text = f.read()
	return eval(text)

u = read_from_file("../data/horizontal0.dat")
v = read_from_file("../data/vertical0.dat")
p = read_from_file("../data/pressure0.dat")

pt.quiver(u, v, p)
pt.show()
