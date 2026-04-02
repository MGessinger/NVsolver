from matplotlib import pyplot as pt

def read_from_file (name):
	with open(name, "r") as f:
		text = f.read()
	return eval(text)

def plot_from_file (number):
	output_file = f'plot{number:08}.jpg'

	u = read_from_file("../data/horizontal{0}.dat".format(number))
	v = read_from_file("../data/vertical{0}.dat".format(number))
	p = read_from_file("../data/pressure{0}.dat".format(number))

	fig = pt.figure()
	ax = fig.add_subplot()

	ax.quiver(u, v, p)
	fig.savefig(output_file)

	pt.close(fig)

	return output_file

n = 0
while True:
	try:
		output_file = plot_from_file(n)
		print(f"Plotted files {output_file}")
		n = n + 1
	except:
		break
