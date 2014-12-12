#!/usr/bin/python

import subprocess
import json
import matplotlib.pyplot as plt

#Run a experiment
#	p  : number of processors
#	px : number of processors in x
#	py : number of processors in y
#	o  : bool to output file
#	w  : omega relaxation parameter 
#	nx : grid width
#	ny : grid height
#	pr : precision_goal
#	mi : max iterations
def run(p = 4, px = 2, py = 2, o = 0, w = 1.95, nx = 100, ny = 100, pr = 0.0001, mi = 5000):
	subprocess.call(["mpirun", "-n", str(p), "./build/ihpc2", str(px), str(py), str(o), str(w), str(nx), str(ny), str(pr), str(mi)])
	return json.load(open('meta.json'))

def experiment_22():
	exp = {"w":[],"it":[],"time":[]}
	for i in range(0,10):
		w = 1.9+i/100.0
		data = run(4,4,1,0,w)
		
		exp["w"].append(w)
		exp["it"].append(data["iterations"])
		exp["time"].append(data["time"])

	plt.plot(exp["w"],exp["it"])
	plt.show()
	plt.plot(exp["w"],exp["time"])
	plt.show()

def experiment_23():
	for c in [[4,1],[2,2]]:
		for g in [[200,200],[400,400],[800,800]]:
			data = run(4,c[0],c[1],0,1.95, g[0], g[1])

experiment_22()
