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
def run(p = 4, px = 2, py = 2, o = 0, w = 1.95, nx = 100, ny = 100, pr = 0.0001, mi = 5000, sweeps = 1, optimized_step = 0):
	subprocess.call(["mpirun", "-n", str(p), "./build/ihpc2", str(px), str(py), str(o), str(w), str(nx), str(ny), str(pr), str(mi), str(sweeps), str(optimized_step)])
	return json.load(open('meta.json'))

def experiment_22():
	exp = {"w":[],"it":[],"time":[]}
	for i in range(0,10):
		w = 1.9+i/100.0
		data = run(4,4,1,0,w)
		
		exp["w"].append(w)
		exp["it"].append(data["iterations"])
		exp["time"].append(data["time"])

	fig = plt.figure()
	fig.suptitle('Exersize 2.2', fontsize=11, fontweight='bold')
	ax1 = fig.add_subplot(211)
	ax2 = fig.add_subplot(212)
	ax1.set_xlabel('Omega')
	ax2.set_xlabel('Omega')
	ax1.set_ylabel('Iterations')
	ax2.set_ylabel('Time')
	ax1.plot(exp["w"],exp["it"],'r')
	ax2.plot(exp["w"],exp["time"],'r')
	plt.show()

def experiment_2_3():
	for c in [[4,1],[2,2]]:
		for g in [[200,200],[400,400],[800,800]]:
			data = run(4,c[0],c[1],0,1.95, g[0], g[1])
			print str(c) + " " + str(g) + " : " + str(data["time"])

def experiment_2_5():
	fig = plt.figure()
	fig.suptitle('Exersize 2.5', fontsize=11, fontweight='bold')
	ax = fig.add_subplot(111)
	ax.set_xlabel('Grid size')
	ax.set_ylabel('Iterations')

	exp = {"grid":[],"it":[]}
	for g in [200,300,400,500,600]:
		data = run(4,2,2,0,1.95, g, g, 0.0001, 999999999)
		exp["grid"].append(g)
		exp["it"].append(data["iterations"])

	ax.plot(exp["grid"],exp["it"],'r')
	plt.show()

def experiment_2_6():
	data = run(4,2,2,0,1.95, 500, 500, 0.0001, 999999999, 2)
	fig = plt.figure()
	fig.suptitle('Exersize 2.6', fontsize=11, fontweight='bold')
	ax = fig.add_subplot(111)
	ax.set_yscale('log')
	ax.set_xlabel('Itterations')
	ax.set_ylabel('Error')
	ax.plot(range(0,data["iterations"]),data["errors"],'r')
	plt.show()

def experiment_2_8():
	fig = plt.figure()
	fig.suptitle('Exersize 2.8', fontsize=11, fontweight='bold')
	ax1 = fig.add_subplot(211)
	ax2 = fig.add_subplot(212)
	exp = {"sweeps":[],"iterations":[],"time":[]}
	for s in range(1,20):	
		data = run(4,2,2,0,1.0, 400, 400, 0.0001, 999999999, s)
		exp["sweeps"].append(s)
		exp["iterations"].append(data["iterations"])
		exp["time"].append(data["time"])

	ax1.set_xlabel('Sweeps')
	ax2.set_xlabel('Sweeps')
	ax1.set_ylabel('Iterations')
	ax2.set_ylabel('Time')

	ax1.plot(exp["sweeps"],exp["iterations"],'r')
	ax2.plot(exp["sweeps"],exp["time"],'r')
	plt.show()

def experiment_2_9():
	fig = plt.figure()
	fig.suptitle('Exersize 2.9', fontsize=11, fontweight='bold')
	ax = fig.add_subplot(111)
	ax.set_xlabel("Grid size")
	ax.set_ylabel("Time")

	exp = {"grid":[],"time1":[],"time2":[]}
	for g in [50,100,150,200,250,300,350,400,450,500,550,600]:
		print g
		data = run(4,2,2,0,1.0, g, g, 0.0001, 999999999, 1, 0)
		exp["grid"].append(g)
		exp["time1"].append(data["time"])
		data = run(4,2,2,0,1.0, g, g, 0.0001, 999999999, 1, 1)
		exp["time2"].append(data["time"])
	ax.plot(exp["grid"],exp["time1"],'b',label = "Not optimized")
	ax.plot(exp["grid"],exp["time2"],'r',label = "optimized")
	ax.legend()
	plt.show()

def experiment_2_11():
	# Data being communicated across borders in bytes excluding overheid
	borderData = lambda nx,ny,px,py: 8*((px-1)*nx+(py-1)*ny)*2

experiment_2_10()