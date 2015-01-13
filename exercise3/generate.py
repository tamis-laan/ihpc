#!/usr/bin/python

import subprocess
import json
import matplotlib.pyplot as plt

#Run a experiment
#	p  : number of processors
#	px : number of processors in x
#	py : number of processors in y
#	nx : grid width
#	ny : grid height

def run(p = 4, px = 2, py = 2, nx = 100, ny = 100, adapt = False):
	if(adapt):
		subprocess.call(["../build/GridDist", str(px), str(py), str(nx), str(ny),"adapt"])
	else:
		subprocess.call(["../build/GridDist", str(px), str(py), str(nx), str(ny)])

	subprocess.call(["mpirun", "-n", str(p), "../build/ihpc3"])
	data = []
	for i in range(0,p):
		data.append(json.load(open('meta'+str(i)+'.json')))
	return data

def experiment_3_2():
	grid = [100,200,400];
	fig = plt.figure()
	
	for p in [[1,4,211],[2,2,212]]:
		ax = fig.add_subplot(p[2])
		ax.set_title(str(p[0]) + "x" + str(p[1]))
		y = {"computation":[],"communication local":[],"communication global":[],"idle":[],"total":[]}
		for g in grid:
			data = run(4,p[0],p[1],g,g)
			y["computation"].append(data[0]["average"]["computation"])
			y["communication local"].append(data[0]["average"]["communication local"])
			y["communication global"].append(data[0]["average"]["communication global"])
			y["idle"].append(data[0]["average"]["idle"])
			y["total"].append(data[0]["average"]["total"])
		ax.set_xlabel('Grid')
		ax.set_ylabel('Time');
		ax.plot(grid,y["computation"],'r', label = "computation")
		ax.plot(grid,y["communication local"],'g', label= "comm. local")
		ax.plot(grid,y["communication global"],'b', label = "comm. global")
		ax.plot(grid,y["idle"],'y', label = "idle")
		ax.plot(grid,y["total"],'m', label = "total")
		ax.legend(loc='upper left')

	plt.show()

def experiment_3_5():
	x = [100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]
	y = []
	for i in x:
		data = run(4,2,2,i,i)
		y.append(data[0]["average"]["computation"]/(data[0]["average"]["communication global"]+data[0]["average"]["computation"]))
	fig = plt.figure()
	ax = fig.add_subplot(211)
	ax.set_xlabel('Grid')
	ax.set_ylabel("Comp. Comm. Ratio")
	ax.plot(x,y)

	p  = [1,2,3,4,5]
	y  = []
	for i in p:
		data = run(i*i,i,i,1000,1000)
		y.append(data[0]["average"]["computation"]/(data[0]["average"]["communication global"]+data[0]["average"]["computation"]))
	ax = fig.add_subplot(212)
	ax.set_xlabel('processors')
	ax.set_ylabel("Comp. Comm. Ratio")
	x = [x**2 for x in p]
	ax.plot(x,y)

	plt.show()

def experiment_3_6():
	
	grid = [100,200,400]
	
	itterations = {"True":[],"False":[]}
	total_time  = {"True":[],"False":[]}
	computation = {"True":[],"False":[]}

	for a in [True,False]:
		for g in grid:
			data = run(4,2,2,g,g,a)
			total_time[str(a)].append(data[0]["average"]["total"])
			computation[str(a)].append(data[0]["average"]["computation"])
			itterations[str(a)].append(data[0]["iterations"])
	
	fig = plt.figure()
	ax = fig.add_subplot(311)
	ax.set_title("Total Time")
	ax.set_xlabel('Grid')
	ax.plot(grid,total_time["True"],'b' , label = "Adaptive")
	ax.plot(grid,total_time["False"],'r' , label = "Not Adaptive")
	ax.legend(loc='upper left')
	ax = fig.add_subplot(312)
	ax.set_title("Computation Time")
	ax.set_xlabel('Grid')
	ax.plot(grid,computation["True"],'b' , label = "True")
	ax.plot(grid,computation["False"],'r' , label = "False")
	ax = fig.add_subplot(313)
	ax.set_title("itterations")
	ax.set_xlabel('Grid')
	ax.plot(grid,itterations["True"],'b' , label = "True")
	ax.plot(grid,itterations["False"],'r' , label = "False")
	plt.show()


experiment_3_6()


