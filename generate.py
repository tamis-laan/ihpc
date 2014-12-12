import subprocess
import json

#Run a experiment
def run(p,px,py,o,w):
	subprocess.call(["mpirun", "-n", str(p), "./build/ihpc2", str(px), str(py), str(o), str(w)])
	return json.load(open('meta.json'))

data = run(4,2,2,0,1.9)
