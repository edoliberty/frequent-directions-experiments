import csv
import ast
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import matplotlib
#path = "output_d_10.csv"
path = "projected_sketch_d_10.csv"
x_ax = range(10,210,10)


f = open(path, "rb")
csvF = csv.reader(f)
csvF.next()
algos = [[] for x in range(6)]
rowCnt = 0
font = {'family' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)
#plt.title("Alpha FD Random Data with Noise ")
for val in csvF:
    data = val[0:6]
    #print 'data=',data
    algoCnt = 0
    for algo_val in data:
        av = ast.literal_eval(algo_val) 
        algos[algoCnt].append(av)
        algoCnt += 1
    
    rowCnt += 1
    if rowCnt == 21:
        break
#print 'algo = ',algos
for val in algos:
    #print len(x_ax)
    #print val
    plt.plot(x_ax, val, linewidth = 3)

plt.xlabel("Sketch Size",fontsize = 18)
plt.ylabel("Projection Error", fontsize = 18)
#plt.ylim(ymin = 10000)
l = plt.legend(['Brute Force', 'Frequent Directions', 'Sampling', 'Hashing', 'Random Projections', 'Naive'], loc = 'upper right',prop = {'size' : 14})
l.draw_frame(False)
#l.fontsize = 'small'
#plt.show()
ofn = "prj_err_sketch_d_10.pdf"

plt.savefig(ofn)
os.system("pdfcrop %s %s" % (ofn, ofn))
