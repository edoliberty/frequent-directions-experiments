import csv
import ast
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import matplotlib
path = "output_time_vs_inputsize.csv"
x_ax = range(10000,110000,10000)

f = open(path, "rb")
csvF = csv.reader(f)
csvF.next()
algos = [[] for x in range(10)]
rowCnt = 0
font = {'family' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)
#plt.title("Alpha FD Random Data with Noise ")
for val in csvF:
    data = val[0:10]
    #print 'data=',data
    #break
    algoCnt = 0
    for algo_val in data:
        av = ast.literal_eval(algo_val) 
        algos[algoCnt].append(av)
        algoCnt += 1
    
    rowCnt += 1
    if rowCnt == 21:
        break
#print 'algos[9] = ',len(algos)
temp = 0
colors = ['r','g','b','k','c','y',(0.4,0,0.4),(1,0.6,0.3),(0.6,0.8,1),(0,0.8,0.4)]
for val in algos:
    #print 'val=',val
    #if temp == 0:
    #    temp += 1
    #    continue
    plt.plot(x_ax, val, linewidth = 3, color = colors[temp])
    temp += 1

plt.xlabel("Input Size", fontsize = 18)
plt.ylabel("Running Time",fontsize = 18)
#plt.ylim(ymin = 10000)
l=plt.legend(['1000', '2000', '3000', '4000', '5000', '6000' , '7000', '8000', '9000', '10000'], loc = 'upper left',prop = {'size' : 14})
l.draw_frame(False)
#plt.show()
ofn = "time_vs_inputsize.pdf"

plt.savefig(ofn)
os.system("pdfcrop %s %s" % (ofn, ofn))
