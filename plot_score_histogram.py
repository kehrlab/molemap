import matplotlib.pyplot as plt
import sys

if len(sys.argv)!=2:
    print("usage: python plot_score_histogram.py bcmap_output.hist")
    quit()
    
hist_file=open(sys.argv[1],"r")
freq=[int(line) for line in hist_file]
scores=[i for i in range(200)]

plt.bar(scores,freq,color='#7070c0')
plt.xlabel("score")
plt.ylabel("frequency")
plt.title("Score Histogram")# for "+sys.argv[1][:-5])
plt.show()
