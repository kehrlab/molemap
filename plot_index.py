import matplotlib.pyplot as plt
import numpy as np



file=open("abundances.txt","r")
Y=[int(x) for x in file.read().split(" ")[:-1]]
X=list(range(0,len(Y)))
plt.plot(X,Y,linewidth=2)
plt.title("Abundance of k-mers")
plt.xlabel('k-mers')
plt.ylabel('abundance of k-mers')
plt.xlim(-round(len(Y)/100),round(len(Y)*1.01))
plt.ylim(-Y[0]/100,Y[0]*1.01)
# plt.yscale("log")
plt.show()

file.close()
