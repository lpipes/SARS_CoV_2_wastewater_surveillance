import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
from matplotlib.pyplot import MultipleLocator
import pylab as pl
import random

basepath = os.path.abspath(os.path.dirname(__file__))

f = open(os.path.join(basepath,'output.csv'),'r') #read the output file

rawloc_p = {}
p = 0
output_data = f.readlines()[1:]
for line in output_data:
	p = eval(line.strip().split(',')[1])*100
	info = line.strip().split(',')[2]
	rawlocation = info.strip('"').split(',')[0]
	# print(rawlocation)
	location = rawlocation.strip('/2020')
	# print(location)
	rawloc_p[p] = location

bar = 0.000726 * 100
loc_p = {}
for p in rawloc_p:
	if p > bar or p == bar:
		loc_p[rawloc_p[p]] = p
# print(loc_p)
d = (sorted(loc_p.items(),key=lambda item:item[1], reverse=True))

final_dic = {}
for i in d:
	final_dic[i[0]] = i[1]
p_list = []
loc_list = []

# print(final_dic)
for loc in final_dic:
	loc_list.append(loc)
	p_list.append(final_dic[loc])

colors = []


plt.barh(loc_list,p_list,color = 'cornflowerblue')

# plt.bar(loc_list,p_list,width = 0.5)
plt.title('Sequence Analysis')

x_major_locator = MultipleLocator(0.001)
ax=plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
plt.xticks(rotation = 90, fontsize = 7)
plt.yticks(fontsize = 7)

plt.ylabel('Sequence', fontsize=12)
plt.xlabel('Percentage(%)', fontsize=12)

for a, b in zip(p_list,loc_list):
	plt.text(a + bar/20, b , (b,'%.6f'%a), ha='center',fontsize=9)
plt.grid(ls='-.')
plt.xlim(min(p_list)*0.9, max(p_list)*1.1)
plt.tight_layout()

plt.savefig(os.path.join(basepath,'pic_plot'))

plt.show()
