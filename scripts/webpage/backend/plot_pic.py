import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from matplotlib.pyplot import MultipleLocator
import pylab as pl
import random


def picture(csv_path, out_path):
	"""
	ourput pictures after R processing
	:param csv_path: CSV path
	:param out_path: Picture path
	:return: None
	"""
	# basepath = os.path.abspath(os.path.dirname(__file__))

	# f = open(os.path.join(basepath,'em_output_test.txt.csv'),'r') #read the output file
	f = open(csv_path, 'r')

	rawloc_p = {}
	p = 0
	output_data = f.readlines()[1:]
	for line in output_data:
		p = eval(line.strip().split(',')[1].strip('"'))*100
		info = line.strip().split(',')[2]
		rawlocation = info.strip('"').split(',')[0]
		# print(rawlocation)
		location = rawlocation.strip('/2020')
		# print(location)
		rawloc_p[p] = location

	bar = 1
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
	plt.title('Analysis')

	# x_major_locator = MultipleLocator(0.001)
	# ax=plt.gca()
	# ax.xaxis.set_major_locator(x_major_locator)
	plt.xticks(rotation = 90, fontsize = 7)
	plt.yticks(fontsize = 7)

	plt.ylabel('Sequence', fontsize=12)
	plt.xlabel('Proportion(%)', fontsize=12)

	for a, b in zip(p_list,loc_list):
		plt.text(max(p_list), b , ('%.4f'%a), ha='center',fontsize=8)
	plt.grid(ls='-.')
	plt.xlim(min(p_list)*0.9, max(p_list)*1.1)
	plt.tight_layout()

	plt.savefig(os.path.join(out_path, 'pic_plot'))

	# plt.show()
