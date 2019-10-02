import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as patches
import matplotlib.lines as lines
from matplotlib.backends.backend_pdf import PdfPages
plt.style.use('ggplot')
plt.rcParams['axes.facecolor'] = 'white'


def plotCNV(matrix_path, output_path, project, percentage=False, aggregate=False):
	N=48
	ticks = np.arange(N)
	width = 0.27
	xticks = []
	fig, ax = plt.subplots(figsize=(20,5))

	df = pd.read_csv(matrix_path, sep='\t')
	#check that matrix is in proper format

	pp = PdfPages(output_path + 'CNV_plot' + project + '.pdf')

	super_class = ['Het', 'LOH', "Hom del"]
	het_sub_class = ['amp+', 'amp', 'gain', 'neut']
	loh_subclass = ['amp+', 'amp', 'gain', 'neut', "del"]
	hom_del_class = ['0-100kb', '100kb-1Mb', '>1Mb']
	x_labels = ['0-100kb', '100kb-1Mb', '1Mb-10Mb', '10Mb-40Mb','>40Mb']
	color_mapping = {'amp+':{'>40Mb':"deeppink", '10Mb-40Mb':"hotpink", '1Mb-10Mb':"palevioletred", '100kb-1Mb':"lightpink", '0-100kb':"lavenderblush"}, 
				 'amp':{'>40Mb':"saddlebrown", '10Mb-40Mb':"sienna", '1Mb-10Mb': "peru", '100kb-1Mb':"sandybrown", '0-100kb':"linen"}, 
				 'gain':{'>40Mb': "rebeccapurple", '10Mb-40Mb':"blueviolet", '1Mb-10Mb':"mediumorchid", '100kb-1Mb':"plum", '0-100kb':"thistle"}, 
				 'neut':{'>40Mb':"olive", '10Mb-40Mb':"olivedrab", '1Mb-10Mb':"yellowgreen", '100kb-1Mb':"lawngreen", '0-100kb':"greenyellow"}, 
				 'del':{'>40Mb':"dimgray", '10Mb-40Mb':"darkgrey", '1Mb-10Mb':"silver", '100kb-1Mb':"lightgray", '0-100kb':"whitesmoke"}}

	hom_del_color_mapping = {'0-100kb':"darkblue" , '100kb-1Mb':"mediumblue", '>1Mb':"cornflowerblue"}


	#ADD BARS
	#each column vector in dataframe contains counts for a specific sample
	i = -1 #used to distinguish first bar from the rest 
	if aggregate:
		df['total_count'] = df.sum(axis=1)
	else: #one plot for each sample all combined in a single pdf		
		for col in df.columns[1:]:
			counts = list(df[col])
			for count, label in zip(counts, df['classification']):
				print (label)
				c = label.split(':')
				c2 = c[0]
				x = c[2]
				print (c)
				print (c2)
				print (x)
				if i == 0: #very first bar
					if count > 0:  				
						ax.bar(ticks[i], count, color=color_mapping[c2][x], edgecolor='black')
					else: #count of 0--don't plot anything but put a tick
						xticks.append(ticks[i])
				else:
					if count > 0: 
						if percentage:
							pass
						else:
							ax.bar(ticks[i] + width, count, color=color_mapping[c2][x], edgecolor='black')
					else: #count of 0--don't plot anything but put a tick
						xticks.append(ticks[i])

		
			ax.set_xticks(xticks);
			ax.set_xticklabels(x_labels * 9 + hom_del_class, rotation=90);


			#ADD PATCHES AND TEXT
			patch_height = 0.09
			patch_width = 0.075
			left_edge = 0.149 #placement of left edge of patch
			y_pos = 0.9 #placement of patch on y-axis
			text_height = 0.91		 
			patch_colors = ['maroon', 'darkorange', 'slateblue', 'green', 'maroon', 'darkorange', 'slateblue', 'green', 'slategrey', 'blue']
			categories = het_sub_class + loh_subclass + ['Hom Del']
			ax.add_patch(plt.Rectangle((left_edge, y_pos), patch_width, patch_height, clip_on=False, facecolor=patch_colors[0], transform=plt.gcf().transFigure))
			plt.text(left_edge, text_height, categories[0], fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
			for i in range(1, 10): #add remaining 9 patches
				left_edge = left_edge + patch_width
				ax.add_patch(plt.Rectangle((left_edge, y_pos), patch_width, patch_height, clip_on=False, facecolor=patch_colors[i], transform=plt.gcf().transFigure))
				plt.text(left_edge, text_height, categories[i], fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)

			#add final patch for hom del
			ax.add_patch(plt.Rectangle((left_edge, y_pos), patch_width, patch_height*2, clip_on=False, facecolor=patch_colors[0], transform=plt.gcf().transFigure))

			#manually add top level patches(het and LOH)
			ax.add_patch(plt.Rectangle((0.149, 0.9+y), 0.3, y, clip_on=False, facecolor='gray', transform=plt.gcf().transFigure))
			ax.add_patch(plt.Rectangle((0.149 + (width * 4), 0.9+y), width*5, y, clip_on=False, facecolor='black', transform=plt.gcf().transFigure))

			#hide y-axis ticks and labels 
			plt.gca().get_yaxis().set_ticks([])
			plt.gca().get_yaxis().set_ticklabels([])



			#ADD SEPARATION LINES

			if percentage:
				plt.ylabel("Percentage", fontsize=35, fontname="Times New Roman", weight = 'bold')
			else:
				plt.ylabel("Frequency", fontsize=35, fontname="Times New Roman", weight = 'bold')

			#SAVE TO PDF
		pp.savefig(plt.gcf())
		plt.close()
	pp.close()


	def plot():
		pass




if __name__ == "__main__":
	matrix_path = '/Users/azhark/Documents/Alexandrov_Lab/CNV/SigProfilerMatrixGenerator/SigProfilerMatrixGenerator/scripts/CNV.matrix.tsv'
	project = 'ASCAT_test'
	output_path = ''
	plotCNV(matrix_path, output_path, project, percentage=False)




	