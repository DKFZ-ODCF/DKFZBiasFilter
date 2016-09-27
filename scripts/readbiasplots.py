"""
        This module provides access to plot functions used for visualizing
        read biases in DNA mutations detected from NGS data.
"""
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#################
# Help Routines #
#################

def calculateRootedSize(size, max_val):
	"""
		Return square root of normalized size.

		Args
		----
		  size: float
		  max_val: float

		Value
		-----
		  float
		    np.sqrt(size/max_val)
	"""
	if(float(size) != 0.0):
		return np.sqrt(size/max_val)
	else:
		return 0.0

#####################
# Plotting Routines #
#####################

def plotSquareSizeLegend(ax, colormap, min_val, max_val, max_mutation_count):
	"""
		Plot legend assigning numbers to the size of squares

		Args
		----
		  ax: matplotlib.aces.Axes instance
		  colormap: string
		    Identifier of colormap being used for plotting (, e.g. "PRGn")
		  min_val: float
		  max_val: float
		  max_mutation_count: int
		    Maximal number of mutations being displayed
	"""
	stepsize = (max_mutation_count-1)/8.

	# Plot square size legend
	freq_list_legend = [[0],[0],[0],[0]]
	error_list_legend = [[0],[0],[0],[0]]
	text_list = [[""]*4, [""]*4, [""]*4, [""]*4]
	for i in range(0, 8):
		if(i==0):
			freq_list_legend[0] += [1.0]
			error_list_legend[0] += [0.0]
			text_list[0][0] = "1"
		elif(i==7):
			freq_list_legend[int(float(i)/2.)] += [float(max_mutation_count)]
			error_list_legend[int(float(i)/2.)] += [0.0]
			text_list[3][3] = str(int(max_mutation_count))
		else:
			if(i<=3):
				freq_list_legend[i] += [int(1.+i*(stepsize))]
				error_list_legend[i] += [0.0]
				text_list[i][0] = str(int(1.+i*(stepsize)))
			else:
				freq_list_legend[i-4] += [int(1.+i*(stepsize))]
				error_list_legend[i-4] += [0.0]
				text_list[i-4][3] = str(int(1.+i*(stepsize)))

	hintonLegend(np.array(freq_list_legend), np.array(error_list_legend), np.array(text_list), colormap, min_val, max_val, max_weight=max_mutation_count, ax=ax)

def hinton(weight_matrix, intensity_matrix, cmap, vmin, vmax, max_weight=None, ax=None):
	"""
		Draw Hinton diagram for visualizing a weight matrix.

		Args
		----
		  weight_matrix: np.array
		  intensity_matrix: np.array
		  cmap: string
		    Identifier of colormap being used for plotting (e.g. "PRGn")
		  vmin: float
		    Minimal value to be displayed in intensity matrix
		  vmax: float
		    Maximal value to be displayed in intensity matrix
		  max_weight: int
		    Force maximal weight of weight matrix
		  ax: matplotlib.axes.Axes instance
	"""
	ax = ax if ax is not None else plt.gca()
	
	# Set colors for intensity matrix
	cm = plt.get_cmap(cmap)
	cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
	scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
	intensity_colors = scalarMap.to_rgba(intensity_matrix)
	
	ax.patch.set_facecolor('gray')
	ax.set_aspect('equal', 'box')
	ax.xaxis.set_major_locator(plt.NullLocator())
	ax.yaxis.set_major_locator(plt.NullLocator())
	
	for (x,y),w in np.ndenumerate(weight_matrix):
		color = intensity_colors[x][y]
		size = 0.
		if(not(w==0)):
			size = calculateRootedSize(float(w), float(weight_matrix.max()))

		if(not(max_weight == None)):
			size = 0.
			if(not(w==0)):
				size = calculateRootedSize(float(w), float(max_weight))
		rect = plt.Rectangle([(3-y) - size / 2, x - size / 2], size, size,
		                     facecolor=color, edgecolor=color)
		ax.add_patch(rect)
	
	plt.ylim([-1,4])
	plt.xlim(-1,4)
	ax.invert_xaxis()
	ax.invert_yaxis()

def hintonLegend(weight_matrix, intensity_matrix, text_matrix, cmap, vmin, vmax, max_weight=None, ax=None):
	"""
		Draw Hinton diagram for visualizing a legend describing the number of mutations corresponding to sizes of squares.

		Args
		----
		  weight_matrix: np.array
		  intensity_matrix: np.array
		  text_matrix: np.array
		  cmap: string
		    Identifier of colormap being used for plotting (e.g. "PRGn")
		  vmin: float
		    Minimal value to be displayed in intensity matrix
		  vmax: float
		    Maximal value to be displayed in intensity matrix
		  max_weight: int
		    Force maximal weight of weight matrix
		  ax: matplotlib.axes.Axes instance
	"""
	ax = ax if ax is not None else plt.gca()
	
	# Set colors for intensity matrix
	cm = plt.get_cmap(cmap)
	cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
	scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
	intensity_colors = scalarMap.to_rgba(intensity_matrix)
	
	ax.patch.set_facecolor('gray')
	ax.set_aspect('equal', 'box')
	ax.xaxis.set_major_locator(plt.NullLocator())
	ax.yaxis.set_major_locator(plt.NullLocator())
	
	for (x,y),w in np.ndenumerate(weight_matrix):
		color = intensity_colors[x][y]
		size = 0.
		if(not(w==0)):
			size = calculateRootedSize(float(w), float(weight_matrix.max()))

		if(not(max_weight == None)):
			size = 0.
			if(not(w==0)):
				size = calculateRootedSize(float(w), float(max_weight))
		rect = plt.Rectangle([(3-y) - size / 2, x - size / 2], size, size,
		                     facecolor=color, edgecolor=color)
		ax.add_patch(rect)

	for (x,y),w in np.ndenumerate(text_matrix):
		ax.add_patch(rect)
		plt.text(3-y, x, w)

	
	plt.ylim([-1,4])
	plt.xlim(-1,4)
	ax.invert_xaxis()
	ax.invert_yaxis()


def plotErrorMatrix(error_matrix, mutation_count_matrix, error_type, pdf, is_bias=False):
	"""
		Draw read bias matrix as hinton plots.

		Args
		----
		  error_matrix: dictionary
		    Contains for each possible mutation and each possible nucleotide triplet read counts
		  mutation_count_matrix: dictionary
		    Contains for each possible mutation and each possible nucleotide triplet the number of mutations
		  error_type: string
		    Title string of plot
		  pdf: matplotlib.backend_pdf.PdfPages instance
		  is_bias: bool
		    Must be set to true, if error_matrix[mut][base_before][base_after] is int, and false if error_matrix[mut][base_before][base_after] is list of int
	"""
	possible_mutations = ["CA", "CG", "CT", "TA", "TC", "TG"]

	# Set figure properties
	figure = plt.figure(figsize=(28,4), dpi=80)
	figure.subplots_adjust(wspace=0.1, bottom=.2, top=.88)
	gs = gridspec.GridSpec(1, 8, height_ratios=[1], width_ratios=[1,1,1,1,1,1,1,0.2])
	colormap="PRGn"
	min_val = -3
	max_val = 3
	number_ticks_colormap = 5.

	# Calculate Maximal mutation count (necessary for normalization)
	max_mutation_count = 0
	for mutation in possible_mutations:
		for base_before in ["A", "C", "G", "T"]:
			for base_after in ["A", "C", "G", "T"]:
				if(mutation_count_matrix[mutation][base_before][base_after] > max_mutation_count):
					max_mutation_count = mutation_count_matrix[mutation][base_before][base_after]

	# Plot square size legend
	ax = plt.subplot(gs[6])
	plotSquareSizeLegend(ax, colormap, min_val, max_val, max_mutation_count)
	plt.title("Number of Mutations")

	counter = 0
	for mutation in possible_mutations:
		counter += 1
		current_error_list = []
		current_frequency_list = []

		for base_before in ["T", "G", "C", "A"]:
			base_before_specific_error_list = []
			base_before_specific_frequency_list = []

			for base_after in ["A", "C", "G", "T"]:
				# Convert counts to log_2 ratios
				if(not(is_bias)):
					base_before_specific_error_list += [math.log(float(error_matrix[mutation][base_before][base_after][0]+1)/float(error_matrix[mutation][base_before][base_after][1]+1),2)]
					base_before_specific_frequency_list += [mutation_count_matrix[mutation][base_before][base_after]]
				else:
					base_before_specific_error_list += [float(error_matrix[mutation][base_before][base_after])]
					base_before_specific_frequency_list += [mutation_count_matrix[mutation][base_before][base_after]]
			current_error_list += [base_before_specific_error_list]
			current_frequency_list += [base_before_specific_frequency_list]

		# Plotting
		ax = plt.subplot(gs[counter-1])
		hinton(np.array(current_frequency_list), np.array(current_error_list), colormap, min_val, max_val, max_weight=max_mutation_count, ax=ax)

#		ax.pcolor(np.array(current_error_list), cmap=plt.get_cmap(colormap), vmin=min_val, vmax=max_val)
		if(mutation[0] == "C" and mutation[1] == "A"):
			plt.title(error_type+"\n"+mutation[0]+"->"+mutation[1])
		else:
			plt.title(mutation[0]+"->"+mutation[1])
		if(counter == 1):
			plt.yticks([0, 1, 2, 3], ("T", "G", "C", "A"))
			plt.ylabel("preceeding")
		else:
			plt.yticks([0, 1, 2, 3], ("", "", "", ""))
		plt.xticks([0, 1, 2, 3], ["T", "G", "C", "A"])
		plt.xlabel("following")

	ax = plt.subplot(gs[7:])
	ax.yaxis.tick_right()
	ax.imshow(np.outer(np.arange(0,1,0.01), np.ones(10)), aspect='auto', cmap=colormap, origin="lower")
	plt.xlim([0,1])
	plt.xticks([],[])
	plt.yticks(np.arange(0.*100,1.*100+((100./number_ticks_colormap)/2.),100./number_ticks_colormap), np.arange(min_val, max_val+((float(max_val)-float(min_val))/number_ticks_colormap)/2., (float(max_val)-float(min_val))/number_ticks_colormap))
	plt.title("values")
	pdf.savefig()

