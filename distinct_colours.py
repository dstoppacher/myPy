# -*- coding: iso-8859-1 -*-

"""
Colour-blind proof distinct colours module, based on work by Paul Tol
Pieter van der Meer, 2011
SRON - Netherlands Institute for Space Research
"""

# colour table in HTML hex format the last two one are grey
hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', 
           '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#225588',
           '#4477AA', '#333333', '#4c4c4c', '#191919', '#AA4466', '#ffab00', '#d2691e', '#000000']

greysafecols = ['#809BC8', '#FF6666', '#FFCC66', '#64C204']

#original Galacticus blue: AA4466, new Galacticus Blue: 225588


xarr = [[11],
	#[6], 
	#[0], 
        #[11, 6],
		[8, 3], 
        #[11, 6, 5],
	#SF-History paper when clustering of 3 sub-samples: red, lowZ, highZ
        [6, 8, 3],  
        [11, 6, 5, 2],
        #[0, 2, 5, 6],
	#for SDSS bands u=blue, y=red 
        [0, 1, 3, 5, 6],
		#[11, 17, 6, 8, 3],
	#Galacticus CMASS Clustering color sample
	#[3, 4, 5, 6, 8], 
        #[12, 6, 5, 2, 0],
        #[0, 1, 2, 5, 6, 8],
        #[0, 1, 2, 3, 5, 6],
		#[0, 1, 2, 3, 17, 6],
		[0, 8, 3, 8, 3, 2], 
        [0, 0, 17, 6, 8, 3, 11],
        [0, 1, 17, 3, 2, 6, 19, 8],
	#for 9 subsamples in the assembly bias paper, customized by Doris 
	  	[0, 1, 2, 3, 5, 6, 12, 8, 9], 
        #[0, 1, 2, 3, 4, 5, 6, 7, 8],
	  	#[6, 6, 11, 11, 5, 5, 12, 8, 9], 
	#for SDSS bands double set (10 bands) u=blue, y=red 
        [0, 1, 3, 5, 6, 0, 1, 3, 5, 6], 
        #[0, 1, 2, 3, 4, 5, 9, 6, 7, 8], 
        [0, 10, 1, 2, 3, 4, 5, 9, 6, 7, 8],
	#12 redshifts 
        [0, 10, 1, 2, 3, 4, 5, 9, 6, 11, 7, 8],
	#13 redshifts 
        [0, 10, 1, 2, 3, 4, 5, 9, 6, 7, 13, 14, 15],
	#14 redshifts/catalogs
		[0, 0, 2, 17, 5, 6, 12, 8, 3, 13, 9, 18, 15, 4]]
		#[0, 4, 2, 1, 5, 6, 12, 8, 3, 13, 9, 18, 15, 17]]
		#[0, 1, 2, 4, 5, 6, 12, 8, 3, 13, 9, 18, 15, 17]]


# get specified nr of distinct colours in HTML hex format.
# in: nr - number of colours [1..12]
# returns: list of distinct colours in HTML hex
def get_distinct(nr):

    #
    # check if nr is in correct range
    #
    
    if nr < 1 or nr > 14:
        print "wrong nr of distinct colours!"
        return

    #
    # get list of indices
    #
    
    lst = xarr[nr-1]
    
    #
    # generate colour list by stepping through indices and looking them up
    # in the colour table
    #

    i_col = 0
    col = [0] * nr
    for idx in lst:
        col[i_col] = hexcols[idx]
        i_col+=1
    return col

# gets 4 colours, which also look distinct in black&white
# returns: list of 4 colours in 
#def get_distinct_grey():
    
# displays usage information and produces example plot.
if __name__ == '__main__':
    import numpy as np
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt

    print __doc__
    print "usage examples: "
    print "print distinct_colours.get_distinct(2)"
    print get_distinct(2)
    print "print distinct_colours.greysafecols"
    print greysafecols

    print "\ngenerating example plot: distinct_colours_example.png"
    plt.close()
    t = np.arange(0.0, 2.0, 0.01)
    s = np.sin(2*np.pi*t)
    c = np.cos(2*np.pi*t)
    cols = get_distinct(2)
    plt.plot(t, s, linewidth=1.0, c=cols[0])
    plt.plot(t, c, linewidth=1.0, c=cols[1])

    plt.xlabel('time (s)')
    plt.ylabel('voltage (mV)')
    plt.title('Distinct colours example')
    plt.grid(True)
    plt.show()
    plt.savefig("distinct_colours_example.png")
