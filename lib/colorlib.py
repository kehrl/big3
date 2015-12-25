"""
Colour-blind proof distinct colours module, based on work by Paul Tol
Pieter van der Meer, 2011
SRON - Netherlands Institute for Space Research
"""

# colour table in HTML hex format
hexcols = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', 
           '#CC6677', '#882255', '#AA4499', '#661100', '#6699CC', '#AA4466',
           '#4477AA']

greysafecols = ['#809BC8', '#FF6666', '#FFCC66', '#64C204']

xarr = [[12], 
        [12, 6], 
        [12, 6, 5], 
        [12, 6, 5, 3], 
        [0, 1, 3, 5, 6], 
        [0, 1, 3, 5, 6, 8], 
        [0, 1, 2, 3, 5, 6, 8], 
        [0, 1, 2, 3, 4, 5, 6, 8], 
        [0, 1, 2, 3, 4, 5, 6, 7, 8], 
        [0, 1, 2, 3, 4, 5, 9, 6, 7, 8], 
        [0, 10, 1, 2, 3, 4, 5, 9, 6, 7, 8], 
        [0, 10, 1, 2, 3, 4, 5, 9, 6, 11, 7, 8]]

# get specified nr of distinct colours in HTML hex format.
# in: nr - number of colours [1..12]
# returns: list of distinct colours in HTML hex
def get_distinct(nr):

    #
    # check if nr is in correct range
    #
    
    if nr < 1 or nr > 12:
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