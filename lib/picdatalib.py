from urlparse import urlparse
 
import pygtk
import gtk
 
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
 
import numpy as np
 
def tellme(s):
    print s
    plt.title(s,fontsize=16)
    plt.draw()
 
def pic2data(source=file,straight=True,origin='upper'):
    """ GUI to get data from a XY graph image. Either provide the graph
        as a path to an image in 'source' or copy it to the clipboard.
    """
         
    image = mpimg.imread(source)
 
    ###### DISPLAY THE IMAGE
     
    plt.ion() # interactive mode !
    fig, ax = plt.subplots(1)
    imgplot = ax.imshow(image, origin=origin)
    fig.canvas.draw()
    plt.draw()
     
    ##### PROMPT THE AXES
     
    def promptPoint(text=None):
         
      if text is not None: tellme(text)
      return  np.array(plt.ginput(1,timeout=-1)[0])
        
    def askValue(text):
      value = float(raw_input(text+': '))
      return value
     
    origin = promptPoint('Place the origin')
    origin_value = askValue('X origin'),askValue('Y origin')
                                          
    Xref =  promptPoint('Place the X reference')
    Xref_value = askValue('X reference')
     
    Yref =  promptPoint('Place the Y reference')
    Yref_value = float(askValue('Y reference'))
     
    if straight :
         
        Xref[1] = origin[1]
        Yref[0] = origin[0]
     
    ##### PROMPT THE POINTS
     
    selected_points = []
     
    tellme("Select your points !")
    print "Right-click or press 's' to select"
    print "Left-click or press 'del' to deselect"
    print "Middle-click or press 'Enter' to confirm"
    print "Note that the keyboard may not work."
     
    selected_points = plt.ginput(-1,timeout=-1)
     
    ##### RETURN THE POINTS COORDINATES
     
    #~ selected_points.sort() # sorts the points in increasing x order
     
    # compute the coordinates of the points in the user-defined system
     
    OXref = Xref - origin
    OYref = Yref - origin
    xScale =  (Xref_value - origin_value[0]) / np.linalg.norm(OXref)
    yScale =  (Yref_value - origin_value[1]) / np.linalg.norm(OYref)
     
    ux = OXref / np.linalg.norm(OXref)
    uy = OYref / np.linalg.norm(OYref)
     
    result = [(ux.dot(pt - origin) * xScale + origin_value[0],
               uy.dot(pt - origin) * yScale + origin_value[1])
               for pt in selected_points ]
    
     
    plt.ioff()
     
    return result