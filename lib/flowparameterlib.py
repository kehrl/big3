# This function takes a temperature and outputs the flow parameter 
# following Cuffey and Paterson 2010, pg. 72.

def arrhenius(T):
  import numpy as np
  
  A=np.zeros_like(T)
  for i in range(0,len(T)):
    if (T[i] > 273.15):
      A[i]=2.4E-24
    elif (T[i] > 263.15):
      A[i]=3.5E-25* np.exp(-115.0E03/8.314 * (1/T[i] - 1/263.15))
    else:
      A[i]=3.5E-25 * np.exp(-60.0E03/8.314 * (1/T[i] - 1/263.15))
    	
  return A 
