import os

# Directories
DIRM=os.path.join(os.getenv("HOME"),"Models/Meshes/3D/Helheim/"+MESHNAME+"/")

# Mesh nodes file
nodes_file=DIRM+"/elmer/mesh.nodes"
nodes=np.genfromtxt(nodes_file,delimiter=' ')

# Kristin's model results
kristin_file=os.path.join(os.getenv("HOME"),"Data/Climate/IceTemperature/Helheim/helheim_TA.xyz")
tempdata=np.genfromtxt(kristin_file,delimiter=',')
tempdata=np.delete(tempdata,(0),axis=0)

# Interpolate to nodes
Anodes_lin = griddata(tempdata[:,0:3],tempdata[:,4],nodes[:,2:5],method='linear')
Anodes_near = griddata(tempdata[:,0:3],tempdata[:,4],nodes[:,2:5],method='nearest')

nans=np.isnan(Anodes_lin)
Anodes_lin[nans]=Anodes_near[nans]

# Write out flow law parameter
fid = open(Inputs+"flowparameter.dat","w")
fid.write('{}\n'.format(len(Anodes_lin)))
for i in range(0,len(Anodes_lin)):
  fid.write('{0} {1} {2} {3} {4}\n'.format(int(nodes[i,0]),nodes[i,2],nodes[i,3],nodes[i,4],Anodes_lin[i]))
fid.close() 
del nan, kristin_file, Anodes_near, fid   


y,inds=np.unique(nodes[:,3],return_index=True)
x=nodes[inds,2]
Eta_s=np.zeros_like(y)
for i in range(0,len(y)):
  colinds=np.where(nodes[:,3]==y[i])
  colinds=np.array(colinds)
  surfind=np.argmin(nodes[colinds,4])
  Eta_s[i]=np.mean(Anodes_lin[colinds[0,:]])


pts=[(308482.09293552686, -2576367.5519362916),
 (303157.50809034629, -2575462.482329594)]

for i in range(0,len(tempdata)):
  if 
 
xdata,xinds=np.unique(tempdata[:,0],return_index=True)
ydata,yinds=np.unique(tempdata[:,1],return_index=True)
inds=np.unique(np.hstack([xinds,yinds]))
del xdata, ydata
xdata=tempdata[inds,0]
ydata=tempdata[inds,1] 


Eta_s_orig=np.zeros_like(xdata)
for i in range(0,len(xdata)):
  colinds=np.where(tempdata[:,0]==xdata[i])
  colinds=np.array(colinds)
  surfind=np.argmax(tempdata[colinds,2])
  Eta_s_orig[i]=tempdata[colinds[0,surfind],4] 
	
  