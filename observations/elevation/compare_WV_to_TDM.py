# We ordered two TDM DEMs (one for Kanger and Helheim) for dates when we have Worldview 
# DEMs to compare the two and calculate uncertainties. For Kanger, the Worldview/TDM pair
# was collected on 2011-07-08. For Helheim, the pair was collected on 2011-06-15.

import os
impmort geotifflib

os.chdir(os.path.join(os.getenv("HOME"),"Bigtmp"))

file_wv_H = os.path.join(os.getenv("DATA_HOME"),"Elevation/Worldview/Helheim/20110615_1442_102001001338D000_1020010013B10800-DEM_32m_trans.tif")
file_td_H = os.path.join(os.getenv("DATA_HOME"),"Elevation/TDM/Helheim/20110615_1957_1032905_DEM_32m_trans.tif")
os.system("geodiff "+file_td_H+" "+file_wv_H+" -o helheim")

x_H,y_H,zdiff_H = geotifflib.read("helheim-diff.tif")
xmask_H,ymask_H,zmask_H = masklib.load_grid("Helheim",x_H[0],x_H[-1],y_H[0],y_H[-1],x_H[1]-x_H[0],icefront_time=datelib.date_to_fracyear(2011,6,15))
ind = np.where(zmask_H == 1)
zdiff_masked_H = np.array(zdiff_H)
zdiff_masked_H[ind] = float('nan')
nonnan = np.where(~(np.isnan(zdiff_masked_H.flatten())))[0]
zdiff_masked_H = zdiff_masked_H.flatten()[nonnan]

file_wv_K = os.path.join(os.getenv("DATA_HOME"),"Elevation/Worldview/Kanger/20110708_1403_1020010015567100_10200100146CB000-DEM_32m_trans.tif")
file_td_K = os.path.join(os.getenv("DATA_HOME"),"Elevation/TDM/Kanger/20110708_1940_1032883_DEM_mask_32m_trans.tif")

os.system("geodiff "+file_td_K+" "+file_wv_K+" -o kanger")

x_K,y_K,zdiff_K = geotifflib.read("kanger-diff.tif")
xmask_K,ymask_K,zmask_K = masklib.load_grid("Kanger",x_K[0],x_K[-1],y_K[0],y_K[-1],x_K[1]-x_K[0],icefront_time=datelib.date_to_fracyear(2011,7,8))
ind = np.where(zmask_K == 1)
zdiff_masked_K = np.array(zdiff_K)
zdiff_masked_K[ind] = float('nan')
nonnan = np.where(~(np.isnan(zdiff_masked_K.flatten())))[0]
zdiff_masked_K = zdiff_masked_K.flatten()[nonnan]


plt.figure(figsize=(3,3))
plt.plot([0,0],[0,0.3],'k:')
n,bins,patches=plt.hist(zdiff_masked_H,bins=40,range=[-10,10],normed=True,histtype='step',color='k',lw=1.5,label='Helheim')
plt.hist(zdiff_masked_K,bins=40,range=[-10,10],normed=True,histtype='step',color='r',lw=1.5,label='Kanger')
plt.xlabel('TDM - WV Elevation difference (m)',fontsize=9,fontname='Arial')
plt.xticks(np.arange(-10,15,5),fontsize=8,fontname='Arial')
plt.yticks(np.arange(0,0.4,0.1),fontsize=8,fontname='Arial')
plt.ylabel('Normalized Frequency',fontname='Arial',fontsize=9)
plt.legend(fontsize=9,numpoints=1,handlelength=0.4,handletextpad=0.4)
plt.tight_layout()
plt.savefig(os.path.join(os.getenv("HOME"),"Bigtmp/TDM_WV_comparison.pdf"),format="PDF")

n,bins,patches=plt.hist(abs(zdiff_masked_H),bins=5000,normed=True,cumulative=True,color='k',lw=1.5,label='Helheim')
ind = np.argmin(abs(n-0.95))
print "Helheim: \n","mean is",np.mean(zdiff_masked_H),"m \n","95% is within",bins[ind],"m \n","RMSE is",np.sqrt(np.mean(zdiff_masked_H**2)),"m"

n,bins,patches=plt.hist(abs(zdiff_masked_K),bins=5000,normed=True,cumulative=True,color='k',lw=1.5,label='Kanger')
ind = np.argmin(abs(n-0.95))
print "Helheim: \n","mean is",np.mean(zdiff_masked_K),"m \n","95% is within",bins[ind],"m \n","RMSE is",np.sqrt(np.mean([zdiff_masked_K**2)),"m"
