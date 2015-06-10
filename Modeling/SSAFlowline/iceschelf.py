B = A**(-1/n)
C = (rho_i*g*(rho_sw-rho_i)/(4*B*rho_sw))**(n)

H_o = 500
U_o = 200/yearinsec

H_ana = (((n+1)*C/(H_o*U_o))*input_x+H_o**(-(n+1)))**(-1/(n+1))

# Plot model inputs
fig = plt.figure(figsize=(5,4))

plt.subplot(211)
plt.plot(input_x/1e3,H_ana,'k',linewidth=3,label='Analytical')
plt.plot(norm_nodes[0:-1:5,0]/1e3,H[0:-1:5,0],'b.-',label='Year 0')
plt.plot(norm_nodes[0:-1:5,183]/1e3,H[0:-1:5,183],'c.-',label='Year 5')
plt.plot(norm_nodes[0:-1:5,183*4]/1e3,H[0:-1:5,183*4],'g.-',label='Year 20')
plt.plot(norm_nodes[0:-1:5,183*10]/1e3,H[0:-1:5,183*10],'r.-',label='Year 50')
plt.xticks(np.arange(0,21,2))
plt.yticks([200,300,400,500],fontsize=13)
plt.ylabel('Ice thickness (m)',fontsize=13)
plt.gca().set_xticklabels([])
plt.xlim([0,21])
plt.text(1,200,'a',fontsize=13,fontweight='bold')
plt.gca().set_xticklabels([])

plt.subplot(212)
plt.plot(stag_nodes[0:-1:5,-2]/1e3,u_stag[0:-1:5,-2]*yearinsec/1e3,'r.-')
plt.ylim([0,0.6])
plt.yticks(np.arange(0,0.6,0.1),fontsize=13)
plt.xticks(np.arange(0,21,2))
plt.xlim([0,21])
plt.ylabel('Velocity (km/yr)',fontsize=13)
plt.text(1,0.1,'b',fontsize=13,fontweight='bold')

plt.subplot(211)
plt.legend(bbox_to_anchor=(1.0, 0.65),
           bbox_transform=plt.gcf().transFigure,fontsize=10)

plt.subplots_adjust(left=None, bottom=None, right=0.73, top=None, wspace=0.05, hspace=0.05)
