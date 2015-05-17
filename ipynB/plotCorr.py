'''

ku_iv = ivDf['log ku']
kf_iv = ivDf['log kf']

# Data from kinetic DB and Ivankov
pf = np.polyfit(kf_iv,ku_iv,1)
yfit = np.polyval(pf,kf_iv)

(cr,p) = pearsonr(ku_iv,kf_iv)

# data from Brioom etal paper,2015
ku_br = BrDf['log ku']
kf_br = BrDf['log kf']

bf = np.polyfit(kf_br,ku_br,1)
yfit_bf = np.polyval(bf,kf_br)

(br,pr) = pearsonr(ku_br,kf_br)

plt.figure(figsize=(10,3))
ax1 = plt.subplot(1,2,1)
ax1.plot(ivDf[ivDf["class"]=="a"]["log kf"],ivDf[ivDf["class"]=="a"]["log ku"],'bo')
ax1.plot(ivDf[ivDf["class"]=="b"]["log kf"],ivDf[ivDf["class"]=="b"]["log ku"],'go')
ax1.plot(ivDf[ivDf["class"]=="ab"]["log kf"],ivDf[ivDf["class"]=="ab"]["log ku"],'ks')

ax1.plot(kf,yfit,'r-')
ax1.set_xlabel('log k_{f}')
ax1.set_ylabel('log k_{u}')
ax1.set_title("Data from KineticDB(=%d)"%len(ku_iv))
ax1.text(-1,4,"log ku= %0.3f+%0.3f*log kf"%(pf[1],pf[0]))
ax1.text(-1,2,"R2=%0.3f"%(cr))


ax2 = plt.subplot(1,2,2)
ax2.plot(BrDf[BrDf["class"]=="a"]["log kf"],BrDf[BrDf["class"]=="a"]["log ku"],'bo')
ax2.plot(BrDf[BrDf["class"]=="b"]["log kf"],BrDf[BrDf["class"]=="b"]["log ku"],'go')
ax2.plot(BrDf[BrDf["class"]=="ab"]["log kf"],BrDf[BrDf["class"]=="ab"]["log ku"],'ks')
ax2.plot(kf_br,yfit_bf,'r-')

ax2.set_xlabel('log k_{f}')
ax2.set_ylabel('log k_{u}')
ax2.set_title("Data from protein Paper (=%d)"%len(ku_br))
ax2.text(-2,4,"log ku= %0.3f+%0.3f*log kf"%(bf[1],bf[0]))
ax2.text(-2,2,"R2=%0.3f"%(br))
plt.show()

#print cr,p
'''
