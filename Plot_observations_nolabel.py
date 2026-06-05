import numpy

facecolor = 'white'
edge_color = 'black'
markeredgewidth = 2

def plot_obs_masses_nolabel(obj):
        obj.errorbar([10.6],[10**6.2],yerr=[[10**6.2-10**5.9],[10**6.4-10**6.2]],marker='s',color='black',ms=10,markerfacecolor=facecolor,markeredgecolor=edge_color,markeredgewidth=markeredgewidth,zorder=1000)
        obj.errorbar([8.7],[10**6.95],yerr=[[10**6.95-10**6.45],[10**7.35-10**6.95]],marker='o',color='black',ms=10,markerfacecolor=facecolor,markeredgecolor=edge_color,markeredgewidth=markeredgewidth,zorder=1000)
        min_U,max_U,mid_U = (8-3.2)*1e7, (8+3.7)*1e7,8e7
        obj.errorbar([10],[mid_U],yerr=[[mid_U-min_U],[max_U-mid_U]],marker='*',color='black',ms=15,markerfacecolor=facecolor,markeredgecolor=edge_color,markeredgewidth=markeredgewidth,zorder=1000)
        min_U,max_U,mid_U = 1e7, 1e8,10**7.5
        obj.errorbar([10.3],[mid_U],yerr=[[mid_U-min_U],[max_U-mid_U]],marker='^',color='black',ms=10,markerfacecolor=facecolor,markeredgecolor=edge_color,markeredgewidth=markeredgewidth,zorder=1000)
        min_U,max_U,mid_U = 10**6.65, 10**8.5,10**((6.65+8.5)/2)
        obj.errorbar([9.22],[mid_U],yerr=[[mid_U-min_U],[max_U-mid_U]],marker='v',color='black',ms=10,markerfacecolor=facecolor,markeredgecolor=edge_color,markeredgewidth=markeredgewidth,zorder=1000)
        data=numpy.loadtxt('/home/whm4dg/observational_data/z6_QSO.txt')
        mask = data[:,1] > 1e9
        obj.errorbar(data[:,0][mask], data[:,1][mask], color='grey', alpha=0.4, linewidth=0, ms=5, marker='o', zorder=0, label=r'$z\sim6-7$ QSO')    
        
