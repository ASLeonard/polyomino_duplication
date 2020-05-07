from scripts.interaction_dynamics import formTime

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.cm import get_cmap
from scipy.stats import expon
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy as np

def plotSingleHomology(data,L,ax,user_color='red',smooth_window=0):
    data = np.apply_along_axis(lambda a: a if np.sum(a)==0 else a/np.sum(a),1,data.astype(np.float))
    if smooth_window:
        cumulative_smooth = np.cumsum(data,0)
        data = cumulative_smooth[smooth_window-1::smooth_window]/smooth_window
        data = data[1:] - data[:-1]

    ##reverse homology so 100% is at top of plot
    pop_grid= ma.masked_equal(data[:80,::-1].T,0)
    cm = LinearSegmentedColormap.from_list("", ['gainsboro',user_color])
    ax.pcolormesh(pop_grid,cmap=cm,norm=LogNorm(.005,.25),rasterized=False)

def plotHomologyEvolution(run,L,edges,pathway='',smooth_window=0):
    raw_data = np.fromfile(f'{pathway}Homology_Run{run}.txt',dtype=np.uint16).reshape((-1,L+1))
    data = {pair:raw_data[pair[0]*4+pair[1]::16] for pair in combinations_with_replacement(range(4),2)}
            
    _, axes= plt.subplots(len(edges),1,sharex=True,sharey=True)
    colors= ['forestgreen','royalblue','orangered','darkviolet','gold']*5

    for i,edge in enumerate(edges):
        plotSingleHomology(data[edge],L=L,ax=axes[i],user_color=colors[i],smooth_window=smooth_window)
        axes[i].axis('off')
    plt.show(block=False)

def plotComplexityEvolution(data,renormalise_t0=False,interaction_sets=(10,),N_simulations=None):

    def deWeight(times,N_int=3,N_sim=None):
        upper_limit = len(times)/N_sim if N_sim else N_int

        heights = list(np.linspace(1/len(times),upper_limit,len(times)))
        unique_times = []
        for index,time in enumerate(times):
            if time not in unique_times:
                unique_times.append(time)
            else:
                heights.pop(len(times)-index-1)
        return unique_times[::-1], np.array(heights)

    _, axes = plt.subplots(3,2,sharex=True)
    ax2 = axes.flatten()[-1]
    for index, c,ax in zip(range(0,len(data),2),('g','c','purple','r','brown'),axes.flatten()):
        spline_data = []
        time_window = None

        for sub_I in (1,0):
            raw_times = np.array(sum((data[index+sub_I].discov_times[slice_] for slice_ in interaction_sets),[]))
            time_points, complexity_avg = deWeight(sorted(raw_times,reverse=True),len(interaction_sets),N_sim=N_simulations)
        
            time_samples = np.logspace(0,np.log10(max(time_points)),100)
            ax.plot(time_points, complexity_avg, ls='--' if sub_I else '-', c=c, label=f'{data[index+sub_I].L},{data[index+sub_I].dup_rate}',alpha=1)

            ## generate spline information for complexity gap
            complexity_spline = InterpolatedUnivariateSpline(time_points, complexity_avg)

            spline_data.append(complexity_spline)

        time_window = np.argmax(spline_data[0](time_samples)>.95*max(complexity_avg))
        new_start=formTime(data[index+sub_I].L//2,data[index+sub_I].S_c)/75
        if renormalise_t0:
            time_window = time_window[time_window>=new_start] - new_start
        else:
            ax2.axvline(new_start,c=c)

        complexity_gap = spline_data[0](time_samples[:time_window])-spline_data[1](time_samples[:time_window])
        max_gap = np.argmax(complexity_gap)

        ax.plot([time_samples[max_gap]]*2,[spline_data[0](time_samples[max_gap]),spline_data[1](time_samples[max_gap])],mfc='none',mec=c,marker='o',ms=12,c=c)
        time_samples=time_samples[:time_window]
        ax2.plot(time_samples,spline_data[0](time_samples)-spline_data[1](time_samples), c=c)

        ax.set_xscale('log')
        ax2.set_xscale('log')
        ax.legend()
    plt.show(block=False)


def plotTimex(*datas,fit_func=expon,renormalise=True,full_renorm=False,row2=False):
    if full_renorm:
        assert renormalise, 'conflicting settings'
    N = 2
    _, ax = plt.subplots(N)
    pop, markers = 100, {60:'o',80:'s',100:'d',120:'X',140:'*'}

    asyms = [formTime(data.L,data.S_c)/formTime(data.L//2,data.S_c) for data in datas]
    scaler = LogNorm(min(asyms),max(asyms)*2)

    for stage in range(N):
        for asym_val, data in zip(asyms,datas):
            L_scaler = 2-stage
            ##all homomers
            if stage == 0:
                mutate_rate_adjust = .33
                combinatoric_adjust = 4
            ##heteromers
            else:
                mutate_rate_adjust = 3/4
                combinatoric_adjust = 1

            gamma_factor = 1
            if full_renorm:
                if stage == 1 and data.dup_rate > 0:
                    gamma_factor = getGammas()[data.L]/100
                    L_scaler = 2

            raw_data = np.array(data.discov_times[10 if stage == 0 else 1]+([] if (stage == 0 or not row2) else data.discov_times[2]))
            renorm_factor = (formTime(data.L//L_scaler,data.S_c) / pop / mutate_rate_adjust / gamma_factor * combinatoric_adjust) if renormalise else 1
            data_scaled = raw_data/renorm_factor

            cmap = get_cmap('Blues_r' if data.dup_rate>0 else 'Oranges_r')
            asym_color = cmap(scaler(asym_val))

            if fit_func:
                fit_p = fit_func.fit(data_scaled,floc=0)
                x_points = np.linspace(0,max(data_scaled),300)
                ax[stage].plot(x_points,fit_func(*fit_p).pdf(x_points),marker=markers[data.L],markevery=[0,-1],c=asym_color,ls='-' if data.dup_rate ==0 else '--',lw=4,mew=5,mfc='none',ms=30,alpha=1,label=f'L:{data.L}, S_c:{data.S_c}' if data.dup_rate !=-1 else None)

        ax[stage].set_yscale('log',nonposy='clip')
    ax[0].legend()
    ax[0].tick_params(axis='both', which='major', labelsize=16)
    ax[1].tick_params(axis='both', which='major', labelsize=16)
    plt.show(block=False)

def printTransitionTable(data):
    for stage in range(max(data.discov_types['stage'])+1):
        print('Incoming at stage ',stage)
        for t_class in ('homodimeric','heterodimeric'):
            print(f'\t{t_class}: ',sum((data.discov_types['class']==t_class) & (data.discov_types['stage']==stage)))
        print('Compositions')
        for t_class in ('homodimeric','Du-Sp','heterodimeric'):
            print(f'\t{t_class}: ',sum((data.composition_types['class']==t_class) & (data.composition_types['stage']==stage)))
    return [sum((data.composition_types['class']==T) & (data.composition_types['stage']==2)) for T in ('Du-Sp','heterodimeric')]
