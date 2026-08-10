# import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from pathlib import Path
from .utilities import rt_quasi_deuteron
from .presets import *

# plot lists
EXP_QVPLOT_LIST = ['Yamaguchi', 'Barreau', 'Jourdan', 'Goldemberg', 'Buki', 'Photo-production']
EXP_Q2PLOT_LIST = ['Yamaguchi', 'Baran', 'Sheren', 'Photo-production']
THEORY_QVPLOT_LIST = ['SuSAv2', 'GFMC','ED-RMF', 'STA-QMC', 'CFG']
THEORY_Q2PLOT_LIST = ['SuSAv2', 'ED-RMF']
MC_QVPLOT_LIST = ['NuWro-SF', 'NuWro-SF-FSI', 'ACHILLES']
MC_Q2PLOT_LIST = ['NuWro-SF', 'NuWro-SF-FSI']

# scatter plotting styles
SCATTER_STYLES = {
    'Yamaguchi':{'color':'cyan', 's':10, 'marker':'o', 'edgecolors':'gray', 'linewidth':0.3, 'zorder':2},
    'Photo-production':{'s':30,'color':'lime','marker':'^','edgecolors':'gray','linewidth':0.5, 'zorder':2},
    'Goldemberg':{'color':'blue','marker':'d','s':15, 'edgecolors':'gray','linewidth':0.5, 'zorder':2},
    'Jourdan':{'color':'limegreen','marker':'v','s':15, 'edgecolors':'gray','linewidth':0.5, 'zorder':2},
    'Barreau':{'color':'yellow', 'marker':'H','s':15, 'edgecolors':'gray','linewidth':0.5, 'zorder':2},
    'Buki':{'color':'blue','marker':'P','s':15, 'edgecolors':'gray','linewidth':0.5, 'zorder':2},
    'Sheren':{'color':'deepskyblue','s':15, 'edgecolors':'gray','linewidth':0.5, 'zorder':2},
    'Baran':{'color':'orange','s':15, 'edgecolors':'gray','linewidth':0.5, 'zorder':2}
}

# error bar plotting styles
ERRORBAR_STYLES = {
    'Goldemberg':{'color':'skyblue','ecolor':'blue','markersize':0, 'capsize':0, 'lw':1, 'fmt':'D', 'elinewidth':1, 'alpha':0.5, 'zorder':1},
    'Yamaguchi':{'color':'gray','markersize':0, 'capsize':0, 'lw':0.5, 'elinewidth':1, 'alpha':0.5, 'zorder':1},
    'Jourdan':{'color':'limegreen','markersize':0, 'capsize':0, 'lw':1, 'fmt':'D', 'elinewidth':1, 'alpha':0.5, 'zorder':1},
    'Barreau':{'color':'black','markersize':0, 'capsize':0, 'lw':1, 'fmt':'D', 'elinewidth':1, 'alpha':0.5, 'zorder':1},
    'Buki':{'color':'blue','markersize':0, 'capsize':0, 'lw':1, 'fmt':'D', 'elinewidth':1, 'alpha':0.5, 'zorder':1},
    'Sheren':{'ecolor':'deepskyblue','markersize':0, 'capsize':0, 'lw':1, 'fmt':'D', 'elinewidth':1, 'alpha':0.5, 'zorder':1},
    'Baran':{'ecolor':'gray','markersize':0, 'capsize':0, 'lw':1, 'fmt':'D', 'elinewidth':1, 'alpha':0.5, 'zorder':1},
    'Photo-production':{'color':'lime','markersize':0, 'capsize':0, 'lw':1, 'fmt':'D', 'elinewidth':1, 'alpha':0.5, 'zorder':1}
}

# line plotting styles
LINE_STYLES = {
    'GFMC':{'color':'violet', 'linestyle':'solid','lw':2,'zorder':-1},
    'ED-RMF':{'color':'cornflowerblue', 'linestyle':'solid','lw':3, 'zorder':-1},
    'CFG':{'color':'brown','linestyle':'-','lw':2, 'zorder':-1},
    'STA-QMC':{'color':'lime', 'linestyle':'solid','lw':2, 'zorder':-1},
    'SuSAv2':{'zorder':-1},
    'NuWro-SF':{'color':'violet','lw':2,'linestyle':'-', 'zorder':-1},
    'NuWro-SF-FSI':{'color':'green','lw':2,'linestyle':'-', 'zorder':-1},
    'ACHILLES':{'color':'blue','linestyle':'dotted','lw':2.5, 'zorder':-1}
}

# plot vertical dash lines at invariant mass w = 0.93, 1.07, 1.23 GeV/c^2
W_PLOT_LIST = [0.93, 1.07, 1.23]
W_PLOT_COLORS = {0.93:'darkorange', 1.07:'turquoise', 1.23:'slateblue'}

# plot hights
RLRT_QVPLOT_HEIGHTS = {
    0.1:[50.0,12],
    0.148:[90.0,27],
    0.167:[80.0,30],
    0.205:[80.0,35],
    0.24:[60.0,41],
    0.3:[50.0,44],
    0.38:[30.0,35],
    0.475:[17.5,45],
    0.57:[10.0,35],
    0.649:[10.0,40],
    0.756:[10.0,42],
    0.991:[10.0,25],
    1.619:[100.0,45],
    1.921:[120.0,35],
    2.213:[80.0,35],
    2.5:[50.0,26],
    2.783:[40.0,15],
    3.5:[12.5,12]
}

RLRT_QVINSET_HEIGHTS = {
    1.619:[0.6, 2.5],
    1.921:[0.3, 1.3],
    2.213:[0.15, 0.8],
    2.5:[0.1, 0.5],
    2.783:[0.06, 0.3],
    3.5:[0.013, 0.07]
}

RLRT_QVINSET_XLIMS = {
    1.619:[0.86, 1.23],
    1.921:[1.10, 1.50],
    2.213:[1.34, 1.78],
    2.5:[1.58, 2.06],
    2.783:[1.82, 2.34],
    3.5:[2.44, 3.08]
}

RLRT_Q2PLOT_HEIGHTS = {
    0.01:[80,20],
    0.02:[65,33],
    0.026:[80,35],
    0.04:[70,35],
    0.056:[60,41],
    0.093:[40,40],
    0.12:[32,35],
    0.16:[22,30],
    0.265:[100,30],
    0.38:[60,25],
    0.5:[50,25],
    0.8:[30,20],
    1.25:[15,15],
    1.75:[12,9],
    2.25:[8,8],
    2.75:[4,6],
    3.25:[3,4],
    3.75:[2,3]
}

RLRT_Q2INSET_HEIGHTS = {
    0.01:[20.0,13.0],
    0.02:[30.0,27.0],
    0.026:[40.0,30.0],
    0.04:[500.0,35.0],
    0.056:[500.0,41.0],
    0.093:[20.0,44.0],
    0.12:[20.0,45.0],
    0.16:[80.0,45.0],
    0.265:[11.0,49.0],
    0.38:[6.0,30.0],
    0.5:[5.0,30.0],
    0.8:[2.0,15.0],
    1.25:[0.6,5.0],
    1.75:[0.2,3.0],
    2.25:[0.2,2.0],
    2.75:[0.05,1.0],
    3.25:[0.03,0.7],
    3.75:[0.02,0.22]
}

RLRT_Q2INSET_XLIMS = {
    0.265:[0,0.7],
    0.38:[0.1,0.67],
    0.5:[0.05,0.8],
    0.8:[0.2,0.9],
    1.25:[0.4,1.2],
    1.75:[0.7,1.5],
    2.25:[0.9,1.8],
    2.75:[1.2,2.0],
    3.25:[1.4,2.4],
    3.75:[1.7,2.7]
}

def plot_response_qvbin(df_this_analysis : pd.DataFrame, qvcenters : list[float] = [0.3, 0.38, 0.57], figsize_per_row : tuple[float, float] = (11, 1.8),
        sharex : bool = False, figshow : bool = False, theory_plot_list : list[str] = THEORY_QVPLOT_LIST, sheet_CBfit : pd.DataFrame = None,
        exp_plot_list : list[str] = EXP_QVPLOT_LIST, mc_plot_list : list[str] = MC_QVPLOT_LIST) -> Figure:
    # TODO: add comments

    fig, axs = plt.subplots(nrows=len(qvcenters), ncols=2, figsize=(figsize_per_row[0], figsize_per_row[1] * len(qvcenters) + 1.5), dpi=200, sharex=sharex, squeeze=False)
    handles = []
    labels = []

    repo_path = Path(__file__).resolve().parent.parent
    plot_xlsx = repo_path/'Carbon/one_sheet_to_rule_them_all.xlsx'
    if sheet_CBfit is None:
        sheet_CBfit = pd.read_excel(plot_xlsx, sheet_name='CBfit_qvbin')
    sheet_mc_rl = pd.read_excel(plot_xlsx, sheet_name='mc_rl_qvbin')
    sheet_mc_rt = pd.read_excel(plot_xlsx, sheet_name='mc_rt_qvbin')
    sheet_exp_rl = pd.read_excel(plot_xlsx, sheet_name='exp_rl_qvbin')
    sheet_exp_rt = pd.read_excel(plot_xlsx, sheet_name='exp_rt_qvbin')
    sheet_theory_rl = pd.read_excel(plot_xlsx, sheet_name='theory_rl_qvbin')
    sheet_theory_rt = pd.read_excel(plot_xlsx, sheet_name='theory_rt_qvbin')

    # helper function to plot on given axes
    def plot_on_rlrt_subplots(ax_rl, ax_rt, qvcenter):
        # plot vertical dash lines at w = 0.93, 1.07, 1.23 GeV/c^2
        for w in W_PLOT_LIST:
            nu_w = 0.025 - MASS_NUCLEON + np.sqrt(qvcenter**2 + w**2)
            if nu_w < qvcenter:
                ax_rl.axvline(x = nu_w, color = W_PLOT_COLORS[w], linestyle = 'dashdot', lw = 1, label = f'$W={w}$ GeV/c$^2$', zorder=-2)
                ax_rt.axvline(x = nu_w, color = W_PLOT_COLORS[w], linestyle = 'dashdot', lw = 1, zorder=-2)

        # plot q2 = 0 vertical dash line
        ax_rl.axvline(x=qvcenter, color = 'brown', linestyle='dashdot',lw=1,label=f'$Q^2$ = 0 (GeV/c)$^2$')
        ax_rt.axvline(x=qvcenter, color = 'brown', linestyle='dashdot',lw=1)

        # plot Christy Bodek fit
        # 25 Sep 25: shift inelastic to left by 18 MeV; 25 July 18: RT quasi deuteron added
        #  FIXME: 26 May 7 inelastic peak shift is removed for now. address later.
        ChristyBodekFit = sheet_CBfit.loc[sheet_CBfit['qv']==qvcenter].copy()
        ax_rl.plot(ChristyBodekFit['nu'],ChristyBodekFit["rltot"],color='black',label="$R_L$(total), $R_T$(total) Christy-Bodek-2024", linestyle='solid',lw=0.8, zorder=0)
        ax_rt.plot(ChristyBodekFit['nu'],ChristyBodekFit["rttot"], color='black', linestyle='solid',lw=0.8, zorder=0)
        ax_rl.plot(ChristyBodekFit['nu'],ChristyBodekFit["rlqe"], color='black',label="$R_L$(QE), $R_T$(QE+TE) Christy-Bodek-2024", linestyle='dotted',lw=0.8, zorder=0)
        ax_rt.plot(ChristyBodekFit['nu'],ChristyBodekFit["rtqe"] + ChristyBodekFit["rte"], color='black', linestyle='dotted',lw=0.8, zorder=0)

        # plot mc as lines
        for mc in mc_plot_list:
            df_rl = sheet_mc_rl.loc[(sheet_mc_rl['qv']==qvcenter) & (sheet_mc_rl['mc']==mc)]
            df_rt = sheet_mc_rt.loc[(sheet_mc_rt['qv']==qvcenter) & (sheet_mc_rt['mc']==mc)]
            if len(df_rl) > 0:
                ax_rl.plot(df_rl['nu'], df_rl['rl'], **LINE_STYLES[mc], label=f'$R_L$, $R_T$ {mc}')
            if len(df_rt) > 0:
                ax_rt.plot(df_rt['nu'], df_rt['rt'], **LINE_STYLES[mc])

        # plot theory as lines
        for theory in theory_plot_list:
            df_rl = sheet_theory_rl.loc[(sheet_theory_rl['qv']==qvcenter) & (sheet_theory_rl['theory']==theory)]
            df_rt = sheet_theory_rt.loc[(sheet_theory_rt['qv']==qvcenter) & (sheet_theory_rt['theory']==theory)]
            if len(df_rl) > 0:
                ax_rl.plot(df_rl['nu'], df_rl['rl'], **LINE_STYLES[theory], label=f'$R_L$, $R_T$ {theory}')
            if len(df_rt) > 0:
                ax_rt.plot(df_rt['nu'], df_rt['rt'], **LINE_STYLES[theory])
        
        # plot experiment data as scatter and error bars
        for exp in exp_plot_list:
            df_rl = sheet_exp_rl.loc[(sheet_exp_rl['qv']==qvcenter) & (sheet_exp_rl['experiment']==exp)]
            df_rt = sheet_exp_rt.loc[(sheet_exp_rt['qv']==qvcenter) & (sheet_exp_rt['experiment']==exp)]
            if len(df_rl) > 0:
                ax_rl.errorbar(df_rl['nu'], df_rl['rl'], yerr=df_rl['rlerr'], **ERRORBAR_STYLES[exp])
                ax_rl.scatter(df_rl['nu'], df_rl['rl'], **SCATTER_STYLES[exp])
            if len(df_rt) > 0:
                ax_rt.errorbar(df_rt['nu'], df_rt['rt'], yerr=df_rt['rterr'], **ERRORBAR_STYLES[exp])
                ax_rt.scatter(df_rt['nu'], df_rt['rt'], **SCATTER_STYLES[exp], label=f'$R_L$, $R_T$ {exp}')

        # plot this analysis
        df_qvbin = df_this_analysis.loc[df_this_analysis['qvcenter'] == qvcenter]
        if qvcenter in [0.148,0.167,0.205,0.24,0.3]:
            # avoid overlapping with Yamaguchi
            nu_yam = sheet_exp_rl.loc[(sheet_exp_rl['qv']==qvcenter) & (sheet_exp_rl['experiment']=='Yamaguchi')]['nu'].max()
            df_qvbin = df_qvbin.loc[df_qvbin['nu'] >= nu_yam]
        ax_rl.errorbar(df_qvbin['nu'], df_qvbin['rl'], yerr = df_qvbin['rlerr'], ecolor='red', markersize=0, capsize=0, lw=1, fmt='D', elinewidth=1, alpha=0.5, zorder=2) 
        ax_rl.scatter(df_qvbin['nu'], df_qvbin['rl'], color='red', label = '$R_L$, $R_T$ this analysis', marker='D', s=15, edgecolors='gray', linewidth=0.5, zorder=2)
        ax_rt.errorbar(df_qvbin['nu'], df_qvbin['rt'], yerr = df_qvbin['rterr'], ecolor='red', markersize=0, capsize=0, lw=1, fmt='D', elinewidth=1, alpha=0.5, zorder=2)
        ax_rt.scatter(df_qvbin['nu'], df_qvbin['rt'], color='red', marker='D', s=15, edgecolors='gray', linewidth=0.5, zorder=2)

    # loop over given qv centers
    for i, qvcenter in enumerate(qvcenters):
        if qvcenter not in QVCENTERS:
            print(f'q={qvcenter} not in spreadsheet')
            continue

        plot_on_rlrt_subplots(axs[i, 0], axs[i, 1], qvcenter)

        # handle legend labels
        h, l = axs[i, 0].get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
        h, l = axs[i, 1].get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)

        # main axs styling
        axs[i, 0].text(0.95, 0.95,'$R_L$ ($|\mathbf{q}|$ = '+f'{qvcenter} GeV/c)',transform=axs[i, 0].transAxes,ha='right',va='top',color='gray')
        axs[i, 0].tick_params(which='both', direction='in', top=True, right=True)
        axs[i, 0].minorticks_on()
        axs[i, 0].set_ylabel(r'$R_L$ (GeV$^{-1}$)')
        axs[i, 0].set_xlim(0, qvcenter*1.05)
        axs[i, 0].set_ylim(0, RLRT_QVPLOT_HEIGHTS[qvcenter][0])
        axs[i, 1].text(0.95, 0.95,'$R_T$ ($|\mathbf{q}|$ = '+f'{qvcenter} GeV/c)',transform=axs[i, 1].transAxes,ha='right',va='top',color='gray')
        axs[i, 1].tick_params(which='both', direction='in', top=True, right=True)
        axs[i, 1].minorticks_on()
        axs[i, 1].set_ylabel(r'$R_T$ (GeV$^{-1}$)')
        axs[i, 1].set_xlim(0, qvcenter*1.05)
        axs[i, 1].set_ylim(0, RLRT_QVPLOT_HEIGHTS[qvcenter][1])

         # put inset plot at top left corner
        if qvcenter >= 1.619:
            ax_inset0 = inset_axes(axs[i, 0], width="45%", height="50%", loc='upper left')
            ax_inset1 = inset_axes(axs[i, 1], width="45%", height="50%", loc='upper left')

            # plot on inset plot
            plot_on_rlrt_subplots(ax_inset0, ax_inset1, qvcenter)

            # inset axs styling
            ax_inset0.yaxis.tick_right()
            ax_inset0.set_xlim(RLRT_QVINSET_XLIMS[qvcenter][0], RLRT_QVINSET_XLIMS[qvcenter][1])
            ax_inset0.set_ylim(0, RLRT_QVINSET_HEIGHTS[qvcenter][0])
            ax_inset0.tick_params(which='both', direction='in', top=True, right=True)
            ax_inset0.minorticks_on()
            ax_inset1.yaxis.tick_right()
            ax_inset1.set_xlim(RLRT_QVINSET_XLIMS[qvcenter][0], RLRT_QVINSET_XLIMS[qvcenter][1])
            ax_inset1.set_ylim(0, RLRT_QVINSET_HEIGHTS[qvcenter][1])
            ax_inset1.tick_params(which='both', direction='in', top=True, right=True)
            ax_inset1.minorticks_on()

    # add x-axis label only to the last subplot
    axs[-1, 0].set_xlabel(r'$\nu$ (GeV)')
    axs[-1, 1].set_xlabel(r'$\nu$ (GeV)')

    # find unique legend labels, put it at bottom
    unique = dict(zip(labels, handles))
    handles = list(unique.values())
    labels = list(unique.keys())
    if len(qvcenters) == 3:
        fig.legend(handles, labels, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.08), frameon=False)
    else:
        fig.legend(handles, labels, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.02), frameon=False)
    # fig.tight_layout()
    if figshow:
        plt.show()
    plt.close()

    return fig

def plot_response_q2bin(df_this_analysis : pd.DataFrame, q2centers : list[float] = [0.093, 0.12, 0.16], figsize_per_row : tuple[float, float] = (11, 1.8),
        sharex : bool = False, figshow : bool = False, theory_plot_list : list[str] = THEORY_Q2PLOT_LIST, sheet_CBfit : pd.DataFrame = None,
        exp_plot_list : list[str] = EXP_Q2PLOT_LIST, mc_plot_list : list[str] = MC_Q2PLOT_LIST) -> Figure:
    
    fig, axs = plt.subplots(nrows=len(q2centers), ncols=2, figsize=(figsize_per_row[0], figsize_per_row[1] * len(q2centers) + 1.5), dpi=200, sharex=sharex, squeeze=False)
    handles = []
    labels = []

    repo_path = Path(__file__).resolve().parent.parent
    plot_xlsx = repo_path/'Carbon/one_sheet_to_rule_them_all.xlsx'
    if sheet_CBfit is None:
        sheet_CBfit = pd.read_excel(plot_xlsx, sheet_name='CBfit_qvbin')
    sheet_mc_rl = pd.read_excel(plot_xlsx, sheet_name='mc_rl_q2bin')
    sheet_mc_rt = pd.read_excel(plot_xlsx, sheet_name='mc_rt_q2bin')
    sheet_exp_rl = pd.read_excel(plot_xlsx, sheet_name='exp_rl_q2bin')
    sheet_exp_rt = pd.read_excel(plot_xlsx, sheet_name='exp_rt_q2bin')
    sheet_theory_rl = pd.read_excel(plot_xlsx, sheet_name='theory_rl_q2bin')
    sheet_theory_rt = pd.read_excel(plot_xlsx, sheet_name='theory_rt_q2bin')

    # helper function to plot on given axes
    def plot_on_rl_subplots(ax_rl, q2center):
        # plot vertical dash lines at w = 0.93, 1.07, 1.23 GeV/c^2
        for w in W_PLOT_LIST:
            nu_w = 0.025 + (w**2 + q2center - MASS_NUCLEON**2)/(2 * MASS_NUCLEON)
            ax_rl.axvline(x = nu_w, color = W_PLOT_COLORS[w], linestyle = 'dashdot', lw = 1, label = f'$W={w}$ GeV/c$^2$', zorder=-2)
        nu_wmax = 0.025 + (W_PLOT_LIST[2]**2 + q2center - MASS_NUCLEON**2)/(2 * MASS_NUCLEON)

        # plot Christy Bodek fit
        # 25 Sep 25: shift inelastic to left by 18 MeV; 25 July 18: RT quasi deuteron added
        #  FIXME: 26 May 7 inelastic peak shift is removed for now. address later.
        ChristyBodekFit = sheet_CBfit.loc[sheet_CBfit['q2']==q2center].copy()
        ax_rl.plot(ChristyBodekFit['nu'],ChristyBodekFit["rltot"],color='black',label="$R_L$(total), $R_T$(total) Christy-Bodek-2024", linestyle='solid',lw=0.8, zorder=0)
        ax_rl.plot(ChristyBodekFit['nu'],ChristyBodekFit["rlqe"], color='black',label="$R_L$(QE), $R_T$(QE+TE) Christy-Bodek-2024", linestyle='dotted',lw=0.8, zorder=0)

        # plot mc as lines
        for mc in mc_plot_list:
            df_rl = sheet_mc_rl.loc[(sheet_mc_rl['q2']==q2center) & (sheet_mc_rl['mc']==mc)]
            if len(df_rl) > 0:
                ax_rl.plot(df_rl['nu'], df_rl['rl'], **LINE_STYLES[mc], label=f'$R_L$, $R_T$ {mc}')

        # plot theory as lines
        for theory in theory_plot_list:
            df_rl = sheet_theory_rl.loc[(sheet_theory_rl['q2']==q2center) & (sheet_theory_rl['theory']==theory)]
            if len(df_rl) > 0:
                ax_rl.plot(df_rl['nu'], df_rl['rl'], **LINE_STYLES[theory], label=f'$R_L$, $R_T$ {theory}')
        
        # plot experiment data as scatter and error bars
        for exp in exp_plot_list:
            df_rl = sheet_exp_rl.loc[(sheet_exp_rl['q2']==q2center) & (sheet_exp_rl['experiment']==exp)]
            if len(df_rl) > 0:
                ax_rl.errorbar(df_rl['nu'], df_rl['rl'], yerr=df_rl['rlerr'], **ERRORBAR_STYLES[exp])
                ax_rl.scatter(df_rl['nu'], df_rl['rl'], **SCATTER_STYLES[exp])

        # plot this analysis
        df_q2bin = df_this_analysis.loc[df_this_analysis['q2center'] == q2center]
        if q2center in [0.020, 0.026, 0.040, 0.056, 0.093]:
            # avoid overlapping with Yamaguchi
            nu_yam = sheet_exp_rl.loc[(sheet_exp_rl['q2']==q2center) & (sheet_exp_rl['experiment']=='Yamaguchi')]['nu'].max()
            df_q2bin = df_q2bin.loc[df_q2bin['nu'] >= nu_yam]
        ax_rl.errorbar(df_q2bin['nu'], df_q2bin['rl'], yerr = df_q2bin['rlerr'], ecolor='red', markersize=0, capsize=0, lw=1, fmt='D', elinewidth=1, alpha=0.5, zorder=2) 
        ax_rl.scatter(df_q2bin['nu'], df_q2bin['rl'], color='red', label = '$R_L$, $R_T$ this analysis', marker='D', s=15, edgecolors='gray', linewidth=0.5, zorder=2)
        if len(df_q2bin) > 0 and df_q2bin['nu'].max() > nu_wmax:
            xmax = df_q2bin['nu'].max()*1.1
        else:
            xmax = nu_wmax * 1.1
        axs[i, 0].set_xlim(0, xmax)
        

    # helper function to plot on given axes
    def plot_on_rt_subplots(ax_rt, q2center):
        # plot vertical dash lines at w = 0.93, 1.07, 1.23 GeV/c^2
        for w in W_PLOT_LIST:
            nu_w = 0.025 + (w**2 + q2center - MASS_NUCLEON**2)/(2 * MASS_NUCLEON)
            ax_rt.axvline(x = nu_w, color = W_PLOT_COLORS[w], linestyle = 'dashdot', lw = 1, zorder=-2)
        nu_wmax = 0.025 + (W_PLOT_LIST[2]**2 + q2center - MASS_NUCLEON**2)/(2 * MASS_NUCLEON)

        # plot Christy Bodek fit
        # 25 Sep 25: shift inelastic to left by 18 MeV; 25 July 18: RT quasi deuteron added
        #  FIXME: 26 May 7 inelastic peak shift is removed for now. address later.
        ChristyBodekFit = sheet_CBfit.loc[sheet_CBfit['q2']==q2center].copy()
        ax_rt.plot(ChristyBodekFit['nu'],ChristyBodekFit["rttot"], color='black', linestyle='solid',lw=0.8, zorder=0)
        ax_rt.plot(ChristyBodekFit['nu'],ChristyBodekFit["rtqe"] + ChristyBodekFit["rte"], color='black', linestyle='dotted',lw=0.8, zorder=0)

        # plot mc as lines
        for mc in mc_plot_list:
            df_rt = sheet_mc_rt.loc[(sheet_mc_rt['q2']==q2center) & (sheet_mc_rt['mc']==mc)]
            if len(df_rt) > 0:
                ax_rt.plot(df_rt['nu'], df_rt['rt'], **LINE_STYLES[mc])

        # plot theory as lines
        for theory in theory_plot_list:
            df_rt = sheet_theory_rt.loc[(sheet_theory_rt['q2']==q2center) & (sheet_theory_rt['theory']==theory)]
            if len(df_rt) > 0:
                ax_rt.plot(df_rt['nu'], df_rt['rt'], **LINE_STYLES[theory])
        
        # plot experiment data as scatter and error bars
        for exp in exp_plot_list:
            df_rt = sheet_exp_rt.loc[(sheet_exp_rt['q2']==q2center) & (sheet_exp_rt['experiment']==exp)]
            if len(df_rt) > 0:
                ax_rt.errorbar(df_rt['nu'], df_rt['rt'], yerr=df_rt['rterr'], **ERRORBAR_STYLES[exp])
                ax_rt.scatter(df_rt['nu'], df_rt['rt'], **SCATTER_STYLES[exp], label=f'$R_L$, $R_T$ {exp}')

        # plot this analysis
        df_q2bin = df_this_analysis.loc[df_this_analysis['q2center'] == q2center]
        if q2center in [0.02, 0.026, 0.04, 0.056, 0.093]:
            # avoid overlapping with Yamaguchi
            nu_yam = sheet_exp_rt.loc[(sheet_exp_rt['q2']==q2center) & (sheet_exp_rt['experiment']=='Yamaguchi')]['nu'].max()
            df_q2bin = df_q2bin.loc[df_q2bin['nu'] >= nu_yam]
        ax_rt.errorbar(df_q2bin['nu'], df_q2bin['rt'], yerr = df_q2bin['rterr'], ecolor='red', markersize=0, capsize=0, lw=1, fmt='D', elinewidth=1, alpha=0.5, zorder=2)
        ax_rt.scatter(df_q2bin['nu'], df_q2bin['rt'], color='red', marker='D', s=15, edgecolors='gray', linewidth=0.5, zorder=2)
        if len(df_q2bin) > 0 and df_q2bin['nu'].max() > nu_wmax:
            xmax = df_q2bin['nu'].max()*1.1
        else:
            xmax = nu_wmax * 1.1
        axs[i, 1].set_xlim(0, xmax)


    # loop over given q2 centers
    for i, q2center in enumerate(q2centers):
        if q2center not in Q2CENTERS:
            print(f'q2={q2center} not in spreadsheet')
            continue

        plot_on_rl_subplots(axs[i, 0], q2center)
        plot_on_rt_subplots(axs[i, 1], q2center)

        # handle legend labels
        h, l = axs[i, 0].get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
        h, l = axs[i, 1].get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)

        # main axs styling
        axs[i, 0].text(0.95, 0.95,'$R_L$ ($Q^2$ = '+f'{q2center} GeV$^2$/c$^2$)',transform=axs[i, 0].transAxes,ha='right',va='top',color='gray')
        axs[i, 0].tick_params(which='both', direction='in', top=True, right=True)
        axs[i, 0].minorticks_on()
        axs[i, 0].set_ylabel(r'$R_L$ (GeV$^{-1}$)')
        axs[i, 0].set_ylim(0, RLRT_Q2PLOT_HEIGHTS[q2center][0])
        axs[i, 1].text(0.95, 0.95,'$R_T$ ($Q^2$ = '+f'{q2center} GeV$^2$/c$^2$)',transform=axs[i, 1].transAxes,ha='right',va='top',color='gray')
        axs[i, 1].tick_params(which='both', direction='in', top=True, right=True)
        axs[i, 1].minorticks_on()
        axs[i, 1].set_ylabel(r'$R_T$ (GeV$^{-1}$)')
        axs[i, 1].set_ylim(0, RLRT_Q2PLOT_HEIGHTS[q2center][1])

         # put inset plot at top left corner
        if q2center >= 0.256:
            ax_inset0 = inset_axes(axs[i, 0], width="45%", height="50%", loc='upper left')
            plot_on_rl_subplots(ax_inset0, q2center)
            ax_inset0.yaxis.tick_right()
            ax_inset0.set_xlim(RLRT_Q2INSET_XLIMS[q2center][0], RLRT_Q2INSET_XLIMS[q2center][1])
            ax_inset0.set_ylim(0, RLRT_Q2INSET_HEIGHTS[q2center][0])
            ax_inset0.tick_params(which='both', direction='in', top=True, right=True)
            ax_inset0.minorticks_on()

        if q2center >= 2.75:
            ax_inset1 = inset_axes(axs[i, 1], width="45%", height="50%", loc='upper left')
            plot_on_rt_subplots(ax_inset1, q2center)
            ax_inset1.yaxis.tick_right()
            ax_inset1.set_xlim(RLRT_Q2INSET_XLIMS[q2center][0], RLRT_Q2INSET_XLIMS[q2center][1])
            ax_inset1.set_ylim(0, RLRT_Q2INSET_HEIGHTS[q2center][1])
            ax_inset1.tick_params(which='both', direction='in', top=True, right=True)
            ax_inset1.minorticks_on()

    # add x-axis label only to the last subplot
    axs[-1, 0].set_xlabel(r'$\nu$ (GeV)')
    axs[-1, 1].set_xlabel(r'$\nu$ (GeV)')

    # find unique legend labels, put it at bottom
    unique = dict(zip(labels, handles))
    handles = list(unique.values())
    labels = list(unique.keys())
    if len(q2centers) == 3:
        fig.legend(handles, labels, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.08), frameon=False)
    else:
        fig.legend(handles, labels, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.015), frameon=False)

    # fig.tight_layout()
    if figshow:
        plt.show()
    plt.close()

    return fig

