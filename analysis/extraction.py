import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from .presets import *
from .utilities import rt_quasi_deuteron, linear_model, special_sigmoid
from .christy_bodek_fit import calculate_response_table

# For development: make individual rosenbluth seperation plot
_plot_rosenbluth = False
_plot_rosenbluth_combined = True

def _plot_rosenbluth_xsec_epsilon(df, qvcenter, a_opt, b_opt, chi2_ndf, nuc, w2center):
    plt.figure(figsize=(15,10))
    for ds in df['dataSet'].unique():
        df_ds = df.loc[df['dataSet'] == ds]
        plt.errorbar(df_ds['epsilon'], df_ds['rosenbluth_xsec'], yerr=df_ds['rosenbluth_xsec_err'], fmt='o',alpha=0.5, capsize=3,label=f'{ds}')
        for index, row in df_ds.iterrows():
            row_nu = row['nu']
            plt.text(row['epsilon'], row['rosenbluth_xsec'], f'{ds}:nu={round(row_nu,3)}')

    xs = np.linspace(0,1,10)
    ys = linear_model(xs, a=a_opt, b=b_opt)
    plt.plot(xs, ys, label=f'y = {round(a_opt, 3)} * x + {round(b_opt, 3)}', alpha=0.7, linestyle='--')
    plt.xlabel('epsilon')
    plt.ylabel('rosenbluth xsec')
    plt.title(f'qvc={qvcenter}, nuc={round(nuc,6)}, W2c={round(w2center,3)}, chi2/ndf={round(chi2_ndf,3)}')
    plt.legend()
    plt.show()

def _plot_rosenbluth_xsec_epsilon_combined(df, qvcenter, df_rlrt):
    df_qvbin = df.copy()
    df_qvbin = df_qvbin.sort_values(by='epsilon')

    plt.figure(figsize=(40,30), dpi=150)
    for w2center in np.sort(df_qvbin['W2center_qv'].unique()):
        nuc = np.sqrt(qvcenter**2+w2center)-MASS_NUCLEON
        if nuc in df_rlrt['nu'].unique() and nuc>0.35:
            qvc2 = qvcenter**2
            q2center = qvc2-nuc**2
            df_result = df_rlrt.loc[df_rlrt['nu']==nuc]
            rl = df_result['rl'].iloc[0]
            rt = df_result['rt'].iloc[0]
            a_opt = rl
            b_opt = rt*qvc2 / (2*q2center)
            df_points = df_qvbin.loc[(df_qvbin['W2center_qv']==w2center)]
            plt.errorbar(df_points['epsilon'], df_points['rosenbluth_xsec'], yerr=df_points['rosenbluth_xsec_err'], 
                fmt='o', alpha=0.5)
            plt.errorbar(df_points['epsilon'], df_points['Hcc_Sig(GeV)'], yerr=df_points['Hcc_error(GeV)'], 
                capsize=3, alpha=0.5)
            # plt.scatter(df_points['epsilon'], df_points['rosenbluth_xsec'],alpha=0.5)
            xs = np.linspace(0,1,10)
            ys = linear_model(xs, a=a_opt, b=b_opt)
            plt.plot(xs, ys, label=f'y = {round(a_opt, 3)} * x + {round(b_opt, 3)}', alpha=0.7, linestyle='--')
            plt.text(xs[0]-0.1, ys[0], f'y = {round(a_opt, 3)} * x + {round(b_opt, 3)}, nuc={round(nuc,6)}')
            for index, row in df_points.iterrows():
                row_ds = row['dataSet']
                row_qv = row['qv']
                row_nu = row['nu']
                row_bc = row['bc_qv_w2']
                # plt.text(df_ds['epsilon']+0.01, df_ds['rosenbluth_xsec']-0.01, f'{ds}:nu={round(row_nu,3)}',alpha=0.7)
                plt.text(row['epsilon'], 
                        # row['rosenbluth_xsec'],
                        row['Hcc_Sig(GeV)'],
                        # f'{row_ds}:nu={round(row_nu,6)}'
                        f'{row_ds}:nu={round(row_nu,6)},    q={round(row_qv,4)}, bc={round(row_bc,3)}'
                    )
    plt.xlabel('epsilon')
    plt.xlabel('rosenbluth xsec')

def prepare_dataframe(df_data : pd.DataFrame, vcoul : float = 0.0031, syst_err : bool = True, mass_nucleus : float = MASS_C12) -> pd.DataFrame:
    """
    Prepare a table (pd.DataFrame) for Rosenbluth RL RT seperation.

    Parameters
    ----------
    df_data : pd.DataFrame
        Data table consisting of columns:
        A, Z, nu(GeV), e0(GeV), theta(degrees), xsec, xsec_err, dataset
    vcoul : float, optional
        Value for nucleus coulomb correction. Default: 0.0031 (carbon)
    syst_err : bool, optional
        If True, add systematic error to add to cross section error.
        Default: True
    mass_nucleus : float, optional
        Mass of nucleus in GeV.
        Default: mass of Carbon nucleus, 11.178.

    Returns
    ----------
    pd.DataFrame
        Table that contains necessary information to extract RL RT.    
    """
    ##############################################################################################################
    ##  A       Z       nu          e0          theta           xsec                xsec_err            dataset ##
    ##  int     int     float(GeV)  float(GeV)  float(degrees)  float(nb/sr/GeV)    float(nb/sr/GeV)    int     ##
    ##############################################################################################################
    df = df_data.copy()
    if 'error' not in df.columns:
        df['error']=0.0

    # calculate normalized cross section:
    if 'dataSet' not in df.columns:
        df['dataSet']=-1
        df['normalization']=1
        df['normError']=0
        df['systError']=0
    else:
        df["normalization"]=df["dataSet"].map(NORMALIZATIONS)
        df["normError"]=df["dataSet"].map(NORMALIZATION_ERRORS)
        df["systError"]=df["dataSet"].map(SYSTEMATIC_ERRORS)

    # add systematics error
    df['error_with_syst'] = np.sqrt(
        df['error']**2
        + ((df['systError']*df['cross'])**2)
    )
    # Barreau systematics error stemming from the calorimeter that varies with E':
    df["Ep"]=df["E0"]-df["nu"]
    barreau_sys2 = 0.025/(1.0 + (df.loc[df['dataSet']==1]['Ep'])/0.1)
    df.loc[df['dataSet']==1, 'error_with_syst'] = np.sqrt(
        df.loc[df['dataSet']==1, 'error_with_syst']**2
        +barreau_sys2 * (df.loc[df['dataSet']==1,'cross']**2)
    )
    if syst_err == False:
        df['error_with_syst'] = df['error']

    df['normCross'] = df['cross'] * df['normalization']
    # df['normCrossError']=df['normCross']*np.sqrt((df['error']/df['cross'])**2+(df['normError']/df['normalization'])**2)    
    df['normCrossError']=df['normCross']*np.sqrt((df['error_with_syst']/df['cross'])**2+(df['normError']/df['normalization'])**2)
    df["ThetaRad"]=df["ThetaDeg"]*np.pi/180
    df["sin2(T/2)"]=(np.sin(df["ThetaRad"]/2))**2
    df["cos2(T/2)"]=(np.cos(df["ThetaRad"]/2))**2
    df["tan2(T/2)"]=(np.tan(df["ThetaRad"]/2))**2
    df["Ex"] = df["nu"] - (df["E0"]-df["E0"]/(1+2*df["E0"]*df["sin2(T/2)"] / mass_nucleus))
    df['Veff'] = vcoul
    df["Ffoc2"]=((df["E0"]+df["Veff"])/df["E0"])**2
    df["E0original"]=df["E0"]
    
    # starting below: kinematics turning into their effective values
    df["E0"]=df["E0"]+df["Veff"]
    df["Ep"]=df["Ep"]+df["Veff"]
    df["Q2"]=4*df["E0"]*(df["Ep"])*df["sin2(T/2)"]
    df["qv2"]=df["nu"]**2+df["Q2"]
    df["qv"]=np.sqrt(df["qv2"])
    df["W2"] = MASS_NUCLEON**2 + 2 * MASS_NUCLEON * df["nu"] - df["Q2"]
    df["epsilon"]=1/(1+2*(1+(df["nu"]**2)/df["Q2"])*df["tan2(T/2)"])
    df["gamma"]=ALPHA_FINE*df["Ep"]*(df["W2"]-MASS_NUCLEON**2)/((4*((np.pi)**2)*df["Q2"]*MASS_NUCLEON*df["E0"])*(1-df["epsilon"]))
    df["Sig_R"]=df["normCross"]/df["gamma"]
    df["D_sig_R"]=df["error"]/df["gamma"]
    df["Sig_mott"]=4*(ALPHA_FINE**2)*(df["Ep"]**2)*df["cos2(T/2)"]/(df["Q2"]**2)
    df["Sig_mott"]=df["Ffoc2"]*ALPHA_FINE**2 *df["cos2(T/2)"]*(2*df["E0"]*df["sin2(T/2)"])**-2

    # Calculate the Rosenbluth quantity:
    df["Hcc"]= ((df["qv"]**4)/(4*(ALPHA_FINE**2)*(df["Ep"]**2)*(df["cos2(T/2)"]+2*(df["qv2"]/df["Q2"])*df["sin2(T/2)"])))/df["Ffoc2"]
    df["Hcc_Sig(nb)"]=df["Hcc"]*df["normCross"]
    df["Hcc_error(nb)"]=df["Hcc"]*df["normCrossError"]
    df["Hcc_Sig(GeV)"]=df["Hcc_Sig(nb)"]/((0.1973269**2)*10000000)
    df["Hcc_error(GeV)"]=df["Hcc_error(nb)"]/((0.1973269**2)*10000000)

    # # RT quasi deuteron added 2025 July 18
    # df["RT_QD_data"] = rt_quasi_deuteron(nus=df['nu'],q2s = df['Q2'],exs = df['Ex'])

    # split into qv, q2 bins:
    df['qvbin'] = 0.0
    df['qvcenter'] = 0.0
    df["qvbin"]=pd.cut(x=df["qv"],bins=QVBINS,labels=QVBIN_NAMES,right=True)
    df["qvcenter"]=pd.cut(x=df["qv"],bins=QVBINS,labels=QVCENTERS,right=True)
    df['qvcenter']=pd.to_numeric(df['qvcenter'])
    df['Q2bin'] = 0.0
    df['Q2center'] = 0.0
    df["Q2bin"]=pd.cut(x=df["Q2"],bins=Q2BINS,labels=Q2BIN_NAMES,right=True)
    df["Q2center"]=pd.cut(x=df["Q2"],bins=Q2BINS,labels=Q2CENTERS,right=True)
    df['Q2center']=pd.to_numeric(df['Q2center'])
    # df = df.dropna()
    df['W2bin_q2']=0.0
    df['W2center_q2']=0.0
    df['nucenter_w2_q2']=0.0
    df['epcenter_w2_q2']=0.0
    df['Exbin_q2']=0.0
    df['Excenter_q2']=0.0
    df['nucenter_ex_q2']=0.0
    df['epcenter_ex_q2']=0.0
    df['W2bin_qv']=0.0
    df['W2center_qv']=0.0
    df['nucenter_w2_qv']=0.0
    df['epcenter_w2_qv']=0.0
    df['Exbin_qv']=0.0
    df['Excenter_qv']=0.0
    df['nucenter_ex_qv']=0.0
    df['epcenter_ex_qv']=0.0

    for Q2center in Q2CENTERS:
        ## ex < EX_CUT: extract at ex center
        mask = (df['Q2center'] == Q2center) & (df['Ex'] < EX_CUT)
        if Q2center == 0.01:
            mask = (df['Q2center'] == Q2center) & (df['Ex'] < EX_CUT_LOWQ)

        if len(df.loc[mask, 'Ex']) > 0:
            Exedges = np.array(Q2CENTER_EXBINEDGES[Q2center])
            Excenters = (Exedges[:-1] + Exedges[1:]) / 2
            df.loc[mask, 'Exbin_q2'] = pd.cut(df.loc[mask, 'Ex'], bins=Exedges, labels=False, include_lowest=True,duplicates='drop')
            df.loc[mask, 'Excenter_q2'] = df.loc[mask, 'Exbin_q2'].map(lambda i: Excenters[int(i)] if pd.notnull(i) else np.nan)
            df.loc[mask, 'nucenter_ex_q2']=df.loc[mask, 'Excenter_q2']+Q2center/(2*mass_nucleus)
            df.loc[mask, 'epcenter_ex_q2']=1/(1+2*(1+(df.loc[mask, 'nucenter_ex_q2']**2)/Q2center)*df.loc[mask, 'tan2(T/2)']) 

        ## Ex >= EX_CUT: extract at W2 center
        mask = (df['Q2center'] == Q2center) & (df['Ex'] >= EX_CUT)
        if Q2center == 0.01:
            mask = (df['Q2center'] == Q2center) & (df['Ex'] > EX_CUT_LOWQ)

        if len(df.loc[mask, 'W2']) > 0:
            nuedges = np.array(Q2CENTER_NUBINEDGES[Q2center])
            W2edges = MASS_NUCLEON**2 + 2*MASS_NUCLEON*nuedges - Q2center
            W2centers = (W2edges[:-1] + W2edges[1:]) / 2
            df.loc[mask, 'W2bin_q2'] = pd.cut(df.loc[mask, 'W2'], bins=W2edges, labels=False, include_lowest=True,duplicates='drop')
            df.loc[mask, 'W2center_q2'] = df.loc[mask, 'W2bin_q2'].map(lambda i: W2centers[int(i)] if pd.notnull(i) else np.nan)
            df.loc[mask, 'nucenter_w2_q2']=(df.loc[mask, 'W2center_q2']-MASS_NUCLEON**2 + Q2center)/(2*MASS_NUCLEON)
            df.loc[mask, 'epcenter_w2_q2']=1/(1+2*(1+(df.loc[mask, 'nucenter_w2_q2']**2)/Q2center)*df.loc[mask, 'tan2(T/2)']) 

    for qvcenter in QVCENTERS:
        ## Ex < EX_CUT: extract at Ex center
        mask = (df['qvcenter'] == qvcenter) & (df['Ex'] < EX_CUT)
        if qvcenter == 0.1:
            mask = (df['qvcenter'] == qvcenter) & (df['Ex'] < EX_CUT_LOWQ)

        if len(df.loc[mask, 'Ex'])>0:
            Exedges = np.array(QVCENTER_EXBINEDGES[qvcenter])
            Excenters = (Exedges[:-1] + Exedges[1:]) / 2
            df.loc[mask, 'Exbin_qv'] = pd.cut(df.loc[mask, 'Ex'], bins=Exedges, labels=False, include_lowest=True,duplicates='drop')
            df.loc[mask, 'Excenter_qv'] = df.loc[mask, 'Exbin_qv'].map(lambda i: Excenters[int(i)] if pd.notnull(i) else np.nan)
            df.loc[mask, 'nucenter_ex_qv'] = np.sqrt(mass_nucleus**2+qvcenter**2+2*mass_nucleus*df.loc[mask, 'Excenter_qv'])-mass_nucleus
            df.loc[mask, 'epcenter_ex_qv']=1/(1+2*(1+(df.loc[mask, 'nucenter_ex_qv']**2)/(qvcenter**2-df.loc[mask, 'nucenter_ex_qv']**2))*df.loc[mask, 'tan2(T/2)']) 

        ## Ex >= EX_CUT: extract at W2 center or (nu center)
        mask = (df['qvcenter'] == qvcenter) & (df['Ex'] >= EX_CUT)        
        # for qvcenter = 0.1, only extract at Ex center, so skip 
        if qvcenter > 0.1:
            nuedges = np.array(QVCENTER_NUBINEDGES[qvcenter])
            if qvcenter in [0.3, 0.38, 0.57]: # use nu bin center, so to overlap with Barreau and Jordan's point
                nucenters = (nuedges[:-1] + nuedges[1:]) / 2
                df.loc[mask, 'W2bin_qv'] = pd.cut(df.loc[mask, 'nu'], bins=nuedges, labels=False, include_lowest=True)        
                df.loc[mask, 'nucenter_w2_qv'] = df.loc[mask, 'W2bin_qv'].map(lambda i: nucenters[int(i)] if pd.notnull(i) else np.nan)
                df.loc[mask, 'epcenter_w2_qv']=1/(1+2*(1+(df.loc[mask, 'nucenter_w2_qv']**2)/(qvcenter**2-df.loc[mask, 'nucenter_w2_qv']**2))*df.loc[mask, 'tan2(T/2)'])
                # df.loc[mask, 'epcenter_w2_qv']=df.loc[mask, 'epsilon']
                
                df.loc[mask, 'W2center_qv'] = MASS_NUCLEON**2 + 2*MASS_NUCLEON*df.loc[mask, 'nucenter_w2_qv'] + df.loc[mask, 'nucenter_w2_qv']**2 - qvcenter**2
            else: # use W2 bin center
                W2edges = MASS_NUCLEON**2 + 2*MASS_NUCLEON*nuedges + nuedges**2 - qvcenter**2
                W2centers = (W2edges[:-1] + W2edges[1:]) / 2
                df.loc[mask, 'W2bin_qv'] = pd.cut(df.loc[mask, 'W2'], bins=W2edges, labels=False, include_lowest=True)
                df.loc[mask, 'W2center_qv'] = df.loc[mask, 'W2bin_qv'].map(lambda i: W2centers[int(i)] if pd.notnull(i) else np.nan)
                df.loc[mask, 'nucenter_w2_qv'] = np.sqrt(qvcenter**2+df.loc[mask, 'W2center_qv'])-MASS_NUCLEON
                df.loc[mask, 'epcenter_w2_qv']=1/(1+2*(1+(df.loc[mask, 'nucenter_w2_qv']**2)/(qvcenter**2-df.loc[mask, 'nucenter_w2_qv']**2))*df.loc[mask, 'tan2(T/2)'])

    # df=df.dropna()
    return df

def calculate_response_table_update_qd_ie(df_qv_nu : pd.DataFrame, a : float = 12.0, z : float = 6.0):
    # FIXME: this function is outdated. Use rtqd calculated in Fortran. - Ziggy Aug 10 2026
    df = calculate_response_table(table = df_qv_nu, a=a, z=z)
    # FIXME: the shift is wrong. Don't do the entire dataframe.
    # # shift the inelastic peak at low q2
    # q2_suppression = special_sigmoid(df['q2'], center = 0.03, width= 0.005) # apply shifts at low q2 only
    # f_rtie = interp1d(df['nu'], df['rtie'], kind="linear", bounds_error=False, fill_value=0.0)
    # rtie_shifts = 1.06*f_rtie(df['nu'] - 0.018) - df['rtie']
    # df['rtie'] = df['rtie'] + rtie_shifts * q2_suppression
    # f_rlie = interp1d(df['nu'], df['rlie'], kind="linear", bounds_error=False, fill_value=0.0)
    # rlie_shifts = 1.06*f_rlie(df['nu'] - 0.018) - df['rlie']
    # df['rlie'] = df['rlie'] + rlie_shifts * q2_suppression

    # quasi-deuteron contribution
    df['rtqd'] = rt_quasi_deuteron(nus = df['nu'], q2s = df['q2'], exs = df['ex'])
    df['rttot'] = df['rtqe'] + df['rte'] + df['rtie'] + df['rtns'] + df['rtqd']
    df['rltot'] = df['rlqe'] + df['rle'] + df['rlie'] + df['rlns']
    for col in ["rttot", "rltot", "rtqe", "rlqe", "rtie", "rlie", "rte", "rle", "rtns", "rlns"]:
        df[col] = df[col] * 1e3 # convert from MeV^-1 to GeV^-1
    return df

def calculate_bin_centering_correction(df_xsec : pd.DataFrame, mass_nucleus : float = MASS_C12, a: float = 12.0, z: float = 6.0) -> pd.DataFrame:
    """
    Calculate bin-centering correction factors for Rosenbluth reduced cross
    section.

    Parameters
    ----------
    df_xsec : pd.DataFrame
        Electron scattering cross section data table consisting of columns:
        A, Z, nu(GeV), e0(GeV), theta(degrees), xsec, xsec_err, dataset
    vcoul : float, optional
        Value for nucleus coulomb correction. Default: 0.0031 (carbon)
    syst_err : float, optional
        Systematic error (0 ~ 1) to add to cross section error.
    mass_nucleus : float, optional
        Mass of nucleus in GeV.
        Default: mass of Carbon nucleus, 11.178.
    a : float
        Nuclear mass number A.
    z : float
        Nuclear charge Z.

    Returns
    ----------
    pd.DataFrame
        Table that contains necessary information to extract RL RT.    
    """
    df = df_xsec.copy()

    # ====== calculate bin centering correction factor to Q2 W2 bin center, bc_q2_w2 ======
    # CBfit response values at Q2 W2 bin center:
    nus = (df['W2center_q2'] + df['Q2center'] - MASS_NUCLEON**2) / (2*MASS_NUCLEON)
    qvs = np.sqrt(df['Q2center'] + nus**2)
    df_response = pd.DataFrame({'qv':qvs, 'nu':nus})
    df_response = calculate_response_table(df_response, a = a, z = z)
    df['RL_q2c_w2'] = df_response['rltot'].values
    df['RT_q2c_w2'] = df_response['rttot'].values
    # df['RT_q2c_w2'] = df['RT_q2c_w2'] + rt_quasi_deuteron(nus=df['nucenter_w2_q2'],q2s = df['Q2center'],
    #     exs = df['nucenter_w2_q2']-df['Q2center']/(2*mass_nucleus)) # RT quasi deuteron added 2025 July 18
    # CBfit response values at data effective Q2 W2:
    nus = (df['W2'] + df['Q2'] - MASS_NUCLEON**2) / (2*MASS_NUCLEON)
    qvs = np.sqrt(df['Q2'] + nus**2)
    df_response = pd.DataFrame({'qv':qvs, 'nu':nus})
    df_response = calculate_response_table(df_response, a = a, z = z)
    df['RL_q2d_w2'] = df_response['rltot'].values
    df['RT_q2d_w2'] = df_response['rttot'].values
    # df['RT_q2d_w2'] = df['RT_q2d_w2'] + df['RT_QD_data'] # RT quasi deuteron added 2025 July 18
    df['bc_q2_w2']=1.0
    for Q2center in Q2CENTERS:
        mask = (df['Q2center'] == Q2center) & (df['Ex'] >= EX_CUT) # use Ex >= EX_CUT
        df.loc[mask, 'bc_q2_w2']=(df.loc[mask, 'epcenter_w2_q2']*df.loc[mask, 'RL_q2c_w2']
                +0.5*((Q2center+df.loc[mask, 'nucenter_w2_q2']**2)/Q2center)*df.loc[mask, 'RT_q2c_w2']
            )/(df.loc[mask, 'epsilon']*df.loc[mask, 'RL_q2d_w2']
            +0.5*(df.loc[mask, 'qv2']/df.loc[mask, 'Q2'])*df.loc[mask, 'RT_q2d_w2'])
    print('RL RT bin centering correction factor bc_q2_w2 done.')


    # ====== calculate bin centering correction factor to Q2 Ex bin center, bc_q2_ex ======
    # CBfit response values at Q2 Ex bin center:
    nus = df['Excenter_q2'] + df['Q2center']/(2*mass_nucleus)
    qvs = np.sqrt(df['Q2center'] + nus**2)
    df_response = pd.DataFrame({'qv':qvs, 'nu':nus})
    df_response = calculate_response_table(df_response, a = a, z = z)
    df['RL_q2c_ex'] = df_response['rltot'].values
    df['RT_q2c_ex'] = df_response['rttot'].values
    # df['RT_q2c_ex'] = df['RT_q2c_ex'] + rt_quasi_deuteron(nus=df['nucenter_ex_q2'],q2s = df['Q2center'],exs = df['Excenter_q2'])# RT quasi deuteron added 2025 July 18
    # # CBfit response values at data effective Q2 Ex:
    nus = df['Ex'] + df['Q2']/(2*mass_nucleus)
    qvs = np.sqrt(df['Q2'] + nus**2)
    df_response = calculate_response_table(df_response, a = a, z = z)
    df['RL_q2d_ex'] = df_response['rltot'].values
    df['RT_q2d_ex'] = df_response['rttot'].values
    # df['RT_q2d_ex'] = df['RT_q2d_ex'] + + df['RT_QD_data'] # RT quasi deuteron added 2025 July 18
    df['bc_q2_ex']=1.0
    for Q2center in Q2CENTERS:
        mask = (df['Q2center'] == Q2center) & (df['Ex'] < EX_CUT) # use Ex < EX_CUT
        if Q2center == 0.01:
            mask = (df['Q2center'] == Q2center) & (df['Ex'] < EX_CUT_LOWQ)
        df.loc[mask, 'bc_q2_ex']=(df.loc[mask, 'epcenter_ex_q2']*df.loc[mask, 'RL_q2c_ex']
                +0.5*((Q2center+df.loc[mask, 'nucenter_ex_q2']**2)/Q2center)*df.loc[mask, 'RT_q2c_ex']
            )/(df.loc[mask, 'epsilon']*df.loc[mask, 'RL_q2d_ex']
            +0.5*(df.loc[mask, 'qv2']/df.loc[mask, 'Q2'])*df.loc[mask, 'RT_q2d_ex'])
    print('RL RT bin centering correction factor bc_q2_ex done.')

    # ====== calculate bin centering correction factor to qv W2 bin center, bc_qv_w2 ======
    # CBfit response values at qv W2 bin center:
    nus = np.sqrt(df['qvcenter']**2 + df['W2center_qv']) - MASS_NUCLEON
    qvs = df['qvcenter']
    df_response = pd.DataFrame({'qv':qvs, 'nu':nus})
    df_response = calculate_response_table(df_response, a = a, z = z)
    df['RL_qvc_w2'] = df_response['rltot'].values
    df['RT_qvc_w2'] = df_response['rttot'].values
    # df['RT_qvc_w2'] = df['RT_qvc_w2'] + rt_quasi_deuteron(nus=df['nucenter_w2_qv'],q2s = df['qvcenter']**2-df['nucenter_w2_qv']**2,
    #                             exs = df['nucenter_w2_qv'] - (df['qvcenter']**2-df['nucenter_w2_qv']**2)/(2*mass_nucleus)) # RT quasi deuteron added 2025 July 18
    # CBfit response values at data effective qv W2:
    nus = np.sqrt(df['qv']**2 + df['W2']) - MASS_NUCLEON
    qvs = df['qv']
    df_response = pd.DataFrame({'qv':qvs, 'nu':nus})
    df_response = calculate_response_table(df_response, a = a, z = z)
    df['RL_qvd_w2'] = df_response['rltot'].values
    df['RT_qvd_w2'] = df_response['rttot'].values
    # df['RT_qvd_w2'] = df['RT_qvd_w2'] + df['RT_QD_data'] # RT quasi deuteron added 2025 July 18
    df['bc_qv_w2']=1.0
    for qvcenter in QVCENTERS:
        ## Ex >= 50MeV:
        mask = (df['qvcenter'] == qvcenter) & (df['Ex'] >= EX_CUT) # use Ex >= EX_CUT
        df.loc[mask, 'bc_qv_w2']=(df.loc[mask, 'epcenter_w2_qv']*df.loc[mask, 'RL_qvc_w2']
                +0.5*((qvcenter**2)/(qvcenter**2-df.loc[mask,'nucenter_w2_qv']**2))*df.loc[mask, 'RT_qvc_w2']
            )/(df.loc[mask, 'epsilon']*df.loc[mask, 'RL_qvd_w2']
            +0.5*(df.loc[mask, 'qv2']/df.loc[mask, 'Q2'])*df.loc[mask, 'RT_qvd_w2'])
    print('RL RT bin centering correction factor bc_qv_w2 done.')

    # ====== calculate bin centering correction factor to qv Ex bin center, bc_qv_ex ======
    # CBfit response values at qv Ex bin center:
    nus = np.sqrt(mass_nucleus**2 + df['qvcenter']**2 + 2*mass_nucleus*df['Excenter_qv']) - mass_nucleus
    qvs = df['qvcenter']
    df_response = pd.DataFrame({'qv':qvs, 'nu':nus})
    df_response = calculate_response_table(df_response, a = a, z = z)
    df['RL_qvc_ex'] = df_response['rltot'].values
    df['RT_qvc_ex'] = df_response['rttot'].values
    # df['RT_qvc_ex'] = df['RT_qvc_ex'] + rt_quasi_deuteron(nus=df['nucenter_ex_qv'],q2s = df['qvcenter']**2-df['nucenter_ex_qv']**2,
    #     exs = df['Excenter_qv']) # RT quasi deuteron added 2025 July 18
    # CBfit response values at data effective qv Ex:
    nus = np.sqrt(mass_nucleus**2 + df['qv']**2 + 2*mass_nucleus*df['Ex']) - mass_nucleus
    qvs = df['qv']
    df_response = pd.DataFrame({'qv':qvs, 'nu':nus})
    df_response = calculate_response_table(df_response, a = a, z = z)
    df['RL_qvd_ex'] = df_response['rltot'].values
    df['RT_qvd_ex'] = df_response['rttot'].values
    # df['RT_qvd_ex'] = df['RT_qvd_ex'] + df['RT_QD_data'] # RT quasi deuteron added 2025 July 18
    df['bc_qv_ex']=1.0
    for qvcenter in QVCENTERS:
        mask = (df['qvcenter'] == qvcenter) & (df['Ex'] < EX_CUT) # use Ex < EX_CUT:
        if qvcenter == 0.1:
            mask = (df['qvcenter'] == qvcenter) & (df['Ex'] < EX_CUT_LOWQ)
        df.loc[mask, 'bc_qv_ex']=(df.loc[mask, 'epcenter_ex_qv']*df.loc[mask, 'RL_qvc_ex']
                +0.5*((qvcenter**2)/(qvcenter**2-df.loc[mask,'nucenter_ex_qv']**2))*df.loc[mask, 'RT_qvc_ex']
            )/(df.loc[mask, 'epsilon']*df.loc[mask, 'RL_qvd_ex']
            +0.5*(df.loc[mask, 'qv2']/df.loc[mask, 'Q2'])*df.loc[mask, 'RT_qvd_ex'])
    print('RL RT bin centering correction factor bc_qv_ex done.')

    return df

def calculate_bc_qv_w2(df_xsec : pd.DataFrame, mass_nucleus : float = MASS_C12, a: float = 12.0, z: float = 6.0) -> pd.DataFrame:
    """
    Calculate bin-centering correction factors for Rosenbluth reduced cross
    section.

    Parameters
    ----------
    df_xsec : pd.DataFrame
        Electron scattering cross section data table consisting of columns:
        A, Z, nu(GeV), e0(GeV), theta(degrees), xsec, xsec_err, dataset
    vcoul : float, optional
        Value for nucleus coulomb correction. Default: 0.0031 (carbon)
    syst_err : float, optional
        Systematic error (0 ~ 1) to add to cross section error.
    mass_nucleus : float, optional
        Mass of nucleus in GeV.
        Default: mass of Carbon nucleus, 11.178.
    a : float
        Nuclear mass number A.
    z : float
        Nuclear charge Z.

    Returns
    ----------
    pd.DataFrame
        Table that contains necessary information to extract RL RT.    
    """
    df = df_xsec.copy()

    # ====== calculate bin centering correction factor to qv W2 bin center, bc_qv_w2 ======
    # CBfit response values at qv W2 bin center:
    nus = np.sqrt(df['qvcenter']**2 + df['W2center_qv']) - MASS_NUCLEON
    qvs = df['qvcenter']
    df_response = pd.DataFrame({'qv':qvs, 'nu':nus})
    df_response = calculate_response_table(df_response, a = a, z = z)
    df['RL_qvc_w2'] = df_response['rltot'].values
    df['RT_qvc_w2'] = df_response['rttot'].values
    # df['RT_qvc_w2'] = df['RT_qvc_w2'] + rt_quasi_deuteron(nus=df['nucenter_w2_qv'],q2s = df['qvcenter']**2-df['nucenter_w2_qv']**2,
    #                             exs = df['nucenter_w2_qv'] - (df['qvcenter']**2-df['nucenter_w2_qv']**2)/(2*mass_nucleus)) # RT quasi deuteron added 2025 July 18
    # CBfit response values at data effective qv W2:
    nus = np.sqrt(df['qv']**2 + df['W2']) - MASS_NUCLEON
    qvs = df['qv']
    df_response = pd.DataFrame({'qv':qvs, 'nu':nus})
    df_response = calculate_response_table(df_response, a = a, z = z)
    df['RL_qvd_w2'] = df_response['rltot'].values
    df['RT_qvd_w2'] = df_response['rttot'].values
    # df['RT_qvd_w2'] = df['RT_qvd_w2'] + df['RT_QD_data'] # RT quasi deuteron added 2025 July 18
    df['bc_qv_w2']=1.0
    for qvcenter in QVCENTERS:
        ## Ex >= 50MeV:
        mask = (df['qvcenter'] == qvcenter) & (df['Ex'] >= EX_CUT) # use Ex >= EX_CUT
        df.loc[mask, 'bc_qv_w2']=(df.loc[mask, 'epcenter_w2_qv']*df.loc[mask, 'RL_qvc_w2']
                +0.5*((qvcenter**2)/(qvcenter**2-df.loc[mask,'nucenter_w2_qv']**2))*df.loc[mask, 'RT_qvc_w2']
            )/(df.loc[mask, 'epsilon']*df.loc[mask, 'RL_qvd_w2']
            +0.5*(df.loc[mask, 'qv2']/df.loc[mask, 'Q2'])*df.loc[mask, 'RT_qvd_w2'])
    print('RL RT bin centering correction factor bc_qv_w2 done.')


    return df

def extract_response_q2bin_w2center(df : pd.DataFrame, q2center : float = 0.01, bin_centering : bool = True,
        absolute_sigma : bool = True, min_epsilon_range : int = 0.3, mass_nucleus : float = MASS_C12) -> pd.DataFrame:
    if q2center == 0.01:
        df_q2bin = df.loc[(df['Q2center']==q2center) & (df['Ex'] >= EX_CUT_LOWQ)].copy()
    else:
        df_q2bin = df.loc[(df['Q2center']==q2center) & (df['Ex'] >= EX_CUT)].copy()
    if bin_centering == True:
        df_q2bin['rosenbluth_xsec']=df_q2bin['Hcc_Sig(GeV)']*df_q2bin['bc_q2_w2']
        df_q2bin['rosenbluth_xsec_err']=df_q2bin['Hcc_error(GeV)']*df_q2bin['bc_q2_w2'] 
    else:
        df_q2bin['rosenbluth_xsec']=df_q2bin['Hcc_Sig(GeV)']
        df_q2bin['rosenbluth_xsec_err']=df_q2bin['Hcc_error(GeV)']
    
    df_rlrt = []
    for w2center in np.sort(df_q2bin['W2center_q2'].unique()):
        df_w2bin = df_q2bin.loc[(df_q2bin['W2center_q2']==w2center)]
        nuc = (w2center - MASS_NUCLEON**2 + q2center)/(2 * MASS_NUCLEON)
        qvc2 = (q2center + nuc**2)
        exc = nuc - q2center/(2 * mass_nucleus)
        x = np.array(df_w2bin["epsilon"])        
        y = np.array(df_w2bin["rosenbluth_xsec"])
        yerr = np.array(df_w2bin["rosenbluth_xsec_err"])

        if (len(y)>2) and ((np.max(x) - np.min(x)) >= min_epsilon_range):
            params, covariance = curve_fit(linear_model, x, y, sigma=yerr, absolute_sigma = absolute_sigma)
            a_opt, b_opt = params
            a_err, b_err = np.sqrt(np.diag(covariance))
            chi2 = np.sum(np.square((y - linear_model(x, a_opt, b_opt)) / yerr))
            # rl, rlerr, rt, rterr in GeV^-1:
            rl = a_opt
            rlerr = a_err
            rt = (2 * b_opt * q2center / qvc2)
            rterr = (2 * b_err * q2center / qvc2)
            df_rlrt.append({'q2center':q2center, 'nu':nuc, 'w2':w2center, 'ex':exc, 'rl':rl, 'rlerr':rlerr, 'rt':rt,
                'rterr':rterr, 'chi2':chi2, 'num_points':len(y)})

    df_rlrt = pd.DataFrame(df_rlrt)
    # if q2center in [0.02,0.026]:
    #     df_rlrt = df_rlrt.loc[df_rlrt['nu']<0.35]
    return df_rlrt

def extract_response_q2bin_excenter(df : pd.DataFrame, q2center : float = 0.01, bin_centering : bool = True, 
        absolute_sigma : bool = True, min_epsilon_range : int = 0.03, mass_nucleus : float = MASS_C12):
    if q2center <= 0.01:
        df_q2bin = df.loc[(df['Q2center'] == q2center) & (df['Ex'] < EX_CUT_LOWQ)].copy()
    else:
        df_q2bin = df.loc[(df['Q2center'] == q2center) & (df['Ex'] < EX_CUT)].copy()
    if bin_centering == True:
        df_q2bin['rosenbluth_xsec']=df_q2bin['Hcc_Sig(GeV)'] * df_q2bin['bc_q2_ex']
        df_q2bin['rosenbluth_xsec_err']=df_q2bin['Hcc_error(GeV)']*df_q2bin['bc_q2_ex']
    else:
        df_q2bin['rosenbluth_xsec']=df_q2bin['Hcc_Sig(GeV)']
        df_q2bin['rosenbluth_xsec_err']=df_q2bin['Hcc_error(GeV)']

    df_rlrt = []
    for excenter in np.sort(df_q2bin['Excenter_q2'].unique()):
        df_exbin = df_q2bin.loc[(df_q2bin['Excenter_q2'] == excenter)]
        nuc = excenter + q2center/(2*mass_nucleus)
        qvc2 = (q2center+nuc**2)
        w2c = MASS_NUCLEON**2 + 2 * MASS_NUCLEON * nuc - q2center
        x = np.array(df_exbin["epsilon"])        
        y = np.array(df_exbin["rosenbluth_xsec"])
        yerr = np.array(df_exbin["rosenbluth_xsec_err"])

        if (nuc > 0) and (len(y) > 2) and ((np.max(x) - np.min(x)) >= min_epsilon_range):
            params, covariance = curve_fit(linear_model, x, y, sigma=yerr, absolute_sigma=absolute_sigma)
            a_opt, b_opt = params
            a_err, b_err = np.sqrt(np.diag(covariance))
            chi2 = np.sum(np.square((y-linear_model(x,a_opt,b_opt))/yerr))
            # rl, rlerr, rt, rterr in GeV^-1:
            rl = a_opt
            rlerr = a_err
            rt = (2*b_opt*q2center/qvc2)
            rterr = (2*b_err*q2center/qvc2)
            df_rlrt.append({'q2center':q2center,'nu':nuc,'w2':w2c,'ex':excenter,'rl':rl,'rlerr':rlerr,'rt':rt,'rterr':rterr,'chi2':chi2,
                    'num_points':len(y)})
    
    df_rlrt = pd.DataFrame(df_rlrt)
    return df_rlrt
    
def extract_response_qvbin_w2center(df : pd.DataFrame, qvcenter : float = 0.01, bin_centering : bool = True,
        absolute_sigma : bool = True, min_epsilon_range : int = 0.3, mass_nucleus : float = MASS_C12) -> pd.DataFrame:
    if qvcenter <= 0.1:
        df_qvbin = df.loc[(df['qvcenter']==qvcenter) & (df['Ex']>=EX_CUT_LOWQ)].copy()
    else:
        df_qvbin = df.loc[(df['qvcenter']==qvcenter) & (df['Ex']>=EX_CUT)].copy()
    if bin_centering == True:
        df_qvbin['rosenbluth_xsec']=df_qvbin['Hcc_Sig(GeV)'] * df_qvbin['bc_qv_w2']
        df_qvbin['rosenbluth_xsec_err']=df_qvbin['Hcc_error(GeV)'] * df_qvbin['bc_qv_w2']
    else:
        df_qvbin['rosenbluth_xsec']=df_qvbin['Hcc_Sig(GeV)']
        df_qvbin['rosenbluth_xsec_err']=df_qvbin['Hcc_error(GeV)']

    df_rlrt = []
    for w2center in np.sort(df_qvbin['W2center_qv'].unique()):
        df_w2bin = df_qvbin.loc[(df_qvbin['W2center_qv']==w2center)]
        nuc = np.sqrt(qvcenter**2+w2center)-MASS_NUCLEON
        qvc2 = qvcenter**2
        q2center = qvc2-nuc**2
        exc = nuc - q2center/(2*mass_nucleus)
        x = np.array(df_w2bin["epsilon"])
        y = np.array(df_w2bin["rosenbluth_xsec"])
        yerr = np.array(df_w2bin["rosenbluth_xsec_err"])

        if (len(y)>2) and ((np.max(x)-np.min(x)) >= min_epsilon_range):
            params, covariance = curve_fit(linear_model, x, y, sigma = yerr, absolute_sigma = absolute_sigma)
            a_opt, b_opt = params
            a_err, b_err = np.sqrt(np.diag(covariance))            
            chi2 = np.sum(np.square((y - linear_model(x,a_opt,b_opt)) / yerr))
            # rl, rlerr, rt, rterr in GeV^-1:
            rl = a_opt
            rlerr = a_err
            rt = (2 * b_opt * q2center / qvc2)
            rterr = (2 * b_err * q2center / qvc2)
            df_rlrt.append({'qvcenter':qvcenter,'nu':nuc,'w2':w2center,'ex':exc,'rl':rl,'rlerr':rlerr,'rt':rt,'rterr':rterr,'chi2':chi2,
                    'num_points':len(y)})             
            if _plot_rosenbluth:
                _plot_rosenbluth_xsec_epsilon(df_w2bin, qvcenter, a_opt, b_opt, chi2/(len(y)-2), nuc, w2center)


    df_rlrt = pd.DataFrame(df_rlrt)
    if _plot_rosenbluth_combined:
        _plot_rosenbluth_xsec_epsilon_combined(df_qvbin, qvcenter, df_rlrt)


    return df_rlrt

def extract_response_qvbin_excenter(df : pd.DataFrame, qvcenter : float = 0.01, bin_centering : bool = True,
        absolute_sigma : bool = True, min_epsilon_range : int = 0.3, mass_nucleus : float = MASS_C12) -> pd.DataFrame:
    if qvcenter <= 0.1:
        df_qvbin = df.loc[(df['qvcenter']==qvcenter) & (df['Ex']<EX_CUT_LOWQ)].copy()
    else:
        df_qvbin = df.loc[(df['qvcenter']==qvcenter) & (df['Ex']<EX_CUT)].copy()
    if bin_centering == True:
        df_qvbin['rosenbluth_xsec']=df_qvbin['Hcc_Sig(GeV)'] * df_qvbin['bc_qv_ex']
        df_qvbin['rosenbluth_xsec_err']=df_qvbin['Hcc_error(GeV)'] * df_qvbin['bc_qv_ex']
    else:
        df_qvbin['rosenbluth_xsec']=df_qvbin['Hcc_Sig(GeV)']
        df_qvbin['rosenbluth_xsec_err']=df_qvbin['Hcc_error(GeV)']

    df_rlrt = []
    for excenter in df_qvbin['Excenter_qv'].unique():
        df_exbin = df_qvbin.loc[(df_qvbin['Excenter_qv']==excenter)]
        nuc = np.sqrt(mass_nucleus**2+qvcenter**2+2*mass_nucleus*excenter)-mass_nucleus
        qvc2 = qvcenter**2
        q2center = qvc2 - nuc**2
        w2c = MASS_NUCLEON**2 + 2 * MASS_NUCLEON * nuc - q2center
        x = np.array(df_exbin["epsilon"])
        y = np.array(df_exbin["rosenbluth_xsec"])
        yerr = np.array(df_exbin["rosenbluth_xsec_err"])

        if (len(y)>2) and (np.max(x)-np.min(x))>=min_epsilon_range:
            params, covariance = curve_fit(linear_model, x, y, sigma=yerr, absolute_sigma=absolute_sigma)
            a_opt, b_opt = params
            a_err, b_err = np.sqrt(np.diag(covariance))
            chi2 = np.sum(np.square((y-linear_model(x,a_opt,b_opt))/yerr))
            # rl, rlerr, rt, rterr in GeV^-1:
            rl = a_opt
            rlerr = a_err
            rt = (2*b_opt*q2center/qvc2)
            rterr = (2*b_err*q2center/qvc2)
            df_rlrt.append({'qvcenter':qvcenter,'nu':nuc,'w2':w2c,'ex':excenter,'rl':rl,'rlerr':rlerr,'rt':rt,'rterr':rterr,'chi2':chi2,
                    'num_points':len(y)})

    df_rlrt = pd.DataFrame(df_rlrt)
    return df_rlrt

def extract_response_qvbins(df : pd.DataFrame, qvcenters : list[float] = QVCENTERS, **kwargs) -> pd.DataFrame:
    df_rlrt = []
    for qvcenter in qvcenters:
        df_rlrt.append(extract_response_qvbin_excenter(df=df, qvcenter=qvcenter, **kwargs))
        df_rlrt.append(extract_response_qvbin_w2center(df=df, qvcenter=qvcenter, **kwargs))
    df_rlrt = pd.concat(df_rlrt)
    return df_rlrt

def extract_response_q2bins(df : pd.DataFrame, q2centers : list[float] = Q2CENTERS, **kwargs) -> pd.DataFrame:
    df_rlrt = []
    for q2center in q2centers:
        df_rlrt.append(extract_response_q2bin_excenter(df=df, q2center=q2center, **kwargs))
        df_rlrt.append(extract_response_q2bin_w2center(df=df, q2center=q2center, **kwargs))
    df_rlrt = pd.concat(df_rlrt)
    return df_rlrt

def extract_photo_production_rt_qvbin(qv : float) -> pd.DataFrame:
    photon = pd.read_excel("Carbon/12C_other_extractions/photo-production/Photon_plotting.xlsx")
    photon = photon.sort_values(by = ['nu'])
    photon = photon.reset_index(drop=True)
    photon_index = (photon['nu'] - qv).abs().idxmin()
    interpolation_rows = photon.iloc[photon_index-5 : photon_index+6]
    nu = np.array(interpolation_rows['nu'])
    rt = np.array(interpolation_rows['RT'])
    rterr = np.array(interpolation_rows['error'])
    params, covariance = curve_fit(linear_model, nu, rt,sigma = rterr, absolute_sigma=True)
    a, b = params
    a_err, b_err = np.sqrt(np.diag(covariance))
    qv_rt = linear_model(qv, a, b)
    qv_rterr = np.sqrt((a_err * qv)**2 + b_err**2)
    df = pd.DataFrame({'qv':[qv], 'nu':[qv], 'rt':[qv_rt], 'rterr':[qv_rterr], 'dataset':['Photo-production']})
    return df
