import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def extract_from_catalogues(samples,path,real_signatures=5,runs=10):
    
    all_detected={}
    all_mse={}
    all_stability={}
    all_min_stability={}

    for sample in samples:    
        detected_sigs=[]
        mses=[]
        mean_stability=[]
        min_stability=[]
        
        for num in range(1,runs+1):
            
            denovo_sig=pd.read_csv(path+str(num)+'_'+str(sample)+'/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt',sep='\t')
            detected_sig=denovo_sig.iloc[:,1:].shape[1]

            denovo_exp=pd.read_csv(path+str(num)+'_'+str(sample)+'/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt',sep='\t')
           
            real_cat=pd.read_csv(path+str(num)+'_'+str(sample)+'/SBS96/Samples.txt',sep='\t')
            pred_cat=pd.DataFrame(denovo_sig.iloc[:,1:].values @ denovo_exp.iloc[:,1:].T.values)
            pred_cat.insert(0,'MutationType',real_cat['MutationType'])
            pred_cat.columns=real_cat.columns
            denovo_sig.to_csv(f'./Signatures_extracted/Sig_{num}_{sample}.csv')
            denovo_exp.to_csv(f'./Signatures_extracted/Exp_{num}_{sample}.csv')
            mean_stability.append(np.mean(pd.read_csv(path+str(num)+'_'+str(sample)+'/SBS96/All_Solutions/SBS96_'+str(detected_sig)+'_Signatures/Solution_Stats/SBS96_S'+str(detected_sig)+'_Signatures_stats.txt',sep='\t')['Stability']))
            min_stability.append(np.min(pd.read_csv(path+str(num)+'_'+str(sample)+'/SBS96/All_Solutions/SBS96_'+str(detected_sig)+'_Signatures/Solution_Stats/SBS96_S'+str(detected_sig)+'_Signatures_stats.txt',sep='\t')['Stability']))       
            mse_tot=np.mean(np.mean(np.square(real_cat.iloc[:,1:] - pred_cat.iloc[:,1:]),axis=None))
            detected_sigs.append(detected_sig)
            mses.append(mse_tot)
            
        all_detected[str(sample)]=detected_sigs
        all_mse[str(sample)]=mses
        all_stability[str(sample)]=mean_stability
        all_min_stability[str(sample)]=min_stability
    
    # mse
    mse=pd.DataFrame(all_mse)
    mse_df=pd.DataFrame(pd.concat([mse[str(sample)] for sample in samples]))
    mse_df.index=sorted(samples*runs)
    mse_df.reset_index(inplace=True)
    mse_df.columns=['Samples','MSE']
    
    #detected signatures
    detected_df=pd.DataFrame(all_detected)
    detected_df=pd.DataFrame(pd.concat([detected_df[str(sample)] for sample in samples]))
    detected_df.index=sorted(samples*runs)
    detected_df.reset_index(inplace=True)
    detected_df['N° of true Signatures']=real_signatures
    detected_df.columns=['Samples','N° of detected Signatures','N° of true Signatures']
    
    # mean stability
    stability=pd.DataFrame(all_stability)
    stability_df=pd.DataFrame(pd.concat([stability[str(sample)] for sample in samples]))
    stability_df.index=sorted(samples*runs)
    stability_df.reset_index(inplace=True)
    stability_df.columns=['Samples','Mean Stability']
    
    # min stability 
    min_stability=pd.DataFrame(all_min_stability)
    min_stability_df=pd.DataFrame(pd.concat([min_stability[str(sample)] for sample in samples]))
    min_stability_df.index=sorted(samples*runs)
    min_stability_df.reset_index(inplace=True)
    min_stability_df.columns=['Samples','Min Stability']    
    return stability_df,min_stability_df,mse_df,detected_df


def statistics_results(detected_df,mse_df,stability_df,min_stability_df):
    
    # percentage
    p=pd.DataFrame(detected_df[detected_df["N° of detected Signatures"]==detected_df["N° of true Signatures"]].groupby('Samples')["N° of true Signatures"].count()/10*100).rename(columns={"N° of true Signatures":'P'})
    
    # median mse + iqr
    mse=pd.DataFrame(mse_df.groupby('Samples').MSE.apply(lambda x : np.round(np.median(x),2)))
    iqr=pd.DataFrame(mse_df.groupby('Samples').MSE.apply(lambda x: [np.round(np.percentile(x,25),2),np.round(np.percentile(x,75),2)])).rename(columns={'MSE':'IQR'}) #- np.percentile(x,[25,75])[0])).rename(columns={'MSE':'IQR'})
    mse_tot=pd.concat([mse,iqr],axis=1)
    
    # median of median stability +iqr
    stab=pd.DataFrame(stability_df.groupby('Samples')['Mean Stability'].apply(lambda x: np.round(np.median(x),2)))
    iqr_stab=pd.DataFrame(stability_df.groupby('Samples')['Mean Stability'].apply(lambda x: [np.round(np.percentile(x,25),2),np.round(np.percentile(x,75),2)])).rename(columns={'Mean Stability':'IQR'})
    stability=pd.concat([stab,iqr_stab],axis=1)
    
    # median of min stability +iqr
    min_stab=pd.DataFrame(min_stability_df.groupby('Samples')['Min Stability'].apply(lambda x: np.round(np.median(x),2)))
    min_iqr_stab=pd.DataFrame(min_stability_df.groupby('Samples')['Min Stability'].apply(lambda x: [np.round(np.percentile(x,25),2),np.round(np.percentile(x,75),2)])).rename(columns={'Min Stability':'IQR'})
    min_stability=pd.concat([min_stab,min_iqr_stab],axis=1)
    
    return pd.concat([p,mse_tot,stability,min_stability],axis=1).sort_index().replace(np.nan,0)



def plot_extraction(scenarios):

    f, axes = plt.subplots(2, 3, figsize=(25, 20), sharex=False)
    for ax, scenario,true_sig in zip(axes.flat,[1,2,3,4,5],[6,5,11,11,20]):
        ax.grid()
        #ax.set_yticks(np.arange(0,11))

        for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(14) 

        for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(14) 

        ax.set_title(f'\n Scenario {scenario} \n (True n° signatures = {true_sig})\n',fontsize=14)
        sns.swarmplot(y='N° of detected Signatures',x='Samples',size=9,data=scenarios[scenarios.Scenario==scenario],ax=ax,edgecolor='black',palette='deep')
        ax.set_ylabel('N° of detected Signatures',fontsize=14)
        ax.set_xlabel('Samples',fontsize=14)
    
    axes.flat[0].set_yticks(np.arange(2,12))
    axes.flat[1].set_yticks(np.arange(2,12))
    axes.flat[2].set_yticks(np.arange(2,12))
    axes.flat[3].set_yticks(np.arange(2,12))
    axes.flat[4].set_yticks(np.arange(4,22,2))
    
    plt.subplots_adjust(hspace=0.35,wspace=0.25)    
    axes[-1, -1].axis('off')
    plt.show()
    
    
def plot_archetypes(array,label,ax, ylim=1,archetypes=0):
    
    color = ((0.196,0.714,0.863),)*16 + ((0.102,0.098,0.098),)*16 + ((0.816,0.180,0.192),)*16 + ((0.777,0.773,0.757),)*16 +   ((0.604,0.777,0.408),)*16 + ((0.902,0.765,0.737),)*16
    color = list(color)

    width = max(array.shape)
    x = np.arange(width)
    if ax == None:
        f,ax = plt.subplots(1,figsize=(20,10))
    bars = ax.bar(x, array)

    for h in range(len(x)):
        bars[h].set_color(color[h])


    plt.title('Archetypes '+str(archetypes+1))
    plt.ylim(0, ylim)
    plt.xlim(0, width)    
