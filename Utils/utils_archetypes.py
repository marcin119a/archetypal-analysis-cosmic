import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


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
            denovo_exp=pd.read_csv(path+str(num)+'_'+str(sample)+'/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt',sep='\t')
            real_cat=pd.read_csv(path+str(num)+'_'+str(sample)+'/SBS96/Samples.txt',sep='\t')
            pred_cat=pd.DataFrame(denovo_sig.iloc[:,1:].values @ denovo_exp.iloc[:,1:].T.values)
            pred_cat.insert(0,'Mutation Types',real_cat['Mutation Types'])
            pred_cat.columns=real_cat.columns
            
            detected_sig=denovo_sig.iloc[:,1:].shape[1]
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
    iqr_stab=pd.DataFrame(stability_df.groupby('Samples')['Mean Stability'].apply(lambda x: [np.round(np.percentile(x,25),2),np.round(np.percentile(x,75),2)])).rename(columns={'Median Stability':'IQR'})
    stability=pd.concat([stab,iqr_stab],axis=1)
    
    # median of min stability +iqr
    min_stab=pd.DataFrame(min_stability_df.groupby('Samples')['Min Stability'].apply(lambda x: np.round(np.median(x),2)))
    min_iqr_stab=pd.DataFrame(min_stability_df.groupby('Samples')['Min Stability'].apply(lambda x: [np.round(np.percentile(x,25),2),np.round(np.percentile(x,75),2)])).rename(columns={'Min Stability':'IQR'})
    min_stability=pd.concat([min_stab,min_iqr_stab],axis=1)
    
    return pd.concat([p,mse_tot,stability,min_stability],axis=1).sort_index().replace(np.nan,0)



def plot_extraction():

    f, axes = plt.subplots(2, 3, figsize=(18, 15), sharex=False)
    for ax, scenario in zip(axes.flat,keys):
        ax.grid()
        ax.set_yticks(np.arange(0,11))

        for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(12) 

        for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(12) 

        ax.set_title('\n Scenario {} \n'.format(scenario))
        sns.countplot(dodge=True,x='N° of detected Signatures',hue='Samples',data=scenarios[scenarios.Scenario==scenario],ax=ax,edgecolor='black',palette='deep')
        ax.set_ylabel('N° of runs',fontsize=12)
        ax.set_xlabel('N° of detected Signatures \n',fontsize=12)

        sns.despine(left=True)

    axes.flat[2].legend(loc='upper right')
    axes.flat[4].legend(loc='upper left')

    #fig1
    patchess=axes.flat[0].patches
    width = patchess[0].get_width()
    new_width = 0.1


    patchess[0].set_width(0.1)
    patchess[0].set_label('4')
    x = patchess[0].get_x()
    patchess[0].set_x(x + (width - new_width) +.05)

    patchess[3].set_width(0.1)
    patchess[3].set_label('6')
    x = patchess[3].get_x()
    patchess[3].set_x(x + (width - 0.5) +.1)

    patchess[1].set_width(0.1)
    x = patchess[1].get_x()
    patchess[1].set_x(x + (width - 0.2) +.1)

    #fig2
    patchess=axes.flat[1].patches
    width = patchess[0].get_width()
    new_width = 0.1

    patchess[0].set_width(0.1)
    x = patchess[0].get_x()
    patchess[0].set_x(x + (width - new_width) +.05)

    patchess[3].set_width(0.1)
    x = patchess[3].get_x()
    patchess[3].set_x(x + (width - 0.5) +0.1)

    patchess[1].set_width(0.1)
    x = patchess[1].get_x()
    patchess[1].set_x(x + (width - 0.2) +0.1)

    #fig3

    patchess=axes.flat[2].patches


    for i in [30,32,33,34]:
        patchess[i].set_width(0.16)
        x = patchess[i].get_x()
        patchess[i].set_x(x-0.3)

    #fig4
    patchess=axes.flat[3].patches

    for i in [0,1,2,4]:
        patchess[i].set_width(0.35)
        x = patchess[i].get_x()
        patchess[i].set_x(x+0.22)

    patchess[-1].set_width(0.35)
    x = patchess[-1].get_x()
    patchess[-1].set_x(x-0.22)

    patchess[12].set_width(0.35)
    x = patchess[12].get_x()
    patchess[12].set_x(x)

    patchess[10].set_width(0.35)
    x = patchess[10].get_x()
    patchess[10].set_x(x)


    #fig5
    patchess=axes.flat[4].patches

    for i in [0,1,2]:
        patchess[i].set_width(0.35)
        x = patchess[i].get_x()
        patchess[i].set_x(x+0.22)
    patchess[-1].set_width(0.35)
    x = patchess[-1].get_x()
    patchess[-1].set_x(x-0.3)

    axes[-1, -1].axis('off')

    axes.flat[1].legend(title="Samples", loc="upper left") 
    axes.flat[2].legend(title="Samples", loc="upper right") 
    axes.flat[3].legend(title="Samples", loc="upper left") 
    axes.flat[4].legend(title="Samples", loc="upper left")
    
    plt.show()
