import pandas as pd
import numpy as np
import SigProfilerExtractor
from SigProfilerExtractor import sigpro
import time

index=['A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T',
       'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T',
       'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T',
       'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T',
       'A[C>G]A', 'A[C>G]C','A[C>G]G', 'A[C>G]T',
       'C[C>G]A', 'C[C>G]C','C[C>G]G', 'C[C>G]T',
       'G[C>G]A', 'G[C>G]C','G[C>G]G', 'G[C>G]T',
       'T[C>G]A', 'T[C>G]C','T[C>G]G', 'T[C>G]T',
       'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T',
       'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T',
       'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T',
       'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T',
       'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T',
       'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T',
       'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T',
       'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T', 
       'A[T>C]A', 'A[T>C]C','A[T>C]G', 'A[T>C]T',  
       'C[T>C]A', 'C[T>C]C','C[T>C]G', 'C[T>C]T',
       'G[T>C]A', 'G[T>C]C','G[T>C]G', 'G[T>C]T',
       'T[T>C]A', 'T[T>C]C','T[T>C]G', 'T[T>C]T',
       'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T',
       'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T',
       'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T',
       'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T']

COSMIC_sig=pd.read_csv('./utils/COSMIC_SBS_GRCh37.txt',sep='\t') 



scenarios={1:{'samples':[200,500],'path_from':'./Scenario_1/flat_','path_to': 'SigProfiler_results/Scenario_1/flat_'},
          2:{'samples':[200,500],'path_from':'./Scenario_2/5_sim_no_flat_','path_to': 'SigProfiler_results/Scenario_2/5_sim_no_flat_'},
          3:{'samples':[200,500,1000],'path_from':'./Scenario_3/eleven_sim_no_flat_','path_to': 'SigProfiler_results/Scenario_3/11_sim_no_flat_'},
          4:{'samples':[200,500,1000,3000,5000],'path_from':'./Scenario_4/both_','path:to': 'SigProfiler_results/Scenario_4/both_'},
          5:{'samples':[1000,3000,5000],'path_from':'./Scenario_5/all_sim_no_flat_','path_to': 'SigProfiler_results/Scenario_5/all_sim_no_flat_'}}          


if __name__ == '__main__':

       for scenario in range(1,6):
              catalogues=scenarios[scenario]['samples']
              path_from=scenarios[scenario]['path_from']
              path_to=scenarios[scenario]['path_to']
              mutations=[5000]

              for run in range(1,11):
                     for catalogue in catalogues:
                            print(f'Extraction with {catalogue} samples')
                            synth_cat=pd.read_csv(path_from + f"{run}_{catalogue}.csv").iloc[:,1:]
                            synth_cat.index=index
                            file_name=path_to + f"{run}_{catalogue}"

                            try:
                                   s=sigpro.sigProfilerExtractor('matrix', file_name , synth_cat.loc[COSMIC_sig['Type']].reset_index(), reference_genome="GRCh37", opportunity_genome = "GRCh37", context_type = "default", exome = False, 
					minimum_signatures=2, maximum_signatures=22, nmf_replicates=30, resample = True, batch_size=1, cpu=-1, gpu=False, nmf_init="random", precision= "single", matrix_normalization= "gmm", 
					seeds= "random", cosmic_version=3.2,min_nmf_iterations= 10000, max_nmf_iterations=1000000, nmf_test_conv= 10000, nmf_tolerance= 1e-8,nnls_add_penalty=0.05,
					nnls_remove_penalty=0.01, initial_remove_penalty=0.05, de_novo_fit_penalty=0.02,get_all_signature_matrices=True)
                            
                            except:
                                   
                                   continue
