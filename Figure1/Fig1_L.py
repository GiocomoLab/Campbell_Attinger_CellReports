import umap
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.stats 
from scipy.ndimage import gaussian_filter1d
import sys
sys.path.append('./python_helpers')
import loadmat as lm
import helpers as helpers
from sklearn import linear_model
# %matplotlib widget
from sklearn.decomposition import PCA
import glob
from sklearn.cluster import DBSCAN,KMeans
import shutil
import traceback



if __name__=='__main__':

    data_dir = 'path_to_data'
    final_files = ['AA1_190727_contrast_1.mat',
    'AA44_190926_gain_1.mat',
    'npF4_1023_gaincontrast_1.mat',
    ]
    files = final_files

   
    

    im_path = os.path.join(data_dir,'umap_baseline_fig')
    
    if not os.path.isdir(im_path):
        os.makedirs(im_path)
    

    shutil.copy2(os.path.abspath(__file__),im_path)
    ds_factor = 5
    for fi in files:
        try:
            fp=os.path.join(data_dir,fi)
            data = lm.loadmat(fp)
              
            
            gain_val = 1
            trial_range = np.arange(5,21)

            try:
                anatomy = data['anatomy']
            except:
                print('no anatomy')
                continue
            
            if 'parent_shifted' in anatomy:
                group = anatomy['parent_shifted']
            else:
                group = anatomy['cluster_parent']
            regions = ('MEC','VISp','RS')
            idx = [str(ss).startswith(regions) for ss in group]
            idxagl = [str(ss).startswith('RSPagl') for ss in group]

            region_idx = np.array(idx) & np.logical_not(np.array(idxagl))

            _,sn = os.path.split(fi)
            
            good_cells = data['sp']['cids'][(data['sp']['cgs']==2) & region_idx]
            if len(good_cells)<5:
                continue

            opt = helpers.options()
            counts,spMapN,stab=helpers.calculateFiringRateMap(data,good_cells=good_cells,trials2extract = trial_range,ops=opt)
            
            X = np.array([])
            tri = np.array([])
            gain = np.array([])
            pos = np.array([])
            for iT in range(counts.shape[2]):
                xs = spMapN[:,:,iT].T
                X = np.vstack([X, xs]) if X.size else xs
                tmp_tri = np.ones((200,1))*iT
                tri = np.vstack([tri,tmp_tri]) if tri.size else tmp_tri
                mult = gain_val if np.isin(iT,range(6,10)) else 1
                tmp_gain = np.ones((200,1))*mult
                gain = np.vstack([gain,tmp_gain]) if gain.size else tmp_gain

                tmp_pos = np.arange(200)
                pos = np.hstack([pos,tmp_pos]) if pos.size else tmp_pos


            m=X.mean(axis=0)
            fr_idx=m>0.1
            Xz=X[:,fr_idx]
            reducer = umap.UMAP(n_components=2,metric='cosine',random_state=42)
            X_um = reducer.fit_transform(Xz)
            

            twoD=True
            fig = plt.figure(figsize=(15,5))
            ax = fig.add_subplot(1,1,1)
            
            bl_idx = np.arange(0,3200)
            
            neg=ax.scatter(X_um[bl_idx,0],X_um[bl_idx,1],c=pos,s=2)
           
            ax.set_aspect('equal', 'box')
            fig.colorbar(neg, ax=ax)
            sn = os.path.join(im_path,sn)
            sn=sn.replace('.mat','.svg')
            fig.savefig(sn)
            fig.clear()
        except ZeroDivisionError:
            print(traceback.format_exc())
            print('did not work for '+sn)
            

