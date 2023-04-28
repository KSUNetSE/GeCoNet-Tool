import numpy as np
import pandas as pd
import csv
from collections import defaultdict
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import networkx as nx 

def save_lcc():
    pass
class GeneCoExp():
    def __init__(self, corr='', paired_items='', filepath='expression.csv'):
        self.filepath = filepath
        self.N = 0
        self.M = 0
        self.corr = corr
        if corr and paired_items:
            print('reading correlation table...')
            self.rho = pd.read_csv(corr, delimiter=',', index_col=0)
            self.paired_items = pd.read_csv(paired_items, delimiter=',', index_col=0)
            self.rho.columns = [int(ele) for ele in self.rho.columns]
            self.paired_items.columns = [int(ele) for ele in self.paired_items.columns]
        
    def read_data(self, zero_removed=False, zscored=False,  rescaled=False, dropna=0, savefile=False):
        '''
        Parameters
        ----------
        zero_removed : required
            remove zeros from the table. The default is False.
        zscored : required
            z-score normalize each column. The default is False.
        rescaled : required
            log2 rescale each value. The default is False.
        savefile : optional
            save the normalized data. The default is ''.
        '''
        print('normalizing data...')
        self.data = pd.read_csv(self.filepath, delimiter=',', index_col=0)
        if self.N != len(self.data): self.N = len(self.data)
        if self.M != len(self.data.columns): self.M = len(self.data.columns)
        
        # self.data[self.data< zero_removed] = np.nan
        if not zero_removed:
            self.data.replace(0, np.nan, inplace=True)#remove 0s 

        if not zscored:  # zscore normalize the data
            if not rescaled: #rescale the expression data, if input is not rescaled
                for col in self.data.columns:
                    self.data[col] = np.log2(self.data[col])
                    
            for col in self.data.columns: # zscore normalize the data
                self.data[col] = (self.data[col] - self.data[col].mean())/self.data[col].std(ddof=0)    
                
                assert abs(self.data[col].mean())<1e-1 and abs(self.data[col].var()-1.)<1e-1
        
        self.data = self.data.dropna(thresh = dropna, axis='columns')
        print(self.data.shape)
        self.savepath = self.filepath.split('.')[0]
        if savefile:
            self.data.to_csv(self.savepath+' (z-scored).csv', sep=',')
            print('normalized data are saved as: ')
            print(self.savepath + ' (z-scored).csv')
        
    def Calculate_PCC(self, min_periods=4,  pcc_path=False, paired_elements_path=False):
        '''
        Parameters
        ----------
        min_periods : int, optional
            Minimum number of observations required per pair of columns to have a valid result. The default is 4.
        gene2index : string, required
            corresponds each gene with a number. The default is 'gene2index.csv'.
        pcc_path : string.csv, optional
            save the PCC matrix table. The default is ''.
        paired_elements_path : string.csv, optional
            save the number of observations between each pair of columns. The default is ''.
        '''
        self.min_periods = min_periods #the min # of paired elements
        if self.corr:
            pass
        else:        
            #calculate PCC

            print('calculating PCC...')
            self.rho = self.data.T.astype('float16').corr(min_periods=self.min_periods)   #compute the correlation coefficient
            self.nodelist = list(self.rho.columns)
            print(self.rho.shape)
            self.rho = self.rho.round(3)
            self.rho = pd.DataFrame(np.triu(self.rho.to_numpy(), 1))
            #save the PCC matrix and paired elements matrix
            if pcc_path:
                self.rho.to_csv(self.savepath + ' PCC_matrix.csv', sep=',')
                print('PCC table is saved as:')
                print(self.savepath + ' PCC_matrix.csv')
                
            self.triu = np.triu_indices(len(self.rho), 1) # the upper triangular matrix
            self.rho =  self.rho.to_numpy()[self.triu] #the upper triangluar of the PCC matrix
            
            
    def paired_elements(self, paired_elements_path = False):
        if self.corr:
            pass
        else:  
            print('determining the # of paired elements for every node pair...')
            #find the # of paired elements between every pair of genes
            val_df = self.data.notnull().astype('int').to_numpy() #replace non-empty with 1, and empty with 0   
            self.data =None
            self.paired_items = np.dot(val_df, val_df.T) #(N, M) * (M, N) -> (N, N)
            self.paired_items = pd.DataFrame(np.triu(self.paired_items, 1))
    
            #save the paired elements matrix
            if paired_elements_path:
                self.paired_items.to_csv(self.savepath + ' paired_elements_matrix.csv', sep =',')
                print('paired element table is saved as:')
                print(self.savepath + ' paired_elements_matrix.csv')
            self.paired_items = self.paired_items.to_numpy()[self.triu] #the upper triangular of paired elements matrix 
            
    def calculate_threshold(self, bin_size=10, cutoff=0.005):
        '''
        Parameters
        ----------
        bin_size : 10 or 5, required
            divide the PCC into different intervals per the # of paired elements. The default is 10.
        cutoff : float, required
            choose the top fraction PCC to construct network. The default is 0.005.
        '''
        #self.rho, self.paired_items
        print('calculating threshold...')
        pair2pcc = defaultdict(list) #{4:[0.3, 0.5,...], 5:[0.1, 0.5, ...],...}
        
        for idx in range(len(self.paired_items)): #
            if self.paired_items[idx] >= self.min_periods:
                pair2pcc[self.paired_items[idx]].append(self.rho[idx]) #
        
        pairs = sorted(list(pair2pcc.keys())) 
        
        print('determining the threshold of each interval...')
        bin2pair = defaultdict() # {4:7, 5:7,...,101:105, 102:105}
        # for p in pairs:
        #     bin2pair[p] = ((p-1)//bin_size)*bin_size + (bin_size+1)/2.
        #     if p <=10:
        #         bin2pair[p] = (min(pairs)+10)/2                         
        #     if max(pairs)//bin_size==((p-1)//bin_size):
        #         bin2pair[p] = (((p-1)//bin_size)*bin_size + max(pairs)+1)/2.
        for p in pairs:
            bin2pair[p] = ((p - self.min_periods)//bin_size)*bin_size + self.min_periods + (bin_size -1)/2.
            
        bin2pcc = defaultdict(list) # {7:[pcc], 15:[pcc], 25:[pcc]}
        for key in pairs:
            if len(pair2pcc[key])>0:
                bin2pcc[bin2pair[key]] += pair2pcc[key]
        
        bin2size={}
        for key,val in bin2pcc.items():
            bin2size[key] = len(val)
        print('The size of each bin: ', bin2size)
        
        self.bin2cut = {}
        for key,val in bin2pcc.items():
            cut_val = sorted(val)[int(len(val)*(1-cutoff))]
            if len(val)>50 and cut_val > 0.3:
                self.bin2cut[key] = cut_val
        print('threshold of each interval: ', self.bin2cut)
        
    def curve_fitting(self, alpha=1,eta=1.5,lam=2,beta=30, para_in= True, edgelist=True, thresholdcurve=True):
        '''
        Parameters
        ----------
        p0 : list, required
            the initial values of the 4 parameter to fit the curve. The default is [1,1.5,2,30].
        edgelist : string, required
            save the edge list. The default is 'edgelist.csv'.
        thresholdcurve : string, optional
            save the fitted threhold curve. The default is 'threshold_curve.png'.
        '''
        if len(self.bin2cut)<4:
            thresh = np.mean([val for key,val in self.bin2cut.items()])
            el = np.where(self.rho >=thresh, self.rho, 0)
        else:
            p0=[alpha,eta,lam,beta]
            xdata=[]
            ydata=[]
            for key,val in self.bin2cut.items():
                xdata.append(key)
                ydata.append(val)
            xdata = np.array(xdata,dtype=np.float64)
            ydata = np.array(ydata,dtype=np.float64)
            self.popt, pcov =  curve_fit(f=func, xdata=xdata, ydata=ydata, p0=p0, maxfev=20000)
            
            if not para_in:
                self.popt = np.array([alpha,eta,lam,beta])
            self.r2 = r2_score(ydata, func(xdata, *self.popt))
            
            if thresholdcurve:
                plt.figure(1)
                plt.scatter(xdata, ydata, color ='red', label = 'raw data')
                plt.plot(xdata, func(xdata, *self.popt), color='blue', label ='fitted curve')
                plt.xlabel('# of paired elements', fontsize =15)
                plt.ylabel('threshold', fontsize =15)
                plt.xticks(fontsize=15)
                plt.yticks(fontsize=15)
                plt.legend(fontsize = 15)
                plt.tight_layout()
                plt.savefig('threshold_curve.png', dpi=300)
                plt.savefig(self.savepath + ' threshold_curve.png', dpi=300)
    
            thres_max = func(min(xdata), *self.popt)
            thres_min = func(max(xdata), *self.popt)
            print('parameters of the curve:')
            print('    alpha: ',self.popt[0])
            print('    eta: ',self.popt[1])
            print('    lambda: ',self.popt[2])
            print('    beta: ',self.popt[3])
            print('generating edge list...')
            thres_matrix = func(self.paired_items, *self.popt)
            el = np.where(self.rho >= thres_matrix, self.rho, 0)
            #el = np.where(self.rho <0.999, self.rho, 0)
            
            thres_matrix = thres_max + thres_min - thres_matrix
            el = np.multiply(el, thres_matrix)
        self.rho = None
        self.paired_items = None
		
        el = np.round(el,decimals = 3)
        
        el_index = np.nonzero(el)
        source = [self.nodelist[self.triu[0][ele]] for ele in el_index[0]]
        target = [self.nodelist[self.triu[1][ele]] for ele in el_index[0]]
        
        self.n_edges = len(source)
        self.n_nodes = len(set(source + target))
        el = map(list, zip(*[source, target, el[el_index[0]]]))
        if edgelist:
            el = pd.DataFrame(el, columns=['source','target','weight'])
            el.to_csv(self.savepath + ' edgelist.csv', sep=',', index=False)
            print('edge list is saved as:')
            print(self.savepath + ' edgelist.csv')


def func(x, alpha, eta, lam, beta):
    return alpha - 1/(eta+lam*np.exp(-x/beta))



