import networkx as nx
import pandas as pd
import community
import matplotlib.pyplot as plt
from collections import Counter

colormap =['black', 'darkgreen', 'royalblue', 'darkorange', 'lime', 'slategrey', 'navy', 'red','blue', 'darkred','teal','saddlebrown',
           'olive','gold','cyan','magenda','darkorchid','pink']
def rn(G,thres):
    for ele in list(G.nodes()):
        if G.degree(ele,weight='weight') <= thres:
            G.remove_node(ele)
    return G

def weightcore(G_in):
    strth=[]
    numnode=[]
    G = G_in.subgraph(list(G_in.nodes())).copy()
    while len(G.nodes())>0:
        deg = [G.degree(ele,weight='weight') for ele in G.nodes()]
        strth.append(min(deg))
        numnode.append(len(deg))
        
        thr = min(deg)
        nl_1 = len(G.nodes())
        
        G_copy = G.subgraph(list(G.nodes())).copy()
        G = rn(G, thr)
        while len(G.nodes())<nl_1:
            nl_1=len(G.nodes())
            G=rn(G,thr)
        if len(G.nodes())>0:
            G = G
        else:
            return thr,G_copy

def save_lcc(el, savepath):
    d_temp = pd.read_csv(el, delimiter=',')
    col = d_temp.columns
    d_temp = d_temp.values.tolist()
    G_temp = nx.Graph()    
    if len(col)>2:
        for e in d_temp:
            G_temp.add_edge(e[0], e[1], weight= round(float(e[2]),3))
            
        gcc = max(nx.connected_components(G_temp), key=lambda x: len(x))
        H_temp = G_temp.subgraph(gcc) 
        if len(H_temp.nodes())<len(G_temp.nodes()):
            save_el = [[u,v,c] for (u, v, c) in H_temp.edges.data('weight')]
            save_el = pd.DataFrame(save_el)
            save_el.columns = ['source', 'target','weight']
            save_el.to_csv( savepath + ' (LCC).csv', sep=',', index=False)  
    else:
        for e in d_temp:
            G_temp.add_edge(e[0], e[1])
            
        gcc = max(nx.connected_components(G_temp), key=lambda x: len(x))
        H_temp = G_temp.subgraph(gcc) 
        if len(H_temp.nodes())<len(G_temp.nodes()):
            save_el = [[u,v] for (u, v) in H_temp.edges()]
            save_el = pd.DataFrame(save_el)
            save_el.columns = ['source', 'target']
            save_el.to_csv( savepath + ' (LCC).csv', sep=',', index=False)      
   
            
class networkanalysis:
    def __init__(self, edgelist='network.csv', weighted=False):
        self.edgelist = edgelist
        self.weighted = weighted
        self.data = pd.read_csv(edgelist, delimiter=',')
        self.data = self.data.values.tolist()
        self.G = nx.Graph()
        
        if not self.weighted:
            for e in self.data:
                self.G.add_edge(e[0],e[1])
        else: 
            for e in self.data:
                self.G.add_edge(e[0], e[1], weight=float(e[2]), inv_weight = round(1/float(e[2]), 3))

        self.corenode =[]
    
        self.savepath = edgelist.split('.')[0]
        
    def find_core(self):
        print('determining core...')
        if not self.weighted:
            coreNet = nx.k_core(self.G)
            deg = [coreNet.degree(nd) for nd in coreNet.nodes()]
            thr = min(deg)
        else:
            thr,coreNet = weightcore(self.G)   
            
        return thr, coreNet
    
    def find_pos(self, lcc=True, core=True, louvain_community=True):
        
        if lcc:
            print('Saving the largest connected component...')
            save_lcc(self.edgelist, self.savepath)
            
        gcc = max(nx.connected_components(self.G), key=lambda x: len(x))
        self.H = self.G.subgraph(gcc)  
        if lcc:
            self.G = self.H  

            print('number of nodes in the largest connected component: ', len(self.G.nodes()))
            print('number of edges in the largest connected component: ', len(self.G.edges())) 
            
        else:    
            print('number of nodes: ', len(self.G.nodes()))
            print('number of edges: ', len(self.G.edges())) 
            
        self.attr=['Id']
        self.nodelist = list(self.G.nodes())
        self.master=[self.nodelist]  
       
        
        if core or louvain_community:
            self.pos= nx.fruchterman_reingold_layout(self.H)
            
        
            
    def netcommu(self,louvain_community = True):
        if louvain_community:
            print('determining communities...')
            comm = community.best_partition(self.G)
            commus = [comm[nd] for nd in self.nodelist]
            self.attr.append('community')
            self.master.append(commus)
            
            Hcommus = [comm[nd] for nd in self.H.nodes()]
            c_cnt = Counter(Hcommus) #{0:1309, 1:1008, 2:12, ..., 20:10, 23: 20}
            commu2cnt = {k: v for k, v in sorted(c_cnt.items(), key=lambda item: item[1], reverse = True)} #{0:1309, 1:1008, 23: 20, 2:12, 20:10, ...}
            print('the sizes of the LCC communities:', commu2cnt)
            
            comm_remap = {c:idx for idx,c in enumerate(commu2cnt.keys())}#{0:0, 1:1, 2:2, 23:3, 2:4, 20:5 ...}
            Hcommus = [comm_remap[e] for e in Hcommus]
            commu2color = []
            # if len(commu2cnt)<=len(colormap):             
            #     commu2color = [colormap[idx] for idx in Hcommus]
            # else:
            keys = {idx:val for idx,val in enumerate(colormap)}
            for ele in Hcommus:
                if ele in keys:
                    commu2color.append(keys[ele])
                else:
                    commu2color.append('lightgray')
                    
            plt.figure(2)
            nx.draw_networkx_nodes(self.H, self.pos, node_size=10,  node_color=commu2color, alpha=0.5)
            # nx.draw_networkx_nodes(H, self.pos, node_size=10, alpha=0.5)
            nx.draw_networkx_edges(self.H, self.pos, edge_color="lightgray")
            plt.axis('off')
            # plt.legend()
            #plt.text(0.1,0.9, str(max(commus)+1)+' communities ')
            plt.tight_layout()
            plt.savefig('community.png', dpi=300)
            
    def netcore(self, core =True):
        if core:
            thr, coreNet = self.find_core()
            self.corenode =  list(coreNet.nodes())
            
            cores =[0]*len(self.nodelist)
            core2pos ={}
            for idx,val in enumerate(self.nodelist):
                if val in self.corenode:
                    cores[idx] = 1 
                    core2pos[val] = self.pos[val]

            self.attr.append('core')
            self.master.append(cores)
            plt.figure(3)
            # pos = nx.spectral_layout(self.G)
            nx.draw_networkx_nodes(self.H, self.pos, node_size=10,  node_color = 'lightgrey')
            nx.draw_networkx_nodes(coreNet, core2pos, node_size=10,  node_color = 'blue')
            nx.draw_networkx_edges(self.H, self.pos, edge_color="lightgray")
            plt.axis('off')
            plt.tight_layout()
            plt.savefig('core.png', dpi=300) 
                
    def netdeg(self, degree=True, eigenvector=True):
           
        if degree:
            print('calculating node degree...')
            if not self.weighted:
                deg = [self.G.degree(nd) for nd in self.nodelist]
                self.attr.append('degree')
            else:
                deg = [self.G.degree(nd, weight = 'weight') for nd in self.nodelist]
                self.attr.append('strength')
            self.master.append(deg)
            
            plt.figure(4, dpi=300)
            plt.hist(deg, bins=100,  edgecolor='black', color= 'blue')
            plt.xlabel(self.attr[-1], fontsize=10)
            plt.ylabel('# of nodes', fontsize=10)
            plt.xticks(fontsize=8)
            plt.yticks(fontsize=8)
            plt.xscale('log')
            plt.yscale('log')
            plt.tight_layout()
            
            plt.savefig('degree_distribution.png', dpi=300)
           
        if eigenvector:
            print('calculating node eigenvector...')
            if not self.weighted:
                eig = nx.eigenvector_centrality(self.G)
            else:
                eig = nx.eigenvector_centrality(self.G, weight = 'weight')
            eigen = [eig[nd] for nd in self.nodelist]
            self.attr.append('eigenvector')
            self.master.append(eigen)
            
            
    def netclose(self, closeness=True):            
        if closeness:
            print('calculating node closeness...')
            if not self.weighted:
                hc = nx.harmonic_centrality(self.G)
            else:
                hc = nx.harmonic_centrality(self.G, distance='inv_weight')
            close = [hc[nd] for nd in self.nodelist]
            self.attr.append('closeness')
            self.master.append(close)

    def netbet(self, betweenness =True):
        if betweenness:
            print('calculating node betweenness...')
            if not self.weighted:
                betc = nx.betweenness_centrality(self.G)
            else:
                betc = nx.betweenness_centrality(self.G, weight= 'inv_weight')
            bet = [betc[nd] for nd in self.nodelist]
            self.attr.append('betweenness')
            self.master.append(bet)
            
    def saveattr(self, attrlist ='attribute_list.csv'):         
        if len(self.master)>1:
            self.master = list(map(list, zip(*self.master)))
            self.master = pd.DataFrame(self.master)
            self.master.columns = self.attr
            
            self.master.to_csv( self.savepath+ ' property list.csv', sep=',', index=False)
            print('network properties are save as:')
            print(self.savepath+ ' property list.csv')
        
        
    def statistics(self):
        if len(self.corenode)>0:
            return len(self.nodelist),len(self.G.edges()),len(self.corenode)
        else:
            return len(self.nodelist),len(self.G.edges())
       
            