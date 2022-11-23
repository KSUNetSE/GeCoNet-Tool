import tkinter.filedialog
import tkinter as tk
import tkinter.ttk as ttk
import dataload
import threading
import os
import network_analysis as na
from PIL import Image as Pil_image, ImageTk as Pil_imageTk

class App(tk.Frame):           
    def __init__(self, master):
        self.window = master
        self.window.config(bg="skyblue")
        self.window.title("GeCoNet_Tool")
        self.window.maxsize(1050,  590)
        # self.window.geometry('1050x600')
        #self.window.iconphoto(False, tk.PhotoImage(file = 'DNA.png'))
        self.gen_win = tk.Frame(self.window,  width=550,  height=200)
        self.gen_win.pack(side='left', fill='both', padx=5,  pady=5,  expand=True)
        self.res_win = tk.Frame(window,  width=500,  height=200)
        self.res_win.pack(side='right', fill='both', padx=5,  pady=5,  expand=True)
        self.res_win.pack_propagate(0)
        
        data_label = tk.Label(self.gen_win, text="Network Generation & Analysis" , bg='light grey')
        data_label.pack( padx=5,  pady=0, fill=tk.X)
        data_label.config(font=("Font", 12))

        self.data_bar = tk.Frame(self.gen_win,  width=260,  height=200, bg='light grey')
        self.data_bar.pack(side='left',  fill='both',  padx=5,  pady=5,  expand=False)
        self.data_bar.pack_propagate(0)
        lbl = tk.Label(self.data_bar, text="Data Processing", bg='grey')
        lbl.pack( fill=tk.X)
        lbl.config(font=("Font", 10))


        self.input_bar  =  tk.Frame(self.data_bar,  width=200,  height=1)
        self.input_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=10,  expand= True)
        # self.input_bar.pack_propagate(0)
        tk.Label(self.input_bar, text="Input data (.csv):                                  ").pack( padx=5,pady=0,fill=tk.X)
        
        # self.data = tk.Entry(self.input_bar, width=20)
        # self.data.insert(0, 'anopheles.csv')
        # self.data.pack(padx=5,  pady=0, fill=tk.X, expand= False)
        
        self.brow_bar  =  tk.Frame(self.input_bar,  width=200,  height=1)
        self.brow_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= True)
        self.data = tk.Entry(self.brow_bar, width=28)
        self.data.insert(0, 'anopheles.csv')
        self.data.pack(side=tk.LEFT, padx=5,  pady=0, fill=tk.X, expand= False)
        tk.Button(self.brow_bar, text="Open...", command=self.read_file).pack(side=tk.RIGHT,padx=5,pady=0, fill=tk.X)

        self.filter_bar  =  tk.Frame(self.gen_win,  width=220,  height=100, bg='light grey')
        self.filter_bar.pack(side='right',  fill='both',  padx=5,  pady=5,  expand=True)

        self.zero_bar  =  tk.Frame(self.data_bar,  width=15,  height=1)
        self.zero_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= False)
        self.log_bar  =  tk.Frame(self.data_bar,  width=15,  height=1)
        self.log_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= False)
        self.zscore_bar  =  tk.Frame(self.data_bar,  width=15,  height=1)
        self.zscore_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= False)

        self.zero_ = tk.BooleanVar(value=True)
        self.rescaled = tk.BooleanVar(value=True)
        self.zscored = tk.BooleanVar(value=False)
        
        # self.zero_label = Label(self.zero_bar, text="Remove values (smaller than): ").pack(side=TOP, padx=5,pady=0, fill=X)
        # self.zero_ = Entry(self.zero_bar, width=15)
        # self.zero_.insert(0, -1000000000)
        # self.zero_.pack(side=TOP, padx=5, pady=0, fill=X)
        
        tk.Checkbutton(self.zero_bar, text="Remove zeros          ", variable=self.zero_, onvalue=False, offvalue=True, height=1, width=15).pack(side=tk.LEFT,padx=5,pady=0, fill=tk.X)
        tk.Checkbutton(self.log_bar, text="re-scale values by log2", variable=self.rescaled, onvalue=False, offvalue=True, height=1, width=17).pack(side=tk.LEFT,padx=5,pady=0, fill=tk.X)
        tk.Checkbutton(self.zscore_bar, text="z-score columns      ", variable=self.zscored, onvalue=False, offvalue=True, height=1, width=15).pack(side=tk.LEFT,padx=5,pady=0, fill=tk.X)
        
        tk.Label(self.data_bar, text="Drop columns (size <):                         ").pack(side=tk.TOP, padx=5,pady=0, fill=tk.X)
        self.dropna = tk.Entry(self.data_bar, width=15)
        self.dropna.insert(0, 200)
        self.dropna.pack(side=tk.TOP, padx=5, pady=0, fill=tk.X)
        
        self.pcc_bar1  =  tk.Frame(self.data_bar,  width=10,  height=1)
        self.pcc_bar1.pack(side='top',  fill=tk.X,  padx=5,  pady=10,  expand= False)
        self.pcc = tk.BooleanVar()
        self.paired = tk.BooleanVar()
        self.save_data = tk.BooleanVar()
        tk.Checkbutton(self.pcc_bar1, text="Save processed data              ", variable=self.save_data, onvalue=True, offvalue=False, height=1, width=15).pack(padx=5,pady=0, fill=tk.X)
        tk.Checkbutton(self.pcc_bar1, text="Save PCC matrix                    ", variable=self.pcc, onvalue=True, offvalue=False, height=1, width=20).pack(padx=5,pady=0, fill=tk.X)
        tk.Checkbutton(self.pcc_bar1, text="Save paired elements matrix", variable=self.paired, onvalue=True, offvalue=False, height=1, width=20).pack(padx=5,pady=0, fill=tk.X)
        
        # ep  =  Frame(data_bar,  width=10,  height=1)
        # ep.pack(side='top',  fill=X,  padx=5,  pady=5,  expand= False)
        lbl = tk.Label(self.data_bar, text="Network Generation", bg='grey')
        lbl.pack( fill=tk.X)
        lbl.config(font=("Font", 10))
        
        self.curve_params = tk.BooleanVar(value=True)
        tk.Checkbutton(self.data_bar, text="Optimize the sliding threshold parameters:", variable=self.curve_params, onvalue=True, offvalue=False, height=1, width=25).pack(padx=5,pady=0, fill=tk.X)
        
        self.param_bar1  =  tk.Frame(self.data_bar,  width=10,  height=1)
        self.param_bar1.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= False)
        tk.Label(self.param_bar1, text="    Alpha:").pack(side=tk.LEFT, fill=tk.X)
        self.alpha = tk.Entry(self.param_bar1, width=8)
        self.alpha.pack(side=tk.LEFT, padx=5, pady=0, fill=tk.X)
        self.alpha.insert(0, 1.28)
        self.eta = tk.Entry(self.param_bar1, width=8)
        self.eta.pack(side=tk.RIGHT, padx=5, pady=0, fill=tk.X)
        self.eta.insert(0, 1.61)
        tk.Label(self.param_bar1, text="Eta:").pack(side=tk.RIGHT, fill=tk.X)
        self.param_bar2  =  tk.Frame(self.data_bar,  width=10,  height=1)
        self.param_bar2.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand=False)
        tk.Label(self.param_bar2, text="Lambda:").pack(side=tk.LEFT, fill=tk.X)
        self.lam = tk.Entry(self.param_bar2, width=8)
        self.lam.insert(0,2.00)
        self.lam.pack(side=tk.LEFT, padx=5, pady=0, fill=tk.X)
        self.beta = tk.Entry(self.param_bar2, width=8)
        self.beta.insert(0,30.00)
        self.beta.pack(side=tk.RIGHT, padx=5, pady=0, fill=tk.X)
        tk.Label(self.param_bar2, text="Beta:").pack(side=tk.RIGHT, fill=tk.X)
        
        self.param_bar2  =  tk.Frame(self.data_bar,  width=10,  height=1)
        self.param_bar2.pack(side='top',  fill=tk.X,  padx=5,  pady=5,  expand= False)
        self.bin_bar = tk.Frame(self.param_bar2,  width=10,  height=1)
        self.bin_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= False)
        tk.Label(self.bin_bar, text="Bin size:").pack(side=tk.LEFT, padx=10,pady=0, fill=tk.X)
        self.bin_size = tk.Entry(self.bin_bar, width=15)
        self.bin_size.insert(0,10)
        self.bin_size.pack(side=tk.LEFT, padx=10,pady=0, fill=tk.X)
        
        self.cut_bar = tk.Frame(self.param_bar2,  width=10,  height=1)
        self.cut_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= False)
        tk.Label(self.cut_bar, text="  Cutoff:").pack(side=tk.LEFT, padx=10,pady=0, fill=tk.X)
        self.cutoff = tk.Entry(self.cut_bar, width=15)
        self.cutoff.insert(0,0.005)
        self.cutoff.pack(side=tk.LEFT, padx=10,pady=0, fill=tk.X)
        
        self.edge_bar = tk.Frame(self.data_bar,  width=10,  height=1)
        self.edge_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=5,  expand= False)
        self.edgelist = tk.BooleanVar(value=True)
        self.thres_curve = tk.BooleanVar(value=True)
        tk.Checkbutton(self.edge_bar, text="Save edge list             ", variable=self.edgelist, onvalue=True, offvalue=False, height=1, width=15).pack(padx=10,pady=0, fill=tk.X)
        tk.Checkbutton(self.edge_bar, text="Save threshold curve", variable=self.thres_curve, onvalue=True, offvalue=False, height=1, width=15).pack(padx=10,pady=0, fill=tk.X)
        
        
        lbl = tk.Label(self.filter_bar, text="Network Analysis", bg='grey')
        lbl.pack( fill=tk.X)
        lbl.config(font=("Font", 10))
        
        self.input_bar  =  tk.Frame(self.filter_bar,  width=15,  height=1)
        self.input_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=10,  expand= False)
        lbl = tk.Label(self.input_bar, text="Input edgelist (.csv):                        ")
        lbl.pack( fill=tk.X)
        
        self.open_bar  =  tk.Frame(self.input_bar,  width=15,  height=1)
        self.open_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= False)
        
        self.readedge = tk.Entry(self.open_bar, width=28)
        self.readedge.insert(0, 'edgelist.csv')
        self.readedge.pack(side=tk.LEFT, padx=5,  pady=0, fill=tk.X, expand= False)
        tk.Button(self.open_bar, text="Open...", command=self.open_edge).pack(side=tk.RIGHT,padx=5,pady=0, fill=tk.X)

        self.set_bar  =  tk.Frame(self.filter_bar,  width=15,  height=1)
        self.set_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= False)
        
        self.lcc = tk.BooleanVar()
        tk.Checkbutton(self.set_bar, text="Save & Analyze Largest Component", variable=self.lcc, onvalue=True, offvalue=False, height=1, width=15).pack(padx=5,pady=0, fill=tk.X)
        
        self.weighted = tk.BooleanVar()
        tk.Checkbutton(self.set_bar, text="weighted (slow mode)", variable=self.weighted, onvalue=True, offvalue=False, height=1, width=15).pack(padx=5,pady=0, fill=tk.X)
        
        self.pro_bar  =  tk.Frame(self.filter_bar,  width=15,  height=1)
        self.pro_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=10,  expand= False)
        
        self.community = tk.BooleanVar(value=True)
        tk.Checkbutton(self.pro_bar, text="community", variable=self.community, onvalue=True, offvalue=False, height=1, width=15).pack(padx=5,pady=0, fill=tk.X)
        self.core = tk.BooleanVar(value=True)
        tk.Checkbutton(self.pro_bar, text="core", variable=self.core, onvalue=True, offvalue=False, height=1, width=15).pack(padx=5,pady=0, fill=tk.X)
        self.degree = tk.BooleanVar(value=True)
        tk.Checkbutton(self.pro_bar, text="degree", variable=self.degree, onvalue=True, offvalue=False, height=1, width=15).pack(padx=5,pady=0, fill=tk.X)
        self.eigenvector = tk.BooleanVar(value=True)
        tk.Checkbutton(self.pro_bar, text="eigenvector", variable=self.eigenvector, onvalue=True, offvalue=False, height=1, width=15).pack(padx=5,pady=0, fill=tk.X)
        self.closeness = tk.BooleanVar()
        tk.Checkbutton(self.pro_bar, text="closeness", variable=self.closeness, onvalue=True, offvalue=False, height=1, width=15).pack(padx=5,pady=0, fill=tk.X)
        self.betweenness = tk.BooleanVar()
        tk.Checkbutton(self.pro_bar, text="betweenness", variable=self.betweenness, onvalue=True, offvalue=False, height=1, width=15).pack(padx=5,pady=0, fill=tk.X)
        
        self.run_bar  =  tk.Frame(self.data_bar,  width=15,  height=1)
        self.run_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=5,  expand= False)    
        self.stopflag = False
        run_code =tk.Button(self.run_bar, text="START", width=15, command=lambda: threading.Thread(target = self.clicked).start()) #threading.Thread(target=clicked).start()
        run_code.pack(side=tk.LEFT,  padx=5,pady=5, fill=tk.X)
        run_code = tk.Button(self.run_bar, text="STOP", width=15, command=self.stop) #lambda: threading.Thread(target = self.stop).start()
        run_code.pack(side=tk.RIGHT,  padx=5,pady=5, fill=tk.X)

        self.run_bar  =  tk.Frame(self.filter_bar,  width=15,  height=1)
        self.run_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=5,  expand= False) 
        run_code = tk.Button(self.run_bar, text="START", width=15, command=lambda: threading.Thread(target = self.NetAttr).start()) #threading.Thread(target=clicked).start()
        run_code.pack(side=tk.LEFT, padx=5,pady=5, fill=tk.X)
        run_code = tk.Button(self.run_bar, text="STOP", width=15, command=self.stop) #lambda: threading.Thread(target = self.stop).start()
        run_code.pack(side=tk.RIGHT,  padx=5,pady=5, fill=tk.X)
        
        self.out_bar  =  tk.Frame(self.filter_bar,  width=15,  height=1)
        self.out_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=10,  expand= False)        
        lbl = tk.Label(self.out_bar, text="Running Status:", bg='grey')
        lbl.pack(fill=tk.X)
        lbl.config(font=("Font", 10))

        self.out_bar  =  tk.Frame(self.filter_bar,  width=15,  height=1)
        self.out_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=10,  expand= False)      
        self.scrollbar = tk.Scrollbar(self.out_bar)
        self.scrollbar.pack( side = tk.RIGHT, fill = tk.Y )
        self.SHBar = tk.Scrollbar(self.out_bar, orient = tk.HORIZONTAL)
        self.SHBar.pack (side = tk.BOTTOM, fill = "x")
        self.mylist = tk.Listbox(self.out_bar, width=35, yscrollcommand = self.scrollbar.set, xscrollcommand = self.SHBar.set )
        self.mylist.pack( side = tk.RIGHT, fill = tk.BOTH ,  expand=True)
        self.scrollbar.config( command = self.mylist.yview )
        self.SHBar.config(command = self.mylist.xview)
        
        res = tk.Label(self.res_win, text="Results", bg ='light grey')
        res.pack(fill=tk.X)
        res.config(font=("Font", 12))
        
        self.dis_bar = tk.Frame(self.res_win,  width=180,  height=150, bg='light grey')
        self.dis_bar.pack(side='left',  fill='both',  padx=5,  pady=5,  expand=True)
        
        self.out_bar  =  tk.Frame(self.dis_bar,  width=15,  height=1)
        self.out_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= False)   
        lbl = tk.Label(self.out_bar, text="Statistics", bg='grey')
        lbl.pack(fill=tk.X)
        lbl.config(font=("Font", 10))
         
        self.out_bar  =  tk.Frame(self.dis_bar,  width=15,  height=1)
        self.out_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= False) 
        tk.Label(self.out_bar, text="number of nodes:").pack(side=tk.LEFT, fill=tk.X)
        self.nodeshow = tk.Label(self.out_bar, text = "N/A")
        self.nodeshow.pack(side=tk.LEFT, fill=tk.X)
        
        self.out_bar  =  tk.Frame(self.dis_bar,  width=15,  height=1)
        self.out_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= False) 
        tk.Label(self.out_bar, text="number of edges:").pack(side=tk.LEFT, fill=tk.X)
        self.edgeshow = tk.Label(self.out_bar, text = "N/A")
        self.edgeshow.pack(side=tk.LEFT, fill=tk.X)

        self.out_bar  =  tk.Frame(self.dis_bar,  width=15,  height=1)
        self.out_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= False) 
        tk.Label(self.out_bar, text="number of cores:").pack(side=tk.LEFT, fill=tk.X)
        self.coreshow = tk.Label(self.out_bar, text = "N/A")
        self.coreshow.pack(side=tk.LEFT, fill=tk.X)

        self.out_bar  =  tk.Frame(self.dis_bar,  width=15,  height=1)
        self.out_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=10,  expand= False)   
        lbl = tk.Label(self.out_bar, text="Figures", bg='grey')
        lbl.pack(fill=tk.X)
        lbl.config(font=("Font", 10))
        
        self.out_bar  =  tk.Frame(self.dis_bar,  width=15,  height=1)
        self.out_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=0,  expand= False)         
        self.options = ['Threshold curve', 'Degree distribution', 'Communities','Cores']
        self.choice = tk.StringVar(value = 'Select a figure')

        dropdown = ttk.OptionMenu( self.out_bar , self.choice, 'Select a figure',*self.options, command = self.showattr) #lambda: threading.Thread(target = self.showattr).start()
        dropdown.pack(side=tk.LEFT, fill=tk.X, padx=5, pady=5, expand = True)
        
        self.fig_bar  =  tk.Frame(self.dis_bar,  width=15,  height=1)
        self.fig_bar.pack(side='top',  fill=tk.X,  padx=5,  pady=5,  expand= False)         

        self.show_fig = tk.Label(self.fig_bar, image = '')
        self.show_fig.pack(side = tk.TOP)
        self.show_fig.image = ''
        self.show_fig.configure(image = '')
        
        #self.event_1.start()
    def stop(self):
        self.stopflag = True
    def stopset(self):
        self.mylist.insert(tk.END, '--Process is terminated.')
        self.stopflag = False    
    
    def read_file(self):
        self.filename = tkinter.filedialog.askopenfilename()
        if self.filename:
            self.data.delete(0,tk.END)
            self.data.insert(0, self.filename)
        print('Read data file:', self.filename)
        
    def open_edge(self):
        self.opend_edge_file = tkinter.filedialog.askopenfilename()
        if self.opend_edge_file:
            self.readedge.delete(0,tk.END)
            self.readedge.insert(0, self.opend_edge_file)
        print('Read edge list:', self.opend_edge_file)        
        
    def clicked(self):

        self.stopflag = False  
        self.mylist.insert(tk.END, '*********************************')
        self.mylist.insert(tk.END, 'Network Generation:')
        
        if not os.path.exists(self.data.get()):
            self.mylist.insert(tk.END, '--Data not found, check the input.')
            return
        
        self.mylist.insert(tk.END, '--Reading data:')
        self.mylist.insert(tk.END, '    '+ self.data.get())
        gene = dataload.GeneCoExp( filepath = self.data.get())
        
        self.mylist.insert(tk.END, '--Processing data...')
            
        gene.read_data( zero_removed = self.zero_.get(), 
                       zscored = self.zscored.get(), 
                       rescaled = self.rescaled.get(), 
                       dropna = int(self.dropna.get()),
                       savefile = self.save_data.get())  
        if gene.M <= 4:
            self.mylist.insert(tk.END, '    check the input data.')
            
        if self.save_data.get():
            self.mylist.insert(tk.END, '--The normalzied dataset is saved as:')
            self.mylist.insert(tk.END, '    ' + gene.savepath + ' (z-scored).csv')
        
        self.mylist.insert(tk.END, '--Calculating PCCs...')
        gene.Calculate_PCC(pcc_path = self.pcc.get() )          
        if self.pcc.get():
            self.mylist.insert(tk.END, '--The PCC matrix is saved as:')
            self.mylist.insert(tk.END, '    ' + gene.savepath + ' PCC_matrix.csv')
        if self.stopflag:
            self.stopset()
            return  
        
        self.mylist.insert(tk.END, '--Determining the number of paired elements for every node pair...')
        gene.paired_elements(paired_elements_path = self.paired.get())
        if self.paired.get():
            self.mylist.insert(tk.END, '--The paired elements matrix is saved as: paired_elements_matrix.csv')        
        if self.stopflag:
            self.stopset()
            return  
              
        self.mylist.insert(tk.END, '--Determining the threshold of each bin...')
        gene.calculate_threshold(bin_size = int(self.bin_size.get()), cutoff = float(self.cutoff.get()))     
        if self.stopflag:
            self.stopset()
            return          
        
        self.mylist.insert(tk.END, '--Fitting Curve & thresholding node pairs...')
        gene.curve_fitting(alpha = float(self.alpha.get()),
                           eta = float(self.eta.get()),
                           lam = float(self.lam.get()),
                           beta = float(self.beta.get()), #assign initial values to the 4 parameters
                           para_in= self.curve_params.get(), 
                           edgelist = self.edgelist.get(), 
                           thresholdcurve = self.thres_curve.get())
        self.alpha.delete(0, tk.END)
        self.alpha.insert(0, round(gene.popt[0],2))
        self.eta.delete(0, tk.END)
        self.eta.insert(0, round(gene.popt[1],2))
        self.lam.delete(0, tk.END)
        self.lam.insert(0, round(gene.popt[2],2))
        self.beta.delete(0, tk.END)
        self.beta.insert(0, round(gene.popt[3],2))

        if self.edgelist.get():
            self.mylist.insert(tk.END, '--Edge list is saved as:')
            self.mylist.insert(tk.END, '    ' + gene.savepath + ' edgelist.csv')
        if self.thres_curve.get():
            self.mylist.insert(tk.END, '--The threshold curve is saved as: threshold_curve.png')
            
            img = Pil_image.open("threshold_curve.png")
            img = img.resize((480, 360)) #, resample =Pil_image.Resampling.LANCZOS
            img = Pil_imageTk.PhotoImage(img)   
            self.show_fig.image = img
            self.show_fig.configure(image = img)  
            
        self.mylist.insert(tk.END, '--The R-suqared of the fitted curve is:')
        self.mylist.insert(tk.END, '    R2 = ' + str(round(gene.r2,3)))
        
        self.nodeshow.config(text = gene.n_nodes)
        self.edgeshow.config(text = gene.n_edges) 
        
        self.mylist.insert(tk.END, '--Completed!')
        self.stopflag = False
        
    def NetAttr(self):
        self.stopflag = False  
        
        self.mylist.insert(tk.END, '*********************************')
        self.mylist.insert(tk.END, 'Network Analysis:')
        
        if not os.path.exists(self.readedge.get()):
            self.mylist.insert(tk.END, '--Edge list not found...')
            return
        
        self.mylist.insert(tk.END, '--Reading edge list:')
        self.mylist.insert(tk.END, '    '+ self.readedge.get())
        netw = na.networkanalysis(self.readedge.get(), weighted = self.weighted.get())
        
        self.mylist.insert(tk.END, '--Pre-processing the network...')
        netw.find_pos(lcc= self.lcc.get(), core = self.core.get(), louvain_community = self.community.get())   
        if self.lcc.get():
            self.mylist.insert(tk.END, '--The largest connected component is saved as:')
            self.mylist.insert(tk.END, '    ' + netw.savepath + ' (LCC).csv')
            
        self.nodeshow.config(text = len(netw.G.nodes()))
        self.edgeshow.config(text = len(netw.G.edges())) 

        if self.stopflag:
            self.stopset()
            return   
        
        if self.community.get():
            self.mylist.insert(tk.END, '--Detecting communities...')
            netw.netcommu()
            if self.stopflag:
                self.stopset()
                return
            
        if self.core.get():
            self.mylist.insert(tk.END, '--Determining cores...')
            netw.netcore()
            self.coreshow.config(text = len(netw.corenode))  
            if self.stopflag:
                self.stopset()
                return         
            
        if self.degree.get() or self.eigenvector.get():
            self.mylist.insert(tk.END, '--Determining degree and (or) eigenvector centrality...')
            netw.netdeg(degree = self.degree.get(), eigenvector = self.eigenvector.get())
            if self.stopflag:
                self.stopset()
                return 
            
        if self.closeness.get():
            self.mylist.insert(tk.END, '--Determining closeness centrality...')
            netw.netclose(closeness = self.closeness.get())            
            if self.stopflag:
                self.stopset()
                return 

        if self.betweenness.get():
            self.mylist.insert(tk.END, '--Determining betweenness centrality...')
            netw.netbet(betweenness = self.betweenness.get())            
            if self.stopflag:
                self.stopset()
                return 
            
        if self.community.get() or self.degree.get() or self.core.get() or self.eigenvector.get() or self.betweenness.get() or self.closeness.get():
            netw.saveattr()
            self.mylist.insert(tk.END, '--Property Table is saved as:')
            self.mylist.insert(tk.END, netw.savepath+ ' property list.csv')
            
        self.mylist.insert(tk.END, '--Completed!')
        self.stopflag = False
        
    def showattr(self, *args):
        imgpath = 'check figure path'
        if self.choice.get() == self.options[0]:
            imgpath = "threshold_curve.png"
        elif self.choice.get() == self.options[1]:
            imgpath = "degree_distribution.png"
        elif self.choice.get() == self.options[2]:
            imgpath = "community.png"
        elif self.choice.get() == self.options[3]:
            imgpath = "core.png"

        if os.path.exists(imgpath):
            img = Pil_image.open(imgpath)
            img = img.resize((480, 360)) #, resample =Pil_image.Resampling.LANCZOS
            img = Pil_imageTk.PhotoImage(img)   
        else:
            img = ''
        self.show_fig.image = img
        self.show_fig.configure(image = img)   
        self.stopflag = False
        # print(args, self.choice.get())
        # if self.choice.get() == self.options[0]:
        #     print(self.choice.get())
        #     self.thr_img.config( text='1') #image = tk.PhotoImage(file="threshold_curve.png")
        # self.thr_img.config( image = tk.PhotoImage(file="threshold_curve.png")) #
if __name__ == "__main__":
    window = tk.Tk()
    myapp = App(window)
    window.mainloop()
    
    

    

    
    


