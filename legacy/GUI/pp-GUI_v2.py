# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 14:27:21 2020

@author: OYCX
"""
import h5py
import time
from threading import Thread
from tkinter import *
from tkinter import messagebox
from tkinter import ttk
from tkinter import filedialog
from PIL import Image, ImageTk
import matplotlib.pyplot as plt
import numpy as np
import os, sys


class preprocess(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.master.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.master.title("pre-processing script")
        self.frame_master = Frame(master, bd=1)
        self.frame_master.pack()

        # Menu bar
        self.menu_filemenu = Menu(self.master, tearoff=0)
        self.master.config(menu=self.menu_filemenu)
        self.menu_filemenu.add_command(label='Save Path', command=lambda: self.save_path())

        # Path & name variables
        self.var_file_path = StringVar(self.master, value="Z:/Lab Data/Rosetta CRUK GC" )
        self.var_savepath = self.var_file_path
        self.var_filename = StringVar(self.master)
        self.var_filename.set(self.var_file_path.get() + "/" + time.strftime("%Y_%m_%d") + ".raw")
        self.var_author = StringVar(self.master)
        self.var_note = StringVar(self.master)
        self.var_specname = StringVar(self.master, value="DESI")


        # File/exp labels
        self.frame_filepath =ttk.Frame(self.master)
        self.frame_filepath.pack(side=TOP, fill=X)
        Label(self.frame_filepath, text="Current data file: ", padx=2, pady=2).pack(side=LEFT)
        Label(self.frame_filepath, textvariable=self.var_filename, bg='gray').pack(side=LEFT)
        Label(self.frame_filepath, textvariable=self.var_specname, bg='gray').pack(side=RIGHT)
        Button(self.frame_filepath,text='Browse',command=lambda: self.choose_path()).pack(side=LEFT)
        Label(self.frame_filepath, text = 'Current exp type: ', padx=2, pady=2).pack(side=RIGHT)

        # Tabs (PAOLO for old script and ALT for new alternative approach, coming soon)
        global filename, specname
        filename = self.var_filename
        specname = self.var_specname
        # tab configuration
        self.notebook = ttk.Notebook(master)
        self.notebook_frame1 = ttk.Frame(self.notebook)
        self.notebook_frame2 = ttk.Frame(self.notebook)
        self.notebook_frame3 = ttk.Frame(self.notebook)
        self.notebook.add(self.notebook_frame1, text="PAOLO")
        self.notebook.add(self.notebook_frame2, text="ALT")
        self.notebook.pack(fill=BOTH)
        
        self.start_switch1 = IntVar()
        self.start_switch2 = IntVar()
        
        ## buttons to run script
        
        # peak picking
        self.labelframe1 = LabelFrame(master, text="Peak picking")
        self.labelframe1.pack(fill="both", expand="yes")
        b1=Button(self.labelframe1, command=lambda: self.threading(self.peak_pick),
    text="start peak picking",
    width=30,
    height=10,
    font=12,
    bg="red",
    fg="blue"
    )
        b1.pack()
        
        # ROI selection buttons, with method selection
        self.labelframe2 = LabelFrame(master, text="ROI selection (choose a method)")
        self.labelframe2.pack(fill="both", expand="yes")
        b2=Button(self.labelframe2, command=lambda: self.threading(self.roi_select),
    text="start ROI selection",
    width=30,
    height=10,
    font=12,
    bg="red",
    fg="blue")
        b2.pack()
        
        self.method = IntVar()

        Radiobutton(root, 
                      text="Manual",
                      padx = 20,
                      justify=LEFT,
                      font=10,
                      variable=self.method, 
                      value=1).pack()
        Radiobutton(root, 
                      text="SVM",
                      padx = 20,
                      justify=LEFT,
                      font=10,
                      variable=self.method, 
                      value=2).pack()
        Radiobutton(root, 
                      text="skip",
                      padx = 20,
                      justify=LEFT,
                      font=10,
                      variable=self.method, 
                      value=3).pack()
        
        self.labelframe3 = LabelFrame(master, text="Peak re-alignment")
        self.labelframe3.pack(fill="both", expand="yes")
        b3=Button(self.labelframe3, command=lambda: self.threading(self.diagnose_match),
    text="start peak matching",
    width=30,
    height=10,
    font=12,
    bg="red",
    fg="blue")
        b3.pack()
        
    #             self.labelframe1 = LabelFrame(master, text="Peak picking")
    #     self.labelframe1.pack(fill="both", expand="yes")
    #     b1=Button(self.labelframe1, text="start peak picking",
    # width=25,
    # height=5,
    # bg="red",
    # fg="blue")
    #     b1.pack()
        
    #             self.labelframe1 = LabelFrame(master, text="Peak picking")
    #     self.labelframe1.pack(fill="both", expand="yes")
    #     b1=Button(self.labelframe1, text="start peak picking",
    # width=25,
    # height=5,
    # bg="red",
    # fg="blue")
    #     b1.pack()
        
    def on_closing(self):
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.master.destroy()

    def choose_path(self, *args):
        self.var_filename.set(filedialog.askdirectory())
        self.var_file_path.set(os.path.dirname(os.path.abspath(self.var_filename.get())))
        # self.var_filename.set(self.var_file_path.get() + "/" + time.strftime("%Y_%m_%d") + ".hdf5")
        
    def choose_savepath(self, *args):
        self.var_savepath.set(filedialog.askdirectory())

    def set_filename(self, *args):
        self.var_filename.set(self.var_file_path.get() + "/" + time.strftime("%Y_%m_%d") + ".hdf5")
        print(self.var_filename.get())
        h5 = h5py.File(self.var_filename.get(), "a")
        h5brillouin = h5.create_group("Brillouin data")
        h5bf = h5.create_group("Bright field")
        h5brillouin.create_group(self.var_specname.get())
        h5bf.create_group(self.var_specname.get())
        timestamp = time.strftime("%Y_%m_%d__%H_%M")
        h5.attrs['time'] = timestamp
        h5.attrs['Author'] = self.var_author.get()
        h5.close()

    # popup window to set saving options
    def save_path(self, *args):
        self.top = Toplevel()
        self.top.title('Where to save?')
        label_csp = Label(self.top, text="Path: ", padx=2, pady=2)
        entry_csp = Entry(self.top, textvariable=self.var_savepath, width=80)
        button_cp = Button(self.top, text="Choose path", command=lambda: self.choose_savepath(), padx=2, pady=2)
        label_expname = Label(self.top, text = "Experiment name: ", padx=2, pady=2)
        entry_expname = Entry(self.top, textvariable=self.var_specname, width=10)
        label_an = Label(self.top, text="Author: ", padx=2, pady=2)
        entry_author = Entry(self.top, textvariable=self.var_author, width=10)
        button_setfp = Button(self.top, text="SET", command=self.top.destroy, padx=5, pady=5)
        label_note = Label(self.top, text="Notes:", padx=2, pady=2)
        entry_note = Entry(self.top, textvariable=self.var_note, width=10)

        label_csp.grid(row=0, column=0,sticky=E+W,)
        entry_csp.grid(row=0, column=1,sticky=E+W,)
        button_cp.grid(row=0, column=2,sticky=E+W,)
        label_expname.grid(row=1, column=0,sticky=E+W,)
        entry_expname.grid(row=1, column=1,sticky=E+W,)
        label_an.grid(row=2, column=0,sticky=E+W,)
        entry_author.grid(row=2, column=1,sticky=E+W,)
        button_setfp.grid(row=2, column=2, rowspan=3,sticky=E+W,)
        label_note.grid(row=3,column=0, stick=E+W)
        entry_note.grid(row=3, column=1, sticky=E+W)
        self.top.attributes('-topmost','TRUE')
        
    def waitbar (self, script):
        self.window1 = Toplevel(width=100,height=60)
        self.window1.title('running script')
        x = root.winfo_x()
        y = root.winfo_y()
        self.window1.geometry("+%d+%d" % (x + 300, y + 200))
        self.progress = ttk.Progressbar(self.window1, orient = HORIZONTAL, 
            length = 300, mode = 'indeterminate') 
        self.progress.pack()
        self.progress.start(interval=50)
        
        script()
        
        self.progress.stop()
        # # disable buttons during
        # for btn in buttons:
        #     btn['state'] = 'disabled'
        w = Label(self.window1, text='Done!', width=30, height=10,font=12)
        w.pack()
        b = Button(self.window1, text="OK", command=self.window1.destroy, width=10)
        b.pack()
        self.window1.attributes('-topmost','TRUE')

        # # Enable all buttons
        # for btn in buttons:
        #     btn['state'] = 'normal'
    
    def threading (self, script):
        Thread(target=self.waitbar, kwargs={'script':script}).start()
        
    # function to run peak_detect.py
    def peak_pick (self):
        detect = open("./peak_detect.py").read()
        input_path=os.path.dirname(os.path.abspath(self.var_filename.get()))
        output_path=self.var_savepath.get()
        sys.argv = ["peak_detect.py", input_path, output_path]     
        exec(detect)
        
    # function to run select_roi.py
    def roi_select (self):
        select = open("./select_roi.py").read()
        input_path=self.var_savepath.get()
        output_path=self.var_savepath.get()
        if self.method.get() == 1:
            sys.argv = ["select_roi.py", input_path]   
        if self.method.get() == 2:
            sys.argv = ["select_roi.py", '--supervised', input_path]
        if self.method.get() == 3:
            sys.argv = ["select_roi.py", '--skip', input_path]    
        exec(select)
             
    # function to run match_peaks.py
    def peak_match (self):
        align = open("./match_peaks.py").read()
        input_path=self.var_savepath.get()
        if self.reference.get()==1:
            if self.method.get() == 2:
                sys.argv = ["match_peaks.py", input_path, '--roi-type', 'auto']
            else:
                sys.argv = ["match_peaks.py", input_path]
        if self.reference.get()==2:
            if self.method.get() == 2:
                sys.argv = ["match_peaks.py", input_path, '--roi-type', 'auto', '--first-ref-mz',self.custom_masses[0],
                           '--sec-ref-mz-1',self.custom_masses[1],'--sec-ref-mz-2',self.custom_masses[2]]
            else:
                sys.argv = ["match_peaks.py", input_path, '--first-ref-mz',self.custom_masses[0],
                           '--sec-ref-mz-1',self.custom_masses[1],'--sec-ref-mz-2',self.custom_masses[2]]
        exec(align)
        
    # function to perform drift diagnostics and then re-alignment of peaks
    def diagnose_match (self):
        if self.start_switch1.get() == 1:
            self.peak_match()
        else:
            if self.start_switch2.get() == 1:
                diagnose = open("./shift_diagnostics.py").read()
                input_path=self.var_savepath.get()
                output_path=self.var_savepath.get()
                if self.reference.get()==1:
                    sys.argv = ["shift_diagnostics.py", input_path]
                if self.reference.get()==2:
                    if self.refs_entry1.get()=='pos':
                        sys.argv = ["shift_diagnostics.py", input_path, "--ref-mz-pos", 
                                    self.custom_masses[0], self.custom_masses[1], self.custom_masses[2]]
                    if self.refs_entry1.get()=='neg':
                        sys.argv = ["shift_diagnostics.py", input_path, "--ref-mz-neg", 
                                    self.custom_masses[0], self.custom_masses[1], self.custom_masses[2]]
                exec(diagnose)
                param1 = TOL[sel_ind]
                # param1 = 10
                param2 = 50
                param3 = 50
                param4 = param1
                self.params = np.array([param1,param2,param3,param4])
                self.show_diagnostics(self.params)
            else:
                self.diagnose_options()

    # functions to set shift_diagnostics options (to set custom references)
    def diagnose_options(self):
        self.reference = IntVar()
        self.progress.stop()
        self.window1.destroy()
        self.window2 = Toplevel(height = 100, width = 300)
        self.window2.title('select options')
        x = root.winfo_x()
        y = root.winfo_y()
        self.window2.geometry("+%d+%d" % (x + 300, y + 200))
        self.b1 = Button(self.window2, text="default", command=lambda:[self.reference.set(1),
                                                                       self.start_switch2.set(1),
                                                                       self.window2.destroy(),
                                                                       self.threading(self.diagnose_match)])
        
        self.b2 = Button(self.window2, text="custom references", command=lambda:[self.reference.set(2),
                                                                                 self.window2.destroy(),
                                                                                 self.custom_refs()])
        self.b1.place( x=125, y = 25)
        self.b2.place( x=100, y = 70)

        self.window2.attributes('-topmost','TRUE')
        
    def custom_refs (self):
        self.refs_window = Toplevel()
        self.refs_window.title('define references masses to use:')
        x = root.winfo_x()
        y = root.winfo_y()
        self.refs_window.geometry("+%d+%d" % (x + 300, y + 200))
        self.refs_entry1=StringVar()
        self.refs_entry2=StringVar()
        label_refs1 = Label(self.refs_window, text="polarity (pos or neg): ", padx=2, pady=2)
        entry_refs1 = Entry(self.refs_window, textvariable=self.refs_entry1, width=50)
        label_refs2 = Label(self.refs_window, text="masses (seprate with ,): ", padx=2, pady=2)
        entry_refs2 = Entry(self.refs_window, textvariable=self.refs_entry2, width=50)
        
        b = Button(self.refs_window, text="Set & Start", command=lambda:[self.set_refs(),
                                                            self.start_switch2.set(1),
                                                            self.refs_window.destroy(),
                                                            self.threading(self.diagnose_match),
                                                            print(self.custom_masses)], width=10)
        label_refs1.grid(row=0, column=0,sticky=E+W,)
        entry_refs1.grid(row=0, column=1,sticky=E+W,)
        label_refs2.grid(row=1, column=0,sticky=E+W,)
        entry_refs2.grid(row=1, column=1,sticky=E+W,)
        b.grid(row=2, column=1, rowspan=2,sticky=E+W,)
        self.refs_window.attributes('-topmost','TRUE')
        
    def set_refs(self):
        lst = self.refs_entry2.get().split(",")
        self.custom_masses = lst           
        
    # functions to show diagnostic results and set alignment parameters    
    def show_diagnostics(self, params):
        self.progress.stop()
        self.window1.destroy()
        self.img_window = Toplevel()
        self.img_window.title('showing images')
        x = root.winfo_x()
        y = root.winfo_y()
        self.img_window.geometry("+%d+%d" % (x + 300, y + 200))
        self.img_window.attributes('-topmost','TRUE')
        self.display = Canvas(self.img_window,width = 1000, height = 800)
        self.display.pack(fill="both", expand=True)
        
        input_dir = self.var_savepath.get()
        filename = os.path.basename(self.var_filename.get())
        filename = os.path.splitext(filename)[0]
        file = os.path.join(input_dir, filename, "diag", "perc_pix_ref_mz.png") 
        while not os.path.exists(file):
            time.sleep(5)
        # self.vsb = Scrollbar(window, orient="vertical", command=self.display.yview)
        # self.vsb.grid(row=0, column=1, sticky="ns")
        # self.display.configure(yscrollcommand= vsb.set)
        # canvas.grid(row=0, column=0, sticky= "news")
        # self.display = Canvas(self, bd=0, highlightthickness=0)
        # self.image1 = Image.open('D:/BOX/Box Sync/RA/codes/DESI/pre-processing/GUI/test data/2019_09_01_porklivertest_100umpixel_neg Analyte 1/diag/255.233_ref_im.png')
        # self.image2 = Image.open('D:/BOX/Box Sync/RA/codes/DESI/pre-processing/GUI/test data/2019_09_01_porklivertest_100umpixel_neg Analyte 1/diag/303.233_ref_im.png')
        # self.image3 = Image.open('D:/BOX/Box Sync/RA/codes/DESI/pre-processing/GUI/test data/2019_09_01_porklivertest_100umpixel_neg Analyte 1/diag/885.5499_ref_im.png')
    
        self.image4 = Image.open(file)
        # self.resized1 = self.image1.resize((800, 500),Image.ANTIALIAS)
        # self.resized2 = self.image2.resize((800, 500),Image.ANTIALIAS)
        # self.resized3 = self.image3.resize((800, 500),Image.ANTIALIAS)
        self.resized4 = self.image4.resize((800, 600),Image.ANTIALIAS)
        # self.render1 = ImageTk.PhotoImage(self.resized1)
        # self.render2 = ImageTk.PhotoImage(self.resized2)
        # self.render3 = ImageTk.PhotoImage(self.resized3)
        self.render4 = ImageTk.PhotoImage(self.resized4)
        
        # self.display.image_names = self.render1
        # self.display.image_names = self.render2
        # self.display.image_names = self.render3
        self.display.image_names = self.render4
        # img = Label(self.img_window, image=render)
        # self.display.create_image((400, 250), image=self.render1)
        # self.display.create_image((1200, 250), image=self.render2)
        # self.display.create_image((400, 750), image=self.render3)
        self.display.create_image((500, 300), image=self.render4)
        
        self.message = ('The recommended match_params are:\n\n' + str(params) )
        self.display.create_text((500,650),text=self.message, font=15)
        
        self.b1 = Button(self.display, text="Accept", command = lambda:[self.generate_params(self.params), 
                                                                        self.start_switch1.set(1), 
                                                                        self.img_window.destroy(), 
                                                                        self.threading(self.peak_match)], width=10)
        self.b2 = Button(self.display, text="Change", command =self.change_params, width=10)
        # b1.configure(width = 10, activebackground = "#33B5E5", relief = FLAT)
        self.display.create_window((400,700), window=self.b1, anchor=CENTER)
        self.display.create_window((500,700), window=self.b2, anchor=CENTER)
        

    def set_params (self):
        self.params[0] = int(self.input1.get())
        self.params[1] = int(self.input2.get())
        self.params[2] = int(self.input3.get())
        self.params[3] = int(self.input4.get())
        
    
    def generate_params (self, params):
        output_dir = self.var_savepath.get()
        filename = os.path.basename(self.var_filename.get())
        output_dir = os.path.join(output_dir, filename[:-4])
        file = os.path.join(output_dir, 'match_params.txt') 
        np.savetxt(file, params, fmt='%-4d', delimiter='\n')
        
    def change_params (self):
        self.params_window = Toplevel()
        self.params_window.title('set peak_match parameters')
        x = root.winfo_x()
        y = root.winfo_y()
        self.params_window.geometry("+%d+%d" % (x + 300, y + 200))
        self.input1=StringVar() 
        self.input2=StringVar() 
        self.input3=StringVar() 
        self.input4=StringVar() 
        label_param1 = Label(self.params_window, text="Parameter 1: ", padx=2, pady=2)
        entry_param1 = Entry(self.params_window, textvariable=self.input1, width=10)
        label_param2 = Label(self.params_window, text="Parameter 2: ", padx=2, pady=2)
        entry_param2 = Entry(self.params_window, textvariable=self.input2, width=10)
        label_param3 = Label(self.params_window, text="Parameter 3: ", padx=2, pady=2)
        entry_param3 = Entry(self.params_window, textvariable=self.input3, width=10)
        label_param4 = Label(self.params_window, text="Parameter 4: ", padx=2, pady=2)
        entry_param4 = Entry(self.params_window, textvariable=self.input4, width=10)
        
        
        b = Button(self.params_window, text="Set & Start", command=lambda:[self.set_params(),
                                                            self.generate_params(self.params), 
                                                            self.start_switch1.set(1),
                                                            self.params_window.destroy(), 
                                                            self.img_window.destroy(), 
                                                            self.threading(self.peak_match)], width=10)
        
        label_param1.grid(row=0, column=0,sticky=E+W,)
        entry_param1.grid(row=0, column=1,sticky=E+W,)
        label_param2.grid(row=1, column=0,sticky=E+W,)
        entry_param2.grid(row=1, column=1,sticky=E+W,)
        label_param3.grid(row=2, column=0,sticky=E+W,)
        entry_param3.grid(row=2, column=1,sticky=E+W,)
        label_param4.grid(row=3, column=0,sticky=E+W,)
        entry_param4.grid(row=3, column=1,sticky=E+W,)
        b.grid(row=4, column=1, rowspan=3,sticky=E+W,)
        self.params_window.attributes('-topmost','TRUE')
        
root = Tk()



# root.geometry("1200x850")

# centre window
ws = root.winfo_screenwidth()
hs = root.winfo_screenheight() 
root.geometry('%dx%d+%d+%d' % (1200, 850, (ws/2) - (1200/2), (hs/2) - (850/2)))


# creation of an instance
app = preprocess(master=root)
# mainloop
root.mainloop()