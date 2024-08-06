#!/usr/bin/python3
import os
import os.path
import tkinter as tk
from tkinter import filedialog
from tkinter import StringVar
from tkinter import messagebox

probefile=None
bedfile=None



def get_probefile():
    global probefile
    probefile_userprovided=tk.filedialog.askopenfilename(defaultextension=".txt",filetypes=[("Probe Files",".txt"),("All Files", "*.*")])
    if os.path.isfile(probefile_userprovided):
        probefile=probefile_userprovided
        #update text field
        iprobe_text.config(state="normal")
        iprobe_text.delete('1.0','end')
        iprobe_text.insert("end",probefile)
        iprobe_text.config(state="disabled")

    
    

def get_bedfile():
    global bedfile
    bedfile_userprovided=tk.filedialog.askopenfilename(defaultextension=".bed",filetypes=[("Bed File",".bed"),("All Files", "*.*")])
    if os.path.isfile(bedfile_userprovided):
        bedfile=bedfile_userprovided
        #update text field
        ibed_text.config(state='normal')
        ibed_text.delete('1.0','end')
        ibed_text.insert('end',bedfile)
        ibed_text.config(state='disabled')


def apply_probefiltering():
    if probefile==None:
        tk.messagebox.showinfo('Error','Please select probes file first')
        return 0
    if bedfile==None:
        tk.messagebox.showinfo('Error','Please select bed file first')
        return 0
    if not os.path.isfile(probefile):
        tk.messagebox.showinfo("Error",f"Error accessing probefile {probefile}\nPlease reimport")
        return 0
    if not os.path.isfile(bedfile):
        tk.messagebox.showinfo("Error",f"Error accessing bedfile {bedfile}\nPlease reimport")
        return 0
    #All parameters seem to be set. Apply filter
    filterprobesbybed(bedfile,probefile)

def bed2dict(bedfile):
    """
        This function reads a 3 column bed file and returns it as a bed object!
    """
    #get beds
    fi=open(bedfile,'r')
    beds={} # {chr1:[[start,end],[],[]],chr.:[]}
    for line in fi:
        if line.strip()=="":continue
        if line.startswith("track"):continue
        if line.startswith("browser"):continue
        sline=line.strip().split('\t')
        cchr=sline[0]
        start=int(sline[1])
        end=int(sline[2])
        if not cchr in beds:beds[cchr]=[]
        beds[cchr].append((start,end))
    fi.close()
    #sort the bed object
    for chrom in beds:
        beds[chrom].sort()
    return(beds)

def getcoordinates(sline):
    coords=[]
    try:
        pend=int(sline[-1])
        pstart=int(sline[-2])
        pchrom=sline[-3]
        coords.append([pchrom,pstart,pend])
    except:
        if "|" in sline[-1]:
            pcords=sline[-1].split("|") # multiple coordinates collapses with pipe for DNA
        elif "," in sline[-1]:
            pcords=sline[-1].split(",") # multiple coordinates collapses with comman for RNA
        else:
            pcords=[sline[-1]]
        for pcord in pcords:
            coord=pcord.replace('-',':').split(':')
            pchr=coord[0]
            pstart=int(coord[1])-1 #fix for SD6 2 SD8
            pend=int(coord[2])
            coords.append([pchr,pstart,pend])
    return (coords)

def filterprobesbybed(bed,probes):
    """
        This script filters the probes out that overlaps with a region from probes.txt files for Sureselect.
    """
    inputbedfile=bed
    if not type(bed) is dict:
        beds=bed2dict(bed) # {chr1:[[start,end],[],[]],chr.:[]}
    else:
        beds=bed
    #read probes
    intarget=0
    outtarget=0
    fi=open(probes,'r')
    fo=open(probes[:-4]+'_filtered_outBed.txt','w')
    fo2=open(probes[:-4]+'_filtered_inBed.txt','w')
    for line in fi:
        if line.strip()=="":continue
        if line.startswith('TargetID\tProbeID'):
            fo.write(line)
            fo2.write(line)
            continue
        sline=line.strip().split('\t')
        istarget=False
        coords=getcoordinates(sline)
        for coord in coords:
            pchr,pstart,pend=coord
            if pchr in beds:
                for bed in beds[pchr]:
                    if bed [0] > pend:break
                    dist=max(bed[0],pstart)-min(bed[1],pend)
                    if dist<=0:
                        istarget=True
                        break
        if istarget:
            fo2.write(line)
            intarget+=1
        else:
            fo.write(line)
            outtarget+=1
    fo.close()
    fo2.close()
    fi.close()
    print("File: {}\ton-target: {}\toff-target: {}\ttotal: {}".format(probes,intarget,outtarget,intarget+outtarget))
    tk.messagebox.showinfo("Completed!","Probe File: {}\nBedFile:{}\non-target: {}\noff-target: {}\ntotal: {}".format(probes,inputbedfile,intarget,outtarget,intarget+outtarget))
#UI 

root=tk.Tk()
root.title("Filter probes by bed")

topframe=tk.Frame(root,height=30, pady=5,padx=5)
topframe.pack(expand="YES",fill='both')


probeframe=tk.Frame(root,height=30,pady=5,padx=5)
probeframe.pack(expand="YES",fill='both')

iprobe_button=tk.Button(probeframe,text="Import probes file",command=get_probefile)
iprobe_button.pack(side="left")

iprobe_text=tk.Text(probeframe,height=1,width=100)
iprobe_text.pack(side="left",expand=True)
iprobe_text.insert("end","None")
iprobe_text.config(state='disabled')

bedframe=tk.Frame(root,height=30,pady=5,padx=5)
bedframe.pack()

ibed_button=tk.Button(bedframe,text="Import bed file",command=get_bedfile)
ibed_button.pack(side='left')
ibed_text=tk.Text(bedframe,height=1,width=100)
ibed_text.pack(side='left')
ibed_text.insert("end","None")
ibed_text.config(state='disabled')

apply_filterbutton=tk.Button(root,text="Apply Filter",command=apply_probefiltering)
apply_filterbutton.pack()

root.mainloop()
