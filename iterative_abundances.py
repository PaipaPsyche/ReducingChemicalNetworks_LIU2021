import numpy as np
import itertools as itr
from decimal import Decimal
import os
from shutil import copyfile
import subprocess
#def gen_random_conditions(templates_path,):
    
def eformat(f, prec, exp_digits):
    s = "%.*E"%(prec, f)
    mantissa, exp = s.split('E')
    # add 1 to digits as 1 is taken by sign +/-
    return "%sE%+0*d"%(mantissa, exp_digits+1, int(exp))


    
    
def change_speci(path,newvals,new_path):
    lines = []
    with open(path,"r") as f:
        while(True):
            line_ = f.readline()
            if not line_:
                break
            if(line_.startswith("#")):
                lines.append(line_)
                continue
            line = line_.split()
            if(len(line)==1):
                lines.append(line_)
                continue
            
            if line[1] in newvals.keys():
                val = newvals[line[1]]
                #print(line[1],val)
                newline =""
                newline = newline+line_[:44]
                newline = newline+'%.3E' % Decimal(val)
                newline = newline+line_[53:]
                lines.append(newline)
            else:
                lines.append(line_)
    with open(new_path,"w") as f:
        for l in lines:
            f.write(l)
                
def change_deple(path,newvals,new_path):
    lines = []
    with open(path,"r") as f:
        while(True):
            line_ = f.readline()
            if not line_:
                break
            if(line_.startswith("#")):
                lines.append(line_)
                continue
            line = line_.split()
            if(len(line)==1):
                lines.append(line_)
                continue
            
            if line[2] in newvals.keys():
                val = newvals[line[2]]
                newline =""
                newline = newline+line_[:4]
                newline = newline+eformat(Decimal(val), 4, 3)
                newline = newline+line_[15:]
                lines.append(newline)
            else:
                lines.append(line_)
    with open(new_path,"w") as f:
        for l in lines:
            f.write(l)
            
            
def generate_cco_abundances(lims=[-3,-6]):
    num1 = np.random.uniform(lims[0],lims[1])
    num2 = np.random.uniform(lims[0],lims[1])

    return [10**num1,10**num2]
                
def generate_scenario(paramspath,outdir,name):
    exn = paramspath.split(".")[-1]
    root = os.path.basename(paramspath).split("_")[0]
    vals = generate_cco_abundances()
    indir = os.path.dirname(paramspath)
    
    
    
    file_speci = indir+"/"+root+"_Speci."+exn
    editspeci = {"c":vals[0],"c+":vals[1],"o":vals[0]+vals[1]}
    new_file_speci = outdir+"/"+name+"_Speci."+exn
    change_speci(file_speci, editspeci, new_file_speci)
    
    
    file_deple = indir+"/"+root+"_Deple."+exn
    editdeple = {"C":vals[0]+vals[1],"O":vals[0]+vals[1]}
    new_file_deple = outdir+"/"+name+"_Deple."+exn
    change_deple(file_deple, editdeple, new_file_deple)
    
    
    
    
    
    file_chemi = indir+"/"+root+"_Chemi."+exn
    new_file_chemi = outdir+"/"+name+"_Chemi."+exn
    try:
        #print(file_add,new_file_add)
        copyfile(file_chemi,new_file_chemi)
    except:
        print("Seems like there was a problem trying to copy ",file_chemi)

    new_file_params = outdir+"/"+name+"_Param."+exn
    try:
        #print(file_add,new_file_add)
        copyfile(paramspath,new_file_params)
    except:
        print("Seems like there was a problem trying to copy ",paramspath)


    
def gen_execute(paramspath,outdir,name):
    
    generate_scenario(paramspath,outdir,name)
    
    #bashCommand = 'ls -l' #.format(name)
    with open('{}{}_exec.txt'.format(outdir,name), 'w') as f:
        
        subprocess.call(['./Chem_Evol {}'.format(name)],cwd="CHIMES_0.6/",shell=True,stdout=f)
    
    subprocess.call(['echo "$(tail -5 {}_exec.txt)" > {}_exect.txt'.format(name,name)],cwd=outdir,shell=True)
    subprocess.call(['rm {}_exec.txt'.format(name)],cwd=outdir,shell=True)
        
        #subprocess.call('ls -l',shell=True)
#        process = subprocess.Popen(['cd','CHIMES_0.6/'], stdout=f,cwd="/")
        #process = subprocess.Popen(['./Chem_Evol', name], stdout=f)
  #  process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
  #  output, error = process.communicate()
    

gen_execute("CHIMES_0.6/Data/iterCorr1_Param.dat","CHIMES_0.6/Data/","testout")              