import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import os
from shutil import copyfile
import subprocess
import time
import pandas as pd
from datetime import datetime,timedelta

CHIMESPATH ="CHIMES_0.6/"    
excluded = ['zeta','nh','T(K)','bthb','ynn','bnn','bee','bchim','gamh2','gamgr','gampeg','zlambd','xlamoh','xlamco','xlah2o','xlamc','xlamcp','xlamo','xlambn','wthb','xlambe','gamrc']  




def plot_info(path_info,path_img):
    lines =[]
    with open(path_info,"r") as f:
        lines = f.readlines()
    
    lines = [l for  l in lines if not l.startswith("#")]
    lines = [l.split()[1:] for l in lines]
    x,y = np.array([[float(l[0].split(":")[1]),float(l[1].split(":")[1])] for l in lines]).T
    
    
    #PLOTTING
    plt.figure(figsize=(8,8))
    plt.scatter(x,y,s=5)
    plt.xlabel("$\chi_C$",fontsize=15)
    plt.ylabel("$\chi_{C+}$",fontsize=15)
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig(path_img,bbox_inches="tight")
    
    return x,y
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

    return vals

def gen_execute(paramspath,outdir,name):

    vals = generate_scenario(paramspath,outdir,name)

    #bashCommand = 'ls -l' #.format(name)
    with open('{}{}_exec.txt'.format(outdir,name), 'w') as f:

        #print(os.system("dir"))
        subprocess.call(['./Chem_Evol {}'.format(name)],cwd="CHIMES_0.6/",shell=True,stdout=f)
        #
        subprocess.call(['echo "$(tail -5 {}_exec.txt)" > ../Out/{}/{}_depldata.txt'.format(name,name,name)],cwd=outdir,shell=True)
        subprocess.call(['rm {}_exec.txt'.format(name)],cwd=outdir,shell=True)
        
        
        #subprocess.call('ls -l',shell=True)
        #        process = subprocess.Popen(['cd','CHIMES_0.6/'], stdout=f,cwd="/")
                #process = subprocess.Popen(['./Chem_Evol', name], stdout=f)
          #  process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
          #  output, error = process.communicate()
    return vals



def generate_train_sample(N,paramspath,paramdir,outdir,name):
    timestart = time.time()
    
    paramdir = CHIMESPATH+paramdir
    outdir = CHIMESPATH+outdir
    
    print("Creating {} samples with name {} :".format(N,name))
    
    startdate = datetime.today().strftime('%m%d')
    name = name+"_"+startdate
    names = []
    description=[]
    for i in range(1,N+1):
        print("generating sample ",i,"...")
        newname = name+"_"+str(i)
        names.append(outdir+name+"/"+newname+"/"+newname)
        vals = gen_execute(paramspath,paramdir,newname)
        new_desc = "{}  c:{:.3e}  c+:{:.3e}  o:{:.3e}".format(newname,vals[0],vals[1],vals[0]+vals[1])
        description.append(new_desc)
    print("cleaning...")
    subprocess.call(['mkdir {}'.format(name)],cwd=paramdir,shell=True)
    subprocess.call(['mv {}_* {}'.format(name,name)],cwd=paramdir,shell=True)
    subprocess.call(['mkdir {}'.format(name)],cwd=outdir,shell=True)
    subprocess.call(['mv {}_* {}'.format(name,name)],cwd=outdir,shell=True)
    
    print("generating log file...")
    elapsed = time.time() - timestart
    
   
    with open("{}{}/{}.info".format(outdir,name,name),"a+") as ff:
        subprocess.call(['var="$(grep "Steady state has been computed" {}/*/*_depldata.txt | wc -l)"; echo "# $var models reached a stready state" >> {}/{}.info'.format(name,name,name)],cwd=outdir,shell=True)
        ff.write("# Elapsed_time  hh:mm:ss  {}\n".format(timedelta(seconds=elapsed)))
        for l in description:
            ff.write(l+"\n")
            
    plot_info("{}{}/{}.info".format(outdir,name,name), "{}{}/{}.png".format(outdir,name,name))
    
    print("Done.")
    return {"filenames":names,"dirname":name}



# build the database
def reformat_graph_file(name_in,name_out):
    with open(name_out,"w") as k:
        with open(name_in,"r") as f:
            while True:
                line = f.readline()
                if(not line):
                    break
                line=line.strip()
                k.write(",".join(line.split()).strip()+"\n")


def get_graph_data(name,ab_list=None,include_t = True,exclude=None):
    reformat_graph_file(name+".graph",name+"_rf.graph")
    
    if exclude:    
        headers = [*pd.read_csv(name+"_rf.graph",sep=",", nrows=1)]
        data = pd.read_csv(name+"_rf.graph",sep=",", usecols=[c for c in headers if c not in exclude])
    else:
        data = pd.read_csv(name+"_rf.graph",sep=",",engine="python")
    #if exclude:
    #    data = data.drop(exclude,axis=1)
    col_list = []
    if not ab_list:
        return data
    elif include_t:
        col_list.append("t(Myrs)")
    col_list = ["t(Myrs)"]+["t(yrs)"] + ab_list
    data.drop(data.tail(1).index,inplace=True) # drop last n=1 rows (supressing final drop (???) )
    data.drop(data.head(1).index,inplace=True) # drop first n=1 rows (supressing initial peak (???) )
    data["t(yrs)"] = data["t(Myrs)"]*1e6
    return data[col_list]


#convfactor from Myrs to seconds (s/Myrs)
def get_deriv_data(name,ab_list=None,exclude = None,include_t = True,convfactor =31556925216000):
    reformat_graph_file(name+".deriv",name+"_rf.deriv")
       
    if exclude:    
        headers = [*pd.read_csv(name+"_rf.deriv",sep=",", nrows=1)]
        data = pd.read_csv(name+"_rf.deriv",sep=",", usecols=[c for c in headers if c not in exclude])
    else:
        data = pd.read_csv(name+"_rf.deriv",sep=",",engine="python")

    #data = pd.read_csv(name+"_rf.deriv",sep=",",engine="python")
    #if exclude:
    #    data = data.drop(exclude,axis=1)
    col_list = []
    if not ab_list:
        return data
    elif include_t:
        col_list.append("t(Myrs)")
    col_list = ["t(Myrs)"]+["t(yrs)"] + ab_list
    data.drop(data.tail(1).index,inplace=True) # drop last n=1 rows (supressing final drop (???) )
    data.drop(data.head(1).index,inplace=True) # drop first n=1 rows (supressing initial peak (???) )
    data["t(yrs)"] = data["t(Myrs)"]*1e6

    for n in ab_list:
        data[n] = data[n].apply(lambda x : x*convfactor)

    return data[col_list]



def compile_database(dirname,N,pref=None,ab_list=None,exclude = None,include_t = True):
    #dfs = []
    outname = CHIMESPATH + "Out/" +dirname +"/"+dirname
    subprocess.call(['mkdir {}_csv'.format(outname)],shell=True)
    names = [CHIMESPATH+"Out/"+dirname+"/"+x+"/"+x for x in os.listdir(CHIMESPATH+"/Out/"+dirname)  if not x.endswith("info") and not x.endswith("csv") and not  x.endswith("png")]


    ("compiling database of ",len(names)," samples")
    for n in names:
        iter = n.split("_")[-1]
        print("compiling sample "+iter+"...")
        name = n if not pref else pref+n
        derivdata = get_deriv_data(name,ab_list=ab_list,exclude=exclude,include_t=include_t)
        graphdata = get_graph_data(name,ab_list=ab_list,exclude=exclude,include_t=include_t)
        derivdata["datatype"]="deriv"
        graphdata["datatype"]="graph"

        newdat = pd.concat([graphdata,derivdata])
        del derivdata
        del graphdata
        #print(newdat.columns)
        newdat["name"]=iter
        #dfs.append(newdat)
        fileout = '{}_csv'.format(outname)+"/"+n.split("/")[-1]+".csv.gz"
        #print(fileout)
        newdat.to_csv(fileout, index=False,compression="gzip")

    #dfs = [df.set_index('id') for df in dfs]
    #print("Merging dataframe...")
    #finaldf = pd.concat(dfs)


    #print(finaldf.describe())
    #print(finaldf.colmns)
    #print(np.shape(finaldf))


    print("done")


