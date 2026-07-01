import os, sys, subprocess

# Script by G.Marozzo to create jobs based on CMS samples. Usage: python3 CreateJobs.py NAMEOFMACRO NAMEOFSAMPLE

def createdirs(num,dataset,identifier,macro):
  os.system("mkdir "+outputdir+dataset+"_Chunk"+str(counter))
  os.system("cp -r JobTemplate "+dataset+"_Chunk"+str(counter))
  with open(dataset+"_Chunk"+str(counter)+"/batchScript.sh", 'r') as batchscript:
    script = batchscript.read()
  script = script.replace("MACRO",macro.strip())
  script = script.replace("FILELOC",line.strip())
  script = script.replace("IDENTIFIER",sample.split(",")[3])
  script = script.replace("XSEC",sample.split(",")[4].strip())
  script = script.replace("NUMEVENTS",numevents)
  with open(dataset+"_Chunk"+str(counter)+"/batchScript.sh", 'w') as batchscript:
    batchscript.write(script)
  print("Created "+dataset+"_Chunk"+str(counter))


macro = sys.argv[1].strip("macros/")
samplelist = sys.argv[2]

outputdir = "/eos/user/g/gmarozzo/JobOutput/" #TO BE CHANGED TO YOUR AREA

if not os.path.exists(outputdir):
  os.mkdir(outputdir)
  
samples = open(samplelist,"r")
for sample in samples:
  if sample.startswith("#") or sample.startswith("-"): continue
  os.system("dasgoclient --query='file dataset="+sample.split(",")[2]+"' > temp.txt")
  result = subprocess.check_output("dasgoclient --query='dataset="+sample.split(",")[2]+" | grep dataset.nevents'",shell=True,text=True)
  os.system("dasgoclient --query='dataset="+sample.split(",")[2]+" | grep dataset.nfiles'")
  numevents = result.strip()
  print(sample.split(",")[1]+"_"+sample.split(",")[0]+" : "+numevents+" events")
  filelist = open("temp.txt","r")
  for counter, line in enumerate(filelist):
    createdirs(counter,sample.split(",")[1]+"_"+sample.split(",")[0],sample.split(",")[3],macro)
    #if counter>0: break
  filelist.close()
  os.system("rm temp.txt")

samples.close()
