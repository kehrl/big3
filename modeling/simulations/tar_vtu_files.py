import argparse, os, sys

##########
# Inputs #
##########

# Get inputs to file
parser = argparse.ArgumentParser()
parser.add_argument("-glacier",dest="glacier",required = True, 
        help = "Name of glacier (Kanger or Helheim)")
parser.add_argument("-mesh", dest="mesh", required = True,
        help = "Name of mesh directory") 
parser.add_argument("-runname",dest="runname",required = True, 
        help = "Model runname for deciding which files to tar.")

args, _ = parser.parse_known_args(sys.argv)

RES = args.mesh
glacier = args.glacier
runname = args.runname

DIRM = os.path.join(os.getenv("MODEL_HOME"),glacier+"/3D/"+RES+"/")

#############
# Tar Files #
#############

files = os.listdir(DIRM+"mesh2d")
os.chdir(DIRM+"mesh2d")
for file in files:
  if file.startswith(runname) and file.endswith('.pvtu'):
    i = int(file[len(runname):len(runname)+4])
    os.system('tar -czvf '+runname+'{0:04d}.pvtu.tar.gz '.format(i)+\
                runname+'*{0:04d}.'.format(i)+'*vtu')
    os.system('rm '+runname+'*{0:04d}.'.format(i)+'*vtu')
