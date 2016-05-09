import smtplib
from os.path import expanduser
from sys import path
import socket
import subprocess,os
from time import strftime
   
def send_email(text,subject='Elmer has finished',to=['kehrl@uw.edu']):
    '''
    Send email.
    '''
    
    from perpython import FROM,password

    message = """\From: %s\nTo: %s\nSubject: %s\n\n%s
            """ % (FROM, ", ".join(to), subject, text)
    try:
        server=smtplib.SMTP('smtp.gmail.com:587')
        server.ehlo()
        server.starttls()
        server.login(FROM,password)
        server.sendmail(FROM,to, message)
        server.close()
    except:
        print "failed to send email"

def run_elmer(sif_file,n=20,to=['kehrl@uw.edu']):
  os.system('echo '+sif_file+'>ELMERSOLVER_STARTINFO')
  
  cmd=['mpiexec', 'ElmerSolver_mpi']
    
  outfile=open(sif_file+'.log','w')
  try:
    p=subprocess.Popen(cmd,stdout=outfile,stderr=subprocess.STDOUT)
    returncode=p.wait()
    outfile.close()
  except KeyboardInterrupt:
    returncode = 'interrupt'

  try:
    os.remove('ELMERSOLVER_STARTINFO')
  except:
    pass
  now=strftime("%Y-%m-%d %H:%M:%S")
  if returncode==0:
	  tail=subprocess.Popen(['tail','-n','2',sif_file+'.log'],stdout=subprocess.PIPE).stdout.read()
	  subject=sif_file+' finished successfully'
	  text=str(sif_file)+' finished running on '+socket.gethostname()+' at '+now+'.\n Last two lines of the log file are:\n'+tail
	  send_email(text,subject=subject,to=to)
  elif returncode=='interrupt':
	  print 'You chose to stop running '+sif_file
  else:
	  tail=subprocess.Popen(['tail','-n','10',sif_file+'.log'],stdout=subprocess.PIPE).stdout.read()
	  subject=sif_file+' returned '+str(returncode)+' at '+now
	  text=str(sif_file)+' exited on '+socket.gethostname()+' at '+now+'.\n Last ten lines of the error log file are:\n'+tail
  print subject
  print text
  return returncode