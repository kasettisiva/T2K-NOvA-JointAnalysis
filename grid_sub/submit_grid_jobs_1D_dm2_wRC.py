from DIRAC.Core.Base import Script
Script.parseCommandLine(ignoreErrors=False)

from DIRAC.Interfaces.API.Dirac import Dirac
from DIRAC.Interfaces.API.Job import Job

from DIRAC import gLogger

import os
import tarfile

job_name = "p-theta compilation"
#tar_sourcefolder_fname = "P-theta.gzip"
tar_sourcefolder_fname = "Packages.gzip"
tar_archive_script = "untar_script.sh"
main_script = "main_script.sh"
runFit_script = "margTemplates_1D_dm2_wRC_30June2022.sh"
test_script = "submit_1D_dm2_wRC.sh"
mode = "wms"
input_files = ["/home/kasetti/T2K_on_grid/P-theta"]
#dest = "LCG.IN2P3-CC.fr"
dest = "LCG.UKI-LT2-QMUL.uk"
#dest = "LCG.UKI-SCOTGRID-GLASGOW.uk"

gLogger.info('creating archive {}'.format(tar_sourcefolder_fname))
with tarfile.open(tar_sourcefolder_fname, "w:gz") as tar:
        for item in input_files:
                gLogger.info("adding {} as {}".format(item, item))
                tar.add(item, arcname=item)
"""
with open(tar_archive_script, 'w') as stream:
        stream.write("#!/bin/bash\n")
        stream.write("echo 'Untar step:'\n")
        stream.write("ls -atr\n")
        stream.write("python -c 'import tarfile; tar = tarfile.open(\"Packages.gzip\"); tar.extractall(); tar.close()'\n".format(tar_sourcefolder_fname))
        stream.write("python -c 'import tarfile; tar = tarfile.open(\"inputs_Spring2020_0325.gzip\"); tar.extractall(); tar.close()'\n".format(tar_sourcefolder_fname))
        stream.write("echo 'Untar step done'\n")
"""
gLogger.info('Job submitter: {}'.format(job_name))

#list=[1,2,3,5,6,7,8,9,11,12]
#list=[2,3,5,6,8,11,12]

#for i in range(200):
#for i in range(25):
#for i in range(len(list)):
#for i in range(30):
for i in range(200):

  j = Job()
  j.setCPUTime(3000)
  j.setName(job_name)
  j.setLogLevel('debug')
  j.setDestination(dest)
  #j.setExecutable(tar_archive_script)
  j.setExecutable(main_script)
  
  lfn_path='LFN:/t2k.org/user/p/pserv.phys.lsu.edu'
  content=''
  with open(runFit_script, "r") as stream, open(test_script,"w+") as writer:
      count=0
      for line in stream:
          count += 1
          if count == 5:
              content += 'start_bin={}'.format(i)
          content += str(line)
      writer.write(content)

  j.setName('API_%d' % i)
  #j.setInputSandbox([tar_sourcefolder_fname, tar_archive_script, main_script, lfn_path+'/'+'subset{:04d}.root'.format(i*100), lfn_path+'/Validation_ToyXP_AsimovA_2018_02May2021.root', lfn_path+'/nova-sl7-novat2k_latest.sif'])
  #j.setInputSandbox([lfn_path+'/'+tar_sourcefolder_fname, tar_archive_script, main_script, lfn_path+'/Validation_ToyXP_AsimovA_2018_16Feb2022.root', lfn_path+'/P-theta_nova-submodules_gitBranchCommit_Aug2021/MakeTemplates_OA2019_30k_t2k-nova_syscorr_set1_16Feb2022/'+'subset{:04d}.root'.format(i*1000), lfn_path+'/nova-sl7-novat2k_v4_fixcosmicsrock.sif', lfn_path+'/inputs_Spring2020_0325.gzip'])
  #j.setInputSandbox([lfn_path+'/'+tar_sourcefolder_fname, tar_archive_script, main_script, lfn_path+'/Validation_ToyXP_AsimovA_2018_16Feb2022.root', lfn_path+'/nova-sl7-novat2k_v4_fixcosmicsrock.sif', lfn_path+'/inputs_Spring2020_0325.gzip'])
  #j.setInputSandbox([lfn_path+'/'+tar_sourcefolder_fname, tar_archive_script, main_script, lfn_path+'/Validation_ToyXP_Asimov0_2018_15July2022.root', lfn_path+'/nova-sl7-novat2k_v4_fixcosmicsrock.sif', lfn_path+'/inputs_Spring2020_0325.gzip'])
  j.setInputSandbox([lfn_path+'/'+tar_sourcefolder_fname, main_script, lfn_path+'/Validation_ToyXP_Asimov4_2018_03Aug2022.root', lfn_path+'/nova-sl7-novat2k_v4_fixcosmicsrock.tar.bz2', lfn_path+'/inputs_Spring2020_0325.gzip'])

  j.setOutputSandbox(['bins{:04d}.root'.format(i*1000)]);
  #j.setOutputSandbox(['bins{:04d}.root'.format(i*1000), 'subset{:04d}.root'.format(i*1000), 'Syst_1k_t2k-nova_syscorr_{:04d}.root'.format(i*1000)])
  #j.setOutputSandbox(['std.err', 'std.out', 'P-theta_install.gzip', 'bins{:04d}.root'.format(i*100)])
  j.setExecutable(test_script)
  dirac = Dirac()
  jobID = dirac.submitJob(j)
  print(dirac.getJobLoggingInfo(jobID['JobID']))
  print('Submission Result: ',jobID)
  if jobID['OK'] == False:
      error_msg = "".join(jobID["CallStack"])
      gLogger.always("Job not submitted:\n{}".format(error_msg))
      DIRAC_exit(1)
      gLogger.always("Job submitted <{}>: \n{}".format(jobID['JobID'], j._toJDL()))
      gLogger.info("Cleaning up!")
      os.remove(tar_sourcefolder_fname)
      #os.remove(tar_archive_script)
