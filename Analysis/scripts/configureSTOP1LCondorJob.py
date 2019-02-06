#!/bin/env python

import sys
import imp
import copy
import os
import filecmp
import shutil
import pickle
import math
import pprint
import subprocess
from datetime import date
from optparse import OptionParser


class BatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_option("--batchqueue", type="string", help="Batch queue")
      self.parser.add_option("--batchscript", type="string", help="Name of the HTCondor script")
      self.parser.add_option("--outdir", type="string", help="Name of the output directory")
      self.parser.add_option("--outlog", type="string", help="Name of the output log file")
      self.parser.add_option("--errlog", type="string", help="Name of the output error file")

      self.parser.add_option("--script", type="string", help="Name of the script to run")
      self.parser.add_option("--fcn", type="string", help="Name of the function in the script")
      self.parser.add_option("--fcnargs", type="string", help="The arguments of the function")

      self.parser.add_option("--dry", dest="dryRun", action="store_true", default=False, help="Do not submit jobs, just set up the files")
      self.parser.add_option("--interactive", dest="interactive", action="store_true", default=False, help="Do not submit jobs; run them interactively")


      (self.opt,self.args) = self.parser.parse_args()

      optchecks=[
         "batchqueue",
         "batchscript",
         "outdir",
         "outlog",
         "errlog",
         "script",
         "fcn",
         "fcnargs"
      ]
      for theOpt in optchecks:
         if not hasattr(self.opt, theOpt) or getattr(self.opt, theOpt) is None:
            sys.exit("Need to set --{} option".format(theOpt))

      if not os.path.isfile(self.opt.script):
         sys.exit("Script {} does not exist. Exiting...".format(self.opt.script))

      self.submitJobs()


   def produceCondorScript(self):
      currentdir = os.getcwd()
      scriptcontents = """executable              = {batchScript}
         arguments               = {mainDir} $(ClusterId)$(ProcId)
         output                  = {outLog}
         error                   = {errLog}
         log                     = $(ClusterId).$(ProcId).log
         Initialdir              = (outDir)
         request_memory          = 4000M
         +JobFlavour             = "tomorrow"
         x509userproxy           = {home}/x509up_u{uid}
         #https://www-auth.cs.wisc.edu/lists/htcondor-users/2010-September/msg00009.shtml
         periodic_remove         = JobStatus == 5
         WhenToTransferOutput    = ON_EXIT_OR_EVICT
         """
      scriptcontents = scriptcontents.format(
         home=os.path.expanduser("~"), uid=os.getuid(),
         mainDir=currentdir, batchScript=self.opt.batchscript,
         outLog=self.opt.outlog, errLog=self.opt.errlog,
         outDir=self.opt.outdir
         )
      self.condorScriptName = "condor.sub"
      condorScriptFile = open(self.condorScriptName,'w')
      condorScriptFile.write(scriptcontents)
      condorScriptFile.close()



   def submitJobs(self):
      self.produceCondorScript()

      strscrcmd = argstr.format(channel=ch,category=cat,achypothesis=hypo,systematic=syst,frmethod=frm)
      strscrcmd = strscrcmd.replace(' ','') # The command passed to bash script should not contain whitespace itself
      jobcmd = "submitHiggsWidthTemplateStageGeneric.sh {} \({}\) {}".format(self.fcnname, strscrcmd, self.opt.batchqueue)
      if self.opt.interactive:
         jobcmd = "root -l -b -q -e \"gROOT->ProcessLine(\\\".x loadLib.C\\\");gROOT->ProcessLine(\\\".x {}.c+({})\\\");\"".format(self.fcnname, strscrcmd)
      if self.opt.dryRun:
         jobcmd = "echo " + jobcmd
      ret = os.system( jobcmd )



if __name__ == '__main__':
   batchManager = BatchManager()
