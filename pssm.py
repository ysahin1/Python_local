import os
import logging
import numpy as np
from Bio import SeqIO
from pssm_lib import config as cfg
import sys
from pssm_lib import workEnv
from pssm_lib.workEnv import TemporaryEnv
from pssm_lib.blast  import runPsiBlast
from pssm_lib.datacache import DataCache

fastafile = "ex.fasta"
dbfile = "db/uniport"

class PSSM():

    def __init__(self,fastafile,dbfile):
        try:
            self.aaOrder = "ARNDCQEGHILKMFPSTWYV"
            self.dbfile = dbfile
            self.pbniter = 3
            self.pbnalign = 5000
            self.pbeval = 0.001
            self.threads = 1
            self.fastafile = fastafile
            self.wim = os.getcwd()
            print ("The instance variables have been assinged")
        except:
            raise print("The instance variables have not been generated")

    def IsFasta(fastafile):
        with open(fastafile, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)

    def GetDataCache(self, cache_dir):
        ret = None
        if cache_dir is not None:
            if os.path.isdir(cache_dir):
                ret = DataCache(cache_dir)
        return ret

    def runBlast(self):
        data_cache = self.GetDataCache("/".join(self.wim+"deneme"))
        try:        
            with open(self.fastafile, "r") as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    print (self.wim)
                    workEnv = TemporaryEnv()
                    prefix = record.id
                    fastaSeq  = workEnv.createFile(prefix+".",".fasta")
                    print (fastaSeq)
                    SeqIO.write([record], fastaSeq, 'fasta')
                    pssmFile = runPsiBlast(prefix,
                                            self.dbfile, 
                                            fastaSeq, 
                                            workEnv=workEnv, 
                                            data_cache=data_cache,
                                            num_alignments= self.pbnalign, 
                                            num_iterations=self.pbniter, 
                                            evalue=self.pbeval,
                                            threads=self.threads)
                    if os.path.isfile(fastaSeq):
                        logging.info("PSSM file has been constructed")
                        workEnv.destroy()
            handle.close()        
        except: 
            logging.error("PSSM file has not been constructed")   
pssm = PSSM(fastafile,dbfile)
run = pssm.runBlast()