import os
import numpy as np
from Bio import SeqIO
from pssm_lib import config as cfg
import sys
from pssm_lib import workEnv
from pssm_lib.blast  import runPsiBlast

class PSSM(self, fastafile, dbfile):
    print ("Starting")
    def __init__:

        if is_fasta:
            self.fastafile = fastafile
            self.aaOrder = "ARNDCQEGHILKMFPSTWYV"
            self.dbfile = dbfile
            self.pbniter = 3
            self.pbnalign = 5000
            self.pbeval = 0.001
            self.threads = 1
            print ("the variables have been constructed")
        else:
            print ("something is wrong")
    def is_fasta(fastafile):
        with open(fastafile, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)
    try:
        with open(self.fastafile, "r") as handle:
            for record in SeqIO.read(handle, 'fasta'):
                workEnv = TemporaryEnv()
                prefix = record.id.replace("|","_")
                seq = record.seq
                fastaSeq  = workEnv.createFile(prefix+".", ".fasta")
                SeqIO.write([record], fastaSeq, 'fasta')
                pssmFile = runPsiBlast(prefix,
                                        self.dbfile, 
                                        fastaSeq, 
                                        workEnv, 
                                        data_cache=data_cache,
                                        num_alignments=ns.pbnalign, num_iterations=ns.pbniter, evalue=ns.pbeval,
                                        threads=ns.threads)
    except: 
        raise    
    
    def _check(self,line): # it checks if the file is PSSM file
        import re
        if not re.search('Last position-specific scoring matrix computed', line):
            raise InvalidCheckpointFileError

    def logistic(x):
        # return math.tanh(x)
        return 1 / (1 + numpy.exp(-x))

    class InvalidCheckpointFileError(Exception):
        def __init__(self):
            pass

    def seq_to_pssm(self):
        X = []
        for seq in self.seq_dic:
            for aa in self.seq_dic[seq]:
                aa == aa.upper()
                x = [0.0]*len(self.aaOrder)
                try:
                    x[self.aaOrder.index(aa)]=1.0
                except:
                    pass
                X.append(x)
        return np.array(X)

    def _pssmParseNew(self,checkpoint, transform):
        
        headerSize = 3
        footerSize = -6

        try:
            _check(checkpoint[1])
        except:
            raise

        pssm = []
        for line in checkpoint[headerSize:footerSize]:
            line = line.split()[2:22]
            pos = numpy.zeros(20)
            for j in range(20):
                if transform:
                    pos[j] = logistic(float(line[j]))  # 1.0 / (1.0 + math.exp(-1 * float(line[j + shift])))
                else:
                    pos[j] = float(line[j])
            pssm.append(pos)
        return numpy.array(pssm)

    def BlastCheckPointPSSM(self,checkpointFile,newFormat = True, transform = True):
        try:
            checkpointFile = open(checkpointFile).readlines()
        except IOError:
            print("Error while open/reading checkpoint file.")
            raise
        pssm = None
        if newFormat:
            try:
                pssm = _pssmParseNew(checkpointFile, transform)
            except:
                raise
        return pssm

fastafile = "/media/yunus/TOSHIBA1TB/Python_local/datasets/example_dataset.fasta"
dbfile = "/home/yunus/uniport_database/uniprot_sprot.pep" 
pssm = PSSM(fastafile,dbfile)
print (pssm)