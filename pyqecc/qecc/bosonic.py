import numpy as np
from .abstruct import *
from ..util import *

from stabilizer import SC

class StabilizerWithGKP(SC):
    NAME = "GKP qubit"
    USAGE = 'GKP_STABILIZER_CLASS'

    def __init__(
        self,
        n,
        k,
        H=None,
        T=None,
        L=None,
        P=None,
        mode="HD",
        BITWISE=True,
        ANALOG_INFORMATION="No",
    ):
        super().__init__(n,k,H=None,T=None,L=None,P=None,mode="HD",BITWISE=True,)
        # create analog infomation calc matrix
        self.matrix_for_genarating_LLR = util.gaussjordan(np.concatenate([H,L],0),change=True)
        self.llr = np.zeros(self.n)
    def get_digital_information(self,E):
        for i in range(self.n):
            count=0
            for j in range(self.n):
                if 1==self.matrix_for_genarating_LLR[i,j]:
                    self.llr[i]+=((-1)**(count%2))*E[i][j]
        llr = np.log(self.llr)
        return llr

    def get_syndrome(self,E):
        '''
        引数: E: アナログ雑音もしくは誤り
        '''
        if E.dtype == "i1":
            return 
        else:
             analog_information = np.zeros(2 * self._n)#ここに，アナログ値を得る方法を記載
             digital_information = self.get_digital_information(E)
             syndrome = self._codeInstance.get_syndrome(digital_information)
             return syndrome,analog_information

    def gkp_cx(a,b):
        return a-b
    def analog_decode(analog_information):
        

    def decode(self,syndrome):
        # do not use the analog information
        #if len(syndrome)==1:
        #    return super().decode(syndrome, analog=False)
        # use the analog information
        #elif len(syndrome)==2:
        return self.analog_decode(syndrome)
    @property
    def codeInstance(self):
        return self._codeInstance

    @property
    def n(self):
        return self._n

class StabilizerWithGKP_old(SC):
    NAME = "GKP qubit"
    USAGE = 'GKP_STABILIZER_CLASS'
    def __init__(self,n,codeInstance):
        self._n = n
        self._codeInstance = codeInstance
    
    def get_digital_information(self,E):
        pass

    def get_syndrome(self,E):
        '''
        引数: E: アナログ雑音もしくは誤り
        '''
        if E.dtype == "i1":
            return 
        else:
             analog_information = np.zeros(2 * self._n)#ここに，アナログ値を得る方法を記載
             digital_information = self.get_digital_information(E)
             syndrome = self._codeInstance.get_syndrome(digital_information)
             return syndrome,analog_information

    def decode(self,syndrome):
        # do not use the analog information
        if len(syndrome)==1:
            return self.codeInstance.decode(syndrome, analog=False)
        # use the analog information
        elif len(syndrome)==2:
            return self.codeInstance.decode(syndrome[0], mode="ANALOG")

    @property
    def codeInstance(self):
        return self._codeInstance

    @property
    def n(self):
        return self._n