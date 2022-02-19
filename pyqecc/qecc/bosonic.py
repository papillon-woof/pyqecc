import numpy as np
from .abstruct import *
from ..util import *

from .stabilizer import SC

class GKP(CODE):
    NAME = "GKP qubit"
    USAGE = 'GKP_STABILIZER_CLASS'

    def __init__(
        self,
        code_instance,
        sigma = 0.1
    ):
        self.decoder_output = {
            "LT": None,
            "LOGICAL_ERROR_PROBABILITY": None,
        }
        self._code_instance = code_instance
        self._sigma = sigma
        super().__init__(self.code_instance.n,self.code_instance.k)
        # create analog infomation calc matrix
        #self.matrix_for_genarating_LLR = util.gaussjordan(np.concatenate([self.code_instance.H,self.code_instance.L],0),change=True)
        self.llr = np.zeros(self.n)

    def set_channel_param(self, sigma):
        self._sigma = sigma

    # def get_analog_information(self,E):
    #     for i in range(self.n):
    #         count=0
    #         for j in range(self.n):
    #             if 1==self.matrix_for_genarating_LLR[i,j]:
    #                 self.llr[i]+=((-1)**(count%2))*E[i][j]
    #     llr = np.log(self.llr)
    #     return llr

    def get_syndrome(self,channel_output,digital=False):
        '''
        引数: E: アナログ雑音もしくは誤り
        Eは，実際に生じたシフトである．
        '''
        err = np.zeros(2*self.n,dtype="i1")
        
        # 2√π>|E|>√π => error if 2√π < or < -2√π => +-2√π 
        delta = dmods(channel_output["DELTA"],2*np.sqrt(np.pi))
        e_pos = np.where(delta>np.sqrt(np.pi)/2)[0]
        err[e_pos]=1
        syndrome = self.code_instance.get_syndrome({"E": err})
        if digital:
            return syndrome

        # calculation Δm
        delta_m = np.abs(delta.copy())
        delta_m[e_pos] = (np.sqrt(np.pi)-np.abs(delta))[e_pos]
        return syndrome, delta_m

    def in_S(self,s):
        return self.code_instance.in_S(s)

    def calc_llr(self,val,sigma):
        return np.log(np.exp((val**2)/(2*sigma**2))/np.exp(((np.sqrt(np.pi)-np.abs(val))**2)/(2*sigma**2)))

    def analog_syndrome_decode(self,information):
        syndrome = information[0]
        delta_m = information[1]
        most_likely_error = np.zeros(self.n,dtype='i1')
        mi = 100000

        for i in range(2 * 2 ** self.k):
            lt = self.code_instance.get_T(syndrome)
            for ii in range(2 * self.k):
                lt ^= (1&(i>>ii))*self.code_instance.get_L(i)
            for j in range(2 * 2 ** (self.n-self.k)):
                c = lt.copy()
                for jj in range(2 * self.nk):
                    c ^= (1&(j>>jj))*self.code_instance.get_S(j)
                llr = 0
                for k in range(len(c)):
                    llr += (((-1)**c[k])*self.calc_llr(delta_m[k],self.sigma))
                if mi > llr:
                    mi = llr
                    most_likely_error = c
        self.decoder_output["LT"] = most_likely_error
        #print(most_likely_error)
        return self.decoder_output
    
    def analog_decode(self,information):
        syndrome = information[0]
        delta_m = information[1]
        most_likely_error = np.zeros(self.n,dtype='i1')
        mi = 100000

        for i in range(2 * 2 ** self.k):
            lt = self.code_instance.get_T(syndrome)
            for ii in range(2 * self.k):
                lt ^= (1&(i>>ii))*self.code_instance.get_L(i)
            llr = 0
            for k in range(len(lt)):
                llr += (((-1)**lt[k])*self.calc_llr(delta_m[k],self.sigma))
            if mi > llr:
                mi = llr
                most_likely_error = lt
        self.decoder_output["LT"] = most_likely_error
        #print(most_likely_error)
        return self.decoder_output

    def decode(self,syndrome):
        
        # do not use the analog information
        if len(syndrome) == 1:
            return super().decode(syndrome, analog=False)
        
        #use the analog information
        else:
            return self.analog_decode(syndrome)

    @property
    def code_instance(self):
        return self._code_instance

    @property
    def sigma(self):
        return self._sigma