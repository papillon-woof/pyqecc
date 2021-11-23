from .stabilizer import *
class CONC(SC):
    def __init__(self,Codeinstances):
        # codelength chack
        self._l = len(Codeinstances)
        for i in range(1,len(Codeinstances)):
            if Codeinstances[i].k != Codeinstances[i-1].n:
                ValueError("Error:codelength mismatch")
        self._Codeinstances = Codeinstances
        super.__init__(self.Codeinstance[self.l-1].n,self.Codeinstance[0].k)

    def hard_decoding(self,s):
        pass
    
    def belief_propagation(self,s):
        pass

    @property
    def l(self):
        return self._l
    
    @property
    def Codeinstances(self):
        return self._Codeinstances