from pyqecc import *
NUM_OF_CONCATENATE = 3
for num_of_concatenate in range(1,NUM_OF_CONCATENATE+1):
    conc_code = [FiveCode()]
    for i in range(1,num_of_concatenate):
        conc_code += [ParaCode([FiveCode() for i in range(5 ** i)])]
    my_code = ConcCode(conc_code)
    print(my_code)
    dec_sim(my_code,channel_instance=DepolarizingChannel(my_code.n,p=[0.13, 0.15, 0.17, 0.18, 0.1885, 0.19]),MONTE=5000)