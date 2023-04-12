import sys
import pymol
import numpy as np
import scipy
import matplotlib
import Bio
import pickle

print('Python Version:', sys.version)
print(pymol.get_version_message())
print('numpy version:', np.__version__)
print('scipy version:', scipy.__version__)
print('matplotlib version:', matplotlib.__version__)
print('biopython version', Bio.__version__)


def reader(loc):
    filo = open(loc,'rb')
    data = pickle.load(filo)
    filo.close()
    return data

loc_info_dict_me = r'D:\Work\DNA Methylation patterns\Results\Steric_Clash\info_nbr_dict'
loc_info_dict_kerstin = r'D:\Work\Kerstin work\Methylation_nucleosome\Methylation_nucleosome\1_steric_clash\info_nbr_dict'

print(' ')
print('New = ', reader(loc_info_dict_me))
print('Old = ', reader(loc_info_dict_kerstin))
print(' ')

loc_clash_dict_me = r'D:\Work\DNA Methylation patterns\Results\Steric_Clash\clash_dict'
loc_clash_dict_kerstin = r'D:\Work\Kerstin work\Methylation_nucleosome\Methylation_nucleosome\1_steric_clash\clash_dict'

print(' ')
print('New = ', reader(loc_clash_dict_me)[1]['steric_clash_list']['X_ALA669'])
# print('Old = ', reader(loc_clash_dict_kerstin)[1].keys())
print(' ')
filo = open(loc_clash_dict_kerstin,'rb')
data = pickle.load(filo, encoding='bytes')
print('Old = ', data[1][b'steric_clash_list'][b'X_ALA669'])
filo.close()
print(' ')