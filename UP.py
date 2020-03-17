from astropy.table import Table
import numpy as np

def main(file_name):
    ID, memb_prob = read_data(file_name)
    for i in range(len(ID)):
        if memb_prob[i]==1:
            memb_prob[i]=0.99
        if memb_prob[i]==0:
            memb_prob[i]=0.01
    return (ID, memb_prob)


def read_data(file_name):
    data = Table.read(file_name, format='ascii')
    return data['ID'], data['probability']

if __name__ == '__main__':
    main()
