from astropy.table import Table
import numpy as np

def main(file_name):
    ID = read_data(file_name)
    memb_prob = np.zeros(len(ID))
    for i in range(len(ID)):
        memb_prob[i] = 0.5
    return (ID, memb_prob)

def read_data(file_name):
    data = Table.read(file_name, format='ascii')
    return data['ID']

if __name__ == '__main__':
    main()