
from astropy.table import Table


def main(file_name):
    data = Table.read(file_name, format='ascii')
    return data['ID'], data['probability']
