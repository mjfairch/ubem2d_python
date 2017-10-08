import numpy as np

__all__ = ['read_data']

def read_data(path):
    '''
    This returns a numpy array from text data arranged in columns of a file.
    '''
    with open(path,'r') as fid:
        line_count = sum(1 for line in fid)
        fid.seek(0)
        column_count = sum(1 for x in fid.readline().split())
        fid.seek(0)
        data = np.zeros((line_count,column_count))
        for i in range(line_count):
            data[i,:] = fid.readline().split()
    return data
