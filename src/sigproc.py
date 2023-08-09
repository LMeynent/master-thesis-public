import numpy as np
import scipy as sp


def shift(x, k, fill=np.nan):
    if fill == 'edge':
        fill = x[-1] if k < 0 else x[0]

    if k > 0:
        return np.concatenate((np.ones((k,)) * fill, x[:-k]))
    else:
        return np.concatenate((x[-k:], np.ones((-k,)) * fill))



def moving_avg(x, w=7, centre=True):
    if centre:
        if w%2 == 0:
            raise ValueError(f'Cannot centre window if window size is even (w={w})')
        return np.convolve(np.pad(x, (w//2, w//2), 'edge'), np.ones(w)/w, mode='valid')
    else:
        return np.convolve(np.pad(x, (w-1, 0), 'edge'), np.ones(w)/w, mode='valid')


def filter(x, h):
    return np.convolve(x, h, 'same')[:len(x)]