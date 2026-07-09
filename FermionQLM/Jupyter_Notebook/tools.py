import numpy as np
import math

def order_of_magnitude2(a_value):
    #return 2
    if np.abs(a_value) < 1.0 and a_value != 0:
        m = np.abs(np.log10(np.abs(a_value)))
        return int(max(math.ceil(m) + 1., 2.))
    else: 
        return 2

def order_of_magnitude(a_value):
    x = a_value - int(a_value)
    x = np.round(x, 8)
    num_str = f"{x}"
    num_str = num_str[2:]
    _size = len(num_str)
    if num_str == "0":
        _size = 0;
    
    return _size
    
def info_raw(Lx, Ly, Nu, Nd, tau, mass):
    arr = [tau, mass]
    names = ['tau', 'm']
    if Nd < 0: 
        Nd = Lx*Ly - Nu
    if Nu < 0: 
        Nu = Lx*Ly - Nd
    
    info = "_Lx=%d,Ly=%d,Nu=%d,Nd=%d"%(Lx,Ly,Nu,Nd)
    for i, var in enumerate(arr):
        # print(names[i], var)
        n = order_of_magnitude( np.abs(var) )
        if np.abs(var) < 1e-15:
            info += str(",%s=0"%(names[i]))
        elif var < 0:
            info += str(",%s=-{:.%df}"%(names[i], n)).format(round(np.abs(var), n))     
        else:
            info += str(",%s={:.%df}"%(names[i], n)).format(round(np.abs(var), n))
    return info

def info(Lx, Ly, Nu, Nd, tau, mass, ext = '.hdf5'):
    return info_raw(Lx, Ly, Nu, Nd, tau, mass) + ext


def legend_label(name, value):
    n = order_of_magnitude(value)
    return str("$%s={:.%df}$"%(name, n)).format(round(value, n))


def diff_central(x, y):
    x0 = x[:-2]
    x1 = x[1:-1]
    x2 = x[2:]
    y0 = y[:-2]
    y1 = y[1:-1]
    y2 = y[2:]
    f = (x2 - x1)/(x2 - x0)
    return (1-f)*(y2 - y1)/(x2 - x1) + f*(y1 - y0)/(x1 - x0)
