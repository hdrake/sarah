import numpy as np
import struct

def f(ϕ): return 2. * (2*np.pi)/(60**2 * 24.) * np.sin(np.deg2rad(ϕ))

def header(file):
    struct.unpack('i',file.read(4))
    
def unpack_float(file, dims=(1)):
    if dims==(1):
        return struct.unpack('d',file.read(8))[0]
    else:
        dtype = np.dtype(np.float64)
        dtype = dtype.newbyteorder('<')
        return np.frombuffer(file.read(8*np.product(dims))).reshape(dims)

def unpack_int(file):
    return struct.unpack('i',file.read(4))[0]

def load_channel(filename):
    class channel():
        pass

    with open(filename, mode="rb") as file:
        # Transfrom (x,y) = (y, -x) to convert from Helfrich to Pratt coordinates
        header(file)
        channel.nout = unpack_int(file)
        channel.tend = unpack_float(file)
        header(file)
        header(file)
        channel.my = unpack_int(file)
        channel.mx = unpack_int(file)
        channel.meqn = unpack_int(file)
        channel.mbc = unpack_int(file)
        header(file)
        header(file)
        channel.dt = unpack_float(file)
        channel.dy = unpack_float(file)
        channel.dx = unpack_float(file)
        gamma = unpack_float(file)
        channel.cf = unpack_float(file)
        channel.cb = unpack_float(file)
        header(file)
        header(file)
        channel.iybc = unpack_int(file)
        channel.ixbc = unpack_int(file)
        header(file)
        header(file)
        channel.y = unpack_float(file, dims=(channel.my,1))
        header(file)
        header(file)
        channel.x = -unpack_float(file, dims=(channel.mx,1))[::-1]
        header(file)
        header(file)
        channel.h = unpack_float(file, dims=(channel.mx,channel.my))[::-1, :]
        header(file)
        header(file)
        channel.cbxy = unpack_float(file, dims=(channel.mx,channel.my))[::-1, :]
        header(file)

        channel.t = np.zeros((channel.nout))
        channel.d = np.zeros((channel.nout, channel.mx, channel.my))
        channel.V = np.zeros((channel.nout, channel.mx, channel.my))
        channel.U = np.zeros((channel.nout, channel.mx, channel.my))

        for i in range(channel.nout):
            try:
                header(file)
                channel.t[i] = unpack_float(file)
                header(file)
                header(file)
                channel.d[i,:,:] = unpack_float(file, dims=(channel.mx,channel.my))[::-1, :]
                header(file)
                header(file)
                channel.V[i,:,:] = unpack_float(file, dims=(channel.mx,channel.my))[::-1, :]
                header(file)
                header(file)
                channel.U[i,:,:] = -unpack_float(file, dims=(channel.mx,channel.my))[::-1, :]
                header(file)
            except:
                print("Broken file")
                break
        
        channel.y = channel.y[:,0][:,np.newaxis]
        channel.x = channel.x[:,0][np.newaxis,:]

        channel.v = channel.V/channel.d
        channel.u = channel.U/channel.d
        
        channel.d_um = np.copy(channel.d)
        
        mask_depth = 1e-5
        channel.v[channel.d <= mask_depth] = np.nan
        channel.u[channel.d <= mask_depth] = np.nan
        channel.V[channel.d <= mask_depth] = np.nan
        channel.U[channel.d <= mask_depth] = np.nan
        channel.d[channel.d <= mask_depth] = np.nan
        channel.η = channel.d + channel.h
        channel.η[channel.d <= mask_depth] = np.nan
        
    return channel