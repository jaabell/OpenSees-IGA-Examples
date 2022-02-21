from matplotlib.pylab import *
import glob
from itertools import product as iterprod
import scipy.signal as sig


from math import log10, floor
def round_to_1(x):
    return round(x, -int(floor(log10(abs(x)))))

def leer_encabezado(fname, maxlines = 38, datadir="./datos"):
    header = {}
    with open(datadir+"/"+fname, "r") as datafile:
        for i, line in zip(range(1, maxlines+1), datafile):
            line = line.replace("\r","")
            line = line.replace("\n","")
            if i == 3:
                header['date'] = line
            if i == 4:
                header['time'] = line
            if i == 10:
                header['units'] = line.split("\t")
            if i == 9:
                header['channel_names'] = line.split("\t")
    return header
        
def get_datos(fname, datadir="./datos"):    
    header = leer_encabezado(fname, datadir=datadir)
    data = genfromtxt(datadir+"/"+fname, skip_header=38)

    i_acc =  [u == 'g' for u in header['units']]
    names = [name for name, i in zip(header['channel_names'], i_acc) if i]

    Nt = data.shape[0]
    t = data[:,0]
    acc = data[:,i_acc]

    dt = t[1] - t[0]
    t = linspace(0, (Nt-1)*dt, Nt)

    return t, acc, header, names


def ver_datos(fname,
        Ncols=3,
        Nrows=6,
        save_this_fig=True,
        show_this_fig=False,
        save_this_input=False,
        f1=2.5,
        f2=125.,
        filter_order=20,
        tmin=-1,
        tmax=-1,
        filter_signal="acc",
        filter_type='sos',
        casename="igainput"
    ):

    print ("Working in: ", fname)

    header = leer_encabezado(fname)
    # data = genfromtxt(fname, skip_header=38)


    # t = data[:,0]
    # acc = data[:,i_acc]
    # t = linspace(0, (Nt-1)*dt, Nt)

    t, acc, header, names = get_datos(fname)
    Nt = len(t)
    i_acc =  [u == 'g' for u in header['units']]
    names = [name for name, i in zip(header['channel_names'], i_acc) if i]


    if tmin == -1:
        tmin = t[0]
    if tmax == -1:
        tmax = t[-1]

    dt = t[1] - t[0]

    difft = diff(t)
    Nchannels = acc.shape[1]
    fnyq = 0.5/dt
    Wn = [f1/fnyq, f2/fnyq]

    print ("  Date: ", header["date"])
    print ("  Time: ", header["time"])
    print ("    dt: ", dt)
    print ("  1/dt: ", 1/dt)
    print ("    Wn: ", Wn)
    print ("  t[0]: ", t[0])
    print (" t[-1]: ", t[-1])
    
    #Setup BA filter
    b, a = sig.butter(filter_order, Wn, btype='bandpass')
    filter_ba = lambda x : sig.filtfilt(b, a, x, axis=0,  method='pad', padlen=acc.shape[0]/2, padtype='even')
    
    #Setup SOS filter
    sos = sig.butter(filter_order, [f1, f2], btype='bandpass', fs=1/dt, output='sos')
    filter_sos = lambda x : sig.sosfilt(sos, x, axis=0)

    filter_none = lambda x : x

    if filter_type == "sos":
        used_filter = filter_sos
    elif filter_type == "ba":
        used_filter = filter_ba
    else:
        used_filter = filter_none
        filter_type = 'none'

    if filter_signal == 'acc':
        acc = used_filter(acc) 

    vel =  0*acc
    dis =  0*acc

    acc = 9.81*acc
    vel[1:,:] = cumsum(acc[1:,:] + acc[0:-1,:], axis=0)*(dt/2)
    dis[1:,:] = cumsum(vel[1:,:] + vel[0:-1,:], axis=0)*(dt/2)

    if filter_signal == 'dis':
        acc = used_filter(acc)
        vel = used_filter(vel)
        dis = used_filter(dis)

    acc = acc/9.81

    fig1, axs1 = subplots(Nrows,Ncols, sharex=True, sharey=True)
    fig1.set_size_inches([13,7], forward=True)

    fig2, axs2 = subplots(Nrows,Ncols, sharex=True, sharey=True)
    fig2.set_size_inches([13,7], forward=True)


    plotstuff = [
        (fig1, axs1, acc, "acc"),
        (fig2, axs2, dis, "dis"),
    ]

    yl = {}
    yl['acc'] = "Acc (g)"
    yl['dis'] = "Dis (m)"

    for fig, axs, signal, name  in plotstuff:
        amax = 0
        chmax = 0
        for i, j in iterprod(range(Nrows), range(Ncols)): 
            channel = i*Ncols + j

            if channel == Nchannels:
                break


            trange = logical_and(t > tmin, t < tmax)
            ax = axs[i][j]
            sca(ax)
            plot(t[trange], signal[trange, channel])
            if i < Nrows-1:
                ax.set_xticklabels([])
            else:
                xlabel("Time (s)")
                xtk_dt = (tmax-tmin)/10.
                xtk = arange(tmin, tmax+xtk_dt, xtk_dt)
                xtk_str = ["{0:4.1f}".format(s) for s in xtk]
                ax.set_xticks(xtk)
                ax.set_xticklabels(xtk_str)
            if j > 0:
                ax.set_yticklabels([])
            else:
                ylabel(yl[name])

            amax_this_ch = abs(signal[trange, channel]).max()
            if amax_this_ch > amax:
                chmax = channel
                amax = amax_this_ch

            print ("    Ch {0:4.0f}  {1} max = {2:4.3f}".format(channel, name, amax_this_ch))
            amax_1sd = round_to_1(amax)
            ylim([-amax_1sd,amax_1sd])
            xlim([tmin, tmax])
            grid(True)

        print ("{0} max = {1:4.3f}  at  channel # {2}".format(name, amax, chmax))

        yticks = [-amax_1sd,-amax_1sd/2,0,amax_1sd/2,amax_1sd]
        for i, j in iterprod(range(Nrows), range(Ncols)): 
            channel = i*Ncols + j            
            if channel == Nchannels:
                break
            axs[i][j].set_ylim([-amax_1sd,amax_1sd])
            axs[i][j].text(tmin,0.85*amax," {0:2.0f}: ".format(channel) + "{}".format(names[channel]), 
                verticalalignment="top", 
                bbox=dict(facecolor="white", alpha=0.3, edgecolor="white"))
        axs[i][0].set_yticks(yticks)
        axs[i][0].set_yticklabels(yticks)
        titlestring = fname.replace("./datos/","") + " Butter Filter on {}.  order = {} type = {} f1 = {} f2 = {}".format(filter_signal, filter_order, filter_type, f1, f2)
        suptitle(titlestring)
        tight_layout(rect=[0, 0.03, 1, 0.95])

        if save_this_fig:
            input_name = "./img/"+casename+"_"+fname.replace(".txt", "")
            savefig(input_name + "_{}.png".format(name), dpi=300)

        if save_this_input:
            input_name = "./iga_inputs/"+casename+"_"+fname.replace(".txt", "")
            print ("Saving ", input_name)
            with  open(input_name+"_y.txt", "w") as fid_y, \
                open(input_name+"_yp.txt", "w") as fid_yp, \
                open(input_name+"_ypp.txt", "w") as fid_ypp:

                trange = logical_and(t > tmin, t < tmax)

                # accel_to_read = 16 
                accel_to_read = 11

                y = dis[trange,accel_to_read]
                yp = vel[trange,accel_to_read]
                ypp = acc[trange,accel_to_read]

                Nt = len(y)

                fid_y.write("{} {}\n".format(Nt, dt))
                fid_yp.write("{} {}\n".format(Nt, dt))
                fid_ypp.write("{} {}\n".format(Nt, dt))

                savetxt(fid_y, y)
                savetxt(fid_yp, yp)
                savetxt(fid_ypp, ypp)
    if show_this_fig:
        show()
    else:
        close("all")



def main():
    # files = glob.glob("./datos/*.txt")
    files = glob.glob("./datos/190612_17_59_58.txt")


    for fname in files:
        ver_datos(fname,
            save_this_fig=True,
            show_this_fig=True,
            tmin=8.,
            tmax=16.,
            f1=2.,
            # filter_type='none',
            filter_type='sos',
            )
        break

if __name__ == "__main__":
    main()