import numpy as np
import matplotlib.pyplot as plt

__all__ = ['plot_wake', 'plot_circulation', 'plot_kinematics',
    'plot_aerodynamics']

def plot_wake(foil, wake, k=None, St=None, CT=None, CL=None, CM=None, eta=None,
    save_path=None):
    fig = plt.figure()
    plt.fill(foil.x/foil.chord, foil.y/foil.chord, color='.9')
    plt.plot(foil.x/foil.chord, foil.y/foil.chord, 'k-')

    plt.scatter(wake.x, wake.y, c=wake.gam, cmap='bwr', edgecolors='none')
    # Uncomment the following lines to plot the vortex cores as well
    #i = np.where(wake.gam >= 0)
    #plt.plot(wake.x[i]/foil.chord, wake.y[i]/foil.chord, 'ro', markersize=1)
    #i = np.where(wake.gam < 0)
    #plt.plot(wake.x[i]/foil.chord, wake.y[i]/foil.chord, 'bo', markersize=1)
    #cs, cx, cy = wake.vortex_cores()
    #plt.plot(cx, cy, 'ks', mfc='none', markersize=10)

    plt.axis('equal')
    plt.xlabel(r'$x/c$')
    plt.ylabel(r'$y/c$')
    plt.title('Wake structure')
    str_list = []
    if (k is not None):
        str_list.append('Reduced frequency: {:.3f}'.format(k))
    if (St is not None):
        str_list.append('Strouhal number: {:.3f}'.format(St))
    if (CT is not None):
        str_list.append('Thrust coefficient: {:.3f}'.format(CT))
    if (CL is not None):
        str_list.append('Lift coefficient: {:.3f}'.format(CL))
    if (CM is not None):
        str_list.append('Moment coefficient: {:.3f}'.format(CM))
    if (eta is not None):
        str_list.append('Efficiency: {:.1f}%'.format(100*eta))
    if (len(str_list) >= 1):
        plt.gca().text(.05, .7, '\n'.join(str_list), transform=\
            plt.gca().transAxes)
    if (save_path is not None):
        fig.savefig(save_path)

def plot_circulation(t, xlabel, wake, bcirc, save_path=None):
    fig = plt.figure()
    plt.subplot(2, 1, 1)
    i = np.where(wake.gam > 0)[0]
    plt.plot(t[i+1], wake.gam[i], 'r.', markersize=3)
    i = np.where(wake.gam < 0)[0]
    plt.plot(t[i+1], wake.gam[i], 'b.', markersize=3)
    plt.grid(True)
    plt.legend(['Pos. vortices','Neg. vortices'])
    plt.title('Circulation')
    plt.subplot(2, 1, 2)
    wcirc = np.concatenate([[0], np.cumsum(wake.gam)])
    plt.plot(t, bcirc)
    plt.plot(t, wcirc)
    plt.plot(t, bcirc + wcirc)
    plt.grid(True)
    plt.legend(['Bound', 'Wake', 'Total'])
    plt.xlabel(xlabel)
    if (save_path is not None):
        fig.savefig(save_path)

def plot_kinematics(t, xlabel, alp_deg=None, yoc=None, save_path=None):
    fig = plt.figure()
    n = sum(1 for x in [alp_deg, yoc] if x is not None)
    i = 1
    if (alp_deg is not None):
        plt.subplot(n,1,i)
        plt.plot(t, alp_deg)
        plt.grid(True)
        if (n > 1):
            plt.gca().get_yaxis().set_label_coords(-.1, .5)
            plt.ylabel('Pitch angle (deg)')
            plt.title('Kinematics')
        else:
            plt.title('Pitch angle (deg)')
        i += 1
    if (yoc is not None):
        plt.subplot(n,1,i)
        plt.plot(t, yoc)
        plt.grid(True)
        if (n > 1):
            plt.gca().get_yaxis().set_label_coords(-.1, .5)
            plt.ylabel('Heave/chord')
        else:
            plt.title('Heave/chord')
        i += 1
    plt.xlabel(xlabel)
    if (save_path is not None):
        fig.savefig(save_path)

def plot_aerodynamics(t, xlabel, CT=None, CL=None, CM=None, save_path=None,
    i0=0):
    fig = plt.figure()
    legend = []
    if (CT is not None):
        plt.plot(t[i0:], CT[i0:], 'g-')
        legend.append('CT')
    if (CL is not None):
        plt.plot(t[i0:], CL[i0:], 'b-')
        legend.append('CL')
    if (CM is not None):
        plt.plot(t[i0:], CM[i0:], 'r-')
        legend.append('CM')
    if (len(legend) == 1):
        plt.ylabel(legend[0])
        plt.title(legend[0])
    elif (len(legend) > 1):
        plt.legend(legend)
        plt.title('Aerodynamic coefficients')
    plt.xlabel(xlabel)
    plt.grid(True)
    if (save_path is not None):
        fig.savefig(save_path)
