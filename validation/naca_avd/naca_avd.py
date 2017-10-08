import os
import glob
import numpy as np
import numpy.linalg as nla
import matplotlib.pyplot as plt
import ubem2d as ubem

this_dir = os.path.dirname(os.path.abspath(__file__))
if __name__ == '__main__':
    '''
    Validate the lift curve (CL vs. angle of attack) for various NACA airfoils
    against experimental data compiled in Abbot & von Doenhoff.
    '''
    pp = 0
    legend = []
    fnames = glob.glob(os.path.join(this_dir,'avd_naca_*_CL.txt'))
    for fname in fnames:
        print('Processing file: {}'.format(fname))
        code = fname.split('_')[-2]
        plt.figure()

        # Plot experimental lift data
        aoa_deg, CL_exp = ubem.read_data(fname).T
        plt.plot(aoa_deg, CL_exp, 'bo:', mfc='none')

        # Compute lift data at same angles of attack
        ubem_CL = np.zeros(len(aoa_deg))
        foil = ubem.naca4(code, 100)
        for i, aoa in enumerate(aoa_deg):
            uinf = (np.cos(aoa*np.pi/180), np.sin(aoa*np.pi/180))
            soln = ubem.solve_hess_smith_body(uinf, foil)
            ubem_CL[i] = ubem.airfoil_cdclcm(uinf, foil, soln.cp, pp)[1]
        plt.plot(aoa_deg, ubem_CL, 'ks-', mfc='none')

        # Compute lift at zero degrees and then thin-airfoil lift
        soln = ubem.solve_hess_smith_body((1,0), foil)
        CL_0 = ubem.airfoil_cdclcm(uinf, foil, soln.cp, pp)[1]
        plt.plot(aoa_deg, CL_0 + aoa_deg*np.pi**2/90, 'r^-.', mfc='none')
        plt.legend([r'NACA Experiment ($Re\sim 10^6$)', r'UBEM2D',
            r'Thin airfoil theory'])
        plt.grid(True)
        plt.xlabel(r'Angle of attack (deg)')
        plt.ylabel(r'Lift coefficient ($C_L$)')
        plt.title('NACA {}'.format(code))
        plt.savefig(os.path.join(this_dir, 'lift_curve_{}.pdf'.format(code)))

    # Show all plots and block
    print('Done')
    plt.show()
