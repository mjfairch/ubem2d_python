import os
import numpy as np
import matplotlib.pyplot as plt
import ubem2d as ubem

fs = ubem.FourierSeries(1.,[1.,1.,1.])
print('Base frequency: ', fs.base_frequency)
print('Period: ', fs.period)
print('Amplitudes: ', fs.amplitudes)
print('Phases: ', fs.phases)

plt.figure()
t = np.linspace(0,2*fs.period,100)
for i in range(len(fs)):
    plt.plot(t/fs.period, fs(t,range(i+1)),'.-')
plt.xlabel('t/T')
plt.ylabel('y(t/T)')
legend = ['1 mode'] + ['{} modes'.format(i+1) for i in range(1,len(fs))]
plt.legend(legend)
plt.title('Base frequency = {:.3f} Hz, period = {:3f} sec'.format(
    fs.base_frequency, fs.period))
plt.grid(True)
plt.savefig(os.path.join(ubem.__plot_dir, 'fourier_series.pdf'))
plt.show()
