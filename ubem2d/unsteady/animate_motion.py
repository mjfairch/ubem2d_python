import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

__all__ = ['animate_motion']

def animate_motion(foil, pp, motion, fps=50, save_path = None, save_count=100,
    bitrate=57600):
    '''
    This is a helper function to animate a given motion.  The resulting
    animation may be either displayed to the screen or saved to the
    filesystem (but not both).
    '''
    def frames(motion):
        t0, alp0, y0 = 0, 0, 0
        for (t, (alp, y)) in motion:
            dt, dalp, dy = t-t0, alp-alp0, (y-y0)*foil.chord
            t0, alp0, y0 = t, alp, y
            yield dt, (dalp, dy)

    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-.5*foil.chord, 1.5*foil.chord)
    ax.set_ylim(-foil.chord, foil.chord)
    ax.set_axis_off()
    foil_line, = ax.plot([], [], 'k-')

    def draw(frame):
        dt, (dalp, dy) = frame
        foil.pitch(dalp*np.pi/180, pp)
        foil.heave(dy)
        foil_line.set_data(foil.x, foil.y)
        return foil_line,

    ani = animation.FuncAnimation(fig, draw, frames(motion), 
        interval=round(1000/fps), save_count=save_count, repeat=False,
        blit=True)
    if (save_path is None):
        plt.show()
    else:
        writer = animation.writers['ffmpeg'](fps=fps, bitrate=bitrate)
        ani.save(save_path, writer=writer)
