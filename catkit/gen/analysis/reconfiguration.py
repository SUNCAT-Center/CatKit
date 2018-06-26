from ..utils import running_mean
import numpy as np
import matplotlib.pyplot as plt


def id_reconstruction(images, save=False):
    """Identify a reconstruction event by analyzing changes in
    the forces.

    Parameters
    ----------
    images : list of ASE atoms objects
        Relaxation trajectory.
    show : bool
        Create a figure to display the events located.

    Returns
    -------
    predicted_events: list of int
        index of images predicted before the event occurs.
    """
    forces = []
    for i, atoms in enumerate(images):
        forces += [np.sqrt((atoms.get_forces()**2).sum())]
    forces = np.array(forces)

    frm = running_mean(forces)
    fdiff = np.diff(frm)
    fterm = np.array([fdiff > 0.25 * frm[:-1]]).astype(int)[0]
    predicted_events = np.where(fterm[:-1] < fterm[1:])[0]

    if save:
        fig, ax = plt.subplots(figsize=(6, 4))
        l, = plt.plot(range(1, len(images) + 1), frm)
        ax.fill_between(
            range(1,
                  len(images) + 1),
            np.zeros(len(images)),
            frm,
            facecolor=l.get_color(),
            alpha=0.5,
            interpolate=True)
        for i in predicted_events:
            plt.text(i - 1, 0.9, i)
            plt.axvline(i, ls='--', color='0.4')

        ylim = ax.get_ylim()
        plt.xlim(4, len(images))
        plt.ylim(0, ylim[1])
        plt.xlabel('Relaxation step')
        plt.ylabel('Force running mean (eV/$\AA$)')
        plt.savefig(save)
        plt.close()

    return predicted_events
