import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

def pixDist(imglist):
    """ Plots the pixel distribution of the individual FLCs in ``img_list``
    for both chips.

    Parameters
    ----------
    imglist : list
        A list of images (FLCs or FLTs)

    Returns
    -------
    fig : obj
        The figure. You can run ``fig.show()`` on an iPython terminal to
        display the plots.
    """
    # ~~~~~~~~~~~~~~~~~~~~~~~ Initializations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fig, axs = plt.subplots(2, 1, figsize=(12,12), sharey=True)
    totalpix = 8400896 # total amt of pixels on 1 UVIS chip w/o overscan regions

    n = len(imglist)
    colors = plt.cm.jet(np.linspace(0., 0.5, n))
    alphas = np.arange(0.4, 1., n)

    
    for img, c, a in zip(imglist, colors, alphas):

        dat1 = fits.getdata(img, 1).flatten()
        dat2 = fits.getdata(img, 4).flatten()
        date = fits.getheader(img)['DATE-OBS']
        # ~~~~~~~~~~~~~~~~~~~~ Setting up Histogram ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # CHIP 1
        counts1, bedges1 = np.histogram(dat1, bins=100, range=(-10, 40))
        counts1 = (counts1/totalpix)*100
        bcenters1 = (bedges1[:-1] + bedges1[1:])/2.

        # CHIP 2
        counts2, bedges2 = np.histogram(dat2, bins=100, range=(-10, 40))
        counts2 = (counts2/totalpix)*100
        bcenters2 = (bedges2[:-1] + bedges2[1:])/2.

        # ~~~~~~~~~~~~~~~~~~~~~~~~~ Plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # CHIP 1
        axs[0].plot(bcenters1, counts1, '-', color=c, alpha=a, label=str(date))
        axs[0].set_ylabel('Amount of Pixels'+'\n'+'(% of Chip 1)')
        axs[0].set_xlabel(r'Pixel Value ($e^-$)')
        # adding the fill bw
        bcenters1 = np.asarray(bcenters1)
        counts1 = np.asarray(counts1)
        axs[0].fill_between(bcenters1[bcenters1>=13.5],counts1[bcenters1>=13.5], \
                            color='k', alpha=0.05)
        axs[0].set_title('Pixel Distribution for January FLCs'+'\n'+\
                        r'(Hot Pixels: $\geq 13.5 e^-$)')

        # ~aesthetics~
        axs[0].spines['right'].set_visible(False)
        axs[0].spines['top'].set_visible(False)
        axs[0].spines['left'].set_visible(False)
        axs[0].yaxis.grid(True)
        axs[0].legend()

        # CHIP 2
        axs[1].plot(bcenters2, counts2, '-', color=c, alpha=a, label=str(date))
        axs[1].set_ylabel('Amount of Pixels'+'\n'+ '(% of Chip 2)')
        axs[1].set_xlabel(r'Pixel Value ($e^-$)')
        # adding the fill bw
        bcenters2 = np.asarray(bcenters2)
        counts2 = np.asarray(counts2)
        axs[1].fill_between(bcenters2[bcenters2>=13.5],counts2[bcenters2>=13.5], \
                            color='k', alpha=0.05)

        # ~aesthetics~
        axs[1].spines['right'].set_visible(False)
        axs[1].spines['top'].set_visible(False)
        axs[1].spines['left'].set_visible(False)
        axs[1].yaxis.grid(True)
        axs[1].legend()

        # final adjustments
        plt.ylim(0)
        plt.subplots_adjust(hspace=0.2)



    plt.show()
    #return figs
