import argparse
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.wcs import WCS
from astropy.visualization import (ImageNormalize, AsymmetricPercentileInterval)
from ligo.skymap.plot.marker import sun
from ligo.skymap.postprocess import (find_greedy_credible_levels, interp_greedy_credible_levels)
from ligo.skymap.io.fits import read_sky_map
from ligo.skymap import (plot, moc)
from ligo.skymap.plot import mellinger
from reproject import reproject_interp
   
def _format_area(area):
    """Format area values for plot annotations.

    Values are printed to at least three significant figures. Values with three
    or more figures to the left of the decimal point are printed as integers.
    A comma separator is added for the thousands place.

    >>> _format_area(0.345678)
    '0.346'
    >>> _format_area(3.45678)
    '3.46'
    >>> _format_area(34.5678)
    '34.6'
    >>> _format_area(345.678)
    '346'
    >>> _format_area(3456.78)
    '3,457'
    >>> _format_area(34567.8)
    '34,568'
    """
    if area <= 100:
        return np.format_float_positional(
            area, precision=3, fractional=False, trim='-')
    else:
        return f'{np.round(area).astype(int):,d}'
        
def plot_skymap_with_ra_dec(skymapFILE,title,ra,dec,color):
    """ Process the given sky map file and extract relevant data. """
    
    # Read HEALPix sky map  
    hpx = read_sky_map(skymapFILE, moc=True) 
    #print(hpx)
    plt.figure(figsize=(10,6))
    # Retrieve axes and add marker (ra,dec)
    axes_args = {}
    axes_args['projection'] = 'astro mollweide'
    ax = plt.axes(**axes_args)
    ax.tick_params(colors='white')
    ax.grid()
    ax.plot_coord(SkyCoord(ra, dec, unit=u.deg), '*', markerfacecolor='white', markeredgecolor='black', markersize=10)
    plot.outline_text(ax)
   
    # Add annotations
    text = []
    text.append(f'Event: {title}')
    try:
       distmean = hpx.meta['distmean']
       diststd = hpx.meta['diststd']
    except KeyError:
        pass
    else:
        text.append(f'Distance: {round(distmean)}±{round(diststd)} Mpc')
    # Add contour on plot and annotations
    dA = moc.uniq2pixarea(hpx['UNIQ'])
    dP = hpx['PROBDENSITY'] * dA
    cls = 100 * find_greedy_credible_levels(dP, hpx['PROBDENSITY'])
    contour = np.array([50, 90])
    
    #Add in annotation
    i = np.flipud(np.argsort(hpx['PROBDENSITY']))
    areas = interp_greedy_credible_levels(contour, cls[i], np.cumsum(dA[i]), right=4*np.pi)
    pp = np.round(contour).astype(int)
    sr_to_deg2 = u.sr.to(u.deg**2)
    ii = areas * sr_to_deg2
    for i, p in zip(ii, pp):
        # FIXME: use Unicode symbol instead of TeX '$^2$'
        # because of broken fonts on Scientific Linux 7.
        text.append('{:d}% area: {} deg²'.format(p, _format_area(i)))
    ax.text(1, 1, '\n'.join(text), transform=ax.transAxes, ha='right')
    # Plot sky map
    center = SkyCoord(ra,dec,unit=u.deg)
    ax_inset = plt.axes([0.75, 0.2, 0.3, 0.3],projection='astro zoom',center=center,radius=10*u.deg)
    for key in ['ra', 'dec']:
        ax_inset.coords[key].set_ticklabel_visible(False)
        ax_inset.coords[key].set_ticks_visible(False)
    ax.mark_inset_axes(ax_inset, color=color)
    ax.connect_inset_axes(ax_inset, 'upper left').set_color(color)
    ax.connect_inset_axes(ax_inset, 'lower left').set_color(color)
    ax.connect_inset_axes(ax_inset, 'upper right').set_color(color)
    ax.connect_inset_axes(ax_inset, 'lower right').set_color(color)
    
    ax_inset.scalebar((0.1, 0.1), 5 * u.deg).label()
    ax_inset.compass(0.9, 0.1, 0.2)
    
    table = Table({'UNIQ': hpx['UNIQ'], 'CLS': cls})
    cs = ax_inset.contour_hpx((table, 'ICRS'), colors='white', linewidths=1,levels=contour, order='nearest-neighbor')
    fmt = r'%g\%%' if rcParams['text.usetex'] else '%g%%'
    plt.clabel(cs, fmt=fmt, fontsize=8, inline=True)
    hpx['PROBDENSITY'] *= 1 / sr_to_deg2
    img = ax.imshow_hpx((hpx, 'ICRS'), vmin=0, order='nearest-neighbor',  cmap='hot')
    img_inset = ax_inset.imshow_hpx((hpx, 'ICRS'), order='nearest-neighbor', cmap='hot')
    backdrop = mellinger()
    backdrop_wcs = WCS(backdrop.header).dropaxis(-1)
    interval = AsymmetricPercentileInterval(45, 98)
    norm = ImageNormalize(backdrop.data, interval)
    backdrop_reprojected = np.asarray([
        reproject_interp((layer, backdrop_wcs), ax.header)[0]
        for layer in norm(backdrop.data)])
    backdrop_reprojected = np.rollaxis(backdrop_reprojected, 0, 3)
    backdrop_reprojected = np.clip(backdrop_reprojected, 0, 1)
    img = ax.imshow(backdrop_reprojected,alpha=0.5)
    #img_inset = ax_inset.imshow(backdrop_reprojected,alpha=0.5)
    ax_inset.plot(center.ra.deg, center.dec.deg,transform=ax_inset.get_transform('world'), marker='*', markerfacecolor='white', markeredgecolor='white', markersize=10,markeredgewidth=1)
    cb = plot.colorbar(img_inset)
    cb.set_label(r'prob. per deg$^2$')
    plotname = f'plot_{title}'
    plt.savefig(plotname)
    plt.show()
    
    return plotname

    
def main():
    """ Main function to parse arguments and call processing function. """

    parser = argparse.ArgumentParser(description="Process a HEALPix sky map and extract relevant information.")
    parser.add_argument("skymapFILE", type=str, help="Path to the sky map FITS file")
    args = parser.parse_args()
    #print(args)
    title = "GW"
    plot_skymap_with_ra_dec(args.skymapFILE,title,47,-44,'grey')


if __name__ == "__main__":
    main()

