import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys


###################
# GEOCLIM OUTPUT: #
###################

# Get argument pass to script

output_files = sys.argv
# remove 1st argument: script name
del(output_files[0])



######################
# DRAWING FUNCTIONS: #
######################


def plot_flux(ncdat, tindex=-1, figsize=(4,4), outfname=None, netflux=True):


    fig = plt.figure(figsize=figsize)
    ax = []
    # ax[0]: Carbon fluxes
    ax += [fig.add_axes([.04, .60, .92, .35])]
    # ax[1]: Carbonate fluxes
    ax += [fig.add_axes([.04, .32, .92, .15])]

    axtxt = fig.add_axes([.01, .01, .98, .27])
    axtxt.set_xlim([0, 1])
    axtxt.set_ylim([0, 1])
    axtxt.axis('off')

    def resizeheight(ax):
        h = ax.get_ylim()[1]
        ax.set_ylim([0, 1.03*h])

    def annotate_flux(ax, flxname, flxval, inout, ypos, fontsize=7, xpad=0.02):
        '''
        put name and flux value on horizontal barplot (all on a same axis).
        inout = +1 if name is wanted outside, -1 if wanted inside.
        xpad is relative to x width (i.e., =1 for entire x axis)
        both fontsize and xpad may be an array
        '''
        # prelim
        n = np.size(flxname)
        if np.size(fontsize)==1:
            fontsize = n*[fontsize]

        if np.size(xpad)==1:
            xpad = xpad*np.ones(n)

        xpad = xpad*np.diff(ax.get_xlim())

        # Annotation
        for name, val, io, y, fsz, xp in zip(flxname, flxval, inout, ypos, fontsize, xpad):
            absval = abs(val)
            if absval < 10:
                valtxt = '{:.2f}'.format(absval)
            elif absval < 100:
                valtxt = '{:.1f}'.format(absval)
            else:
                valtxt = '{:.0f}'.format(absval)

            orient = np.sign(val*io)
            if orient >= 0:
                ax.annotate(name,
                            (val+xp, y+0.17),
                            va='center', ha='left', fontsize=fsz)
                ax.annotate(valtxt,
                            (val+xp, y-0.17),
                            va='center', ha='left', fontsize=fsz)
            else:
                ax.annotate(name,
                            (val-xp, y+0.17),
                            va='center', ha='right', fontsize=fsz)
                ax.annotate(valtxt,
                            (val-xp, y-0.17),
                            va='center', ha='right', fontsize=fsz)



    # Carbon fluxes
    #--------------

    flxname = ['degass', 'sil wth', 'bas wth', 'ker wth', 'OC bur']
    flx = ['tot_CO2_degassing', 'sil_wth_C_flux', 'bas_wth_C_flux', 'ker_wth_C_flux', 'org_C_tot_dep_flux']
    col = ['crimson', 'steelblue', 'rebeccapurple', 'indianred', 'seagreen']
    y = [0, 0, 0, 1, 1]

    flxval = [1e-12*ncdat[fname][tindex] for fname in flx]

    # sources: +; sinks: -
    for k in [1, 2, 4]:
        flxval[k] = -flxval[k]

    #<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>#
    ax[0].barh(y, flxval, color=col, edgecolor='black', linewidth=0.5)
    #<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>#
    if netflux:
        net = [round(sum(flxval[:2]), 1), round(sum(flxval[3:]), 1)]
        ax[0].vlines(net[0], y[0]+0.4, y[0]-0.4, linestyles='--', colors='black', linewidths=1) 
        ax[0].vlines(net[1], y[3]+0.4, y[3]-0.4, linestyles='--', colors='black', linewidths=1) 
        annotate_flux(ax[0], flxname+['', ''], flxval+net, [-1, -1, +1, -1, -1, +1, +1], y+[y[0], y[3]])
    else:
        annotate_flux(ax[0], flxname, flxval, [-1, -1, +1, -1, -1], y)

    ax[0].set_xlabel('Tmol(C)/a', fontsize=8)


    # Carbonate weathering and precipitation
    #---------------------------------------

    flxname = ['+sw', 'carb wth', 'tot', 'carb prec (ner.)']
    flx = ['sil_wth_C_flux', 'carb_wth_C_flux', 'carb_pel_tot_dep_flux', 'carb_ner_tot_dep_flux']
    col = ['seashell', 'orchid', 'mediumblue', 'skyblue']
    y = len(col)*[0]

    flxval = [1e-12*ncdat[fname][tindex] for fname in flx]
    # sources: +; sinks: -
    for k in [-2, -1]:
        flxval[k] = -flxval[k]

    flxval[-2] = flxval[-2] + flxval[-1] # sum pelagic + neritic carbonate deposition
    flxval[0] += flxval[1] # sum carbonate and silicate weathering
    #<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>#
    ax[1].barh(y, flxval, color=col, edgecolor='black', linewidth=0.5)
    #<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>#
    if netflux:
        net = round(flxval[0]+flxval[2], 1)
        ax[1].vlines(net, y[0]+0.4, y[0]-0.4, linestyles='--', colors='black', linewidths=1) 
        annotate_flux(ax[1], flxname+[''], flxval+[net], [-1, -1, -1, -1], y+[y[0]])#,
        #            xpad=[0.005, 0.005, 0.02, 0.02, 0.02, 0.02],
        #            fontsize=[6, 6, 7, 7, 7, 7])
    else:
        annotate_flux(ax[1], flxname, flxval, [-1, -1, -1, -1], y)#,
        #            xpad=[0.005, 0.005, 0.02, 0.02, 0.02],
        #            fontsize=[6, 6, 7, 7, 7])

    ax[1].set_xlabel('Tmol(Ca)/a', fontsize=8)


    # Customization
    #--------------

    for axi in ax:
        ylim = axi.get_ylim()
        axi.vlines(0, *ylim, colors='grey', linewidths=0.5)
        axi.set_ylim(ylim)
        axi.tick_params(axis='x', labelsize=7)
        axi.get_yaxis().set_visible(False)
        for wch in ['top', 'left', 'right']:
            axi.spines[wch].set_visible(False)


    # water discharge, sediment discharge, sedimentation flux
    #--------------------------------------------------------

    axtxt.annotate('Disch: {:.3f} Sv'.format((1e-6/(365.2422*24*60*60))*ncdat['discharge'][tindex]),
                   (.05, .75), ha='left', va='top')
    axtxt.annotate('TSS: {:.3f} Gt/yr'.format(1e-12*ncdat['TSS'][tindex]),
                   (.05, .50), ha='left', va='top')
    axtxt.annotate('sedim: {:.3f} Gt/yr'.format(1e-12*ncdat['sedim_flux'][tindex,:].sum()),
                   (.05, .25), ha='left', va='top')


    # CO2 and O2
    #-----------

    axtxt.annotate('bioprod: {:.1f} Tmol/yr'.format(1e-12*ncdat['org_C_bio_prod'][tindex,:].sum()),
                   (.55, .75), ha='left', va='top')
    axtxt.annotate('pCO2: {:.3f} PAL'.format(ncdat['CO2_atm_level'][tindex]),
                   (.55, .50), ha='left', va='top')
    axtxt.annotate('pO2: {:.3f} PAL'.format(ncdat['O2_atm_level'][tindex]),
                   (.55, .25), ha='left', va='top')


    # Printing/Saving
    #----------------

    if outfname is not None:
        fig.savefig(outfname+'.pdf', format='pdf', transparent=True)




############################################################################
#--------------------------------------------------------------------------#
############################################################################



def plot_ocechem(ncdat, tindex=-1, figsize=(7,4), outfname=None):

    colors = ['navy', 'firebrick', 'forestgreen']

    if outfname is not None:
        pdf = PdfPages(outfname+'.pdf')
    else:
        # Trick: create an object with attributes `savefig' and `close': methods that do nothing
        class dummypdfpages():
            def savefig(self):
                pass
            close = savefig

        pdf = dummypdfpages()



    ########################################################
    # Oceanic profiles (polar, Mid-lat and epicontinental) #
    ########################################################

    for specie in ['alkalinity', 'DIC', ['H2CO3', 'HCO3', 'CO3'], 'pH', 'temperature', 'Ca', 'POC', 'O2', 'PO4', 'POP', 'lysoc_depth_carb', 'lysoc_depth_arag', 'DIC_d13C', 'Sr', 'Sr_iso_ratio']:

        if type(specie) is not list:
            listspe = [specie]
        else:
            listspe = specie

        fig = plt.figure(figsize=figsize)
        ax1 = fig.add_axes([.10, .3, .26, .67])
        ax2 = fig.add_axes([.41, .3, .26, .67], sharex=ax1)
        ax3 = fig.add_axes([.72, .3, .26, .67], sharex=ax1)
        axtxt = fig.add_axes([.01, .01, .98, .12])
        axtxt.set_xlim([0, 1])
        axtxt.set_ylim([0, 1])
        axtxt.axis('off')

        for k,spe in enumerate(listspe):
            units = ncdat[spe].units
            if 'd13C' in spe and units == '-':
                fact = 1e3
                units = '‰'
            elif spe in ['PO4', 'POP', 'Sr']:
                fact = 1e3
                units = 'm'+units
            else:
                fact = 1

            ax1.plot(fact*ncdat[spe][tindex,0:2], [0,-2],    '-', marker='o', ms=5, color=colors[k])
            ax1.plot(fact*ncdat[spe][tindex,7:9], [0,-2],   '--', marker='o', ms=5, color=colors[k])

            ax2.plot(fact*ncdat[spe][tindex,2:5], [0,-1,-2], '-', marker='o', ms=5, color=colors[k])

            ax3.plot(fact*ncdat[spe][tindex,5:7], [0,-1],    '-', marker='o', ms=5, color=colors[k])

        for ax in [ax1, ax2, ax3]:
            ax.set_ylim([-2.2, 0.2])
            ax.set_yticks([-2,-1,0])

        ax1.set_yticklabels(['deep', 'thermo', 'surf'], rotation=40)
        for ax in [ax2, ax3]:
            ax.set_yticklabels([])

        for ax in [ax1, ax2, ax3]:
            ax.tick_params(axis='x', labelsize=7)

        if listspe[0] == 'H2CO3':
            ax1.set_xscale('log')

        ax1.set_xlabel('POLAR')
        ax2.set_xlabel('MID-LAT')
        ax3.set_xlabel('EPICONT')

        ax1.legend(['N', 'S'])

        if len(listspe) > 1:
            ax3.legend(listspe)


        # Mean value:
        #------------

        if specie == 'alkalinity':
            globname = 'alk_glob'
        elif listspe[0] in ['lysoc_depth_carb', 'lysoc_depth_arag', 'H2CO3']:
            globname = None
        else:
            globname = specie+'_glob'

        name = listspe[0]
        for spe in listspe[1:]:
            name += ', '+spe

        axtxt.annotate(name+' ('+units+')', (0.1, 0.5), fontsize=12, fontweight='bold', ha='left', va='bottom')
        if globname is not None:
            meanval = fact*ncdat[globname][tindex]

            axtxt.annotate('Mean value: {:.5e}'.format(meanval), (0.5, 0.5), ha='left', va='bottom')


        # Save page
        #----------

        pdf.savefig()



    ###########################
    # Sedimentation variables #
    ###########################

    color = ['goldenrod', 'maroon', 'navy', 'turquoise', 'plum']
    color_list = ['maroon', 'goldenrod', 'slateblue']
    index = [5, 6, 4, 1, 8]
    legnd = ['shelf', 'ep deep', 'open ML', 'open NP', 'open SP']
    xpos  = [1, 3, 5, 6.5, 8]

    for var in ['sedim_rate', 'sedim_flux', 'org_C_dep_flux', 'burial_efficiency', 'P_dep']:

        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([.09, .3, .88, .65])
        axtxt = fig.add_axes([.01, .01, .98, .11])
        axtxt.set_xlim([0, 1])
        axtxt.set_ylim([0, 1])
        axtxt.axis('off')


        if var == 'P_dep':
            value = [1e-9*f[v][tindex,index] for v in ['P_dep_flux_orgC', 'P_dep_flux_phosph', 'P_dep_flux_hydro']]
            value[1] += value[2]
            value[0] += value[1]
            unit = 'Gmol/a'
            name = 'Phosphorus burial fluxes'
            legend = ['orgC-bound', 'phosphorite', 'Hydro Fe-bound']
        elif var == 'sedim_rate':
            value = 1e3 * f[var][tindex,index]
            unit = 'mm/a'
            name = f[var].long_name
        elif var in ['sedim_flux', 'org_C_dep_flux']:
            value = 100 * f[var][tindex,index] / f[var][tindex,:].sum()
            unit = '% of total'
            name = f[var].long_name
        else:
            value = f[var][tindex,index]
            unit = f[var].units
            name = f[var].long_name



        # Bar plot
        #---------

        if type(value) is list:
            for val,col in zip(value, color_list):
                ax.bar(xpos, val, width=1, color=col, edgecolor='black', linewidth=0.5)

            ax.legend(legend)
        else:
            ax.bar(xpos, value, width=1, color=color, edgecolor='black', linewidth=0.5)

        ax.set_xticks(xpos)
        ax.set_xticklabels(legnd, rotation=30, ha='right')
        ax.set_ylabel(unit)

        axtxt.annotate(name, (0.1, 0.5), fontsize=12, fontweight='bold', ha='left', va='bottom')


        # Save page
        #----------

        pdf.savefig()


    # Close pdf file
    ################

    pdf.close()


############################################################################


for name in output_files:

    f = nc.Dataset(name)

    # remove "geoclim_output."
    k = name.find('geoclim_output')
    name = name[k+15:]
    # remove ".nc"
    name = name[:-3]

    #it = int(np.argwhere(f['time'][:]==700000))
    it = -1

    # Fluxes histograms
    plot_flux(f, tindex=it, outfname='final_fluxes--'+name)

    # Ocean chemistry
    plot_ocechem(f, tindex=it, outfname='final_ocean_chemistry--'+name)

    f.close()
