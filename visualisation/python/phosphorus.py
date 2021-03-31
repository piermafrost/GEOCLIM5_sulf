import netCDF4 as nc
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from math import sqrt

TFACT = 1e-6
TUNITS = 'Ma'
FFACT = 1e-12
FUNITS = 'Tmol/a'

#flist = [nc.Dataset('../../OUTPUT/geoclim_output.ERA5_PI_equilibrium.nc')]
#f = nc.Dataset('../../OUTPUT/geoclim_output.ERA5_nocontOC.nc')
#f = nc.Dataset('../../OUTPUT/geoclim_output.ERA5_test_2xeros.nc')
#f = nc.Dataset('../../OUTPUT/geoclim_output.GFDL_PD_equilibrium.nc')
#flist = [nc.Dataset('../../OUTPUT/geoclim_output.GFDL_noIA_equilibrium.nc',
#         nc.Dataset('../../OUTPUT/geoclim_output.GFDL_noIA_equilibrium_2.nc')]
flist = [nc.Dataset('../../OUTPUT/geoclim_output.GFDL_abruptIA.nc'),
         nc.Dataset('../../OUTPUT/geoclim_output.GFDL_abruptIA_2.nc')]
#flist2 = []
flist2 = [nc.Dataset('../../OUTPUT/geoclim_output.GFDL_abruptIA_Pwth-noIA.nc'),
          nc.Dataset('../../OUTPUT/geoclim_output.GFDL_abruptIA_Pwth-noIA_2.nc')]

# Saving params:
FIGSIZE=(8,3)
ROOT='./'
RUNNAME = 'abruptIA_Pwth-noIA'


def figure(twinx=False):

    fig = plt.figure(figsize=FIGSIZE)

    SC = sqrt(FIGSIZE[0]*FIGSIZE[1]) / 5
    AR = FIGSIZE[0]/FIGSIZE[1]
    hmrg = 0.15/sqrt(AR)/SC
    vmrg = 0.11*sqrt(AR)/SC

    if twinx:
        ax = fig.add_axes([hmrg, vmrg, 1-2*hmrg, 1-vmrg-0.04])
        axb = ax.twinx()
        return fig, ax, axb

    else:
        ax = fig.add_axes([hmrg, vmrg, 1-hmrg-0.02, 1-vmrg-0.04])
        return fig, ax
        

pdf = PdfPages(ROOT+RUNNAME+'_phosphorus.pdf')


# oceanic PO4 and POP time evolution
fig, ax, axb = figure(twinx=True)
legid = []
for f in flist:
    pid, = ax.plot(TFACT*f['time'][:], 1e3*f['PO4_glob'][:], color='slateblue')
    axb.plot(TFACT*f['time'][:], 1e6*f['POP_glob'][:], color='tomato')
legid.append(pid)
for f in flist2:
    pid, = ax.plot(TFACT*f['time'][:], 1e3*f['PO4_glob'][:], '--', color='slateblue')
    axb.plot(TFACT*f['time'][:], 1e6*f['POP_glob'][:], '--', color='tomato')
legid.append(pid)
ax.set_xlabel('time ('+TUNITS+')')
ax.set_ylabel('oceanic PO$_4^{2-}$ (mmol/m$^3$)', color='slateblue')
axb.set_ylabel('ocn Part Org P ($\mu$mol/m$^3$)', color='tomato')
if flist2 != []:
    ax.legend(legid, ['regular', 'no Pwth from SEAI'])
pdf.savefig()

# phosphorus weathering and bioproductivity
fig, ax, axb = figure(twinx=True)
legid = []
for f in flist:
    pid, = ax.plot(TFACT*f['time'][:], 1e-9*f['P_wth_flux'][:], color='indigo')
    axb.plot(TFACT*f['time'][1:], FFACT*f['org_C_bio_prod'][1:,:].sum(1), color='teal')
legid.append(pid)
for f in flist2:
    pid, = ax.plot(TFACT*f['time'][:], 1e-9*f['P_wth_flux'][:], '--', color='indigo')
    axb.plot(TFACT*f['time'][1:], FFACT*f['org_C_bio_prod'][1:,:].sum(1), '--', color='teal')
legid.append(pid)
ax.set_xlabel('time ('+TUNITS+')')
ax.set_ylabel('phosphorus weathering (Gmol/a)', color='indigo')
axb.set_ylabel('ocn bioproductivity ('+FUNITS+')', color='teal')
if flist2 != []:
    ax.legend(legid, ['regular', 'no Pwth from SEAI'])
pdf.savefig()


pdf.close()
plt.show()
