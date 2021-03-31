import netCDF4 as nc
from matplotlib import pyplot as plt
from math import sqrt


############################################################################

TFACT = 1e-6
TUNITS = 'Ma'
FFACT = 1e-12
FUNITS = 'Tmol/a'

# List of run to plot. Will be plotted continuously, as if it was a single run
#flist = [nc.Dataset('../../OUTPUT/geoclim_output.test.nc')]
flist = [nc.Dataset('../../OUTPUT/geoclim_output.test-ERA5.bis.nc')]
# Can be used to compare with a second list of runs (dashed)
flist2 = []
#flist = [nc.Dataset('../../OUTPUT/geoclim_output.GFDL_noIA_equilibrium.nc',
#         nc.Dataset('../../OUTPUT/geoclim_output.GFDL_noIA_equilibrium_2.nc')]
#flist = [nc.Dataset('../../OUTPUT/geoclim_output.GFDL_abruptIA.nc'),
#         nc.Dataset('../../OUTPUT/geoclim_output.GFDL_abruptIA_2.nc')]
#flist2 = [nc.Dataset('../../OUTPUT/geoclim_output.GFDL_abruptIA_Pwth-noIA.nc'),
#          nc.Dataset('../../OUTPUT/geoclim_output.GFDL_abruptIA_Pwth-noIA_2.nc')]

# Saving params:
SAVE = False
FIGSIZE=(8,2.5)
ROOT='./'
RUNNAME = 'ERA5_PI_eq'

############################################################################


def figure(twinx=False):

    fig = plt.figure(figsize=FIGSIZE)

    SC = sqrt(FIGSIZE[0]*FIGSIZE[1]) / 5
    AR = FIGSIZE[0]/FIGSIZE[1]
    hmrg = 0.15/sqrt(AR)/SC
    vmrg = 0.11*sqrt(AR)/SC

    if twinx:
        ax = fig.add_axes([hmrg, vmrg, 1-2*hmrg, 1-vmrg-0.02])
        axb = ax.twinx()
        return fig, ax, axb

    else:
        ax = fig.add_axes([hmrg, vmrg, 1-hmrg-0.02, 1-vmrg-0.02])
        return fig, ax
        


# CO2 and O2 time evolution
fig, ax = figure()
for f in flist:
    ax.plot(TFACT*f['time'][:], f['CO2_atm_level'], color='slateblue')
    ax.plot(TFACT*f['time'][:], f['O2_atm_level'], color='tomato')
for f in flist2:
    ax.plot(TFACT*f['time'][:], f['CO2_atm_level'], '--', color='slateblue')
    ax.plot(TFACT*f['time'][:], f['O2_atm_level'], '--', color='tomato')
ax.legend(['CO2', 'O2'])
ax.set_xlabel('time ('+TUNITS+')')
ax.set_ylabel('atmospheric level (PAL)')
if SAVE:
    fig.savefig(ROOT+RUNNAME+'_CO2-O2.pdf', format='pdf', transparent=True)

# main carbon fluxes
fig, ax = figure()
for f in flist:
    ax.plot(TFACT*f['time'][:], FFACT*f['tot_CO2_degassing'][:], color='crimson')
    ax.plot(TFACT*f['time'][:], FFACT*(f['gran_wth_C_flux'][:]+f['bas_wth_C_flux'][:]), color='indigo')
    ax.plot(TFACT*f['time'][:], FFACT*f['ker_wth_C_flux'][:], color='goldenrod')
    ax.plot(TFACT*f['time'][:], FFACT*f['org_C_tot_dep_flux'][:], color='forestgreen')
for f in flist2:
    ax.plot(TFACT*f['time'][:], FFACT*f['tot_CO2_degassing'][:], '--', color='crimson')
    ax.plot(TFACT*f['time'][:], FFACT*(f['gran_wth_C_flux'][:]+f['bas_wth_C_flux'][:]), '--', color='indigo')
    ax.plot(TFACT*f['time'][:], FFACT*f['ker_wth_C_flux'][:], '--', color='goldenrod')
    ax.plot(TFACT*f['time'][:], FFACT*f['org_C_tot_dep_flux'][:], '--', color='forestgreen')
ax.legend(['degass', 'sil wth', 'ker wth', 'OC bur'])
ax.set_xlabel('time ('+TUNITS+')')
ax.set_ylabel('fluxes ('+FUNITS+')')
if SAVE:
    fig.savefig(ROOT+RUNNAME+'_C-fluxes.pdf', format='pdf', transparent=True)

# carbonates fluxes
fig, ax = figure()
for f in flist:
    ax.plot(TFACT*f['time'][:], FFACT*f['carb_wth_C_flux'][:], color='magenta')
    ax.plot(TFACT*f['time'][:], FFACT*(f['carb_wth_C_flux'][:]+f['gran_wth_C_flux'][:]+f['bas_wth_C_flux'][:]), color='darkmagenta')
    ax.plot(TFACT*f['time'][:], FFACT*(f['carb_ner_tot_dep_flux'][:]+f['carb_pel_tot_dep_flux'][:]), color='slateblue')
for f in flist2:
    ax.plot(TFACT*f['time'][:], FFACT*f['carb_wth_C_flux'][:], '--', color='magenta')
    ax.plot(TFACT*f['time'][:], FFACT*(f['carb_wth_C_flux'][:]+f['gran_wth_C_flux'][:]+f['bas_wth_C_flux'][:]), '--', color='darkmagenta')
    ax.plot(TFACT*f['time'][:], FFACT*(f['carb_ner_tot_dep_flux'][:]+f['carb_pel_tot_dep_flux'][:]), '--', color='slateblue')
ax.legend(['carb wth', 'tot wth', 'carb precip'])
ax.set_xlabel('time ('+TUNITS+')')
ax.set_ylabel('fluxes ('+FUNITS+')')
if SAVE:
    fig.savefig(ROOT+RUNNAME+'_C+Carb-fluxes.pdf', format='pdf', transparent=True)

# OC production (continental and oceanic)
fig, ax, axb = figure(twinx=True)
for f in flist:
    ax.plot(TFACT*f['time'][1:], FFACT*f['total_cont_POC_export'][1:], color='indigo')
    axb.plot(TFACT*f['time'][1:], FFACT*f['org_C_bio_prod'][1:,:].sum(1), color='teal')
for f in flist2:
    ax.plot(TFACT*f['time'][1:], FFACT*f['total_cont_POC_export'][1:], '--', color='indigo')
    axb.plot(TFACT*f['time'][1:], FFACT*f['org_C_bio_prod'][1:,:].sum(1), '--', color='teal')
ax.set_xlabel('time ('+TUNITS+')')
ax.set_ylabel('cont OC export flux ('+FUNITS+')', color='indigo')
axb.set_ylabel('ocn bioproductivity ('+FUNITS+')', color='teal')
if SAVE:
    fig.savefig(ROOT+RUNNAME+'_OC-prod-flux.pdf', format='pdf', transparent=True)

# phosphorus weathering
fig, ax = figure()
for f in flist:
    ax.plot(TFACT*f['time'][:], 1e-9*f['P_wth_flux'][:], color='indigo')
for f in flist2:
    ax.plot(TFACT*f['time'][:], 1e-9*f['P_wth_flux'][:], '--', color='indigo')
ax.set_xlabel('time ('+TUNITS+')')
ax.set_ylabel('phosphorus weathering (Gmol/a)')
if SAVE:
    fig.savefig(ROOT+RUNNAME+'_P-wth.pdf', format='pdf', transparent=True)

# erosion
fig, ax = figure()
for f in flist:
    ax.plot(TFACT*f['time'][:], 1e-12*f['TSS'][:], color='indigo')
for f in flist2:
    ax.plot(TFACT*f['time'][:], 1e-12*f['TSS'][:], '--', color='indigo')
ax.set_xlabel('time ('+TUNITS+')')
ax.set_ylabel('TSS (Gt/a)')
if SAVE:
    fig.savefig(ROOT+RUNNAME+'_erosion.pdf', format='pdf', transparent=True)


if not SAVE:
    plt.show()
