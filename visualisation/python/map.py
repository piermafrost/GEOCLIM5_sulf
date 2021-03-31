from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np

from matplotlib.backends.backend_pdf import PdfPages
pdf_pages=PdfPages('map.pdf')


def plot_z(lat, lon, z):
    """
        Plot data.
        
        Parameters
        ----------
        lat : latitude vector
        lon : longitude vector
        z : data array
        """
    #wrap the longitude vector around an arbitary point in order to fix non-plotting issue on date line
    wrapped_lon = np.concatenate([lon[5:],lon[:5]])
    wrapped_z = np.concatenate([z[:,5:],z[:,:5]],1)
    
    #create figure handle
    fig = plt.figure(figsize=(15,6))
    fig.suptitle('Regolith thickness')
    
    #create axis handle and add gridlines
    ax = plt.subplot(1,1,1, projection=ccrs.PlateCarree())
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                 color='gray', alpha=0.5, linestyle='--')
        
    # plot the data (plotted twice because of the date line bug)
    im = ax.pcolormesh(wrapped_lon, lat, wrapped_z, color='lightblue', transform=ccrs.PlateCarree())
    im = ax.pcolormesh(lon, lat, z, cmap='jet',transform=ccrs.PlateCarree())
    fig.colorbar(im, ax=ax)
    
    # show the plot
    # plt.show()
    pdf_pages.savefig(fig)




#dataset
mep = Dataset('/Users/yves/fortran/GEOCLIM4_ber/OUTPUT/dynsoil_output.ref.nc')

#map
plot_z(mep.variables['lat'][:], mep.variables['lon'][:], (mep.variables['runoff'][3,:,:])*100)
pdf_pages.close()
