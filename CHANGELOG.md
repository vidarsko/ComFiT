# Changelog

All notable changes to this project will be documented in this file.
Full documentation can be found here: [https://comfitlib.com/](https://comfitlib.com/).

## [1.9.5] - 2025-04-16
- Bug fix for tool_configure_axis function.

## [1.9.4] - 2025-04-16
- Removed the need to pass both `ax` and `fig` to plot functions.

## [1.9.3] - 2025-04-16
- Changed definition of `xlim`. It is now the limit of the x-axis, not the size `xmax` of the domain. See comfitlib.com for more details.
- Added shadows to 3D complex plots with plotly.
- Changed default plotting setting for fourier space from 2pi/a0 to 1/a0. 
- Bug fixes.

## [1.9.2] - 2025-03-25
- Fixed bug related to trying to plot a constant field.

## [1.9.1] - 2025-03-24
- Bug fix: setting ylim with plotly now works for 2D plots. 
- Added `width` and `height` parameters to plot_save-function.
- Added phase angle plot to `plotly` library.
- Changing the requirements to work with numpy>=2.2.4

## [1.9.0] - 2025-03-18
- Added the option to pass `fourier=True` to the plot functions to plot a field that is in Fourier space.

## [1.8.7] - 2025-03-13
- Changed default colormap to 'angle' for plot_complex_field_in_plane function.

## [1.8.6] - 2025-03-13
- Fixed bug in plot_save-function.

## [1.8.5] - 2025-03-11
- Added 'opacity' to 3D plot field functions.

## [1.8.4] - 2025-03-03
- Changed plotting convention cfi.plot_save(n,fig) -> cfi.plot_save(fig,n)

## [1.8.3] - 2025-02-28
- Changed calculation of integrating factors in base system class remove division by zero errors.

## [1.8.2] - 2025-02-21
- Fixed plotly colorbar placement.

## [1.8.1] - 2025-02-21
- Plotly now returns a fig, ax object, like the matplotlib functions.

## [1.8.0] - 2025-02-04
- Created more seamless transition between plotly and matplotlib. 

## [1.7.0] - 2025-01-08
- Changed default plotting library to plotly
- Changes to PFC class

## [1.6.2] - 2024-09-07
- Bug fixes for plotly vector plots.

## [1.6.1] - 2024-09-05
- Added more plotly features for alternative plotting.

## [1.6.0] - 2024-06-30
- Added plotly feature integration for alternative plotting (under development).
- Added polycrystalline methods for the PFC models.

## [1.5.0] - 2024-05-30
- Finalized version of software after JOSS software review. 

## [1.4.2] - 2024-03-18
- Bug fixes and stability improvements.

## [1.4.1] - 2024-03-15
- Fixed error in plot_complex_field function.

## [1.4.0] - 2024-03-15
- Completed writing all the plot functions. 
- Added stress tensor calculations for all the PFC models.
- Completed much of the documentation.

## [1.3.1] - 2024-02-24
- Added stress divergence calculation method to all PFC models.

## [1.3.0] - 2024-02-24
- All mayavi functionality removed and placed in optional extension (ComFiT-mayavi)

## [1.2.2] - 2024-02-24
- Relaxed vtk requirement

## [1.2.1] - 2024-02-24
- Made a wheel for the distribution

## [1.2.0] - 2024-02-24
- Bug fixes and stability imrpovements.
- Included the Mayavi package for 3D plotting.

## [1.1.0] - 2023-12-07
### Fixed
- Many bug fixes and stability improvements, in particular for the BoseEinsteinCondensate class.

## [1.0.0] - 2023-11-22
- Initial release of the package.
- Containing most functions, but lacking in stability.