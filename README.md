# About

A pipeline for the reduction of Liverpool Telescope SPRAT spectra using IRAF.

Requires ```Pyraf```, see Astroconda.

# Installation

To run, the code just needs to be cloned or downloaded.

Download the required spectrun and arc ```.fits``` files from the LT data archive.

# Commands

Fom the command line, navigate to the pipeline folder and  run

```python SPRATManualExtraction.py```

The code will automatically identify the spectrum file and the arc file.

### Arc line identification

Follow standard ```IRAF``` commands to move through the windows. 

The lines to be identified, and an example of their positions can be found at http://telescope.livjm.ac.uk/TelInst/Inst/SPRAT/sprat_xenon_arc.png

When the line identification window is reached, press "m" in the ```IRAF``` window over the strongest line, enter "7642" and press enter Press "m" over the second strongest line, type "4671" and press enter. Press "l" to automatically find the lines.

_Note: keep an eye on the terminal window output, if it says "no monotonic solution found", try again and extract the arc file in a different position. If the line fit is obviously incorrect, try again and add more lines from. If this still fails, copy the database and arc.ms from a successful run._

### Spectrum extraction

Follow standard ```IRAF``` commands to move through the windows. 

When reaching the spectrum extraction window, for bright targets the object will be correctly identified. In some cases, another source is selected for extraction. To delete the aperture press "d" over the designated aperture and then press "n" over the correct location (approximately 126).

To manually set the background if required -- press "b", then "z" to delete the background regions and "s" to set them. Finish by pressing "f" to fit the background and "q" to exit.

### Flux calibration

The spectrum undergoes an initial flux calibration using ```newsen.fits```. The final flux calibration applies a correction to the spectrum using the standard star used in ```newsen.fits```. This compensates for some of the calibration and throughput issues with SPRAT. Another consequence that it applies a telluric correction.

# Output

The output spectrum is ```<name>_<observation date>.w.txt```