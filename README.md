# LT-SPRAT-pipeline
A quick pipeline for manual reductions of SPRAT spectra

Requires Pyraf, see Astroconda

To run, enter the pyraf conda environment of your choice, navigate to the folder with the files and type 'python manualextraction.py'

Download both the arc and the the spectrum files from the LT repository. Rename arc file to 'arc.fits'

Pipeline will ask to extract arc first. Follow standard iraf commands, usually lots of 'y', occasionally type 'q' to quit the window

On line identification. Press 'n' . (may be 'm'...) in the iraf window over the strongest line, enter 7642 and press enter
Press 'n' over the second strongest line, type 4671 and press enter. Press 'l' to automatically find the lines.

Next, extract the spectrum as per standard IRAF procedures

The rest of the pipeline should run automatically.

The final flux calibration applys a correction to the spectrum using the standard star used in the sensitivity function. This compensates for some of the calibration and throughput issues with SPRAT. Another consequence that it applies a telluric correction. Future updates will see corrections for different airmasses.
