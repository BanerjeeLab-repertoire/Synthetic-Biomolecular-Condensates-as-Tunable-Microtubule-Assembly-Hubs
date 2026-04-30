# Synthetic-Biomolecular-Condensates-as-Tunable-Microtubule-Assembly-Hubs
The repository contains the codes used for processing and analyzing the nanorheology data for the paper, "Synthetic Biomolecular Condensates as Tunable Microtubule Assembly Hubs".

# Software: Anaconda distribution v2024.10-1

# System requirements

Windows 10 or later, 64-bit macOS 10.15+ (for Intel) or 64-bit macOS 11.1+ (for Apple Silicon), or Linux, including Ubuntu, RedHat, CentOS 7+, and others.

# Installation guide

1. Download the installer from the Anaconda website
2. Go to your Downloads folder (or Home folder if downloaded via CLI) and double-click the installer to launch.
3. Click Next and then agree to Anaconda’s Terms of Service (TOS).
4. Select an installation option - Just me/ All users
5. Click Next.
6. Select a destination folder to install Anaconda, then click Next.
7. Click Install. The installation might take a few minutes to complete. Click Show details to view the packages being installed.
8. Click Next twice, then click Finish to close the installer.

# Demo

# For the determination of the mean square displacement (MSD) of the beads inside the condensate system, follow the steps below;
1. Launch the Jupyter notebook using Anaconda Navigator
2. In the Jupyter notebook, open the custom Python codes - "Video particle tracking nanorheology.ipynb"
3. Run the code by providing the .xml file as input. A file named "Bead trajectories after tracking using Fiji.xml" is provided in the folder - Demo/MSD-determination/. 

# Expected outcome
MSDs of the beads tracked in the condensates. The codes also allow you to estimate the viscosity of the condensate microenvironment.

# For estimation of the complex moduli from MSDs 
1. Open the custom Python codes - "Determination of complex moduli from MSDs using Evans method.ipynb" in the Jupyter notebook.
2. The .xml file provided in the folder - /Demo/Complex-moduli-determination should be used as an input by providing the address of the file in the code

# Expected outcome 
Complex moduli (Storage and loss moduli)

# Expected run time for demo on a "normal" desktop computer
For the entire demo (both MSD determination and estimation of complex moduli), it should take ~10 minutes. 

# Instructions for use
Videos acquired with known acquisition parameters (number of frames and exposure time) can be processed to get the particle trajectory data using Fiji. The custom Python codes provided here can be used on trajectory data with the input file format .xml to compute the particle MSDs and the complex moduli from the estimated MSDs.
