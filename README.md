# oycxms
miscellaneous scripts/codes useful for MS analysis / ZTGROUP

## CWTPP - new continuous wavelet transform (CWT)-based preprocessing & visualisation workflow for ambient 1D/2D (imaging) MS data

## Instructions for installation
1. Install Matlab, recommended versions are anything >R2020a (preferably not the latest beta version)

1.	Install Windows SDK7 for Matlab as per instructions here:
https://uk.mathworks.com/matlabcentral/answers/101105-how-do-i-install-microsoft-windows-sdk-7-1
*Note that the second answer (not the ‘accepted’) by Andre Silva is the more useful in most cases, in summary:
•	Remove programs > uninstall Visual studios C++ (all of them)
•	Remove programs > uninstall NET. (all of them)
•	Restart PC
•	Install the SDK (**winsdk_web.exe**)
press ok on the warning, make sure you untick "Visual C++ Compilers" and "Microsoft Visual C++ 2010" components.
•	Install 
NDP452-KB2901907-x86-x64-AllOS-ENU.(**.net installation**)* _if the installation states that you already have this installed (most likely for modern Windows), then you can skip this_
VC-Compiler-KB2519277.exe (SDK patches)
vcredist_x86.exe (c++ compiler package)
vcredist_x64.exe (c++ compiler package)
•	Open Matlab, try typing ‘mex -setup’ in the command window, if the installation was completed correctly then you will see this:

'MEX configured to use '**Microsoft Windows SDK 7.1(C)**' for C language compilation.'

•	If you did not see this, then you many have to copy some configuration files (xml files, **please email me**) to put in C: Users > your username > AppData (you may need to reveal this hidden folder) > Roaming > MathWorks > MATLAB > your MATLAB version >

3. You will also need the following toolboxes:
'Signal Processing Toolbox'	'8.4'
'Statistics and Machine Learning Toolbox'	'11.7'
'Wavelet Toolbox'	'5.4'
'Parallel Computing Toolbox'	'7.2'

You can find & intall these within Matlab by HOME > Add-Ons > Get Add-Ons


11/04/2023 v2 release  notes:

- all known bugs fixed
- added visualisation functions (TIC_image, RGB_image etc. refer to wiki (in progress))
- now compatible with .imzml format raw files as well
- inter-/intra-data alignment functions can now be run independently and save intermediate results
- option to take metadata input coming soon
- recalibration functions (with and without metadata input) coming soon


## spectral processing 

- v2 update 05/10/2022, most cells functionalised & added enhanced graphic output utilities; worked examples included for most cells
- contains relevant scripts for pre-/post-processing of 2D spectral data
- for generalised script for the unsupervised/supervised analysis of MS spectral data, download both 'classifiers_n_features.ipynb' with functions included in 'classifiers_functions.ipynb' under 'MS classification'
- v1 update 11/02/2022, wrapped installation file/worked example to come

## GUI - the graphical user interface element for the old py_DESI_MSI workflow for imaging data analysis (REDUNDANT and UNSERVICED)
