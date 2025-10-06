# oycxms
miscellaneous scripts/codes useful for MS analysis / ZTGROUP

for external access, please fill in some basic information here:
https://uniregensburg-my.sharepoint.com/:x:/g/personal/xiy23208_ads_uni-regensburg_de/EaBSLNacqtFMkqUZjDYqyFgBRRM2Yoe36JffvONkpLTxPA?e=oWm5SP

and email me on 
yuchen.xiang11@imperial.ac.uk
yuchen.xiang@ur.de

## CWTPP - new continuous wavelet transform (CWT)-based preprocessing & visualisation workflow for ambient 1D/2D (imaging) MS data

## Instructions for installation Oct 2025 - YOU NO LONGER NEED THE WINDOWS SDK 7
1. Install Matlab, recommended versions are anything >R2020a (preferably not the latest beta version)

2. It is highly recommended that you first convert any raw files into the more efficient, faster (~40X for preprocessing) format MZ5 through the use of ProteoWizard (MSconvertGUI):

https://proteowizard.sourceforge.io/download.html
select 'mz5' format and use all default parameters.

(Should you for some reason want to work with .raw files, you will stil need to install Windows SDK7 for Matlab as per instructions here:)

https://uk.mathworks.com/matlabcentral/answers/101105-how-do-i-install-microsoft-windows-sdk-7-1 (**no longer available, please download from link below**)

[https://www.dropbox.com/scl/fo/p27kwayp9mrz3tvqvnh1i/AI3-P76BoJus89o8-6hPlMg?rlkey=0fxxkjjqwoa7uxufmzhn1bhrr&st=kx5rbzsc&dl=0](https://www.dropbox.com/scl/fo/p27kwayp9mrz3tvqvnh1i/AI3-P76BoJus89o8-6hPlMg?rlkey=0fxxkjjqwoa7uxufmzhn1bhrr&st=mxznmi56&dl=0)

3. You will also need at least the following toolboxes:
'Signal Processing Toolbox'	'8.4'
'Statistics and Machine Learning Toolbox'	'11.7'
'Wavelet Toolbox'	'5.4'
'Parallel Computing Toolbox'	'7.2'

You can find & intall these within Matlab by HOME > Add-Ons > Get Add-Ons


06/10/2025 update notes:
- compatibility with mz5 added
- parameter tuning ongoing, optimised parameters+validation/test data coming soon
- detailed documentation coming soon



## spectral processing 

- v2 update 05/10/2022, most cells functionalised & added enhanced graphic output utilities; worked examples included for most cells
- contains relevant scripts for pre-/post-processing of 2D spectral data
- for generalised script for the unsupervised/supervised analysis of MS spectral data, download both 'classifiers_n_features.ipynb' with functions included in 'classifiers_functions.ipynb' under 'MS classification'
- v1 update 11/02/2022, wrapped installation file/worked example to come

## GUI - the graphical user interface element for the old py_DESI_MSI workflow for imaging data analysis (REDUNDANT and UNSERVICED)
