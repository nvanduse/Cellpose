Installation:

Install Python for Windows (v3.12) and pip
Open command prompt
Install Cellpose with the following command: pip install cellpose[gui]
Start the Cellpose GUI by entering the following command: cellpose

Make a .bat file on your desktop called: "run_cellpose.bat". Open it with a text editor and save the shell script into it.
Make a folder on your desktop named: input_cellpose
Make another folder on your desktop called: "scripts_and_notes_cellpose" and put the python script (cellpose_quantify.py) into it


Put your images into the input folder and double click the .bat file. A command line window will pop up and 
give you some options for running the segmentation.

If the segmentation isn't working well it's probably because the auto settings for cell size aren't right.
To figure out the right cell size setting, open up the cellpose GUI and drag & drop an image in. The red
dot in the bottom left of the image is the expected cell size. If it's too big/small then adjust the size using
the menu in the left panel, until the red dot is about the same size as a cell or nuclei; then run the script,
select "M" for manual specification, and enter the correct size.

I'm not sure what the difference is between the cyto and nuclei models, but the cyto model seems to work better
for nuclei of isolated cardiomyocytes than the nuclei model (I suspect due to their elongated non-spherical shape).


I think the GUI has options for fine-tuning a model by manually adding or removing annotations, but the standard models
seem pretty good and I haven't gotten around to figuring that out and adding compatibility with the batch script yet.
