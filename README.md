# ThesisProject
A pipeline that shuffles mutations in two homologous sequences then uses AlphaFold on the outputs to search for never-born proteins.


Requirements:
Local colab fold: https://github.com/YoshitakaMo/localcolabfold
Python 3.11.5
PyQt6 
Perl 5.31.1 (including List::Util and Time:HiRes modules)

So far tested only on MacOS


Directions:

1. Open the GUI from the command line using "python NMThesisGUI.py"



To use the mutation shuffling program:

1. Use fasta files that contain two homologous sequences 

2. Provide the path to the folder containing the fasta files, either starting from the folder containing the program files or from the root directory. Then select whether every file in the folder should be shuffled or provide names of specific fasta files

3. Do the same for the folder where the shuffle files will be created, can be a folder that does or doesn't already exist.

4. Provide the window size and the number of shuffles to create per homologous pair

5. Select whether the input files are DNA or Protein sequences, and if they are DNA sequences select whether they should be translated before or after shuffling.

6. Select whether insertions or deletions should be shuffled and select whether the input files contain only the two sequences to be used or provide the two sequence headers to use.

7. Press the "Run shuffler" button


To run mutation shuffles in AlphaFold

1. Provide the input folder, input files, and output folder in the same manner as the mutation shuffler

2. Provide the number of shuffled sequences to be run from each file and the number of recycles for each sequence that is run.



Important notes:

- A GPU is highly recommended when using AlphaFold

- Keep the three programs in the same folder and make sure the folder does not have a file called temporaryInput.fasta

- Input files can not have spaces in their names

- An error will not be shown if the incorrect fasta headers are given when "Select species" is checked

- An error will not be shown if "Entire folder" is selected and incompatible files and present in the folder for either the shuffler or AlphaFold

- Do not edit the mutation shuffle files in any way before running them in AlphaFold using the GUI


