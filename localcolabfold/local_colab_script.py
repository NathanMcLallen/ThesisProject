import os
import sys
from Bio import SeqIO #type: ignore

# Read command line arguments
if not len(sys.argv) == 7:
    raise Exception("Wrong number of args, should have 6 args but got " + str(len(sys.argv) - 1))

inputFile = sys.argv[1]
outputFolder = sys.argv[2]
foldRefs = sys.argv[3]
numRecycles = sys.argv[4]
templates  = sys.argv[5]
relaxation  = sys.argv[6]

# Various errors for bad arguments
if not os.path.exists(inputFile):
    raise FileNotFoundError(inputFile + " doesn't exist")
if (not foldRefs == "y") and (not foldRefs == "n"):
    raise ValueError("Argument for whether or not to fold refs isn't y or n")
try:
    numRecycles = int(numRecycles)
except:
    raise TypeError("Argument for number of recycles isn't an integer")
if (not templates == "y") and (not templates == "n"):
    raise ValueError("Argument for whether or not to use templates isn't y or n")
if (not relaxation == "y") and (not relaxation == "n"):
    raise ValueError("Argument for whether or not to use relaxation isn't y or n")

# Define the command that will be reused for each sequence
generalCommand = "colabfold_batch"
if relaxation == "y":
    generalCommand = generalCommand + " --amber"
if templates == "y":
    generalCommand = generalCommand + " --templates"
generalCommand = generalCommand + " --num-recycle " + str(numRecycles) + " "

# For each sequence in the input file, create a temporary fasta with just one sequence and call local colab fold
i = 0
for record in SeqIO.parse(inputFile, "fasta"):
    i = i + 1
    if (i == 1 or i == 2) and (not foldRefs):
        continue

    tempFilePath = os.path.join(outputFolder, "temporary.fasta")
    with open(tempFilePath, "w") as tempFile:
        tempFile.write(">" + record.description + "\n")
        tempFile.write(str(record.seq))
    
    specificOutput = os.path.join(outputFolder, record.description)

    os.system(generalCommand + " " + tempFilePath + " " + specificOutput)

    os.remove(tempFilePath)