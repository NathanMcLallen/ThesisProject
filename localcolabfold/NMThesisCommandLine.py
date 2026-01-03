import argparse
import sys
import os

parser = argparse.ArgumentParser()
parser.add_argument("input_path", type=str, help="Alignment file, FASTA file, or folder depending on step")
parser.add_argument("-o", "--output", type=str, help="FASTA file, empty folder, or none depending on step")
parser.add_argument("-x", "--cx_path", type=str)
parser.add_argument("-s", "--step", type=int, choices=[0, 1, 2, 3], help="0 for all steps, 1 to generate dropouts, 2 for AF predictions, 3 for analysis")
parser.add_argument("-rc", "--num_recycles", type=int, choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
parser.add_argument("-c", "--include_conservatives", type=str, choices=["y", "n"])
parser.add_argument("-t", "--use_templates", type=str, choices=["y", "n"])
parser.add_argument("-rx", "--use_relaxation", type=str, choices=["y", "n"])
parser.add_argument("-l", "--same_length_refs", type=str, choices=["y", "n"])
parser.add_argument("-a", "--align_algo", type=str, choices=["matchmake", "align"])
parser.add_argument("-d", "num_dropouts", type=int)
args = parser.parse_args()

if not os.path.exists(args.input_path):
    raise Exception("Input path does not exist")

if args.step:
    stepToDo = args.step
else:
    stepToDo = 0

# Errors for invalid input arguments or argument combinations
if not args.cx_path and (stepToDo == 0 or stepToDo == 3):
    raise Exception("ChimeraX path argument is required for step 3 or all steps")
elif not os.path.exists(args.cx_path):
    raise Exception("ChimeraX path is incorrect")

if not stepToDo == 3 and not args.output:
    raise Exception("Output path is needed but not provided")
#elif stepToDo == 1 and (args.output[-1] == "/" 

if args.num_recycles:
    numRecycles = args.num_recycles
else:
    numRecycles = 3

if not args.include_conservatives:
    useConservatives = "y"

if not args.use_templates:
    useTemplates = "n"

if not args.use_relaxation:
    useRelaxation = "y"

if not args.same_length_refs:
    sameLength = "n"

if args.align_algo == "align":
    alignAlgo = "a"
else:
    alignAlgo = "mm"

if alignAlgo == "a" and sameLength == "n":
    raise Exception("Refs aren't same length but align algorithm was chosen")


if stepToDo == 1:
    outputFolder = os.path.dirname(args.output)
    dropoutsPath = args.output
else:
    outputFolder = args.output

if not os.path.exists(outputFolder):
    os.makedirs(outputFolder)

if stepToDo == 0 or stepToDo == 1:
    dropoutCommand = ["python", "Dropout_generator.py", args.input_path, dropoutsPath, useConservatives, 0]
    os.system(" ".join(dropoutCommand))

if stepToDo == 0 or stepToDo == 2:
    foldingCommand = ["python", "local_colab_script.py", dropoutsPath, outputFolder, "y", numRecycles, useTemplates, useRelaxation]
    os.system(" ".join(foldingCommand))

if stepToDo == 0 or stepToDo == 3:
    with open("cxScriptSettings.txt", "w") as toWrite:
        toWrite.write(outputFolder + "\n")
        toWrite.write(alignAlgo + "\n")
        toWrite.write(sameLength + "\n")
        toWrite.write("n" + "\n")
        toWrite.write("NA" + "\n")
        toWrite.write("NA")

    thisPath = sys.path[0]
    os.system(args.cx_path + " " + os.path.join(thisPath, "chimeraxScript.py"))
    os.remove("cxScriptSettings.txt")
