# Usage:
# /path/ucsf-chimerax/bin/ChimeraX /absolute/path/to/chimeraxScript.py 
# With arguments in a 5 line cxScriptSettings.txt file with this program
#<folder path> <use align (a) or matchmaker (mm)> <Use uniprot (y) or in folder (n)>  <Uni ID or header 1> <Uni ID or header 2>

from chimerax.core.commands import run #type: ignore
import sys
import requests #type: ignore
import os
import json
import shutil
from PDAnalysis import Protein, AverageProtein, Deformation #type: ignore

def main():
    # If running in ChimeraX app
    #thisProgramPath = os.path.dirname(sys.argv[0])
    # If running from command line
    thisProgramPath = os.path.dirname(sys.argv[1])
    with open(os.path.join(thisProgramPath, "cxScriptSettings.txt"), "r") as file:
        lines = file.readlines()
        lines = [line.strip() for line in lines]

    if not len(lines) == 7:
        print("Error: Wrong number of arguments")
        return

    mainFolderPath = lines[0]
    rmsdFunction = lines[1] # string of "mm" or "a"
    refsSameLength = lines[3] # String of "y" or "n"
    refSource = lines[4] # String of "y" or "n"
    uniprotIDOne = lines[5]
    uniprotIDTwo = lines[6]

    # Define stuff for PDanalysis library
    min_plddt = 0
    neigh_cut = 13.0
    exitApp = False
    protein_kwargs = {"neigh_cut":neigh_cut, "min_plddt":min_plddt}

    if (not refsSameLength == "y") and (not refsSameLength == "n"):
        print("Error: Same length refs argument isn't y or n")
        return
    if (rmsdFunction == "mm"):
        rmsdFunction = "matchmake"
    elif (rmsdFunction == "a"):
        rmsdFunction = "align"
    else:
        print("Error: align/matchmaker function isn't mm or a")
        return
    

    if (not refSource == "y") and (not refSource == "n"):
        print("Error: Ref source argument isn't y or n")
        return
    

    # If uniprot (currently unused)
    if refSource == "y":

        # Get first part of headers from uniprot webserver, print error and exit program if needed
        refNameOne = getUniprotName(uniprotIDOne)
        if "failed" in refNameOne:
            print(refNameOne)
            return
        refNameTwo = getUniprotName(uniprotIDTwo)
        if "failed" in refNameTwo:
            print(refNameTwo)
            return
        # Fetch alphafold structure
        run(session, "alphafold fetch " + uniprotIDOne)
        run(session, "alphafold fetch " + uniprotIDTwo)
        # Rename structures
        run(session, "rename #1 " + refNameOne)
        run(session, "rename #2 " + refNameTwo)


    # If reference structures are in the folder, open them in ChimeraX
    elif refSource == "n":
        # This code for finding the folder name is from ChatGPT
        refNameOne = next((d for d in os.listdir(mainFolderPath) if "first_ref_" in d and os.path.isdir(os.path.join(mainFolderPath, d))), "missing")
        refNameTwo = next((d for d in os.listdir(mainFolderPath) if "second_ref_" in d and os.path.isdir(os.path.join(mainFolderPath, d))), "missing")
        if refNameOne == "missing":
            print("Error: No subfolder for first reference")
            return
        if refNameTwo == "missing":
            print("Error: No subfolder for second reference")
            return
        
        refFolderOne = os.path.join(mainFolderPath, refNameOne)
        refFileOne = getCorrectFilename(refFolderOne)
        if "could not find" in refFileOne:
            print(refFileOne)
            return
        
        refFolderTwo = os.path.join(mainFolderPath, refNameTwo)
        refFileTwo = getCorrectFilename(refFolderTwo)
        if "could not find" in refFileTwo:
            print(refFileTwo)
            return

        run(session, "open " + refFileOne)
        run(session, "open " + refFileTwo)

        run(session, "rename #1 " + refNameOne)
        run(session, "rename #2 " + refNameTwo)


    # Color the reference structures
    run(session, "color #1 pale goldenrod")
    run(session, "color #2 aquamarine")
    modelList = []
    modelList.append(refNameOne)
    modelList.append(refNameTwo)

    # Open shuffle structures
    numDropouts = 0
    i = 1
    while True:
        dropoutNameStart = "dropout" + str(i)
        dropoutName = next((d for d in os.listdir(mainFolderPath) if dropoutNameStart in d and os.path.isdir(os.path.join(mainFolderPath, d))), "missing")
        if dropoutName == "missing":
            if i == 1:
                print("Error: could not find a folder for first dropout")
                return
            else:
                break
        modelList.append(dropoutName)
        dropoutFolder = os.path.join(mainFolderPath, dropoutName)
        correctModel = getCorrectFilename(dropoutFolder)
        if "could not find" in correctModel:
            print(correctModel)
            return

        run(session, "open " + correctModel)
        modelNumber = i + 2
        run(session, "rename #" + str(modelNumber) + " " + dropoutName)

        numDropouts = numDropouts + 1
        i = i + 1
        
    # Align the reference structures and save RMSD
    refsRMSD = cxRMSD(1, 2, rmsdFunction)

    # Creates folders where aligned structures of each dropout to both refs will be saved
    individualFolderOne = os.path.join(mainFolderPath, "individual_models_aligned_" + refNameOne)
    individualFolderTwo = os.path.join(mainFolderPath, "individual_models_aligned_" + refNameTwo)
    if os.path.exists(individualFolderOne):
        shutil.rmtree(individualFolderOne)
    if os.path.exists(individualFolderTwo):
        shutil.rmtree(individualFolderTwo)
    os.makedirs(individualFolderOne)
    os.makedirs(individualFolderTwo)

    # Align and color based on first reference
    firstCXRMSDList = []
    for i in range(3, numDropouts + 3):
        firstCXRMSDList.append(cxRMSD(1, i, rmsdFunction))

    maxRMSD = max(firstCXRMSDList)
    minRMSD = min(firstCXRMSDList)

    for i in range(0, len(firstCXRMSDList)):
        normalized = (firstCXRMSDList[i] - minRMSD) / (maxRMSD - minRMSD)
        green = (1 - normalized) * 100
        red = normalized * 100
        modelNum = i + 3
        run(session, "color #" + str(modelNum) + " rgb(" + str(red) + ", " + str(green) + ", 0)")



    # Save as pdb and glb (aligned to first reference)
    pdbPath = os.path.join(mainFolderPath, "aligned_ref1_" + refNameOne + ".pdb")
    glbPath = os.path.join(mainFolderPath, "aligned_ref1_" + refNameOne + ".glb")

    if os.path.exists(pdbPath):
        os.remove(pdbPath)

    if os.path.exists(glbPath):
        os.remove(glbPath)

    run(session, "save " + pdbPath + " format pdb")
    run(session, "save " +  glbPath + " format gltf")

    for i in range(numDropouts + 2):
        run(session, "save " + os.path.join(individualFolderOne, modelList[i]) + ".pdb #" + str(i + 1))


    # Align and color based on second reference
    secondCXRMSDList = []
    for i in range(3, numDropouts + 3):
        secondCXRMSDList.append(cxRMSD(2, i, rmsdFunction))

    maxRMSD = max(secondCXRMSDList)
    minRMSD = min(secondCXRMSDList)

    for i in range(0, len(secondCXRMSDList)):
        normalized = (secondCXRMSDList[i] - minRMSD) / (maxRMSD - minRMSD)
        green = (1 - normalized) * 100
        red = normalized * 100
        modelNum = i + 3
        run(session, "color #" + str(modelNum) + " rgb(" + str(red) + ", " + str(green) + ", 0)")



    # Save as pdb and glb (aligned to first reference)
    pdbPath = os.path.join(mainFolderPath, "aligned_ref2_" + refNameOne + ".pdb")
    glbPath = os.path.join(mainFolderPath, "aligned_ref2_" + refNameTwo + ".glb")

    if os.path.exists(pdbPath):
        os.remove(pdbPath)

    if os.path.exists(glbPath):
        os.remove(glbPath)

    run(session, "save " + pdbPath + " format pdb")
    run(session, "save " +  glbPath + " format gltf")

    for i in range(numDropouts + 2):
        run(session, "save " + os.path.join(individualFolderTwo, modelList[i]) + ".pdb #" + str(i + 1))


    # ES/RMSD of refs and each shuffle to each ref
    if (refsSameLength == "Yes"):
        refsESRMSD = ESRMSD(os.path.join(individualFolderTwo, modelList[0]) + ".pdb", os.path.join(individualFolderTwo, modelList[1]) + ".pdb", protein_kwargs)
    else:
        refsESRMSD = ["NA", "NA"]

    firstESList = []
    secondESList = []
    firstOtherRMSDList = []
    secondOtherRMSDList = []
    for i in range(numDropouts):
        theESRMSD = ESRMSD(os.path.join(individualFolderOne, modelList[0]) + ".pdb", os.path.join(individualFolderOne, modelList[i+2]) + ".pdb", protein_kwargs)
        firstESList.append(theESRMSD[0])
        firstOtherRMSDList.append(theESRMSD[1])
        if (refsSameLength == "Yes"):
            theESRMSD = ESRMSD(os.path.join(individualFolderTwo, modelList[1]) + ".pdb", os.path.join(individualFolderTwo, modelList[i+2]) + ".pdb", protein_kwargs)
            secondESList.append(theESRMSD[0])
            secondOtherRMSDList.append(theESRMSD[1])
        else:
            secondESList.append("NA")
            secondOtherRMSDList.append("NA")

    # Save RMSD and ES values in a spreadsheet
    csvPath = os.path.join(mainFolderPath, "RMSD_values.csv")
    with open(csvPath, "w") as toWrite:
        if rmsdFunction == "y":
            toWrite.write("Alignment tool,align\n\n")
        elif rmsdFunction == "n":
            toWrite.write("Alignment tool,matchmaker\n\n")
        toWrite.write("Ref one to ref two\n")
        toWrite.write("ChimeraX RMSD,Effective strain,Other RMSD\n")
        toWrite.write(str(refsRMSD)+ "," + str(refsESRMSD[0]) + "," + str(refsESRMSD[1]) + "\n\n")

        toWrite.write(" ," + refNameOne + ", , , ," + refNameTwo + "\n")
        toWrite.write(" ,ChimeraX RMSD,Effective strain,Other RMSD, ,ChimeraX RMSD,Effective strain,Other RMSD")
        for i in range(1, numDropouts + 1):
            toWrite.write("\ndropout" + str(i) + "," + str(firstCXRMSDList[i]) + "," + str(firstESList[i]) + "," + str(firstOtherRMSDList[i]) + ", ," + str(secondCXRMSDList[i]) + "," + str(secondESList[i]) + "," + str(secondOtherRMSDList[i]))





    if exitApp:
        run(session, "close all")
        run(session, "exit")




# Unused function for getting the name of a protein given the uniprot ID
def getUniprotName(ID):
    url = f"https://rest.uniprot.org/uniprotkb/{ID}.fasta"
    response = requests.get(url)

    if response.status_code == 200:
        fastaData = response.text
        headerLine = fastaData.split("\n", 1)[0]
        headerLine = headerLine.split(" ")[0]
        toReturn = headerLine[1:]
        toReturn = toReturn.replace("|", "_")
        return toReturn
    else:
        return "Error: Uniprot server failed for " + ID

# Function for finding the top ranked structure given a localcolabfold output folder
def getCorrectFilename(folderPath):
    # Search each file for relaxed rank 1 or unrelaxed rank 1 and .pdb
    foundRelaxed = False
    relaxedName = ""
    unrelaxedName = ""
    for file in os.listdir(folderPath):
        filename = os.fsdecode(file)
        if ("relaxed_rank_001" in filename) and (".pdb" in filename):
            foundRelaxed = True
            relaxedName = filename
        if ("unrelaxed_rank_001" in filename) and (".pdb" in filename):
            unrelaxedName = filename
    if (relaxedName == "") and (unrelaxedName == ""):
        return "Error: Could not find a rank 1 fold in " + folderPath
    if foundRelaxed:
        return(os.path.join(folderPath, relaxedName))
    else:
        return(os.path.join(folderPath, unrelaxedName))

    # If found a relaxed rank 1 return that attached to folder path, otherwise return the unrelaxed

def cxRMSD(firstModelNum, secondModelNum, rmsdFunction):
    firstBackbone = "#" + str(firstModelNum) + "@ca"
    secondBackbone = "#" + str(secondModelNum) + "@ca"
    modelPair = firstBackbone + " to " + secondBackbone
    # Align function
    if rmsdFunction == "y":
        run(session, "align " + modelPair)
        return run(session, "rmsd " + modelPair)
    elif rmsdFunction == "n":
        jsonData = run(session, "matchmake " + modelPair)
        jsonData = jsonData[0]
        return jsonData["full RMSD"]


def ESRMSD(firstPath, secondPath, protein_kwargs):
    deform_kwargs = {'method':'all'}
    protA = Protein(firstPath, **protein_kwargs)
    protB = Protein(secondPath, **protein_kwargs)
    deform = Deformation(protA, protB, **deform_kwargs)
    deform.run()
    metrics = ["Effective Strain", "RMSD"]
    deformation_metrics = [deform.strain, deform.rmsd_per_residue]
    toReturn = []
    for entry in deformation_metrics:
        toReturn.append(sum(entry)/len(entry))
    return toReturn




main()
