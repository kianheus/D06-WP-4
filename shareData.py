# script to export data to and to update a file

# DATA SHOULD BE FORMATTED AS FOLLOWS!
data = {"T": 99999, "W": 555555, "lambda": 4444444}


# ---------- copy code below ----------
def updateFileList(data):
    print()
    print("Updating variables.dat...")
    file = open('variables.dat')  # open the file
    print("Opening the file...")

    fileList = []  # making file a list
    for line in file:
        fileList.append(line.rstrip().split())

    print("Making backup...")
    oldFile = open('oldvariables.dat', 'w')
    for line in fileList:
        oldFile.write(line[0] + ' ' + line[1] + '\n')
    oldFile.close()

    for variableName in data:  # check if variable exists
        print("Checking if", variableName, "is already a variable...")
        existsAlready = False
        sameValue = False
        for i in range(len(fileList)):
            if str(variableName) == str(fileList[i][0]):
                existsAlready = True
                foundIndex = i
                print(variableName, "already exists")
                print("Checking if the value is the same...")
                if float(data[variableName]) == float(fileList[i][1]):
                    sameValue = True
                    print(variableName, "has the same value")

        if not existsAlready:
            fileList.append([variableName, data[variableName]])
            print(variableName, "does not exists yet, appending to the list...")

        elif not sameValue:
            print()
            print("Do you want me to replace", variableName, "=", fileList[foundIndex][1], "with",
                  fileList[foundIndex][0], "=", data[fileList[foundIndex][0]], "?")
            choice = "a"
            while choice.lower() != "y" and choice.lower() != "n":
                choice = input("Enter 'Y' for yes and 'N' for no: ")

            if choice == "y":
                print("Appending your choice...")
                fileList[foundIndex][1] = data[fileList[foundIndex][0]]

    print("Writing to new file...")
    newFile = open('variables.dat', 'w')
    for i in range(len(fileList)):
        line = fileList[i][0] + " " + str(fileList[i][1]) + "\n"
        newFile.write(line)

    file.close()

    test


updateFileList(data)
