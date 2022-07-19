function newFolder=createFolder(pathFolder,newFolderPath)
    newFolder = [pathFolder,newFolderPath]; %Macbook
    if exist(newFolder) == 0, mkdir(newFolder); end
end