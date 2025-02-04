clear all; close all;

fnames = ["T7", "T7"]

for i = 1:length(fnames)
    fname = strcat("C:\Users\Zi Wang\Desktop\",fnames(i), ".txt")
    data = importXPfile(fname,2, inf);
    save([fname ".mat"],"data");
end