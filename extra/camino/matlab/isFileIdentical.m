function p = isFileIdentical(file1, file2)

fid1 = fopen(file1); a = fread(fid1, 'float', 'b'); fclose(fid1);
fid2 = fopen(file2); b = fread(fid2, 'float', 'b'); fclose(fid2);

p = isequalwithequalnans(a,b);

if nargout==0 && p==1
   display('Files are identical')
elseif nargout==0 && p==0
   display('Files are different')
end
    