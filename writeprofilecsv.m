function writeprofilecsv(outputdir,selnames,suffix,profiledatalist)

for m = 1:size(profiledatalist,2)
     fp = fopen([outputdir '/' selnames{m} '-' suffix '.csv'],'wt');
     fprintf(fp,'r1,c1,r2,c2,len\n');
     fclose(fp); % clear files if exist
end

for k=1:size(profiledatalist,1) % lines
    for m = 1:size(profiledatalist,2) % regions
        pdata = profiledatalist(k,m);
        if ~isempty(pdata.len)
            fp = fopen([outputdir '/' selnames{m} '-' suffix '.csv'],'at');
            fprintf(fp,'%.2f,%.2f,%.2f,%.2f,%.2f\n',pdata.st(1,1),pdata.st(1,2),pdata.en(1,1),pdata.en(1,2),pdata.len(1));
            fclose(fp);
        end
    end
end