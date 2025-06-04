% ---------- pretty loading-table to Command Window ---------------------
function listLoadings(coeff, rowNames, hdr)
    fprintf('\n--- %s ---\n', hdr);
    disp(array2table(coeff, RowNames=rowNames, ...
                     VariableNames={'PC1','PC2','PC3'}));
end
