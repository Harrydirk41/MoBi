function sb_csv_to_xlsx(csvPath, xlsxPath)
%SB_CSV_TO_XLSX Convert a CSV to XLSX with column names preserved verbatim, so a Python-written
%   Vpop table (row 1 = exact parameter names) becomes the .xlsx sb_run_vpop expects - using only
%   base MATLAB (readtable/writetable), no Python spreadsheet package needed.
    t = readtable(char(csvPath), 'VariableNamingRule', 'preserve');
    writetable(t, char(xlsxPath));
end
