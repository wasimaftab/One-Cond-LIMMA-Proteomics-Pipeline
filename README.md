# This is a pipeline to analyze proteomic data (SILAC-style ratios) that is already log2 transformed with one condition.

## The accepted data structure is as follows:
- The data file should be in xlsx format.
- First column is the names of genes/proteins followed by data replicates.
- The column names for replicates must have some part in common starting from the beginning of the string, for example,  "Log2R_1"    "Log2R_2"    "Log2R_3"    "Log2R_4". Here "Log2R_" part is common in all the replicates.

## Note:
- Make sure to change the line `myFilePath <- "Parabolic_volcano/1110_data_before.xlsx"` for your data and path.
- Comment the temporary line `colnames(data)[2:ncol(data)] <- c("Log2R_1", "Log2R_2", "Log2R_3", "Log2R_4")`.

    
