In the file LULU_bug_corrected.R, there is a corrected version of the LULU program that fixes the bug reported by Mahe when running MUMU (see https://github.com/tobiasgf/lulu/issues/8).
The bug causes MOTUs to be merged only when co-occurrence is equal to 1, even if a different co-occurrence value is specified. 
The input must be a table of MOTUs in CSV format with their sequence information.
