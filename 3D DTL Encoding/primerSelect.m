function [F_primers, R_primers] = primerSelect(filename)

if strcmpi(filename,"Jikji1")
    F_primers = "TTCGTTCGTCGTTGATTGGT"; 
    R_primers = "AAACGGAGCCATGAGTTTGT";
    
elseif strcmpi(filename,"Jikji2")
    F_primers = "AAATCCTTTGTGCCTGCCAT"; 
    R_primers = "AAACTCAAGGCCGACCAATT";
    
elseif strcmpi(filename,"Jikji3")
    F_primers = "AATCATGGCCTTCAAACCGT"; 
    R_primers = "AACGCTCCGAAAGTCTTGTT";
    
elseif strcmpi(filename,"Jikji4")
    F_primers = "CTGTCCATAGCCTTGTTCGT"; 
    R_primers = "GCGGAAACGTAGTGAAGGTA";
    
elseif strcmpi(filename,"Jikji5")
    F_primers = "TGTATTTCCTTCGGTGCTCC"; 
    R_primers = "TTTCGACAACGGTCTGGTTT"; 
    
elseif strcmpi(filename,"Jikji6")
    F_primers = "AGCCTTGTGTCCATCAATCC"; 
    R_primers = "TGCGCTATGGTTTGGCTAAT";
    
elseif strcmpi(filename,"Jikji7")
    F_primers = "GTCCAGGCAAAGATCCAGTT"; 
    R_primers = "ACCACCGTTAGGCTAAAGTG";
    
elseif strcmpi(filename,"tripitaka")
    F_primers = "TAGCCTCCAGAATGAAACGG"; 
    R_primers = "TTCAAGCCAAACCGTGTGTA"; 
    
elseif strcmpi(filename,"liberty")
    F_primers = "AAGGCAAGTTGTTACCAGCA"; 
    R_primers = "TGCGACCGTAATCAAACCAA";
    
elseif strcmpi(filename,"seokga")
    F_primers = "ATTCGCGTCGCCTAATTTGT";
    R_primers = "AAACTGGAGGCGGCAAATTA";
    
elseif strcmpi(filename,"pisa")
    F_primers = "AATGGACGTTCCGCAATCAT"; 
    R_primers = "AGAGCCGTGGCAATGTAAAT";
end

F_primers = char(F_primers);
R_primers = char(R_primers);