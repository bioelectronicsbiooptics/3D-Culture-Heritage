function [F,R] = PrimerSelect(filename)
if strcmp(filename,'Jikji1')
    F = 'TTCGTTCGTCGTTGATTGGT';
    R = 'AAACGGAGCCATGAGTTTGT';
    
elseif strcmp(filename,'Jikji2')
    F = 'AAATCCTTTGTGCCTGCCAT';
    R = 'AAACTCAAGGCCGACCAATT';
    
elseif strcmp(filename,'Jikji3')
    F = 'AATCATGGCCTTCAAACCGT';
    R = 'AACGCTCCGAAAGTCTTGTT';
    
elseif strcmp(filename,'Jikji4')
    F = 'CTGTCCATAGCCTTGTTCGT';
    R = 'GCGGAAACGTAGTGAAGGTA';
    
elseif strcmp(filename,'Jikji5')
    F = 'TGTATTTCCTTCGGTGCTCC';
    R = 'TTTCGACAACGGTCTGGTTT';
    
elseif strcmp(filename,'Jikji6')
    F = 'AGCCTTGTGTCCATCAATCC';
    R = 'TGCGCTATGGTTTGGCTAAT';
    
elseif strcmp(filename,'Jikji7')
    F = 'GTCCAGGCAAAGATCCAGTT';
    R = 'ACCACCGTTAGGCTAAAGTG';
   
elseif strcmp(filename,'tripitaka')
    F = 'TAGCCTCCAGAATGAAACGG';
    R = 'TTCAAGCCAAACCGTGTGTA';
    
elseif strcmp(filename,'liberty')
    F = 'AAGGCAAGTTGTTACCAGCA';
    R = 'TGCGACCGTAATCAAACCAA';
    
elseif strcmp(filename,'seokga')
    F = 'ATTCGCGTCGCCTAATTTGT';
    R = 'AAACTGGAGGCGGCAAATTA';
    
elseif strcmp(filename,'pisa')
    F = 'AATGGACGTTCCGCAATCAT';
    R = 'AGAGCCGTGGCAATGTAAAT';
    
elseif strcmp(filename,'venus')
    F = 'ATCCTGCAAACGCATTTCCT';
    R = 'ATGCCTTTCCGAAGTTTCCA';
    
elseif strcmp(filename,'goldengate')
    F = 'TCCACCGTTCCTTGATTTCG';
    R = 'AATCCGTTTGCCTGCCTTTA';
end