% GB comments
1.	100
2a. 100
2b. 100
2c. 100
3a 100 
3b. 100
3c. 100  	
Overall: 100


%HW3

%% Problem 1 - Smith-Waterman alignment
% Consider two sequences 'GTAATCC' and 'GTATCCG'

% Give sequences and sequence length
S1 = 'GTAATCC';
S2 = 'GTATCCG';
len = length(S1);

% Construct the scoring matrix for this with the parameters:
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 

% Define values and construct score matrix
matchval = 2;
mismatchval = -1;
gap_pen = 1;
score_matrix = eye(len) * matchval + (ones(len) - eye(len)) * mismatchval;

% Evaluate the optimal alignment with Smith-Waterman algorithm
[~, align, ~] = swalign(S1, S2, 'Alphabet', 'nt', 'ScoringMatrix', score_matrix, 'GapOpen', gap_pen);

disp(' ')
disp(align)
disp(' ')

%% Problem 2 - using the NCBI databases and sequence alignments

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 

erk1_acc = 'NM_002746';
erk2_acc = 'NM_002745';

% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 

% Get sequences from genbank for ERK1, ERK2
genbank_erk1 = getgenbank(erk1_acc);
erk1_seq = genbank_erk1.Sequence;
genbank_erk2 = getgenbank(erk2_acc);
erk2_seq = genbank_erk2.Sequence;

% Get Smith-Waterman alignment for the two sequences
[~, align, ~] = swalign(erk1_seq, erk2_seq);

% Count the aligned base pairs in the alignment and report the fraction of
% aligned base pairs.
num_align = length(find(align == '|'));
frac = num_align / length(erk1_seq);

disp(' ')
fprintf('The fraction of base pairs in ERK1 that can align to ERK2 is %1.4f. \n', frac)

% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?

% Get translated amino acid sequence for ERK1 and ERK2
erk1_aa = genbank_erk1.CDS.translation;
erk2_aa = genbank_erk2.CDS.translation;

% Get Smith-Waterman alignment for the two amino acid sequences
[~, align, ~] = swalign(erk1_aa, erk2_aa);

% Count the aligned base pairs in the alignment and report the fraction of
% aligned base pairs.
num_align = length(find(align == '|'));
frac = num_align / length(erk1_aa);

disp(' ')
fprintf('The fraction of amino acids in ERK1 that can align to ERK2 is %1.4f. \n', frac)
disp(' ')

% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 

mouse_erk1_acc = 'NM_011952';
mouse_erk2_acc = 'NM_001038663';

% Get DNA base pair and amino acid sequence info on each erk for mice
mouse_erk1 = getgenbank(mouse_erk1_acc);
mouse_erk2 = getgenbank(mouse_erk2_acc);
mouse_erk1_seq = mouse_erk1.Sequence;
mouse_erk2_seq = mouse_erk2.Sequence;
mouse_erk1_aa = mouse_erk1.CDS.translation;
mouse_erk2_aa = mouse_erk2.CDS.translation;

% Align DNA and protein sequences
[~, align_erk1_seq, ~] = swalign(erk1_seq, mouse_erk1_seq);
[~, align_erk2_seq, ~] = swalign(erk2_seq, mouse_erk2_seq);
[~, align_erk1_aa, ~] = swalign(erk1_aa, mouse_erk1_aa);
[~, align_erk2_aa, ~] = swalign(erk2_aa, mouse_erk1_aa);

% Get fraction of matching DNA base pairs and amino acids for each
frac_seq1 = length(find(align_erk1_seq == '|')) / length(mouse_erk1_seq);
frac_seq2 = length(find(align_erk2_seq == '|')) / length(mouse_erk2_seq);
frac_aa1 = length(find(align_erk1_aa == '|')) / length(mouse_erk1_aa);
frac_aa2 = length(find(align_erk2_aa == '|')) / length(mouse_erk2_aa);

% Each alignment has over 75% matches of mouse ERK1/2 base pairs or amino
% acids to human base pairs or amino acids.

%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 

% See function written below.
acc_nums = top_blast_hits('NM_002746', 3)

% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 

% See function written below.
[ acc_human, acc_nonhuman ] = species_top_hits('NM_002746')

% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 

% The accession number for the human NKX2 gene is:
nkx2 = 'NM_001079668';
% This gene plays a role in human lung development

[ nkx2_human, nkx2_nonhuman ] = species_top_hits(nkx2)

% The accession number for the rhesus monkey gene Pax6 is:
pax6 = 'NM_001266257';
% This gene is involved in the development of the eye of this species.

[ pax6_human, pax6_nonhuman ] = species_top_hits(pax6)

% The results show that the gene for each species had a homolog in both
% human and nonhuman species. This shows that many genes are conserved and
% only slightly altered by the evolutionary process.




% FUNCTIONS %

function [ acc_human, acc_nonhuman ] = species_top_hits( accession )

N = 15; % To try to ensure both a human and non-human match are found

% Use the function written for Part 1
acc_nums = top_blast_hits( accession, N );

% For first accession number, store for appropriate category
% Store data
data = getgenbank(char(acc_nums(1)));

% Check species and store
if contains(data.Source, 'Homo sapiens (human)')
    acc_human = char(acc_nums(1));
    
else
    acc_nonhuman = char(acc_nums(1));
    
end

% Iterate through all other accession numbers until both human and 
% non-human are found
for iter = 2:length(acc_nums)
    
    % Store data for iterated gene
    data = getgenbank(char(acc_nums(iter)));
    
    % Check species and store accession for category not occupied by 1st
    % accession
    if ~contains(data.Source, 'Homo sapiens (human)') && exist('acc_human', 'var')
        
        acc_nonhuman = char(acc_nums(iter));
        break
        
    elseif contains(data.Source, 'Homo sapiens (human)') && exist('acc_nonhuman', 'var')
        
        acc_human = char(acc_nums(iter));
        break
        
    end
    
end

% Return warning if no accession is found for humans
if ~exist('acc_human', 'var')
   
    acc_human = 'None';
    disp(' ')
    disp('No human gene found matching the given accession number.')   
    disp(' ')
    
end

end

function [ acc_nums ] = top_blast_hits( accession, N )

% Get data from genbank for given accession number
genbank_dat = getgenbank(accession);

% Get BLAST data for the sequence returned by genbank
[requestID, requestTime] = blastncbi(genbank_dat.Sequence, 'blastn');
blast_data = getblast(requestID, 'WaitTime', requestTime);

% Generate cell array to hold accession numbers
acc_nums = {};

% Iterate through number of hits
for i = 1:N
    
    % Add accession number to cell array
    acc_nums{i} = blast_data.Hits(i).Accession;
    
end

end


