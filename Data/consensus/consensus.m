clear all

%Import Clustal data
readtable clustal.txt;
clustal = ans;
clustal = clustal{:,1};
clustal = char(clustal);
clustal = clustal(:,35:end);

%Now we must remove *'s 
[r, c] = size(clustal);
row = [];
%for i=r:-1:1;
%    annoyance = ((clustal(i,1)=='*')||(clustal(i,1)=='*'));
%    if annoyance == 1
%        clustal(i,:)=[];
%    end
%end

%Sequences are now in 61 columns across rows, with every sequence repeating
%every four columns. We want to make it one giant string, so 4x(length)

for i=1:r
    row(i) = mod(i,4);
end

%Correct for modular making rows 4, 0
rowshift = row + 1;
for i = length(row):-1:2
    rowshift(i)=rowshift(i-1);
end
rowshift(1)=1;

%We now have clustal, which lists all sequences and another vector that
%lists the sequence identity of each row. 
alignment = [];
for i = 1:(r/4)
    alignment = [alignment,clustal(i*4-3:i*4,:)];
end

%Now to construct the consensus sequence
singlesequence = char(0);
for i = 1:length(alignment)
    singlesequence(i) = mode(alignment(:,i));
end
singlesequence = cellstr(singlesequence);
T = cell2table(singlesequence,'VariableNames',{'sequence'});
writetable(T,'tabledata.dat');

