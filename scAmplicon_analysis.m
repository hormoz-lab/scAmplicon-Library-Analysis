clear all

GeneName = ''; %Add gene name, for example: 0114_Jak2


fileName = ''; %Add file name, for exmpale: 0114_Jak2_Miseq
%Parse barcodes from read1

%Read all 4 lanes
bc_reads1 = fastqread([fileName, '']);%Finish file name, for example: _L001_R1_001.fastq
bc_reads2 = fastqread([fileName, '']);
%bc_reads3 = fastqread([fileName, '']);
%bc_reads4 = fastqread([fileName, '']);
%Aggregate all the reads
bc_reads = [bc_reads1, bc_reads2];%, bc_reads3, bc_reads4]; %, bc_reads3, bc_reads4];

bcsOrig = {bc_reads.Sequence};
bcsQCOrig = cellfun(@(x) double(x)-33, {bc_reads.Quality}, 'UniformOutput', false); %in integer format

%Read in the sequences R2
reads1 = fastqread([fileName, '']);%Finish file name, for example: _L001_R2_001.fastq
reads2 = fastqread([fileName, '']);
%reads3 = fastqread([fileName, '']);
%reads4 = fastqread([fileName, '']);
%Aggregate all the reads
reads = [reads1, reads2];%, reads3, reads4];

seqsOrig = {reads.Sequence};
seqsQCOrig = cellfun(@(x) double(x)-33, {reads.Quality}, 'UniformOutput', false); %in integer format

disp('Number of reads:')
length(reads)

%NextSeq version only reads 98 bp into the Gene2 sequence

%%

%Correct Gene2 sequence, 91 bp downstream of and including the
%locus_specific primer 3, for example, the correct gene 2 for jak2 is 
%GCAGCAAGTATGATGAGCAAGCTTTCTCACAAGCATTTGGTTTTAAATTATGGAGTATGTGTCTGTGGAGACGAGAATATTCTGGTTCAGG
corrGene2 = ''
            
posMut = ; %Location of the mutation, for example: 61 for jak2
baseMut = '' %Mutated sequence, for example: T for jak2

%Look +- 5 bp of the mutation for the fragment libraries

wtCompSeq = [corrGene2(posMut-5:posMut-1) corrGene2(posMut) corrGene2(posMut+1:posMut+5)]
mutCompSeq = [corrGene2(posMut-5:posMut-1) baseMut corrGene2(posMut+1:posMut+5)]

%%
%Go through all sequences and compute their distance to corrJake

qc_reads = [];
qc_read1 = [];
qc_bcs = [];
qc_umis = [];
qc_bc_umis = [];
qc_index = [];
selected_index = [];
Gene2BCs = [];
dnaBase = ['A', 'C', 'G', 'T'];


bcQCTot = zeros(1,100);
seqQCTot = zeros(1,100);
ctrQC = 0;
ctrIndex = 0;

tmpTotal = zeros(1,26);

for i = 1:length(reads)
    
    if mod(i,round(length(reads)/100)) == 0
        perc = round(i/length(reads)*100)
    end
    
    thisRead = reads(i).Sequence;
    thisIndex = bc_reads(i).Header(end-7:end);
    ctrIndex = ctrIndex + 1;
    
    if mean(bcsQCOrig{i}(1:26)) >= 30 && mean(seqsQCOrig{i}) >= 30 && length(thisRead) >= 80 %&& indexMatch
  
        %Number of reads that pass the QC threshold
        ctrQC = ctrQC + 1;
        
    
        qc_reads{ctrQC} = thisRead;
        qc_bcs{ctrQC} = bc_reads(i).Sequence(1:16);
        qc_umis{ctrQC} = bc_reads(i).Sequence(17:28);
        qc_bc_umis{ctrQC} = bc_reads(i).Sequence(1:28);
        qc_index{ctrQC} = thisIndex;
        
    end
    
        
end


disp('Fraction that pass QC threshold:')
ctrQC / ctrIndex

%%
%Plot the distrubtion of indices

[uniqueIndices,~,idx] = unique(qc_index);
numOccurrences = histcounts(idx,numel(uniqueIndices));

fig = figure
bar(numOccurrences)
set(gca,'fontsize', 18)
xlabel('Unique Library Index','fontsize', 18);
ylabel('Number of reads','fontsize', 18);
savefig('fig1.fig')


[rankOfOccurrences,rankIndex] = sort(numOccurrences,'descend');
%Display the top 4 indices
for k = 1:4
    uniqueIndices{rankIndex(k)}
end
%%
%Correct barcodes

%Two options
% 1: Take the detected Miseq barcodes as sequenced
% 2: Collapse the detected Miseq barcodes to the list of known NovaSeq
% barcodes
option = 2;

clear uniqueBCs idx
[uniqueBCs,~,idx] = unique(qc_bcs);

selected_reads = [];
selected_read1 = [];
selected_bcs =[];
original_bcs = [];
selected_umis = [];
selected_bc_umis = [];

%List of barcodes
barcode_inputfile = '0114barcodes.tsv';

%List of known barcodes
if option == 2 %load the NovaSeq barcodes

    fid = fopen(barcode_inputfile, 'r');

    clear bc10X

    nLines = 0;
    thisLine = fgetl(fid);
    while (thisLine ~= -1)

      nLines = nLines + 1;
      bc10X{nLines} = thisLine(1:16);
      %Read next line
      thisLine = fgetl(fid);

    end

    nLines
end


ctr = 0;

%Go through the list
for i = 1:numel(uniqueBCs)
    
    if mod(i,round(numel(uniqueBCs)/100)) == 0
        perc = round(i/numel(uniqueBCs)*100)
    end

    thisBC = uniqueBCs{i};
    
    if option == 2 %NovaSeq 10x barcodes
    
        %Match against NovaSeq 10x barcodes
        %Correct it
        %Find the BC in the list of known barcodes
        found = 0;
        myBCdist = zeros(nLines,1);

        for q = 1:nLines
            myBCdist(q) = sum(thisBC ~= bc10X{q});
        end

        [qval, qpos] = min(myBCdist);

        %**Key parameter. At what distance should we collapse the barcodes.
        if  qval <= 2 && sum(myBCdist == qval) == 1
            found = 1;
            corrBC = bc10X{qpos};
        end
   
    elseif option == 1 %Just take the detected barcode
        corrBC = thisBC; 
        found = 1; 
    end
    
    %If barcode match found
    if found
       tmpPos = find(idx == i);
        
       for k = 1:length(tmpPos)
            ctr = ctr + 1;  
            selected_reads{ctr} = qc_reads{tmpPos(k)};
            selected_bcs{ctr} = corrBC;
            original_bcs{ctr} = qc_bcs{tmpPos(k)}; %Store original barcode
            selected_umis{ctr} = qc_umis{tmpPos(k)};
            selected_bc_umis{ctr} = [corrBC,qc_umis{tmpPos(k)}];

       end
    end
    
    
end

ctr/ctrIndex


%%

%Number of unique molecules
clear uniqueMols idx
[uniqueMols,~,idx] = unique(selected_bc_umis);

numOccurrences = histcounts(idx,numel(uniqueMols));

[rankOfOccurrences,rankIndex] = sort(numOccurrences,'descend');
fig = figure
loglog(rankOfOccurrences,'bo');
set(gca,'fontsize', 18)
xlabel('Rank of molecule','fontsize', 18);
ylabel('Number of reads','fontsize', 18);
print(fig,[GeneName,'_ReadCounts_Mols'],'-dpng')
savefig('fig2.fig')

%%


%Number of reads
disp('Number of reads:')
ctr
disp('Number of unique molecules:')
length(uniqueMols)
disp('Number of molecules with >= 1000 reads')
sum(numOccurrences>=1000)
disp('Number of molecules with >= 100 reads')
sum(numOccurrences>=100)
disp('Number of single-cell barcodes:')
uniqueBCs = unique(selected_bcs);
length(uniqueBCs)

%%
%Plot hte distribution of occurrence of single cell barcodes
clear uniqueBCs idx
[uniqueBCs,~,idx] = unique(selected_bcs);
numOccurrences = histcounts(idx,numel(uniqueBCs));

[rankOfOccurrences,rankIndex] = sort(numOccurrences,'descend');
fig = figure
loglog(rankOfOccurrences,'bo');
set(gca,'fontsize', 18)
xlabel('Rank of sc barcode','fontsize', 18);
ylabel('Number of reads','fontsize', 18);
print(fig,[GeneName,'_ReadCounts_SC_Barcodes'],'-dpng')
savefig('fig3.fig')

%%
%Get the most common R2
[uniqueReads,~,idx] = unique(selected_reads);
numOccurrences = histcounts(idx,numel(uniqueReads));

[rankOfOccurrences,rankIndex] = sort(numOccurrences,'descend');
fig = figure
loglog(rankOfOccurrences,'bo');
set(gca,'fontsize', 18)
xlabel('Unique R2','fontsize', 18);
ylabel('Number of reads','fontsize', 18);
print(fig,'ReadCounts_R2','-dpng')
savefig('fig4.fig')

%Write this informatio to a file
fileID = fopen([GeneName,'_TopReads.txt'],'w');

[rankOfOccurrences,rankIndex] = sort(numOccurrences,'descend');
%Display the top reads
for k = 1:200
    [Score2, Alignment2] = nwalign(corrGene2,uniqueReads{rankIndex(k)});
    fprintf(fileID,'Rank: %d\n',k);
    fprintf(fileID, '%s\n',Alignment2(1,:));
    fprintf(fileID, '%s\n',Alignment2(2,:));
    fprintf(fileID, '%s\n',Alignment2(3,:));
    fprintf(fileID,'Alignmennt score: %d\n',Score2);
    fprintf(fileID,'Number of reads: %d\n',rankOfOccurrences(k));
    fprintf(fileID,'\n');
end

fclose(fileID)

%%

%Go through the barcodes
clc
clear uniqueBCs idx

[uniqueBCs,~,idx] = unique(selected_bcs);

%For each unique barcode count the number of unique UMIs

numOccurrences = histcounts(idx,numel(uniqueBCs));

size(uniqueBCs)


%Initialize
scGene2Data = [];
scBCSeq = [];
scBCSeqOrig = [];
scUMI = [];

ctrCell = 0;


for i = 1:numel(uniqueBCs) %go through all distinct barcodes
    
    if numOccurrences(i) >= 1 %** Pre-selection / if at least a minimum number of reads detected for that cell
        
        ctrUMI = 0;

        if mod(i,round(numel(uniqueBCs)/100)) == 0
            perc = round(i/numel(uniqueBCs)*100)
        end

        clear tmPos
        tmpPos = find(idx == i);

        clear thisUMIs
        clear thisUMIPos

       for k = 1:length(tmpPos)
           thisUMIs{k} = selected_umis{tmpPos(k)};
           thisUMIPos(k) = tmpPos(k);
        end
        %Number of unique UMIs for this barcode
        [uniqueUMIsBC,~,idxUMIs] = unique(thisUMIs);

        foundGene2Cell = 0;

        for j = 1:length(uniqueUMIsBC) %Go through all unique UMIs

            posUMI = find(idxUMIs==j);
            ctrReads = 0;

            foundGene2UMI = 0;

            for l = 1:length(posUMI) %Go through all reads of the same UMI
                %Get the Gene2 sequence
                thisR2Seq = selected_reads{thisUMIPos(posUMI(l))};
                
                %Map read2
                [Score2, Alignment2] = nwalign(thisR2Seq,corrGene2);
                
            
                if  Score2 > 150 %Record the read if either read mapped to Gene2

                    %If the frist read detected, record the cell
                    if ~foundGene2Cell
                        ctrCell = ctrCell + 1;
                        scBCSeq{ctrCell} = uniqueBCs{i};
                        foundGene2Cell = 1;
                    end

                    if ~foundGene2UMI
                        ctrUMI = ctrUMI + 1;
                        scUMI{ctrCell}{ctrUMI}.sequence = uniqueUMIsBC{j};
                        scUMI{ctrCell}{ctrUMI}.totnumReads = length(posUMI);
                        foundGene2UMI = 1;
                    end

                    %Record
                    ctrReads = ctrReads + 1;
                    scGene2Data{ctrCell}{ctrUMI}{ctrReads}.R2Seq = thisR2Seq;
                    scGene2Data{ctrCell}{ctrUMI}{ctrReads}.isR2Gene2 = Score2;
       
                end
        
            end %End going through all reads of the same UMI

    end %End going through all UMIs
    
    

    end %pre-selection crtierion
    
end %go through all distinct barcodes

%%

%Correct UMIs

%Correct the UMIs by collapsing UMIs at certain hamming distance


cellUMISeqCorr = [];
cellUMIReadsCorr = [];
cellUMIRead2Corr = [];
cellUMIReadPos = [];
scUMICounts2 = [];
scGene2Counts2 = [];


scGene2DataCorr = [];
scUMICorr = [];

thresholdUMI= 15; %threshold of similariy to combine UMIs

clc

for icell = 1:length(scBCSeq)
    
    if mod(icell,round(length(scBCSeq)/100)) == 0
        perc = round(icell/length(scBCSeq)*100)
    end

    clear umiSeq
    umiNum = [];
   
    for k = 1:length(scUMI{icell})
        umiSeq(k,:) = scUMI{icell}{k}.sequence;
        umiNum(k) = scUMI{icell}{k}.totnumReads;
    end

 
    myCell = [];
    myUMI = [];
    myOrigPos = [1:length(umiNum)];
    
    
 
    ctrUMI = 0;
    while ~isempty(umiNum)

        [maxNum, maxPos] = max(umiNum);
        maxSeq = umiSeq(maxPos,:);

        
        %Keep the UMI with the maximum number of reads
        ctrUMI = ctrUMI + 1;
        myCell{ctrUMI} = scGene2Data{icell}{myOrigPos(maxPos)};
        myUMI{ctrUMI} = scUMI{icell}{myOrigPos(maxPos)};
      
        %find the most common UMI and compute all distances to it
        udist = [];
        for k = 1:length(umiNum)
             %Many of the UMIs are shifted versions. To correct for this allign
             %the sequences insteed
             udist(k) = sum(umiSeq(k,:)~=maxSeq);
             
            if udist(k) <= 2  && k ~= maxPos
                %Merge with the max UMI
                 myCell{ctrUMI} = [myCell{ctrUMI} scGene2Data{icell}{k}];
                 myUMI{ctrUMI}.totnumReads = myUMI{ctrUMI}.totnumReads + scUMI{icell}{k}.totnumReads;
             end
                 
        end
        
   
        remUMIs = udist <= 2; %use the alignment score
        %Remove the redundant UMIs
        umiSeq(remUMIs,:) = [];
        umiNum(remUMIs) = [];
        myOrigPos(remUMIs) = [];
    
    end
    
    scGene2DataCorr{icell} = myCell;
    scUMICorr{icell} = myUMI;
    

end


%%

%Analyze the results

%Go through all cells

clc

numCells = length(scBCSeq)

scUMI_count = [];
UMI_read_count = [];

sc_Gene2Detected = []; %vector of 3, total Gene2, WT Gene2, Cancer Gene2
sc_Gene2NumReads = [];
sc_Gene2R1Detected = [];
umi_Gene2Detected = [];
umi_Gene2R1Detected = [];
umi_Gene2R1pos = [];
umi_fracGene2 = [];
umi_numReads = [];

clear scInfo

ctrUMI = 0; 
for icell = 1:numCells
    
    if mod(icell,round(numCells/100)) == 0
        perc = round(icell/numCells*100)
    end
    
    scUMI_count(icell) = length(scUMICorr{icell});
    sc_Gene2Detected(icell,:) = [0 0 0];
    sc_Gene2NumReads(icell,:) = [0 0 0];
    sc_Gene2R1Detected(icell) = 0;
    
    %Go through each UMI
    ctrx = 0;
    for j = 1:length(scUMICorr{icell})
        
        ctrUMI = ctrUMI + 1;
        UMI_read_count(ctrUMI) = scUMICorr{icell}{j}.totnumReads;
        
        Gene2R2reads = 0;
        Gene2R1reads = 0;
        Gene2R2pos = [];
        Gene2R2call = []; %Store if WT or cancer
        Gene2R1pos = [];
        Gene2R2AlgnScore = [];
        
        %Go through each read
        for k = 1:length(scGene2DataCorr{icell}{j})
            
            if scGene2DataCorr{icell}{j}{k}.isR2Gene2 >= 150
                
                Gene2R2reads = Gene2R2reads + 1;
                Gene2R2AlgnScore(Gene2R2reads) = scGene2DataCorr{icell}{j}{k}.isR2Gene2;
                
                %Check for the cancer mutation 0. undetermined; 1. WT; 2.
                %cancer
                if ~isempty(strfind(scGene2DataCorr{icell}{j}{k}.R2Seq, wtCompSeq)) %WT
                    Gene2R2call(Gene2R2reads) = 1;
                elseif ~isempty(strfind(scGene2DataCorr{icell}{j}{k}.R2Seq, mutCompSeq)) %Cancer
                    Gene2R2call(Gene2R2reads) = 2;
                else
                    Gene2R2call(Gene2R2reads) = 0;
                end
                
                
            end
            
            
        end %reads
        
        %Store the information for each UMI
        umi_Gene2Detected(ctrUMI,:) = [0 0 0];
        umi_fracGene2(ctrUMI) = 0;
        
        if Gene2R2reads >= 1000 %Please pick this threshold based on the knee plot
            ctrx = ctrx + 1;
            scInfo{icell}{ctrx}.isGene2 = 1;
            scInfo{icell}{ctrx}.Gene2R2reads = Gene2R2reads;
            scInfo{icell}{ctrx}.numWT = sum(Gene2R2call == 1);
            scInfo{icell}{ctrx}.numCancer = sum(Gene2R2call == 2);
            scInfo{icell}{ctrx}.algnScores = Gene2R2AlgnScore;
            scInfo{icell}{ctrx}.UMISeq =  scUMICorr{icell}{j}.sequence;
            scInfo{icell}{ctrx}.cellBC =  scBCSeq{icell};
            
            %Record the number of reads
            umi_Gene2Detected(ctrUMI,1) = Gene2R2reads; %all Gene2 reads
            umi_Gene2Detected(ctrUMI,2) = sum(Gene2R2call == 1); %all WT reads
            umi_Gene2Detected(ctrUMI,3) = sum(Gene2R2call == 2); %all cancer reads
            scUMICorr{icell}{j}.fracReadsGene2 = Gene2R2reads / length(scGene2DataCorr{icell}{j});
            umi_fracGene2(ctrUMI) = Gene2R2reads / length(scGene2DataCorr{icell}{j});
            umi_numReads(ctrUMI) = length(scGene2DataCorr{icell}{j});
            sc_Gene2Detected(icell,1) = sc_Gene2Detected(icell,1) + 1;
            sc_Gene2NumReads(icell,1) = sc_Gene2NumReads(icell,1) + Gene2R2reads;
            %Mutation call
            scUMICorr{icell}{j}.call = mode(Gene2R2call);
            %if mode(Gene2R2call) == 1 %WT
            if sum(Gene2R2call==1) >= 0.5*Gene2R2reads %WT
                sc_Gene2Detected(icell,2) = sc_Gene2Detected(icell,2) + 1;
                sc_Gene2NumReads(icell,2) = sc_Gene2NumReads(icell,2) + sum(Gene2R2call==1);
            %elseif mode(Gene2R2call) == 2 %Cancer
            elseif sum(Gene2R2call==2) >= 0.5*Gene2R2reads %Cancer
                sc_Gene2Detected(icell,3) = sc_Gene2Detected(icell,3) + 1;
                sc_Gene2NumReads(icell,3) = sc_Gene2NumReads(icell,3) + sum(Gene2R2call==2);
            end
            
            %*Hack
            %if sum(sc_Gene2NumReads(icell,2:3)) > 20000
                %scInfo{icell}{ctrx}.cellBC
            %    scInfo{icell}{ctrx}.UMISeq
            %end
                
        end
        
    
    end %UMIs
    
    
end %Cells


%%

%Diplay the results

disp(['Number of unique cells detected: ', num2str(numCells)])
disp(['Number of unique cells with correct Gene2 R2: ', num2str(sum(sc_Gene2Detected(:,1)>0))])
disp(['Number of unique cells with correct WT Gene2 R2: ', num2str(sum(sc_Gene2Detected(:,2)>0))])
disp(['Number of unique cells with correct Cancer Gene2 R2: ', num2str(sum(sc_Gene2Detected(:,3)>0))])
disp(['Number of unique heterozygous cells: ', num2str(sum(sc_Gene2Detected(:,3)>0 & sc_Gene2Detected(:,2)>0))])
disp(['Fraction of cells with corr Gene2 R2: ', num2str(sum(sc_Gene2Detected(:,1)>0)/numCells)])

disp(['Total number of Gene2 R2 molecules: ', num2str(sum(umi_Gene2Detected(:,1)>0))])
disp(['Total number of WT Gene2 R2 molecules: ', num2str(sum(umi_Gene2Detected(:,2)>0.9*umi_Gene2Detected(:,3)))])
disp(['Total number of Cancer Gene2 R2 molecules: ', num2str(sum(umi_Gene2Detected(:,3)>0.9*umi_Gene2Detected(:,2)))])
disp(['Fraction of molecules that are Gene2 R2 molecules: ', num2str(sum(umi_Gene2Detected(:,1)>0) / ctrUMI)])
disp(['Fraction of reads that are Gene2 R2 molecules: ', num2str(sum(umi_Gene2Detected(:,1)) / length(reads))])

%%
%Map to the NextSeq barcodes

%List of known barcodes
fid = fopen(barcode_inputfile, 'r');

clear bcNS

nLines = 0;
thisLine = fgetl(fid);
while (thisLine ~= -1)
    
  nLines = nLines + 1;
  bcNS{nLines} = thisLine(1:16);
  %Read next line
  thisLine = fgetl(fid);
  
end

nLines

numGene2WT = zeros(nLines,1);
numGene2Cancer = zeros(nLines,1);



mapped = zeros(numCells,1);

%Store the data
clear matchedBCs matchedUMIs
ctrmatch = 0;

%Go through the cells and assign ground truth
for i = 1:numCells
    
    if mod(i,round(numCells/100)) == 0
        perc = round(i/numCells*100)
        %waitbar(perc/100,h,sprintf('%d%% along...',perc))
    end
    
    thisBC = scBCSeq{i};
    %Find this barcode in the clusterMaps
    bcDist = ones(nLines,1).*16;
    for k = 1:nLines
        bcDist(k) = sum(bcNS{k} ~= thisBC);
    end
    
    [minVal, minPos] = min(bcDist);
    
    if minVal == 0 && sum(sc_Gene2Detected(i,2:3))>0 %If unique match found at a distance of 2 bps
        numGene2WT(minPos) = sc_Gene2Detected(i,2);
        numGene2Cancer(minPos) = sc_Gene2Detected(i,3);
        mapped(i) = minPos;
        %Store the barcode
        ctrmatch = ctrmatch + 1;
        matchedBCs{ctrmatch}.bcSeq = thisBC;
        matchedBCs{ctrmatch}.numGene2 = sc_Gene2Detected(i,:);
        matchedUMIs{ctrmatch} = scUMICorr{i};
        %Get the UMIs also
        
    end
    
end

%Write to file for Python visualization
fid = fopen(['Python_AnnData_',GeneName,'.csv'], 'w');
fprintf(fid, [GeneName,'WT,',GeneName,'Cancer\n']);
for i = 1:nLines
    fprintf(fid, '%s, %s\n', num2str(numGene2WT(i)), num2str(numGene2Cancer(i)));
end
fclose(fid);
    


disp(['Number of nextseq cells: ', num2str(nLines)]);
disp(['Number of nextseq cells with at least 1 wt Gene2 transcript: ', num2str(sum(numGene2WT>0))]);
disp(['Number of nextseq cells with at least 1 cancer Gene2 transcript: ', num2str(sum(numGene2Cancer>0))]);



