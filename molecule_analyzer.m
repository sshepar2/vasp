%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                   %%%%%%%%%%%%%%
%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%
%%%%%%%%%%%   %                                            %%%%%   %%%%%%%%
%%%%%%%%%%%   %                                              %%%%%   %%%%%%
%%%%%%%%%%%%   %                                             %%%%%%%   %%%%
%%%%%%%%%%%%   %                                             %%%%%%%%   %%%
%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%    %%%%%%%%%%%%            %%%%%%%%%   %%%
%%%%%%%%%%%%%%%    %%%%%%%         %%%%%%%%%%              %%%%%%%%%   %%%%
%%%%%%%%%%%%%%%%%%   %            %%%%%%%                 %%%%%%%%%   %%%%%
%%%%%%%%%%%%%%%%%   %            %%%%                    %%%%%%%%%   %%%%%%
%%%%%%%%%%%%%%%%   %            %                     %%%%%%%%%%%   %%%%%%%
%%%%%%%%%%%%%%%   %                                %%%%%%%%%%%   %%%%%%%%%%
%%%%%%%%%%%%%%   %                              %%%%%%%%%%%   %%%%%%%%%%%%%
%%%%%%%%%%%%%   %             %                    %%%%   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   %             %%%%%                  %%%%%%   %%%%%%%%%%%%%%
%%%%%%%%%%%   %             %%%%%%%                 %%%%%%%%   %%%%%%%%%%%%
%%%%%%%%%%   %             %%%%%%%%%               %%%%%%%%%%   %%%%%%%%%%%
%%%%%%%%%   %             %%%%%%%%%%              %%%%%%%%%%%   %%%%%%%%%%%
%%%%%%%%   %             %%%%%%%%%               %%%%%%%%%%%   %%%%%%%%%%%%
%%%%%%%   %                                     %%%%%%%%%%%   %%%%%%%%%%%%%
%%%%%%   %                                     %%%%%%%%%%%   %%%%%%%%%%%%%%
%%%%%   %                                     %%%%%%%%%%%   %%%%%%%%%%%%%%%
%%%%   %                                   %%%%%%%%%%%%%   %%%%%%%%%%%%%%%%
%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%
%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%
%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%
%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                            %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuart Shepard %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Binghamton University %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2015 - 2020 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOLECULE ANALYZER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% This code reads in a file formatted like the VASP output file XDATCAR
% resulting from a molecular dynamics calculation. The numbering within the 
% headers for each time step are not important, as long as the sets of
% atomic coordinates are in chronological order and each set of coordinates
% has a header line starting with 'Direct..'. This means XDATCAR files from
% consecutive molecular dynamics runs can be concatinated as long as the
% 'top header' of all XDATCARS, besides the first one, are deleted prior to
% concatenation. 'top header' here refers to the first 7 lines of the
% standard XDATCAR/POSCAR/CONTCAR VASP out put files: title line, scale
% factor, basis vector 1, basis vector 2, basis vector 3, atom type, atom
% type numbers. Note XDATCAR output does not include a 'Selective Dynamics' 
% line even if it was used in the initial POSCAR file and calculation.
% In this script such an extra line is not anticipated in the 'top header'
% since it would not be there.
%
% This code is broken into two parts,
%
% 1. Molecular Cluster Algorithm    2. Molecular Analysis
%
% The first part will determine which atoms make up the different molecules
% by looking at the first time step of data only and save the necessary
% data. Once the algorithm has been performed once it can be switched off
% and the saved data is loaded instead. If the atom cluster data is known
% prior it can simply be loaded by providing the data in the same format
% and the alogrithm switched off.
% 
% The second part calculates various properties of the molecules (bond
% lengths, breaking, time of breaking, location of molecule when breaking,
% etc.) over the entire trajectory specified and the time samplings
% specified. Data is saved to different variables for plotting and further
% analysis.
%
% As of 2020-03-20 the code is very specific to a system of ethylene 
% carbonate molecules in parts 1 and 2. The code is fairly close to being
% general. For instance one can write their own cluter algorithm, or just
% provide data to load. Part 2 could easily made to accept general cluster
% data but is currently not set up to do so.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% HOW TO USE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% There are TWO user input sections with 'User Input Required' written and
% surrounded by '%%...' so they are easy to spot. The first one in located
% in the 'Initial User Input' section and the second located in the
% 'Molecular Analysis' section. The input variables have explanations of 
% what they mean next to them.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% CHANGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Make it so the lines_keep or atype variables can be loaded and not
% required to be calculated prior. So add a section in the loading portion
% which calculates lines_keep variable from a loaded atype variable. 
% Essentially repeat whats in the Initial user input section.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Initial User Input %%

tic;
fprintf('\n...reading user input\n\n\t...time elapsed: %f seconds\n\n\n',toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% User Input Required %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change variable below to path of XDATCAR type file
xdatfile='XDATCAR_2_mol_relax';

% Time info
ts=1; % time step used in MD calculation (in femtoseconds)
sample = 1; % sampling rate in units of ts. will keep ts 1 and ts 1+sample..
time_start = 1; % analyze starting at this time step in units of time step (ts)
time_end = 5820; % analyze up to and including this time step

% Algo info
atype = {'C','H','O'}; % atom types which make up molecule
algo_atom = {'C'}; % atom which the molecule clustering algorithm is built around
cut = 4; % cutoff distance from algo_atom to exclude in search for atoms in that cluster (Å)
num_mol = 30; % number of molecules
atom_mol = 10; % number of atoms per molecule
algo = 'no'; % if yes will use built in algorithm to determine molecule clusters
save_mol_data = 'yes'; % set to yes if you want the algo to save cluster data
data_load = 'mol_data.mat'; % as mol_data.m so subsequent calculations do not
                            % have to run the algorithm. To use the stored
                            % data set algo = 'no'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thank you %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('###################################################################\n');
fprintf('##### INPUT INFO - START ##########################################\n');
fprintf('###################################################################\n\n');
fprintf('Description             (variable)   =                        input\n');
fprintf('===================================================================\n');
fprintf('input data file         (xdatfile)   = %+28s\n',xdatfile);
fprintf('time step               (ts)         = %25.1f fs\n',ts); start_fs = time_start*ts; end_fs = time_end*ts; sample_fs = sample*ts;
fprintf('starting time step      (time_start) = %4i ts units,\t %7.1f fs\n',time_start,start_fs);
fprintf('ending time step        (time_end)   = %4i ts units,\t %7.1f fs\n',time_end,end_fs);
fprintf('sampling time           (sample)     = %4i ts units,\t %7.1f fs\n',sample,sample_fs); tot_time_points = floor((time_end - time_start + 1) / sample);
fprintf('-------------------------------------------------------------------\n');
fprintf('                   total time points = %4i\n',tot_time_points);
fprintf('-------------------------------------------------------------------\n');
fprintf('Molecule Info\n');
fprintf('-------------------------------------------------------------------\n');
fprintf('number of molecules     (num_mol)    = %i\n',num_mol);
fprintf('atoms per molecule      (atom_mol)   = %i\n',atom_mol);
fprintf('atom types in molecule  (atype)      = '); fprintf(1,'%s ',atype{:});
fprintf('\nuse cluster algo?       (algo)       = %s\n',algo);
fprintf('===================================================================\n');


%% Getting Preliminary data from xdatfile %% 

fprintf('\n\ngetting preliminary info from "%s"\n\n',xdatfile);

% Header Data
fid = fopen(xdatfile);
fgetl(fid);
fgetl(fid);
% getting lattice vectors
a1 = fscanf(fid,'%f %f %f',3);
a2 = fscanf(fid,'%f %f %f',3);
a3 = fscanf(fid,'%f %f %f',3);
A = [a1 a2 a3]; clear a1 a2 a3;
% getting atom type and number
fgetl(fid);
line = fgetl(fid);
atoms = regexp(line, ' ', 'split'); atoms = atoms(~cellfun('isempty',atoms));
atom_nums = fscanf(fid,'%i',size(atoms,2));
tot_atoms = sum(atom_nums);

fprintf('Atom Information\n');
fprintf('===================================================================\n');
fprintf('Element\t#Atoms\tLine#s\n');
fprintf('-------------------------------------------------------------------\n');
counter=0;
num_atoms_to_keep = 0;
for i=1:length(atom_nums)
    if any(strcmp(atype,atoms{i}))
        num_atoms_to_keep = num_atoms_to_keep + atom_nums(i);
    end
end
atom_lines = zeros(length(atom_nums),2);
for i=1:length(atoms)
    initial=counter + 1;
    atom_lines(i,1) = initial;
    counter=counter+atom_nums(i);
    atom_lines(i,2) = counter;
    fprintf('%+6s:\t%6.0f\t%3.0f-%-3.0f\n', atoms{i},atom_nums(i),initial,counter);
end
lines_keep = zeros(num_atoms_to_keep,1);
counter = 1;
atoms_left = cell(length(atype),1);
atom_nums_left = zeros(length(atype),1);
for i=1:length(atom_nums)
    if any(strcmp(atype,atoms{i}))
        atoms_left{find(cellfun('isempty', atoms_left),1)} = atoms{i};
        atom_nums_left(find(~atom_nums_left,1)) = atom_nums(i);
        lines_keep(counter:counter+atom_nums(i)-1) = atom_lines(i,1):atom_lines(i,2);
        counter = counter + atom_nums(i);
    end
end
tot_atoms_left = sum(atom_nums_left);

fprintf('===================================================================\n\n');
fprintf('###################################################################\n');
fprintf('#####  INPUT INFO - END  ##########################################\n');
fprintf('###################################################################\n\n');

fprintf('\t...time elapsed: %f seconds\n\n',toc);

%% Molecule Cluster Input / Algorithm %%

if strcmp(algo,'yes')
    
    
    
    % Retreiving Position Data for molecule clustering algorithm (first ts)
    content=fileread(xdatfile);
    blockStarts = strfind(content, '=')+7;
    ts_found = length(blockStarts);

    % check if time_end input is valid
    if time_end > ts_found
        fprintf('=======================\n')
        fprintf('===== ERROR PANIC =====\n')
        fprintf('=======================\n')
        fprintf('The ending time step you have chosen is larger than the amount of steps found in "%s".\ntime_end: %i\nxdatfile steps: %i\nCheck variable: time_end.\n', xdatfile, ts, ts_found);
        clear content
        return
    end
    fprintf('###################################################################\n');
    fprintf('##### ALGO - START ################################################\n');
    fprintf('###################################################################\n\n');
    fprintf('You have chosen to use the built-in algorithm to determine molecule clusters.\n\n');
    fprintf('Description             (variable)   = input\n');
    fprintf('===================================================================\n');
    fprintf('atoms in molecule       (atype)      = {%s %s %s}\n',atype{1},atype{2},atype{3});
    fprintf('algo start atom         (algo_atom)  = {%s}\n',algo_atom{1});
    fprintf('cluster cutoff radius   (cut)        = %d\n',cut);
    fprintf('===================================================================\n\n');
    fprintf('reading initial configuration from %s...',xdatfile);
    
    
    blockEnds = strfind(content, 'Direct');

    % keeping only the first time step
    blockStart1 = blockStarts(1);
    blockEnd1 = blockEnds(2);

    clear blockStarts blockEnds ts_found;

    data = reshape(sscanf(content(blockStart1:blockEnd1), '%f'),3,[])';
    
    clear blockStart1 blockEnd1 content;

    % keeping atoms types in molecule based on user input
    data = data(lines_keep,:);
    
    fprintf('success.\n\n');
    fprintf('\t...time elapsed %f seconds\n\n',toc);
    
    % convert to cartesian
    fprintf('converting to cartesian coordinates...');
    data = (A*data')';
    fprintf('done \n\n');
    fprintf('\t...time elapsed %f seconds\n\n',toc);
    
    % adjacent cells to consider in distance calculation
    cls = [perms([0 0 0]); perms([1 0 0]); perms([-1 0 0]); perms([1 1 0]);...
        perms([1 -1 0]); perms([-1 -1 0]); perms([1 1 1]); perms([-1 1 1]);...
        perms([-1 -1 1]); perms([-1 -1 -1])];
    cls = unique(cls,'rows');
    cls = (A*cls')';
    
    % re-define atom type ranges
    %algo_atom_id = find(strcmp(atoms_left,algo_atom)); % Index of Carbon atoms
    %H_atom_id = find(strcmp(atoms_left,'H'));
    %O_atom_id = find(strcmp(atoms_left,'O'));
    
    % will find the general way to do this later
    C_range = [1 90]; 
    H_range = [91 90+120];
    O_range = [H_range(2)+1 H_range(2)+90]; 
    C = C_range(1):C_range(2);
    
    H = H_range(1):H_range(2);
    O = O_range(1):O_range(2);
    % changeables
    C_left = C;
    all_lns_left = (1:tot_atoms_left)';
    
    
    % goal is to cluster atom line numbers and associate them with a molecule
    % columns:
    %      1              2         3          4      5      6   7  8  9  10
    % O surrounded C, H bonded C, H bonded C, O end, O in, O in, H, H, H, H
    mol_data = zeros(num_mol,atom_mol);
    mol_count = 0;
    fprintf('algorithm starting...\n\n\t...time elapsed %f seconds\n\n',toc);
    
    % choose first algo atom (Carbon atom)
    fprintf('choosing first %s atom...\n',algo_atom{1});
    for i=1:length(C) % do not change anything about C
        if ~ismember(C(i),C_left)
            continue
        else % remove this C from list now
            all_lns_left = setdiff(all_lns_left,C(i));
        end
        % we're gonna find a new molecule
        mol_count = mol_count + 1;
        %remove current C atom from list
        %algo_atom_lns(algo_atom_lns==i)=0;
        %all_lns(all_lns==i)=0;
        %all_lns_loop = all_lns;
        %all_lns_loop(all_lns_loop==0) = [];
        % remove all candidates further than cutoff (cut) depends on
        % molecule
        mol_lns = all_lns_left;
        dists = zeros(length(mol_lns),1);
        for j=1:length(mol_lns)
            d_min = min_dist(data(C(i),:),data(mol_lns(j),:),cls);
            if  d_min > cut
                mol_lns(mol_lns==mol_lns(j))=0;
            else
                dists(j,1) = d_min;
            end
        end
        mol_lns(mol_lns==0)=[]; % now have lns of relevant atoms
        dists(dists==0)=[];
        cut_set = sortrows([mol_lns dists],2); % array of line number and distance to C(i) sorted from shortest to furthest distance
        
        
        % section for C(i) type 1 (3 nn are Oxygens)
        if ge(cut_set(1,1),O_range(1)) &&  le(cut_set(1,1),O_range(2)) && ge(cut_set(2,1),O_range(1)) &&  le(cut_set(2,1),O_range(2)) && ge(cut_set(3,1),O_range(1)) &&  le(cut_set(3,1),O_range(2))
            mol_data(mol_count,1) = C(i); % C down
            mol_data(mol_count,4) = cut_set(1,1); mol_lns(mol_lns==cut_set(1,1))=[]; % O type 1 down
            mol_data(mol_count,5) = cut_set(2,1); mol_lns(mol_lns==cut_set(2,1))=[]; % O type 2 down
            mol_data(mol_count,6) = cut_set(3,1); mol_lns(mol_lns==cut_set(3,1))=[]; % O type 2 down
            
            % if only 2 carbons remain add them to the mol_data
            CI = and(ge(mol_lns',C_range(1)),le(mol_lns',C_range(2)));
            if sum(CI) == 2
                givens = zeros(2,2);
                givens(:,1) = mol_lns(CI);
                % check which is closest to O5
                givens(1,2) = min_dist(data(mol_data(mol_count,5),:),data(givens(1,1),:),cls); % O5-gC
                givens(2,2) = min_dist(data(mol_data(mol_count,5),:),data(givens(2,1),:),cls); % O5-gC
                givens = sortrows(givens,2); % put closest gC in slot 2 corresponds to O in slot 5, farther gC goes in slot 3 corresopnds to O in slot 6
                mol_data(mol_count,2) = givens(1,1); mol_lns(mol_lns==givens(1,1))=[]; % C down
                mol_data(mol_count,3) = givens(2,1); mol_lns(mol_lns==givens(2,1))=[]; % C down
                clear givens
            else % find C closest to each O type 2. should weed out all but C's within unk_set before
                O2_dists = zeros(length(mol_lns),2);
                O2_dists(:,1) = mol_lns;
                for l=1:length(mol_lns)
                    O2_dists(l,2) = min_dist(data(mol_data(mol_count,5),:),data(mol_lns(l),:),cls);
                end
                O2_dists = sortrows(O2_dists,2); % the first element will be the C type 2 atom closest to this O type 2
                mol_data(mol_count,2) = O2_dists(1,1);
                
                
                mol_lns(mol_lns==O2_dists(1,1))=[];
                
                
                O3_dists = zeros(length(mol_lns),2);
                O3_dists(:,1) = mol_lns;
                for l=1:length(mol_lns)
                    O3_dists(l,2) = min_dist(data(mol_data(mol_count,6),:),data(mol_lns(l),:),cls);
                end
                O3_dists = sortrows(O3_dists,2); % the first element will be the C type 2 atom closest to this O type 2
                mol_data(mol_count,3) = O3_dists(1,1);
                mol_lns(mol_lns==O3_dists(1,1))=[];
            end
            % if only 4 hydrogen remain add them to mol data
            HI = and(ge(mol_lns',H_range(1)), le(mol_lns',H_range(2)));
            if sum(HI) == 4
                givens = mol_lns(HI);
                mol_data(mol_count,7) = givens(1); mol_lns(mol_lns==givens(1))=[]; % H down
                mol_data(mol_count,8) = givens(2); mol_lns(mol_lns==givens(2))=[]; % H down
                mol_data(mol_count,9) = givens(3); mol_lns(mol_lns==givens(3))=[]; % H down
                mol_data(mol_count,10) = givens(4); mol_lns(mol_lns==givens(4))=[]; % H down
            else % find 2 closest H's to the 2 C type 2's you just found
                C2_dists = zeros(length(mol_lns),2);
                C2_dists(:,1) = mol_lns;
                for l=1:length(mol_lns)
                    C2_dists(l,2) = min_dist(data(mol_data(mol_count,2),:),data(mol_lns(l),:),cls);
                end
                C2_dists = sortrows(C2_dists,2); % the first 2 elements will be the H's
                mol_data(mol_count,7) = C2_dists(1,1);
                mol_data(mol_count,8) = C2_dists(2,1);
                
                mol_lns(mol_lns==C2_dists(1,1))=[];
                mol_lns(mol_lns==C2_dists(2,1))=[];
                
                C3_dists = zeros(length(mol_lns),2);
                C3_dists(:,1) = mol_lns;
                for l=1:length(mol_lns)
                    C3_dists(l,2) = min_dist(data(mol_data(mol_count,3),:),data(mol_lns(l),:),cls);
                end
                C3_dists = sortrows(C3_dists,2); % the first 2 elements will be the H's
                mol_data(mol_count,9) = C3_dists(1,1);
                mol_data(mol_count,10) = C3_dists(2,1);
            end
        % section for C(i) type 2 (2nn are Hydrogens)      
        elseif ge(cut_set(1,1),H_range(1)) &&  le(cut_set(1,1),H_range(2)) && ge(cut_set(2,1),H_range(1)) &&  le(cut_set(2,1),H_range(2))
            mol_data(mol_count,2) = C(i); % C down
            mol_data(mol_count,7) = cut_set(1,1); mol_lns(mol_lns==cut_set(1,1))=[]; % H down
            mol_data(mol_count,8) = cut_set(2,1); mol_lns(mol_lns==cut_set(2,1))=[]; % H down
            if ge(cut_set(3,1),C_range(1)) &&  le(cut_set(3,1),C_range(2))
                mol_data(mol_count,3) = cut_set(3,1); % C down
                mol_data(mol_count,5) = cut_set(4,1); % O down
            elseif  ge(cut_set(3,1),O_range(1)) &&  le(cut_set(3,1),O_range(2))
                mol_data(mol_count,5) = cut_set(3,1); % O down
                mol_data(mol_count,3) = cut_set(4,1); % C down
            end
            mol_lns(mol_lns==cut_set(3,1))=[];
            mol_lns(mol_lns==cut_set(4,1))=[];
            % if only 2 hydrogen remain add them to mol data
            HI = and(ge(mol_lns,H_range(1)), le(mol_lns,H_range(2)));
            if sum(HI) == 2
                givens = mol_lns(HI);
                mol_data(mol_count,9) = givens(1); mol_lns(mol_lns==givens(1))=[]; % H down
                mol_data(mol_count,10) = givens(2); mol_lns(mol_lns==givens(2))=[]; % H down
            else % find 2 closest hydrogens to C type 2 mol_data(mol_count,3)
                C3_dists = zeros(length(mol_lns),2);
                C3_dists(:,1) = mol_lns;
                for l=1:length(mol_lns)
                    C3_dists(l,2) = min_dist(data(mol_data(mol_count,3),:),data(mol_lns(l),:),cls);
                end
                C3_dists = sortrows(C3_dists,2); % the first 2 elements will be the H's
                mol_data(mol_count,9) = C3_dists(1,1);
                mol_data(mol_count,10) = C3_dists(2,1);
                
                mol_lns(mol_lns==C3_dists(1,1))=[];
                mol_lns(mol_lns==C3_dists(2,1))=[];
            end
            % just recalculate those distances and the closest should be
            % the other O type 2
            C3_dists = zeros(length(mol_lns),2);
            C3_dists(:,1) = mol_lns;
            for l=1:length(mol_lns)
                C3_dists(l,2) = min_dist(data(mol_data(mol_count,3),:),data(mol_lns(l),:),cls);
            end
            C3_dists = sortrows(C3_dists,2); % the first element will be the remaining O type 2
            mol_data(mol_count,6) = C3_dists(1,1);
            
            mol_lns(mol_lns==C3_dists(1,1))=[];
            
            % now getting C of type 1
            O6_dists = zeros(length(mol_lns),2);
            O6_dists(:,1) = mol_lns;
            for l=1:length(mol_lns)
                O6_dists(l,2) = min_dist(data(mol_data(mol_count,6),:),data(mol_lns(l),:),cls);
            end
            O6_dists = sortrows(O6_dists,2); % the first element will be the C type 1
            mol_data(mol_count,1) = O6_dists(1,1);
            
            mol_lns(mol_lns==O6_dists(1,1))=[];
            
            % finally the last O type 1
            C1_dists = zeros(length(mol_lns),2);
            C1_dists(:,1) = mol_lns;
            for l=1:length(mol_lns)
                C1_dists(l,2) = min_dist(data(mol_data(mol_count,1),:),data(mol_lns(l),:),cls);
            end
            C1_dists = sortrows(C1_dists,2); % the first element will be the O type 1
            mol_data(mol_count,4) = C1_dists(1,1);
        end
        % rm the ten atoms making up molecule
        all_lns_left = setdiff(all_lns_left,mol_data(mol_count,:));
        C_left = setdiff(C_left,mol_data(mol_count,1:3));
        if mol_count == num_mol
            break;
        end
    end
    
    % checking results
    if ~isequal(size(mol_data),[num_mol atom_mol])
        fprintf('=======================\n')
        fprintf('===== ERROR PANIC =====\n')
        fprintf('=======================\n')
        if num_mol*atom_mol ~= sum(atom_nums_left)
            fprintf('Check variables:\nnum_mol, atom_mol\nnum_mol * atom_mol should add up to the total atoms you are providing.')
        else
            fprintf('Something weird has happened. Time to debug algorithm.')
        end
        return
    end
    fprintf('algorithm completed.\n\n\t...time elapsed: %f seconds\n\n%i clusters of %i atoms found from specified atoms.\nchecking results...\n',toc,num_mol,atom_mol)
    all_used = all_lns_left == [];
    C_used = C_left == [];
    c_data = mol_data(:,1:3); c_data = c_data(:);
    c_check1 = all(and(ge(c_data,C_range(1)),le(c_data,C_range(2))));
    c_check2 = length(c_data) == length(unique(c_data));
    o_data = mol_data(:,4:6); o_data = o_data(:);
    o_check1 = all(and(ge(o_data,O_range(1)),le(o_data,O_range(2))));
    o_check2 = length(o_data) == length(unique(o_data));
    h_data = mol_data(:,7:10); h_data = h_data(:);
    h_check1 = all(and(ge(h_data,H_range(1)),le(h_data,H_range(2))));
    h_check2 = length(h_data) == length(unique(h_data));
    if all([C_used all_used c_check1 c_check2 o_check1 o_check2 h_check1 h_check2])
        fprintf('Wow it appears to have worked, good job.\n\n');
        fprintf('\t...time elapsed: %f seconds\n\n',toc);
    else
        fprintf('=======================\n')
        fprintf('===== ERROR PANIC =====\n')
        fprintf('=======================\n')
        fprintf('Either an atom has been clustered more than once,\nor the wrong atom was assigned to a molecular position.\nTime to debug algorithm.')
        return
    end
    
    if strcmp(save_mol_data,'yes')
        save('mol_data.mat','mol_data','num_mol','atom_mol','lines_keep','atype');
    end
    
    
    clear all_lns_left c_data o_data h_data all_used c_check1 c_check2...
        o_check1 o_check2 h_check1 h_check2 cut C1_dists C2_dists C3_dists...
        C_used C_left CI counter d_min dists fid givens HI i j line mol_count...
        O2_dists O3_dists O6_dists cut_set data C H O;
    
    fprintf('###################################################################\n');
    fprintf('#####  ALGO - END  ################################################\n');
    fprintf('###################################################################\n\n');
    fprintf('\t...time elapsed: %f seconds\n\n',toc);
elseif strcmp(algo,'no')
    fprintf('You have opted to provide the molecular cluster data manually\n\tbased on the setting of the variable "algo".\n\n...attempting to load cluster data\n\n')
        
    % the molecule cluster array should have dimension number of molecules
    % by number of atoms in the molecule by 1 (containing atom line number)
    % mol_data = zeros(num_mol,atom_mol);
    % enter data below
    
    load(data_load);
    
    % checking the appropriate data has been provided
    if ~exist('mol_data','var')
        fprintf('=======================\n')
        fprintf('===== ERROR PANIC =====\n')
        fprintf('=======================\n')
        fprintf('Did not find variable ''mol_data'' cannot continue.\n')
        return
    elseif ~ndims(mol_data) == 2
        fprintf('=======================\n')
        fprintf('===== ERROR PANIC =====\n')
        fprintf('=======================\n')
        fprintf('Variable ''mol_data'' does not have correct dimension/nIt should be a 2-D integer array with size [# of molecules x # atoms/molecule].\n')
        return
    elseif ~exist('lines_keep','var')
        fprintf('=======================\n')
        fprintf('===== ERROR PANIC =====\n')
        fprintf('=======================\n')
        fprintf('Did not find variable ''lines_keep'' cannot comtinue.\nEither provide a array of length [number of lines in molecules],\nor fill in the variable ''atype'' in the initial user input section\n')
        fprintf('The variable ''lines_keep'' is an array containing which atom lines (those which make up molecules) to keep from the XDATCAR data.\nThis variable will be built from the variable ''atype'' if provided prior to loading data\n')
        return
    end
    
    num_mol = size(mol_data,1);
    atom_mol = size(mol_data,2);
    if num_mol*atom_mol == length(lines_keep)
        fprintf('...cluster data successfully loaded.\n\n')
        fprintf('...found %i molecules with %i atoms each.\n\n',num_mol,atom_mol)
    else
        fprintf('=======================\n')
        fprintf('===== ERROR PANIC =====\n')
        fprintf('=======================\n')
        fprintf('The number of molecules * the number of atoms per molecule does not match the total number of atoms you wish to analyze!.\n\n')
        return
    end
else
    fprintf('=======================\n')
    fprintf('===== ERROR PANIC =====\n')
    fprintf('=======================\n')
    fprintf('Input string for variable ''algo'' not recognized.\n\n')
end

%% Molecular Analysis %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% User Input Required %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bond distances of interest and their cutoffs use graphic below for
% guidance.

%bond#  eq.
%#1     1.33               #2
C1_O5 = 1.80;              C1_O6 = C1_O5; % C type 1 - O type 2's
%#3     1.42               #4  
C2_O5 = 1.80;              C3_O6 = C2_O5; % C type 2 - O type 2's
%#5     1.53
C2_C3 = 1.90; % C type 2 - C type 2

% create array of maximum bond lengths from above, the order matters here
% one can change the code from here if changes are made to the bonds and
% ordering
bond_max = [C1_O5 C1_O6 C2_O5 C3_O6 C2_C3];

% also want to know which electrode (left of right) the molecule broke on,
% below. Provide a z coordinate for which if the z coordinate of the C 
% atom in the breaking bond of the molecule is less than the coordinate is
% is assigned to the left electrode and if greater assigned to the right
% electrode. This could also be chosen as half of the min and max z coordinate
% of the atoms in all molecules. At this point in the code that information
% is not stored yet. For now it seems easy enough to enter manually a good
% enough guess.
z_boundary = 20;

% data to store
mol_bonds = 0; % saves all bond lengths if == 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      %                                %
%          O4          %                                %
%          ||          %   bond#  atoms  equ.  break    %
%          C1          %                                %
%        /    \        %     #1   C1-O5  1.33  >1.80    %
%      O5      O6      %     #2   C1-O6  1.33  >1.80    %
%       \      /       %     #3   C2-O5  1.42  >1.80    %
%        C2--C3        %     #4   C2-O6  1.42  >1.80    %
%       / \  / \       %     #5   C2-C3  1.53  >1.90    %
%      H7 H8H9 H10     %                                %
%                      %                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Ethylene Carbonate  %  Table of conditions for bond  %
%       Molecule       %  breaking in a molecule        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thank you %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('###################################################################\n');
fprintf('#####  MOLECULE ANALYSIS - START  #################################\n');
fprintf('###################################################################\n\n');
fprintf('\t...time elapsed: %f seconds\n\n',toc);

fprintf('Molecule Breakdown\n');

fprintf('|=======================================================|\n');
fprintf('|                      |                                |\n');
fprintf('|          O4          |                                |\n');
fprintf('|          ||          |   bond#  atoms  equ.  break    |\n');
fprintf('|          C1          |                                |\n');
fprintf('|        /    \\        |     #1   C1-O5  1.33  >%4.2f    |\n',C1_O5);
fprintf('|      O5      O6      |     #2   C1-O6  1.33  >%4.2f    |\n',C1_O6);
fprintf('|       \\      /       |     #3   C2-O5  1.42  >%4.2f    |\n',C2_O5);
fprintf('|        C2--C3        |     #4   C3-O6  1.42  >%4.2f    |\n',C3_O6);
fprintf('|       / \\  / \\       |     #5   C2-C3  1.53  >%4.2f    |\n',C2_C3);
fprintf('|      H7 H8H9 H10     |                                |\n');
fprintf('|                      |                                |\n');
fprintf('|=======================================================|\n');
fprintf('|  Ethylene Carbonate  |  Table of conditions for bond  |\n');
fprintf('|       Molecule       |  breaking in a molecule        |\n');
fprintf('|=======================================================|\n\n\n');



% Retreiving Position Data
content=fileread(xdatfile);

blockStarts = strfind(content, '=')+7;
ts_found = length(blockStarts);


% check if time_end input is valid
if time_end > ts_found
    fprintf('=======================\n')
    fprintf('===== ERROR PANIC =====\n')
    fprintf('=======================\n')
    fprintf('The ending time step you have chosen is larger than\nthe amount of steps found in "%s".\ntime_end: %i\nxdatfile steps: %i\nCheck variable: time_end.\n', xdatfile, ts, ts_found);
    clear content
    return
end

% reimplementing: adjacent cells to consider in distance calculation
cls = [perms([0 0 0]); perms([1 0 0]); perms([-1 0 0]); perms([1 1 0]);...
    perms([1 -1 0]); perms([-1 -1 0]); perms([1 1 1]); perms([-1 1 1]);...
    perms([-1 -1 1]); perms([-1 -1 -1])];
cls = unique(cls,'rows');
cls = (A*cls')';




fprintf('...storing position data of {%s %s %s} atoms\n\tfrom %4.1f fs to %4.1f fs every %4.1f fs\n',atype{1},atype{2},atype{3},start_fs,end_fs,sample_fs);


blockEnds = strfind(content, 'Direct');
blockEndsLast = length(content);
blockEnds = [blockEnds(2:end) blockEndsLast];

% restricting to user input request
blockStarts =  blockStarts(time_start:sample:time_end);
blockEnds =  blockEnds(time_start:sample:time_end);

nBlocks = numel(blockStarts);

data = zeros(nBlocks,tot_atoms,3);

for i=1:nBlocks
    data(i,:,:) = reshape(sscanf(content(blockStarts(i):blockEnds(i)), '%f'),3,[])';
end

if mol_bonds == 0
    fprintf('\n...expect array size of\n\t%i x %i x %i\n\t(time slices x atoms x monitored bonds)\n\n',nBlocks,length(lines_keep),length(bond_max));
elseif mol_bonds == 1
    fprintf('\n...expect two array sizes of\n\t%i x %i x %i\n\t(time slices x atoms x monitored bonds)\n\n',nBlocks,length(lines_keep),length(bond_max));
end

clear content blockStarts blockEnds blockEndsLast i ts_found;

% keeping only specified atoms
data = data(:,lines_keep,:);

clear lines_keep;

fprintf('\t...time elapsed: %f seconds\n\n',toc);
fprintf('...converting to cartesian coordinates\n\n');
for i=1:nBlocks
    data(i,:,:) = (A*squeeze(data(i,:,:))')';
end

fprintf('\t...time elapsed: %f seconds\n\n',toc);
fprintf('...rearranging data\n\n');

% rearranging atoms into molecule clusters
mol_time_data = zeros(nBlocks,num_mol,atom_mol,3);
for i=1:nBlocks
    for j=1:num_mol
        mol_time_data(i,j,:,:) = data(i,mol_data(j,:),:);
    end 
end
clear data
fprintf('\t...time elapsed: %f seconds\n\n',toc);

% mol_time_data now stores coordinates in the fashion below
% time 1:
    % mol 1:
        % atom 1: x y z
        % atom 2: x y z
        % ...
        % atom 10: x y z
    % mol 2:
        % atom 1: x y z
        % ...
if mol_bonds == 0
    fprintf('...calculating bond distances to determine breaks\n\n');
elseif mol_bonds == 1
    fprintf('...calculating bond distances to determine breaks\n...storing bind distances\n\n');
end
        
 
% OPTION 1: section that stores what bonds break and when
if mol_bonds == 0
    bond_breaks = zeros(nBlocks,num_mol,length(bond_max)); % set to 1 if broken
    for i=1:nBlocks
        for j=1:num_mol
            bond_breaks(i,j,:) = ...
                [min_dist(squeeze(mol_time_data(i,j,1,:))',squeeze(mol_time_data(i,j,5,:))',cls)...
                min_dist(squeeze(mol_time_data(i,j,1,:))',squeeze(mol_time_data(i,j,6,:))',cls)...
                min_dist(squeeze(mol_time_data(i,j,2,:))',squeeze(mol_time_data(i,j,5,:))',cls)...
                min_dist(squeeze(mol_time_data(i,j,3,:))',squeeze(mol_time_data(i,j,6,:))',cls)...
                min_dist(squeeze(mol_time_data(i,j,2,:))',squeeze(mol_time_data(i,j,3,:))',cls)]...
                > bond_max;
        end
    end
end

% OPTION 2: section that stores what bonds break, when, and the bond distances
if mol_bonds == 1
    bond_lengths = zeros(nBlocks,num_mol,length(bond_max)); % 5 bond lengths of interest
    bond_breaks = zeros(nBlocks,num_mol,length(bond_max)); % set to 1 if broken
    for i=1:nBlocks
        for j=1:num_mol
            bond_lengths(i,j,1) = min_dist(squeeze(mol_time_data(i,j,1,:))',squeeze(mol_time_data(i,j,5,:))',cls);
            bond_lengths(i,j,2) = min_dist(squeeze(mol_time_data(i,j,1,:))',squeeze(mol_time_data(i,j,6,:))',cls);
            bond_lengths(i,j,3) = min_dist(squeeze(mol_time_data(i,j,2,:))',squeeze(mol_time_data(i,j,5,:))',cls);
            bond_lengths(i,j,4) = min_dist(squeeze(mol_time_data(i,j,3,:))',squeeze(mol_time_data(i,j,6,:))',cls);
            bond_lengths(i,j,5) = min_dist(squeeze(mol_time_data(i,j,2,:))',squeeze(mol_time_data(i,j,3,:))',cls);
            bond_breaks(i,j,:) = squeeze(bond_lengths(i,j,:)) > bond_max;
        end
    end
end

fprintf('...calculations complete\n\n');
fprintf('\t...time elapsed: %f seconds\n\n',toc);

% for output purposes
fprintf('Molecule breaking summary\n');
fprintf('===================================================================\n');
mol_left = 1:num_mol;
mol_break_plot = zeros(nBlocks,1);
mol_break_plot_L = zeros(nBlocks,1);
mol_break_plot_R = zeros(nBlocks,1);
broke_pos = zeros(50,4); % seems like enough
for i=1:nBlocks
    for j=1:length(mol_left)
        if any(bond_breaks(i,mol_left(j),:)) == 1
            bond = find(bond_breaks(i,mol_left(j),:),1);
            fprintf('Molecule %3i broke at %5i fs via bond %2i\n',mol_left(j),start_fs-1 + i*sample_fs,bond);
            mol_break_plot(i:end) = mol_break_plot(i:end) + 1;
            if bond == 1 || 2
                broke_pos(mol_break_plot(end),:) = [squeeze(mol_time_data(i,mol_left(j),1,:)); bond];
            elseif bond == 3 || 5
                broke_pos(mol_break_plot(end),:) = [squeeze(mol_time_data(i,mol_left(j),2,:)); bond];
            elseif bond == 4
                broke_pos(mol_break_plot(end),:) = [squeeze(mol_time_data(i,mol_left(j),3,:)); bond];
            end
            if broke_pos(mol_break_plot(end),3) < z_boundary
                mol_break_plot_L(i:end) = mol_break_plot_L(i:end) + 1;
            elseif broke_pos(mol_break_plot(end),3) > z_boundary
                mol_break_plot_R(i:end) = mol_break_plot_R(i:end) + 1;
            end
            mol_left(j) = 0;
        end
    end
    mol_left(mol_left==0) = [];
end
broke_pos = broke_pos(any(broke_pos,2),:);
fprintf('-------------------------------------------------------------------\n');
fprintf('\tTotal Molecules Broken: %3i\n',max(mol_break_plot));
fprintf('===================================================================\n');

fprintf('\n...molecular analysis complete\n\n');


fprintf('Output data stored\n\n');
fprintf('Variable\t Data Type\n');
fprintf('===================================================================\n');
fprintf('bond_breaks\t %ix%ix%i logical array\n',nBlocks,num_mol,length(bond_max));
fprintf('-------------------------------------------------------------------\n');
fprintf('Description:\n');
fprintf('[time slices x number molecules x number of bonds monitored]\n');
fprintf('An entry of ''1'' signifies a broken bond.\n');
fprintf('===================================================================\n');
fprintf('mol_break_plot\t %ix1 integer vector\n',nBlocks);
fprintf('-------------------------------------------------------------------\n');
fprintf('Description:\n');
fprintf('[time slices x 1]\n');
fprintf('Contains the cumulative number of broken molecules up to that\n');
fprintf('time slice. Convenient for 2D plot versus time.\n');
fprintf('===================================================================\n');
fprintf('mol_break_plot_L\t %ix1 integer vector\n',nBlocks);
fprintf('-------------------------------------------------------------------\n');
fprintf('Description:\n');
fprintf('[time slices x 1]\n');
fprintf('The same data as the mol_break_plot variable but only on the left\n');
fprintf('electrode. Determined using user specified variable ''z_boundary''.\n');
fprintf('===================================================================\n');
fprintf('mol_break_plot_R\t %ix1 integer vector\n',nBlocks);
fprintf('-------------------------------------------------------------------\n');
fprintf('Description:\n');
fprintf('[time slices x 1]\n');
fprintf('The same data as the mol_break_plot variable but only on the right\n');
fprintf('electrode. Determined using user specified variable ''z_boundary''.\n');
fprintf('===================================================================\n');
fprintf('broke_pos\t %ix4 double array\n',size(broke_pos,1));
fprintf('-------------------------------------------------------------------\n');
fprintf('Description:\n');
fprintf('[number of broken molecules x spatial dimension+label broken bond]\n');
fprintf('Contains the x,y,z coordinates (cartesian) of the C atom in the\n');
fprintf('broken bond of the molecule at the first time slice the molecule\n');
fprintf('broke. Element 4 in the array contains the bond type that has\nbroken (1-5). Can use to plot histogram plot along z.\n');
if mol_bonds == 1
    fprintf('===================================================================\n');
    fprintf('bond_length\t %ix%ix%i double array\n',nBlocks,num_mol,mength(bond_max));
    fprintf('-------------------------------------------------------------------\n');
    fprintf('Description:\n');
    fprintf('[time slices x number molecules x number of bonds monitored]\n');
    fprintf('Contains the bond lengths at all time slices of all monitored bonds.\n');
end
fprintf('===================================================================\n\n\n');


clearvars -except bond_breaks mol_break_plot broke_pos bond_length A mol_break_plot_L mol_break_plot_R


fprintf('###################################################################\n');
fprintf('#####  MOLECULE ANALYSIS - END  ###################################\n');
fprintf('###################################################################\n\n');
fprintf('\t...time elapsed: %f seconds\n\n\n\n',toc);


fprintf('###########################################################################\n')
fprintf('###########################################################################\n')
fprintf('##########                                                   ##############\n')
fprintf('##########   ################################################   ###########\n')
fprintf('###########   #                                            #####   ########\n')
fprintf('###########   #                                              #####   ######\n')
fprintf('############   #                                             #######   ####\n')
fprintf('############   #                                             ########   ###\n')
fprintf('#############   ################    ############            #########   ###\n')
fprintf('###############    #######         ##########              #########   ####\n')
fprintf('##################   #            #######                 #########   #####\n')
fprintf('#################   #            ####                    #########   ######\n')
fprintf('################   #            #                     ###########   #######\n')
fprintf('###############   #                                ###########   ##########\n')
fprintf('##############   #                              ###########   #############\n')
fprintf('#############   #             #                    ####   #################\n')
fprintf('############   #             #####                  ######   ##############\n')
fprintf('###########   #             #######                 ########   ############\n')
fprintf('##########   #             #########               ##########   ###########\n')
fprintf('#########   #             ##########              ###########   ###########\n')
fprintf('########   #             #########               ###########   ############\n')
fprintf('#######   #                                     ###########   #############\n')
fprintf('######   #                                     ###########   ##############\n')
fprintf('#####   #                                     ###########   ###############\n')
fprintf('####   #                                   #############   ################\n')
fprintf('###   #################################################   #################\n')
fprintf('####   ##############################################   ###################\n')
fprintf('#####   ############################################   ####################\n')
fprintf('######   #########################################   ######################\n')
fprintf('#######                                            ########################\n')
fprintf('###########################################################################\n')
fprintf('###########################################################################\n\n\n')









%% Define in-script functions %%

function a = min_dist(p1,p2,adj)
% returns minimum distance between two atoms positioned at p1 and p2 within 
% the unit cell where p1 and p2 are 3x1 arrays of their cartesian
% coordinates. The variable adj should be a 3xN (N = 27 n.n. cells) array
% consisting of different Real space lattice vecter in cartesian coordinates.
% p2 is translated by all R's provided in adj and the distance between those
% positions and p1 are calculated. This function returns the minimum distance
% of the set of distances.
    a = min(vecnorm(repmat(p2,length(adj),1) + adj - repmat(p1,length(adj),1),2,2));
end