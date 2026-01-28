clear; clc; close all;
file_name = 'iJO1366.mat';
load(file_name);
disp(['Successfully loaded file: ', file_name]);
disp('------------------------------');

if exist('Model', 'var')
    stoichiometric_matrix = Model.US;
    [num_nodes, num_total_reactions] = size(stoichiometric_matrix);
    
    fprintf('Parsing directed hypergraph from stoichiometric matrix...\n');
    fprintf('  - Number of nodes (metabolites): %d\n', num_nodes);
    fprintf('  - Total reactions: %d\n', num_total_reactions);
    disp('------------------------------');
    
    directed_hyperedges = struct('inputs', {}, 'outputs', {});
    
    for j = 1:num_total_reactions
        current_reaction_vector = stoichiometric_matrix(:, j);
        input_nodes = find(current_reaction_vector < 0);
        output_nodes = find(current_reaction_vector > 0);
        
        is_self_loop = ~isempty(intersect(input_nodes, output_nodes));
        
        if ~isempty(input_nodes) && ~isempty(output_nodes) && ~is_self_loop
            idx = length(directed_hyperedges) + 1;
            directed_hyperedges(idx).inputs = input_nodes;
            directed_hyperedges(idx).outputs = output_nodes;
        end
    end
    
    num_valid_hyperedges = length(directed_hyperedges);
    disp('Hypergraph structure successfully filtered and converted.');
    fprintf('Number of valid hyperedges (reactions): %d\n', num_valid_hyperedges);
    disp('------------------------------');
    
    output_filename = 'iJO1366_hypergraph.txt';
    fprintf('\n--- Writing hypergraph to file : %s ---\n', output_filename);
    try
        fileID = fopen(output_filename, 'wt');
        fprintf(fileID, '%d %d\n', num_nodes, num_valid_hyperedges);
        
        for i = 1:num_valid_hyperedges
            he = directed_hyperedges(i);
            m_in = numel(he.inputs);
            m_out = numel(he.outputs);
            
            fprintf(fileID, '%d %d\n', m_in, m_out);
            
            fprintf(fileID, '%d ', he.inputs - 1);
            fprintf(fileID, '\n');
            
            fprintf(fileID, '%d ', he.outputs - 1);
            fprintf(fileID, '\n');
        end
        
        fclose(fileID);
        disp('File written successfully!');
        
    catch ME
        if exist('fileID', 'var') && fileID ~= -1
            fclose(fileID);
        end
        error('Error writing file: %s', ME.message);
    end
    
else
    warning("Warning: 'Model' variable not found in .mat file.");
end
