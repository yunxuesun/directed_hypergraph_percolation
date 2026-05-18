#include "globals.h"

void load_hypergraph_from_file(const char* filepath) {
    FILE *file = fopen(filepath, "r");
    if (!file) {
        perror("Error opening hypergraph file");
        exit(EXIT_FAILURE);
    }

    if (fscanf(file, "%d %d\n", &N, &Q) != 2) {
        fprintf(stderr, "Error: Could not read N and Q from file.\n");
        exit(EXIT_FAILURE);
    }
    printf("Loading graph from file: N=%d nodes, Q=%d hyperedges.\n", N, Q);

    hyperedge_T = (int**)calloc(Q, sizeof(int*));
    hyperedge_H = (int**)calloc(Q, sizeof(int*));
    size_T = (int*)calloc(Q, sizeof(int));
    size_H = (int*)calloc(Q, sizeof(int));

    for (int mu = 0; mu < Q; mu++) {
        if (fscanf(file, "%d %d\n", &size_T[mu], &size_H[mu]) != 2) exit(EXIT_FAILURE);
        
        hyperedge_T[mu] = (int*)calloc(size_T[mu], sizeof(int));
        hyperedge_H[mu] = (int*)calloc(size_H[mu], sizeof(int));

        for (int n = 0; n < size_T[mu]; n++) {
            if (fscanf(file, "%d", &hyperedge_T[mu][n]) != 1) exit(EXIT_FAILURE);
        }
        for (int n = 0; n < size_H[mu]; n++) {
            if (fscanf(file, "%d", &hyperedge_H[mu][n]) != 1) exit(EXIT_FAILURE);
        }
    }
    fclose(file);
    printf("Hypergraph data loaded successfully.\n");
}

void generate_synthetic_hypergraph() {
    N = N_SYNTHETIC;
    Q = Q_SYNTHETIC;

    hyperedge_T = malloc(Q * sizeof(int*));
    hyperedge_H = malloc(Q * sizeof(int*));
    size_T = malloc(Q * sizeof(int));
    size_H = malloc(Q * sizeof(int));

    for (int mu = 0; mu < Q; mu++) {
        size_T[mu] = M_IN_SYNTHETIC;
        size_H[mu] = M_OUT_SYNTHETIC;
        hyperedge_T[mu] = malloc(size_T[mu] * sizeof(int));
        hyperedge_H[mu] = malloc(size_H[mu] * sizeof(int));
        
        for (int t = 0; t < size_T[mu]; t++) hyperedge_T[mu][t] = rand() % N;
        for (int h = 0; h < size_H[mu]; h++) hyperedge_H[mu][h] = rand() % N;
    }
    printf("Synthetic hypergraph generated: N=%d, Q=%d.\n", N, Q);
}

void build_mappings() {
    // Build a directed mapping
    k_in = calloc(N, sizeof(int));
    k_out = calloc(N, sizeof(int));

    for (int mu = 0; mu < Q; mu++) {
        for (int t = 0; t < size_T[mu]; t++) k_out[hyperedge_T[mu][t]]++;
        for (int h = 0; h < size_H[mu]; h++) k_in[hyperedge_H[mu][h]]++;
    }

    node_as_input_of = malloc(N * sizeof(int*));
    node_as_input_of_pos = malloc(N * sizeof(int*));
    node_as_output_of = malloc(N * sizeof(int*));
    node_as_output_of_pos = malloc(N * sizeof(int*));

    for (int i = 0; i < N; i++) {
        node_as_input_of[i] = k_out[i] > 0 ? malloc(k_out[i] * sizeof(int)) : NULL;
        node_as_input_of_pos[i] = k_out[i] > 0 ? malloc(k_out[i] * sizeof(int)) : NULL;
        node_as_output_of[i] = k_in[i] > 0 ? malloc(k_in[i] * sizeof(int)) : NULL;
        node_as_output_of_pos[i] = k_in[i] > 0 ? malloc(k_in[i] * sizeof(int)) : NULL;
    }
    
    int *cur_in_idx = calloc(N, sizeof(int));
    int *cur_out_idx = calloc(N, sizeof(int));
    for (int mu = 0; mu < Q; mu++) {
        for (int t = 0; t < size_T[mu]; t++) {
            int node = hyperedge_T[mu][t];
            int idx = cur_out_idx[node]++;
            node_as_input_of[node][idx] = mu;
            node_as_input_of_pos[node][idx] = t;
        }
        for (int h = 0; h < size_H[mu]; h++) {
            int node = hyperedge_H[mu][h];
            int idx = cur_in_idx[node]++;
            node_as_output_of[node][idx] = mu;
            node_as_output_of_pos[node][idx] = h;
        }
    }
    free(cur_in_idx); free(cur_out_idx);

    // Build an undirected mapping
    hyperedge_all = malloc(Q * sizeof(int*));
    size_all = malloc(Q * sizeof(int));
    k_all = calloc(N, sizeof(int));

    for (int mu = 0; mu < Q; mu++) {
        int *temp = malloc((size_T[mu] + size_H[mu]) * sizeof(int));
        int cnt = 0;
        
        for (int i = 0; i < size_T[mu]; i++) {
            int n = hyperedge_T[mu][i];
            bool dup = false;
            for(int j=0; j<cnt; j++) if(temp[j]==n) {dup=true; break;}
            if(!dup) temp[cnt++] = n;
        }
        for (int i = 0; i < size_H[mu]; i++) {
            int n = hyperedge_H[mu][i];
            bool dup = false;
            for(int j=0; j<cnt; j++) if(temp[j]==n) {dup=true; break;}
            if(!dup) temp[cnt++] = n;
        }
        size_all[mu] = cnt;
        hyperedge_all[mu] = malloc(cnt * sizeof(int));
        memcpy(hyperedge_all[mu], temp, cnt * sizeof(int));
        free(temp);

        for (int i = 0; i < cnt; i++) k_all[hyperedge_all[mu][i]]++;
    }

    node_in_hyperedge = malloc(N * sizeof(int*));
    node_in_hyperedge_pos = malloc(N * sizeof(int*));
    int *cur_all_idx = calloc(N, sizeof(int));
    
    for (int i = 0; i < N; i++) {
        node_in_hyperedge[i] = k_all[i] > 0 ? malloc(k_all[i] * sizeof(int)) : NULL;
        node_in_hyperedge_pos[i] = k_all[i] > 0 ? malloc(k_all[i] * sizeof(int)) : NULL;
    }

    for (int mu = 0; mu < Q; mu++) {
        for (int pos = 0; pos < size_all[mu]; pos++) {
            int node = hyperedge_all[mu][pos];
            int idx = cur_all_idx[node]++;
            node_in_hyperedge[node][idx] = mu;
            node_in_hyperedge_pos[node][idx] = pos;
        }
    }
    free(cur_all_idx);
}

void assign_random_anchors() {
    hyperedge_anchors = (int**)malloc(Q * sizeof(int*));
    num_anchors = (int*)calloc(Q, sizeof(int));

    for (int mu = 0; mu < Q; mu++) {
        int max_members = size_T[mu] + size_H[mu];
        if (max_members == 0) {
            hyperedge_anchors[mu] = NULL;
            num_anchors[mu] = 0;
            continue;
        }

        int *members = (int*)malloc(max_members * sizeof(int));
        int members_count = 0;

        if (ANCHOR_SELECTION_MODE == ANCHOR_MODE_INPUT || ANCHOR_SELECTION_MODE == ANCHOR_MODE_GLOBAL) {
            memcpy(members + members_count, hyperedge_T[mu], size_T[mu] * sizeof(int));
            members_count += size_T[mu];
        }
        if (ANCHOR_SELECTION_MODE == ANCHOR_MODE_OUTPUT || ANCHOR_SELECTION_MODE == ANCHOR_MODE_GLOBAL) {
            memcpy(members + members_count, hyperedge_H[mu], size_H[mu] * sizeof(int));
            members_count += size_H[mu];
        }

        if (members_count == 0) {
            hyperedge_anchors[mu] = NULL;
            num_anchors[mu] = 0;
            free(members);
            continue;
        }

        int *unique_members = (int*)malloc(members_count * sizeof(int));
        int unique_count = 0;

        for (int i = 0; i < members_count; i++) {
            bool is_duplicate = false;
            for (int j = 0; j < unique_count; j++) {
                if (unique_members[j] == members[i]) { is_duplicate = true; break; }
            }
            if (!is_duplicate) unique_members[unique_count++] = members[i];
        }
        
        int *current_anchors = (int*)malloc(unique_count * sizeof(int));
        int anchor_count = 0;
        for (int i = 0; i < unique_count; i++) {
            if ((double)rand() / RAND_MAX <= THETA) {
                current_anchors[anchor_count++] = unique_members[i];
            }
        }

        num_anchors[mu] = anchor_count;
        if (anchor_count > 0) {
            hyperedge_anchors[mu] = (int*)malloc(anchor_count * sizeof(int));
            memcpy(hyperedge_anchors[mu], current_anchors, anchor_count * sizeof(int));
        } else {
            hyperedge_anchors[mu] = NULL;
        }

        free(members); 
        free(unique_members); 
        free(current_anchors);
    }
}

void save_hypergraph_structure(const char* filename) {
    FILE *outfile = fopen(filename, "w");
    if (!outfile) return;

    fprintf(outfile, "N=%d\nQ=%d\n\n", N, Q);
    for (int mu = 0; mu < Q; mu++) {
        fprintf(outfile, "H_T %d", mu);
        for (int t = 0; t < size_T[mu]; t++) fprintf(outfile, " %d", hyperedge_T[mu][t]);
        fprintf(outfile, "\nH_H %d", mu);
        for (int h = 0; h < size_H[mu]; h++) fprintf(outfile, " %d", hyperedge_H[mu][h]);
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
    for (int mu = 0; mu < Q; mu++) {
        if (num_anchors[mu] > 0) {
            fprintf(outfile, "ANCHOR %d", mu);
            for (int a = 0; a < num_anchors[mu]; a++) fprintf(outfile, " %d", hyperedge_anchors[mu][a]);
            fprintf(outfile, "\n");
        }
    }
    fclose(outfile);
}