#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#define USE_REAL_DATA 0
const char* HYPERGRAPH_FILE_PATH = "ijo1366_hypergraph.txt";
const char * output_filename = "percolation.txt";

#define MONTE_CARLO_RUNS 50
#define NETWORK_GENERATION_MODE 1

#define GAMMA_Q 3.5  
#define GAMMA_M 3.5 
#define MIN_DEGREE 1
#define MAX_DEGREE 1000 

#define N_SYNTHETIC 5000
#define Q_SYNTHETIC 1000
#define M_IN_SYNTHETIC 2
#define M_OUT_SYNTHETIC 2

#define MAX_ITERS 100
#define THETA 0.2  
#define CONVERGENCE_THRESHOLD 1e-9

int N; 
int Q;

int **hyperedge_T, **hyperedge_H;
int *size_T, *size_H;
int **node_as_input_of, **node_as_output_of;
int **node_as_input_of_pos, **node_as_output_of_pos;
int *k_in, *k_out;

int **hyperedge_anchors;
int *num_anchors;

int *dam_nodes;

int *vis;
int *size_cluster;
int *stack;
int stack_ptr;
int **scc_list;
int scc_count;
int *current_scc;
int current_scc_size;

void load_hypergraph_from_file(const char* filepath);
void generate_synthetic_hypergraph_fixed();
void generate_synthetic_hypergraph_uncorrelated();
void generate_synthetic_hypergraph_correlated();
void assign_random_anchors();
void build_reverse_mappings();
void save_hypergraph_structure(const char* filename);
void run_single_experiment(double damage_rate, FILE* outfile);
void shuffle_array(int *array, long n);
int* generate_power_law_sequence(int num_elements, double gamma, int min_degree, int max_degree, long* out_sum);
void rescale_sequence(int* seq, int num_elements, long current_sum, long target_sum);
void populate_hypergraph_from_stubs(long L_q_in, long L_q_out, int* q_in_seq, int* q_out_seq, int* m_in_seq, int* m_out_seq);
void kosaraju_dfs1(int u);
void kosaraju_dfs2(int u);
void find_all_sccs_simulation();
void DFS_GOC(int start_node);
void DFS_GIC(int start_node);
void run_message_passing_probabilistic(double p_N, double *pred_gout, double *pred_gin, double *pred_gscc);

int main() {
    srand(2025); 
    
    if (USE_REAL_DATA) {
        printf("--- Mode: Loading REAL hypergraph data ---\n");
        load_hypergraph_from_file(HYPERGRAPH_FILE_PATH);
    } else {
        printf("--- Mode: Generating SYNTHETIC hypergraph data ---\n");
        
        if (NETWORK_GENERATION_MODE == 0) {
            generate_synthetic_hypergraph_fixed();
        } else if (NETWORK_GENERATION_MODE == 1) {
            generate_synthetic_hypergraph_uncorrelated();
        } else if (NETWORK_GENERATION_MODE == 2) {
            generate_synthetic_hypergraph_correlated();
        } else {
            fprintf(stderr, "Error: Unknown NETWORK_GENERATION_MODE %d\n", NETWORK_GENERATION_MODE);
            return 1;
        }
    }

    assign_random_anchors();
    build_reverse_mappings();

    printf("Hypergraph setup complete. Starting simulation loop.\n");
    save_hypergraph_structure("result/hypergraph_structure_anchor.txt");
    
    dam_nodes = malloc(N * sizeof(int));
    vis = calloc(N, sizeof(int));
    size_cluster = calloc(N, sizeof(int));
    stack = malloc(N * sizeof(int));
    scc_list = calloc(N, sizeof(int*));
    current_scc = malloc(N * sizeof(int));
    
    FILE *outfile = fopen(output_filename, "w");
    if (outfile == NULL) {
        fprintf(stderr, "Error: Could not open file for writing.\n");
        return 1;
    }

    fprintf(outfile, "node_survival,mp_gout,mp_gin,mp_gscc,sim_gout,sim_gin,sim_gscc\n");

    for (int nc = 0; nc <= 100; nc++) {
        double p_damage = 1.0 - nc * 0.01;
        if(p_damage < 0.0) p_damage = 0.0;
        run_single_experiment(p_damage, outfile);
        fprintf(stderr, "\rFinished damage rate: %.4f (p_survival=%.4f)", p_damage, 1.0-p_damage);
        fflush(stderr);
    }

    fclose(outfile);
    fprintf(stderr, "\nSimulation finished. Data saved to %s.\n", output_filename);

    free(dam_nodes);
    for (int mu = 0; mu < Q; mu++) { 
        free(hyperedge_T[mu]); 
        free(hyperedge_H[mu]);
        if (hyperedge_anchors[mu]) free(hyperedge_anchors[mu]);
    }
    free(hyperedge_T); free(hyperedge_H);
    free(size_T); free(size_H);
    free(hyperedge_anchors);
    free(num_anchors);
    for (int i = 0; i < N; i++) {
        if (node_as_input_of[i]) free(node_as_input_of[i]);
        if (node_as_input_of_pos[i]) free(node_as_input_of_pos[i]);
        if (node_as_output_of[i]) free(node_as_output_of[i]);
        if (node_as_output_of_pos[i]) free(node_as_output_of_pos[i]);
    }
    free(node_as_input_of); free(node_as_input_of_pos);
    free(node_as_output_of); free(node_as_output_of_pos);
    free(k_in); free(k_out);
    free(vis);
    free(size_cluster);
    free(stack);
    for(int i=0; i<scc_count; i++) if(scc_list[i]) free(scc_list[i]);
    free(scc_list);
    free(current_scc);
    
    return 0;
}

void shuffle_array(int *array, long n) {
    if (n > 1) {
        for (long i = 0; i < n - 1; i++) {
            long j = i + rand() / (RAND_MAX / (n - i) + 1);
            int t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

int* generate_power_law_sequence(int num_elements, double gamma, int min_degree, int max_degree, long* out_sum) {
    int* seq = (int*)malloc(num_elements * sizeof(int));
    if (!seq) {
        fprintf(stderr, "Memory allocation failed for power-law sequence.\n");
        exit(EXIT_FAILURE);
    }
    
    *out_sum = 0;
    for (int i = 0; i < num_elements; i++) {
        double r = (double)rand() / (RAND_MAX + 1.0);
        double val = min_degree * pow(1.0 - r, -1.0 / (gamma - 1.0));
        int d = (int)round(val);
        
        if (d < min_degree) d = min_degree;
        if (d > max_degree) d = max_degree;
        
        seq[i] = d;
        *out_sum += d;
    }
    return seq;
}

void rescale_sequence(int* seq, int num_elements, long current_sum, long target_sum) {
    if (current_sum == 0) {
        fprintf(stderr, "Warning: sequence sum is 0, cannot rescale.\n");
        return;
    }

    long new_sum = 0;
    double scale = (double)target_sum / (double)current_sum;
    
    for (int i = 0; i < num_elements; i++) {
        seq[i] = (int)round(seq[i] * scale);

        if (seq[i] == 0 && (int)round(seq[i] / scale) > 0) seq[i] = 1;
        new_sum += seq[i];
    }
    
    long diff = target_sum - new_sum;
    for (long i = 0; i < labs(diff); i++) {
        int node = rand() % num_elements;
        if (diff > 0) {
            seq[node]++;
        } else if (seq[node] > MIN_DEGREE) { 
            seq[node]--;
        }
    }
    
    long final_sum = 0;
    for(int i=0; i<num_elements; i++) final_sum += seq[i];
    if(final_sum != target_sum) {
         fprintf(stderr, "Warning: Rescale final sum %ld != target %ld. Retrying...\n", final_sum, target_sum);

         if (final_sum < target_sum) seq[rand() % num_elements] += (target_sum - final_sum);
         if (final_sum > target_sum) seq[rand() % num_elements] -= (final_sum - target_sum);
    }
}

void populate_hypergraph_from_stubs(long L_q_in, long L_q_out, int* q_in_seq, int* q_out_seq, int* m_in_seq, int* m_out_seq) {
    int *q_in_stubs = (int*)malloc(L_q_in * sizeof(int));   
    int *m_out_stubs = (int*)malloc(L_q_in * sizeof(int)); 
    
    int *q_out_stubs = (int*)malloc(L_q_out * sizeof(int)); 
    int *m_in_stubs = (int*)malloc(L_q_out * sizeof(int));
    
    long idx = 0;
    for(int i=0; i<N; i++) for(int j=0; j<q_in_seq[i]; j++) q_in_stubs[idx++] = i;
    idx = 0;
    for(int i=0; i<Q; i++) for(int j=0; j<m_out_seq[i]; j++) m_out_stubs[idx++] = i;
    
    idx = 0;
    for(int i=0; i<N; i++) for(int j=0; j<q_out_seq[i]; j++) q_out_stubs[idx++] = i;
    idx = 0;
    for(int i=0; i<Q; i++) for(int j=0; j<m_in_seq[i]; j++) m_in_stubs[idx++] = i;

    shuffle_array(q_in_stubs, L_q_in);
    shuffle_array(m_out_stubs, L_q_in);
    shuffle_array(q_out_stubs, L_q_out);
    shuffle_array(m_in_stubs, L_q_out);
    
    hyperedge_T = (int**)malloc(Q * sizeof(int*));
    hyperedge_H = (int**)malloc(Q * sizeof(int*));
    size_T = (int*)calloc(Q, sizeof(int)); 
    size_H = (int*)calloc(Q, sizeof(int));
    
    for(int mu=0; mu<Q; mu++) {
        hyperedge_T[mu] = (int*)malloc(m_in_seq[mu] * sizeof(int));
        size_T[mu] = m_in_seq[mu]; 
        hyperedge_H[mu] = (int*)malloc(m_out_seq[mu] * sizeof(int));
        size_H[mu] = m_out_seq[mu]; 
    }
    
    int *current_T_idx = (int*)calloc(Q, sizeof(int));
    int *current_H_idx = (int*)calloc(Q, sizeof(int));

    for(long i=0; i<L_q_out; i++) {
        int node_i = q_out_stubs[i];      
        int hyperedge_mu = m_in_stubs[i]; 
        if (current_T_idx[hyperedge_mu] < size_T[hyperedge_mu]) {
            hyperedge_T[hyperedge_mu][current_T_idx[hyperedge_mu]++] = node_i;
        }
    }
    
    for(long i=0; i<L_q_in; i++) {
        int node_i = q_in_stubs[i];       
        int hyperedge_mu = m_out_stubs[i]; 
        if (current_H_idx[hyperedge_mu] < size_H[hyperedge_mu]) {
            hyperedge_H[hyperedge_mu][current_H_idx[hyperedge_mu]++] = node_i;
        }
    }
    
    free(q_in_stubs); free(m_in_stubs); free(q_out_stubs); free(m_out_stubs);
    free(current_T_idx); free(current_H_idx);
}

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

    if (!hyperedge_T || !hyperedge_H || !size_T || !size_H) {
        fprintf(stderr, "Memory allocation failed for hypergraph structure.\n");
        exit(EXIT_FAILURE);
    }

    for (int mu = 0; mu < Q; mu++) {
        if (fscanf(file, "%d %d\n", &size_T[mu], &size_H[mu]) != 2) {
            fprintf(stderr, "Error reading size for hyperedge %d\n", mu);
            exit(EXIT_FAILURE);
        }
        
        hyperedge_T[mu] = (int*)calloc(size_T[mu], sizeof(int));
        hyperedge_H[mu] = (int*)calloc(size_H[mu], sizeof(int));
        if (!hyperedge_T[mu] || !hyperedge_H[mu]) {
            fprintf(stderr, "Memory allocation failed for hyperedge %d details.\n", mu);
            exit(EXIT_FAILURE);
        }

        for (int n = 0; n < size_T[mu]; n++) {
            if (fscanf(file, "%d", &hyperedge_T[mu][n]) != 1) {
                fprintf(stderr, "Error reading tail node for hyperedge %d\n", mu);
                exit(EXIT_FAILURE);
            }
        }
        for (int n = 0; n < size_H[mu]; n++) {
            if (fscanf(file, "%d", &hyperedge_H[mu][n]) != 1) {
                fprintf(stderr, "Error reading head node for hyperedge %d\n", mu);
                exit(EXIT_FAILURE);
            }
        }
    }
    fclose(file);
    printf("Hypergraph data loaded successfully from file.\n");
}

void assign_random_anchors() {
    hyperedge_anchors = malloc(Q * sizeof(int*));
    num_anchors = calloc(Q, sizeof(int));

    for (int mu = 0; mu < Q; mu++) {
        int max_members = size_T[mu] + size_H[mu];
        if (max_members == 0) {
            hyperedge_anchors[mu] = NULL;
            num_anchors[mu] = 0;
            continue;
        }
        
        int *members = malloc(max_members * sizeof(int));
        int *unique_members = malloc(max_members * sizeof(int));
        int unique_count = 0;

        memcpy(members, hyperedge_T[mu], size_T[mu] * sizeof(int));
        memcpy(members + size_T[mu], hyperedge_H[mu], size_H[mu] * sizeof(int));

        for (int i = 0; i < max_members; i++) {
            bool is_duplicate = false;
            for (int j = 0; j < unique_count; j++) {
                if (unique_members[j] == members[i]) {
                    is_duplicate = true;
                    break;
                }
            }
            if (!is_duplicate) {
                unique_members[unique_count++] = members[i];
            }
        }
        
        int current_anchors[unique_count];
        int anchor_count = 0;
        for (int i = 0; i < unique_count; i++) {
            if ((double)rand() / RAND_MAX < THETA) {
                current_anchors[anchor_count++] = unique_members[i];
            }
        }

        num_anchors[mu] = anchor_count;
        if (anchor_count > 0) {
            hyperedge_anchors[mu] = malloc(anchor_count * sizeof(int));
            memcpy(hyperedge_anchors[mu], current_anchors, anchor_count * sizeof(int));
        } else {
            hyperedge_anchors[mu] = NULL;
        }

        free(members);
        free(unique_members);
    }
    printf("Random anchors (THETA=%.2f) assigned to all hyperedges.\n", THETA);
}

void generate_synthetic_hypergraph_fixed() {
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
        
        for (int t = 0; t < size_T[mu]; t++) {
            hyperedge_T[mu][t] = rand() % N;
        }
        for (int h = 0; h < size_H[mu]; h++) {
            hyperedge_H[mu][h] = rand() % N;
        }
    }
    printf("Synthetic (Fixed Cardinality) hypergraph generated: N=%d, Q=%d, M_in=%d, M_out=%d.\n", N, Q, M_IN_SYNTHETIC, M_OUT_SYNTHETIC);
}

void generate_synthetic_hypergraph_uncorrelated() {
    N = N_SYNTHETIC;
    Q = Q_SYNTHETIC;
    
    long sum_qin, sum_qout, sum_min, sum_mout;
    
    int* q_in_seq = generate_power_law_sequence(N, GAMMA_Q, MIN_DEGREE, MAX_DEGREE, &sum_qin);
    int* q_out_seq = generate_power_law_sequence(N, GAMMA_Q, MIN_DEGREE, MAX_DEGREE, &sum_qout);
    int* m_in_seq = generate_power_law_sequence(Q, GAMMA_M, MIN_DEGREE, MAX_DEGREE, &sum_min);
    int* m_out_seq = generate_power_law_sequence(Q, GAMMA_M, MIN_DEGREE, MAX_DEGREE, &sum_mout);

    rescale_sequence(m_out_seq, Q, sum_mout, sum_qin); 
    rescale_sequence(m_in_seq, Q, sum_min, sum_qout); 

    long L_q_in = 0; for(int i=0; i<N; i++) L_q_in += q_in_seq[i]; 
    long L_q_out = 0; for(int i=0; i<N; i++) L_q_out += q_out_seq[i]; 

    printf("Generating Uncorrelated PL Graph: L_q_in/m_out = %ld, L_q_out/m_in = %ld\n", L_q_in, L_q_out);

    populate_hypergraph_from_stubs(L_q_in, L_q_out, q_in_seq, q_out_seq, m_in_seq, m_out_seq);
    
    printf("Synthetic UNCORRELATED power-law hypergraph generated (g_q=%.2f, g_m=%.2f).\n", GAMMA_Q, GAMMA_M);
    
    free(q_in_seq); free(q_out_seq); free(m_in_seq); free(m_out_seq);
}

void generate_synthetic_hypergraph_correlated() {
    N = N_SYNTHETIC;
    Q = Q_SYNTHETIC;
    
    long sum_q, sum_m;
    
    int* q_seq = generate_power_law_sequence(N, GAMMA_Q, MIN_DEGREE, MAX_DEGREE, &sum_q);
    int* m_seq = generate_power_law_sequence(Q, GAMMA_M, MIN_DEGREE, MAX_DEGREE, &sum_m);
    
    rescale_sequence(m_seq, Q, sum_m, sum_q);
    
    long L = 0; for(int i=0; i<N; i++) L += q_seq[i];
    printf("Generating Correlated PL Graph: L_total = %ld\n", L);

    populate_hypergraph_from_stubs(L, L, q_seq, q_seq, m_seq, m_seq);
    
    printf("Synthetic MAXIMALLY CORRELATED power-law hypergraph generated (g_q=%.2f, g_m=%.2f).\n", GAMMA_Q, GAMMA_M);
    
    free(q_seq); free(m_seq);
}

void build_reverse_mappings() {
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
    free(cur_in_idx); 
    free(cur_out_idx);
    printf("Reverse mappings (k_in, k_out, etc.) built.\n");
}

void save_hypergraph_structure(const char* filename) {
    FILE *outfile = fopen(filename, "w");
    if (outfile == NULL) {
        fprintf(stderr, "错误：无法打开文件 %s 进行写入。\n", filename);
        return;
    }

    fprintf(stderr, "正在将网络结构保存到 %s...\n", filename);

    fprintf(outfile, "N=%d\n", N);
    fprintf(outfile, "Q=%d\n", Q);
    fprintf(outfile, "\n");

    for (int mu = 0; mu < Q; mu++) {
        fprintf(outfile, "H_T %d", mu);
        for (int t = 0; t < size_T[mu]; t++) {
            fprintf(outfile, " %d", hyperedge_T[mu][t]);
        }
        fprintf(outfile, "\n");

        fprintf(outfile, "H_H %d", mu);
        for (int h = 0; h < size_H[mu]; h++) {
            fprintf(outfile, " %d", hyperedge_H[mu][h]);
        }
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");

    for (int mu = 0; mu < Q; mu++) {
        if (num_anchors[mu] > 0) {
            fprintf(outfile, "ANCHOR %d", mu);
            for (int a = 0; a < num_anchors[mu]; a++) {
                fprintf(outfile, " %d", hyperedge_anchors[mu][a]);
            }
            fprintf(outfile, "\n");
        }
    }

    fclose(outfile);
    fprintf(stderr, "网络结构保存完毕。\n");
}

void kosaraju_dfs1(int u) {
    if (vis[u] || dam_nodes[u] == 0) return;
    vis[u] = 1;
    for (int i = 0; i < k_out[u]; i++) {
        int hyperedge_idx = node_as_input_of[u][i];
        int is_active = 1;
        for (int j = 0; j < num_anchors[hyperedge_idx]; j++) {
            if (dam_nodes[hyperedge_anchors[hyperedge_idx][j]] == 0) {
                is_active = 0;
                break;
            }
        }
        if (is_active) {
            for (int j = 0; j < size_H[hyperedge_idx]; j++) {
                if(dam_nodes[hyperedge_H[hyperedge_idx][j]]==0) continue;
                kosaraju_dfs1(hyperedge_H[hyperedge_idx][j]);
            }
        }
    }
    stack[stack_ptr++] = u;
}

void kosaraju_dfs2(int u) {
    if (vis[u] || dam_nodes[u] == 0) return;
    vis[u] = 1;
    current_scc[current_scc_size++] = u;
    for (int i = 0; i < k_in[u]; i++) {
        int hyperedge_idx = node_as_output_of[u][i];
        int is_active = 1;
        for (int j = 0; j < num_anchors[hyperedge_idx]; j++) {
            if (dam_nodes[hyperedge_anchors[hyperedge_idx][j]] == 0) {
                is_active = 0;
                break;
            }
        }
        if (is_active) {
            for (int j = 0; j < size_T[hyperedge_idx]; j++) {
                if(dam_nodes[hyperedge_T[hyperedge_idx][j]]==0) continue;
                kosaraju_dfs2(hyperedge_T[hyperedge_idx][j]);
            }
        }
    }
}

void find_all_sccs_simulation() {
    stack_ptr = 0;
    memset(vis, 0, N * sizeof(int));
    for (int i = 0; i < N; i++) {
        if (dam_nodes[i] == 1 && !vis[i]) kosaraju_dfs1(i);
    }
    scc_count = 0;
    memset(vis, 0, N * sizeof(int));
    while (stack_ptr > 0) {
        int u = stack[--stack_ptr];
        if (dam_nodes[u] == 1 && !vis[u]) {
            current_scc_size = 0;
            kosaraju_dfs2(u);
            if (current_scc_size > 0) {
                scc_list[scc_count] = (int*)malloc(current_scc_size * sizeof(int));
                memcpy(scc_list[scc_count], current_scc, current_scc_size * sizeof(int));
                size_cluster[scc_count] = current_scc_size;
                scc_count++;
            }
        }
    }
}

void DFS_GOC(int start_node) {
    if (vis[start_node] || dam_nodes[start_node] == 0) return;
    vis[start_node] = 1;
    for (int i = 0; i < k_out[start_node]; i++) {
        int hyperedge_idx = node_as_input_of[start_node][i];
        int is_active = 1;
        for (int j = 0; j < num_anchors[hyperedge_idx]; j++) {
            if (dam_nodes[hyperedge_anchors[hyperedge_idx][j]] == 0) {
                is_active = 0;
                break;
            }
        }
        if (is_active) {
            for (int j = 0; j < size_H[hyperedge_idx]; j++) {
                if(dam_nodes[hyperedge_H[hyperedge_idx][j]]==0) continue;
                DFS_GOC(hyperedge_H[hyperedge_idx][j]);
            }
        }
    }
}

void DFS_GIC(int start_node) {
    if (vis[start_node] || dam_nodes[start_node] == 0) return;
    vis[start_node] = 1;
    for (int i = 0; i < k_in[start_node]; i++) {
        int hyperedge_idx = node_as_output_of[start_node][i];
        int is_active = 1;
        for (int j = 0; j < num_anchors[hyperedge_idx]; j++) {
            if (dam_nodes[hyperedge_anchors[hyperedge_idx][j]] == 0) {
                is_active = 0;
                break;
            }
        }
        if (is_active) {
            for (int j = 0; j < size_T[hyperedge_idx]; j++) {
                if(dam_nodes[hyperedge_T[hyperedge_idx][j]]==0) continue;
                DFS_GIC(hyperedge_T[hyperedge_idx][j]);
            }
        }
    }
}

void run_message_passing_probabilistic(double p_N, double *pred_gout, double *pred_gin, double *pred_gscc) {
    double p_H = 1.0;

    double **w_gout = (double**)malloc(N * sizeof(double*));
    double **v_gout = (double**)malloc(Q * sizeof(double*));
    double **w_gin = (double**)malloc(N * sizeof(double*));
    double **v_gin = (double**)malloc(Q * sizeof(double*));
    double **w_gout_new = (double**)malloc(N * sizeof(double*));
    double **v_gout_new = (double**)malloc(Q * sizeof(double*));
    double **w_gin_new = (double**)malloc(N * sizeof(double*));
    double **v_gin_new = (double**)malloc(Q * sizeof(double*));

    for(int i=0; i<N; i++) {
        w_gout[i] = calloc(k_out[i], sizeof(double));
        w_gout_new[i] = calloc(k_out[i], sizeof(double));
        w_gin[i] = calloc(k_in[i], sizeof(double));
        w_gin_new[i] = calloc(k_in[i], sizeof(double));
    }
    for(int mu=0; mu<Q; mu++) {
        v_gout[mu] = calloc(size_H[mu], sizeof(double));
        v_gout_new[mu] = calloc(size_H[mu], sizeof(double));
        v_gin[mu] = calloc(size_T[mu], sizeof(double));
        v_gin_new[mu] = calloc(size_T[mu], sizeof(double));
    }

    for(int i=0; i<N; i++) {
        for(int k=0; k<k_out[i]; k++) w_gout[i][k] = p_N; 
        for(int k=0; k<k_in[i]; k++) w_gin[i][k] = p_N;  
    }
    for(int mu=0; mu<Q; mu++) {
        for(int h=0; h<size_H[mu]; h++) v_gout[mu][h] = 1.0;
        for(int t=0; t<size_T[mu]; t++) v_gin[mu][t] = 1.0;
    }
    
    int iter = 0;
    while(iter < MAX_ITERS) {
        double max_diff = 0.0;

        for (int mu = 0; mu < Q; mu++) {
            for (int h = 0; h < size_H[mu]; h++) {
                int i = hyperedge_H[mu][h];
                
                int Z_ia = 0; 
                for(int j = 0; j < num_anchors[mu]; j++) {
                    if (hyperedge_anchors[mu][j] != i) Z_ia++;
                }
                double prob_anchors_ok = pow(p_N, Z_ia);

                double prod_term = 1.0;
                for (int t = 0; t < size_T[mu]; t++) {
                    int j = hyperedge_T[mu][t];
                    int k_idx = -1;
                    for(int k=0;k<k_out[j];++k) if(node_as_input_of[j][k]==mu) { k_idx=k; break;}
                    if(k_idx != -1) prod_term *= (1.0 - w_gout[j][k_idx]);
                }
                v_gout_new[mu][h] = p_H * prob_anchors_ok * (1.0 - prod_term);
            }
        }
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < k_out[i]; k++) {
                int mu = node_as_input_of[i][k];
                
                bool is_anchor = false;
                for(int j=0; j<num_anchors[mu]; j++) if(hyperedge_anchors[mu][j] == i) {is_anchor = true; break;}
                double prob_node_modifier = pow(p_N, 1 - (is_anchor ? 1 : 0));
                
                double prod_term = 1.0;
                for (int h_idx = 0; h_idx < k_in[i]; h_idx++) {
                    int beta = node_as_output_of[i][h_idx];
                    if (beta == mu) continue; 
                    int beta_pos = node_as_output_of_pos[i][h_idx];
                    prod_term *= (1.0 - v_gout[beta][beta_pos]);
                }
                w_gout_new[i][k] = prob_node_modifier * (1.0 - prod_term);
            }
        }

        for (int mu = 0; mu < Q; mu++) {
            for (int t = 0; t < size_T[mu]; t++) {
                int i = hyperedge_T[mu][t];
                
                int Z_ia = 0;
                for(int j = 0; j < num_anchors[mu]; j++) {
                    if (hyperedge_anchors[mu][j] != i) Z_ia++;
                }
                double prob_anchors_ok = pow(p_N, Z_ia);

                double prod_term = 1.0;
                for (int h = 0; h < size_H[mu]; h++) {
                    int j = hyperedge_H[mu][h];
                    int k_idx = -1;
                    for(int k=0; k<k_in[j]; ++k) if(node_as_output_of[j][k]==mu) { k_idx=k; break;}
                    if(k_idx != -1) prod_term *= (1.0 - w_gin[j][k_idx]);
                }
                v_gin_new[mu][t] = p_H * prob_anchors_ok * (1.0 - prod_term);
            }
        }
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < k_in[i]; k++) {
                int mu = node_as_output_of[i][k];

                bool is_anchor = false;
                for(int j=0; j<num_anchors[mu]; j++) if(hyperedge_anchors[mu][j] == i) {is_anchor = true; break;}
                double prob_node_modifier = pow(p_N, 1 - (is_anchor ? 1 : 0));
                
                double prod_term = 1.0;
                for (int t_idx = 0; t_idx < k_out[i]; t_idx++) {
                    int beta = node_as_input_of[i][t_idx];
                    if (beta == mu) continue; 
                    int beta_pos = node_as_input_of_pos[i][t_idx];
                    prod_term *= (1.0 - v_gin[beta][beta_pos]);
                }
                w_gin_new[i][k] = prob_node_modifier * (1.0 - prod_term);
            }
        }
        
        for(int i=0; i<N; i++) for(int k=0; k<k_out[i]; k++) { max_diff=fmax(max_diff, fabs(w_gout[i][k]-w_gout_new[i][k])); w_gout[i][k]=w_gout_new[i][k]; }
        for(int mu=0; mu<Q; mu++) for(int h=0; h<size_H[mu]; h++) { max_diff=fmax(max_diff, fabs(v_gout[mu][h]-v_gout_new[mu][h])); v_gout[mu][h]=v_gout_new[mu][h]; }
        for(int i=0; i<N; i++) for(int k=0; k<k_in[i]; k++) { max_diff=fmax(max_diff, fabs(w_gin[i][k]-w_gin_new[i][k])); w_gin[i][k]=w_gin_new[i][k]; }
        for(int mu=0; mu<Q; mu++) for(int t=0; t<size_T[mu]; t++) { max_diff=fmax(max_diff, fabs(v_gin[mu][t]-v_gin_new[mu][t])); v_gin[mu][t]=v_gin_new[mu][t]; }

        if (max_diff < CONVERGENCE_THRESHOLD) break;
        iter++;
    }

    double gout_total = 0, gin_total = 0, gscc_total = 0;
    for (int i = 0; i < N; i++) {
        double prod_gout = 1.0;
        for (int k = 0; k < k_in[i]; k++) {
            int mu = node_as_output_of[i][k];
            int h = node_as_output_of_pos[i][k];
            prod_gout *= (1.0 - v_gout[mu][h]);
        }
        double prob_in_gout = p_N * (1.0 - prod_gout); 
        
        double prod_gin = 1.0;
        for (int k = 0; k < k_out[i]; k++) {
            int mu = node_as_input_of[i][k];
            int t = node_as_input_of_pos[i][k];
            prod_gin *= (1.0 - v_gin[mu][t]);
        }
        double prob_in_gin = p_N * (1.0 - prod_gin); 
        
        gout_total += prob_in_gout;
        gin_total += prob_in_gin;
        gscc_total += p_N * (1.0 - prod_gout) * (1.0 - prod_gin); 
    }
    *pred_gout = gout_total;
    *pred_gin = gin_total;
    *pred_gscc = gscc_total;

    for(int i=0; i<N; i++) { free(w_gout[i]); free(w_gout_new[i]); free(w_gin[i]); free(w_gin_new[i]); }
    for(int mu=0; mu<Q; mu++) { free(v_gout[mu]); free(v_gout_new[mu]); free(v_gin[mu]); free(v_gin_new[mu]); }
    free(w_gout); free(w_gout_new); free(w_gin); free(w_gin_new);
    free(v_gout); free(v_gout_new); free(v_gin); free(v_gin_new);
}

void run_single_experiment(double damage_rate, FILE* outfile) {
    double p_N = 1.0 - damage_rate;
    double gout_size_mp, gin_size_mp, gscc_size_mp;
    run_message_passing_probabilistic(p_N, &gout_size_mp, &gin_size_mp, &gscc_size_mp);

    double total_gout_sim = 0, total_gin_sim = 0, total_gscc_sim = 0;

    for (int run = 0; run < MONTE_CARLO_RUNS; ++run) {
        for (int i = 0; i < N; i++) {
            dam_nodes[i] = ((double)rand() / RAND_MAX < damage_rate) ? 0 : 1;
        }

        for(int i = 0; i < scc_count; i++) {
            if(scc_list[i] != NULL) { free(scc_list[i]); scc_list[i] = NULL; }
        }
        memset(size_cluster, 0, N * sizeof(int));
    
        find_all_sccs_simulation();
        int gscc_size_sim = 0, gscc_idx = -1;
        for (int i = 0; i < scc_count; i++) {
            if (size_cluster[i] > gscc_size_sim) {
                gscc_size_sim = size_cluster[i];
                gscc_idx = i;
            }
        }

        int gout_size_sim = 0, gin_size_sim = 0;
        if (gscc_idx != -1) {
            memset(vis, 0, N * sizeof(int));
            for (int j = 0; j < gscc_size_sim; j++) DFS_GOC(scc_list[gscc_idx][j]);
            for (int i = 0; i < N; i++) if (vis[i]) gout_size_sim++;

            memset(vis, 0, N * sizeof(int));
            for (int j = 0; j < gscc_size_sim; j++) DFS_GIC(scc_list[gscc_idx][j]);
            for (int i = 0; i < N; i++) if (vis[i]) gin_size_sim++;
        }
        
        total_gout_sim += gout_size_sim;
        total_gin_sim += gin_size_sim;
        total_gscc_sim += gscc_size_sim;
    }

    double avg_gout_sim = total_gout_sim / MONTE_CARLO_RUNS;
    double avg_gin_sim = total_gin_sim / MONTE_CARLO_RUNS;
    double avg_gscc_sim = total_gscc_sim / MONTE_CARLO_RUNS;

    fprintf(outfile, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
           1 - damage_rate,
           gout_size_mp / N, gin_size_mp / N, gscc_size_mp / N,
           avg_gout_sim / N, avg_gin_sim / N, avg_gscc_sim / N);
}