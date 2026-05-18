#include "globals.h"

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

int find_max_hwcc_simulation() {
    memset(vis_node_hwcc, 0, N * sizeof(int));
    memset(vis_edge_hwcc, 0, Q * sizeof(int));
    int max_wcc_size = 0;

    for (int i = 0; i < N; i++) {
        if (dam_nodes[i] == 1 && !vis_node_hwcc[i]) {
            int head = 0, tail = 0;
            q_nodes_hwcc[tail++] = i;
            vis_node_hwcc[i] = 1;
            int current_wcc_size = 0;

            while (head < tail) {
                int u = q_nodes_hwcc[head++];
                current_wcc_size++;

                for (int k = 0; k < k_all[u]; k++) {
                    int mu = node_in_hyperedge[u][k];
                    if (vis_edge_hwcc[mu]) continue;

                    int is_active = dam_edges[mu];
                    if (is_active) {
                        for (int j = 0; j < num_anchors[mu]; j++) {
                            if (dam_nodes[hyperedge_anchors[mu][j]] == 0) {
                                is_active = 0; break;
                            }
                        }
                    }

                    if (is_active) {
                        vis_edge_hwcc[mu] = 1;
                        for (int pos = 0; pos < size_all[mu]; pos++) {
                            int v = hyperedge_all[mu][pos];
                            if (dam_nodes[v] == 1 && !vis_node_hwcc[v]) {
                                vis_node_hwcc[v] = 1;
                                q_nodes_hwcc[tail++] = v;
                            }
                        }
                    }
                }
            }
            if (current_wcc_size > max_wcc_size) max_wcc_size = current_wcc_size;
        }
    }
    return max_wcc_size;
}

void run_message_passing_hwcc(double p_N, double p_H, double *pred_hwcc) {
    double **w_hwcc = malloc(N * sizeof(double*));
    double **v_hwcc = malloc(Q * sizeof(double*));
    double **w_hwcc_new = malloc(N * sizeof(double*));
    double **v_hwcc_new = malloc(Q * sizeof(double*));

    for(int i=0; i<N; i++) {
        w_hwcc[i] = calloc(k_all[i], sizeof(double));
        w_hwcc_new[i] = calloc(k_all[i], sizeof(double));
        for(int k=0; k<k_all[i]; k++) w_hwcc[i][k] = p_N;
    }
    for(int mu=0; mu<Q; mu++) {
        v_hwcc[mu] = calloc(size_all[mu], sizeof(double));
        v_hwcc_new[mu] = calloc(size_all[mu], sizeof(double));
        for(int pos=0; pos<size_all[mu]; pos++) v_hwcc[mu][pos] = p_H;
    }

    int iter = 0;
    while(iter < MAX_ITERS) {
        double max_diff = 0.0;

        for (int mu = 0; mu < Q; mu++) {
            for (int pos = 0; pos < size_all[mu]; pos++) {
                int i = hyperedge_all[mu][pos];
                int Z_ia = 0;
                for(int j = 0; j < num_anchors[mu]; j++) if (hyperedge_anchors[mu][j] != i) Z_ia++;
                double prob_anchors_ok = pow(p_N, Z_ia);

                double prod_term = 1.0;
                for (int t = 0; t < size_all[mu]; t++) {
                    if (t == pos) continue; 
                    int j = hyperedge_all[mu][t];
                    int k_idx = -1;
                    for(int k=0; k<k_all[j]; ++k) if(node_in_hyperedge[j][k] == mu) { k_idx=k; break;}
                    if(k_idx != -1) prod_term *= (1.0 - w_hwcc[j][k_idx]);
                }
                v_hwcc_new[mu][pos] = p_H * prob_anchors_ok * (1.0 - prod_term);
            }
        }

        for (int i = 0; i < N; i++) {
            for (int k = 0; k < k_all[i]; k++) {
                int mu = node_in_hyperedge[i][k];
                bool is_anchor = false;
                for(int j=0; j<num_anchors[mu]; j++) if(hyperedge_anchors[mu][j] == i) {is_anchor = true; break;}
                double prob_node_modifier = pow(p_N, 1 - (is_anchor ? 1 : 0));

                double prod_term = 1.0;
                for (int h_idx = 0; h_idx < k_all[i]; h_idx++) {
                    int beta = node_in_hyperedge[i][h_idx];
                    if (beta == mu) continue;
                    int beta_pos = node_in_hyperedge_pos[i][h_idx];
                    prod_term *= (1.0 - v_hwcc[beta][beta_pos]);
                }
                w_hwcc_new[i][k] = prob_node_modifier * (1.0 - prod_term);
            }
        }

        for(int i=0; i<N; i++) {
            for(int k=0; k<k_all[i]; k++) { 
                max_diff=fmax(max_diff, fabs(w_hwcc[i][k]-w_hwcc_new[i][k])); 
                w_hwcc[i][k]=w_hwcc_new[i][k]; 
            }
        }
        for(int mu=0; mu<Q; mu++) {
            for(int pos=0; pos<size_all[mu]; pos++) { 
                max_diff=fmax(max_diff, fabs(v_hwcc[mu][pos]-v_hwcc_new[mu][pos])); 
                v_hwcc[mu][pos]=v_hwcc_new[mu][pos]; 
            }
        }
        if (max_diff < CONVERGENCE_THRESHOLD) break;
        iter++;
    }

    double hwcc_total = 0;
    for (int i = 0; i < N; i++) {
        double prod_hwcc = 1.0;
        for (int k = 0; k < k_all[i]; k++) {
            int mu = node_in_hyperedge[i][k];
            int pos = node_in_hyperedge_pos[i][k];
            prod_hwcc *= (1.0 - v_hwcc[mu][pos]);
        }
        hwcc_total += p_N * (1.0 - prod_hwcc);
    }
    *pred_hwcc = hwcc_total;

    for(int i=0; i<N; i++) { free(w_hwcc[i]); free(w_hwcc_new[i]); }
    for(int mu=0; mu<Q; mu++) { free(v_hwcc[mu]); free(v_hwcc_new[mu]); }
    free(w_hwcc); free(w_hwcc_new); free(v_hwcc); free(v_hwcc_new);
}


void run_message_passing_probabilistic(double p_N, double p_H, double *pred_gout, double *pred_gin, double *pred_gscc) {
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
        for(int h=0; h<size_H[mu]; h++) v_gout[mu][h] = p_H;
        for(int t=0; t<size_T[mu]; t++) v_gin[mu][t] = p_H;
    }
    int iter = 0;
    while(iter < MAX_ITERS) {
        double max_diff = 0.0;

        // --- GOUT (+) MESSAGE UPDATE ---
        // v_alpha -> i (hyperedge to node)
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
        // w_i -> alpha (node to hyperedge)
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < k_out[i]; k++) {
                int mu = node_as_input_of[i][k];
                
                bool is_anchor = false;
                for(int j=0; j<num_anchors[mu]; j++) if(hyperedge_anchors[mu][j] == i) {is_anchor = true; break;}
                double prob_node_modifier = pow(p_N, 1 - (is_anchor ? 1 : 0));
                
                double prod_term = 1.0;
                for (int h_idx = 0; h_idx < k_in[i]; h_idx++) {
                    int beta = node_as_output_of[i][h_idx];
                    if (beta == mu) continue; // Cavity
                    int beta_pos = node_as_output_of_pos[i][h_idx];
                    prod_term *= (1.0 - v_gout[beta][beta_pos]);
                }
                w_gout_new[i][k] = prob_node_modifier * (1.0 - prod_term);
            }
        }

        // --- GIN (-) MESSAGE UPDATE ---
        // v_alpha -> i (backwards)
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
        // w_i -> alpha (backwards)
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
        
        // --- Convergence Check and Swap ---
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

void run_single_experiment(double p_N, double p_H, FILE* outfile) {
    double hwcc_size_mp, gout_size_mp, gin_size_mp, gscc_size_mp;

    run_message_passing_probabilistic(p_N, p_H, &gout_size_mp, &gin_size_mp, &gscc_size_mp);

    run_message_passing_hwcc(p_N, p_H, &hwcc_size_mp);

    double total_hwcc_sim = 0, total_gout_sim = 0, total_gin_sim = 0, total_gscc_sim = 0;

    for (int run = 0; run < MONTE_CARLO_RUNS; ++run) {
        
        for (int i = 0; i < N; i++) dam_nodes[i] = ((double)rand() / RAND_MAX < p_N) ? 1 : 0;
        for (int mu = 0; mu < Q; mu++) dam_edges[mu] = ((double)rand() / RAND_MAX < p_H) ? 1 : 0;

        int hwcc_size_sim = find_max_hwcc_simulation();
        total_hwcc_sim += hwcc_size_sim;
        
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

    double avg_hwcc_sim = total_hwcc_sim / MONTE_CARLO_RUNS;
    double avg_gout_sim = total_gout_sim / MONTE_CARLO_RUNS;
    double avg_gin_sim = total_gin_sim / MONTE_CARLO_RUNS;
    double avg_gscc_sim = total_gscc_sim / MONTE_CARLO_RUNS;


    fprintf(outfile, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
           p_N, p_H,
           hwcc_size_mp / N, gout_size_mp / N, gin_size_mp / N, gscc_size_mp / N,
           avg_hwcc_sim / N, avg_gout_sim / N, avg_gin_sim / N, avg_gscc_sim / N);
}