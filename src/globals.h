#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <unistd.h>

#define USE_REAL_DATA 0

#define ANCHOR_MODE_INPUT  1  
#define ANCHOR_MODE_OUTPUT 2  
#define ANCHOR_MODE_GLOBAL 3  


#define ANCHOR_SELECTION_MODE 1
#define MONTE_CARLO_RUNS 50

#define N_SYNTHETIC 10000
#define Q_SYNTHETIC 10000
#define M_IN_SYNTHETIC 3
#define M_OUT_SYNTHETIC 2

#define MAX_ITERS 100
#define THETA 0.2
#define CONVERGENCE_THRESHOLD 1e-9

#define SWEEP_MODE 1
#define FIXED_P_H 1.0 
#define FIXED_P_N 1.0 

extern const char* HYPERGRAPH_FILE_PATH;
extern const char* output_filename;


extern int N; 
extern int Q; 

extern int **hyperedge_T, **hyperedge_H;
extern int *size_T, *size_H;
extern int **node_as_input_of, **node_as_output_of;
extern int **node_as_input_of_pos, **node_as_output_of_pos;
extern int *k_in, *k_out;

extern int **hyperedge_all;
extern int *size_all;
extern int **node_in_hyperedge;
extern int **node_in_hyperedge_pos;
extern int *k_all;

extern int **hyperedge_anchors;
extern int *num_anchors;

extern int *dam_nodes;
extern int *dam_edges; 
extern int *vis;
extern int *size_cluster;
extern int *stack;
extern int stack_ptr;
extern int **scc_list;
extern int scc_count;
extern int *current_scc;
extern int current_scc_size;

extern int *vis_node_hwcc;
extern int *vis_edge_hwcc;
extern int *q_nodes_hwcc;
// hypergraph.c
void load_hypergraph_from_file(const char* filepath);
void generate_synthetic_hypergraph(void);
void build_mappings(void);
void assign_random_anchors(void);
void save_hypergraph_structure(const char* filename);

// percolation.c
void kosaraju_dfs1(int u);
void kosaraju_dfs2(int u);
void find_all_sccs_simulation(void);
void DFS_GOC(int start_node);
void DFS_GIC(int start_node);
void run_message_passing_probabilistic(double p_N, double p_H, double *pred_gout, double *pred_gin, double *pred_gscc);
void run_single_experiment(double p_N, double p_H, FILE* outfile);
int find_max_hwcc_simulation(void);
void run_message_passing_hwcc(double p_N, double p_H, double *pred_hwcc);

#endif // GLOBALS_H