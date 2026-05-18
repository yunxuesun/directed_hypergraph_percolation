#include "globals.h"

const char* HYPERGRAPH_FILE_PATH = "./data_processing/iJO1366_hypergraph.txt";
const char* output_filename = "result/ph_data.txt";

int N; 
int Q; 


int **hyperedge_all;
int *size_all;
int **node_in_hyperedge;
int **node_in_hyperedge_pos;
int *k_all;


int *dam_edges; 


int *vis_node_hwcc;
int *vis_edge_hwcc;
int *q_nodes_hwcc;
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
int stack_ptr = 0;
int **scc_list;
int scc_count = 0;
int *current_scc;
int current_scc_size = 0;