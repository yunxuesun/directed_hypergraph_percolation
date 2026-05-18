#include "globals.h"

int main() {
    setvbuf(stdout, NULL, _IONBF, 0);
    srand((unsigned)time(NULL));

    if (USE_REAL_DATA) {
        printf("--- Mode: Loading REAL hypergraph data ---\n");
        load_hypergraph_from_file(HYPERGRAPH_FILE_PATH);
    } else {
        printf("--- Mode: Generating SYNTHETIC hypergraph data ---\n");
        generate_synthetic_hypergraph();
    }

    assign_random_anchors();
    build_mappings();

    printf("Hypergraph setup complete. Starting simulation loop.\n");
    fflush(stdout);
    save_hypergraph_structure("result/hypergraph_structure_anchor.txt");
    
    // Allocate global arrays for simulation
    dam_nodes = malloc(N * sizeof(int));
    dam_edges = malloc(Q * sizeof(int)); 
    vis = calloc(N, sizeof(int));
    size_cluster = calloc(N, sizeof(int));
    stack = malloc(N * sizeof(int));
    scc_list = calloc(N, sizeof(int*));
    current_scc = malloc(N * sizeof(int));

    vis_node_hwcc = calloc(N, sizeof(int));
    vis_edge_hwcc = calloc(Q, sizeof(int));
    q_nodes_hwcc = malloc(N * sizeof(int));
    
    FILE *outfile = fopen(output_filename, "w");
    if (outfile == NULL) {
        fprintf(stderr, "Error: Could not open file for writing.\n");
        return 1;
    }

    fprintf(outfile, "p_N,p_H,mp_hwcc,mp_gout,mp_gin,mp_gscc,sim_hwcc,sim_gout,sim_gin,sim_gscc\n");

    for (int nc = 0; nc <= 100; nc++) {
        double p_N, p_H;
        if (SWEEP_MODE == 1) {
            p_N = nc * 0.01;
            p_H = FIXED_P_H;
        } else if (SWEEP_MODE == 2){
            p_N = FIXED_P_N;
            p_H = nc * 0.01;
        } else {
            p_N = nc * 0.01;
            p_H = p_N;
        }

        run_single_experiment(p_N, p_H, outfile);
        fprintf(stderr, "\rFinished p_N=%.2f, p_H=%.2f", p_N, p_H);
        fflush(stderr);
    }

    fclose(outfile);
    fprintf(stderr, "\nSimulation finished. Data saved to %s.\n", output_filename);
    /* --- Memory Cleanup --- */
    free(dam_nodes);
    free(dam_edges); 
    for (int mu = 0; mu < Q; mu++) { 
        free(hyperedge_T[mu]); free(hyperedge_H[mu]);
        if (hyperedge_all[mu]) free(hyperedge_all[mu]);
        if (hyperedge_anchors[mu]) free(hyperedge_anchors[mu]);
    }
    free(hyperedge_T); free(hyperedge_H);
    free(size_T); free(size_H);
    free(hyperedge_all); free(size_all);
    free(hyperedge_anchors); free(num_anchors);

    for (int i = 0; i < N; i++) {
        if (node_as_input_of[i]) free(node_as_input_of[i]);
        if (node_as_input_of_pos[i]) free(node_as_input_of_pos[i]);
        if (node_as_output_of[i]) free(node_as_output_of[i]);
        if (node_as_output_of_pos[i]) free(node_as_output_of_pos[i]);
        
        if (node_in_hyperedge[i]) free(node_in_hyperedge[i]);
        if (node_in_hyperedge_pos[i]) free(node_in_hyperedge_pos[i]);
    }
    free(node_as_input_of); free(node_as_input_of_pos);
    free(node_as_output_of); free(node_as_output_of_pos);
    free(node_in_hyperedge); free(node_in_hyperedge_pos);
    free(k_in); free(k_out); free(k_all);
    
    free(vis); free(size_cluster); free(stack);
    for(int i=0; i<scc_count; i++) if(scc_list[i]) free(scc_list[i]);
    free(scc_list); free(current_scc);
    
    free(vis_node_hwcc); free(vis_edge_hwcc); free(q_nodes_hwcc);
    return 0;
}