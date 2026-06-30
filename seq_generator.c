#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

#define DEFAULT_N_SYNTHETIC 10000
#define M_SYNTHETIC 3 

#define TOPO_POISSON_UNCOR 0
#define TOPO_POISSON_COR 1
#define TOPO_SF_UNCOR 2
#define TOPO_SF_COR 3

int N, Q;
int N_synthetic = DEFAULT_N_SYNTHETIC;
int **hyperedge_T, **hyperedge_H;
int *size_T, *size_H;
int *k_in, *k_out;

void generate_degree_seq(int *q, int num_elements, int target_sum, int is_sf, double gamma, int max_allowed) {
    if (num_elements <= 0) {
        return;
    }

    if (is_sf == 0) { 
        for(int i = 0; i < num_elements; i++) q[i] = 1;
        int remain = target_sum - num_elements;
        while(remain > 0) { 
            int idx = rand() % num_elements;
            if (q[idx] < max_allowed) { q[idx]++; remain--; }
        }
    } else { 
        double alpha = gamma - 1.0; 
        int cutoff = (int)pow((double)num_elements, 1.0 / alpha);
        if (cutoff < 100) cutoff = 100; 
        if (max_allowed > cutoff) max_allowed = cutoff; 
        
        double *weights = malloc(num_elements * sizeof(double));
        if (weights == NULL) {
            fprintf(stderr, "Failed to allocate degree weights.\n");
            exit(EXIT_FAILURE);
        }
        double w_sum = 0;
        for(int i = 0; i < num_elements; i++) {
            double u = (double)rand() / RAND_MAX;
            if (u >= 0.999999) u = 0.999999;
            weights[i] = pow(1.0 - u, -1.0 / alpha);
            w_sum += weights[i];
        }

        double *cdf = malloc(num_elements * sizeof(double));
        if (cdf == NULL) {
            fprintf(stderr, "Failed to allocate degree CDF.\n");
            free(weights);
            exit(EXIT_FAILURE);
        }
        cdf[0] = weights[0] / w_sum;
        for(int i = 1; i < num_elements; i++) cdf[i] = cdf[i-1] + (weights[i] / w_sum);

        for(int i = 0; i < num_elements; i++) q[i] = 1;
        int remain = target_sum - num_elements;
        int consecutive_fails = 0;

        while(remain > 0) {
            double r = (double)rand() / RAND_MAX;
            int left = 0, right = num_elements - 1, chosen = right;
            while(left <= right) {
                int mid = left + (right - left) / 2;
                if(cdf[mid] >= r) { chosen = mid; right = mid - 1; } 
                else { left = mid + 1; }
            }
            if (q[chosen] < max_allowed) {
                q[chosen]++; remain--; consecutive_fails = 0;
            } else {
                consecutive_fails++;
                if (consecutive_fails > 100) {
                    for (int k = 0; k < num_elements; k++) {
                        if (q[k] < max_allowed) { q[k]++; remain--; break; }
                    }
                }
            }
        }
        free(weights); free(cdf);
    }
}

void generate_synthetic_hypergraph(int topo_mode, double lam, double gamma) {
    N = N_synthetic;
    Q = (int)((lam * N) / M_SYNTHETIC); 
    
    hyperedge_T = malloc(Q * sizeof(int*));
    hyperedge_H = malloc(Q * sizeof(int*));
    size_T = malloc(Q * sizeof(int));
    size_H = malloc(Q * sizeof(int));

    int total_stubs = Q * M_SYNTHETIC;
    int *stubs_in = malloc(total_stubs * sizeof(int));
    int *stubs_out = malloc(total_stubs * sizeof(int));
    int *q_in_seq = calloc(N, sizeof(int));
    int *q_out_seq = calloc(N, sizeof(int));
    if (!hyperedge_T || !hyperedge_H || !size_T || !size_H ||
        !stubs_in || !stubs_out || !q_in_seq || !q_out_seq) {
        fprintf(stderr, "Failed to allocate synthetic hypergraph buffers.\n");
        exit(EXIT_FAILURE);
    }

    int is_sf = (topo_mode == TOPO_SF_UNCOR || topo_mode == TOPO_SF_COR) ? 1 : 0;

    generate_degree_seq(q_in_seq, N, total_stubs, is_sf, gamma, Q);
    if (topo_mode == TOPO_POISSON_COR || topo_mode == TOPO_SF_COR) {
        memcpy(q_out_seq, q_in_seq, N * sizeof(int)); 
    } else {
        generate_degree_seq(q_out_seq, N, total_stubs, is_sf, gamma, Q);
    }

    generate_degree_seq(size_T, Q, total_stubs, is_sf, gamma, N); 
    if (topo_mode == TOPO_POISSON_COR || topo_mode == TOPO_SF_COR) {
        memcpy(size_H, size_T, Q * sizeof(int)); 
    } else {
        generate_degree_seq(size_H, Q, total_stubs, is_sf, gamma, N);
    }

    int idx_in = 0, idx_out = 0;
    for (int i = 0; i < N; i++) {
        for(int k=0; k<q_in_seq[i]; k++) stubs_in[idx_in++] = i;
        for(int k=0; k<q_out_seq[i]; k++) stubs_out[idx_out++] = i;
    }

    int idx_T = 0;
    for (int mu = 0; mu < Q; mu++) {
        hyperedge_T[mu] = malloc(size_T[mu] * sizeof(int));
        for (int t = 0; t < size_T[mu]; t++) {
            int selected_idx = idx_T;
            bool duplicate;
            int attempts = 0;
            do {
                duplicate = false;
                selected_idx = idx_T + rand() % (total_stubs - idx_T);
                int candidate_node = stubs_in[selected_idx];
                for (int k = 0; k < t; k++) {
                    if (hyperedge_T[mu][k] == candidate_node) { duplicate = true; break; }
                }
                attempts++;
                if (attempts > 50 && duplicate) {
                    bool found_any = false;
                    for (int scan = idx_T; scan < total_stubs; scan++) {
                        bool scan_dup = false;
                        for (int k = 0; k < t; k++) {
                            if (hyperedge_T[mu][k] == stubs_in[scan]) { scan_dup = true; break; }
                        }
                        if (!scan_dup) { selected_idx = scan; duplicate = false; found_any = true; break; }
                    }
                    if (!found_any) break; 
                }
            } while (duplicate);
            
            if (duplicate) { size_T[mu] = t; break; }
            int temp = stubs_in[idx_T];
            stubs_in[idx_T] = stubs_in[selected_idx];
            stubs_in[selected_idx] = temp;
            hyperedge_T[mu][t] = stubs_in[idx_T];
            idx_T++;
        }
    }

    int idx_H = 0;
    for (int mu = 0; mu < Q; mu++) {
        hyperedge_H[mu] = malloc(size_H[mu] * sizeof(int));
        for (int h = 0; h < size_H[mu]; h++) {
            int selected_idx = idx_H;
            bool duplicate;
            int attempts = 0;
            do {
                duplicate = false;
                selected_idx = idx_H + rand() % (total_stubs - idx_H);
                int candidate_node = stubs_out[selected_idx];
                for (int k = 0; k < h; k++) {
                    if (hyperedge_H[mu][k] == candidate_node) { duplicate = true; break; }
                }
                attempts++;
                if (attempts > 50 && duplicate) {
                    bool found_any = false;
                    for (int scan = idx_H; scan < total_stubs; scan++) {
                        bool scan_dup = false;
                        for (int k = 0; k < h; k++) {
                            if (hyperedge_H[mu][k] == stubs_out[scan]) { scan_dup = true; break; }
                        }
                        if (!scan_dup) { selected_idx = scan; duplicate = false; found_any = true; break; }
                    }
                    if (!found_any) break; 
                }
            } while (duplicate);
            
            if (duplicate) { size_H[mu] = h; break; }
            int temp = stubs_out[idx_H];
            stubs_out[idx_H] = stubs_out[selected_idx];
            stubs_out[selected_idx] = temp;
            hyperedge_H[mu][h] = stubs_out[idx_H];
            idx_H++;
        }
    }

    // 计算实际的 k_in 和 k_out
    k_in = calloc(N, sizeof(int));
    k_out = calloc(N, sizeof(int));
    for (int mu = 0; mu < Q; mu++) {
        for (int t = 0; t < size_T[mu]; t++) k_out[hyperedge_T[mu][t]]++;
        for (int h = 0; h < size_H[mu]; h++) k_in[hyperedge_H[mu][h]]++;
    }

    free(stubs_in); free(stubs_out); free(q_in_seq); free(q_out_seq);
}

void generate_synthetic_sequences(int topo_mode, double lam, double gamma) {
    N = N_synthetic;
    Q = (int)((lam * N) / M_SYNTHETIC);

    int total_stubs = Q * M_SYNTHETIC;
    k_in = calloc(N, sizeof(int));
    k_out = calloc(N, sizeof(int));
    size_T = calloc(Q, sizeof(int));
    size_H = calloc(Q, sizeof(int));
    if (!k_in || !k_out || !size_T || !size_H) {
        fprintf(stderr, "Failed to allocate sequence buffers.\n");
        exit(EXIT_FAILURE);
    }

    int is_sf = (topo_mode == TOPO_SF_UNCOR || topo_mode == TOPO_SF_COR) ? 1 : 0;

    generate_degree_seq(k_out, N, total_stubs, is_sf, gamma, Q);
    if (topo_mode == TOPO_POISSON_COR || topo_mode == TOPO_SF_COR) {
        memcpy(k_in, k_out, N * sizeof(int));
    } else {
        generate_degree_seq(k_in, N, total_stubs, is_sf, gamma, Q);
    }

    generate_degree_seq(size_T, Q, total_stubs, is_sf, gamma, N);
    if (topo_mode == TOPO_POISSON_COR || topo_mode == TOPO_SF_COR) {
        memcpy(size_H, size_T, Q * sizeof(int));
    } else {
        generate_degree_seq(size_H, Q, total_stubs, is_sf, gamma, N);
    }
}

int main(int argc, char *argv[]) {
    if (argc != 4 && argc != 5) return 1;
    int topo = atoi(argv[1]);
    double lam = atof(argv[2]);
    double gamma = atof(argv[3]);
    if (argc == 5) {
        N_synthetic = atoi(argv[4]);
        if (N_synthetic <= 0) {
            fprintf(stderr, "N must be positive.\n");
            return 1;
        }
    }

    // 防止在短时间内密集调用时产生相同种子
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    srand((unsigned)(ts.tv_nsec ^ ts.tv_sec));

    generate_synthetic_sequences(topo, lam, gamma);

    // 格式化输出到 stdout 给 Python 读取
    printf("%d %d\n", N, Q);
    for(int i = 0; i < N; i++) printf("%d %d\n", k_in[i], k_out[i]);
    for(int mu = 0; mu < Q; mu++) printf("%d %d\n", size_T[mu], size_H[mu]);

    // 释放内存
    free(size_T); free(size_H);
    free(k_in); free(k_out);
    
    return 0;
}
