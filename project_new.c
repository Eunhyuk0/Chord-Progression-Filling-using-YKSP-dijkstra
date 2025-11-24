#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VERTICES 24
#define SCALENUM 7 //nodes for each graph
#define MAX_LAYERS 16
#define MAX_NODE_LAYER (SCALENUM * MAX_LAYERS)
#define SRC_NODE (MAX_NODE_LAYER)
#define SINK_NODE (MAX_NODE_LAYER + 1)
#define NODE_NUM (MAX_NODE_LAYER + 2) // + source + sink
#define MAX_EDGES 2000
#define TOP_K 10
#define INF_COST 100000000


typedef struct {
    int to;
    int weight; // inverted
    int orig_from; // absolute chord id
    int orig_to;
    int active; // 1 if active
} Edge;

Edge edges[MAX_EDGES];
int edge_count = 0;

// adjacency list for layered graph
int adj_count[NODE_NUM]; //number of adjacent edges for each nodes
int adj_list[NODE_NUM][SCALENUM * MAX_LAYERS];

// path struct
typedef struct {
    int nodes[NODE_NUM]; // (only first 'len' of whole array)
    int len;
    int total_cost; // sum of converted weights
} Path;

Path A[TOP_K+5]; // found shortest paths (with dijkstra)
int A_count = 0;
Path B[MAX_EDGES]; // all possible paths in Yen's
int B_count = 0;

struct valuegraph {
    int num;
    int matrix[VERTICES][VERTICES];
} WeightGraph;


int Layers = 4; // default

void init_valuegraph(int num) {
    WeightGraph.num = num;
    for(int i=0;i<num;i++) 
        for(int j=0;j<num;j++) 
            WeightGraph.matrix[i][j] = 0;
}

void valuegraph_insert_edge(int v1,int v2,int weight) { 
    WeightGraph.matrix[v1][v2] = weight; 
}

void get_weight(int quality, int style)
{
     //constants
    int maj_chords[7] = {0, 14, 16, 5, 7, 21, 23};
    int min_chords[7] = {12, 14, 4, 17, 7, 9, 11};
    int min_natural_chords[7] = {12, 14, 4, 17, 19, 9, 11};

    int maj_weight[7][7] =
    {
        {5, 10, 5, 30, 25, 20, 5},
        {10, 5, 5, 20, 45, 10, 5},
        {10, 15, 5, 20, 5, 40, 5},
        {30, 15, 4, 3, 35, 8, 5},
        {60, 5, 3, 10, 2, 15, 5},
        {10, 30, 5, 25, 20, 5, 5},
        {55, 5, 10, 5, 20, 3, 2},
    };
    int min_weight[7][7] =
    {
        {5, 5, 10, 30, 20, 20, 10},
        {25, 5, 10, 10, 40, 5, 5},
        {35, 5, 5, 20, 5, 20, 10},
        {30, 5, 5, 10, 30, 10, 10},
        {50, 5, 5, 15, 5, 15, 5},
        {25, 5, 10, 25, 0, 5, 30},
        {40, 5, 20, 15, 10, 10, 0},
    };
    int jazz_maj_weight[7][7] =
    {
        {10, 25, 10, 40, 35, 30, 5},  // I
        {60, 10, 15, 10, 90, 10, 10}, // ii
        {15, 30, 10, 20, 60, 50, 10}, // iii
        {70, 40, 20, 10, 80, 15, 15}, // IV
        {95, 40, 10, 30, 10, 10, 20}, // V
        {40, 75, 20, 30, 60, 10, 10}, // vi
        {85, 10, 10, 10, 90, 15, 5}   // vii°
    };
    int jazz_min_weight[7][7] =
    {
        {10, 25, 15, 40, 50, 35, 10}, // i
        {60, 10, 15, 10, 90, 10, 10}, // iiø
        {20, 20, 10, 30, 40, 60, 10}, // III
        {75, 20, 15, 10, 85, 20, 10}, // iv
        {95, 30, 10, 30, 10, 20, 15}, // V7
        {50, 60, 20, 20, 70, 10, 10}, // VI
        {80, 10, 10, 10, 85, 15, 5}   // vii°
    };

    if(quality) //maj
    {
        if(style)
        {
            for(int i=0;i<7;i++) 
            {
                for(int j=0;j<7;j++)
                    valuegraph_insert_edge(maj_chords[i], maj_chords[j], jazz_maj_weight[i][j]);
            }
        }
        else
        {
            for(int i=0;i<7;i++) 
            {
                for(int j=0;j<7;j++)
                    valuegraph_insert_edge(maj_chords[i], maj_chords[j], maj_weight[i][j]);
            }
        }
        
    } 
    else  //min
    {
        if(style)
        {
            for(int i=0;i<7;i++) 
            {
                for(int j=0;j<7;j++)
                    valuegraph_insert_edge(min_chords[i], min_chords[j], jazz_min_weight[i][j]);
            }
        }
        else
        {
            for(int i=0;i<7;i++) 
            {
                for(int j=0;j<7;j++)
                    valuegraph_insert_edge(min_chords[i], min_chords[j], min_weight[i][j]);
            }
        }
        
    }
}

void init_graph() {
    edge_count = 0;
    for(int i=0;i<NODE_NUM;i++) {
        adj_count[i] = 0;
        for(int j = 0;j<SCALENUM*MAX_LAYERS;j++)
            adj_list[i][j] = 0;
    }
}

void insert_edge(int u, int v, int w, int orig_from, int orig_to) {
    if(edge_count >= MAX_EDGES) {
        printf("\nerror - edge storage overflow\n");
    }
    edges[edge_count].to = v;
    edges[edge_count].weight = w;
    edges[edge_count].orig_from = orig_from;
    edges[edge_count].orig_to = orig_to;
    edges[edge_count].active = 1;
    adj_list[u][adj_count[u]++] = edge_count;
    edge_count++;
}

int dijkstra_path(int source, int sink, int prev_out[], int dist_out[]) {
    //simple djikstra's
    int N = NODE_NUM;
    int visited[NODE_NUM];
    int dist[NODE_NUM];//distance(cost) from SRC to i
    int prev[NODE_NUM];//nodes passed from SRC to i
    int eidx;
    int u,v,w;
    int best = INF_COST;

    for(int i=0;i<N;i++){ visited[i]=0; dist[i]=INF_COST; prev[i]=-1; }
    dist[source] = 0;
    //dijkstra's
    while(1) {
        u=-1; 
        best = INF_COST;
        for(int i=0;i<N;i++){
            if(!visited[i] && dist[i] < best){ best = dist[i]; u=i; }
        }
        if(u==-1 || u == sink) 
            break; //no reachable nodes or no paths
        visited[u] = 1;
        for(int ei_idx=0; ei_idx < adj_count[u]; ei_idx++){
            eidx = adj_list[u][ei_idx];
            if(!edges[eidx].active) continue; //pass if deactivated
            v = edges[eidx].to;
            w = edges[eidx].weight;
            if(dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                prev[v] = u;
            }
        }
    }
    if(dist[sink] >= INF_COST) 
        return 0; //fail

    for(int i=0;i<N;i++)
    { 
        prev_out[i] = prev[i]; //outputs
        dist_out[i] = dist[i]; 
    }
    return 1;
} //run dijkstra's and fill path_nodes and path_len

int build_path_from_prev(int source, int sink, int prev[], int path_nodes[]) {
    int tmp[NODE_NUM];
    int tlen = 0;
    int cur = sink;
    int len;
    while(cur != -1 && cur != source) {
        tmp[tlen++] = cur;
        cur = prev[cur];
    }
    if(cur == -1) 
        return 0; //no path
    tmp[tlen++] = source;
    //reverse
    len = tlen;
    for(int i=0;i<len;i++) path_nodes[i] = tmp[len-1-i];
    return len;
}

int find_edge_index_between(int u, int v) {
    int eidx;
    for(int ei=0; ei < adj_count[u]; ei++){
        eidx = adj_list[u][ei];
        if(edges[eidx].to == v) return eidx;
    }
    return -1;
} //number of edges between 2 elements

void add_candidate_to_B(Path *p) {
    int same;
    //without duplicates
    for(int i=0;i<B_count;i++){
        if(B[i].len == p->len){
            same = 1;
            for(int j=0;j<p->len;j++) 
                if(B[i].nodes[j] != p->nodes[j]) { same = 0; break; }
            if(same) return;
        }
    }
    B[B_count++] = *p;
}

int pop_best_from_B(Path *out) {
    int best_i;
    if(B_count==0) return 0;
    best_i = 0;
    for(int i=1;i<B_count;i++) 
        if(B[i].total_cost < B[best_i].total_cost) best_i = i;
    *out = B[best_i];
    for(int k=best_i+1;k<B_count;k++) B[k-1] = B[k];
    B_count--;
    return 1;
} //remove the shortest path from B[]

void build_layered_graph(int diatonic_chords[], int chordinput[]) {
    int node;
    int u,v;
    int to, from, w;
    int conv;
    init_graph();
    // SRC to L0
    for(int i=0;i<SCALENUM;i++){
        int node = 0 * SCALENUM + i;
        if(chordinput[0] != 24 && diatonic_chords[i] != chordinput[0]) 
            continue; //first one fixed
        insert_edge(SRC_NODE, node, 0, -1, -1);
    }
    // Edges between layers
    for(int layer=0; layer<Layers-1; layer++){
        for(int i=0;i<SCALENUM;i++){
            u = layer * SCALENUM + i;
            //skip if chord is fixed (input 24)
            if(chordinput[layer] != 24 && diatonic_chords[i] != chordinput[layer]) 
                continue;
            for(int j=0;j<SCALENUM;j++){
                v = (layer+1) * SCALENUM + j;
                if(chordinput[layer+1] != 24 && diatonic_chords[j] != chordinput[layer+1]) 
                    continue;
                from = diatonic_chords[i];
                to = diatonic_chords[j];
                //algorithmic cost : 100-weight
                w = WeightGraph.matrix[from][to];
                conv = 100 - w;
                if(conv < 0) conv = 0;
                insert_edge(u, v, conv, from, to);
            }
        }
    } 
    // Last layer to sink
    for(int i=0;i<SCALENUM;i++){
        int node = (Layers-1)*SCALENUM + i;
        if(chordinput[Layers-1] != 24 && diatonic_chords[i] != chordinput[Layers-1]) 
            continue;
        insert_edge(node, SINK_NODE, 0, -1, -1);
    }
} //implement edges[], adj_list[][] for each layers (if input != 24)

int prefix_equal(int *a, int *b, int len) {
    for(int i=0;i<len;i++) 
    {
        if(a[i] != b[i]) 
            return 0;
    }
    return 1;
} //for 'len' length if 2 arrays are equal return 1 and if not 0


int yen_k_shortest(int K) {
    A_count = 0; B_count = 0;
    int dis_count = 0; //num of disabled edges
    int rootPath[NODE_NUM], spur_nodes[NODE_NUM];
    int prev[NODE_NUM], dists[NODE_NUM]; //result of dijkstra's
    int path_nodes[NODE_NUM]; //temp
    int plen,spurNode,rootPathLen,splen;
    int u,v,eidx,eidx2,idx,total_cost;
    int disabled_edges_idx[MAX_EDGES];
    Path best;
    Path candidate;
    int p_spur[NODE_NUM], d_spur[NODE_NUM]; //result of dijkstra from spurNode

    //find shortest path from SRC to SINK
    if(!dijkstra_path(SRC_NODE, SINK_NODE, prev, dists)) 
        return 0;
    plen = build_path_from_prev(SRC_NODE, SINK_NODE, prev, path_nodes);
    if(plen == 0) 
        return 0;
    //store in A
    A[0].len = plen;
    for(int i=0;i<plen;i++) 
        A[0].nodes[i] = path_nodes[i];
    A[0].total_cost = dists[SINK_NODE];
    A_count = 1;

    //main loop : find k-1 paths
    for(int k=1; k<K; k++){
        //for each node in previous best path (without SINK)
        for(int i=0;i<A[k-1].len - 1; i++)
        {
            spurNode = A[k-1].nodes[i];
            rootPathLen = i+1;
            for(int r=0;r<rootPathLen;r++) 
                rootPath[r] = A[k-1].nodes[r];

            //for each path p in A (with paths with same rootPath)
            for(int p_i=0; p_i<A_count; p_i++)
            {
                if(A[p_i].len >= rootPathLen && prefix_equal(A[p_i].nodes, rootPath, rootPathLen))
                {
                    //remove the edge that follows root path
                    u = A[p_i].nodes[rootPathLen-1];
                    v = A[p_i].nodes[rootPathLen];
                    eidx = find_edge_index_between(u,v);
                    if(eidx>=0 && edges[eidx].active){
                        edges[eidx].active = 0;
                        disabled_edges_idx[dis_count++] = eidx;
                    }
                }
            }
            //compute spur path from spurNode to SINK (with edge removed)
            if(dijkstra_path(spurNode,SINK_NODE,p_spur,d_spur))
            {
                splen = build_path_from_prev(spurNode,SINK_NODE,p_spur,spur_nodes);
                if(splen>0){
                    //total path : rootPath(w/o last node) + spur_nodes
                    idx=0;
                    for(int r=0;r<rootPathLen-1;r++) 
                        candidate.nodes[idx++]=rootPath[r];
                    for(int r=0;r<splen;r++) 
                        candidate.nodes[idx++]=spur_nodes[r];
                    candidate.len=idx;

                    total_cost = 0;
                    for(int rr=0; rr<candidate.len-1; rr++){
                        u=candidate.nodes[rr];
                        v=candidate.nodes[rr+1];
                        eidx2=find_edge_index_between(u,v);
                        if (eidx2 >= 0)
                            total_cost+=edges[eidx2].weight;
                        else 
                            total_cost += INF_COST / 10; //error case
                    }
                    candidate.total_cost=total_cost;
                    add_candidate_to_B(&candidate);
                }
            }
            //restore removed edges
            for(int di=0;di<dis_count;di++) 
                edges[disabled_edges_idx[di]].active = 1;
        } //end for each spur node
        // pick best from B
        if(!pop_best_from_B(&best)) 
            break;
        //add best to A
        A[A_count++] = best;
    }//end for k
    return A_count;
} // Yen’s algorithm

void print_top_results(int diatonic_chords[]) {
    int limit; //number of results
    int seq[MAX_LAYERS];
    int seq_idx;
    int weightsum;
    int node;
    int idx;

    if(A_count == 0) {
        printf("\nerror - No paths found.\n");
        return;
    }
    // compute musical weight sums and print top k paths
    if(A_count < TOP_K)
        limit = A_count;
    else
        limit = TOP_K;
    printf("\nTop %d progressions (format: each chord per layer):\n\n", limit);
    for(int i=0;i<limit;i++){
        Path *p=&A[i]; //absolute chord id seq by mapping
       seq_idx=0;
        for(int j=0;j<p->len;j++){
            node = p->nodes[j];
            if(node==SRC_NODE || node==SINK_NODE) 
                continue;
            idx = node % SCALENUM;
            seq[seq_idx++] = diatonic_chords[idx];
        }
        weightsum=0;
        for(int t=0;t<seq_idx-1;t++)
            weightsum += WeightGraph.matrix[seq[t]][seq[t+1]];
        printf("%2d) ", i+1);
        for(int t=0;t<seq_idx;t++) printf("%2d ", seq[t]);
            printf("   (weight sum=%d, alg cost=%d)\n", weightsum, p->total_cost);
    }
} //print top K results

int main()
{
    int key;
    char tempinput;
    int quality; //1:maj, 0:min
    int style; //1 for jazz
    int diatonic_chords[7];
    int chordinput[MAX_LAYERS];
    int flag_diatonic = 0; //flag for checking if user's input chord fits the selected diatonic scale

    init_graph();
    //constants
    int maj_chords[7] = {0, 14, 16, 5, 7, 21, 23};
    int min_chords[7] = {12, 14, 4, 17, 7, 9, 11};

    printf("\nPick key (0=C ... 11=B):\t");
    scanf("%d", &key);
    if(key<0 || key>=12)
    { 
        printf("\nerror - invalid key\n"); 
        return 0; 
    }

    printf("\nMajor for M, minor for m:\t");
    scanf(" %c", &tempinput);
    if(tempinput == 'M') 
        quality = 1;
    else if(tempinput == 'm')
    { 
        quality = 0;
    }
    else
    {
        printf("\nerror - type in 'M' or 'm'");
        return 0;
    }

    printf("\nwhich harmony theory to follow: c for classical, j for jazz/blues\t");
    scanf(" %c", &tempinput);
    if(tempinput == 'c')
        style = 0;
    else if(tempinput == 'j')
        style = 1;
    else
    {
        printf("\nerror - type in c or j");
        return 0;
    }

    printf("\nChord progression length (2~%d):\t", MAX_LAYERS);
    scanf("%d", &Layers);
    if(Layers < 2 || Layers > MAX_LAYERS){
        printf("\nerror - invalid length value\n");
        return 0;
    }

    // build diatonic chords array with absolute chord ids
    if(quality) 
    {
        for(int i=0;i<7;i++) 
            diatonic_chords[i]=maj_chords[i];
    }
    else 
    {
        for(int i=0;i<7;i++) 
            diatonic_chords[i]=min_chords[i];
    }
    get_weight(quality, style);

    printf("\nType chords (0:I, 14:ii, ..., 24=blank) - total %d:\n\t", Layers);
    for(int i=0;i<Layers;i++)
    {
        scanf("%d", &chordinput[i]);
        flag_diatonic=0;
        for(int j=0;j<7;j++)
        {
            if(chordinput[i]==24 || diatonic_chords[j]==chordinput[i])
            {
                flag_diatonic=1; 
                break;
            }
        }
        if(flag_diatonic==0)
        {
            printf("\nerror - chord %d not in scale", i+1);
            return 0;
        }
    }

    build_layered_graph(diatonic_chords, chordinput);

    //run yen's
    int found = yen_k_shortest(TOP_K);
    if(found==0){
        printf("\nNo valid paths found.\n");
        return 0;
    }

    print_top_results(diatonic_chords);
    return 0;
}
