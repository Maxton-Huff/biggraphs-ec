#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1

void vertex_set_clear(vertex_set* list) {
    list->count = 0;
}

void vertex_set_init(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

void top_down_step_parallel(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{
    // Use OpenMP to parallelize the outer loop
    #pragma omp parallel for
    for (int i = 0; i < frontier->count; i++) {

        int node = frontier->vertices[i];

        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->outgoing_starts[node + 1];

        // Iterate over neighbors of the current node
        for (int neighbor = start_edge; neighbor < end_edge; neighbor++) {
            int outgoing = g->outgoing_edges[neighbor];

            // Use compare-and-swap to safely update the distance array
            if (__sync_bool_compare_and_swap(&distances[outgoing], NOT_VISITED_MARKER, distances[node] + 1)) {
                // Atomically add the node to the new frontier
                int index;
                #pragma omp atomic capture
                index = new_frontier->count++;
                new_frontier->vertices[index] = outgoing;
            }
        }
    }
}

void bfs_top_down(Graph graph, solution* sol) {

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // Initialize all nodes to NOT_VISITED
    for (int i = 0; i < graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // Setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif
    vertex_set_clear(new_frontier);
    // Perform a parallel step
    top_down_step_parallel(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // Swap pointers for the next step
        vertex_set* temp = frontier;
        frontier = new_frontier;
        new_frontier = temp;
    }
}

void bfs_bottom_up(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.
}

void bfs_hybrid(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
}
