#include "Graph.h"
#include "Utility.h"
#include "Timer.h"
#include "popl.hpp"

using namespace std;
using namespace popl; 

void print_usage()
{
    printf("Example usage: ./SVC-TD -g [path_to_graph] -l [size_lower_bound] -u [size_upper_bound] -q [query_vertex]\n");
}

int main(int argc, char *argv[])
{
#ifndef NDEBUG
    printf("**** SVC-TD (Debug) build at %s %s ***\n", __TIME__, __DATE__);
    printf("!!! You may want to define NDEBUG in Utility.h to get better performance!\n");
#else
    printf("**** SVC-TD (Release) build at %s %s ***\n", __TIME__, __DATE__);
#endif

    OptionParser op("Allowed options");
    auto help_option = op.add<Switch>("h", "help", "\'produce help message\'");
    auto graph_option = op.add<Value<string>>("g", "graph", "\'path to input graph file\'");
    auto size_lb_option = op.add<Value<int>>("l", "size lower bound", "\'the size lower bound for SVC\'");
    auto size_ub_option = op.add<Value<int>>("u", "size upper bound", "\'the size upper bound for SVC\'");
    auto query_vertex_option = op.add<Value<int>>("q", "query vertex", "\'the input query vertex for SVC\'");
    op.parse(argc, argv);

    if (help_option->is_set() || argc <= 4)
    {
        cout << op << endl;
        if (argc <= 4)
        {
            print_usage();
            return 0;
        }
    }

    if (!graph_option->is_set())
    {
        printf("!!! The argument -g is required! Exit !!!\n");
        print_usage();
        return 0;
    }

    if (!size_lb_option->is_set())
    {
        printf("!!! The argument -l is required! Exit !!!\n");
        print_usage();
        return 0;
    }

    if (!size_ub_option->is_set())
    {
        printf("!!! The argument -u is required! Exit !!!\n");
        print_usage();
        return 0;
    }

    if (!query_vertex_option->is_set())
    {
        printf("!!! The argument -q is required! Exit !!!\n");
        print_usage();
        return 0;
    }

    int LB_size = size_lb_option->value();
    int UB_size = size_ub_option->value();
    int query_vertex = query_vertex_option->value();
    printf("The size lower bound LB is: %d ,the upper bound UB is: %d ,and the query vertex is %d\n", LB_size, UB_size, query_vertex);
    Graph *graph = new Graph(graph_option->value().c_str());
    graph->readBinFile();
    graph->SVC_TD(query_vertex, LB_size, UB_size);

    printf("Everything is OK!!!");
    delete graph;
    return 0;
}
