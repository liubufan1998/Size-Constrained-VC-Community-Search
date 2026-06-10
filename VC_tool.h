#ifndef _VC_TOOL_H
#define _VC_TOOL_H

#include <stdio.h>
#include <utility>
#include <fstream>
#include <bits/stdc++.h>

const unsigned int INF = 0x3f3f3f3f;
using ui = unsigned int;  // vertex type
using ept = unsigned int; // edge pointer type; unsigned int can be used to process upto two billion undirected edges
using namespace std;

class MaximumFlow 
{
private:
    struct Node
    {
        ui from, to, next;
        ui cap;
    };
    Node *edge;
    ui *cap_backup;
    ui tol;
    ui *head;
    ui *dep;
    ui *gap, *que;
    ui *cur;
    ui *tau;
    ui n;
    ui m;
    ui *nV, *oID;
    bool flag; 
    vector<bool> visited;
    vector<ui> real_minCut; 

public:
    MaximumFlow();
    ~MaximumFlow();
    ui SAP(ui start, ui end, ui nn, vector<ui> &real_minCut);
    void BFS(ui start, ui end);
    void reserve(ui n, ui m);
    void init(ui _n);
    void addedge(ui u, ui v, ui w);
    void BFSMinCut(ui start, vector<bool> &visited);
    bool satisfy_necessary_requirement(ui *ids, ui C_end, char *matrix, long long row, int tau, bool is_global);                                     // 在邻接矩阵中判断一个图是否是s-bundle，但是图顶点数组ids是ui*形式
    bool satisfy_necessary_requirement(ui *ids, ui C_end, ui *pstart, ui *pend, ui *edges, int tau);                                                 // 在邻接表中判断一个图是否是s-bundle，参数pointer用来指示是判断S_end还是R_end                       // 在邻接表中判断一个图是否是s-bundle，参数pointer用来指示是判断S_end还是R_end
    ui get_exact_vertex_connectivity(vector<ui> &ids, ui size, ui *pstart, ui *edges, ui current_VC);                                                  // 用来精确的获取一个图的顶点连通度
    void get_branching_case_based_on_OSP(ui *ids, ui C_end, char *matrix, vector<ui> &s_connected, long long row, ui tau, ui &u, ui &v);                         // 用来在BB_matrix分支过程中得到具体的分支Scenario                                                              // verify_Sbundle_by_maxFlowAlg2的变体
    void get_branching_case_based_on_OSP(ui *SC, ui &u1, ui &u2, vector<char> &s_connected, ui S_end, ui C_end, ui *pstart, ui *pend, ui *edges, int tau, ui n); // 用来在邻接表形式中找到基于局部顶点连通度的多分支顶点对
};
#endif