#include <stdio.h>
#include "VC_tool.h"

MaximumFlow::MaximumFlow()
{
    edge = nullptr;
    cap_backup = nullptr;
    head = nullptr;
    dep = nullptr;
    gap = nullptr;
    que = nullptr;
    cur = nullptr;
    tau = nullptr;
    nV = nullptr; 
    oID = nullptr;
}

MaximumFlow::~MaximumFlow()
{
    if (edge != nullptr)
    {
        delete[] edge;
        edge = nullptr;
    }
    if (cap_backup != nullptr)
    {
        delete[] cap_backup;
        cap_backup = nullptr;
    }
    if (head != nullptr)
    {
        delete[] head;
        head = nullptr;
    }
    if (dep != nullptr)
    {
        delete[] dep;
        dep = nullptr;
    }
    if (gap != nullptr)
    {
        delete[] gap;
        gap = nullptr;
    }
    if (que != nullptr)
    {
        delete[] que;
        que = nullptr;
    }
    if (cur != nullptr)
    {
        delete[] cur;
        cur = nullptr;
    }
    if (tau != nullptr)
    {
        delete[] tau;
        tau = nullptr;
    }
    if (nV != nullptr)
    {
        delete[] nV;
        nV = nullptr;
    }
    if (oID != nullptr)
    {
        delete[] oID;
        oID = nullptr;
    }
}

void MaximumFlow::BFS(ui start, ui end)
{
    memset(dep, -1, n * sizeof(ui));
    memset(gap, 0, n * sizeof(ui));
    gap[0] = 1;
    ui front, rear;
    front = rear = 0;
    dep[end] = 0;
    que[rear++] = end;
    while (front != rear)
    {
        ui u = que[front++];
        for (ui i = head[u]; i != -1; i = edge[i].next)
        {
            ui v = edge[i].to;
            if (dep[v] != -1)
                continue;
            que[rear++] = v;
            dep[v] = dep[u] + 1;
            ++gap[dep[v]];
        }
    }
}

void MaximumFlow::reserve(ui n, ui m)
{
    n += 3;
    edge = new Node[m * 4 + n * 2];
    cap_backup = new ui[m * 4 + n * 2];
    head = new ui[n + n];
    dep = new ui[n + n];
    gap = new ui[n + n];
    que = new ui[n + n];
    cur = new ui[n + n];
    tau = new ui[n + n];
}

void MaximumFlow::init(ui _n)
{
    tol = 0;
    n = _n + 2;
    memset(head, -1, n * sizeof(ui));
}

void MaximumFlow::addedge(ui u, ui v, ui w)
{

    edge[tol].from = u;
    edge[tol].to = v;
    edge[tol].cap = w;
    cap_backup[tol] = w;
    edge[tol].next = head[u];
    head[u] = tol++;

    edge[tol].from = v;
    edge[tol].to = u;
    edge[tol].cap = 0;
    cap_backup[tol] = 0;
    edge[tol].next = head[v];
    head[v] = tol++;
}

ui MaximumFlow::SAP(ui start, ui end, ui new_size, vector<ui> &real_minCut)
{
    real_minCut.clear();

    for (ui i = 0; i < tol; ++i)
        edge[i].cap = cap_backup[i];
    ui res = 0;
    BFS(start, end);

    ui top = 0;
    memcpy(cur, head, n * sizeof(ui));
    ui u = start;
    ui i;
    while (dep[start] < n)
    {
        if (u == end)
        {
            ui temp = INF + 1;
            ui inser = 0;
            for (i = 0; i < top; ++i)
                if (temp > edge[tau[i]].cap)
                {
                    temp = edge[tau[i]].cap;
                    inser = i;
                }

            for (i = 0; i < top; ++i)
            {
                edge[tau[i]].cap -= temp;
                edge[tau[i] ^ 1].cap += temp;
            }
            res += temp;
            top = inser;
            u = edge[tau[top]].from;
        }
        if (u != end && gap[dep[u] - 1] == 0)
            break;
        for (i = cur[u]; i != -1; i = edge[i].next)
            if (edge[i].cap != 0 && dep[u] == dep[edge[i].to] + 1)
                break;

        if (i != -1)
        {
            cur[u] = i;
            tau[top++] = i;
            u = edge[i].to;
        }
        else
        {
            ui min = n;
            for (i = head[u]; i != -1; i = edge[i].next)
            {
                if (edge[i].cap == 0)
                    continue;
                if (min > dep[edge[i].to])
                {
                    min = dep[edge[i].to];
                    cur[u] = i;
                }
            }
            --gap[dep[u]];
            dep[u] = min + 1;
            ++gap[dep[u]];

            if (u != start)
                u = edge[tau[--top]].from;
        }
    }

    vector<bool> visited(n, false);
    BFSMinCut(start, visited);
    for (ui i = 0; i < tol; i += 2)
    {
        ui u = edge[i].from;
        ui v = edge[i].to;

        if (visited[u] && !visited[v])
        {
            assert(u < v);
            real_minCut.push_back(u);
        }
    }
    assert(real_minCut.size() == res);
    return res;
}

void MaximumFlow::BFSMinCut(ui start, vector<bool> &visited)
{
    queue<ui> q;
    q.push(start);
    visited[start] = true;
    while (!q.empty())
    {
        ui u = q.front();
        q.pop();
        for (ui i = head[u]; i != -1; i = edge[i].next)
        {
            ui v = edge[i].to;
            if (!visited[v] && edge[i].cap > 0)
            {
                visited[v] = true;
                q.push(v);
            }
        }
    }
}

bool MaximumFlow::satisfy_necessary_requirement(ui *SC, ui C_end, char *matrix, long long row, int tau, bool is_global)
{
    if (C_end <= tau)
        return true;

    assert(C_end > tau);

    MaximumFlow mf;
    mf.reserve(C_end, C_end * C_end);
    nV = new ui[C_end + 3];
    oID = new ui[C_end + 3];

    int new_size = 0;
    for (ui i = 0; i < C_end; i++)
        nV[i] = new_size, oID[new_size++] = SC[i];

    mf.init(new_size + new_size);

    for (ui i = 0; i < C_end; i++)
        mf.addedge(nV[i], nV[i] + new_size, 1);

    for (ui i = 0; i < C_end; i++)
    {
        ui u = SC[i];
        char *t_matrix = matrix + u * row;
        for (ui j = i + 1; j < C_end; j++)
        {
            ui v = SC[j];
            if (!t_matrix[v])
                continue;
            mf.addedge(nV[i] + new_size, nV[j], INF);
            mf.addedge(nV[j] + new_size, nV[i], INF);
        }
    }

    bool is_sbundle = true;
    if (is_global)
    {
        for (ui i = 0; i < new_size; ++i)
        {
            if (is_sbundle)
            {
                char *t_matrix = matrix + oID[i] * row;
                for (ui j = i + 1; j < new_size; ++j)
                    if (!t_matrix[oID[j]] && mf.SAP(new_size + i, j, new_size, real_minCut) < C_end - tau)
                    {
                        is_sbundle = false;
                        break;
                    }
            }
        }
    }
    else
    {
        ui u = oID[new_size - 1];
        assert(u == SC[C_end - 1]);
        char *t_matrix = matrix + u * row;
        for (ui i = 0; i < new_size - 1; ++i)
        {
            assert(oID[i] != u);
            if (!t_matrix[oID[i]] && mf.SAP(new_size + new_size - 1, i, new_size, real_minCut) < C_end - tau)
            {
                is_sbundle = false;
                break;
            }
        }
    }

    if (nV != nullptr)
    {
        delete[] nV;
        nV = nullptr;
    }
    if (oID != nullptr)
    {
        delete[] oID;
        oID = nullptr;
    }
    return is_sbundle;
}

ui MaximumFlow::get_exact_vertex_connectivity(vector<ui> &SC, ui C_end, ui *pstart, ui *edges, ui current_VC)
{

    MaximumFlow mf;
    mf.reserve(C_end, C_end * C_end);
    nV = new ui[C_end + 3];
    oID = new ui[C_end + 3];

    int new_size = 0;
    for (ui i = 0; i < C_end; i++)
        nV[i] = new_size, oID[new_size++] = SC[i];

    mf.init(new_size + new_size);
    for (ui i = 0; i < C_end; i++)
        mf.addedge(nV[i], nV[i] + new_size, 1);

    for (ui i = 0; i < C_end; i++)
    {
        ui u = SC[i];
        for (ui j = i + 1; j < C_end; j++)
        {
            ui v = SC[j];

            bool is_adj = false;
            for (ui k = pstart[u]; k < pstart[u + 1]; k++)
            {
                if (edges[k] == v)
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
                continue;
            mf.addedge(nV[i] + new_size, nV[j], INF);
            mf.addedge(nV[j] + new_size, nV[i], INF);
        }
    }

    ui min_VC = C_end - 1;
    for (ui i = 0; i < new_size; ++i)
        for (ui j = i + 1; j < new_size; ++j)
        {

            bool is_adj = false;
            for (ui k = pstart[oID[i]]; k < pstart[oID[i] + 1]; k++)
            {
                if (edges[k] == oID[j])
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
            {
                ui t_VC = mf.SAP(new_size + i, j, new_size, real_minCut);
                if (t_VC < min_VC)
                    min_VC = t_VC;
                if (min_VC <= current_VC)
                    break;
            }
        }

    if (nV != nullptr)
    {
        delete[] nV;
        nV = nullptr;
    }
    if (oID != nullptr)
    {
        delete[] oID;
        oID = nullptr;
    }
    return min_VC;
}

bool MaximumFlow::satisfy_necessary_requirement(ui *SC, ui C_end, ui *pstart, ui *pend, ui *edges, int tau)
{
    if (C_end <= tau)
        return true;

    MaximumFlow mf;
    mf.reserve(C_end, C_end * C_end);
    nV = new ui[C_end + 3];
    oID = new ui[C_end + 3];

    int new_size = 0;
    for (ui i = 0; i < C_end; i++)
        nV[i] = new_size, oID[new_size++] = SC[i];

    mf.init(new_size + new_size);
    for (ui i = 0; i < C_end; i++)
        mf.addedge(nV[i], nV[i] + new_size, 1);

    for (ui i = 0; i < C_end; i++)
    {
        ui u = SC[i];
        for (ui j = i + 1; j < C_end; j++)
        {
            ui v = SC[j];

            bool is_adj = false;
            for (ui k = pstart[u]; k < pend[u]; k++)
            {
                if (edges[k] == v)
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
                continue;
            mf.addedge(nV[i] + new_size, nV[j], INF);
            mf.addedge(nV[j] + new_size, nV[i], INF);
        }
    }

    bool is_sbundle = true;
    for (ui i = 0; i < new_size; ++i)
        for (ui j = i + 1; j < new_size; ++j)
        {
            if (is_sbundle)
            {
                bool is_adj = false;
                for (ui k = pstart[oID[i]]; k < pend[oID[i]]; k++)
                {
                    if (edges[k] == oID[j])
                    {
                        is_adj = true;
                        break;
                    }
                }
                if (!is_adj && mf.SAP(new_size + i, j, new_size, real_minCut) < C_end - tau)
                {
                    is_sbundle = false;
                    break;
                }
            }
        }

    if (nV != nullptr)
    {
        delete[] nV;
        nV = nullptr;
    }
    if (oID != nullptr)
    {
        delete[] oID;
        oID = nullptr;
    }
    return is_sbundle;
}

void MaximumFlow::get_branching_case_based_on_OSP(ui *SC, ui C_end, char *matrix, vector<ui> &HG_mark, long long row, ui tau, ui &v1, ui &v2)
{
    if (C_end <= tau)
        return;

    assert(C_end > tau);

    MaximumFlow mf;
    mf.reserve(C_end, C_end * C_end);
    nV = new ui[C_end + 3];
    oID = new ui[C_end + 3];

    int new_size = 0;
    for (ui i = 0; i < C_end; i++)
        nV[i] = new_size, oID[new_size++] = SC[i];

    mf.init(new_size + new_size);

    for (ui i = 0; i < C_end; i++)
        mf.addedge(nV[i], nV[i] + new_size, 1);

    for (ui i = 0; i < C_end; i++)
    {
        ui u = SC[i];
        char *t_matrix = matrix + u * row;
        for (ui j = i + 1; j < C_end; j++)
        {
            ui v = SC[j];
            if (!t_matrix[v])
                continue;
            mf.addedge(nV[i] + new_size, nV[j], INF);
            mf.addedge(nV[j] + new_size, nV[i], INF);
        }
    }

    for (ui i = 0; i < new_size; ++i)
    {
        char *t_matrix = matrix + oID[i] * row;
        for (ui j = i + 1; j < new_size; ++j)
            if (!t_matrix[oID[j]] && mf.SAP(new_size + i, j, new_size, real_minCut) < C_end - tau)
            {
                v1 = oID[i];
                v2 = oID[j];
                HG_mark.resize(row);
                for (ui k = 0; k < real_minCut.size(); k++)
                    HG_mark[oID[real_minCut[k]]] = 1;
                if (nV != nullptr)
                {
                    delete[] nV;
                    nV = nullptr;
                }
                if (oID != nullptr)
                {
                    delete[] oID;
                    oID = nullptr;
                }
                return;
            }
    }
}

void MaximumFlow::get_branching_case_based_on_OSP(ui *SC, ui &u1, ui &u2, vector<char> &HG_mark, ui S_end, ui C_end, ui *pstart, ui *pend, ui *edges, int tau, ui n)
{

    MaximumFlow mf;
    mf.reserve(C_end, C_end * C_end);
    nV = new ui[C_end + 3];
    oID = new ui[C_end + 3];

    int new_size = 0;
    for (ui i = 0; i < C_end; i++)
        nV[i] = new_size, oID[new_size++] = SC[i];
    mf.init(new_size + new_size);
    for (ui i = 0; i < C_end; i++)
        mf.addedge(nV[i], nV[i] + new_size, 1);
    for (ui i = 0; i < C_end; i++)
    {
        ui u = SC[i];
        for (ui j = i + 1; j < C_end; j++)
        {
            ui v = SC[j];
            bool is_adj = false;
            for (ui k = pstart[u]; k < pend[u]; k++)
            {
                if (edges[k] == v)
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
                continue;
            mf.addedge(nV[i] + new_size, nV[j], INF);
            mf.addedge(nV[j] + new_size, nV[i], INF);
        }
    }
    assert(new_size == C_end);
    for (ui i = 0; i < S_end; ++i)
        for (ui j = i + 1; j < S_end; ++j)
        {

            bool is_adj = false;
            for (ui k = pstart[oID[i]]; k < pend[oID[i]]; k++)
                if (edges[k] == oID[j])
                {
                    is_adj = true;
                    break;
                }
            if (!is_adj)
            {
                vector<ui> real_minCut;
                ui temp_vc = mf.SAP(new_size + i, j, new_size, real_minCut);
                if (temp_vc < C_end - tau)
                {
                    assert(!is_adj);
                    u1 = oID[i];
                    u2 = oID[j];
                    vector<char> value(n, 0);
                    for (ui k = 0; k < real_minCut.size(); k++)
                    {
                        assert(oID[real_minCut[k]] < n);
                        value[oID[real_minCut[k]]] = 1;
                    }
                    HG_mark = value;
                    assert(!HG_mark.empty());
                    value.clear();
                    return;
                }
            }
        }

    for (ui i = 0; i < S_end; ++i)
    {
        for (ui j = S_end; j < C_end; ++j)
        {

            bool is_adj = false;
            for (ui k = pstart[oID[i]]; k < pend[oID[i]]; k++)
                if (edges[k] == oID[j])
                {
                    is_adj = true;
                    break;
                }
            if (!is_adj)
            {
                vector<ui> real_minCut;
                ui temp_vc = mf.SAP(new_size + i, j, new_size, real_minCut);
                if (temp_vc < C_end - tau)
                {
                    assert(!is_adj);
                    u1 = oID[i];
                    u2 = oID[j];
                    vector<char> value(n, 0);
                    for (ui k = 0; k < real_minCut.size(); k++)
                    {
                        assert(oID[real_minCut[k]] < n);
                        value[oID[real_minCut[k]]] = 1;
                    }
                    HG_mark = value;
                    assert(!HG_mark.empty());
                    value.clear();
                    return;
                }
            }
        }
    }

    for (ui i = S_end; i < C_end; ++i)
    {
        for (ui j = i + 1; j < C_end; ++j)
        {
            bool is_adj = false;
            for (ui k = pstart[oID[i]]; k < pend[oID[i]]; k++)
            {
                if (edges[k] == oID[j])
                {
                    is_adj = true;
                    break;
                }
            }
            if (!is_adj)
            {
                vector<ui> real_minCut;
                ui temp_vc = mf.SAP(new_size + i, j, new_size, real_minCut);
                if (temp_vc < C_end - tau)
                {
                    assert(!is_adj);
                    u1 = oID[i];
                    u2 = oID[j];
                    vector<char> value(n, 0);
                    for (ui k = 0; k < real_minCut.size(); k++)
                    {
                        assert(oID[real_minCut[k]] < n);
                        value[oID[real_minCut[k]]] = 1;
                    }
                    HG_mark = value;
                    assert(!HG_mark.empty());
                    value.clear();
                    return;
                }
            }
        }
    }
}
