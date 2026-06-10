#include "Graph.h"
#include "Sym_SE_BB_matrix.h"
#include "Sym_SE_BB.h"
#include "CTPrune.h"
#include <fstream>

using namespace std;
ui pre_n;
ui pre_m; 

Graph::Graph(const char *_dir)
{
	dir = string(_dir);

	n = m = 0;

	pstart = nullptr;
	pend = pend_buf = nullptr;
	edges = nullptr;
	edgelist_pointer = nullptr;

	candidate_community.clear();
}

Graph::~Graph()
{
	if (pstart != nullptr)
	{
		delete[] pstart;
		pstart = nullptr;
	}
	if (pend != nullptr)
	{
		delete[] pend;
		pend = nullptr;
	}
	if (pend_buf != nullptr)
	{
		delete[] pend_buf;
		pend_buf = nullptr;
	}
	if (edges != nullptr)
	{
		delete[] edges;
		edges = nullptr;
	}
	if (edgelist_pointer != nullptr)
	{
		delete[] edgelist_pointer;
		edgelist_pointer = nullptr;
	}
}

void Graph::readBinFile()
{
	ifstream in(dir);
	if (!in.is_open())
	{
		printf("Failed to open %s \n", dir.c_str());
		fflush(stdout);
		exit(1);
	}
	string suffix = "bin";
	if (suffix == "bin")
	{
		FILE *in = fopen(dir.c_str(), "rb");
		if (in == nullptr)
		{
			printf("Failed to open %s \n", dir.c_str());
			exit(1);
		}
		ui size_int;
		fread(&size_int, sizeof(ui), 1, in);
		if (size_int != sizeof(ui))
		{
			printf("sizeof ui is different: graph_file(%u), machine(%u)\n", size_int, (ui)sizeof(ui));
			exit(1);
		}
		fread(&n, sizeof(ui), 1, in);
		fread(&m, sizeof(ui), 1, in);
		cout << "File: " << get_file_name_without_suffix(dir) << " n= " << n << " m= " << m / 2 << endl;
		ui *degree = new ui[n];
		if (pstart == nullptr)
			pstart = new ept[n + 1];
		if (edges == nullptr)
			edges = new ui[m];
		fread(degree, sizeof(ui), n, in);
		fread(edges, sizeof(ui), m, in);
		pstart[0] = 0;
		for (ui i = 1; i <= n; i++)
			pstart[i] = pstart[i - 1] + degree[i - 1];
		delete[] degree;
	}

	printf("Graph init ok\n");
	fflush(stdout);
	in.close();
}

std::string Graph::get_file_name_without_suffix(const std::string &file_path)
{
	size_t last_slash = file_path.find_last_of("/\\");
	size_t last_dot = file_path.find_last_of(".");
	return file_path.substr(last_slash + 1, last_dot - last_slash - 1);
}

void Graph::SVC_TD(ui query_vertex, const ui LB_size, const ui UB_size)
{
	if (query_vertex > n - 1)
	{
		printf("The number of the query vertex exceeds the permitted maximum value %d !\n", n - 1);
		return;
	}

	Timer t;
	ui vc_min = 0, vc_max = UB_size - 1, core_q = n;

	char *vis = new char[n];
	ui *peel_sequence = new ui[n], *core = new ui[n], *degree = new ui[n];
	ListLinearHeap *heap = new ListLinearHeap(n, n - 1);
	ui UB = degen(n, peel_sequence, core, pstart, edges, degree, vis, heap, true, core_q, query_vertex);
	vc_max = (vc_max > core_q) ? core_q : vc_max;

	ui original_query_vertex = query_vertex;
	ui t_vc_min = 0;
	vector<ui> optimal_community;
	bool find_better_community = false, find_optimal_community = false;
	vector<ui> heu_candidate;

	SVC_heu(heu_candidate, t_vc_min, vc_max, query_vertex, LB_size, UB_size);
	vc_min = vc_min > t_vc_min ? vc_min : t_vc_min;

	bool heu_is_optimal = false;
	if (vc_min == vc_max)
	{
		find_optimal_community = true;
		heu_is_optimal = true;
	}

	if (!find_optimal_community)
	{
		tau = 1;
		candidate_community.resize(tau + vc_min);
		ui old_size = candidate_community.size(), *out_mapping = new ui[n], *rid = new ui[n], pre_n = n;
		core_shrink_graph(n, m, peel_sequence, core, out_mapping, nullptr, rid, pstart, edges, false);
		if (n != pre_n)
		{
			update_query_vertex(n, out_mapping, query_vertex, original_query_vertex);
			if (query_vertex == n)
			{
				printf("Sorry! But the query vertex is deleted from the graph! There is no better community!!\n");
				find_optimal_community = true;
			}
		}

		if (!find_optimal_community)
		{
			vector<ui> ids;
			vector<pair<ui, ui>> vp;
			memset(vis, 0, sizeof(char) * n);
			extract_subgraph_within_2hops(query_vertex, ids, rid, vp, vis, pstart, edges);

			tau = LB_size > vc_max ? LB_size - vc_max : 1;
			ui size_threshold = (tau > LB_size - 1) ? tau : LB_size - 1;
			size_threshold = (size_threshold > tau + vc_min) ? size_threshold : (tau + vc_min);
			candidate_community.resize(size_threshold);
			while (vc_min < vc_max)
			{
				old_size = candidate_community.size();
				UB = (UB_size > vc_max + tau) ? vc_max + tau : UB_size;
				assert(tau + vc_max <= UB_size && tau + vc_max >= LB_size);

				Sym_SE_BB_matrix *optimal_solver_m = new Sym_SE_BB_matrix();
				optimal_solver_m->allocateMemory(n);
				optimal_solver_m->load_graph(ids.size(), vp);
				optimal_solver_m->get_optimal_community(tau, UB, candidate_community, true);
				if (candidate_community.size() >= LB_size && candidate_community.size() > vc_min + tau)
				{
					if (candidate_community.size() > old_size)
					{
						old_size = candidate_community.size();
						for (ui j = 0; j < candidate_community.size(); j++)
							candidate_community[j] = ids[candidate_community[j]];
					}

					if (candidate_community.size() >= vc_max + tau)
					{

						vc_min = vc_max;
						update_candidate_community(original_query_vertex, out_mapping, optimal_community, find_better_community);
						find_optimal_community = true;
						break;
					}

					ui t_size = candidate_community.size() > UB_size ? UB_size : candidate_community.size();
					ui t_vc_min = t_size - tau, old_vc_min = vc_min;
					vc_min = (vc_min < t_vc_min) ? t_vc_min : vc_min;
					if (vc_min > old_vc_min)
						update_candidate_community(original_query_vertex, out_mapping, optimal_community, find_better_community);
				}
				delete optimal_solver_m;

				if (n > candidate_community.size() && candidate_community.size() < 2 * tau - 2 && 2 * tau - 2 >= LB_size)
				{
					UB = (2 * tau - 2 < UB) ? 2 * tau - 2 : UB;
					Sym_SE_BB *optimal_solver = new Sym_SE_BB();
					optimal_solver->allocateMemory(n);
					optimal_solver->load_graph(n, pstart, pstart + 1, edges);
					optimal_solver->get_optimal_community(tau, UB, candidate_community, false, query_vertex, UB);

					if (candidate_community.size() > 2 * tau - 2)
						printf("!!! WA in SVC_exact!\n");

					if (candidate_community.size() >= LB_size && candidate_community.size() > vc_min + tau)
					{

						if (candidate_community.size() >= vc_max + tau - 1)
						{
							if (candidate_community.size() >= vc_max + tau)
								vc_min = vc_max;
							else
								vc_min = vc_max = vc_max - 1;
							assert(candidate_community.size() >= vc_max + tau - 1);
							update_candidate_community(original_query_vertex, out_mapping, optimal_community, find_better_community);
							find_optimal_community = true;
							break;
						}

						ui t_size = candidate_community.size() > UB_size ? UB_size : candidate_community.size();
						ui t_vc_min = t_size - tau, old_vc_min = vc_min;
						vc_min = (vc_min < t_vc_min) ? t_vc_min : vc_min;
						if (vc_min > old_vc_min)
							update_candidate_community(original_query_vertex, out_mapping, optimal_community, find_better_community);
					}
					delete optimal_solver;
				}

				assert(candidate_community.size() < UB_size);

				if (tau + vc_max >= UB_size)
					vc_max--;
				tau++;
				assert(tau < UB_size);
				size_threshold = (tau > LB_size - 1) ? tau : LB_size - 1;
				size_threshold = (size_threshold > tau + vc_min) ? size_threshold : tau + vc_min;
				candidate_community.resize(size_threshold);
			}
		}

		delete[] out_mapping;
		delete[] rid;
	}
	delete heap;
	delete[] core;
	delete[] peel_sequence;
	delete[] vis;
	delete[] degree;

	printf("\tThe total running time is: %s (microseconds)\n", Utility::integer_to_string(t.elapsed()).c_str());
	printf("\tThe optimal vertex connectivity is: %d\n", vc_min);
	printf("The members of the found optimal community are shown as follow:\n");

	if (!find_optimal_community && !find_better_community)
		heu_is_optimal = true;

	if (!heu_is_optimal)
	{
		for (auto result : optimal_community)
			printf("%u, ", result);
		printf("\n");
	}
	else
	{
		for (auto result : heu_candidate)
			printf("%u, ", result);
		printf("\n");
	}
	fflush(stdout);
}

void Graph::SVC_heu(vector<ui> &heu_solution, ui &vc_min, ui vc_max, int QID, int LB_size, int UB_size)
{

	heu_solution.clear();
	MaximumFlow *vc = new MaximumFlow();

	vector<ui> t_solution;
	vector<ui> vis(n, 0);
	vector<ui> deg(n, 0);
	vector<ui> deg_in_S(n, 0);
	priority_queue<pair<ui, ui>> Q;

	t_solution.push_back(QID);
	vis[QID] = 1;
	for (ui i = pstart[QID]; i < pstart[QID + 1]; ++i)
	{
		++deg_in_S[edges[i]];
		if (!vis[edges[i]])
			Q.push(make_pair(deg_in_S[edges[i]], edges[i]));
	}

	while (t_solution.size() < UB_size)
	{

		ui v = n;
		while (true)
		{
			if (Q.empty())
				break;
			ui t_v = Q.top().second;
			Q.pop();
			assert(deg_in_S[t_v] > 0);
			if (!vis[t_v])
			{
				v = t_v;
				break;
			}
		}

		if (v == n)
			break;
		assert(v != n && !vis[v]);

		t_solution.push_back(v);
		vis[v] = 1;
		for (ui i = pstart[v]; i < pstart[v + 1]; i++)
		{
			ui t_v = edges[i];
			deg_in_S[t_v]++;
			if (!vis[t_v])
				Q.push(make_pair(deg_in_S[t_v], t_v));

			if (vis[t_v])
			{
				++deg[v];
				++deg[t_v];
			}
		}

		assert(t_solution.size() <= UB_size);
		if (t_solution.size() >= LB_size)
		{
			ui min_d = n;
			for (ui i = 0; i < t_solution.size(); ++i)
				if (deg[t_solution[i]] < min_d)
					min_d = deg[t_solution[i]];
			assert(min_d < n);

			if (min_d > vc_min)
			{

				ui t_VC = vc->get_exact_vertex_connectivity(t_solution, t_solution.size(), pstart, edges, vc_min);
				if (t_VC > vc_min)
				{
					vc_min = t_VC;
					heu_solution = t_solution;
					if (vc_min == vc_max)
						break;
				}
			}
		}
	}

	printf("The value of vc_min identified by our heuristic method is %d\n", vc_min);
	delete vc;
}

void Graph::update_query_vertex(ui n, ui *out_mapping, ui &query_vertex, const ui original_q_ID)
{
	ui new_query_vertex = n;
	for (ui i = 0; i < n; i++)
		if (out_mapping[i] == original_q_ID)
		{
			new_query_vertex = i;
			break;
		}
	query_vertex = new_query_vertex;
}

void Graph::update_candidate_community(ui original_query_vertex, ui *out_mapping, vector<ui> &optimal_community, bool &find_better_community)
{
	find_better_community = true;
	optimal_community.resize(candidate_community.size());
	if (out_mapping[candidate_community[0]] != original_query_vertex)
	{
		printf("当前的out_mapping[candidate_community[0]]的值是%u\n", out_mapping[candidate_community[0]]);
		printf("当前community中的成员分别是: \n");
		for (ui vertex : candidate_community)
			printf("%u, ", vertex);
		printf("\n");
	}
	assert(out_mapping[candidate_community[0]] == original_query_vertex);
	for (ui i = 0; i < candidate_community.size(); i++)
		optimal_community[i] = out_mapping[candidate_community[i]];
}

void Graph::extract_subgraph_within_2hops(ui u, vector<ui> &ids, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ui *edges)
{
	ids.clear();
	vp.clear();
	ids.push_back(u);
	exists[u] = 1;
	rid[u] = 0;

	for (ept i = pstart[u]; i < pstart[u + 1]; i++)
	{
		ids.push_back(edges[i]);
		exists[edges[i]] = 1;
		rid[edges[i]] = ids.size() - 1;
	}

	ui old_size = ids.size();
	for (ui i = 1; i < old_size; i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pstart[v + 1]; j++)
		{
			ui w = edges[j];
			if (exists[w])
				continue;
			ids.push_back(w);
			exists[w] = 1;
			rid[w] = ids.size() - 1;
		}
	}

	for (ui i = 0; i < ids.size(); i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pstart[v + 1]; j++)
			if (exists[edges[j]] && rid[v] < rid[edges[j]])
				vp.push_back(make_pair(rid[v], rid[edges[j]]));
	}

	for (ui i = 0; i < ids.size(); i++)
		exists[ids[i]] = 0;
}

void Graph::reorganize_adjacency_lists(ui n, ui *peel_sequence, ui *rid, ui *pstart, ui *pend, ui *edges)
{
	for (ui i = 0; i < n; i++)
		rid[peel_sequence[i]] = i;

	for (ui i = 0; i < n; i++)
	{
		ui &end = pend[i] = pstart[i];
		for (ui j = pstart[i]; j < pstart[i + 1]; j++)
			if (rid[edges[j]] > rid[i])
				edges[end++] = edges[j];
	}

	for (ui i = n; i > 0; i--)
	{
		ui u = peel_sequence[i - 1];
		for (ui j = pstart[u]; j < pend[u] && rid[edges[j]] > rid[u]; j++)
		{
			ui v = edges[j];
			edges[pend[v]++] = u;
			assert(pend[v] <= pstart[v + 1]);
		}
	}
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(pend[i] == pstart[i + 1]);
#endif

	for (ui i = 0; i < n; i++)
	{
		ui &end = pend[i] = pstart[i];
		while (end < pstart[i + 1] && rid[edges[end]] > rid[i])
			++end;
	}
}

void Graph::extract_subgraph_with_prune(ui u, ui degree_threshold, ui triangle_threshold, ui cn_threshold, const ui *p_rid, ui *degree, vector<ui> &ids, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges)
{
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(!exists[i]);
#endif

	ids.clear();
	vp.clear();
	ids.push_back(u);
	exists[u] = 1;
	for (ept i = pstart[u]; i < pend[u]; i++)
	{
		assert(p_rid[edges[i]] > p_rid[u]);
		ids.push_back(edges[i]);
		exists[edges[i]] = 2;
	}
	assert(pend[u] >= pstart[u + 1] || p_rid[edges[pend[u]]] < p_rid[u]);

	ui *Q = rid;
	ui Q_n = 0;
	for (ui i = 1; i < ids.size(); i++)
	{
		ui v = ids[i];
		degree[v] = 0;
		for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++)
		{
			if (exists[edges[j]])
				++degree[v];
		}
		if (degree[v] < triangle_threshold)
			Q[Q_n++] = v;
	}

	for (ui i = 0; i < Q_n; i++)
	{
		ui v = Q[i];
		exists[v] = 3;
		for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++)
			if (exists[edges[j]] == 2)
			{
				if (degree[edges[j]] == triangle_threshold)
					Q[Q_n++] = edges[j];
				--degree[edges[j]];
			}
	}
	assert(Q_n < ids.size());

	if (ids.size() - Q_n - 1 < degree_threshold)
	{
		for (ui i = 0; i < ids.size(); i++)
			exists[ids[i]] = 0;
		ids.clear();
		return;
	}

	ui old_size = ids.size();
	for (ui i = 1; i < old_size; i++)
		if (exists[ids[i]] == 2)
		{
			ui v = ids[i];
			for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++)
			{
				if (!exists[edges[j]])
				{
					ids.push_back(edges[j]);
					exists[edges[j]] = 1;
					degree[edges[j]] = 1;
				}
				else
					++degree[edges[j]];
			}
		}

	ui new_size = 1;
	for (ui i = 1; i < old_size; i++)
	{
		if (exists[ids[i]] == 3)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
	assert(new_size + Q_n == old_size);

	for (ui i = old_size; i < ids.size(); i++)
	{
		if (degree[ids[i]] < cn_threshold)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
	ids.resize(new_size);

	for (ui i = 0; i < ids.size(); i++)
		rid[ids[i]] = i;

	for (ui i = 0; i < ids.size(); i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pend[v]; j++)
			if (exists[edges[j]])
			{
				assert(rid[v] < ids.size() && rid[edges[j]] < ids.size());
				vp.push_back(make_pair(rid[v], rid[edges[j]]));
			}
	}

	for (ui i = 0; i < ids.size(); i++)
		exists[ids[i]] = 0;
}

void Graph::extract_subgraph_wo_prune(ui u, const ui *p_rid, vector<ui> &ids, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges)
{
	ids.clear();
	vp.clear();
	ids.push_back(u);
	exists[u] = 1;
	rid[u] = 0;

	for (ept i = pstart[u]; i < pend[u]; i++)
	{
		assert(p_rid[edges[i]] > p_rid[u]);
		ids.push_back(edges[i]);
		exists[edges[i]] = 1;
		rid[edges[i]] = ids.size() - 1;
	}
	assert(pend[u] >= pstart[u + 1] || p_rid[edges[pend[u]]] < p_rid[u]);
	ui old_size = ids.size();
	for (ui i = 1; i < old_size; i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++)
		{
			ui w = edges[j];
			if (exists[w])
				continue;
			ids.push_back(w);
			exists[w] = 1;
			rid[w] = ids.size() - 1;
		}
	}
	for (ui i = 0; i < ids.size(); i++)
	{
		ui v = ids[i];
		for (ept j = pstart[v]; j < pend[v]; j++)
			if (exists[edges[j]])
				vp.push_back(make_pair(rid[v], rid[edges[j]]));
	}
	for (ui i = 0; i < ids.size(); i++)
		exists[ids[i]] = 0;
}

void Graph::write_subgraph(ui n, const vector<pair<ui, ui>> &edge_list)
{
	FILE *fout = Utility::open_file("edges.txt", "w");

	fprintf(fout, "%u %lu\n", n, edge_list.size());
	for (ui i = 0; i < edge_list.size(); i++)
		fprintf(fout, "%d %d\n", edge_list[i].first, edge_list[i].second);

	fclose(fout);
}

void Graph::extract_subgraph(ui u, ui *ids, ui &ids_n, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer)
{
	ids_n = 0;
	vp.clear();
	ids[ids_n++] = u;
	exists[u] = 1;
	rid[u] = 0;
	ui u_n = pstart[u];
	for (ept i = pstart[u]; i < pend[u]; i++)
		if (!deleted[edgelist_pointer[i]])
		{
			edges[u_n] = edges[i];
			edgelist_pointer[u_n++] = edgelist_pointer[i];
			ui v = edges[i];
			rid[v] = ids_n;
			ids[ids_n++] = v;
			exists[v] = 1;
		}
	pend[u] = u_n;
	ui old_size = ids_n;
	for (ui i = 1; i < old_size; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				ui v = edges[j];
				if (exists[v])
					continue;
				rid[v] = ids_n;
				ids[ids_n++] = v;
				exists[v] = 1;
			}
		pend[u] = u_n;
	}
	for (ui i = 0; i < old_size; i++)
	{
		u = ids[i];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (edges[j] > u)
			{
				vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
	}
	for (ui i = old_size; i < ids_n; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (edges[j] > u && exists[edges[j]])
					vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
		pend[u] = u_n;
	}
	for (ui i = 0; i < ids_n; i++)
		exists[ids[i]] = 0;
}

void Graph::extract_subgraph_full(const ui *ids, ui ids_n, ui *rid, vector<pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer)
{
	vp.clear();
	for (ui i = 0; i < ids_n; i++)
	{
		rid[ids[i]] = i;
		exists[ids[i]] = 1;
	}
	for (ui i = 0; i < ids_n; i++)
	{
		ui u = ids[i];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]] && edges[j] > u && exists[edges[j]])
			{
				vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
	}
	for (ui i = 0; i < ids_n; i++)
		exists[ids[i]] = 0;
}

void Graph::extract_subgraph_and_prune(ui u, ui *ids, ui &ids_n, ui *rid, vector<pair<ui, ui>> &vp, ui *Q, ui *degree, char *exists, ept *pend, char *deleted, ui *edgelist_pointer)
{
	vp.clear();
	ids_n = 0;
	ids[ids_n++] = u;
	exists[u] = 1;
	ui u_n = pstart[u];
	for (ept i = pstart[u]; i < pend[u]; i++)
		if (!deleted[edgelist_pointer[i]])
		{
			edges[u_n] = edges[i];
			edgelist_pointer[u_n++] = edgelist_pointer[i];
			ui v = edges[i];
			ids[ids_n++] = v;
			exists[v] = 2;
		}
	pend[u] = u_n;

	ui Q_n = 0;
	for (ui i = 1; i < ids_n; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		degree[u] = 0;
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (exists[edges[j]] == 2)
					++degree[u];
			}
		pend[u] = u_n;
		if (degree[u] + 2 * tau <= candidate_community.size())
			Q[Q_n++] = u;
	}
	for (ui i = 0; i < Q_n; i++)
	{
		u = Q[i];
		exists[u] = 10;
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (exists[edges[j]] == 2)
			{
				if ((degree[edges[j]]--) + 2 * tau == candidate_community.size() + 1)
				{
					assert(Q_n < m / 2);
					Q[Q_n++] = edges[j];
				}
			}
	}
	assert(Q_n <= ids_n);
	if (ids_n - 1 - Q_n + tau <= candidate_community.size())
	{
		for (ui i = 0; i < ids_n; i++)
			exists[ids[i]] = 0;
		ids_n = 0;
		return;
	}

	ui nr_size = ids_n;
	for (ui i = 1; i < nr_size; i++)
		if (exists[ids[i]] == 2)
		{
			u = ids[i];
			for (ept j = pstart[u]; j < pend[u]; j++)
			{
				if (!exists[edges[j]])
				{
					ids[ids_n++] = edges[j];
					exists[edges[j]] = 3;
					degree[edges[j]] = 1;
				}
				else if (exists[edges[j]] == 3)
					++degree[edges[j]];
			}
		}

#ifndef NDEBUG

#endif

	ui new_size = 1;
	for (ui i = 1; i < nr_size; i++)
	{
		if (exists[ids[i]] == 10)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
#ifndef NDEBUG
	if (new_size + Q_n != nr_size)
	{
		printf("new_size: %u, Q_n: %u, nr_size: %u\n", new_size, Q_n, nr_size);
		printf("New list: ");
		for (ui i = 0; i < new_size; i++)
			printf(" %u", ids[i]);
		printf("\n");
		printf("Pruned list: ");
		for (ui i = 0; i < Q_n; i++)
			printf(" %u", Q[i]);
		printf("\n");
	}
#endif
	assert(new_size + Q_n == nr_size);
	ui old_nr_size = nr_size;
	nr_size = new_size;
	for (ui i = old_nr_size; i < ids_n; i++)
	{
		if (degree[ids[i]] + 2 * tau <= candidate_community.size() + 2)
			exists[ids[i]] = 0;
		else
			ids[new_size++] = ids[i];
	}
	ids_n = new_size;
#ifndef NDEBUG
	assert(exists[ids[0]] == 1);
	for (ui i = 1; i < nr_size; i++)
		assert(exists[ids[i]] == 2);
	for (ui i = nr_size; i < ids_n; i++)
		assert(exists[ids[i]] == 3);
#endif

	for (ui i = 0; i < ids_n; i++)
	{
		assert(exists[ids[i]]);
		rid[ids[i]] = i;
	}

	for (ui i = 0; i < nr_size; i++)
	{
		u = ids[i];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (exists[edges[j]] && edges[j] > u)
			{
				assert(!deleted[edgelist_pointer[j]]);
				vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
	}
	for (ui i = nr_size; i < ids_n; i++)
	{
		u = ids[i];
		u_n = pstart[u];
		for (ept j = pstart[u]; j < pend[u]; j++)
			if (!deleted[edgelist_pointer[j]])
			{
				edges[u_n] = edges[j];
				edgelist_pointer[u_n++] = edgelist_pointer[j];
				if (edges[j] > u && exists[edges[j]])
					vp.push_back(make_pair(rid[u], rid[edges[j]]));
			}
		pend[u] = u_n;
	}
	for (ui i = 0; i < ids_n; i++)
		exists[ids[i]] = 0;
#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(exists[i] == 0);
#endif
}

ui Graph::degen(ui n, ui *peel_sequence, ui *core, ept *pstart, ui *edges, ui *degree, char *vis, ListLinearHeap *heap, bool output, ui &core_q, ui query_vertex)
{
	Timer t;

	ui threshold = (candidate_community.size() + 1 > tau ? candidate_community.size() + 1 - tau : 0);

	for (ui i = 0; i < n; i++)
		degree[i] = pstart[i + 1] - pstart[i];

	ui queue_n = 0, new_size = 0;
	for (ui i = 0; i < n; i++)
		if (degree[i] < threshold)
			peel_sequence[queue_n++] = i;

	for (ui i = 0; i < queue_n; i++)
	{
		ui u = peel_sequence[i];
		degree[u] = 0;
		for (ept j = pstart[u]; j < pstart[u + 1]; j++)
			if (degree[edges[j]] > 0)
				if ((degree[edges[j]]--) == threshold)
					peel_sequence[queue_n++] = edges[j];
	}
	ui UB = n;
	if (queue_n == n)
		UB = candidate_community.size();

	memset(vis, 0, sizeof(char) * n);
	for (ui i = 0; i < n; i++)
	{
		if (degree[i] >= threshold)
			peel_sequence[queue_n + (new_size++)] = i;
		else
		{
			vis[i] = 1;
			core[i] = 0;
			if (i == query_vertex)
				core_q = 0;
		}
	}
	assert(queue_n + new_size == n);

	if (new_size != 0)
	{
		heap->init(new_size, new_size - 1, peel_sequence + queue_n, degree);
		ui max_core = 0;
		ui idx = n;
		UB = 0;
		for (ui i = 0; i < new_size; i++)
		{
			ui u, key;
			heap->pop_min(u, key);
			if (key > max_core)
				max_core = key;
			core[u] = max_core;

			if (u == query_vertex)
				core_q = core[u];

			peel_sequence[queue_n + i] = u;

			ui t_UB = core[u] + tau;
			if (new_size - i < t_UB)
				t_UB = new_size - i;
			if (t_UB > UB)
				UB = t_UB;

			if (idx == n && key + tau >= new_size - i)
				idx = i;
			vis[u] = 1;

			for (ept j = pstart[u]; j < pstart[u + 1]; j++)
				if (vis[edges[j]] == 0)
					heap->decrement(edges[j], 1);
		}
	}

	return UB;
}

void Graph::core_shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *&pstart, ui *&edges, bool output)
{

	ui cnt = 0;
	for (ui i = 0; i < n; i++)
		if (core[i] + tau > candidate_community.size())
		{
			rid[i] = cnt;
			if (in_mapping == nullptr)
				out_mapping[cnt] = i;
			else
				out_mapping[cnt] = in_mapping[i];
			++cnt;
		}

	if (cnt != n)
	{
		cnt = 0;
		ept pos = 0;
		for (ui i = 0; i < n; i++)
			if (core[i] + tau > candidate_community.size())
			{
				ept t_start = pstart[i];
				pstart[cnt] = pos;
				for (ept j = t_start; j < pstart[i + 1]; j++)
					if (core[edges[j]] + tau > candidate_community.size())
					{
						edges[pos++] = rid[edges[j]];
					}
				++cnt;
			}
		pstart[cnt] = pos;

		assert(core[peel_sequence[n - cnt - 1]] == 0 || core[peel_sequence[n - cnt - 1]] + tau <= candidate_community.size());
		assert(cnt == 0 || core[peel_sequence[n - cnt]] + tau > candidate_community.size());

		for (ui i = 0; i < cnt; i++)
		{
			peel_sequence[i] = rid[peel_sequence[n - cnt + i]];
			core[i] = core[out_mapping[i]];
		}

		if (pos > 0 && pos < m / 2)
		{

			ept *pstart_new = new ept[cnt + 1];
			ui *edges_new = new ui[pos];
			memcpy(pstart_new, pstart, sizeof(ept) * (cnt + 1));
			memcpy(edges_new, edges, sizeof(ui) * pos);

			delete[] pstart;
			pstart = pstart_new;
			delete[] edges;
			edges = edges_new;
		}
		n = cnt;
		m = pos;
	}

	if (output)
		printf("*** After core shrink: n = %s, m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m / 2).c_str());
}

void Graph::orient_graph(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *rid)
{
	for (ui i = 0; i < n; i++)
		rid[peel_sequence[i]] = i;
	for (ui i = 0; i < n; i++)
	{
		ept &end = pend[i] = pstart[i];
		for (ept j = pstart[i]; j < pstart[i + 1]; j++)
			if (rid[edges[j]] > rid[i])
				edges[end++] = edges[j];
	}

#ifndef NDEBUG
	long long sum = 0;
	for (ui i = 0; i < n; i++)
		sum += pend[i] - pstart[i];
	assert(sum * 2 == m);
#endif
}

void Graph::oriented_triangle_counting(ui n, ui m, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj)
{
	memset(adj, 0, sizeof(ui) * n);
	long long cnt = 0;
	memset(tri_cnt, 0, sizeof(ui) * m);
	for (ui u = 0; u < n; u++)
	{
		for (ept j = pstart[u]; j < pend[u]; j++)
			adj[edges[j]] = j + 1;

		for (ept j = pstart[u]; j < pend[u]; j++)
		{
			ui v = edges[j];
			for (ept k = pstart[v]; k < pend[v]; k++)
				if (adj[edges[k]])
				{
					++tri_cnt[j];
					++tri_cnt[k];
					++tri_cnt[adj[edges[k]] - 1];
					++cnt;
				}
		}

		for (ept j = pstart[u]; j < pend[u]; j++)
			adj[edges[j]] = 0;
	}

#ifndef NDEBUG

#endif
}

void Graph::reorganize_oriented_graph(ui n, ui *tri_cnt, ui *edge_list, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edgelist_pointer, ui *buf)
{
	for (ui i = 0; i < n; i++)
		pend2[i] = pend[i];
	ept pos = 0;
	for (ui i = 0; i < n; i++)
	{
		for (ept j = pstart[i]; j < pend[i]; j++)
		{
			tri_cnt[pos >> 1] = edgelist_pointer[j];
			edge_list[pos++] = i;
			edge_list[pos++] = edges[j];

			ept &k = pend2[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j] = (pos >> 1) - 1;
			edges[k++] = i;
		}
	}

#ifndef NDEBUG
	for (ui i = 0; i < n; i++)
		assert(pend2[i] == pstart[i + 1]);
#endif

	for (ui i = 0; i < n; i++)
	{
		pend2[i] = pend[i];
		pend[i] = pstart[i];
	}
	for (ui i = 0; i < n; i++)
	{
		for (ept j = pend2[i]; j < pstart[i + 1]; j++)
		{
			ept &k = pend[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j];
			edges[k++] = i;
		}
	}

	ept *ids = pend2;
	for (ui i = 0; i < n; i++)
	{
		if (pend[i] == pstart[i] || pend[i] == pstart[i + 1])
			continue;
		ept j = pstart[i], k = pend[i], pos = 0;
		while (j < pend[i] && k < pstart[i + 1])
		{
			if (edges[j] < edges[k])
			{
				ids[pos] = edges[j];
				buf[pos++] = edgelist_pointer[j++];
			}
			else
			{
				ids[pos] = edges[k];
				buf[pos++] = edgelist_pointer[k++];
			}
		}
		while (j < pend[i])
		{
			ids[pos] = edges[j];
			buf[pos++] = edgelist_pointer[j++];
		}
		while (k < pstart[i + 1])
		{
			ids[pos] = edges[k];
			buf[pos++] = edgelist_pointer[k++];
		}
		for (ept j = 0; j < pos; j++)
		{
			edges[pstart[i] + j] = ids[j];
			edgelist_pointer[pstart[i] + j] = buf[j];
		}
	}
}

/*
		该函数接受以下参数:
		u: 顶点 u
		w: 顶点 w
		b: 边列表的起始索引
		e: 边列表的结束索引
		deleted: 一个标记被删除边的数组
		idx: 返回找到的边在边列表中的索引
		edgelist_pointer: 边列表的指针
		edges: 边列表本身
*/
char Graph::find(ui u, ui w, ept &b, ept e, char *deleted, ept &idx, ui *edgelist_pointer, ui *edges)
{
	if (b >= e)
		return 0;

	while (b + 1 < e)
	{
		idx = b + (e - b) / 2;
		if (edges[idx] > w)
			e = idx;
		else
			b = idx;
	}

	if (edges[b] == w)
	{
		idx = edgelist_pointer[b];
		if (!deleted[idx])
			return 1;
	}

	return 0;
}

ept Graph::peeling(ListLinearHeap *linear_heap, ui *Qv, ui &Qv_n, ui d_threshold, ui *Qe, bool initialize_Qe, ui t_threshold, ui *tri_cnt, ui *active_edgelist, ui &active_edgelist_n, ui *edge_list, ui *edgelist_pointer, char *deleted, ui *degree, ept *pstart, ept *pend, ui *edges, char *exists)
{
	ept Qe_n = 0;
#ifndef NO_TRUSS_PRUNE
	if (initialize_Qe)
	{
		ept active_edgelist_newn = 0;
		for (ept j = 0; j < active_edgelist_n; j++)
			if (!deleted[active_edgelist[j]])
			{
				if (tri_cnt[active_edgelist[j]] < t_threshold)
					Qe[Qe_n++] = active_edgelist[j];
				else
					active_edgelist[active_edgelist_newn++] = active_edgelist[j];
			}
		active_edgelist_n = active_edgelist_newn;
	}
#endif

	ept deleted_edges_n = 0;
	ui Qv_idx = 0;
	while (Qv_idx < Qv_n || Qe_n)
	{
		if (Qe_n == 0)
		{

			ui u = Qv[Qv_idx++];
			ept u_n = pstart[u];
			for (ept k = pstart[u]; k < pend[u]; k++)
				if (!deleted[edgelist_pointer[k]])
				{
					edges[u_n] = edges[k];
					edgelist_pointer[u_n++] = edgelist_pointer[k];
					exists[edges[k]] = 1;
				}
			pend[u] = u_n;

			for (ept k = pstart[u]; k < pend[u]; k++)
				deleted[edgelist_pointer[k]] = 1;
			deleted_edges_n += pend[u] - pstart[u];
			linear_heap->decrement(u, degree[u]);
			degree[u] = 0;

			for (ept k = pstart[u]; k < pend[u]; k++)
			{
				ui v = edges[k];
#ifndef NO_TRUSS_PRUNE
				ept v_n = pstart[v];
				for (ept x = pstart[v]; x < pend[v]; x++)
					if (!deleted[edgelist_pointer[x]])
					{
						edges[v_n] = edges[x];
						edgelist_pointer[v_n++] = edgelist_pointer[x];
						if (edges[x] > v && exists[edges[x]])
						{
							if ((tri_cnt[edgelist_pointer[x]]--) == t_threshold)
								Qe[Qe_n++] = edgelist_pointer[x];
						}
					}
				pend[v] = v_n;
#endif

				linear_heap->decrement(v, 1);
				if ((degree[v]--) == d_threshold)
					Qv[Qv_n++] = v;
			}

			for (ept k = pstart[u]; k < pend[u]; k++)
				exists[edges[k]] = 0;
		}
#ifdef NO_TRUSS_PRUNE
		Qe_n = 0;
#endif
		for (ept j = 0; j < Qe_n; j++)
		{
			ept idx = Qe[j];
			ui u = edge_list[idx << 1], v = edge_list[(idx << 1) + 1];
			ui tri_n = tri_cnt[idx];

			deleted[idx] = 1;
			linear_heap->decrement(u, 1);
			linear_heap->decrement(v, 1);
			if ((degree[u]--) == d_threshold)
				Qv[Qv_n++] = u;
			if ((degree[v]--) == d_threshold)
				Qv[Qv_n++] = v;
			deleted_edges_n++;

			if (degree[u] < degree[v])
				swap(u, v);

			if (degree[u] > degree[v] * 2)
			{

				ept v_n = pstart[v], start = pstart[u];
				for (ept k = pstart[v]; k < pend[v]; k++)
					if (!deleted[edgelist_pointer[k]])
					{
						edges[v_n] = edges[k];
						edgelist_pointer[v_n++] = edgelist_pointer[k];

						if (tri_n && find(u, edges[k], start, pend[u], deleted, idx, edgelist_pointer, edges))
						{
							--tri_n;
							if ((tri_cnt[idx]--) == t_threshold)
								Qe[Qe_n++] = idx;
							if ((tri_cnt[edgelist_pointer[k]]--) == t_threshold)
								Qe[Qe_n++] = edgelist_pointer[k];
						}
					}
				pend[v] = v_n;
				assert(tri_n == 0);
			}
			else
			{
				ept ii = pstart[u], jj = pstart[v];
				ept u_n = pstart[u], v_n = pstart[v];

				while (true)
				{
					while (ii < pend[u] && deleted[edgelist_pointer[ii]])
						++ii;
					while (jj < pend[v] && deleted[edgelist_pointer[jj]])
						++jj;
					if (ii >= pend[u] || jj >= pend[v])
						break;

					if (edges[ii] == edges[jj])
					{
						edges[u_n] = edges[ii];
						edgelist_pointer[u_n++] = edgelist_pointer[ii];
						edges[v_n] = edges[jj];
						edgelist_pointer[v_n++] = edgelist_pointer[jj];

						if ((tri_cnt[edgelist_pointer[ii]]--) == t_threshold)
							Qe[Qe_n++] = edgelist_pointer[ii];
						if ((tri_cnt[edgelist_pointer[jj]]--) == t_threshold)
							Qe[Qe_n++] = edgelist_pointer[jj];

						++ii;
						++jj;
					}
					else if (edges[ii] < edges[jj])
					{
						edges[u_n] = edges[ii];
						edgelist_pointer[u_n++] = edgelist_pointer[ii];
						++ii;
					}
					else
					{
						edges[v_n] = edges[jj];
						edgelist_pointer[v_n++] = edgelist_pointer[jj];
						++jj;
					}
				}
				while (ii < pend[u])
				{
					if (!deleted[edgelist_pointer[ii]])
					{
						edges[u_n] = edges[ii];
						edgelist_pointer[u_n++] = edgelist_pointer[ii];
					}
					++ii;
				}
				while (jj < pend[v])
				{
					if (!deleted[edgelist_pointer[jj]])
					{
						edges[v_n] = edges[jj];
						edgelist_pointer[v_n++] = edgelist_pointer[jj];
					}
					++jj;
				}
				pend[u] = u_n;
				pend[v] = v_n;
			}
		}
		Qe_n = 0;
	}
#ifndef NDEBUG
	printf("*** Truss removed %s undirected edges\n", Utility::integer_to_string(deleted_edges_n).c_str());
#endif
	return deleted_edges_n;
}
