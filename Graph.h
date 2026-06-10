#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"

class Graph 
{
private:
	std::string dir; // input graph directory
	ui n;			 // number of nodes of the graph
	ept m;			 // number of edges of the graph
	ui tau;			 // the value of s in s-bundle

	ept *pstart; // offset of neighbors of nodes
	ept *pend;	 // used in search
	ept *pend_buf;
	ui *edges; // adjacent ids of edges
	ui *edgelist_pointer;

	std::vector<ui> candidate_community;

public:
	Graph(const char *_dir);
	~Graph();

	void readBinFile();													
	std::string get_file_name_without_suffix(const std::string &file_path);
	void SVC_TD(ui query_vertex, const ui LB_size, const ui UB_size);

private:
	void reorganize_adjacency_lists(ui n, ui *peel_sequence, ui *rid, ui *pstart, ui *pend, ui *edges);
	void extract_subgraph_with_prune(ui u, ui degree_threshold, ui triangle_threshold, ui cn_threshold, const ui *p_rid, ui *degree, std::vector<ui> &ids, ui *rid, std::vector<std::pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges);
	void extract_subgraph_wo_prune(ui u, const ui *p_rid, std::vector<ui> &ids, ui *rid, std::vector<std::pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges);

	void write_subgraph(ui n, const std::vector<std::pair<ui, ui>> &edge_list);
	void extract_subgraph(ui u, ui *ids, ui &ids_n, ui *rid, std::vector<std::pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer);
	void extract_subgraph_and_prune(ui u, ui *ids, ui &ids_n, ui *rid, std::vector<std::pair<ui, ui>> &vp, ui *Q, ui *degree, char *exists, ept *pend, char *deleted, ui *edgelist_pointer);
	void extract_subgraph_full(const ui *ids, ui ids_n, ui *rid, std::vector<std::pair<ui, ui>> &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer);

	ui degen(ui n, ui *peel_sequence, ui *core, ept *pstart, ui *edges, ui *degree, char *vis, ListLinearHeap *heap, bool output, ui &core_q, ui query_vertex);
	void core_shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *&pstart, ui *&edges, bool output);
	void orient_graph(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *rid);
	void oriented_triangle_counting(ui n, ui m, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj);
	void reorganize_oriented_graph(ui n, ui *tri_cnt, ui *edge_list, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edgelist_pointer, ui *buf);
	ept peeling(ListLinearHeap *linear_heap, ui *Qv, ui &Qv_n, ui d_threshold, ui *Qe, bool initialize_Qe, ui t_threshold, ui *tri_cnt, ui *active_edgelist, ui &active_edgelist_n, ui *edge_list, ui *edgelist_pointer, char *deleted, ui *degree, ept *pstart, ept *pend, ui *edges, char *exists);
	char find(ui u, ui w, ept &b, ept e, char *deleted, ept &idx, ui *edgelist_pointer, ui *edges);

	void update_query_vertex(ui n, ui *out_mapping, ui &query_vertex, const ui original_q_ID);
	void extract_subgraph_within_2hops(ui u, std::vector<ui> &ids, ui *rid, std::vector<std::pair<ui, ui>> &vp, char *exists, ept *pstart, ui *edges);
	void update_candidate_community(ui original_vertex, ui *out_mapping, std::vector<ui> &record_solution, bool &find_feasible_solution);
	void SVC_heu(std::vector<ui> &heu_solution, ui &vc_min, ui vc_max, int QID, int LB_size, int UB_size);

};
#endif
