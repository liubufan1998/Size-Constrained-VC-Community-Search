#ifndef _SYM_SE_BB_MATRIX_
#define _SYM_SE_BB_MATRIX_

#include "Utility.h"
#include "Timer.h"
#include "VC_tool.h"

class Sym_SE_BB_matrix
{
private:
	long long n;
	char *matrix;
	long long matrix_size; 
	vector<ui> temp_deleted;
	MaximumFlow *vc;
	long long branch_count;

	ui *degree;
	ui *degree_in_S;

	ui tau;
	ui *best_solution;
	ui best_solution_size;
	ui _UB_;
	ui *neighbors;
	ui *nonneighbors;
	ui *SC;
	ui *SC_rid;
	queue<ui> Qv;
	ui *level_id;
	vector<pair<ui, ui>> vp;

public:
	Sym_SE_BB_matrix()
	{
		n = 0;
		matrix = nullptr;
		matrix_size = 0;

		degree = degree_in_S = nullptr;

		best_solution = nullptr;
		tau = best_solution_size = _UB_ = 0;

		neighbors = nonneighbors = nullptr;

		SC = SC_rid = nullptr;
		level_id = nullptr;
	}

	~Sym_SE_BB_matrix()
	{
		if (matrix != NULL)
		{
			delete[] matrix;
			matrix = NULL;
		}
		if (degree != NULL)
		{
			delete[] degree;
			degree = NULL;
		}
		if (degree_in_S != NULL)
		{
			delete[] degree_in_S;
			degree_in_S = NULL;
		}
		if (best_solution != NULL)
		{
			delete[] best_solution;
			best_solution = NULL;
		}
		if (SC != NULL)
		{
			delete[] SC;
			SC = NULL;
		}
		if (SC_rid != NULL)
		{
			delete[] SC_rid;
			SC_rid = NULL;
		}
		if (neighbors != NULL)
		{
			delete[] neighbors;
			neighbors = NULL;
		}
		if (nonneighbors != NULL)
		{
			delete[] nonneighbors;
			nonneighbors = NULL;
		}
		if (level_id != NULL)
		{
			delete[] level_id;
			level_id = NULL;
		}
	}

	void allocateMemory(ui n)
	{
		if (n <= 0)
			return;

		matrix_size = 1;
		matrix = new char[matrix_size];
		degree = new ui[n];
		degree_in_S = new ui[n];
		best_solution = new ui[n];
		SC = new ui[n];
		SC_rid = new ui[n];
		neighbors = new ui[n];
		nonneighbors = new ui[n];
		level_id = new ui[n];
	}

	void load_graph(ui _n, const vector<pair<ui, ui>> &vp)
	{
		n = _n;
		if (((long long)n) * n > matrix_size)
		{
			do
			{
				matrix_size *= 2;
			} while (((long long)n) * n > matrix_size);
			delete[] matrix;
			matrix = new char[matrix_size];
		}
		memset(matrix, 0, sizeof(char) * ((long long)n) * n);
		for (ui i = 0; i < n; i++)
			degree[i] = 0;
		for (ui i = 0; i < vp.size(); i++)
		{
			assert(vp[i].first >= 0 && vp[i].first < n && vp[i].second >= 0 && vp[i].second < n);
			ui a = vp[i].first, b = vp[i].second;
			degree[a]++;
			degree[b]++;
			if (matrix[a * n + b])
				printf("Duplicate edge in KPLEX_BB_matrix.load_graph()\n");

			matrix[a * n + b] = matrix[b * n + a] = 1;
		}
	}

	void get_optimal_community(ui tau_, ui UB_, vector<ui> &candidate_community, bool must_include_0)
	{
		tau = tau_;
		_UB_ = UB_;
		best_solution_size = candidate_community.size();
		ui C_end;
		initialization(C_end, must_include_0, true);
		vc = new MaximumFlow();

		if (C_end)
			branching_with_fast_VCL(0, C_end, 1, must_include_0);

		if (vc != nullptr)
		{
			delete vc;
			vc = nullptr;
		}
		if (best_solution_size > _UB_)
			best_solution_size = _UB_;
		if (best_solution_size > candidate_community.size())
		{
			candidate_community.clear();
			for (int i = 0; i < best_solution_size; i++)
				candidate_community.push_back(best_solution[i]);
		}
	}

	void initialization(ui &C_end, bool must_include_0, bool is_precise)
	{

		ui *peel_sequence = neighbors;
		ui *core = nonneighbors;
		ui *vis = SC;
		memset(vis, 0, sizeof(ui) * n);
		ui max_core = 0, idx = n;
		ui idx_splex = n;
		for (ui i = 0; i < n; i++)
		{
			ui u, min_degree = n;
			for (ui j = 0; j < n; j++)
				if (!vis[j] && degree[j] < min_degree)
				{
					u = j;
					min_degree = degree[j];
				}
			if (min_degree > max_core)
				max_core = min_degree;
			core[u] = max_core;
			peel_sequence[i] = u;
			vis[u] = 1;

			for (ui j = 0; j < n; j++)
				if (!vis[j] && matrix[u * n + j])
					--degree[j];
		}

		C_end = 0;
		for (ui i = 0; i < n; i++)
			SC_rid[i] = n;
		for (ui i = 0; i < n; i++)
			if (core[i] + tau > best_solution_size)
			{
				SC[C_end] = i;
				SC_rid[i] = C_end;
				++C_end;
			}

		if ((must_include_0 && SC_rid[0] == n) || best_solution_size >= _UB_)
		{
			C_end = 0;
			return;
		}

		for (ui i = 0; i < C_end; i++)
		{
			ui u = SC[i];
			degree[u] = degree_in_S[u] = 0;
			for (ui j = 0; j < C_end; j++)
				if (matrix[u * n + SC[j]])
					++degree[u];
		}

		memset(level_id, 0, sizeof(ui) * n);
		for (ui i = 0; i < C_end; i++)
			level_id[SC[i]] = n;

		if (!Qv.empty())
			printf("!!! Something wrong. Qv must be empty in initialization\n");
	}

	void store_better_community(ui size)
	{
		if (size <= best_solution_size)
		{
			printf("!!! the solution to store is no larger than the current best solution!");
			return;
		}

		best_solution_size = size;
		for (ui i = 0; i < best_solution_size; i++)
			best_solution[i] = SC[i];
	}

	bool satisfy_degree_requirement(ui C_end)
	{
		for (ui i = 0; i < C_end; i++)
			if (degree[SC[i]] + tau < C_end)
				return false;
		return true;
	}

	void branching_with_fast_VCL(ui S_end, ui C_end, ui level, bool choose_zero)
	{

		if (S_end > best_solution_size && vc->satisfy_necessary_requirement(SC, S_end, matrix, n, tau, true))
			store_better_community(S_end);
		if (C_end > best_solution_size && satisfy_degree_requirement(C_end) && vc->satisfy_necessary_requirement(SC, C_end, matrix, n, tau, true))
			store_better_community(C_end);
		if (C_end <= best_solution_size + 1 || best_solution_size >= _UB_)
			return;

#ifndef NDEBUG
		for (ui i = 0; i < C_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SC[i] * n + SC[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < C_end; j++)
				if (matrix[SC[i] * n + SC[j]])
					++d2;
			assert(d1 == degree_in_S[SC[i]]);
			assert(d2 == degree[SC[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SC[i]] + tau >= S_end);
		for (ui i = S_end; i < C_end; i++)
			assert(degree_in_S[SC[i]] + tau > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SC[i]] + tau == S_end)
			{
				char *t_matrix = matrix + SC[i] * n;
				for (ui j = S_end; j < C_end; j++)
					assert(t_matrix[SC[j]]);
			}
		for (ui i = 0; i < C_end; i++)
			assert(level_id[SC[i]] > level);
#endif
		ui old_S_end = S_end, old_R_end = C_end;
		assert(Qv.empty());

		if (choose_zero && SC_rid[0] < C_end && !move_u_to_S_with_prune(0, S_end, C_end, level))
		{
			restore_SC(S_end, C_end, old_S_end, old_R_end, level);
			return;
		}

#ifdef _RR3_
		Reduction_Rule_3(S_end, C_end, level);
		if (!remove_vertices_and_prune(S_end, C_end, level))
		{
			restore_SC(S_end, C_end, old_S_end, old_R_end, level);
			return;
		}
#endif

		if (C_end > best_solution_size && satisfy_degree_requirement(C_end) && vc->satisfy_necessary_requirement(SC, C_end, matrix, n, tau, true))
		{
			store_better_community(C_end);
			restore_SC(S_end, C_end, old_S_end, old_R_end, level);
			return;
		}

		if (C_end <= best_solution_size + 1 || best_solution_size >= _UB_)
		{
			restore_SC(S_end, C_end, old_S_end, old_R_end, level);
			return;
		}

#ifndef NDEBUG
		for (ui i = 0; i < C_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SC[i] * n + SC[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < C_end; j++)
				if (matrix[SC[i] * n + SC[j]])
					++d2;
			assert(d1 == degree_in_S[SC[i]]);
			assert(d2 == degree[SC[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SC[i]] + tau >= S_end);
		for (ui i = S_end; i < C_end; i++)
			assert(degree_in_S[SC[i]] + tau > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SC[i]] + tau == S_end)
			{
				char *t_matrix = matrix + SC[i] * n;
				for (ui j = S_end; j < C_end; j++)
					assert(t_matrix[SC[j]]);
			}
		for (ui i = 0; i < C_end; i++)
			assert(level_id[SC[i]] > level);
#endif

		ui u = n;
#ifdef _FAST_VCL_
		ui min_deg = n;

		for (ui i = 0; i < C_end; i++)
			if (degree[SC[i]] < min_deg)
			{
				u = SC[i];
				min_deg = degree[SC[i]];
			}
		assert(u != n && SC_rid[u] < C_end && min_deg != n);

		if (min_deg >= C_end - tau)
		{

#ifndef NDEBUG
			for (ui i = 0; i < C_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SC[i] * n + SC[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < C_end; j++)
					if (matrix[SC[i] * n + SC[j]])
						++d2;
				assert(d1 == degree_in_S[SC[i]]);
				assert(d2 == degree[SC[i]]);
			}
#endif
			ui pre_S_end = S_end, pre_R_end = C_end;
			branching_with_VCL(S_end, C_end, level + 1, false);
			assert(S_end == pre_S_end && C_end == pre_R_end);
#ifndef NDEBUG
			for (ui i = 0; i < C_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SC[i] * n + SC[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < C_end; j++)
					if (matrix[SC[i] * n + SC[j]])
						++d2;
				assert(d1 == degree_in_S[SC[i]]);
				assert(d2 == degree[SC[i]]);
			}
#endif
			restore_SC(S_end, C_end, old_S_end, old_R_end, level);
			return;
		}

		if (SC_rid[u] < S_end)
			u = choose_branch_vertex_based_on_non_neighbors(S_end, C_end, u);
#else

		ui u2 = n;
		vector<ui> HG_mark;
		vc->get_branching_case_based_on_OSP(SC, C_end, matrix, HG_mark, n, tau, u, u2);
		assert(u != u2 && SC_rid[u] >= 0 && SC_rid[u] < C_end && SC_rid[u2] >= 0 && SC_rid[u2] < C_end);

		if (SC_rid[u] < S_end && SC_rid[u2] < S_end)
		{
			for (ui i = S_end; i < C_end; i++)
				if (HG_mark[SC[i]] == 0)
				{
					u = SC[i];
					break;
				}
		}
		else if (SC_rid[u] < S_end)
		{
			assert(SC_rid[u2] >= S_end && SC_rid[u2] < C_end);
			u = u2;
		}
#endif
		assert(SC_rid[u] >= S_end && SC_rid[u] < C_end);
		assert(degree[u] + tau > best_solution_size && degree[u] + tau > S_end);

		swap_pos(S_end, SC_rid[u]);
		bool can_add_u = vc->satisfy_necessary_requirement(SC, S_end + 1, matrix, n, tau, false);

		if (can_add_u)
		{

			ui pre_best_solution_size = best_solution_size, t_old_S_end = S_end, t_old_R_end = C_end;
			assert(Qv.empty());

#ifndef NDEBUG
			for (ui i = 0; i < C_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SC[i] * n + SC[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < C_end; j++)
					if (matrix[SC[i] * n + SC[j]])
						++d2;
				assert(d1 == degree_in_S[SC[i]]);
				assert(d2 == degree[SC[i]]);
			}
			for (ui i = 0; i < S_end; i++)
				assert(degree_in_S[SC[i]] + tau >= S_end);
			for (ui i = S_end; i < C_end; i++)
				assert(degree_in_S[SC[i]] + tau > S_end);
			for (ui i = 0; i < S_end; i++)
				if (degree_in_S[SC[i]] + tau == S_end)
				{
					char *t_matrix = matrix + SC[i] * n;
					for (ui j = S_end; j < C_end; j++)
						assert(t_matrix[SC[j]]);
				}
			for (ui i = 0; i < C_end; i++)
				assert(level_id[SC[i]] > level);
#endif
			if (move_u_to_S_with_prune(u, S_end, C_end, level))
				branching_with_fast_VCL(S_end, C_end, level + 1, false);

			if (best_solution_size >= _UB_)
				return;
			assert(S_end == t_old_S_end + 1 && SC[S_end - 1] == u);
			restore_SC(S_end, C_end, S_end, t_old_R_end, level);

#ifndef NDEBUG
			for (ui i = 0; i < C_end; i++)
			{
				ui d1 = 0, d2 = 0;
				for (ui j = 0; j < S_end; j++)
					if (matrix[SC[i] * n + SC[j]])
						++d1;
				d2 = d1;
				for (ui j = S_end; j < C_end; j++)
					if (matrix[SC[i] * n + SC[j]])
						++d2;
				assert(d1 == degree_in_S[SC[i]]);
				assert(d2 == degree[SC[i]]);
			}
			for (ui i = 0; i < C_end; i++)
				assert(level_id[SC[i]] > level);
#endif

			assert(Qv.empty());
			bool succeed = remove_u_from_S_with_prune(S_end, C_end, level);
			if (succeed && best_solution_size > pre_best_solution_size)
				succeed = collect_removable_vertices(S_end, C_end, level);
			if (succeed)
				succeed = remove_vertices_and_prune(S_end, C_end, level);

			if (succeed)
				branching_with_fast_VCL(S_end, C_end, level + 1, false);

			if (best_solution_size >= _UB_)
				return;
			assert(S_end >= old_S_end && C_end <= old_R_end);
			restore_SC(S_end, C_end, old_S_end, old_R_end, level);
		}
		else
		{
			ui pre_best_solution_size = best_solution_size;
			assert(Qv.empty());
			bool succeed = remove_u_from_C_with_prune(u, S_end, C_end, level);
			if (succeed && best_solution_size > pre_best_solution_size)
				succeed = collect_removable_vertices(S_end, C_end, level);
			if (succeed)
				succeed = remove_vertices_and_prune(S_end, C_end, level);

			if (succeed)
				branching_with_fast_VCL(S_end, C_end, level + 1, false);

			if (best_solution_size >= _UB_)
				return;
			assert(S_end >= old_S_end && C_end <= old_R_end);
			restore_SC(S_end, C_end, old_S_end, old_R_end, level);
		}
	}

	void branching_with_VCL(ui S_end, ui C_end, ui level, bool choose_zero)
	{

		if (S_end > best_solution_size && vc->satisfy_necessary_requirement(SC, S_end, matrix, n, tau, true))
			store_better_community(S_end);
		if (C_end > best_solution_size && satisfy_degree_requirement(C_end) && vc->satisfy_necessary_requirement(SC, C_end, matrix, n, tau, true))
		{
			store_better_community(C_end);
			return;
		}
		if (C_end <= best_solution_size + 1 || best_solution_size >= _UB_)
			return;

		if (choose_zero && SC_rid[0] < C_end)
		{
			assert(S_end == 0 && SC_rid[0] < C_end);
			addtoS(0, S_end);
		}

		assert(C_end > best_solution_size);

		ui u = n;
		ui u2 = n;
		vector<ui> HG_mark;
		vc->get_branching_case_based_on_OSP(SC, C_end, matrix, HG_mark, n, tau, u, u2);
		assert(u != u2 && SC_rid[u] >= 0 && SC_rid[u] < C_end && SC_rid[u2] >= 0 && SC_rid[u2] < C_end);
		u = SC_rid[u], u2 = SC_rid[u2];
		if (u < S_end && u2 < S_end)
			branching_scenario_1(u, u2, HG_mark, S_end, C_end, level);
		else
		{
			if (u >= S_end)
			{
				ui temp = u;
				u = u2;
				u2 = temp;
			}
			branching_scenarios_2_and_3(u, u2, S_end, C_end, level);
		}
	}

	void restore_deleted_vertices(ui &C_end, ui level, ui &count_restore)
	{
		while (!temp_deleted.empty())
		{
			ui last_deleted = temp_deleted.back();
			if (level_id[last_deleted] == level)
			{
				addtoC(last_deleted, C_end, level);
				temp_deleted.pop_back();
				count_restore++;
			}
			else
				break;
		}
	}

	void branching_scenario_1(ui branching_vertex_1, ui branching_vertex_2, vector<ui> &HG_mark, ui S_end, ui C_end, ui level)
	{
		ui pre_Rend = C_end, u = branching_vertex_1, u2 = branching_vertex_2;
		assert(u < S_end && u2 < S_end);

		ui HG_intersect_S = 0;
		vector<ui> HG_intersect_C;
		get_HG_intersect_S_and_C(HG_mark, HG_intersect_S, HG_intersect_C, S_end, C_end, level);
		assert(!HG_intersect_C.empty() && HG_intersect_S <= tau);

		assert(tau >= HG_intersect_S);
		ui canselect = tau - HG_intersect_S;
		ui pos = n;
		assert(canselect >= 0 && canselect < HG_intersect_C.size());
		bool all = true;
		generate_canselect_branches(all, canselect, S_end, C_end, level, pos, HG_intersect_C);
		generate_the_last_branch(all, canselect, S_end, C_end, level, HG_intersect_C);

		if (pos != n)
			for (ui i = pos; i >= 0 && i <= HG_intersect_C.size(); i--)
				delfrS(HG_intersect_C[i], S_end);
		ui after_R_end = C_end;
		assert(pre_Rend == after_R_end);
	}

	void branching_scenarios_2_and_3(ui u, ui u2, ui S_end, ui C_end, ui level)
	{
		ui pre_Rend = C_end;
		assert((u < S_end && u2 >= S_end) || (u >= S_end && u2 >= S_end));
		ui record_index = SC[u2];

		assert(level_id[record_index] > level);
		if (canadd(record_index, S_end))
		{
			addtoS(record_index, S_end);
			ui count_deleted = 0, count_restore = 0;
			delete_vertices_from_C(S_end, C_end, level, count_deleted);
			branching_with_VCL(S_end, C_end, level + 1, false);
			restore_deleted_vertices(C_end, level, count_restore);
			assert(count_restore == count_deleted);
			delfrS(record_index, S_end);
		}

		assert(SC_rid[record_index] < C_end && SC_rid[record_index] >= S_end);
		delfrC(record_index, C_end, level);
		branching_with_VCL(S_end, C_end, level + 1, false);
		addtoC(record_index, C_end, level);
		assert(pre_Rend == C_end);
	}

	void generate_canselect_branches(bool &all, ui canselect, ui &S_end, ui &C_end, ui level, ui &pos, vector<ui> HG_intersect_C)
	{
		for (ui i = 0; i < canselect; ++i)
		{
			assert(SC_rid[HG_intersect_C[i]] < C_end && SC_rid[HG_intersect_C[i]] >= S_end);
			delfrC(HG_intersect_C[i], C_end, level);

			if (i && !canadd(HG_intersect_C[i - 1], S_end))
			{
				addtoC(HG_intersect_C[i], C_end, level);
				all = false;
				break;
			}
			ui count_deleted = 0;
			if (i)
			{
				assert(canadd(HG_intersect_C[i - 1], S_end));
				addtoS(HG_intersect_C[i - 1], S_end);
				delete_vertices_from_C(S_end, C_end, level, count_deleted);

				pos = i - 1;
			}

			branching_with_VCL(S_end, C_end, level + 1, false);

			ui count_restore = 0;
			if (i)
				restore_deleted_vertices(C_end, level, count_restore);
			assert(count_restore == count_deleted);
			addtoC(HG_intersect_C[i], C_end, level);
		}
	}

	void generate_the_last_branch(bool all, ui canselect, ui &S_end, ui &C_end, ui level, vector<ui> HG_intersect_C)
	{
		if (all && (canselect == 0 || canadd(HG_intersect_C[canselect - 1], S_end)))
		{
			bool t_canadd = false;
			if (canselect)
				t_canadd = canadd(HG_intersect_C[canselect - 1], S_end);
			ui pre_dd = 0;
			for (ui i = canselect; i < (ui)HG_intersect_C.size(); ++i)
			{

				assert(level_id[HG_intersect_C[i]] > level);
				if (level_id[HG_intersect_C[i]] > level)
				{
					assert(SC_rid[HG_intersect_C[i]] < C_end && SC_rid[HG_intersect_C[i]] >= S_end);
					delfrC(HG_intersect_C[i], C_end, level);
					pre_dd++;
				}
			}

			ui count_deleted = 0;
			if (canselect && t_canadd)
			{
				addtoS(HG_intersect_C[canselect - 1], S_end);
				delete_vertices_from_C(S_end, C_end, level, count_deleted);
			}

			branching_with_VCL(S_end, C_end, level + 1, false);

			ui count_restore = 0;
			if (canselect && t_canadd)
			{
				restore_deleted_vertices(C_end, level, count_restore);
				delfrS(HG_intersect_C[canselect - 1], S_end);
			}
			assert(count_restore == count_deleted);
			ui after_dd = 0;
			for (ui i = HG_intersect_C.size() - 1; i >= canselect && i <= HG_intersect_C.size(); i--)
			{
				if (level_id[HG_intersect_C[i]] == level)
				{
					addtoC(HG_intersect_C[i], C_end, level);
					after_dd++;
				}
			}
			assert(pre_dd == after_dd);
		}
	}

	void get_HG_intersect_S_and_C(vector<ui> &HG_mark, ui &HG_intersect_S, vector<ui> &HG_intersect_C, ui S_end, ui C_end, ui level)
	{
		assert(HG_intersect_C.empty());
		for (ui i = 0; i < C_end; i++)
			if (HG_mark[SC[i]] == 0)
				if (i < S_end)
					HG_intersect_S++;
				else
				{
					assert(level_id[SC[i]] > level);
					HG_intersect_C.push_back(SC[i]);
				}

		assert(HG_intersect_S <= tau);
		assert(HG_intersect_C.size() > tau - HG_intersect_S);
	}

	void delete_vertices_from_C(ui S_end, ui &C_end, ui level, ui &count_deleted)
	{
		for (ui i = S_end; i < C_end; i++)
		{
			ui temp_record = SC[i];
			if (level_id[temp_record] > level)
			{
				bool can_delete = false;
				if (degree[temp_record] < best_solution_size - tau)
					can_delete = true;
				if (can_delete)
				{
					temp_deleted.push_back(temp_record);
					assert(level_id[temp_record] > level);
					assert(SC_rid[temp_record] < C_end && SC_rid[temp_record] >= S_end);
					delfrC(temp_record, C_end, level);
					count_deleted++;
				}
			}
		}
	}

	bool canadd(ui u, ui S_end)
	{
		if (SC_rid[u] != S_end)
			swap_pos(S_end, SC_rid[u]);
		return vc->satisfy_necessary_requirement(SC, S_end + 1, matrix, n, tau, false);
	}

	void delfrC(ui u, ui &C_end, ui level)
	{
		assert(level_id[u] == n);
		level_id[u] = level;
		--C_end;

		swap_pos(C_end, SC_rid[u]);

		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < C_end; i++)
			if (t_matrix[SC[i]])
			{
				assert(degree[SC[i]] > 0);
				--degree[SC[i]];
			}
	}

	void addtoC(ui u, ui &C_end, ui level)
	{

		assert(level_id[u] != n && level_id[u] == level);
		level_id[u] = n;
		assert(SC_rid[u] == C_end);
		SC_rid[u] = C_end;
		SC[C_end] = u;
		++C_end;

		char *t_matrix = matrix + u * n;
		degree[u] = 0;
		for (ui i = 0; i < C_end; i++)
			if (t_matrix[SC[i]])
			{
				++degree[SC[i]];
				degree[u]++;
			}
	}

	void addtoS(ui u, ui &S_end)
	{
		swap_pos(S_end, SC_rid[u]);
		++S_end;
	}

	void delfrS(ui u, ui &S_end)
	{
		--S_end;
		assert(S_end == SC_rid[u]);
	}

	void Reduction_Rule_3(ui S_end, ui C_end, ui level)
	{
		vp.resize(C_end - S_end);
		ui max_loss = 0;
		for (ui i = S_end; i < C_end; i++)
		{
			ui loss = 0;
			loss = S_end - degree_in_S[SC[i]];
			if (loss > max_loss)
				max_loss = loss;
			vp[i - S_end].first = SC[i];
			vp[i - S_end].second = loss;
		}
		ui *cnt = neighbors;
		for (ui i = 0; i <= max_loss; i++)
			cnt[i] = 0;
		for (ui i = 0; i < vp.size(); i++)
			++cnt[vp[i].second];
		for (ui i = 0; i < max_loss; i++)
			cnt[i + 1] += cnt[i];
		for (ui i = max_loss; i > 0; i--)
			cnt[i] = cnt[i - 1];
		cnt[0] = 0;
		ui *ids = nonneighbors;
		for (ui i = 0; i < vp.size(); i++)
			ids[cnt[vp[i].second]++] = i;

		ui rest = 0;
		for (ui i = 0; i < S_end; i++)
			rest += tau - S_end + degree_in_S[SC[i]];

		ui new_n = 0;
		while (!Qv.empty())
			Qv.pop();
		for (ui i = 0; i < vp.size(); i++)
		{
			ui idx = ids[i], v = vp[ids[i]].first;
			ui individual_rest = rest - vp[idx].second;
			char *t_matrix = matrix + v * n;
			ui j = 0, v_support = tau - 1 - S_end + degree_in_S[v], ub = S_end + 1;
			while (true)
			{
				if (j == new_n)
					j = i + 1;
				if (j >= vp.size() || ub > best_solution_size || ub + vp.size() - j <= best_solution_size)
					break;
				ui u = vp[ids[j]].first, loss = vp[ids[j]].second;
				if (individual_rest < loss)
					break;
				if (t_matrix[u])
				{
					individual_rest -= loss;
					++ub;
				}
				else if (v_support > 0)
				{
					--v_support;
					individual_rest -= loss;
					++ub;
				}
				++j;
			}
			if (ub <= best_solution_size)
			{
				level_id[v] = level;
				Qv.push(v);
			}
			else
				ids[new_n++] = ids[i];
		}
	}

	void get_neighbors_and_nonneighbors(ui u, ui C_end, ui &neighbors_n, ui &nonneighbors_n)
	{
		neighbors_n = 0;
		nonneighbors_n = 0;
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < C_end; i++)
			if (SC[i] != u)
			{
				if (t_matrix[SC[i]])
					neighbors[neighbors_n++] = SC[i];
				else
					nonneighbors[nonneighbors_n++] = SC[i];
			}
	}

	bool move_u_to_S_with_prune(ui u, ui &S_end, ui &C_end, ui level)
	{
		assert(SC_rid[u] >= S_end && SC_rid[u] < C_end && SC[SC_rid[u]] == u);
		assert(degree_in_S[u] + tau > S_end);
#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SC[i]] + tau >= S_end);
		for (ui i = S_end; i < C_end; i++)
			assert(degree_in_S[SC[i]] + tau > S_end);
		for (ui i = 0; i < C_end; i++)
			assert(level_id[SC[i]] > level);
#endif
		if (SC_rid[u] != S_end)
			swap_pos(S_end, SC_rid[u]);
		++S_end;

		ui neighbors_n = 0, nonneighbors_n = 0;
		get_neighbors_and_nonneighbors(u, C_end, neighbors_n, nonneighbors_n);
		assert(neighbors_n + nonneighbors_n == C_end - 1);
		for (ui i = 0; i < neighbors_n; i++)
			++degree_in_S[neighbors[i]];
#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SC[i]] + tau >= S_end);
#endif

		while (!Qv.empty())
			Qv.pop();

		for (ui i = 0; i < nonneighbors_n; i++)
		{
			ui v = nonneighbors[i];
			if (SC_rid[v] >= S_end)
			{
				if (level_id[v] == level)
					continue;
				bool is_found = false;
#ifdef _RR1_
				if (S_end - degree_in_S[v] >= tau)
				{
					level_id[v] = level;
					Qv.push(v);
					is_found = true;
				}
#endif
#ifdef _RR2_
				if (!is_found && S_end - degree_in_S[u] == tau)
				{
					level_id[v] = level;
					Qv.push(v);
				}
#endif
			}
			else if (S_end - degree_in_S[v] == tau)
			{
#ifdef _RR2_
				char *tt_matrix = matrix + v * n;
				for (ui j = S_end; j < C_end; j++)
					if (level_id[SC[j]] > level && !tt_matrix[SC[j]])
					{
						level_id[SC[j]] = level;
						Qv.push(SC[j]);
					}
#endif
			}
		}

#ifndef NDEBUG
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SC[i]] + tau == S_end)
			{
				char *t_matrix = matrix + SC[i] * n;
				for (ui j = S_end; j < C_end; j++)
					assert(level_id[SC[j]] == level || t_matrix[SC[j]]);
			}
#endif

		return remove_vertices_and_prune(S_end, C_end, level);
	}

	bool remove_vertices_and_prune(ui S_end, ui &C_end, ui level)
	{
		while (!Qv.empty())
		{
			ui u = Qv.front();
			Qv.pop();
			assert(SC[SC_rid[u]] == u && SC_rid[u] >= S_end && SC_rid[u] < C_end);
			--C_end;
			swap_pos(SC_rid[u], C_end);

			bool terminate = false;
			ui neighbors_n = 0;
			char *t_matrix = matrix + u * n;
			for (ui i = 0; i < C_end; i++)
				if (t_matrix[SC[i]])
				{
					ui w = SC[i];
					neighbors[neighbors_n++] = w;
					--degree[w];
					if (degree[w] + tau <= best_solution_size)
					{
						if (i < S_end)
							terminate = true;
						else if (level_id[w] > level)
						{
							level_id[w] = level;
							Qv.push(w);
						}
					}
				}

			if (terminate)
			{

				for (ui i = 0; i < neighbors_n; i++)
					++degree[neighbors[i]];
				level_id[u] = n;
				++C_end;
				return false;
			}
		}
		return true;
	}

	void restore_SC(ui &S_end, ui &C_end, ui old_S_end, ui old_R_end, ui level)
	{
#ifndef NDEBUG
		for (ui i = 0; i < C_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SC[i] * n + SC[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < C_end; j++)
				if (matrix[SC[i] * n + SC[j]])
					++d2;
			assert(d1 == degree_in_S[SC[i]]);
			assert(d2 == degree[SC[i]]);
		}
#endif

		while (!Qv.empty())
		{
			ui u = Qv.front();
			Qv.pop();
			assert(level_id[u] == level);
			assert(SC_rid[u] < C_end);
			level_id[u] = n;
		}

#ifndef NDEBUG
		for (ui i = 0; i < C_end; i++)
		{
			if (level_id[SC[i]] <= level)
				printf("level_id[%u] = %u, level = %u, n = %u\n", SC[i], level_id[SC[i]], level, n);
			assert(level_id[SC[i]] > level);
		}
#endif

		for (; C_end < old_R_end; C_end++)
		{
			ui u = SC[C_end];
			assert(level_id[u] == level && SC_rid[u] == C_end);
			level_id[u] = n;

			ui neighbors_n = 0;
			char *t_matrix = matrix + u * n;
			degree[u] = degree_in_S[u] = 0;
			for (ui i = 0; i < C_end; i++)
				if (t_matrix[SC[i]])
				{
					ui w = SC[i];
					neighbors[neighbors_n++] = w;
					++degree[w];
					++degree[u];
					if (i < S_end)
						++degree_in_S[u];
				}
		}

		for (; S_end > old_S_end; S_end--)
		{
			ui u = SC[S_end - 1];
			assert(SC_rid[u] == S_end - 1);

			ui neighbors_n = 0;
			char *t_matrix = matrix + u * n;
			for (ui i = 0; i < C_end; i++)
				if (t_matrix[SC[i]])
				{
					ui w = SC[i];
					neighbors[neighbors_n++] = w;
					--degree_in_S[w];
				}
		}

#ifndef NDEBUG
		for (ui i = 0; i < C_end; i++)
		{
			ui d1 = 0, d2 = 0;
			for (ui j = 0; j < S_end; j++)
				if (matrix[SC[i] * n + SC[j]])
					++d1;
			d2 = d1;
			for (ui j = S_end; j < C_end; j++)
				if (matrix[SC[i] * n + SC[j]])
					++d2;
			assert(d1 == degree_in_S[SC[i]]);
			assert(d2 == degree[SC[i]]);
		}
		for (ui i = 0; i < S_end; i++)
			assert(degree_in_S[SC[i]] + tau >= S_end);
		for (ui i = S_end; i < C_end; i++)
			assert(old_S_end == S_end || degree_in_S[SC[i]] + tau > S_end);
		for (ui i = 0; i < S_end; i++)
			if (degree_in_S[SC[i]] + tau == S_end)
			{
				char *t_matrix = matrix + SC[i] * n;
				for (ui j = S_end; j < C_end; j++)
					assert(old_S_end == S_end || t_matrix[SC[j]]);
			}
		for (ui i = 0; i < C_end; i++)
			assert(level_id[SC[i]] > level);
#endif
	}

	bool remove_u_from_S_with_prune(ui &S_end, ui &C_end, ui level)
	{
		assert(S_end);
		ui u = SC[S_end - 1];
		--S_end;
		--C_end;
		swap_pos(S_end, C_end);
		level_id[u] = level;

		bool terminate = false;
		ui neighbors_n = 0;
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < C_end; i++)
			if (t_matrix[SC[i]])
				neighbors[neighbors_n++] = SC[i];
		for (ui i = 0; i < neighbors_n; i++)
		{
			ui v = neighbors[i];
			--degree_in_S[v];
			--degree[v];
			if (degree[v] + tau <= best_solution_size)
			{
				if (SC_rid[v] < S_end)
					terminate = true;
				else
				{
					assert(level_id[v] > level);
					level_id[v] = level;
					Qv.push(v);
				}
			}
		}
		if (terminate)
			return false;
		return true;
	}

	bool remove_u_from_C_with_prune(ui deleted_vertex, ui S_end, ui &C_end, ui level)
	{
		assert(S_end);
		ui u = deleted_vertex;
		--C_end;
		swap_pos(SC_rid[u], C_end);
		level_id[u] = level;

		bool terminate = false;
		ui neighbors_n = 0;
		char *t_matrix = matrix + u * n;
		for (ui i = 0; i < C_end; i++)
			if (t_matrix[SC[i]])
				neighbors[neighbors_n++] = SC[i];
		for (ui i = 0; i < neighbors_n; i++)
		{
			ui v = neighbors[i];
			--degree[v];
			if (degree[v] + tau <= best_solution_size)
			{
				if (SC_rid[v] < S_end)
					terminate = true;
				else
				{
					assert(level_id[v] > level);
					level_id[v] = level;
					Qv.push(v);
				}
			}
		}
		if (terminate)
			return false;
		return true;
	}

	bool collect_removable_vertices(ui S_end, ui C_end, ui level)
	{
		for (ui i = 0; i < S_end; i++)
			if (degree[SC[i]] + tau <= best_solution_size)
				return false;

		for (ui i = S_end; i < C_end; i++)
			if (level_id[SC[i]] > level)
			{
				ui v = SC[i];
				if (S_end - degree_in_S[v] >= tau || degree[v] + tau <= best_solution_size)
				{
					level_id[v] = level;
					Qv.push(v);
					continue;
				}
				char *t_matrix = matrix + v * n;
				for (ui j = 0; j < S_end; j++)
					if (S_end - degree_in_S[SC[j]] == tau && !t_matrix[SC[j]])
					{
						level_id[v] = level;
						Qv.push(v);
						break;
					}
			}
		return true;
	}

	void swap_pos(ui i, ui j)
	{
		swap(SC[i], SC[j]);
		SC_rid[SC[i]] = i;
		SC_rid[SC[j]] = j;
	}

	ui choose_branch_vertex_based_on_non_neighbors(ui S_end, ui C_end, ui vertex_in_S)
	{
		assert(SC_rid[vertex_in_S] < S_end);
		ui u = n, min_degree_in_S = n;
		for (ui i = S_end; i < C_end; i++)
		{
			ui v = SC[i];
			if (!matrix[vertex_in_S * n + v])
			{
				if (degree_in_S[v] < min_degree_in_S)
				{
					u = v;
					min_degree_in_S = degree_in_S[v];
				}
				else if (degree_in_S[v] == min_degree_in_S && degree[v] > degree[u])
				{
					u = v;
					min_degree_in_S = degree_in_S[v];
				}
			}
		}
		assert(u != n && !matrix[vertex_in_S * n + u]);
		return u;
	}
};
#endif
