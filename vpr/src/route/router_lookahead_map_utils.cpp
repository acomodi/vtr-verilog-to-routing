#include "router_lookahead_map_utils.h"

#include "globals.h"
#include "vpr_context.h"
#include "vtr_math.h"
#include "route_common.h"

constexpr int DIRECT_CONNECT_SPECIAL_SEG_TYPE = -1;

#define MAX_EXPANSION_LEVEL 1

namespace util {

PQ_Entry::PQ_Entry(
    int set_rr_node_ind,
    int switch_ind,
    float parent_delay,
    float parent_R_upstream,
    float parent_congestion_upstream,
    bool starting_node,
    float Tsw_adjust) {
    this->rr_node_ind = set_rr_node_ind;

    auto& device_ctx = g_vpr_ctx.device();
    this->delay = parent_delay;
    this->congestion_upstream = parent_congestion_upstream;
    this->R_upstream = parent_R_upstream;
    if (!starting_node) {
        float Tsw = device_ctx.rr_switch_inf[switch_ind].Tdel;
        Tsw += Tsw_adjust;
        VTR_ASSERT(Tsw >= 0.f);
        float Rsw = device_ctx.rr_switch_inf[switch_ind].R;
        float Cnode = device_ctx.rr_nodes[set_rr_node_ind].C();
        float Rnode = device_ctx.rr_nodes[set_rr_node_ind].R();

        float T_linear = 0.f;
        if (device_ctx.rr_switch_inf[switch_ind].buffered()) {
            T_linear = Tsw + Rsw * Cnode + 0.5 * Rnode * Cnode;
        } else { /* Pass transistor */
            T_linear = Tsw + 0.5 * Rsw * Cnode;
        }

        float base_cost = 0.f;
        if (device_ctx.rr_switch_inf[switch_ind].configurable()) {
            base_cost = get_single_rr_cong_base_cost(set_rr_node_ind);
        }

        VTR_ASSERT(T_linear >= 0.);
        VTR_ASSERT(base_cost >= 0.);
        this->delay += T_linear;

        this->congestion_upstream += base_cost;
    }

    /* set the cost of this node */
    this->cost = this->delay;
}

util::PQ_Entry_Delay::PQ_Entry_Delay(
    int set_rr_node_ind,
    int switch_ind,
    const util::PQ_Entry_Delay* parent) {
    this->rr_node_ind = set_rr_node_ind;

    if (parent != nullptr) {
        auto& device_ctx = g_vpr_ctx.device();
        float Tsw = device_ctx.rr_switch_inf[switch_ind].Tdel;
        float Rsw = device_ctx.rr_switch_inf[switch_ind].R;
        float Cnode = device_ctx.rr_nodes[set_rr_node_ind].C();
        float Rnode = device_ctx.rr_nodes[set_rr_node_ind].R();

        float T_linear = 0.f;
        if (device_ctx.rr_switch_inf[switch_ind].buffered()) {
            T_linear = Tsw + Rsw * Cnode + 0.5 * Rnode * Cnode;
        } else { /* Pass transistor */
            T_linear = Tsw + 0.5 * Rsw * Cnode;
        }

        VTR_ASSERT(T_linear >= 0.);
        this->delay_cost = parent->delay_cost + T_linear;
    } else {
        this->delay_cost = 0.f;
    }
}

util::PQ_Entry_Base_Cost::PQ_Entry_Base_Cost(
    int set_rr_node_ind,
    int switch_ind,
    const util::PQ_Entry_Base_Cost* parent) {
    this->rr_node_ind = set_rr_node_ind;

    if (parent != nullptr) {
        auto& device_ctx = g_vpr_ctx.device();
        if (device_ctx.rr_switch_inf[switch_ind].configurable()) {
            this->base_cost = parent->base_cost + get_single_rr_cong_base_cost(set_rr_node_ind);
        } else {
            this->base_cost = parent->base_cost;
        }
    } else {
        this->base_cost = 0.f;
    }
}

/* returns cost entry with the smallest delay */
util::Cost_Entry util::Expansion_Cost_Entry::get_smallest_entry() const {
    util::Cost_Entry smallest_entry;

    for (auto entry : this->cost_vector) {
        if (!smallest_entry.valid() || entry.delay < smallest_entry.delay) {
            smallest_entry = entry;
        }
    }

    return smallest_entry;
}

/* returns a cost entry that represents the average of all the recorded entries */
util::Cost_Entry util::Expansion_Cost_Entry::get_average_entry() const {
    float avg_delay = 0;
    float avg_congestion = 0;

    for (auto cost_entry : this->cost_vector) {
        avg_delay += cost_entry.delay;
        avg_congestion += cost_entry.congestion;
    }

    avg_delay /= (float)this->cost_vector.size();
    avg_congestion /= (float)this->cost_vector.size();

    return util::Cost_Entry(avg_delay, avg_congestion);
}

/* returns a cost entry that represents the geomean of all the recorded entries */
util::Cost_Entry util::Expansion_Cost_Entry::get_geomean_entry() const {
    float geomean_delay = 0;
    float geomean_cong = 0;
    for (auto cost_entry : this->cost_vector) {
        geomean_delay += log(cost_entry.delay);
        geomean_cong += log(cost_entry.congestion);
    }

    geomean_delay = exp(geomean_delay / (float)this->cost_vector.size());
    geomean_cong = exp(geomean_cong / (float)this->cost_vector.size());

    return util::Cost_Entry(geomean_delay, geomean_cong);
}

/* returns a cost entry that represents the medial of all recorded entries */
util::Cost_Entry util::Expansion_Cost_Entry::get_median_entry() const {
    /* find median by binning the delays of all entries and then chosing the bin
     * with the largest number of entries */

    int num_bins = 10;

    /* find entries with smallest and largest delays */
    util::Cost_Entry min_del_entry;
    util::Cost_Entry max_del_entry;
    for (auto entry : this->cost_vector) {
        if (!min_del_entry.valid() || entry.delay < min_del_entry.delay) {
            min_del_entry = entry;
        }
        if (!max_del_entry.valid() || entry.delay > max_del_entry.delay) {
            max_del_entry = entry;
        }
    }

    /* get the bin size */
    float delay_diff = max_del_entry.delay - min_del_entry.delay;
    float bin_size = delay_diff / (float)num_bins;

    /* sort the cost entries into bins */
    std::vector<std::vector<util::Cost_Entry>> entry_bins(num_bins, std::vector<util::Cost_Entry>());
    for (auto entry : this->cost_vector) {
        float bin_num = floor((entry.delay - min_del_entry.delay) / bin_size);

        VTR_ASSERT(vtr::nint(bin_num) >= 0 && vtr::nint(bin_num) <= num_bins);
        if (vtr::nint(bin_num) == num_bins) {
            /* largest entry will otherwise have an out-of-bounds bin number */
            bin_num -= 1;
        }
        entry_bins[vtr::nint(bin_num)].push_back(entry);
    }

    /* find the bin with the largest number of elements */
    int largest_bin = 0;
    int largest_size = 0;
    for (int ibin = 0; ibin < num_bins; ibin++) {
        if (entry_bins[ibin].size() > (unsigned)largest_size) {
            largest_bin = ibin;
            largest_size = (unsigned)entry_bins[ibin].size();
        }
    }

    /* get the representative delay of the largest bin */
    util::Cost_Entry representative_entry = entry_bins[largest_bin][0];

    return representative_entry;
}

template<typename Entry>
void expand_dijkstra_neighbours(const t_rr_graph_storage& rr_nodes,
                                const Entry& parent_entry,
                                std::vector<util::Search_Path>* paths,
                                std::vector<bool>* node_expanded,
                                std::priority_queue<Entry,
                                                    std::vector<Entry>,
                                                    std::greater<Entry>>* pq) {
    int parent_ind = size_t(parent_entry.rr_node_ind);

    auto& parent_node = rr_nodes[parent_ind];

    for (int iedge = 0; iedge < parent_node.num_edges(); iedge++) {
        int child_node_ind = parent_node.edge_sink_node(iedge);
        int switch_ind = parent_node.edge_switch(iedge);

        /* skip this child if it has already been expanded from */
        if ((*node_expanded)[child_node_ind]) {
            continue;
        }

        Entry child_entry(child_node_ind, switch_ind, &parent_entry);
        VTR_ASSERT(child_entry.cost() >= 0);

        /* Create (if it doesn't exist) or update (if the new cost is lower)
         * to specified node */
        Search_Path path_entry = {child_entry.cost(), parent_ind, iedge};
        auto& path = (*paths)[child_node_ind];
        if (path_entry.cost < path.cost) {
            pq->push(child_entry);
            path = path_entry;
        }
    }
}

template void expand_dijkstra_neighbours(const t_rr_graph_storage& rr_nodes,
                                         const PQ_Entry_Delay& parent_entry,
                                         std::vector<Search_Path>* paths,
                                         std::vector<bool>* node_expanded,
                                         std::priority_queue<PQ_Entry_Delay,
                                                             std::vector<PQ_Entry_Delay>,
                                                             std::greater<PQ_Entry_Delay>>* pq);
template void expand_dijkstra_neighbours(const t_rr_graph_storage& rr_nodes,
                                         const PQ_Entry_Base_Cost& parent_entry,
                                         std::vector<Search_Path>* paths,
                                         std::vector<bool>* node_expanded,
                                         std::priority_queue<PQ_Entry_Base_Cost,
                                                             std::vector<PQ_Entry_Base_Cost>,
                                                             std::greater<PQ_Entry_Base_Cost>>* pq);

template<typename Entry>
std::pair<float, int> run_dijkstra(int start_node_ind,
                                   std::vector<bool>* node_expanded,
                                   std::vector<util::Search_Path>* paths,
                                   RoutingCosts* routing_costs);

vtr::Point<int> pick_sample_tile(t_physical_tile_type_ptr tile_type, vtr::Point<int> prev) {
    //Very simple for now, just pick the fist matching tile found
    vtr::Point<int> loc(OPEN, OPEN);

    //VTR_LOG("Prev: %d,%d\n", prev.x(), prev.y());

    auto& device_ctx = g_vpr_ctx.device();
    auto& grid = device_ctx.grid;

    int y_init = prev.y() + 1; //Start searching next element above prev

    for (int x = prev.x(); x < int(grid.width()); ++x) {
        if (x < 0) continue;

        //VTR_LOG("  x: %d\n", x);

        for (int y = y_init; y < int(grid.height()); ++y) {
            if (y < 0) continue;

            //VTR_LOG("   y: %d\n", y);
            if (grid[x][y].type->index == tile_type->index) {
                loc.set_x(x);
                loc.set_y(y);
                return loc;
            }
        }

        if (loc.x() != OPEN && loc.y() != OPEN) {
            break;
        } else {
            y_init = 0; //Prepare to search next column
        }
    }
    //VTR_LOG("Next: %d,%d\n", loc.x(), loc.y());

    return loc;
}

void dijkstra_flood_to_wires(int itile, RRNodeId node, t_src_opin_delays& src_opin_delays) {
    auto& device_ctx = g_vpr_ctx.device();
    auto& rr_graph = device_ctx.rr_nodes;

    struct t_pq_entry {
        float delay;
        float congestion;
        RRNodeId node;

        bool operator<(const t_pq_entry& rhs) const {
            return this->delay < rhs.delay;
        }
    };

    std::priority_queue<t_pq_entry> pq;

    t_pq_entry root;
    root.congestion = 0.;
    root.delay = 0.;
    root.node = node;

    int ptc = rr_graph.node_ptc_num(node);

    /*
     * Perform Djikstra from the SOURCE/OPIN of interest, stopping at the the first
     * reachable wires (i.e until we hit the inter-block routing network), or a SINK
     * (via a direct-connect).
     *
     * Note that typical RR graphs are structured :
     *
     *      SOURCE ---> OPIN --> CHANX/CHANY
     *              |
     *              --> OPIN --> CHANX/CHANY
     *              |
     *             ...
     *
     *   possibly with direct-connects of the form:
     *
     *      SOURCE --> OPIN --> IPIN --> SINK
     *
     * and there is a small number of fixed hops from SOURCE/OPIN to CHANX/CHANY or
     * to a SINK (via a direct-connect), so this runs very fast (i.e. O(1))
     */
    pq.push(root);
    while (!pq.empty()) {
        t_pq_entry curr = pq.top();
        pq.pop();

        e_rr_type curr_rr_type = rr_graph.node_type(curr.node);
        if (curr_rr_type == CHANX || curr_rr_type == CHANY || curr_rr_type == SINK) {
            //We stop expansion at any CHANX/CHANY/SINK
            int seg_index;
            if (curr_rr_type != SINK) {
                //It's a wire, figure out its type
                int cost_index = rr_graph.node_cost_index(curr.node);
                seg_index = device_ctx.rr_indexed_data[cost_index].seg_index;
            } else {
                //This is a direct-connect path between an IPIN and OPIN,
                //which terminated at a SINK.
                //
                //We treat this as a 'special' wire type
                seg_index = DIRECT_CONNECT_SPECIAL_SEG_TYPE;
            }

            //Keep costs of the best path to reach each wire type
            if (!src_opin_delays[itile][ptc].count(seg_index)
                || curr.delay < src_opin_delays[itile][ptc][seg_index].delay) {
                src_opin_delays[itile][ptc][seg_index].wire_rr_type = curr_rr_type;
                src_opin_delays[itile][ptc][seg_index].wire_seg_index = seg_index;
                src_opin_delays[itile][ptc][seg_index].delay = curr.delay;
                src_opin_delays[itile][ptc][seg_index].congestion = curr.congestion;
            }

        } else if (curr_rr_type == SOURCE || curr_rr_type == OPIN || curr_rr_type == IPIN) {
            //We allow expansion through SOURCE/OPIN/IPIN types
            int cost_index = rr_graph.node_cost_index(curr.node);
            float incr_cong = device_ctx.rr_indexed_data[cost_index].base_cost; //Current nodes congestion cost

            for (RREdgeId edge : rr_graph.edge_range(curr.node)) {
                int iswitch = rr_graph.edge_switch(edge);
                float incr_delay = device_ctx.rr_switch_inf[iswitch].Tdel;

                RRNodeId next_node = rr_graph.edge_sink_node(edge);

                t_pq_entry next;
                next.congestion = curr.congestion + incr_cong; //Of current node
                next.delay = curr.delay + incr_delay;          //To reach next node
                next.node = next_node;

                pq.push(next);
            }
        } else {
            VPR_ERROR(VPR_ERROR_ROUTE, "Unrecognized RR type");
        }
    }
}

void dijkstra_flood_to_ipins(RRNodeId node, t_chan_ipins_delays& chan_ipins_delays) {
    auto& device_ctx = g_vpr_ctx.device();
    auto& rr_graph = device_ctx.rr_nodes;

    struct t_pq_entry {
        float delay;
        float congestion;
        RRNodeId node;
        int level;

        bool operator<(const t_pq_entry& rhs) const {
            return this->delay < rhs.delay;
        }
    };

    std::priority_queue<t_pq_entry> pq;

    t_pq_entry root;
    root.congestion = 0.;
    root.delay = 0.;
    root.node = node;
    root.level = 0;

    /*
     * Perform Djikstra from the CHAN of interest, stopping at the the first
     * reachable IPIN
     *
     * Note that typical RR graphs are structured :
     *
     *      CHANX/CHANY --> CHANX/CHANY --> ... --> CHANX/CHANY --> IPIN --> SINK
     *                  |
     *                  --> CHANX/CHANY --> ... --> CHANX/CHANY --> IPIN --> SINK
     *                  |
     *                  ...
     *
     * and there is a variable number of hops from a given CHANX/CHANY  to IPIN.
     * To avoid impacting on run-time, a fixed number of hops is performed. This
     * should be enough to find the delay from the last CAHN to IPIN connection.
     */
    pq.push(root);

    float site_pin_delay = std::numeric_limits<float>::infinity();

    while (!pq.empty()) {
        t_pq_entry curr = pq.top();
        pq.pop();

        e_rr_type curr_rr_type = rr_graph.node_type(curr.node);
        if (curr_rr_type == IPIN) {
            int node_x = rr_graph.node_xlow(curr.node);
            int node_y = rr_graph.node_ylow(curr.node);

            auto tile_type = device_ctx.grid[node_x][node_y].type;
            int itile = tile_type->index;

            int ptc = rr_graph.node_ptc_num(curr.node);

            if (ptc >= int(chan_ipins_delays[itile].size())) {
                chan_ipins_delays[itile].resize(ptc + 1); //Inefficient but functional...
            }

            site_pin_delay = std::min(curr.delay, site_pin_delay);
            //Keep costs of the best path to reach each wire type
            chan_ipins_delays[itile][ptc].wire_rr_type = curr_rr_type;
            chan_ipins_delays[itile][ptc].delay = site_pin_delay;
            chan_ipins_delays[itile][ptc].congestion = curr.congestion;
        } else if (curr_rr_type == CHANX || curr_rr_type == CHANY) {
            if (curr.level >= MAX_EXPANSION_LEVEL) {
                continue;
            }

            //We allow expansion through SOURCE/OPIN/IPIN types
            int cost_index = rr_graph.node_cost_index(curr.node);
            float new_cong = device_ctx.rr_indexed_data[cost_index].base_cost; //Current nodes congestion cost

            for (RREdgeId edge : rr_graph.edge_range(curr.node)) {
                int iswitch = rr_graph.edge_switch(edge);
                float new_delay = device_ctx.rr_switch_inf[iswitch].Tdel;

                RRNodeId next_node = rr_graph.edge_sink_node(edge);

                t_pq_entry next;
                next.congestion = new_cong; //Of current node
                next.delay = new_delay;     //To reach next node
                next.node = next_node;
                next.level = curr.level + 1;

                pq.push(next);
            }
        } else {
            VPR_ERROR(VPR_ERROR_ROUTE, "Unrecognized RR type");
        }
    }
}

} // namespace util
