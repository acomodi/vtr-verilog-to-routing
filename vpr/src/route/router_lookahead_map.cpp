#include "router_lookahead_map.h"

#include <vector>
#include <queue>

#include "rr_node.h"
#include "router_lookahead_map_utils.h"
#include "router_lookahead_sampling.h"
#include "globals.h"
#include "vtr_math.h"
#include "vtr_time.h"
#include "vtr_geometry.h"
#include "echo_files.h"
#include "rr_graph.h"

#include "route_timing.h"
#include "route_common.h"

#if defined(VPR_USE_TBB)
#    include <tbb/parallel_for_each.h>
#    include <tbb/mutex.h>
#endif

/* we're profiling routing cost over many tracks for each wire type, so we'll
 * have many cost entries at each |dx|,|dy| offset. There are many ways to
 * "boil down" the many costs at each offset to a single entry for a given
 * (wire type, chan_type) combination we can take the smallest cost, the
 * average, median, etc. This define selects the method we use.
 *
 * NOTE: Currently, only SMALLEST is supported.
 *
 * See e_representative_entry_method */
#define REPRESENTATIVE_ENTRY_METHOD util::SMALLEST

#define CONNECTION_BOX_LOOKAHEAD_MAP_PRINT_COST_MAPS

#define X_OFFSET 2
#define Y_OFFSET 2

/* we will profile delay/congestion using this many tracks for each wire type */
#define MAX_TRACK_OFFSET 16

// Don't continue storing a path after hitting a lower-or-same cost entry.
static constexpr bool BREAK_ON_MISS = false;

static constexpr int MIN_PATH_COUNT = 1000;

static RRNodeId get_start_node(int start_x, int start_y, int target_x, int target_y, t_rr_type rr_type, int seg_index, int track_offset);

//RR position adjustments
static void adjust_rr_position(const RRNodeId rr, int& x, int& y);
static void adjust_rr_pin_position(const RRNodeId rr, int& x, int& y);
static void adjust_rr_wire_position(const RRNodeId rr, int& x, int& y);
static void adjust_rr_src_sink_position(const RRNodeId rr, int& x, int& y);

static void adjust_rr_position(const RRNodeId rr, int& x, int& y) {
    auto& device_ctx = g_vpr_ctx.device();
    auto& rr_graph = device_ctx.rr_nodes;

    e_rr_type rr_type = rr_graph.node_type(rr);

    if (is_chan(rr_type)) {
        adjust_rr_wire_position(rr, x, y);
    } else if (is_pin(rr_type)) {
        adjust_rr_pin_position(rr, x, y);
    } else {
        VTR_ASSERT_SAFE(is_src_sink(rr_type));
        adjust_rr_src_sink_position(rr, x, y);
    }
}

static void adjust_rr_pin_position(const RRNodeId rr, int& x, int& y) {
    /*
     * VPR uses a co-ordinate system where wires above and to the right of a block
     * are at the same location as the block:
     *
     *
     *       <-----------C
     *    D
     *    |  +---------+  ^
     *    |  |         |  |
     *    |  |  (1,1)  |  |
     *    |  |         |  |
     *    V  +---------+  |
     *                    A
     *     B----------->
     *
     * So wires are located as follows:
     *
     *      A: (1, 1) CHANY
     *      B: (1, 0) CHANX
     *      C: (1, 1) CHANX
     *      D: (0, 1) CHANY
     *
     * But all pins incident on the surrounding channels
     * would be at (1,1) along with a relevant side.
     *
     * In the following, we adjust the positions of pins to
     * account for the channel they are incident too.
     *
     * Note that blocks at (0,*) such as IOs, may have (unconnected)
     * pins on the left side, so we also clip the minimum x to zero.
     * Similarly for blocks at (*,0) we clip the minimum y to zero.
     */
    auto& device_ctx = g_vpr_ctx.device();
    auto& rr_graph = device_ctx.rr_nodes;

    VTR_ASSERT_SAFE(is_pin(rr_graph.node_type(rr)));
    VTR_ASSERT_SAFE(rr_graph.node_xlow(rr) == rr_graph.node_xhigh(rr));
    VTR_ASSERT_SAFE(rr_graph.node_ylow(rr) == rr_graph.node_yhigh(rr));

    x = rr_graph.node_xlow(rr);
    y = rr_graph.node_ylow(rr);

    e_side rr_side = rr_graph.node_side(rr);

    if (rr_side == LEFT) {
        x -= 1;
        x = std::max(x, 0);
    } else if (rr_side == BOTTOM && y > 0) {
        y -= 1;
        y = std::max(y, 0);
    }
}

static void adjust_rr_wire_position(const RRNodeId rr, int& x, int& y) {
    auto& device_ctx = g_vpr_ctx.device();
    auto& rr_graph = device_ctx.rr_nodes;

    VTR_ASSERT_SAFE(is_chan(rr_graph.node_type(rr)));

    e_direction rr_dir = rr_graph.node_direction(rr);

    if (rr_dir == DEC_DIRECTION) {
        x = rr_graph.node_xhigh(rr);
        y = rr_graph.node_yhigh(rr);
    } else if (rr_dir == INC_DIRECTION) {
        x = rr_graph.node_xlow(rr);
        y = rr_graph.node_ylow(rr);
    } else {
        VTR_ASSERT_SAFE(rr_dir == BI_DIRECTION);
        //Not sure what to do here...
        //Try average for now.
        x = vtr::nint((rr_graph.node_xlow(rr) + rr_graph.node_xhigh(rr)) / 2.);
        y = vtr::nint((rr_graph.node_ylow(rr) + rr_graph.node_yhigh(rr)) / 2.);
    }
}

static void adjust_rr_src_sink_position(const RRNodeId rr, int& x, int& y) {
    //SOURCE/SINK nodes assume the full dimensions of their
    //associated block
    //
    //Use the average position.
    auto& device_ctx = g_vpr_ctx.device();
    auto& rr_graph = device_ctx.rr_nodes;

    VTR_ASSERT_SAFE(is_src_sink(rr_graph.node_type(rr)));

    x = vtr::nint((rr_graph.node_xlow(rr) + rr_graph.node_xhigh(rr)) / 2.);
    y = vtr::nint((rr_graph.node_ylow(rr) + rr_graph.node_yhigh(rr)) / 2.);
}

std::pair<float, float> MapLookahead::get_src_opin_delays(RRNodeId from_node, int delta_x, int delta_y, float criticality_fac) const {
    auto& device_ctx = g_vpr_ctx.device();
    auto& rr_graph = device_ctx.rr_nodes;

    e_rr_type from_type = rr_graph.node_type(from_node);
    if (from_type == SOURCE || from_type == OPIN) {
        //When estimating costs from a SOURCE/OPIN we look-up to find which wire types (and the
        //cost to reach them) in f_src_opin_delays. Once we know what wire types are
        //reachable, we query the f_wire_cost_map (i.e. the wire lookahead) to get the final
        //delay to reach the sink.

        t_physical_tile_type_ptr tile_type = device_ctx.grid[rr_graph.node_xlow(from_node)][rr_graph.node_ylow(from_node)].type;
        auto tile_index = std::distance(&device_ctx.physical_tile_types[0], tile_type);

        auto from_ptc = rr_graph.node_ptc_num(from_node);

        if (this->src_opin_delays[tile_index][from_ptc].empty()) {
            //During lookahead profiling we were unable to find any wires which connected
            //to this PTC.
            //
            //This can sometimes occur at very low channel widths (e.g. during min W search on
            //small designs) where W discretization combined with fraction Fc may cause some
            //pins/sources to be left disconnected.
            //
            //Such RR graphs are of course unroutable, but that should be determined by the
            //router. So just return an arbitrary value here rather than error.

            return std::pair<float, float>(0.f, 0.f);
        } else {
            //From the current SOURCE/OPIN we look-up the wiretypes which are reachable
            //and then add the estimates from those wire types for the distance of interest.
            //If there are multiple options we use the minimum value.

            float delay = 0;
            float congestion = 0;
            float expected_cost = std::numeric_limits<float>::infinity();

            for (const auto& kv : this->src_opin_delays[tile_index][from_ptc]) {
                const util::t_reachable_wire_inf& reachable_wire_inf = kv.second;

                util::Cost_Entry cost_entry;
                if (reachable_wire_inf.wire_rr_type == SINK) {
                    //Some pins maybe reachable via a direct (OPIN -> IPIN) connection.
                    //In the lookahead, we treat such connections as 'special' wire types
                    //with no delay/congestion cost
                    cost_entry.delay = 0;
                    cost_entry.congestion = 0;
                } else {
                    //For an actual accessible wire, we query the wire look-up to get it's
                    //delay and congestion cost estimates
                    cost_entry = cost_map_.find_cost(reachable_wire_inf.wire_seg_index, reachable_wire_inf.wire_rr_type, delta_x, delta_y);
                }

                float this_cost = (criticality_fac) * (reachable_wire_inf.delay + cost_entry.delay)
                                  + (1. - criticality_fac) * (reachable_wire_inf.congestion + cost_entry.congestion);

                if (this_cost < expected_cost) {
                    delay = reachable_wire_inf.delay;
                    congestion = reachable_wire_inf.congestion;
                }
            }

            return std::pair<float, float>(delay, congestion);
        }

        VTR_ASSERT_SAFE_MSG(false,
                            vtr::string_fmt("Lookahead failed to estimate cost from %s: %s",
                                            rr_node_arch_name(size_t(from_node)).c_str(),
                                            describe_rr_node(size_t(from_node)).c_str())
                                .c_str());
    }

    return std::pair<float, float>(0.f, 0.f);
}

float MapLookahead::get_chan_ipin_delays(RRNodeId to_node) const {
    auto& device_ctx = g_vpr_ctx.device();
    auto& rr_graph = device_ctx.rr_nodes;

    auto to_tile_type = device_ctx.grid[rr_graph.node_xlow(to_node)][rr_graph.node_ylow(to_node)].type;
    auto to_tile_index = to_tile_type->index;

    auto to_ptc = rr_graph.node_ptc_num(to_node);

    if (this->chan_ipins_delays[to_tile_index].size() != 0) {
        auto delay = this->chan_ipins_delays[to_tile_index][to_ptc].delay;

        if (std::isfinite(delay)) {
            return delay;
        }
    }

    return 0.f;
}

// derive a cost from the map between two nodes
float MapLookahead::get_map_cost(int from_node_ind,
                                 int to_node_ind,
                                 float criticality_fac) const {
    if (from_node_ind == to_node_ind) {
        return 0.f;
    }

    auto& device_ctx = g_vpr_ctx.device();

    auto from_node_type = device_ctx.rr_nodes[to_node_ind].type();

    RRNodeId to_node(to_node_ind);
    RRNodeId from_node(from_node_ind);

    int dx, dy;
    std::tie(dx, dy) = this->get_xy_deltas(from_node, to_node);

    float extra_delay, extra_congestion;
    std::tie(extra_delay, extra_congestion) = this->get_src_opin_delays(from_node, dx, dy, criticality_fac);

    int from_seg_index = cost_map_.node_to_segment(from_node_ind);
    util::Cost_Entry cost_entry = cost_map_.find_cost(from_seg_index, from_node_type, dx, dy);

    if (!cost_entry.valid()) {
        // there is no route
        VTR_LOGV_DEBUG(f_router_debug,
                       "Not connected %d (%s, %d) -> %d (%s)\n",
                       from_node_ind, device_ctx.rr_nodes[from_node_ind].type_string(), from_seg_index,
                       to_node_ind, device_ctx.rr_nodes[to_node_ind].type_string());
        return std::numeric_limits<float>::infinity();
    }

    float expected_delay = cost_entry.delay + extra_delay;
    float expected_congestion = cost_entry.congestion + extra_congestion;

    float site_pin_delay = this->get_chan_ipin_delays(to_node);
    expected_delay += site_pin_delay;

    float expected_cost = criticality_fac * expected_delay + (1.0 - criticality_fac) * expected_congestion;

    VTR_LOGV_DEBUG(f_router_debug, "Requested lookahead from node %d to %d\n", from_node_ind, to_node_ind);
    const std::string& segment_name = device_ctx.segment_inf[from_seg_index].name;
    VTR_LOGV_DEBUG(f_router_debug, "Lookahead returned %s (%d) with distance (%zd, %zd)\n",
                   segment_name.c_str(), from_seg_index,
                   dx, dy);
    VTR_LOGV_DEBUG(f_router_debug, "Lookahead delay: %g\n", expected_delay);
    VTR_LOGV_DEBUG(f_router_debug, "Lookahead congestion: %g\n", expected_congestion);
    VTR_LOGV_DEBUG(f_router_debug, "Criticality: %g\n", criticality_fac);
    VTR_LOGV_DEBUG(f_router_debug, "Lookahead cost: %g\n", expected_cost);
    VTR_LOGV_DEBUG(f_router_debug, "Site pin delay: %g\n", site_pin_delay);

    if (!std::isfinite(expected_cost)) {
        VTR_LOG_ERROR("infinite cost for segment %d at (%d, %d)\n", from_seg_index, (int)dx, (int)dy);
        VTR_ASSERT(0);
    }

    if (expected_cost < 0.f) {
        VTR_LOG_ERROR("negative cost for segment %d at (%d, %d)\n", from_seg_index, (int)dx, (int)dy);
        VTR_ASSERT(0);
    }

    return expected_cost;
}

// add a best cost routing path from start_node_ind to node_ind to routing costs
template<typename Entry>
bool MapLookahead::add_paths(int start_node_ind,
                             Entry current,
                             const std::vector<util::Search_Path>& paths,
                             util::RoutingCosts* routing_costs) const {
    auto& device_ctx = g_vpr_ctx.device();
    auto& rr_graph = device_ctx.rr_nodes;

    int node_ind = current.rr_node_ind;
    RRNodeId node(node_ind);

    bool new_sample_found = false;

    // reconstruct the path
    std::vector<int> path;
    for (int i = paths[node_ind].parent; i != start_node_ind; i = paths[i].parent) {
        VTR_ASSERT(i != -1);
        path.push_back(i);
    }
    path.push_back(start_node_ind);

    float site_pin_delay = this->get_chan_ipin_delays(node);
    current.adjust_Tsw(-site_pin_delay);

    // add each node along the path subtracting the incremental costs from the current costs
    Entry start_to_here(start_node_ind, UNDEFINED, nullptr);
    int parent = start_node_ind;
    for (auto it = path.rbegin(); it != path.rend(); it++) {
        RRNodeId this_node(*it);
        auto& here = device_ctx.rr_nodes[*it];
        int seg_index = device_ctx.rr_indexed_data[here.cost_index()].seg_index;

        auto chan_type = rr_graph.node_type(node);

        int ichan = 0;
        if (chan_type == CHANY) {
            ichan = 1;
        }

        int delta_x, delta_y;
        std::tie(delta_x, delta_y) = this->get_xy_deltas(this_node, node);
        vtr::Point<int> delta(delta_x, delta_y);

        util::RoutingCostKey key = {
            ichan,
            seg_index,
            delta};

        if (*it != start_node_ind) {
            auto& parent_node = device_ctx.rr_nodes[parent];
            start_to_here = Entry(*it, parent_node.edge_switch(paths[*it].edge), &start_to_here);
            parent = *it;
        }

        float cost = current.cost() - start_to_here.cost();
        if (cost < 0.f && cost > -10e-15 /* 10 femtosecond */) {
            cost = 0.f;
        }

        VTR_ASSERT(cost >= 0.f);

        // NOTE: implements REPRESENTATIVE_ENTRY_METHOD == SMALLEST
        auto result = routing_costs->insert(std::make_pair(key, cost));
        if (!result.second) {
            if (cost < result.first->second) {
                result.first->second = cost;
                new_sample_found = true;
            } else if (BREAK_ON_MISS) {
                break;
            }
        } else {
            new_sample_found = true;
        }
    }
    return new_sample_found;
}

/* Runs Dijkstra's algorithm from specified node until all nodes have been
 * visited. Each time a pin is visited, the delay/congestion information
 * to that pin is stored to an entry in the routing_cost_map.
 *
 * Returns the maximum (last) minimum cost path stored, and
 * the number of paths from start_node_ind stored. */
template<typename Entry>
std::pair<float, int> MapLookahead::run_dijkstra(int start_node_ind,
                                                 std::vector<bool>* node_expanded,
                                                 std::vector<util::Search_Path>* paths,
                                                 util::RoutingCosts* routing_costs) const {
    auto& device_ctx = g_vpr_ctx.device();
    int path_count = 0;

    /* a list of boolean flags (one for each rr node) to figure out if a
     * certain node has already been expanded */
    std::fill(node_expanded->begin(), node_expanded->end(), false);
    /* For each node keep a list of the cost with which that node has been
     * visited (used to determine whether to push a candidate node onto the
     * expansion queue.
     * Also store the parent node so we can reconstruct a specific path. */
    std::fill(paths->begin(), paths->end(), util::Search_Path{std::numeric_limits<float>::infinity(), -1, -1});
    /* a priority queue for expansion */
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;

    /* first entry has no upstream delay or congestion */
    Entry first_entry(start_node_ind, UNDEFINED, nullptr);
    float max_cost = first_entry.cost();

    pq.push(first_entry);

    /* now do routing */
    while (!pq.empty()) {
        auto current = pq.top();
        pq.pop();

        int node_ind = current.rr_node_ind;

        /* check that we haven't already expanded from this node */
        if ((*node_expanded)[node_ind]) {
            continue;
        }

        /* if this node is an ipin record its congestion/delay in the routing_cost_map */
        if (device_ctx.rr_nodes[node_ind].type() == IPIN) {
            // the last cost should be the highest
            max_cost = current.cost();

            path_count++;
            this->add_paths<Entry>(start_node_ind, current, *paths, routing_costs);
        } else {
            util::expand_dijkstra_neighbours(device_ctx.rr_nodes,
                                             current, paths, node_expanded, &pq);
            (*node_expanded)[node_ind] = true;
        }
    }
    return std::make_pair(max_cost, path_count);
}

/* returns index of a node from which to start routing */
static RRNodeId get_start_node(int start_x, int start_y, int target_x, int target_y, t_rr_type rr_type, int seg_index, int track_offset) {
    auto& device_ctx = g_vpr_ctx.device();
    auto& rr_graph = device_ctx.rr_nodes;

    int result = UNDEFINED;

    if (rr_type != CHANX && rr_type != CHANY) {
        VPR_FATAL_ERROR(VPR_ERROR_ROUTE, "Must start lookahead routing from CHANX or CHANY node\n");
    }

    /* determine which direction the wire should go in based on the start & target coordinates */
    e_direction direction = INC_DIRECTION;
    if ((rr_type == CHANX && target_x < start_x) || (rr_type == CHANY && target_y < start_y)) {
        direction = DEC_DIRECTION;
    }

    int start_lookup_x = start_x;
    int start_lookup_y = start_y;
    if (rr_type == CHANX) {
        //Bizarely, rr_node_indices stores CHANX with swapped x/y...
        std::swap(start_lookup_x, start_lookup_y);
    }

    const std::vector<int>& channel_node_list = device_ctx.rr_node_indices[rr_type][start_lookup_x][start_lookup_y][0];

    /* find first node in channel that has specified segment index and goes in the desired direction */
    for (unsigned itrack = 0; itrack < channel_node_list.size(); itrack++) {
        int node_ind = channel_node_list[itrack];
        if (node_ind < 0) continue;

        RRNodeId node_id(node_ind);

        VTR_ASSERT(rr_graph.node_type(node_id) == rr_type);

        e_direction node_direction = rr_graph.node_direction(node_id);
        int node_cost_ind = rr_graph.node_cost_index(node_id);
        int node_seg_ind = device_ctx.rr_indexed_data[node_cost_ind].seg_index;

        if ((node_direction == direction || node_direction == BI_DIRECTION) && node_seg_ind == seg_index) {
            /* found first track that has the specified segment index and goes in the desired direction */
            result = node_ind;
            if (track_offset == 0) {
                break;
            }
            track_offset -= 2;
        }
    }

    return RRNodeId(result);
}

/* returns the absolute delta_x and delta_y offset required to reach to_node from from_node */
std::pair<int, int> MapLookahead::get_xy_deltas(const RRNodeId from_node, const RRNodeId to_node) const {
    auto& device_ctx = g_vpr_ctx.device();
    auto& rr_graph = device_ctx.rr_nodes;

    e_rr_type from_type = rr_graph.node_type(from_node);
    e_rr_type to_type = rr_graph.node_type(to_node);

    int delta_x, delta_y;

    if (!is_chan(from_type) && !is_chan(to_type)) {
        //Alternate formulation for non-channel types
        int from_x = 0;
        int from_y = 0;
        adjust_rr_position(from_node, from_x, from_y);

        int to_x = 0;
        int to_y = 0;
        adjust_rr_position(to_node, to_x, to_y);

        delta_x = to_x - from_x;
        delta_y = to_y - from_y;
    } else {
        //Traditional formulation

        /* get chan/seg coordinates of the from/to nodes. seg coordinate is along the wire,
         * chan coordinate is orthogonal to the wire */
        int from_seg_low;
        int from_seg_high;
        int from_chan;
        int to_seg;
        int to_chan;
        if (from_type == CHANY) {
            from_seg_low = rr_graph.node_ylow(from_node);
            from_seg_high = rr_graph.node_yhigh(from_node);
            from_chan = rr_graph.node_xlow(from_node);
            to_seg = rr_graph.node_ylow(to_node);
            to_chan = rr_graph.node_xlow(to_node);
        } else {
            from_seg_low = rr_graph.node_xlow(from_node);
            from_seg_high = rr_graph.node_xhigh(from_node);
            from_chan = rr_graph.node_ylow(from_node);
            to_seg = rr_graph.node_xlow(to_node);
            to_chan = rr_graph.node_ylow(to_node);
        }

        /* now we want to count the minimum number of *channel segments* between the from and to nodes */
        int delta_seg, delta_chan;

        /* orthogonal to wire */
        int no_need_to_pass_by_clb = 0; //if we need orthogonal wires then we don't need to pass by the target CLB along the current wire direction
        if (to_chan > from_chan + 1) {
            /* above */
            delta_chan = to_chan - from_chan;
            no_need_to_pass_by_clb = 1;
        } else if (to_chan < from_chan) {
            /* below */
            delta_chan = from_chan - to_chan + 1;
            no_need_to_pass_by_clb = 1;
        } else {
            /* adjacent to current channel */
            delta_chan = 0;
            no_need_to_pass_by_clb = 0;
        }

        /* along same direction as wire. */
        if (to_seg > from_seg_high) {
            /* ahead */
            delta_seg = to_seg - from_seg_high - no_need_to_pass_by_clb;
        } else if (to_seg < from_seg_low) {
            /* behind */
            delta_seg = from_seg_low - to_seg - no_need_to_pass_by_clb;
        } else {
            /* along the span of the wire */
            delta_seg = 0;
        }

        /* account for wire direction. lookahead map was computed by looking up and to the right starting at INC wires. for targets
         * that are opposite of the wire direction, let's add 1 to delta_seg */
        e_direction from_dir = rr_graph.node_direction(from_node);
        if (is_chan(from_type)
            && ((to_seg < from_seg_low && from_dir == INC_DIRECTION) || (to_seg > from_seg_high && from_dir == DEC_DIRECTION))) {
            delta_seg++;
        }

        if (from_type == CHANY) {
            delta_x = delta_chan;
            delta_y = delta_seg;
        } else {
            delta_x = delta_seg;
            delta_y = delta_chan;
        }
    }

    VTR_ASSERT_SAFE(std::abs(delta_x) < (int)device_ctx.grid.width());
    VTR_ASSERT_SAFE(std::abs(delta_y) < (int)device_ctx.grid.height());

    return std::make_pair(delta_x, delta_y);
}

// compute the cost maps for lookahead
void MapLookahead::compute(const std::vector<t_segment_inf>& segment_inf) {
    this->compute_router_src_opin_lookahead();
    this->compute_router_chan_ipin_lookahead();

    auto& device_ctx = g_vpr_ctx.device();
    auto& grid = device_ctx.grid;
    auto& rr_graph = device_ctx.rr_nodes;
    vtr::ScopedStartFinishTimer timer("Computing lookahead map");

    // Initialize rr_node_route_inf if not already
    alloc_and_load_rr_node_route_structs();

    size_t num_segments = segment_inf.size();
    //std::vector<SampleRegion> sample_regions = find_sample_regions(num_segments);

    int longest_length = 0;
    for (const auto& seg_inf : segment_inf) {
        longest_length = std::max(longest_length, seg_inf.length);
    }

    //Start sampling at the lower left non-corner
    int ref_x = 1;
    int ref_y = 1;

    //Sample from locations near the reference location (to capture maximum distance paths)
    //Also sample from locations at least the longest wire length away from the edge (to avoid
    //edge effects for shorter distances)
    std::vector<int> ref_increments = {0, 1,
                                       longest_length, longest_length + 1};

    //Uniquify the increments (avoid sampling the same locations repeatedly if they happen to
    //overlap)
    std::sort(ref_increments.begin(), ref_increments.end());
    ref_increments.erase(std::unique(ref_increments.begin(), ref_increments.end()), ref_increments.end());

    //Upper right non-corner
    int target_x = device_ctx.grid.width() - 2;
    int target_y = device_ctx.grid.height() - 2;

    std::vector<std::map<t_rr_type, std::vector<RRNodeId>>> samples;
    //Profile each wire segment type
    for (int iseg = 0; iseg < int(segment_inf.size()); iseg++) {
        //First try to pick good representative sample locations for each type
        std::map<t_rr_type, std::vector<RRNodeId>> sample_nodes;
        for (e_rr_type chan_type : {CHANX, CHANY}) {
            for (int ref_inc : ref_increments) {
                int sample_x = ref_x + ref_inc;
                int sample_y = ref_y + ref_inc;

                if (sample_x >= int(grid.width())) continue;
                if (sample_y >= int(grid.height())) continue;

                for (int track_offset = 0; track_offset < MAX_TRACK_OFFSET; track_offset += 2) {
                    /* get the rr node index from which to start routing */
                    RRNodeId start_node = get_start_node(sample_x, sample_y,
                                                         target_x, target_y, //non-corner upper right
                                                         chan_type, iseg, track_offset);

                    if (!start_node) {
                        continue;
                    }

                    sample_nodes[chan_type].push_back(RRNodeId(start_node));
                }
            }
        }

        //If we failed to find any representative sample locations, search exhaustively
        //
        //This is to ensure we sample 'unusual' wire types which may not exist in all channels
        //(e.g. clock routing)
        for (e_rr_type chan_type : {CHANX, CHANY}) {
            if (!sample_nodes[chan_type].empty()) continue;

            //Try an exhaustive search to find a suitable sample point
            for (int inode = 0; inode < int(device_ctx.rr_nodes.size()); ++inode) {
                auto rr_node = RRNodeId(inode);
                auto rr_type = rr_graph.node_type(rr_node);
                if (rr_type != chan_type) continue;

                int cost_index = rr_graph.node_cost_index(rr_node);
                VTR_ASSERT(cost_index != OPEN);

                int seg_index = device_ctx.rr_indexed_data[cost_index].seg_index;

                if (seg_index == iseg) {
                    sample_nodes[chan_type].push_back(RRNodeId(inode));
                }

                if (sample_nodes[chan_type].size() >= ref_increments.size()) {
                    break;
                }
            }
        }

        samples.push_back(sample_nodes);
    }

    /* free previous delay map and allocate new one */
    cost_map_.set_counts(segment_inf.size());
    cost_map_.build_segment_map();

    VTR_ASSERT(REPRESENTATIVE_ENTRY_METHOD == util::SMALLEST);
    util::RoutingCosts all_delay_costs;
    util::RoutingCosts all_base_costs;

    /* run Dijkstra's algorithm for each segment type & channel type combination */

#if defined(VPR_USE_TBB)
    tbb::mutex all_costs_mutex;
    tbb::parallel_for_each(samples, [&](const std::map<t_rr_type, std::vector<RRNodeId>>& sample) {
#else
    for (const auto& sample : samples) {
#endif
        // holds the cost entries for a run
        util::RoutingCosts delay_costs;
        util::RoutingCosts base_costs;
        int total_path_count = 0;
        std::vector<bool> node_expanded(device_ctx.rr_nodes.size());
        std::vector<util::Search_Path> paths(device_ctx.rr_nodes.size());
        for (auto& entry : sample) {
            // statistics

            vtr::Timer run_timer;
            float max_delay_cost = 0.f;
            float max_base_cost = 0.f;
            int path_count = 0;
            int seg_index = OPEN;
            for (auto node : entry.second) {
                auto node_ind = size_t(node);

                int cost_index = rr_graph.node_cost_index(node);
                seg_index = device_ctx.rr_indexed_data[cost_index].seg_index;
                {
                    auto result = run_dijkstra<util::PQ_Entry_Delay>(node_ind, &node_expanded, &paths, &delay_costs);
                    max_delay_cost = std::max(max_delay_cost, result.first);
                    path_count += result.second;
                }
                {
                    auto result = run_dijkstra<util::PQ_Entry_Base_Cost>(node_ind, &node_expanded, &paths, &base_costs);
                    max_base_cost = std::max(max_base_cost, result.first);
                    path_count += result.second;
                }
            }

            if (path_count > 0) {
                VTR_LOG("Expanded %d paths of segment type %s(%d), max_cost %e %e (%g paths/sec)\n",
                        path_count, segment_inf[seg_index].name.c_str(), seg_index,
                        max_delay_cost, max_base_cost,
                        path_count / run_timer.elapsed_sec());
            }

            if (path_count == 0) {
                for (auto node : entry.second) {
                    VTR_LOG("Expanded node %s\n", describe_rr_node(size_t(node)).c_str());
                }
            }

            total_path_count += path_count;
            if (total_path_count > MIN_PATH_COUNT) {
                break;
            }
        }

#if defined(VPR_USE_TBB)
        all_costs_mutex.lock();
#endif
        // combine the cost map from this run with the final cost maps for each segment
        for (const auto& cost : delay_costs) {
            const auto& val = cost.second;
            auto result = all_delay_costs.insert(std::make_pair(cost.first, val));
            if (!result.second) {
                // implements REPRESENTATIVE_ENTRY_METHOD == SMALLEST
                result.first->second = std::min(result.first->second, val);
            }
        }
        for (const auto& cost : base_costs) {
            const auto& val = cost.second;
            auto result = all_base_costs.insert(std::make_pair(cost.first, val));
            if (!result.second) {
                // implements REPRESENTATIVE_ENTRY_METHOD == SMALLEST
                result.first->second = std::min(result.first->second, val);
            }
        }
#if defined(VPR_USE_TBB)
        all_costs_mutex.unlock();
    });
#else
    }
#endif

    VTR_LOG("Combining results\n");
    /* boil down the cost list in routing_cost_map at each coordinate to a
     * representative cost entry and store it in the lookahead cost map */
    cost_map_.set_cost_map(all_delay_costs, all_base_costs);

// diagnostics
#if defined(CONNECTION_BOX_LOOKAHEAD_MAP_PRINT_COST_ENTRIES)
    for (auto& cost : all_costs) {
        const auto& key = cost.first;
        const auto& val = cost.second;
        VTR_LOG("%d -> %d (%d, %d): %g, %g\n",
                val.from_node, val.to_node,
                key.delta.x(), key.delta.y(),
                val.cost_entry.delay, val.cost_entry.congestion);
    }
#endif

#if defined(CONNECTION_BOX_LOOKAHEAD_MAP_PRINT_COST_MAPS)
    for (int iseg = 0; iseg < (ssize_t)num_segments; iseg++) {
        VTR_LOG("cost map for %s(%d)\n",
                segment_inf[iseg].name.c_str(), iseg);
        cost_map_.print(iseg);
    }
#endif

#if defined(CONNECTION_BOX_LOOKAHEAD_MAP_PRINT_EMPTY_MAPS)
    for (std::pair<int, int> p : cost_map_.list_empty()) {
        int ichan, iseg;
        std::tie(ichan, iseg) = p;
        VTR_LOG("cost map for %s(%d), chan %d EMPTY\n",
                segment_inf[iseg].name.c_str(), iseg, box_id);
    }
#endif
}

void MapLookahead::compute_router_src_opin_lookahead() const {
    vtr::ScopedStartFinishTimer timer("Computing src/opin lookahead");
    auto& device_ctx = g_vpr_ctx.device();
    auto& rr_graph = device_ctx.rr_nodes;

    this->src_opin_delays.clear();

    this->src_opin_delays.resize(device_ctx.physical_tile_types.size());

    std::vector<int> rr_nodes_at_loc;

    //We assume that the routing connectivity of each instance of a physical tile is the same,
    //and so only measure one instance of each type
    for (size_t itile = 0; itile < device_ctx.physical_tile_types.size(); ++itile) {
        for (e_rr_type rr_type : {SOURCE, OPIN}) {
            vtr::Point<int> sample_loc(-1, -1);

            size_t num_sampled_locs = 0;
            bool ptcs_with_no_delays = true;
            while (ptcs_with_no_delays) { //Haven't found wire connected to ptc
                ptcs_with_no_delays = false;

                sample_loc = util::pick_sample_tile(&device_ctx.physical_tile_types[itile], sample_loc);

                if (sample_loc.x() == -1 && sample_loc.y() == -1) {
                    //No untried instances of the current tile type left
                    VTR_LOG_WARN("Found no %ssample locations for %s in %s\n",
                                 (num_sampled_locs == 0) ? "" : "more ",
                                 rr_node_typename[rr_type],
                                 device_ctx.physical_tile_types[itile].name);
                    break;
                }

                //VTR_LOG("Sampling %s at (%d,%d)\n", device_ctx.physical_tile_types[itile].name, sample_loc.x(), sample_loc.y());

                rr_nodes_at_loc.clear();

                get_rr_node_indices(device_ctx.rr_node_indices, sample_loc.x(), sample_loc.y(), rr_type, &rr_nodes_at_loc);
                for (int inode : rr_nodes_at_loc) {
                    if (inode < 0) continue;

                    RRNodeId node_id(inode);

                    int ptc = rr_graph.node_ptc_num(node_id);

                    if (ptc >= int(this->src_opin_delays[itile].size())) {
                        this->src_opin_delays[itile].resize(ptc + 1); //Inefficient but functional...
                    }

                    //Find the wire types which are reachable from inode and record them and
                    //the cost to reach them
                    util::dijkstra_flood_to_wires(itile, node_id, this->src_opin_delays);

                    if (this->src_opin_delays[itile][ptc].empty()) {
                        VTR_LOGV_DEBUG(f_router_debug, "Found no reachable wires from %s (%s) at (%d,%d)\n",
                                       rr_node_typename[rr_type],
                                       rr_node_arch_name(inode).c_str(),
                                       sample_loc.x(),
                                       sample_loc.y());

                        ptcs_with_no_delays = true;
                    }
                }

                ++num_sampled_locs;
            }
            if (ptcs_with_no_delays) {
                VPR_ERROR(VPR_ERROR_ROUTE, "Some SOURCE/OPINs have no reachable wires\n");
            }
        }
    }
}

void MapLookahead::compute_router_chan_ipin_lookahead() const {
    vtr::ScopedStartFinishTimer timer("Computing chan/ipin lookahead");
    auto& device_ctx = g_vpr_ctx.device();

    this->chan_ipins_delays.clear();

    this->chan_ipins_delays.resize(device_ctx.physical_tile_types.size());

    std::vector<int> rr_nodes_at_loc;

    //We assume that the routing connectivity of each instance of a physical tile is the same,
    //and so only measure one instance of each type
    for (auto tile_type : device_ctx.physical_tile_types) {
        vtr::Point<int> sample_loc(-1, -1);

        sample_loc = util::pick_sample_tile(&tile_type, sample_loc);

        if (sample_loc.x() == -1 && sample_loc.y() == -1) {
            //No untried instances of the current tile type left
            VTR_LOG_WARN("Found no sample locations for %s\n",
                         tile_type.name);
            continue;
        }

        int min_x = std::max(0, sample_loc.x() - X_OFFSET);
        int min_y = std::max(0, sample_loc.y() - Y_OFFSET);
        int max_x = std::min(int(device_ctx.grid.width()), sample_loc.x() + X_OFFSET);
        int max_y = std::min(int(device_ctx.grid.height()), sample_loc.y() + Y_OFFSET);

        for (int ix = min_x; ix < max_x; ix++) {
            for (int iy = min_y; iy < max_y; iy++) {
                for (auto rr_type : {CHANX, CHANY}) {
                    rr_nodes_at_loc.clear();

                    get_rr_node_indices(device_ctx.rr_node_indices, ix, iy, rr_type, &rr_nodes_at_loc);
                    for (int inode : rr_nodes_at_loc) {
                        if (inode < 0) continue;

                        RRNodeId node_id(inode);

                        //Find the IPINs which are reachable from the wires within the bounding box
                        //around the selected tile location
                        util::dijkstra_flood_to_ipins(node_id, this->chan_ipins_delays);
                    }
                }
            }
        }
    }
}

// get an expected minimum cost for routing from the current node to the target node
float MapLookahead::get_expected_cost(
    int current_node,
    int target_node,
    const t_conn_cost_params& params,
    float /*R_upstream*/) const {
    auto& device_ctx = g_vpr_ctx.device();

    t_rr_type rr_type = device_ctx.rr_nodes[current_node].type();

    if (rr_type == CHANX || rr_type == CHANY) {
        return get_map_cost(
            current_node, target_node, params.criticality);
    } else if (rr_type == IPIN) { /* Change if you're allowing route-throughs */
        return (device_ctx.rr_indexed_data[SINK_COST_INDEX].base_cost);
    } else { /* Change this if you want to investigate route-throughs */
        return (0.);
    }
}

#ifndef VTR_ENABLE_CAPNPROTO

void MapLookahead::read(const std::string& file) {
    VPR_THROW(VPR_ERROR_ROUTE, "MapLookahead::read not implemented");
}
void MapLookahead::write(const std::string& file) const {
    VPR_THROW(VPR_ERROR_ROUTE, "MapLookahead::write not implemented");
}

#else

void MapLookahead::read(const std::string& file) {
    cost_map_.read(file);

    this->compute_router_src_opin_lookahead();
    this->compute_router_chan_ipin_lookahead();
}

void MapLookahead::write(const std::string& file) const {
    cost_map_.write(file);
}

#endif
