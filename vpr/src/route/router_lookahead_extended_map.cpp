#include "router_lookahead_extended_map.h"

#include <vector>
#include <queue>

#include "rr_node.h"
#include "router_lookahead_map_utils.h"
#include "router_lookahead_sampling.h"
#include "globals.h"
#include "vtr_math.h"
#include "vtr_time.h"
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

#define LOOKAHEAD_EXTENDED_MAP_PRINT_COST_MAPS

// Don't continue storing a path after hitting a lower-or-same cost entry.
static constexpr bool BREAK_ON_MISS = false;

static constexpr int MIN_PATH_COUNT = 1000;

// compute the cost maps for lookahead
void ExtendedMapLookahead::compute(const std::vector<t_segment_inf>& segment_inf) {
    this->compute_router_src_opin_lookahead();
    this->compute_router_chan_ipin_lookahead();

    vtr::ScopedStartFinishTimer timer("Computing connection box lookahead map");

    // Initialize rr_node_route_inf if not already
    alloc_and_load_rr_node_route_structs();

    size_t num_segments = segment_inf.size();
    std::vector<SampleRegion> sample_regions = find_sample_regions(num_segments);

    /* free previous delay map and allocate new one */
    auto& device_ctx = g_vpr_ctx.device();
    cost_map_.set_counts(segment_inf.size());
    cost_map_.build_segment_map();

    VTR_ASSERT(REPRESENTATIVE_ENTRY_METHOD == util::SMALLEST);
    util::RoutingCosts all_delay_costs;
    util::RoutingCosts all_base_costs;

    /* run Dijkstra's algorithm for each segment type & channel type combination */
#if defined(VPR_USE_TBB)
    tbb::mutex all_costs_mutex;
    tbb::parallel_for_each(sample_regions, [&](const SampleRegion& region) {
#else
    for (const auto& region : sample_regions) {
#endif
        // holds the cost entries for a run
        util::RoutingCosts delay_costs;
        util::RoutingCosts base_costs;
        int total_path_count = 0;
        std::vector<bool> node_expanded(device_ctx.rr_nodes.size());
        std::vector<util::Search_Path> paths(device_ctx.rr_nodes.size());

        for (auto& point : region.points) {
            // statistics
            vtr::Timer run_timer;
            float max_delay_cost = 0.f;
            float max_base_cost = 0.f;
            int path_count = 0;
            for (auto node_ind : point.nodes) {
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
                VTR_LOG("Expanded %d paths of segment type %s(%d) starting at (%d, %d) from %d segments, max_cost %e %e (%g paths/sec)\n",
                        path_count, segment_inf[region.segment_type].name.c_str(), region.segment_type,
                        point.location.x(), point.location.y(),
                        (int)point.nodes.size(),
                        max_delay_cost, max_base_cost,
                        path_count / run_timer.elapsed_sec());
            }

            total_path_count += path_count;
            if (total_path_count > MIN_PATH_COUNT) {
                break;
            }
        }

#if defined(VPR_USE_TBB)
        all_costs_mutex.lock();
#endif

        if (total_path_count == 0) {
            VTR_LOG_WARN("No paths found for sample region %s(%d, %d)\n",
                         segment_inf[region.segment_type].name.c_str(), region.grid_location.x(), region.grid_location.y());
        }

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
#if defined(LOOKAHEAD_EXTENDED_MAP_PRINT_COST_ENTRIES)
    for (auto& cost : all_costs) {
        const auto& key = cost.first;
        const auto& val = cost.second;
        VTR_LOG("%d -> %d (%d, %d): %g, %g\n",
                val.from_node, val.to_node,
                key.delta.x(), key.delta.y(),
                val.cost_entry.delay, val.cost_entry.congestion);
    }
#endif

#if defined(LOOKAHEAD_EXTENDED_MAP_PRINT_COST_MAPS)
    for (int iseg = 0; iseg < (ssize_t)num_segments; iseg++) {
        VTR_LOG("cost map for %s(%d)\n",
                segment_inf[iseg].name.c_str(), iseg);
        cost_map_.print(iseg);
    }
#endif

#if defined(LOOKAHEAD_EXTENDED_MAP_PRINT_EMPTY_MAPS)
    for (std::pair<int, int> p : cost_map_.list_empty()) {
        int ichan, iseg;
        std::tie(ichan, iseg) = p;
        VTR_LOG("cost map for %s(%d), chan %d EMPTY\n",
                segment_inf[iseg].name.c_str(), iseg, box_id);
    }
#endif
}

std::pair<int, int> ExtendedMapLookahead::get_xy_deltas(const RRNodeId from_node, const RRNodeId to_node) const {
    auto& device_ctx = g_vpr_ctx.device();
    auto& rr_graph = device_ctx.rr_nodes;

    int from_x = vtr::nint((rr_graph.node_xlow(from_node) + rr_graph.node_xhigh(from_node)) / 2.);
    int from_y = vtr::nint((rr_graph.node_ylow(from_node) + rr_graph.node_yhigh(from_node)) / 2.);

    int to_x = vtr::nint((rr_graph.node_xlow(to_node) + rr_graph.node_xhigh(to_node)) / 2.);
    int to_y = vtr::nint((rr_graph.node_ylow(to_node) + rr_graph.node_yhigh(to_node)) / 2.);

    int dx, dy;
    dx = to_x - from_x;
    dy = to_y - from_y;

    return std::make_pair(dx, dy);
}
