#ifndef ROUTER_LOOKAHEAD_MAP_H_
#define ROUTER_LOOKAHEAD_MAP_H_

#include <vector>
#include "physical_types.h"
#include "router_lookahead.h"
#include "router_lookahead_map_utils.h"
#include "router_lookahead_cost_map.h"

// Implementation of RouterLookahead based on source segment and destination rr node types
class MapLookahead : public RouterLookahead {
  protected:
    //Look-up table from SOURCE/OPIN to CHANX/CHANY of various types
    mutable util::t_src_opin_delays src_opin_delays;

    //Look-up table from CHANX/CHANY to SINK/IPIN of various types
    mutable util::t_chan_ipins_delays chan_ipins_delays;

    //Runs the dijkstra algorithm and stores the paths found during the expansion
    template<typename Entry>
    std::pair<float, int> run_dijkstra(int start_node_ind,
                                       std::vector<bool>* node_expanded,
                                       std::vector<util::Search_Path>* paths,
                                       util::RoutingCosts* routing_costs) const;

    //Returns delay and congestion for the SOURCE/OPIN to CHAN connections.
    //It relies on a pre-computed data structure (t_src_opin_delays) for fast_lookups.
    std::pair<float, float> get_src_opin_delays(RRNodeId from_node, int delta_x, int delta_y, float criticality_fac) const;

    //Returns delay and congestion for the CHAN to SINK/IPIN connections.
    //It relies on a pre-computed data structure (t_chan_ipin_delays) for fast_lookups.
    float get_chan_ipin_delays(RRNodeId to_node) const;

    void compute_router_src_opin_lookahead() const;
    void compute_router_chan_ipin_lookahead() const;

    template<typename Entry>
    bool add_paths(int start_node_ind,
                   Entry current,
                   const std::vector<util::Search_Path>& paths,
                   util::RoutingCosts* routing_costs) const;

    virtual std::pair<int, int> get_xy_deltas(const RRNodeId from_node, const RRNodeId to_node) const;

  public:
    float get_expected_cost(int node, int target_node, const t_conn_cost_params& params, float R_upstream) const override;
    float get_map_cost(int from_node_ind, int to_node_ind, float criticality_fac) const;
    void compute(const std::vector<t_segment_inf>& segment_inf) override;

    void read(const std::string& file) override;
    void write(const std::string& file) const override;

    CostMap cost_map_;
};

#endif
