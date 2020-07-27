#ifndef ROUTER_LOOKAHEAD_EXTENDED_MAP_H_
#define ROUTER_LOOKAHEAD_EXTENDED_MAP_H_

#include <vector>
#include "physical_types.h"
#include "router_lookahead_map.h"
#include "router_lookahead_map_utils.h"
#include "vtr_geometry.h"

// Implementation of RouterLookahead based on source segment and destination connection box types
class ExtendedMapLookahead : public MapLookahead {
  protected:
    std::pair<int, int> get_xy_deltas(const RRNodeId from_node, const RRNodeId to_node) const override;

  public:
    void compute(const std::vector<t_segment_inf>& segment_inf) override;
};

#endif
