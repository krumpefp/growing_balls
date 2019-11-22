/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2017  Filip Krumpe <filip.krumpe@fmi.uni-stuttgart.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef ELIMINATIONORDER_H
#define ELIMINATIONORDER_H

#include <algorithm>
#include <limits>
#include <optional>
#include <queue>
#include <stdint.h>
#include <string.h>
#include <unordered_map>

#include "geofunctions.h"
#include "io.h"
#include "pointofinterest.h"
#include "spatialhelper.h"
#include "timer.h"

namespace {
bool
prefer(const growing_balls::PointOfInterest& p1,
       const growing_balls::PointOfInterest& p2)
{
  return p1.get_priority() > p2.get_priority();
};
}

namespace growing_balls {
using Time = ElimTime;

class EliminationOrder
{
public:
  struct Elimination
  {
    PointOfInterest m_eliminated;
    PointOfInterest m_eliminated_by;
    Time m_elimination_time;

    Elimination(Time elim_t, PointOfInterest elim, PointOfInterest elim_by)
      : m_eliminated(std::move(elim))
      , m_eliminated_by(std::move(elim_by))
      , m_elimination_time(elim_t){};
  };

public:
  EliminationOrder() = default;

  EliminationOrder(EliminationOrder&& other) = default;
  EliminationOrder& operator=(EliminationOrder&& other) = default;

  std::vector<Elimination> compute_elimination_order(std::string file);

private:
  std::vector<double> m_debug_times;
};
}

#endif // ELIMINATIONORDER_H

// BEGIN DEFINITION

namespace growing_balls {

// BEGIN Helpers
namespace {
using PoiList = std::vector<PointOfInterest>;

enum class EventType : int32_t
{
  UPDATE_EVENT = 1,
  COLLISION_EVENT = 2,
};

struct Event
{
  Time m_trigger_time;
  PoiId m_coll1;
  PoiId m_coll2;

  EventType m_evt_type;

  Event(Time t, PoiId c1)
    : m_trigger_time(t)
    , m_coll1(c1)
    , m_coll2(c1)
    , m_evt_type(EventType::UPDATE_EVENT){};

  Event(Time t, PoiId c1, PoiId c2)
    : m_trigger_time(t)
    , m_coll1(c1)
    , m_coll2(c2)
    , m_evt_type(EventType::COLLISION_EVENT){};

  Event(const Event&) = default;

  // In order to use a priority_queue (a max heap) as min heap: overload
  // operator <
  bool operator<(const Event& other) const
  {
    if (m_trigger_time == other.m_trigger_time) {
      if (m_evt_type == EventType::UPDATE_EVENT) {
        // update events are prefered to collision events
        return false;
      } else if (m_evt_type == other.m_evt_type) {
        // make things deterministic: use m_coll1 ids to break ties with equal
        // event type
        return m_coll1 > other.m_coll1;
      } else {
        // other is an update event and we are not -> other has higher priority
        return true;
      }
    } else {
      return m_trigger_time > other.m_trigger_time;
    }
  }
};

Time
compute_collision_time(const PointOfInterest& p1, const PointOfInterest& p2)
{
  auto d = distance_in_centimeters(
    p1.get_lat(), p1.get_lon(), p2.get_lat(), p2.get_lon());

  return d / (p1.get_radius() + p2.get_radius());
}

std::optional<Event>
predict_collision(const PointOfInterest& p,
                  Time t,
                  SpatialHelper& sh,
                  const PoiList& poi_list)
{
  auto id_nn = sh.get_nearest_neighbor(p.get_pid());

  if (id_nn == sh.UNDEFINED_ID) {
    return Event(
      std::numeric_limits<Time>::max(), p.get_pid(), p.get_pid());
  }

  auto& nn = poi_list.at(id_nn);

  auto distance_pnn = distance_in_centimeters(
    p.get_lat(), p.get_lon(), nn.get_lat(), nn.get_lon());
  Time t_upd = distance_pnn / (2 * p.get_radius());

  if (t < t_upd) {
    return Event(t_upd, p.get_pid());
  } else {
    auto p_rad = p.get_radius();

    Time min_coll_t = std::numeric_limits<Time>::max();
    OsmId min_coll_p = 0;
    for (auto id : sh.get_in_range(p.get_pid(), 2 * distance_pnn)) {
      auto& p2 = poi_list.at(id);
      if (p2.get_radius() > p_rad) {
        // p2 is responsible for predicting the collision
        continue;
      }

      Time coll_t = compute_collision_time(p, p2);
      if (coll_t < min_coll_t) {
        min_coll_t = coll_t;
        min_coll_p = p2.get_pid();
      }
    }

    if (min_coll_t != std::numeric_limits<Time>::max()) {
      return Event(min_coll_t, p.get_pid(), min_coll_p);
    } else {
      // return nullopt if no event to a more or equally prioritized point was
      // found
      return {};
    }
  }
}
}
// END Helpers

// BEGIN class EliminationOrder

std::vector<EliminationOrder::Elimination>
EliminationOrder::compute_elimination_order(std::string file)
{
  debug_timer::Timer timer;
  timer.start();
  PoiList pois = IO::import_label(file);
  std::vector<bool> alive;
  alive.resize(pois.size(), true);

  timer.createTimepoint();

  SpatialHelper spatial_helper(pois);

  timer.createTimepoint();

  // remove dupplicated pois
  for (auto& id : spatial_helper.get_dupplicates()) {
    auto dupplicate = pois.at(id);
    std::cout << "Ignoring dupplicated poi: " << dupplicate.print()
              << std::endl;
    dupplicate.set_elimination(0, dupplicate.get_osm_id());
    alive.at(id) = false;
  }

  std::vector<Elimination> result;

  std::priority_queue<Event> Q;

  timer.createTimepoint();
  // initialize
  for (const auto& p : pois) {
    if (!alive.at(p.get_pid())) continue;
    if (auto evt = predict_collision(p, 0., spatial_helper, pois)) {
      assert(evt->m_evt_type == EventType::UPDATE_EVENT);
      Q.push(*evt);
    }
  }

  timer.createTimepoint();
  Time t = 0;
  while (!Q.empty()) {
    auto current_evt = Q.top();
    Q.pop();
    t = current_evt.m_trigger_time;

    if (alive.at(current_evt.m_coll1)) {
      auto p1 = pois.at(current_evt.m_coll1);
      if (current_evt.m_evt_type == EventType::UPDATE_EVENT) {
        if (auto evt = predict_collision(p1, t, spatial_helper, pois))
          Q.push(*evt);
      } else if (current_evt.m_evt_type == EventType::COLLISION_EVENT) {
        if (current_evt.m_coll1 == current_evt.m_coll2) {
          result.emplace_back(t, p1, p1);
          break;
        }

        if (alive.at(current_evt.m_coll2)) {
          auto p2 = pois.at(current_evt.m_coll2);
          // here p1 and p2 are alive
          if (prefer(p1, p2)) {
            result.emplace_back(t, p2, p1);
            spatial_helper.erase(p2.get_pid());
            alive.at(p2.get_pid()) = false;

            if (auto evt =
                  predict_collision(p1, t, spatial_helper, pois))
              Q.push(*evt);
          } else {
            result.emplace_back(t, p1, p2);
            spatial_helper.erase(p1.get_pid());
            alive.at(p1.get_pid()) = false;

            if (auto evt =
                  predict_collision(p2, t, spatial_helper, pois))
              Q.push(*evt);
          }
        } else {
          // p1 is the only collision partner alive => repredict collision
          if (auto evt = predict_collision(p1, t, spatial_helper, pois))
            Q.push(*evt);
        }
      }
    } else {
      // check if p2 is alive and predict its next collision
      // otherwise do nothing
      if (alive.at(current_evt.m_coll2)) {
        auto p2 = pois.at(current_evt.m_coll2);
        if (auto evt = predict_collision(p2, t, spatial_helper, pois))
          Q.push(*evt);
      }
    }
  }

  timer.stop();

  auto times = timer.getTimes();

  std::cout << "Computation of the elimination order finished." << std::endl;
  std::cout << "Computed " << pois.size() << " many eliminations." << std::endl;
  std::cout << "Required times in ms were as follows:" << std::endl;
  std::cout << times[0] << "\t Label import" << std::endl;
  std::cout << times[1] << "\t Initialization of the spatial helper"
            << std::endl;
  std::cout << times[3] << "\t Initialization of the algorithm" << std::endl;
  std::cout << times[4] << "\t Main algorithm loop" << std::endl;

  return result;
}

// END class EliminationOrder
}
// END DEFINITION
