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

#ifndef SPATIALHELPER_H
#define SPATIALHELPER_H

#include <limits>
#include <unordered_map>
#include <vector>

#include "datastorage.h"
#include "geofunctions.h"
#include "io.h"
#include "pointofinterest.h"

namespace growing_balls {

class SpatialHelper
{
public:
  using IO = growing_balls::IO;
  using POI = PointOfInterest;

  using DataStorageOsmId = DataStorage<OsmId>;

  using Distance = growing_balls::Distance;

  static const DataStorageOsmId::ElementId UNDEFINED_ID = DataStorageOsmId::ElementIdFactory::UNDEFINED_ID;

public:
  // Default constructor: initialize with an empty poi set
  SpatialHelper() = default;
  SpatialHelper(const std::vector<POI>& pois);
  SpatialHelper(SpatialHelper&& other) = delete;
  ~SpatialHelper();
  SpatialHelper& operator=(SpatialHelper&& other) = delete;

  OsmId get_nearest_neighbor(OsmId id) const;

  std::vector<OsmId> get_in_range(OsmId id, Distance d) const;

  void erase(OsmId);

  std::vector<OsmId> get_dupplicates() const { return m_dupplicates; };

private:
  std::unordered_map<OsmId, DataStorageOsmId::ElementId> m_id_mapper;
  std::vector<OsmId> m_dupplicates;
  DataStorageOsmId m_storage;
};

} // growing_balls

// BEGIN Declaration

namespace growing_balls {
using Distance = growing_balls::Distance;

using LabelElement = SpatialHelper::DataStorageOsmId::Element;
using LabelElementId = SpatialHelper::DataStorageOsmId::ElementId;
using ElementIdFactory = SpatialHelper::DataStorageOsmId::ElementIdFactory;

// BEGIN helpers and stuff

namespace {
struct ID_Mapper
{
  std::unordered_map<OsmId, LabelElementId> m_map;

  void operator()(const LabelElement& elem)
  {
    m_map.emplace(elem.get_info(), elem.get_id());
  }
};

struct NN_Extractor
{
  Distance m_nn_distance = std::numeric_limits<Distance>::max();
  // TODO: OsmId != UNDEFINED_ID!!!!
  OsmId m_nn_id = ElementIdFactory::UNDEFINED_ID;
  const LabelElement& m_query;

  NN_Extractor(OsmId query_id, const SpatialHelper::DataStorageOsmId& storage)
    : m_query(storage.get(query_id)){};

  // CAUTION this code might not be thread save
  void operator()(const LabelElement& elem)
  {
    Distance d = growing_balls::distance_in_centimeters(m_query.get_coord_1(),
                                                        m_query.get_coord_2(),
                                                        elem.get_coord_1(),
                                                        elem.get_coord_2());
    if (d < m_nn_distance) {
      m_nn_distance = d;
      m_nn_id = elem.get_info();
    }
  };
};

struct RangeExplorer
{
  std::unordered_set<OsmId> m_neighbors;
  Distance m_distance_limit;
  const LabelElement& m_query;
  LabelElementId m_current_id;
  const SpatialHelper::DataStorageOsmId& m_storage;

  RangeExplorer(OsmId query_id,
                Distance max_dist,
                const SpatialHelper::DataStorageOsmId& storage)
    : m_distance_limit(max_dist)
    , m_query(storage.get(query_id))
    , m_current_id(query_id)
    , m_storage(storage){};

  void operator()(const LabelElement& elem)
  {
    LabelElementId elem_id = elem.get_id();
    if (m_neighbors.count(elem.get_info()) != 0 || elem_id == m_query.get_id()) {
      return;
    }

    Distance d = growing_balls::distance_in_centimeters(m_query.get_coord_1(),
                                                        m_query.get_coord_2(),
                                                        elem.get_coord_1(),
                                                        elem.get_coord_2());
    if (d <= m_distance_limit) {
      m_neighbors.insert(elem.get_info());

      m_current_id = elem_id;
      m_storage.visit_neighborhood(elem_id, *this);
    }
  }
};
}

// END helpers and stuff

// BEGIN class SpatialHelper
SpatialHelper::SpatialHelper(const std::vector<POI>& pois)
{
  std::vector<LabelElement> lbl_elems;
  for (auto poi : pois) {
    lbl_elems.emplace_back(poi.get_lat(), poi.get_lon(), poi.get_pid());
  }

  m_dupplicates = m_storage.insert(lbl_elems.begin(), lbl_elems.end());

  ID_Mapper mapper;
  m_storage.visit_all(mapper);
  m_id_mapper = std::move(mapper.m_map);
}

SpatialHelper::~SpatialHelper()
{
  //   delete(m_storage);
}

OsmId
SpatialHelper::get_nearest_neighbor(OsmId id) const
{
  auto i_id = m_id_mapper.at(id);
  NN_Extractor nn(i_id, m_storage);

  m_storage.visit_neighborhood(i_id, nn);

  if (nn.m_nn_id == ElementIdFactory::UNDEFINED_ID) {
    return UNDEFINED_ID;
  }

  return nn.m_nn_id;
}

std::vector<OsmId>
SpatialHelper::get_in_range(OsmId id, Distance d) const
{
  auto i_id = m_id_mapper.at(id);
  RangeExplorer r_exp(i_id, d, m_storage);

  m_storage.visit_neighborhood(i_id, r_exp);

  std::vector<OsmId> res;
  res.insert(res.begin(), r_exp.m_neighbors.begin(), r_exp.m_neighbors.end());

  return res;
}

void
SpatialHelper::erase(OsmId id)
{
  m_storage.remove(m_id_mapper.at(id));
  m_id_mapper.erase(id);
}

// END class SpatialHelper
}

// END Declaration

#endif // SPACIALHELPER_H
