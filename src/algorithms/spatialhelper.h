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
#include "textinput.h"

namespace growing_balls {

class SpatialHelper
{
public:
  using TextInput = growing_balls::TextInput;
  using OsmId = TextInput::OsmId;
  using POI = TextInput::PointOfInterest;

  using DataStorage = DataStorage<OsmId>;

  using Distance = growing_balls::Distance;

public:
  SpatialHelper(std::vector<POI>& pois);
  SpatialHelper(const SpatialHelper& other) = delete;
  ~SpatialHelper();
  SpatialHelper& operator=(const SpatialHelper& other);
  bool operator==(const SpatialHelper& other) const = delete;

  OsmId get_nearest_neighbor(OsmId id);

  std::vector<OsmId> get_in_range(OsmId id, Distance d);

private:
  std::unordered_map<OsmId, DataStorage::ElementId> m_id_mapper;
  DataStorage m_storage;
};

} // growing_balls

// BEGIN Declaration

namespace growing_balls {
using Distance = growing_balls::Distance;

using LabelElement = SpatialHelper::DataStorage::Element;
using LabelElementId = SpatialHelper::DataStorage::ElementId;
using ElementIdFactory = SpatialHelper::DataStorage::ElementIdFactory;

using OsmId = SpatialHelper::OsmId;

// BEGIN helpers and stuff

namespace {
struct ID_Mapper
{
  std::unordered_map<OsmId, LabelElementId> m_map;

  void operator()(LabelElement& elem)
  {
    m_map.emplace(elem.get_info(), elem.get_id());
  }
};

struct NN_Extractor
{
  Distance m_nn_distance = std::numeric_limits<Distance>::max();
  OsmId m_nn_id = ElementIdFactory::UNDEFINED_ID;
  LabelElement& m_query;

  NN_Extractor(OsmId query_id, SpatialHelper::DataStorage& storage)
    : m_query(storage.get(query_id)){};

  // CAUTION this code might not be thread save
  void operator()(LabelElement& elem)
  {
    Distance d = growing_balls::distance_in_centimeters(
      m_query.get_coord_1(), m_query.get_coord_2(), elem.get_coord_1(),
      elem.get_coord_2());
    if (d < m_nn_distance) {
      m_nn_distance = d;
      m_nn_id = elem.get_id();
    }
  };
};

struct RangeExplorer
{
  std::unordered_set<LabelElementId> m_neighbors;
  Distance m_distance_limit;
  LabelElement& m_query;
  LabelElementId m_current_id;
  SpatialHelper::DataStorage& m_storage;

  RangeExplorer(OsmId query_id, Distance max_dist,
                SpatialHelper::DataStorage& storage)
    : m_distance_limit(max_dist)
    , m_query(storage.get(query_id))
    , m_current_id(query_id)
    , m_storage(storage){};

  void operator()(LabelElement& elem)
  {
    LabelElementId elem_id = elem.get_id();
    if (m_neighbors.count(elem_id) != 0 || elem_id == m_query.get_id()) {
      return;
    }

    Distance d = growing_balls::distance_in_centimeters(
      m_query.get_coord_1(), m_query.get_coord_2(), elem.get_coord_1(),
      elem.get_coord_2());
    if (d <= m_distance_limit) {
      m_neighbors.insert(elem_id);

      m_current_id = elem_id;
      m_storage.visit_neighborhood(elem_id, *this);
    }
  }
};
}

// END helpers and stuff

// BEGIN class SpatialHelper
SpatialHelper::SpatialHelper(std::vector<POI>& pois)
{
  std::vector<LabelElement> lbl_elems;
  for (auto poi : pois) {
    lbl_elems.emplace_back(poi.get_lat(), poi.get_lon(), poi.get_osm_id());
  }

  m_storage.insert(lbl_elems.begin(), lbl_elems.end());

  ID_Mapper mapper;
  m_storage.visit_all(mapper);

  m_id_mapper = std::move(mapper.m_map);
}

SpatialHelper::~SpatialHelper()
{
  //   delete(m_storage);
}

SpatialHelper::OsmId
SpatialHelper::get_nearest_neighbor(OsmId id)
{
  NN_Extractor nn(m_id_mapper.at(id), m_storage);

  m_storage.visit_neighborhood(m_id_mapper.at(id), nn);

  return m_storage.get(nn.m_nn_id).get_info();
}

std::vector<OsmId>
SpatialHelper::get_in_range(OsmId id, Distance d)
{
  RangeExplorer r_exp(m_id_mapper.at(id), d, m_storage);

  m_storage.visit_neighborhood(m_id_mapper.at(id), r_exp);

  std::vector<OsmId> res;
  std::transform(r_exp.m_neighbors.begin(), r_exp.m_neighbors.end(),
                 std::back_inserter(res), [this](LabelElementId elem_id) {
                   return m_storage.get(elem_id).get_info();
                 });

  return res;
}

// END class SpatialHelper
}

// END Declaration

#endif // SPACIALHELPER_H
