/*
 * Represents a geographic Point of Interest with additional labeling
 * information.
 *
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

#ifndef POINTOFINTEREST_H
#define POINTOFINTEREST_H

#include <stdint.h>
#include <string>
#include <vector>

namespace growing_balls {
using Coord_Lat = double;
using Coord_Lon = double;

using FontFactor = double;
using Label = std::string;
using OsmId = int64_t;
using PoiId = uint64_t;
using Priority = uint32_t;
using Radius = double;

using ElimTime = double;

class PointOfInterest
{
public:
  PointOfInterest(PoiId pid,
                  Coord_Lat lat,
                  Coord_Lon lon,
                  FontFactor f_fac,
                  Label lbl,
                  OsmId oid,
                  Priority prio,
                  Radius rad)
    : m_pid(pid)
    , m_lat(lat)
    , m_lon(lon)
    , m_font_fac(f_fac)
    , m_label(lbl)
    , m_osmid(oid)
    , m_priority(prio)
    , m_radius(rad){};

  /**
   * parses a string of the following form and initializes the PointOfInterest
   * <lat> <lon> <priority> <radius> <osm_id> '<label>' <font_factor>
   */
  PointOfInterest(PoiId pid, std::string input_str);

  ElimTime get_elim_time() const { return m_elim_t; };
  FontFactor get_font_factor() const { return m_font_fac; };
  Label get_label() const { return m_label; };
  Coord_Lat get_lat() const { return m_lat; };
  Coord_Lon get_lon() const { return m_lon; };
  OsmId get_osm_id() const { return m_osmid; };
  Priority get_priority() const { return m_priority; };
  Radius get_radius() const { return m_radius; };

  std::string print() const;

  void set_elimination(ElimTime et, OsmId e_partner);

private:
  PoiId m_pid;

  // position
  Coord_Lat m_lat;
  Coord_Lon m_lon;

  // label related stuff
  FontFactor m_font_fac;
  Label m_label;
  OsmId m_osmid;
  Priority m_priority;
  Radius m_radius;

  ElimTime m_elim_t;
  OsmId m_elim_partner;
};
}

#endif // POINTOFINTEREST_H
