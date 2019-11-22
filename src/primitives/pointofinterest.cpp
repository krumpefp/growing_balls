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

#include "pointofinterest.h"

#include <sstream>
#include <string>

namespace growing_balls {
// BEGIN class PointOfInterest
PointOfInterest::PointOfInterest(PoiId pid, std::string input_str)
  : m_pid(pid)
  , m_elim_t(0)
  , m_elim_partner()
{
  std::stringstream ss(input_str);
  std::string first, second, last;
  std::getline(ss, first, '\'');
  std::getline(ss, second, '\'');
  std::getline(ss, last);

  m_label = second;
  std::stringstream ss_first(first);
  ss_first >> m_lat >> m_lon >> m_priority >> m_radius >> m_osmid;
  std::stringstream ss_last(last);
  ss_last >> m_font_fac;
}

std::string
PointOfInterest::print() const
{
  std::stringstream ss_result;

  ss_result.precision(22);
  ss_result << "Poi #" << m_pid << ": " << m_lat << " " << m_lon << " "
            << m_priority << " " << m_osmid << " '" << m_label << "' "
            << m_font_fac;

  return ss_result.str();
}

void
PointOfInterest::set_elimination(ElimTime elim_t, OsmId elim_p)
{
  m_elim_t = elim_t;
  m_elim_partner = elim_p;
}
// END class PointOfInterest
}
