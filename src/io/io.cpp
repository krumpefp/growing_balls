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

#include "io.h"

#include <assert.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace growing_balls {

// BEGIN class IO

bool
IO::export_eliminationorder(std::string& export_file,
                            std::vector<PointOfInterest>& pois)
{
  std::ofstream outFile(export_file);
  if (!outFile) {
    std::cerr << "File " << export_file << " could not be opened!" << std::endl;
    std::cerr << "CAUTION! No output file has been generated!" << std::endl;
    return false;
  }

  // write the number of labels
  outFile << pois.size() << std::endl;

  // write the column header
  outFile << "lat lon osm_id priority collision_time radius font_factor label"
          << std::endl;

  outFile << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
  for (auto& poi : pois) {
    outFile << poi.get_lat() << " " << poi.get_lon() << " " << poi.get_osm_id()
            << ' ' << poi.get_priority() << " " << poi.get_elim_time() << " "
            << poi.get_radius() << " " << poi.get_font_factor() << " '"
            << poi.get_label() << "'" << std::endl;
  }
  return true;
}

std::vector<PointOfInterest>
IO::import_label(std::string input_file)
{
  std::vector<PointOfInterest> result;

  std::ifstream inFile(input_file);
  if (!inFile) {
    // if the file could not be opened: return an empty result
    std::cerr << "File " << input_file
              << " could not be opened to import label data" << std::endl;
    return result;
  }

  std::string header;
  std::getline(inFile, header);
  std::size_t total = std::stol(header);

  result.reserve(total);
  PoiId current_id = 0;
  for (std::string line; std::getline(inFile, line);) {
    if (line.size() > 0 && line.at(0) != '#') {
      result.emplace_back(current_id, line);
      ++current_id;
    }
  }

  assert(result.size() == total);

  return result;
}
}

// END class IO
