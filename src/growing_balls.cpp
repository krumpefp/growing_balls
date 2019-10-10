#include <algorithm>
#include <assert.h>
#include <iostream>
#include <string>
#include <vector>

#include "eliminationorder.h"
#include "io.h"
#include "pointofinterest.h"

int
main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << "Please provide a valid input file!"
              << std::endl;
    return 1;
  }

  std::string input_path(argv[1]);
  growing_balls::EliminationOrder eo;
  auto es = eo.compute_elimination_order(input_path);

  std::vector<growing_balls::PointOfInterest> elimination_order;
  for (auto& e : es) {
    auto et = e.m_elimination_time;
    auto poi = e.m_eliminated;
    auto elim_p = e.m_eliminated_by;

    poi.set_elimination(et, elim_p.get_osm_id());
    elimination_order.push_back(std::move(poi));
  }

  std::string output_path = input_path;
  if (input_path.find(".complete.txt") != std::string::npos) {
    // replace the old ending '.complete.txt' by .ce
    auto pos = output_path.rfind(".complete.txt");
    output_path.replace(pos, std::string::npos, ".ce");
  } else {
    output_path = output_path + ".ce";
  }

  growing_balls::IO::export_eliminationorder(output_path, elimination_order);
}
