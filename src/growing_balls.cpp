#include <algorithm>
#include <assert.h>
#include <iostream>
#include <string>
#include <vector>

#include "eliminationorder.h"
#include "io.h"
#include "pointofinterest.h"

growing_balls::EliminationOrder::Heuristic parse_heur(std::string& sz_heur) {
  if (sz_heur == "RAD_INV") return growing_balls::EliminationOrder::Heuristic::RAD_INV;
  if (sz_heur == "ID") return growing_balls::EliminationOrder::Heuristic::ID;
  if (sz_heur == "ID_INV") return growing_balls::EliminationOrder::Heuristic::ID_INV;
  if (sz_heur == "NXE2") return growing_balls::EliminationOrder::Heuristic::NXE2;
  if (sz_heur == "NXE2_INV") return growing_balls::EliminationOrder::Heuristic::NXE2_INV;
  if (sz_heur == "NXE4") return growing_balls::EliminationOrder::Heuristic::NXE4;
  if (sz_heur == "NXE4_INV") return growing_balls::EliminationOrder::Heuristic::NXE4_INV;
  if (sz_heur == "NXE8") return growing_balls::EliminationOrder::Heuristic::NXE8;
  if (sz_heur == "NXE8_INV") return growing_balls::EliminationOrder::Heuristic::NXE8_INV;
  if (sz_heur == "NXE16") return growing_balls::EliminationOrder::Heuristic::NXE16;
  if (sz_heur == "NXE16_INV") return growing_balls::EliminationOrder::Heuristic::NXE16_INV;
  std::cerr << "Could not parse heuristic " << sz_heur << ".\nFalling back to RAD" << std::endl;
  return growing_balls::EliminationOrder::Heuristic::RAD;
}

int
main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << "Please provide a valid input file!"
              << std::endl;
    return 1;
  }
  growing_balls::EliminationOrder::Heuristic heuristic = growing_balls::EliminationOrder::Heuristic::RAD;
  if (argc >= 3) {
    std::string sz_heur(argv[2]);
    heuristic = parse_heur(sz_heur);
    std::cout << "Using heuristic " << sz_heur << std::endl;
  }

  std::string input_path(argv[1]);
  growing_balls::EliminationOrder eo(heuristic);
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
