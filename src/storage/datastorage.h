#ifndef DATASTORAGE_H
#define DATASTORAGE_H

#include <assert.h>
#include <cstdio>
#include <limits>
#include <list>
#include <stdint.h>
#include <unordered_set>
#include <vector>

#include <CGAL/spatial_sort.h>
#include <libdts2/Constrained_delaunay_triangulation_s2.h>

#include "geofunctions.h"

namespace {
struct Helpers;
}

namespace growing_balls {

template<typename Info>
class DataStorage
{
public:
  using Dimension_1 = double;
  using Dimension_2 = double;
  using LessDim1 = std::less<Dimension_1>;
  using LessDim2 = std::less<Dimension_2>;

  using NodeID = std::size_t;
  using Priority = uint32_t;
  using Radius = double;

  using InfoType = Info;

public:
  using ElementId = uint32_t;

protected:
  using CDT = dts2::Delaunay_triangulation_with_info_s2<ElementId, void>;

  using FaceHandle = CDT::Face_handle;
  using LocateType = CDT::Locate_type;

  using Point3 = CDT::Point_3;

  using Vertex = CDT::Vertex;
  using VertexHandle = CDT::Vertex_handle;
  using VertexCirculator = CDT::Vertex_circulator;

  using ProjectOnSphere = CDT::Project_on_sphere;

public:
  class Element
  {
    Dimension_1 m_dim1;
    Dimension_2 m_dim2;

    InfoType m_info;

    ElementId m_id = std::numeric_limits<uint32_t>::max();
    VertexHandle m_cdt_handle = nullptr;

  public:
    Element(Dimension_1 dim1, Dimension_2 dim2, InfoType info);

    Element(Element&& other);
    Element& operator=(Element&& other);

    void set_handle(VertexHandle handle) { m_cdt_handle = handle; };
    void set_id(ElementId id) { m_id = id; };

    // getters
    Dimension_1 get_coord_1() const { return m_dim1; };
    Dimension_2 get_coord_2() const { return m_dim2; };
    InfoType get_info() const { return m_info; };

    ElementId get_id() const { return m_id; };

  protected:
    VertexHandle get_handle() const { return m_cdt_handle; };

    friend class DataStorage;
  };

  /**
   * Ids 0-5 are reserved:
   * max    : UNDEFINED_ID
   * 1 - 4  : auxiliary points
   *
   * further ids map to the vector like follows:
   * 5 -> 0   - first element
   * i -> i-5 - for i > 4
   */
  struct ElementIdFactory
  {
    const static ElementId UNDEFINED_ID = std::numeric_limits<ElementId>::max();

    static std::size_t get_vpos_from_id(ElementId id) { return id - 5; };

    static ElementId get_next_id(const std::vector<Element>& data)
    {
      assert(data.size() < std::numeric_limits<uint32_t>::max());
      return data.size() + 5;
    }

    static bool is_undefined(ElementId id) { return id == UNDEFINED_ID; };
  };

public:
  DataStorage();

  /**
   * Insert the emenent into the data storage.
   *
   * If there exists an element with the same position, the new element will not
   * be inserted and the "undefined id" is returned.
   * If the element was inserted successfully, its internal id is returned.
   */
  ElementId insert(Element&& elem);

  /**
   * Insert the elements in the range begin to end.
   *
   * Elements with a position equals to an already inserted element are skipped.
   * Skipped elements are returned in the returned vector,
   */
  std::vector<InfoType> insert(typename std::vector<Element>::iterator begin,
                               typename std::vector<Element>::iterator end);

  std::size_t remove(ElementId id);

  Element& get(ElementId id)
  {
    return m_elements.at(ElementIdFactory::get_vpos_from_id(id));
  };

  template<typename Visitor>
  void visit_all(Visitor& v);

  /**
   * The neighborhood visitor is applied for each neighbor of the specified
   * element.
   *
   * It is guaranteed that one of these neighbors is the nearest neighbor.
   */
  template<typename Visitor>
  void visit_neighborhood_with_aux(const VertexCirculator& begin,
                                   const VertexCirculator& current,
                                   const VertexCirculator& end,
                                   ElementId query,
                                   Visitor& v);

  /**
   * The neighborhood visitor is applied for each neighbor of the specified
   * element.
   *
   * It is guaranteed that one of these neighbors is the nearest neighbor.
   */
  template<typename Visitor>
  void visit_neighborhood(ElementId elem_id, Visitor& v);

  friend Helpers;

private:
  std::vector<Element> m_elements;

  CDT m_cdt;

  ProjectOnSphere m_proj;
};
}

// BEGIN Declaration

// BEGIN Helpers
namespace growing_balls {
#define TMPL_HDR template<typename Info>
#define TMPL_CLS DataStorage<Info>

TMPL_HDR
struct Helpers
{
  using Element = typename TMPL_CLS::Element;

  struct SpatialSortingTrait
  {
    using Point_2 = Element;
    using Less_x_2 = std::function<bool(const Element&, const Element&)>;
    using Less_y_2 = std::function<bool(const Element&, const Element&)>;

    Less_x_2 less_x_2_object() const
    {
      auto cmp_d1 = typename TMPL_CLS::LessDim1();
      return [cmp_d1](const Element& a, const Element& b) {
        return cmp_d1(a.get_coord_1(), b.get_coord_1());
      };
    }

    Less_y_2 less_y_2_object() const
    {
      auto cmp_d2 = typename TMPL_CLS::LessDim2();
      return [cmp_d2](const Element& a, const Element& b) {
        return cmp_d2(a.get_coord_2(), b.get_coord_2());
      };
    }
  };

  static void spatial_sort(typename std::vector<Element>::iterator begin,
                           typename std::vector<Element>::iterator end)
  {
    SpatialSortingTrait sst;
    CGAL::spatial_sort(begin, end, sst);
  }
};
// END Helpers

// BEGIN DataStorage
TMPL_HDR
TMPL_CLS::DataStorage()
  : m_cdt(64)
  , m_proj(m_cdt.geom_traits().project_on_sphere_object()){};

TMPL_HDR typename TMPL_CLS::ElementId
TMPL_CLS::insert(typename TMPL_CLS::Element&& elem)
{
  LocateType lt;
  int32_t li;

  Point3 pos = m_cdt.project(elem.get_coord_1(), elem.get_coord_2());

  FaceHandle fh = m_cdt.locate(pos, lt, li);
  if (lt == LocateType::VERTEX) {
    // skip insertion: point is already contained
    return ElementIdFactory::UNDEFINED_ID;
  }
  auto vh = m_cdt.insert(pos, fh);

  elem.set_id(ElementIdFactory::get_next_id(m_elements));
  elem.set_handle(vh);
  vh->info() = elem.get_id();

  assert(m_elements.size() ==
         ElementIdFactory::get_vpos_from_id(elem.get_id()));
  m_elements.push_back(std::move(elem));

  return m_elements.back().get_id();
}

TMPL_HDR std::vector<typename TMPL_CLS::InfoType>
TMPL_CLS::insert(typename std::vector<Element>::iterator begin,
                 typename std::vector<Element>::iterator end)
{
  std::vector<typename TMPL_CLS::InfoType> result;

  // ATTENTION: reserving the storage space in advance may lead to performance
  // issues if many smaller chunks of elements are inserted into the data
  // structure!
  m_elements.reserve(m_elements.size() + std::distance(begin, end));

  Helpers<InfoType>::spatial_sort(begin, end);

  FaceHandle fh;
  LocateType lt;
  int32_t li;

  for (auto it = begin; it != end; ++it) {
    Point3 pos = m_cdt.project(it->get_coord_1(), it->get_coord_2());

    fh = m_cdt.locate(pos, lt, li, fh);
    if (lt == LocateType::VERTEX) {
      // skip insertion: point is already contained
      result.push_back(it->get_info());
      continue;
    }
    auto vh = m_cdt.insert(pos, fh);

    it->set_id(ElementIdFactory::get_next_id(m_elements));
    it->set_handle(vh);
    vh->info() = it->get_id();

    assert(m_elements.size() ==
           ElementIdFactory::get_vpos_from_id(it->get_id()));

    m_elements.push_back(std::move(*it));
  }

  return result;
}

TMPL_HDR template<typename Visitor>
void
TMPL_CLS::visit_all(Visitor& v)
{
  for (auto& elem : m_elements) {
    v(elem);
  }
}

TMPL_HDR template<typename Visitor>
void
TMPL_CLS::visit_neighborhood_with_aux(const VertexCirculator& begin,
                                      const VertexCirculator& current,
                                      const VertexCirculator& end,
                                      ElementId query,
                                      Visitor& v)
{
  VertexCirculator it = begin;

  // add the elements that were already visited to skip
  std::unordered_set<ElementId> skip;
  skip.insert(query);
  while (it != current) {
    assert(!m_cdt.is_auxiliary(it));
    skip.insert(it->info());

    ++it;
  }

  // now it == current and should be the first aux point in the Circulator
  assert(it == current && m_cdt.is_auxiliary(it));
  do {
    if (!m_cdt.is_auxiliary(it)) {
      ElementId id = it->info();
      skip.insert(id);
      v(m_elements.at(ElementIdFactory::get_vpos_from_id(id)));
    } else {
      // recursively cycle through the neighbors of the aux point to avoid
      // breaking the guarantee that the real nearest neighbor is visited
      VertexCirculator it2 = m_cdt.incident_vertices(it);
      auto end2 = it2;
      do {
        if (!m_cdt.is_auxiliary(it2)) {
          ElementId id2 = it2->info();
          if (skip.find(id2) == skip.end()) {
            skip.insert(id2);
            v(m_elements.at(ElementIdFactory::get_vpos_from_id(id2)));
          }
        } else {
          // ... and once again to ensure we jumped over the triangle at the
          // north pole
          VertexCirculator it3 = m_cdt.incident_vertices(it2);
          auto end3 = it3;
          do {
            if (!m_cdt.is_auxiliary(it3)) {
              ElementId id3 = it3->info();
              skip.insert(id3);
              if (skip.find(id3) == skip.end()) {
                v(m_elements.at(ElementIdFactory::get_vpos_from_id(id3)));
              }
            }
            ++it3;
          } while (it3 != end3);
        }
        ++it2;
      } while (it2 != end2);
    }

    ++it;
  } while (it != end);
}

TMPL_HDR template<typename Visitor>
void
TMPL_CLS::visit_neighborhood(TMPL_CLS::ElementId elem_id, Visitor& v)
{
  auto& elem = m_elements.at(ElementIdFactory::get_vpos_from_id(elem_id));

  auto elem_hdl = elem.get_handle();

  VertexCirculator begin = m_cdt.incident_vertices(elem_hdl);
  VertexCirculator vc = begin;
  auto end = vc;
  do {
    if (m_cdt.is_auxiliary(vc)) {
      visit_neighborhood_with_aux(begin, vc, end, elem_id, v);
      break;
    }

    ElementId id = vc->info();
    v(m_elements.at(ElementIdFactory::get_vpos_from_id(id)));

    ++vc;
  } while (vc != end);
}

TMPL_HDR std::size_t
TMPL_CLS::remove(ElementId id)
{
  ElementId i_id = ElementIdFactory::get_vpos_from_id(id);
  assert(i_id < m_elements.size());

  VertexHandle vh = m_elements.at(i_id).get_handle();
  if (vh == nullptr) {
    assert(m_elements.at(i_id).get_id() == ElementIdFactory::UNDEFINED_ID);
    return 0;
  }
  m_cdt.remove(vh);

  auto& elem = m_elements.at(i_id);
  elem.set_handle(nullptr);
  elem.set_id(ElementIdFactory::UNDEFINED_ID);

  return 1;
}
// END DataStorage

// BEGIN DataStorage::Element
TMPL_HDR
TMPL_CLS::Element::Element(Dimension_1 dim1, Dimension_2 dim2, InfoType info)
  : m_dim1(dim1)
  , m_dim2(dim2)
  , m_info(info)
  , m_id(ElementIdFactory::UNDEFINED_ID){};

TMPL_HDR
TMPL_CLS::Element::Element(Element&& other)
  : m_dim1(std::move(other.m_dim1))
  , m_dim2(std::move(other.m_dim2))
  , m_info(std::move(other.m_info))
  , m_id(std::move(other.m_id))
  , m_cdt_handle(std::move(other.m_cdt_handle)){};

TMPL_HDR typename TMPL_CLS::Element&
TMPL_CLS::Element::operator=(Element&& other)
{
  m_dim1 = std::move(other.m_dim1);
  m_dim2 = std::move(other.m_dim2);
  m_info = std::move(other.m_info);
  m_id = std::move(other.m_id);
  m_cdt_handle = std::move(other.m_cdt_handle);

  return *this;
};
// END Datastorage::Element
}
// END Declaration
#endif // DATASTORAGE_H
