// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


Point<3>
trans_func(const Point<3> &p)
{
  Point<3> r(p[0] + p[1] * p[1], p[1], p[2]);
  return r;
}



void
test()
{
  Triangulation<3> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.execute_coarsening_and_refinement();


  typename Triangulation<3, 3>::active_cell_iterator cell = tria.begin_active();
  for (unsigned int i = 0; i < 3; i++)
    {
      cell->set_refine_flag(RefinementCase<3>::cut_axis(i));
      cell++;
    }
  tria.execute_coarsening_and_refinement();

  deallog << "Unchanged grid:" << std::endl;
  GridOut().write_gnuplot(tria, deallog.get_file_stream());
  {
    std::ofstream f("grid1");
    GridOut().write_gnuplot(tria, f);
  }

  GridTools::transform(trans_func<3>, tria);
  deallog << "transformed grid:" << std::endl;
  GridOut().write_gnuplot(tria, deallog.get_file_stream());
  {
    std::ofstream f("grid2");
    GridOut().write_gnuplot(tria, f);
  }
}


int
main()
{
  initlog();

  test();

  return 0;
}
