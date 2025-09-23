/*
 * This file belongs to the Galois project, a C++ library for exploiting
 * parallelism. The code is being released under the terms of the 3-Clause BSD
 * License (a copy is located in LICENSE.txt at the top-level directory).
 *
 * Copyright (C) 2018, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 */

#include <iostream>
#include <limits>
#include "DistBench/Start.h"
#include "galois/DistGalois.h"
#include "galois/gstl.h"
#include <fstream>

/******************************************************************************/
/* Declaration of command line arguments */
/******************************************************************************/

namespace cll = llvm::cl;

/******************************************************************************/
/* Graph structure declarations + other initialization */
/******************************************************************************/

struct NodeData {
  uint32_t dummy;
};

typedef galois::graphs::DistGraph<NodeData, void> Graph;
typedef typename Graph::GraphNode GNode;
std::unique_ptr<galois::graphs::GluonSubstrate<Graph>> syncSubstrate;

/******************************************************************************/
/* Main */
/******************************************************************************/

constexpr static const char* const name = "Partition";
constexpr static const char* const desc = "Partitions a normal graph.";
constexpr static const char* const url  = 0;

static cll::opt<unsigned int>
    numParts("numParts", cll::desc("Number of parts to partition into"),
             cll::init(3));

static cll::opt<std::string>
    output_path("output_path", cll::desc("Output file for partitioned graph"),
           cll::init(""));

int main(int argc, char** argv) {
  galois::DistMemSys G;
  DistBenchStart(argc, argv, name, desc, url);

  std::unique_ptr<Graph> hg;
  std::tie(hg, syncSubstrate) = distGraphInitialization<NodeData, void>();

  auto& net = galois::runtime::getSystemNetworkInterface();

  if (net.ID == 0) {

    std::ofstream outfile;
    if (output_path != "") {
      outfile.open(output_path);
    } else {
      outfile.open("part_" + std::to_string(numParts) + "_" +
                   std::to_string(hg->numGlobalNodes) + "_" +
                   std::to_string(hg->numGlobalEdges) + ".txt");
    }

    for (uint64_t i = 0; i < hg->numGlobalNodes; ++i) {
      auto part = hg->getHostID(i);
      outfile << part << "\n";
    }

    outfile.close();
  }

  return 0;
}
