// graph.h <Starter Code>
// < Hsien-Hao Chang >
//
// Basic graph class using adjacency matrix representation.  Currently
// limited to a graph with at most 100 vertices.
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
//

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <algorithm>

using namespace std;

template<typename VertexT, typename WeightT>
class graph {
 private:
	int numEdge;
	typedef unordered_map<VertexT, WeightT> vwMap;
	set<VertexT> vertices;
	unordered_map<VertexT, vwMap> adjList;

 public:
  //
  // constructor:
  //
  // Constructs an empty graph where n is the max # of vertices
  // you expect the graph to contain.
  //
  // NOTE: the graph is implemented using an adjacency matrix.
  // If n exceeds the dimensions of this matrix, an exception
  // will be thrown to let you know that this implementation
  // will not suffice.
  //
  // Deafualt Constructor
  graph() {
    numEdge = 0;
  }

  //
  // NumVertices
  //
  // Returns the # of vertices currently in the graph.
  //
  int NumVertices() const {
    return vertices.size();
  }

  //
  // NumEdges
  //
  // Returns the # of edges currently in the graph.
  //
  int NumEdges() const {
    return numEdge;
  }

  //
  // addVertex
  //
  // Adds the vertex v to the graph if there's room, and if so
  // returns true.  If the graph is full, or the vertex already
  // exists in the graph, then false is returned.
  //
  bool addVertex(VertexT v) {
    //
    // is the vertex already in the graph?  If so, we do not
    // insert again otherwise Vertices may fill with duplicates:
    //
    if (vertices.find(v) != vertices.end()) {
      return false;
    }

    //
    // if we get here, vertex does not exist so insert.  Where
    // we insert becomes the rows and col position for this
    // vertex in the adjacency matrix.
    //
    vertices.insert(v);
    return true;
  }

  //
  // addEdge
  //
  // Adds the edge (from, to, weight) to the graph, and returns
  // true.  If the vertices do not exist or for some reason the
  // graph is full, false is returned.
  //
  // NOTE: if the edge already exists, the existing edge weight
  // is overwritten with the new edge weight.
  //
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
    // If from and to don't exist in the graph, false.
    if (vertices.find(from) == vertices.end()
				|| vertices.find(to) == vertices.end()) {
      return false;
    }
    if (adjList[from].count(to) != 0) {
    // Do nothing...
    } else {
      numEdge++;
    }
    adjList[from][to] = weight;
    return true;
  }

  //
  // getWeight
  //
  // Returns the weight associated with a given edge.  If
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not
  // exist, the weight parameter is unchanged and false is
  // returned.
  //
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
    // check to, if not existed return false
    if (adjList.count(from) == 0) {  // fix seg fault checking from exists
      return false;
    }
    if (adjList.at(from).count(to) == 0) {
      return false;
    } else {
      weight = adjList.at(from).at(to);
      return true;
    }
  }

  //
  // neighbors
  //
  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  //
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT>  S;
    // big map, then small, all keys push into set. 2 for each loop
    if (adjList.count(v) == 0) {  // check if v is in adj or will seg fault
      return S;
    }
    for (auto p : adjList.at(v)) {
      S.insert(p.first);
    }
    return S;
  }

  //
  // getVertices
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  //
  vector<VertexT> getVertices() const {
    vector <VertexT> result;
    // set.... for loop   o(N) ?????
    for (auto p : vertices) {
      result.push_back(p);
    }
    return result;
  }

  //
  // dump
  //
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G(26);
  //    ...
  //    G.dump(cout);  // dump to console
  //
	void dump(ostream& output) const {
//     output << "***************************************************" << endl;
//     output << "********************* GRAPH ***********************" << endl;

//     output << "**Num vertices: " << this->NumVertices() << endl;
//     output << "**Num edges: " << this->NumEdges() << endl;

//     output << endl;
//     output << "**Vertices:" << endl;
//     for (int i = 0; i < this->NumVertices(); ++i) {
//       output << " " << i << ". " << this->Vertices[i] << endl;
//     }

//     output << endl;
//     output << "**Edges:" << endl;
//     for (int row = 0; row < this->NumVertices(); ++row) {
//       output << " row " << row << ": ";

//       for (int col = 0; col < this->NumVertices(); ++col) {
//         if (this->AdjMatrix[row][col].EdgeExists == false) {
//           output << "F ";
//         } else {
//           output << "(T,"
//             << this->AdjMatrix[row][col].Weight
//             << ") ";
//         }
//       }
//       output << endl;
//     }
//     output << "**************************************************" << endl;
	}
};
