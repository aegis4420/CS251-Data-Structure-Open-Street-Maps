// application.cpp <Starter Code>
// <Hsien-Hao Chang>
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
//
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <queue>

#include "tinyxml2.h"
#include "dist.h"
#include "osm.h"
#include "graph.h"

using namespace std;
using namespace tinyxml2;

const double INF = numeric_limits<double>::max();
//
// Implement your creative component application here
// TO DO: add arguments
//
class prioritize {
	public:
		bool operator()(const pair<long long, double> &p1,
										const pair<long long, double> &p2) const{
			return p1.second > p2.second;
		}
};
void creative() {
}

bool searchBuilding(string query,
										vector<BuildingInfo> &Buildings,
										BuildingInfo &building) {
	// Check for person1Building exists in the vector
	for (unsigned int i = 0; i < Buildings.size(); i++) {
			if (query == Buildings.at(i).Abbrev) {
				building = Buildings.at(i);
				return true;
			}
			// If Abbrevs don't match, look for Fullname..
			if (i == Buildings.size() - 1) {
				for (unsigned int j = 0; j < Buildings.size(); j++) {
					if (Buildings.at(j).Fullname.find(query) != string::npos) {
						building = Buildings.at(j);
						return true;
					}
				}
			}
		}
	return false;
}

BuildingInfo nearestBuilding(vector<BuildingInfo> &Buildings,
														 Coordinates &midpoint,
														 set <long long> visited) {
	double min = INF;
	double distance;
	int minIndex;
	for (unsigned int i = 0; i < Buildings.size(); i++) {
		if (visited.count(Buildings[i].Coords.ID) != 0) {
			continue;
		}
		distance = distBetween2Points(midpoint.Lat, midpoint.Lon,
																	Buildings.at(i).Coords.Lat,
																	Buildings.at(i).Coords.Lon);
		if (distance < min) {
			min = distance;
			minIndex = i;  // Recording the min's index.
		}
	}
	return Buildings.at(minIndex);
}

long long nearestNode(map<long long, Coordinates> &Nodes,
											vector<FootwayInfo> &Footways,
											BuildingInfo building) {
	double min = INF;
	double distance;
	long long minID;
	for (auto i : Footways) {
		for (unsigned int e = 0; e < i.Nodes.size(); e++) {
			distance = distBetween2Points(Nodes[i.Nodes.at(e)].Lat,
																		Nodes[i.Nodes.at(e)].Lon,
																		building.Coords.Lat,
																		building.Coords.Lon);
			if (distance < min) {
				min = distance;
				minID = i.Nodes.at(e);  // Recording the min node's index
			}
		}
  }
	return minID;
}

map <long long, long long> Dijkstra(graph <long long, double> &G,
																		long long start,
																		long long dest,
																		map <long long, double> &distance) {
	map <long long, long long> pred;  // will be return for the path
	vector<long long>  visited;
  priority_queue<pair<long long, double>,
	vector<pair<long long, double>>, prioritize> pq;
  vector<long long> allNodes = G.getVertices();
	for (long long vertex : allNodes) {  // Change 1
    distance[vertex] = INF;
    pq.push(make_pair(vertex, INF));
  }
  pq.push(make_pair(start, 0));
  distance[start] = 0;
  while (!pq.empty()) {
    pair<long long, double> thisNode = pq.top();
    pq.pop();
		// If we reach destination, break
    if (thisNode.second == INF || thisNode.first == dest) {
      break;
    }
    bool found = false;
    for (long long node : visited) {
      if (thisNode.first == node) {
        found = true;
      }
    }
    if (found) {
      continue;
    }
    visited.push_back(thisNode.first);
    set<long long> neighbors = G.neighbors(thisNode.first);
    for (long long neighbor : neighbors) {
      double edgeWeight = 0.0;
      if (G.getWeight(thisNode.first, neighbor, edgeWeight)) {
				double pathDist = thisNode.second + edgeWeight;  // this need to be double
				if (pathDist < distance[neighbor]) {
					pred[neighbor] = thisNode.first;  // Why ta added this
					distance[neighbor] = pathDist;
					pq.push(make_pair(neighbor, pathDist));
				}
			}
    }
  }
  return pred;
}
//
// Implement your standard application here
// TO DO: add a parameter for the graph you make.
//
void application(
    map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways,
    vector<BuildingInfo>& Buildings, graph <long long, double> &G) {
  string person1Building, person2Building;
	BuildingInfo building1;
	BuildingInfo building2;
	bool found1 = false;
	bool found2 = false;
  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);
  while (person1Building != "#") {
		cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);
    //
    // TO DO: lookup buildings, find nearest start and dest nodes, find center
    // run Dijkstra's alg from each start,
    // output distances and paths to destination:
    //
		// MS 7, searching for buildings
		found1 = searchBuilding(person1Building, Buildings, building1);
		found2 = searchBuilding(person2Building, Buildings, building2);
		if (found1 != true) {
			cout << "Person 1's building not found" << endl;
			cout << endl;
			cout << "Enter person 1's building (partial name or abbreviation), or #> ";
			getline(cin, person1Building);
			continue;
		}
		if (found2 != true) {
			cout << "Person 2's building not found" << endl;
			cout << endl;
			cout << "Enter person 1's building (partial name or abbreviation), or #> ";
			getline(cin, person1Building);
			continue;
		}
		// print if both are found
		if (found1 == true && found2 == true) {
			// Printing building info
			cout << "Person 1's point:" << endl;
			cout << " " << building1.Fullname << endl;
			cout << " (" << building1.Coords.Lat << ", "
			<< building1.Coords.Lon << ")" << endl;
			// Printing building info
			cout << "Person 2's point:" << endl;
			cout << " " << building2.Fullname << endl;
			cout << " (" << building2.Coords.Lat << ", "
			<< building2.Coords.Lon << ")" << endl;
		}
		Coordinates midpoint;  // move this out of the while loop
		BuildingInfo buildingCenter;
		midpoint = centerBetween2Points(building1.Coords.Lat,
																		building1.Coords.Lon,
																		building2.Coords.Lat,
																		building2.Coords.Lon);
		set <long long> visited;
		long long ID1;
		long long ID2;
		long long IDCenter;
		map <long long, double> distance;
		map <long long, double> distance1;
		map <long long, double> distance2;
		map <long long, long long> path1;
			map <long long, long long> path2;
		// 8 - 10 should be in a loop
		// MS 8, Locating the center building
		bool reachable1 = false;
		bool reachable2 = false;
		bool reachableDest = true;
		while (reachable1 == false && reachable2 == false) {
			buildingCenter = nearestBuilding(Buildings, midpoint, visited);
			visited.insert(buildingCenter.Coords.ID);
			// Printing building info
			if (visited.size() == 1) {
				cout << "Destination Building:" << endl;
				cout << " " << buildingCenter.Fullname << endl;
				cout << " (" << buildingCenter.Coords.Lat << ", "
				<< buildingCenter.Coords.Lon << ")" << endl;
				cout << endl;
			} else {
				cout << "New destination building:" << endl;
				cout << " " << buildingCenter.Fullname << endl;
				cout << " (" << buildingCenter.Coords.Lat << ", "
				<< buildingCenter.Coords.Lon << ")" << endl;
				IDCenter = nearestNode(Nodes, Footways, buildingCenter);
				cout << "Nearest destination node:" << endl;
				cout << " " << IDCenter << endl;
				cout << " (" << Nodes[IDCenter].Lat << ", "
				<< Nodes[IDCenter].Lon << ")" << endl;
				cout << endl;
			}
			// MS 9, get all nearest Nodes' id
			if (visited.size() == 1) {
				ID1 = nearestNode(Nodes, Footways, building1);
				cout << "Nearest P1 node:" << endl;
				cout << " " << ID1 << endl;
				cout << " (" << Nodes[ID1].Lat << ", " << Nodes[ID1].Lon << ")" << endl;
				ID2 = nearestNode(Nodes, Footways, building2);
				cout << "Nearest P2 node:" << endl;
				cout << " " << ID2 << endl;
				cout << " (" << Nodes[ID2].Lat << ", " << Nodes[ID2].Lon << ")" << endl;
				IDCenter = nearestNode(Nodes, Footways, buildingCenter);
				cout << "Nearest destination node:" << endl;
				cout << " " << IDCenter << endl;
				cout << " (" << Nodes[IDCenter].Lat << ", "
				<< Nodes[IDCenter].Lon << ")" << endl;
				cout << endl;
			}
			// MS 10, Dijkstra algorithm
			// checking if 1 to 2 is reachable
			Dijkstra(G, ID1, ID2, distance);
			if (distance[ID2] >= INF) {
				cout << "Sorry, destination unreachable." << endl;
				reachableDest = false;
				break;
			} else {
				reachableDest = true;
			}
			path1 = Dijkstra(G, ID1, IDCenter, distance1);
			path2 = Dijkstra(G, ID2, IDCenter, distance2);
			if (distance1[IDCenter] >= INF || distance2[IDCenter] >= INF) {
				cout << "At least one person was unable to reach the destination building."
				<< " Finding next closest building..." << endl;
				continue;
			} else {
					reachable1 = true;
					reachable2 = true;
				// break;
			}
		}
		if (reachableDest == false) {
			cout << endl;
			cout << "Enter person 1's building (partial name or abbreviation), or #> ";
			getline(cin, person1Building);
			continue;
		}
		// MS 11, create a path stack
		vector <long long> pred1;
		vector <long long> pred2;
		long long p1 = IDCenter;  // Push the path to the vector from the destination
		double totalDist1 = distance1[p1];
		while (p1 != 0) {
			pred1.push_back(p1);
			p1 = path1[p1];
		}
		long long p2 = IDCenter;
		double totalDist2 = distance2[p2];
		while (p2 != 0) {
			pred2.push_back(p2);
			p2 = path2[p2];
		}
		// Printing out the path
		cout << "Person 1's distance to dest: " << totalDist1 << " miles" << endl;
		cout << "Path: ";
		for (int i = pred1.size() - 1; i >= 0; i--) {
			if (i > 0) {
				cout << pred1.at(i) << "->";
			} else if (i == 0) {
				cout << pred1.at(i) << endl;
			}
		}
		cout << endl;
		cout << "Person 2's distance to dest: " << totalDist2 << " miles" << endl;
		cout << "Path: ";
		for (int i = pred2.size()-1; i >= 0; i--) {
			if (i > 0) {
				cout << pred2.at(i) << "->";
			} else if (i == 0) {
				cout << pred2.at(i) << endl;
			}
		}
    // another navigation?
    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  } 
}

int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates>  Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo>          Footways;
  // info about each building, in no particular order
  vector<BuildingInfo>         Buildings;
  XMLDocument                  xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;


  graph <long long, double> G;
  // ms5
  for (auto p : Nodes) {  // adding IDs to the graph
    G.addVertex(p.first);
  }
  // ms6
  for (auto i : Footways) {  // adding edges to the graph
    for (unsigned int e = 0; e < i.Nodes.size()-1; e++) {
      double distance = distBetween2Points(Nodes[i.Nodes[e]].Lat,
																					 Nodes[i.Nodes[e]].Lon,
																					 Nodes[i.Nodes[e + 1]].Lat,
																					 Nodes[i.Nodes[e + 1]].Lon);
      G.addEdge(i.Nodes[e], i.Nodes[e + 1], distance);
			G.addEdge(i.Nodes[e + 1], i.Nodes[e], distance);
    }
  }
  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  //
  // Menu
  //
  string userInput;
  cout << "Enter \"a\" for the standard application or "
        << "\"c\" for the creative component application> ";
  getline(cin, userInput);
  if (userInput == "a") {
    // TO DO: add argument for the graph you make.
    application(Nodes, Footways, Buildings, G);
  } else if (userInput == "c") {
    // TO DO: add arguments
    creative();
  }
  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}
