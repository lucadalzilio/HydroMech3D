//
// Created by Federico Ciardo on 28.07.21. All rights reserved.
//

// Inclusion from InsideLibrary
#include <armadillo>
#include <cassert>
// Import from the project
#include "ElementData.h"

// Import from standard library
#include <algorithm>

#ifndef INC_3DEQSIM_SRC_MESH_H
#define INC_3DEQSIM_SRC_MESH_H

namespace EQSim {

class Mesh {
 private:
  arma::mat coordinates_;
  arma::imat nodes_connectivity_;
  arma::imat edges_connectivity_;
  arma::uword interpolation_order_ = 0;
  arma::mat centroids_;
  arma::ivec indexes_inner_centroids_;
  arma::ivec boundary_elements_;

 public:
  // Constructors
  Mesh() = default;

  Mesh(const arma::mat &Coordinates,
       const arma::imat &Connectivity,
       const arma::uword &interpolationOrder) {
    assert(Coordinates.n_rows > 1 && Coordinates.n_cols == 3);
    assert(interpolationOrder == 0 || interpolationOrder == 4);

    coordinates_ = Coordinates;
    nodes_connectivity_ = Connectivity;
    interpolation_order_ = interpolationOrder;

    arma::uword Nelts = nodes_connectivity_.n_rows;
    arma::mat Centroids(Nelts, 3, arma::fill::zeros);

    for (arma::uword I = 0; I < Nelts; ++I) {
      Centroids(I, 0) = (Coordinates(Connectivity(I, 0), 0) +
                         Coordinates(Connectivity(I, 1), 0) +
                         Coordinates(Connectivity(I, 2), 0) +
                         Coordinates(Connectivity(I, 3), 0)) /
                        4.;
      Centroids(I, 1) = (Coordinates(Connectivity(I, 0), 1) +
                         Coordinates(Connectivity(I, 1), 1) +
                         Coordinates(Connectivity(I, 2), 1) +
                         Coordinates(Connectivity(I, 3), 1)) /
                        4.;
      Centroids(I, 2) = (Coordinates(Connectivity(I, 0), 2) +
                         Coordinates(Connectivity(I, 1), 2) +
                         Coordinates(Connectivity(I, 2), 2) +
                         Coordinates(Connectivity(I, 3), 2)) /
                        4.;
    }
    centroids_ = Centroids;

    arma::imat EdgeConnectivity(Nelts, 8, arma::fill::zeros);
      for (arma::uword I = 0; I < Nelts; ++I) {
        EdgeConnectivity(I, 0) = Connectivity(I, 0);
        EdgeConnectivity(I, 1) = Connectivity(I, 1);
        EdgeConnectivity(I, 2) = Connectivity(I, 1);
        EdgeConnectivity(I, 3) = Connectivity(I, 2);
        EdgeConnectivity(I, 4) = Connectivity(I, 2);
        EdgeConnectivity(I, 5) = Connectivity(I, 3);
        EdgeConnectivity(I, 6) = Connectivity(I, 3);
        EdgeConnectivity(I, 7) = Connectivity(I, 0);
      }
    edges_connectivity_ = EdgeConnectivity;
  }

  // Getter methods
  arma::mat getCoordinates() const { return coordinates_; };
  double getCoordinates(arma::uword k, arma::uword i) const {
    return coordinates_(k, i);
  };

  arma::imat getNodesConnettivity() const {
    return nodes_connectivity_;
  };
  arma::uword getNodesConnettivity(arma::uword k, arma::uword i) const {
    return nodes_connectivity_(k, i);
  };

  arma::imat getEdgesConnettivity() const {
    return edges_connectivity_;
  };
  arma::uword getEdgesConnettivity(arma::uword k, arma::uword i) const {
    return edges_connectivity_(k, i);
  };

  arma::mat getCentroids() const { return centroids_; };
  double getCentroids(arma::uword k, arma::uword i) const {
    return centroids_(k, i);
  };

  arma::mat getInnerCentroids(arma::imat &neigh_elts);

  arma::ivec getIndexesInnerCentroids() const {
    return indexes_inner_centroids_;
  };

  arma::uword getIndexesInnerCentroids(arma::uword centroid_i) const {
    return indexes_inner_centroids_[centroid_i];
  };

  arma::uword getInterpolationOrder() const { return interpolation_order_; };

  arma::uword getNumberOfElts() const { return nodes_connectivity_.n_rows; };

  arma::uword getNumberOfNodes() const { return coordinates_.n_rows; };

  arma::uword getNumberOfDofs() const { return 3 * nodes_connectivity_.n_rows; };

  EQSim::ElementData getElementData(arma::uword ne) const;

  arma::imat getNeighbourElements_NonUniformMesh() const;
  arma::imat getNeighbourElements_UniformMesh(
      EQSim::Mesh &Mesh);

  arma::ivec getBoundaryElements() const {
    return boundary_elements_;
  };

  // Setter methods
  void setCentroids(arma::mat &centroids) { centroids_ = centroids; };
};

}  // namespace EQSim
#endif  // INC_3DEQSIM_SRC_MESH_H
