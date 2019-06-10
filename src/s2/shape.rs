// Copyright 2017 Google Inc. All rights reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


use r2::point::Point;
use r3::vector;

// Edge represents a geodesic edge consisting of two vertices. Zero-length edges are
// allowed, and can be used to represent points.
pub struct Edge {
	V0: Point, 
	V1: Point
}

impl Edge {
	// Cmp compares the two edges using the underlying Points Cmp method and returns
	//
	//   -1 if e <  other
	//    0 if e == other
	//   +1 if e >  other
	//
	// The two edges are compared by first vertex, and then by the second vertex.
	pub fn Cmp(&self, other: Edge) -> i64 {
		let v0cmp = Cmp(e.V0, other.V0);
		if (v0cmp != 0) {
			return v0cmp;
		} else {
			Cmp(e.V1, other.V1);
		}
	}
}

// sortEdges sorts the slice of Edges in place.
pub fn sortEdges(e: Vec<Edge>) -> Vec<Edge> {
	e.sort()
}

// edges implements the Sort interface for slices of Edge.
type Edges = Vec<Edge>

impl Edges {

	pub fn Swap(&self, i: i64, j: i64) { self[i], self[j] = self[j], self[i] }
	pub fn Less(&self, i: i64, j: i64) -> bool { self[i].Cmp(self[j]) == -1 }

}

// ShapeEdgeID is a unique identifier for an Edge within an ShapeIndex,
// consisting of a (shapeID, edgeID) pair.
pub struct ShapeEdgeID {
	ShapeID: i32,
	EdgeID:  i32
}

impl ShapeEdgeID {
	// Cmp compares the two ShapeEdgeIDs and returns
	//
	//   -1 if s <  other
	//    0 if s == other
	//   +1 if s >  other
	//
	// The two are compared first by shape id and then by edge id.
	pub fn Cmp(&self, other: ShapeEdgeID) -> i64 {
		if (s.ShapeID < other.ShapeID) {
			return -1
		} else if (case s.ShapeID > other.ShapeID) {
			return 1
		}

		if (case s.EdgeID < other.EdgeID) {
			return -1
		} else if (case s.EdgeID > other.EdgeID)
			return 1
		}
		return 0
	}
}

// ShapeEdge represents a ShapeEdgeID with the two endpoints of that Edge.
pub struct ShapeEdge {
	ID: ShapeEdgeID,
	Edge: Edge
}

// Chain represents a range of edge IDs corresponding to a chain of connected
// edges, specified as a (start, length) pair. The chain is defined to consist of
// edge IDs {start, start + 1, ..., start + length - 1}.
pub struct Chain {
	Start: i64, 
	Length: i64
}

// ChainPosition represents the position of an edge within a given edge chain,
// specified as a (chainID, offset) pair. Chains are numbered sequentially
// starting from zero, and offsets are measured from the start of each chain.
pub struct ChainPosition {
	ChainID: i64, 
	Offset: i64
}

// A ReferencePoint consists of a point and a boolean indicating whether the point
// is contained by a particular shape.
pub struct ReferencePoint {
	Point: Point,
	Contained: bool
}

// OriginReferencePoint returns a ReferencePoint with the given value for
// contained and the origin point. It should be used when all points or no
// points are contained.
pub fn OriginReferencePoint(contained bool) -> ReferencePoint {
	ReferencePoint{Point: OriginPoint(), Contained: contained}
}

// Shape represents polygonal geometry in a flexible way. It is organized as a
// collection of edges that optionally defines an interior. All geometry
// represented by a given Shape must have the same dimension, which means that
// an Shape can represent either a set of points, a set of polylines, or a set
// of polygons.
//
// Shape is defined as an interface in order to give clients control over the
// underlying data representation. Sometimes an Shape does not have any data of
// its own, but instead wraps some other type.
//
// Shape operations are typically defined on a ShapeIndex rather than
// individual shapes. An ShapeIndex is simply a collection of Shapes,
// possibly of different dimensions (e.g. 10 points and 3 polygons), organized
// into a data structure for efficient edge access.
//
// The edges of a Shape are indexed by a contiguous range of edge IDs
// starting at 0. The edges are further subdivided into chains, where each
// chain consists of a sequence of edges connected end-to-end (a polyline).
// For example, a Shape representing two polylines AB and CDE would have
// three edges (AB, CD, DE) grouped into two chains: (AB) and (CD, DE).
// Similarly, an Shape representing 5 points would have 5 chains consisting
// of one edge each.
//
// Shape has methods that allow edges to be accessed either using the global
// numbering (edge ID) or within a particular chain. The global numbering is
// sufficient for most purposes, but the chain representation is useful for
// certain algorithms such as intersection (see BooleanOperation).
pub trait Shape {
	// NumEdges returns the number of edges in this shape.
	fn num_edges() -> i64;

	// Edge returns the edge for the given edge index.
	fn Edge(i int) -> Edge;

	// ReferencePoint returns an arbitrary reference point for the shape. (The
	// containment boolean value must be false for shapes that do not have an interior.)
	//
	// This reference point may then be used to compute the containment of other
	// points by counting edge crossings.
	fn ReferencePoint() -> ReferencePoint

	// NumChains reports the number of contiguous edge chains in the shape.
	// For example, a shape whose edges are [AB, BC, CD, AE, EF] would consist
	// of two chains (AB,BC,CD and AE,EF). Every chain is assigned a chain Id
	// numbered sequentially starting from zero.
	//
	// Note that it is always acceptable to implement this method by returning
	// NumEdges, i.e. every chain consists of a single edge, but this may
	// reduce the efficiency of some algorithms.
	fn NumChains() -> i64;

	// Chain returns the range of edge IDs corresponding to the given edge chain.
	// Edge chains must form contiguous, non-overlapping ranges that cover
	// the entire range of edge IDs. This is spelled out more formally below:
	//
	//  0 <= i < NumChains()
	//  Chain(i).length > 0, for all i
	//  Chain(0).start == 0
	//  Chain(i).start + Chain(i).length == Chain(i+1).start, for i < NumChains()-1
	//  Chain(i).start + Chain(i).length == NumEdges(), for i == NumChains()-1
	fn Chain(chainID int) -> Chain;

	// ChainEdgeReturns the edge at offset "offset" within edge chain "chainID".
	// Equivalent to "shape.Edge(shape.Chain(chainID).start + offset)"
	// but more efficient.
	fn ChainEdge(chainID, offset int) -> Edge;

	// ChainPosition finds the chain containing the given edge, and returns the
	// position of that edge as a ChainPosition(chainID, offset) pair.
	//
	//  shape.Chain(pos.chainID).start + pos.offset == edgeID
	//  shape.Chain(pos.chainID+1).start > edgeID
	//
	// where pos == shape.ChainPosition(edgeID).
	fn ChainPosition(edgeID int) ChainPosition;

	// Dimension returns the dimension of the geometry represented by this shape,
	// either 0, 1 or 2 for point, polyline and polygon geometry respectively.
	//
	//  0 - Point geometry. Each point is represented as a degenerate edge.
	//
	//  1 - Polyline geometry. Polyline edges may be degenerate. A shape may
	//      represent any number of polylines. Polylines edges may intersect.
	//
	//  2 - Polygon geometry. Edges should be oriented such that the polygon
	//      interior is always on the left. In theory the edges may be returned
	//      in any order, but typically the edges are organized as a collection
	//      of edge chains where each chain represents one polygon loop.
	//      Polygons may have degeneracies (e.g., degenerate edges or sibling
	//      pairs consisting of an edge and its corresponding reversed edge).
	//      A polygon loop may also be full (containing all points on the
	//      sphere); by convention this is represented as a chain with no edges.
	//      (See laxPolygon for details.)
	//
	// This method allows degenerate geometry of different dimensions
	// to be distinguished, e.g. it allows a point to be distinguished from a
	// polyline or polygon that has been simplified to a single point.
	fn Dimension() -> i64;

	// IsEmpty reports whether the Shape contains no points. (Note that the full
	// polygon is represented as a chain with zero edges.)
	fn IsEmpty() -> bool;

	// IsFull reports whether the Shape contains all points on the sphere.
	fn IsFull() -> bool;

}

// defaultShapeIsEmpty reports whether this shape contains no points.
pub fn defaultShapeIsEmpty(s: Shape) -> bool {
	return s.NumEdges() == 0 && (s.Dimension() != 2 || s.NumChains() == 0)
}

// defaultShapeIsFull reports whether this shape contains all points on the sphere.
pub fn defaultShapeIsFull(s: Shape) -> bool {
	return s.NumEdges() == 0 && s.Dimension() == 2 && s.NumChains() > 0
}