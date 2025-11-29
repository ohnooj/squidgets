#ifndef DEFORMER_CURVES_HPP
#define DEFORMER_CURVES_HPP

#include <maya/MDagPath.h>
#include <maya/MObjectArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MStringArray.h>

#include <memory>
#include <set>
#include <stack>
#include <vector>

#include "abstract_squidget/types.hpp"

namespace abstract_squidget {
namespace deformer {

struct VisibleSegment {
  double startParam;
  double endParam;
};

MObjectArray getImplicitSquidgetCurves(MStringArray inputGeom, bool worldSpace,
                                       bool isRigCurves);
bool isPointVisible(const MDagPath &cameraDagPath, const MPoint &point,
                    const MObjectArray &exclusionObjects);

/**
 * @brief Go along the curve and find continuous visible segments.
 *
 * @param curveFn MFnNurbsCurve object
 * @param cameraDagPath MDagPath to the camera
 * @param exclusionObjects MObjectArray of objects to exclude from visibility
 * @param numSamples Number of samples to take along the curve (actually
 * numSamples + 1)
 *
 * @return std::vector<VisibleSegment> List of visible segments
 */
std::vector<VisibleSegment>
getVisibleSegments(const MFnNurbsCurve &curveFn, const MDagPath &cameraDagPath,
                   const MObjectArray &exclusionObjects, int numSamples = 100);

MPointArray findShadowCurve(MPointArray strokePoints2D, MPointArray cvs,
                            RenderMatrices renderMats, double noise);

// Forward declarations
struct HalfEdge;
struct Vertex;

// Vertex Structure
struct Vertex {
  MPoint position;
  std::vector<std::shared_ptr<HalfEdge>> outgoingHalfEdges;
  std::string name;

  Vertex(const MPoint &pos) : position(pos) {
    static int id = 0;
    name = "v" + std::to_string(id++);
  }
};

// HalfEdge Structure
struct HalfEdge {
  std::shared_ptr<Vertex> origin;      // Starting vertex
  std::shared_ptr<HalfEdge> twin;      // Opposite half-edge
  std::shared_ptr<HalfEdge> next;      // Next half-edge around the face
  std::shared_ptr<HalfEdge> prev;      // Previous half-edge around the face
  std::shared_ptr<Vertex> destination; // Ending vertex
  bool visited;                        // Flag to indicate if traversed

  HalfEdge(std::shared_ptr<Vertex> orig, std::shared_ptr<Vertex> dest)
      : origin(orig), destination(dest) {}
};

// Graph Structure
struct Graph {
  // All vertices
  std::vector<std::shared_ptr<Vertex>> vertices;

  // All half-edges
  std::vector<std::shared_ptr<HalfEdge>> halfEdges;

  // Mapping from position to vertex for quick lookup (using a spatial hash or
  // similar) For simplicity, we'll use a vector and linear search with
  // tolerance
};

// Structure to hold traversal results
struct TraversalResult {
  std::vector<std::vector<std::shared_ptr<Vertex>>> cycles; // Detected cycles
  std::vector<std::vector<std::shared_ptr<Vertex>>>
      openPaths; // Detected open paths
};

bool pointsAreEqual(const MPoint &p1, const MPoint &p2,
                    double tolerance = 1e-5);
std::shared_ptr<Vertex> findOrCreateVertex(Graph &graph, const MPoint &point,
                                           double tolerance);
MStatus buildGraphFromCurves(const MObjectArray &curves, Graph &graph,
                             double tolerance = 1e-5);
void linkHalfEdgesForCycles(Graph &graph, double tolerance = 1e-5);

void traverseGraph(const Graph &graph, TraversalResult &result);

MObject createNurbsCurveFromVertices(
    const std::vector<std::shared_ptr<Vertex>> &vertices, bool isClosed,
    MStatus &status);
MStatus createCurvesFromTraversal(const TraversalResult &traversal,
                                  MObjectArray &newCurves);

} // namespace deformer
} // namespace abstract_squidget
#endif // DEFORMER_CURVES_HPP