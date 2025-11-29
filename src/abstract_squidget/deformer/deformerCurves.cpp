#include "abstract_squidget/deformer/deformerCurves.hpp"

#include <maya/MFnNurbsCurve.h>
#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>

#include <maya/MFloatPoint.h>
#include <maya/MFloatVector.h>
#include <maya/MFnCamera.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnMesh.h>
#include <maya/MPlug.h>
#include <maya/MPoint.h>
#include <maya/MVector.h>

#include <cmath>
#include <functional>
#include <maya/MPoint.h>
#include <memory>
#include <unordered_map>
#include <vector>

#include "common/maya/toonRenderer.hpp"
#include "common/maya/factory.hpp"

namespace abstract_squidget {
namespace deformer {
MObjectArray getImplicitSquidgetCurves(MStringArray inputGeom,
                                       bool worldSpace, bool isRigCurves) {
  /**
   * @brief Gets the implicit squidget curves from the input geometry at the
   * current camera view.  These curves are created using the toon renderer in
   * 3D.
   */
  MStatus status;
  MObjectArray curves;
  MSelectionList oldList;
  MGlobal::getActiveSelectionList(oldList);

  if (isRigCurves) {
    // For pre-authored implicit squidgets, use curves from scene, not toon renders.
    for (MString g : inputGeom) {
      MGlobal::selectByName(g, MGlobal::kReplaceList);

      MSelectionList sel_list;
      MGlobal::getActiveSelectionList(sel_list);

      MDagPath geomPath;
      MFnNurbsCurve curveFn;
      MPointArray cvs;
      sel_list.getDagPath(0, geomPath);
      geomPath.extendToShapeDirectlyBelow(0);
      curveFn.setObject(geomPath);
      // curveFn.getCVs(cvs, MSpace::kWorld);

      double dlength = 0.2;
      double currLength = 0.0;
      while (currLength < curveFn.length()) {
        MPoint p;
        curveFn.getPointAtParam(curveFn.findParamFromLength(currLength), p,
                                MSpace::kWorld);
        cvs.append(p);
        currLength += dlength;
      }

      MObject copyCurve =
          common::maya::CreateNurbsCurve(cvs, 3); // Create fake implicit curves
      curves.append(copyCurve);
    }

  } else {
    common::maya::ToonRenderer r;
    curves = r.renderCurves(inputGeom, worldSpace);
  }

  // Remove unseen curves
  // M3dView view = M3dView::active3dView();
  // MDagPath camera;
  // view.getCamera(camera);

  // // Check if entire curve is visible or is partially visible
  // MObjectArray fullyVisibleCurves;
  // MObjectArray partiallyVisibleCurves;
  // for (MObject curve : curves) {
  //   MDagPath objPath = common::util::getNurbsCurveDagPathFromObject(curve);
  //   MFnNurbsCurve curveFn(objPath);

  //   MPointArray cvs;
  //   status = curveFn.getCVs(cvs, MSpace::kWorld);
  //   CHECK_MSTATUS(status);

  //   bool isVisible = true;
  //   for (MPoint p : cvs) {
  //     if (!isPointVisible(camera, p, curves)) {
  //       isVisible = false;
  //       break;
  //     }
  //   }

  //   if (!isVisible) {
  //     partiallyVisibleCurves.append(curve);
  //   } else {
  //     fullyVisibleCurves.append(curve);
  //   }
  // }


  // Combine fully visible curves into a single curves since there are a lot of
  // disconnected curves that share endpoints.
  // Graph graph;
  // buildGraphFromCurves(fullyVisibleCurves, graph);
  // linkHalfEdgesForCycles(graph);

  // TraversalResult traversalResult;
  // traverseGraph(graph, traversalResult);

  // MObjectArray combinedCurves;
  // createCurvesFromTraversal(traversalResult, combinedCurves);

  // Go through partially visible curves and find visible segments.
  // Then create new curves from these segments.
  // MObjectArray partialNewCurves;
  // for (MObject curve : curves) {
  //   MDagPath objPath = common::util::getNurbsCurveDagPathFromObject(curve);
  //   MFnNurbsCurve curveFn(objPath);

  //   MPointArray cvs;
  //   status = curveFn.getCVs(cvs, MSpace::kWorld);
  //   CHECK_MSTATUS(status);

  //   std::vector<VisibleSegment> visibleSegments =
  //       getVisibleSegments(curveFn, camera, curves);

  //   // Create new curves from visible segments
  //   for (VisibleSegment segment : visibleSegments) {
  //     MPointArray newCvs;
  //     for (unsigned int i = 0; i < cvs.length(); ++i) {
  //       double param =
  //           curveFn.findParamFromLength(curveFn.length() * i / cvs.length());
  //       if (param >= segment.startParam && param <= segment.endParam) {
  //         newCvs.append(cvs[i]);
  //       }
  //     }

  //     MObject newCurve = common::maya::CreateNurbsCurve(newCvs, 0);
  //     partialNewCurves.append(newCurve);
  //   }
  // }


  // Add fully visible curves and partially visible curves to new_curves
  MObjectArray new_curves;
  for (MObject curve : curves) {
    new_curves.append(curve);
  }

  // for (MObject curve : combinedCurves) {
  //   new_curves.append(curve);
  // }
  // for (MObject curve : partialNewCurves) {
  //   new_curves.append(curve);
  // }
  // for (MObject curve : curves) {
  //   common::util::deleteCurve(curve);
  // }

  // for (MObject curve : fullyVisibleCurves) {
  //   new_curves.append(curve);
  // }
  // for (MObject curve : fullyVisibleCurves) {
  //   common::util::deleteCurve(curve);
  // }
  // for (MObject curve : partiallyVisibleCurves) {
  //   common::util::deleteCurve(curve);
  // }

  MGlobal::setActiveSelectionList(oldList);
  return new_curves;
}

// Function to check if a point is visible from the camera
bool isPointVisible(const MDagPath &cameraDagPath, const MPoint &point,
                    const MObjectArray &exclusionObjects) {
  MStatus status;
  // Get camera position
  MFnCamera cameraFn(cameraDagPath);
  MPoint cameraPos = cameraFn.eyePoint(MSpace::kWorld);

  // Compute direction vector
  MVector direction = point - cameraPos;
  double distance = direction.length();
  direction.normalize();

  // Define ray origin and direction
  MPoint rayOrigin(cameraPos.x, cameraPos.y, cameraPos.z);
  MVector rayDirection(direction.x, direction.y, direction.z);

  // Prepare for intersection
  bool hit = false;
  MIntArray ids;

  // Iterate through all meshes in the scene
  MItDag itDag(MItDag::kDepthFirst, MFn::kMesh);
  for (; !itDag.isDone(); itDag.next()) {
    MDagPath meshDagPath;
    itDag.getPath(meshDagPath);

    // Exclude specified objects (e.g., the curves themselves)
    bool exclude = false;
    for (unsigned int i = 0; i < exclusionObjects.length(); ++i) {
      if (meshDagPath.node() == exclusionObjects[i]) {
        exclude = true;
        break;
      }
    }
    if (exclude)
      continue;

    MFloatPoint hitPoint;
    float hitParam;
    MFnMesh meshFn(meshDagPath);
    hit = meshFn.closestIntersection(
        rayOrigin, rayDirection, nullptr, nullptr, false, MSpace::kWorld,
        distance - 1.0E-1f, false, nullptr, hitPoint, &hitParam, nullptr,
        nullptr, nullptr, nullptr, 1.0E-6f, &status);

    if (hit) {
      return false; // Point is occluded
    }
  }

  return true; // Point is visible
}

std::vector<VisibleSegment>
getVisibleSegments(const MFnNurbsCurve &curveFn, const MDagPath &cameraDagPath,
                   const MObjectArray &exclusionObjects, int numSamples) {
  MStatus status;
  std::vector<VisibleSegment> visibleSegments;
  std::vector<bool> visibilityFlags(numSamples + 1, false);

  // Sample points along the curve and mark visibility
  for (int i = 0; i <= numSamples; ++i) {
    double param =
        curveFn.findParamFromLength(curveFn.length() * i / numSamples);
    MPoint point;
    status = curveFn.getPointAtParam(param, point, MSpace::kWorld);
    visibilityFlags[i] = isPointVisible(cameraDagPath, point, exclusionObjects);
  }

  // Identify continuous visible segments
  int startIdx = -1;
  for (int i = 0; i <= numSamples; ++i) {
    if (visibilityFlags[i] && startIdx == -1) {
      // Start a new segment if starting index is not set
      startIdx = i;
    }

    if (startIdx != -1 && (!visibilityFlags[i] || i == numSamples)) {
      // End the segment if visibility changes or at the end of the curve
      int endIdx = i - 1;

      // Convert indices to parameter ranges
      double startParam =
          curveFn.findParamFromLength(curveFn.length() * startIdx / numSamples);
      double endParam =
          curveFn.findParamFromLength(curveFn.length() * endIdx / numSamples);
      visibleSegments.push_back(VisibleSegment{startParam, endParam});
      startIdx = -1;
    }
  }

  return visibleSegments;
}

// =======================================================================
// Graph =================================================================
// =======================================================================

// Helper function to compare two points within a tolerance
bool pointsAreEqual(const MPoint &p1, const MPoint &p2, double tolerance) {
  return p1.distanceTo(p2) < tolerance;
}

// Function to find or create a vertex in the graph
std::shared_ptr<Vertex> findOrCreateVertex(Graph &graph, const MPoint &point,
                                           double tolerance) {
  for (auto &vertex : graph.vertices) {
    if (pointsAreEqual(vertex->position, point, tolerance)) {
      return vertex;
    }
  }

  // If not found, create a new vertex
  auto newVertex = std::make_shared<Vertex>(point);
  graph.vertices.push_back(newVertex);
  return newVertex;
}

// Function to build the graph from MObjectArray of NURBS Curves
MStatus buildGraphFromCurves(const MObjectArray &curves, Graph &graph,
                             double tolerance) {
  MStatus status;

  for (unsigned int i = 0; i < curves.length(); ++i) {
    MObject curveObj = curves[i];
    MDagPath curvePath;
    curvePath = common::util::getNurbsCurveDagPathFromObject(curveObj);
    MFnNurbsCurve fnCurve(curvePath, &status);
    CHECK_MSTATUS(status);

    // Get CVs in world space
    MPointArray cvs;
    status = fnCurve.getCVs(cvs, MSpace::kWorld);
    if (cvs.length() < 2) {
      MGlobal::displayWarning("Curve has less than 2 CVs. Skipping.");
      continue;
    }

    // Iterate through CVs and create edges
    unsigned int j = 0;
    MPoint startPt = cvs[j];
    MPoint endPt = cvs[j + 1];
    auto origin = findOrCreateVertex(graph, startPt, tolerance);
    auto destination = findOrCreateVertex(graph, endPt, tolerance);

    while (j < cvs.length() - 1) {
      // Find or create vertices

      // Create half-edges
      auto he1 = std::make_shared<HalfEdge>(origin, destination);
      auto he2 = std::make_shared<HalfEdge>(destination, origin); // Twin

      // Link twins
      he1->twin = he2;
      he2->twin = he1;

      // Add half-edges to the graph
      graph.halfEdges.push_back(he1);
      graph.halfEdges.push_back(he2);

      // Add outgoing half-edge to origin vertex
      origin->outgoingHalfEdges.push_back(he1);
      destination->outgoingHalfEdges.push_back(he2);

      ++j;
      MPoint startPt = cvs[j];
      MPoint endPt = cvs[j + 1];
      origin = destination;
      destination = findOrCreateVertex(graph, endPt, tolerance);
    }
  }

  MGlobal::displayInfo("Graph successfully built from curves.");

  return MS::kSuccess;
}

// Function to link next and prev pointers for half-edges forming cycles
void linkHalfEdgesForCycles(Graph &graph, double tolerance) {
  // Iterate through each vertex and its outgoing half-edges
  for (auto &vertex : graph.vertices) {
    // If a vertex has more than one outgoing half-edge, it might be a junction
    // For simplicity, assume planar graphs where cycles can be easily
    // identified
    if (vertex->outgoingHalfEdges.size() < 1)
      continue;

    // Sort outgoing half-edges based on angle to ensure consistent ordering
    // This is essential for correctly linking next and prev pointers
    // Calculate angles relative to a reference direction (e.g., X-axis)

    struct HalfEdgeWithAngle {
      std::shared_ptr<HalfEdge> halfEdge;
      double angle;
    };

    std::vector<HalfEdgeWithAngle> sortedHalfEdges;

    MVector refDir(1.0, 0.0, 0.0); // Reference direction (X-axis)

    for (auto &he : vertex->outgoingHalfEdges) {
      MVector dir = he->destination->position - he->origin->position;
      dir.normalize();
      double angle = refDir.angle(dir);
      // Determine the sign of the angle using cross product
      double cross = refDir.x * dir.y - refDir.y * dir.x;
      if (cross < 0)
        angle = -angle;
      sortedHalfEdges.push_back({he, angle});
    }

    // Sort the half-edges by angle in ascending order
    std::sort(
        sortedHalfEdges.begin(), sortedHalfEdges.end(),
        [](const HalfEdgeWithAngle &a, const HalfEdgeWithAngle &b) -> bool {
          return a.angle < b.angle;
        });

    // Link the half-edges in a circular manner
    size_t n = sortedHalfEdges.size();
    for (size_t k = 0; k < n; ++k) {
      auto currentHE = sortedHalfEdges[k].halfEdge;
      auto nextHE = sortedHalfEdges[(k + 1) % n].halfEdge;

      // The next half-edge in the cycle is the twin of the previous half-edge
      currentHE->next = nextHE->twin;
      nextHE->twin->prev = currentHE;
    }
  }

  MGlobal::displayInfo("Half-edges linked for cycles.");
}

// Function to detect cycles and open paths using next pointers
void traverseGraph(const Graph &graph, TraversalResult &result) {
  for (const auto &he : graph.halfEdges) {
    if (he->visited)
      continue;

    // Start traversal from this half-edge
    std::shared_ptr<HalfEdge> startHE = he;
    std::shared_ptr<HalfEdge> currentHE = startHE;
    std::vector<std::shared_ptr<Vertex>> path;
    bool isCycle = true;

    while (currentHE != nullptr && !currentHE->visited) {
      currentHE->visited = true;
      path.push_back(currentHE->origin);
      // Move to next half-edge
      currentHE = currentHE->next;
    }

    if (currentHE == startHE) {
      // Completed a cycle
      result.cycles.push_back(path);
    } else {
      // Reached an open end
      result.openPaths.push_back(path);
    }
  }

  // Delete cycles that share the exact same vertices.
  // Check by putting all vertices in a set and checking if other vertices
  // create the same set.
  std::vector<std::set<std::string>> vertexSets;
  for (unsigned int i = 0; i < result.cycles.size(); ++i) {
    std::set<std::string> vertexSet;
    for (unsigned int j = 0; j < result.cycles[i].size(); ++j) {
      vertexSet.insert(result.cycles[i][j]->name);
    }
    vertexSets.push_back(vertexSet);
  }

  for (unsigned int i = 0; i < result.cycles.size(); ++i) {
    for (unsigned int j = i + 1; j < result.cycles.size(); ++j) {
      if (vertexSets[i] == vertexSets[j]) {
        result.cycles.erase(result.cycles.begin() + j);
        vertexSets.erase(vertexSets.begin() + j);
        --j;
      }
    }
  }
}

// Function to create a NURBS curve from a list of vertices
MObject createNurbsCurveFromVertices(
    const std::vector<std::shared_ptr<Vertex>> &vertices, bool isClosed,
    MStatus &status) {
  MPointArray cvs;
  for (const auto &vertex : vertices) {
    cvs.append(vertex->position);
  }

  // Optionally, ensure the curve is closed by adding the first point at the end
  if (isClosed && (cvs[0] != cvs[cvs.length() - 1])) {
    cvs.append(cvs[0]);
  }

  // Create the curve
  MFnNurbsCurve curveFn;
  MObject curve = common::maya::CreateNurbsCurve(cvs, 0);

  if (status != MS::kSuccess) {
    MGlobal::displayError("Failed to create NURBS curve.");
    return MObject::kNullObj;
  }

  return curve;
}

// Function to create curves from traversal results
MStatus createCurvesFromTraversal(const TraversalResult &traversal,
                                  MObjectArray &newCurves) {
  MStatus status;
  MObjectArray tempCurves;

  // Create cycles as closed curves
  for (const auto &cycle : traversal.cycles) {
    if (cycle.size() < 2)
      continue; // Need at least two vertices
    MObject curveObj = createNurbsCurveFromVertices(cycle, true, status);
    if (status == MS::kSuccess) {
      tempCurves.append(curveObj);
    }
  }

  // Create open paths as open curves
  for (const auto &path : traversal.openPaths) {
    if (path.size() < 2)
      continue; // Need at least two vertices
    MObject curveObj = createNurbsCurveFromVertices(path, false, status);
    if (status == MS::kSuccess) {
      tempCurves.append(curveObj);
    }
  }

  newCurves = tempCurves;
  return MS::kSuccess;
}

MPointArray findShadowCurve(MPointArray _strokePoints2D, MPointArray cvs,
                            RenderMatrices renderMats, double noise) {
  /**
   * @brief pojects strokePoints3D onto cvs using curveObj MFnCurve
   parameters
   * in screen space.  Returns the projected points in object space.
   *
   * DOESN"T RETURN PROJECTION.  RETURNS SAMPLED OUTLINE NOW.
   * Why did I do this again?
   *
   * 1) project both strokePoints2D and cvs to screen space
   * 2) Move the stroke points closer to cvs2D using average distance
   * 3) Project stroke points onto the curve
   * 4) Order the parameter points
   * 5) Find the mapping of the parameter points to object space points
   */

  // 0) Get cv points
  MObject curveObj = common::maya::CreateNurbsCurve(cvs, 0);
  MDagPath curvePath = common::util::getNurbsCurveDagPathFromObject(curveObj);
  MFnNurbsCurve curveFn(curvePath);
  MPointArray sampled_cvs;
  curveFn.getCVs(sampled_cvs);

  // 1) Project cvs to screen space
  MPointArray worldCVs;
  MPointArray screenCVs;
  for (MPoint p : sampled_cvs) {
    MPoint worldp = renderMats.worldMatrix * renderMats.localMatrix * p;
    worldCVs.append(worldp);

    MPoint screenp = renderMats.screenMatrix * renderMats.projMatrix *
                     renderMats.cameraMatrix * worldp;
    screenp.cartesianize();
    screenCVs.append(screenp);
  }

  // 1.5) Get the initial closest point from strokePoints2D to cvs
  curveFn.setCVs(screenCVs);
  MPointArray closestPoints;
  for (MPoint p : _strokePoints2D) {
    double param;
    MPoint closestPoint = curveFn.closestPoint(p, &param);
    closestPoints.append(closestPoint);
  }

  // 2) Move the stroke points closer to cvs2D using means of closest points
  // Projecting directly works bad for far away points.
  MPointArray strokePoints2D(_strokePoints2D);
  MPoint strokePoints2DMean = MPoint::origin;
  for (MPoint p : strokePoints2D) {
    strokePoints2DMean += p;
  }
  strokePoints2DMean = strokePoints2DMean / strokePoints2D.length();

  MPointArray closestPointsMean(closestPoints);
  MPoint closestPoints2DMean = MPoint::origin;
  for (MPoint p : closestPointsMean) {
    closestPoints2DMean += p;
  }
  closestPoints2DMean = closestPoints2DMean / closestPointsMean.length();

  MVector diff = closestPoints2DMean - strokePoints2DMean;
  for (size_t i = 0; i < strokePoints2D.length(); ++i) {
    strokePoints2D[i] += diff;
  }

  // 3) Project stroke points onto the curve
  // We use the parameter however to find the corresponding worldCVs point.
  MDoubleArray shadowParam;
  for (MPoint p : strokePoints2D) {
    double param;
    MPoint shadowPoint = curveFn.closestPoint(p, &param);
    shadowParam.append(param);
  }

  // 3.5) Remove parameters that are too close to each other
  for (size_t i = 0; i < shadowParam.length() - 1; ++i) {
    if (abs(shadowParam[i] - shadowParam[i + 1]) < 0.01) {
      shadowParam.remove(i);
      i--;
    }
  }

  // 4) Order the parameter points

  // 5) Find the mapping of the parameter points to object space points
  // Find shadow curve in world space and then convert to object space.
  curveFn.setCVs(worldCVs); // Reset original cvs
  MPointArray screenShadowCurve;
  MPointArray worldShadowCurve;
  MPointArray objSpaceShadowCurve;
  for (double param : shadowParam) {
    MPoint shadowCVPoint;
    curveFn.getPointAtParam(param, shadowCVPoint);

    // Project to screen space
    MPoint shadowCVScreen = renderMats.screenMatrix * renderMats.projMatrix *
                            renderMats.cameraMatrix * shadowCVPoint;
    shadowCVScreen.cartesianize();

    // Convert to object space
    MPoint objp = (renderMats.worldMatrix * renderMats.localMatrix).inverse() *
                  shadowCVPoint;

    screenShadowCurve.append(shadowCVScreen);
    worldShadowCurve.append(shadowCVPoint);
    objSpaceShadowCurve.append(objp);
  }
  
  common::util::deleteCurve(curveObj);
  return objSpaceShadowCurve;
  // return sampled_cvs;
}
} // namespace deformer
} // namespace abstract_squidget
