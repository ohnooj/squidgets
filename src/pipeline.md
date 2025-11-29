# How the Tool Works (Sorta)

We use the Maya C++ API to implement squidgets. This is the abstract squidget pipeline.

0. Plugin Loading and Setup:
   1. The plugin is loaded through the `pluginMain.cpp` file which loads an MPxContext `abstractContext` and MPxToolCommand `abstractSquidgetTool` together through an MPxContextCommand `abstractContextCmd`.

1. Stroke Input: MPxContext
   1. All mouse and keyboard IO (such as mouse press, drag, release, enter key) and overlaying scene UI display (drawing a stroke) is handled in the `abstractContext` class. A 2D stroke in screen space is created and rendered for the user here.
   2. The `abstractSquidgetTool` class is the Maya command that is responsible for doing the squidget operations. All maya actions are counted as commands with arguments, thus we have to follow this pattern as well. The arguments are set up and called in `AbstractContext::callToolCommand`.

2. Maya Tool Command: MPxToolCommand
   1. This class is responsible for the squidget's behavior. There is a selection/query phase followed by a deformation/manipulation phase.
   2. This command also allows for undoing and redoing squidget operations in the Maya software.

3. Squidget Selection: `abstractDeformer::queryClosestObject()`
   1. The squidget object is queried in `abstractDeformer::queryClosestObject()`.
   2. Abstract Squidgets curves are created using the Maya Toon Rendrer at the current camera perspective.
   3. Each squidget curve is projected into the 2D viewport where the 2D stroke is.  Curves are matched in screen space.
      1. We choose to simply project the 2D stroke onto the 2D squidget curve to get a projected shadow curve. You can replace this with a more sophisticated curve matching algorithm.  Iteratively doing this a couple times does get a pretty okay shadow projection curve.
   4. We perform a 2D curve registration algorithm to match the 2D stroke to the squidget curve. The output of this is a 2D transformation matrix to move the 2D stroke to the squidget curve.
      1. The transformation matrix is inverted such that the squidget curve is moved to the 2D stroke.
   5. Rather than performing a transformation on the 2D squidget curve, we need to perform the transformation on the object in the object's local space. We find the transformation matrix B that manipulates object s.t. the 2D projection of the object matches the 2D stroke.
   6. We calculate a score on this matrix using 2 metrics:
      1. Transformation displacement.
      2. Point 2 Projection distance from the stroke to the transformed 2D squidget curve.
   7. The scores are placed into a priority queue and the object-transformation pair with the highest score is selected.

4. Squidget Deformation: `abstractDeformer::transformObject()` and `abstractDeformer::deformObject()`
   1. Our tool currently supports translation, translation + rotation, and 2D vertex deformation.
   2. For object transformation, we apply the transformation matrix B to the object.
   3. For vertex deformation, we project the 2D stroke onto the 2D squidget curve to capture an interval of vertices [vi, vj] to deform.  We then deform the vertices in the object's local space.
      1. There can be fall-off applied to the deformation at the end points.
      2. Note this is a very basic vertex deformation.  Other deformers (wire deformer, skin deformer) can be used to deform the object.

*Notes from the author:*
The original paper was written with bookmark curves in mind back in 2023 for CHI 2024; abstraction curves were added later for CHI 2025 and UIST 2025.  Thus, the bookmark curves code is much older than the abstraction curves code.  I hardly touched the bookmark curves code for the 2025 deadlines, so that code is less robust.  Around this time, I refactored abstraction curves as a separate tool, that is why there are some duplicate code between bookmark curves and abstraction curves.

The bookmark curves work sorta similarly with much of the code base written to organize the canvases to 1) draw squidget strokes and 2) couple strokes to object scene attributes automatically.  This part is very messy.

This code is written for Maya, so a majority of the code is fluff for the Maya API integration.  I don't think I will spend ime cleaning this up much, but feel free to reach out if you have questions.

Overall, the project is mainly focused on translation and deformation.  Theoretically, other attributes can be manipulated as we discuss in the paper, but the practical problem is the generalization and scope of the problem.  I don't think this method will replace simple tasks like translation, but deformation I think is a promising direction, or even animation.

Also, this code is written without AI assistance.  :)
