/*
    headless.h -- headless "GUI parity" processing entry points

    This file adds a batch/headless API layer on top of Instant Meshes that can:
    - load a mesh
    - apply constraint strokes (orientation brush + edge brush)
    - solve orientation/position fields using the same Optimizer kernels as the GUI
    - export/import intermediate state (Serializer PLY)
*/

#pragma once

#include "common.h"
#include <set>

struct HeadlessOptions {
    std::string mode; // orientation | position | extract | applyConstraints
    std::string inputPath; // mesh input (orientation)
    std::string statePath; // serialized state input (position/extract/applyConstraints)
    std::string constraintsPath; // constraints text file
    std::string outputStatePath; // serialized state output
    std::string outputMeshPath; // extracted mesh output (extract)

    int rosy = 4;
    int posy = 4;
    bool extrinsic = true;
    bool deterministic = false;
    bool alignToBoundaries = false;
    bool dominant = false;
    int smoothIter = 2;
    Float creaseAngle = -1;
    Float scale = -1;
    int faceCount = -1;
    int vertexCount = -1;
};

extern int headless_process(const HeadlessOptions &options);
