/*
    headless.cpp -- headless "GUI parity" processing
*/

#include "headless.h"

#include "meshio.h"
#include "meshstats.h"
#include "subdivide.h"
#include "normal.h"
#include "bvh.h"
#include "field.h"
#include "extract.h"
#include "hierarchy.h"
#include "serializer.h"
#include "dedge.h"
#include "adjacency.h"

#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>

/* Defined in field.cpp (used by GUI), but not declared in headers */
bool move_orientation_singularity(MultiResolutionHierarchy &mRes, uint32_t f_src, uint32_t f_target);
bool move_position_singularity(MultiResolutionHierarchy &mRes, uint32_t f_src, uint32_t f_target);

struct ConstraintSample {
    uint32_t f = 0;
    Vector3f p = Vector3f::Zero();
};

struct StrokeConstraint {
    uint32_t type = 0; // 0=orientation brush, 1=edge brush (position constraint)
    std::vector<ConstraintSample> samples;
};

struct AttractorConstraint {
    bool orientations = true;
    std::vector<uint32_t> faces;
};

struct HeadlessConstraints {
    bool alignToBoundaries = false;
    std::vector<StrokeConstraint> strokes;
    std::vector<AttractorConstraint> attractors;
};

static inline std::string trim(const std::string &s) {
    size_t b = 0, e = s.size();
    while (b < e && std::isspace((unsigned char) s[b])) b++;
    while (e > b && std::isspace((unsigned char) s[e-1])) e--;
    return s.substr(b, e-b);
}

static HeadlessConstraints load_constraints_txt(const std::string &filename) {
    HeadlessConstraints out;
    if (filename.empty())
        return out;

    std::ifstream is(filename);
    if (!is)
        throw std::runtime_error("Could not open constraints file: " + filename);

    std::string line;
    if (!std::getline(is, line))
        return out;
    line = trim(line);
    if (line != "IMC1")
        throw std::runtime_error("Invalid constraints file (missing IMC1 header)");

    while (std::getline(is, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream ss(line);
        std::string kind;
        ss >> kind;
        if (kind == "align_to_boundaries") {
            int value = 0;
            ss >> value;
            out.alignToBoundaries = value != 0;
            continue;
        }
        if (kind == "stroke") {
            StrokeConstraint stroke;
            uint32_t count = 0;
            ss >> stroke.type >> count;
            stroke.samples.reserve(count);
            for (uint32_t i = 0; i < count; ++i) {
                std::string sampleLine;
                if (!std::getline(is, sampleLine))
                    throw std::runtime_error("Unexpected EOF while reading stroke samples");
                sampleLine = trim(sampleLine);
                if (sampleLine.empty()) { i--; continue; }
                std::istringstream ss2(sampleLine);
                ConstraintSample sample;
                ss2 >> sample.f >> sample.p.x() >> sample.p.y() >> sample.p.z();
                stroke.samples.push_back(sample);
            }
            out.strokes.push_back(stroke);
            continue;
        }
        if (kind == "attractor") {
            AttractorConstraint attr;
            int orient = 1;
            uint32_t count = 0;
            ss >> orient >> count;
            attr.orientations = orient != 0;
            attr.faces.reserve(count);
            for (uint32_t i = 0; i < count; ++i) {
                std::string faceLine;
                if (!std::getline(is, faceLine))
                    throw std::runtime_error("Unexpected EOF while reading attractor faces");
                faceLine = trim(faceLine);
                if (faceLine.empty()) { i--; continue; }
                std::istringstream ss2(faceLine);
                uint32_t f = 0;
                ss2 >> f;
                attr.faces.push_back(f);
            }
            out.attractors.push_back(attr);
            continue;
        }
        throw std::runtime_error("Unknown constraint directive: " + kind);
    }
    return out;
}

static void apply_stroke_constraints(MultiResolutionHierarchy &mRes, const HeadlessConstraints &constraints, int rosy, int posy, bool alignBoundaries) {
    if (mRes.levels() == 0)
        return;
    if (mRes.F().size() == 0)
        return;

    mRes.clearConstraints();

    if (alignBoundaries) {
        for (uint32_t i = 0; i < 3 * mRes.F().cols(); ++i) {
            if (mRes.E2E()[i] == INVALID) {
                uint32_t i0 = mRes.F()(i % 3, i / 3);
                uint32_t i1 = mRes.F()((i + 1) % 3, i / 3);
                Vector3f p0 = mRes.V().col(i0), p1 = mRes.V().col(i1);
                Vector3f edge = p1 - p0;
                if (edge.squaredNorm() > 0) {
                    edge.normalize();
                    mRes.CO().col(i0) = p0;
                    mRes.CO().col(i1) = p1;
                    mRes.CQ().col(i0) = mRes.CQ().col(i1) = edge;
                    mRes.CQw()[i0] = mRes.CQw()[i1] = mRes.COw()[i0] = mRes.COw()[i1] = 1.0f;
                }
            }
        }
    }

    const MatrixXu &F = mRes.F();
    const MatrixXf &N = mRes.N();

    for (auto const &stroke : constraints.strokes) {
        auto const &samples = stroke.samples;
        if (samples.size() < 2)
            continue;

        for (uint32_t i = 0; i < samples.size(); ++i) {
            Vector3f tangent;
            if (i == 0)
                tangent = samples[1].p - samples[0].p;
            else if (i == samples.size() - 1)
                tangent = samples[samples.size() - 1].p - samples[samples.size() - 2].p;
            else
                tangent = samples[i + 1].p - samples[i - 1].p;
            if (tangent.squaredNorm() <= RCPOVERFLOW)
                continue;
            tangent.normalize();

            uint32_t f = samples[i].f;
            if (f >= (uint32_t) F.cols())
                continue;

            for (int j = 0; j < 3; ++j) {
                uint32_t v = F(j, f);
                Vector3f tlocal = tangent;
                tlocal -= tlocal.dot(N.col(v)) * N.col(v);
                if (tlocal.squaredNorm() <= RCPOVERFLOW)
                    continue;
                tlocal.normalize();

                mRes.CQ().col(v) = tlocal;
                mRes.CQw()[v] = 1.0f;

                if (stroke.type == 1) {
                    mRes.CO().col(v) = samples[i].p;
                    mRes.COw()[v] = 1.0f;
                }
            }
        }
    }
    mRes.propagateConstraints(rosy, posy);
}

static void apply_attractor_constraints(MultiResolutionHierarchy &mRes, const HeadlessConstraints &constraints, bool orientations) {
    for (auto const &attr : constraints.attractors) {
        if (attr.orientations != orientations)
            continue;
        if (attr.faces.size() < 2)
            continue;
        // Move singularities along the face path (from end -> start, like GUI).
        for (int i = (int) attr.faces.size() - 1; i >= 1; --i) {
            uint32_t f_src = attr.faces[(size_t) i];
            uint32_t f_target = attr.faces[(size_t) i - 1];
            bool ok = orientations
                ? move_orientation_singularity(mRes, f_src, f_target)
                : move_position_singularity(mRes, f_src, f_target);
            if (!ok)
                break;
        }
    }
}

static void save_state(
    const std::string &filename,
    MultiResolutionHierarchy &mRes,
    const HeadlessOptions &options,
    const HeadlessConstraints &constraints,
    const std::set<uint32_t> &creaseIn,
    const MeshStats &stats
) {
    Serializer state;
    state.set("version", (uint32_t) 1);
    state.set("rosy", (uint32_t) options.rosy);
    state.set("posy", (uint32_t) options.posy);
    state.set("extrinsic", (uint8_t) (options.extrinsic ? 1 : 0));
    state.set("deterministic", (uint8_t) (options.deterministic ? 1 : 0));
    state.set("creaseAngle", (float) options.creaseAngle);
    state.set("smoothIter", (uint32_t) options.smoothIter);
    state.set("dominant", (uint8_t) (options.dominant ? 1 : 0));
    state.set("alignToBoundaries", (uint8_t) ((options.alignToBoundaries || constraints.alignToBoundaries) ? 1 : 0));

    state.pushPrefix("meshStats");
    state.set("weightedCenter", stats.mWeightedCenter);
    state.set("averageEdgeLength", stats.mAverageEdgeLength);
    state.set("surfaceArea", stats.mSurfaceArea);
    state.set("aabb.min", stats.mAABB.min);
    state.set("aabb.max", stats.mAABB.max);
    state.popPrefix();

    VectorXu creaseVec((uint32_t) creaseIn.size());
    uint32_t idx = 0;
    for (uint32_t v : creaseIn)
        creaseVec[idx++] = v;
    state.set("creaseIn", creaseVec);

    state.pushPrefix("constraints");
    state.set("strokeCount", (uint32_t) constraints.strokes.size());
    for (uint32_t i = 0; i < constraints.strokes.size(); ++i) {
        state.pushPrefix("stroke." + std::to_string(i));
        state.set("type", constraints.strokes[i].type);
        auto const &samples = constraints.strokes[i].samples;
        MatrixXf p(3, samples.size());
        VectorXu f(samples.size());
        for (uint32_t j = 0; j < samples.size(); ++j) {
            p.col(j) = samples[j].p;
            f[j] = samples[j].f;
        }
        state.set("p", p);
        state.set("f", f);
        state.popPrefix();
    }
    state.set("attractorCount", (uint32_t) constraints.attractors.size());
    for (uint32_t i = 0; i < constraints.attractors.size(); ++i) {
        state.pushPrefix("attractor." + std::to_string(i));
        state.set("orientations", (uint8_t) (constraints.attractors[i].orientations ? 1 : 0));
        VectorXu faces((uint32_t) constraints.attractors[i].faces.size());
        for (uint32_t j = 0; j < constraints.attractors[i].faces.size(); ++j)
            faces[j] = constraints.attractors[i].faces[j];
        state.set("faces", faces);
        state.popPrefix();
    }
    state.popPrefix();

    state.pushPrefix("mres");
    mRes.save(state);
    state.popPrefix();

    state.write(filename);
}

static void load_state(
    const std::string &filename,
    MultiResolutionHierarchy &mRes,
    HeadlessOptions &options,
    HeadlessConstraints &constraints,
    std::set<uint32_t> &creaseIn,
    MeshStats &stats
) {
    Serializer state(filename);
    uint32_t version = 0, rosy = 4, posy = 4, smoothIter = 2;
    uint8_t extrinsic = 1, deterministic = 0, dominant = 0, align = 0;
    float creaseAngle = -1;

    state.get("version", version);
    state.get("rosy", rosy);
    state.get("posy", posy);
    state.get("extrinsic", extrinsic);
    state.get("deterministic", deterministic);
    state.get("creaseAngle", creaseAngle);
    state.get("smoothIter", smoothIter);
    state.get("dominant", dominant);
    state.get("alignToBoundaries", align);

    options.rosy = (int) rosy;
    options.posy = (int) posy;
    options.extrinsic = extrinsic != 0;
    options.deterministic = deterministic != 0;
    options.creaseAngle = creaseAngle;
    options.smoothIter = (int) smoothIter;
    options.dominant = dominant != 0;
    options.alignToBoundaries = align != 0;

    state.pushPrefix("meshStats");
    state.get("weightedCenter", stats.mWeightedCenter);
    state.get("averageEdgeLength", stats.mAverageEdgeLength);
    state.get("surfaceArea", stats.mSurfaceArea);
    state.get("aabb.min", stats.mAABB.min);
    state.get("aabb.max", stats.mAABB.max);
    state.popPrefix();

    VectorXu creaseVec;
    state.get("creaseIn", creaseVec);
    creaseIn.clear();
    for (uint32_t i = 0; i < creaseVec.size(); ++i)
        creaseIn.insert(creaseVec[i]);

    state.pushPrefix("constraints");
    uint32_t strokeCount = 0;
    state.get("strokeCount", strokeCount);
    constraints.strokes.clear();
    constraints.strokes.resize(strokeCount);
    for (uint32_t i = 0; i < strokeCount; ++i) {
        state.pushPrefix("stroke." + std::to_string(i));
        uint32_t type = 0;
        MatrixXf p;
        VectorXu f;
        state.get("type", type);
        state.get("p", p);
        state.get("f", f);
        StrokeConstraint stroke;
        stroke.type = type;
        stroke.samples.resize((size_t) p.cols());
        for (uint32_t j = 0; j < stroke.samples.size(); ++j) {
            stroke.samples[j].p = p.col(j);
            stroke.samples[j].f = f[j];
        }
        constraints.strokes[i] = stroke;
        state.popPrefix();
    }
    uint32_t attractorCount = 0;
    state.get("attractorCount", attractorCount);
    constraints.attractors.clear();
    constraints.attractors.resize(attractorCount);
    for (uint32_t i = 0; i < attractorCount; ++i) {
        state.pushPrefix("attractor." + std::to_string(i));
        uint8_t orientations = 1;
        VectorXu faces;
        state.get("orientations", orientations);
        state.get("faces", faces);
        AttractorConstraint attr;
        attr.orientations = orientations != 0;
        attr.faces.resize((size_t) faces.size());
        for (uint32_t j = 0; j < faces.size(); ++j)
            attr.faces[j] = faces[j];
        constraints.attractors[i] = attr;
        state.popPrefix();
    }
    state.popPrefix();

    mRes.free();
    state.pushPrefix("mres");
    mRes.load(state);
    state.popPrefix();
}

static void build_from_mesh(
    const std::string &input,
    MultiResolutionHierarchy &mRes,
    std::set<uint32_t> &creaseIn,
    MeshStats &stats,
    const HeadlessOptions &options
) {
    MatrixXu F;
    MatrixXf V, N;
    VectorXf A;
    AdjacencyMatrix adj = nullptr;

    load_mesh_or_pointcloud(input, F, V, N);
    if (F.size() == 0)
        throw std::runtime_error("Headless parity currently supports triangle meshes only (no point clouds)");

    stats = compute_mesh_stats(F, V, options.deterministic);

    Float scale = options.scale;
    int face_count = options.faceCount;
    int vertex_count = options.vertexCount;

    if (scale < 0 && vertex_count < 0 && face_count < 0) {
        vertex_count = (int) V.cols() / 16;
    }
    if (scale > 0) {
        Float face_area = options.posy == 4 ? (scale * scale) : (std::sqrt(3.f) / 4.f * scale * scale);
        face_count = (int) (stats.mSurfaceArea / face_area);
        vertex_count = options.posy == 4 ? face_count : (face_count / 2);
    } else if (face_count > 0) {
        Float face_area = stats.mSurfaceArea / (Float) face_count;
        vertex_count = options.posy == 4 ? face_count : (face_count / 2);
        scale = options.posy == 4 ? std::sqrt(face_area) : (2 * std::sqrt(face_area * std::sqrt(1.f / 3.f)));
    } else if (vertex_count > 0) {
        face_count = options.posy == 4 ? vertex_count : (vertex_count * 2);
        Float face_area = stats.mSurfaceArea / (Float) face_count;
        scale = options.posy == 4 ? std::sqrt(face_area) : (2 * std::sqrt(face_area * std::sqrt(1.f / 3.f)));
    }

    /* Subdivide the mesh if necessary */
    VectorXu V2E, E2E;
    VectorXb boundary, nonManifold;
    if (stats.mMaximumEdgeLength * 2 > scale || stats.mMaximumEdgeLength > stats.mAverageEdgeLength * 2) {
        build_dedge(F, V, V2E, E2E, boundary, nonManifold);
        subdivide(F, V, V2E, E2E, boundary, nonManifold,
                  std::min(scale / 2, (Float) stats.mAverageEdgeLength * 2),
                  options.deterministic);
        stats = compute_mesh_stats(F, V, options.deterministic);
    }

    /* Compute a directed edge data structure */
    build_dedge(F, V, V2E, E2E, boundary, nonManifold);

    /* Compute adjacency matrix */
    adj = generate_adjacency_matrix_uniform(F, V2E, E2E, nonManifold);

    /* Compute vertex/crease normals */
    if (options.creaseAngle >= 0)
        generate_crease_normals(F, V, V2E, E2E, boundary, nonManifold, options.creaseAngle, N, creaseIn);
    else
        generate_smooth_normals(F, V, V2E, E2E, nonManifold, N);

    /* Compute dual vertex areas */
    compute_dual_vertex_areas(F, V, V2E, E2E, nonManifold, A);

    mRes.setE2E(std::move(E2E));

    mRes.setAdj(std::move(adj));
    mRes.setF(std::move(F));
    mRes.setV(std::move(V));
    mRes.setA(std::move(A));
    mRes.setN(std::move(N));
    mRes.setScale(scale);
    mRes.build(options.deterministic);
    mRes.resetSolution();
}

int headless_process(const HeadlessOptions &options) {
    HeadlessOptions local = options;
    const HeadlessConstraints constraints = load_constraints_txt(options.constraintsPath);

    Timer<> timer;

    if (options.mode == "orientation") {
        if (options.inputPath.empty() || options.outputStatePath.empty())
            throw std::runtime_error("headless orientation requires --input and --output-state");

        MultiResolutionHierarchy mRes;
        std::set<uint32_t> creaseIn;
        MeshStats stats;
        build_from_mesh(options.inputPath, mRes, creaseIn, stats, options);

        apply_stroke_constraints(mRes, constraints, options.rosy, options.posy, options.alignToBoundaries || constraints.alignToBoundaries);

        Optimizer optimizer(mRes, false);
        optimizer.setRoSy(options.rosy);
        optimizer.setPoSy(options.posy);
        optimizer.setExtrinsic(options.extrinsic);

        optimizer.optimizeOrientations(-1);
        optimizer.notify();
        optimizer.wait();
        apply_attractor_constraints(mRes, constraints, true);
        optimizer.optimizeOrientations(-1);
        optimizer.notify();
        optimizer.wait();
        optimizer.shutdown();

        save_state(options.outputStatePath, mRes, options, constraints, creaseIn, stats);
        return 0;
    }

    if (options.mode == "position") {
        if (options.statePath.empty() || options.outputStatePath.empty())
            throw std::runtime_error("headless position requires --state and --output-state");

        MultiResolutionHierarchy mRes;
        HeadlessOptions loadedOptions;
        HeadlessConstraints loadedConstraints;
        std::set<uint32_t> creaseIn;
        MeshStats stats;
        load_state(options.statePath, mRes, loadedOptions, loadedConstraints, creaseIn, stats);

        const HeadlessConstraints mergedConstraints = constraints.strokes.empty() && constraints.attractors.empty() ? loadedConstraints : constraints;
        apply_stroke_constraints(mRes, mergedConstraints, loadedOptions.rosy, loadedOptions.posy, loadedOptions.alignToBoundaries || mergedConstraints.alignToBoundaries);

        Optimizer optimizer(mRes, false);
        optimizer.setRoSy(loadedOptions.rosy);
        optimizer.setPoSy(loadedOptions.posy);
        optimizer.setExtrinsic(loadedOptions.extrinsic);

        optimizer.optimizePositions(-1);
        optimizer.notify();
        optimizer.wait();
        apply_attractor_constraints(mRes, mergedConstraints, false);
        optimizer.optimizePositions(-1);
        optimizer.notify();
        optimizer.wait();
        optimizer.shutdown();

        save_state(options.outputStatePath, mRes, loadedOptions, mergedConstraints, creaseIn, stats);
        return 0;
    }

    if (options.mode == "applyConstraints") {
        if (options.statePath.empty() || options.outputStatePath.empty())
            throw std::runtime_error("headless applyConstraints requires --state and --output-state");

        MultiResolutionHierarchy mRes;
        HeadlessOptions loadedOptions;
        HeadlessConstraints loadedConstraints;
        std::set<uint32_t> creaseIn;
        MeshStats stats;
        load_state(options.statePath, mRes, loadedOptions, loadedConstraints, creaseIn, stats);

        apply_stroke_constraints(mRes, constraints, loadedOptions.rosy, loadedOptions.posy, loadedOptions.alignToBoundaries || constraints.alignToBoundaries);
        save_state(options.outputStatePath, mRes, loadedOptions, constraints, creaseIn, stats);
        return 0;
    }

    if (options.mode == "extract") {
        if (options.statePath.empty() || options.outputMeshPath.empty())
            throw std::runtime_error("headless extract requires --state and --output-mesh");

        MultiResolutionHierarchy mRes;
        HeadlessOptions loadedOptions;
        HeadlessConstraints loadedConstraints;
        std::set<uint32_t> creaseIn;
        MeshStats stats;
        load_state(options.statePath, mRes, loadedOptions, loadedConstraints, creaseIn, stats);

        BVH *bvh = nullptr;
        if (loadedOptions.smoothIter > 0) {
            bvh = new BVH(&mRes.F(), &mRes.V(), &mRes.N(), stats.mAABB);
            bvh->build();
        }

        MatrixXf O_extr, N_extr, Nf_extr;
        std::vector<std::vector<TaggedLink>> adj_extr;
        std::set<uint32_t> creaseOut;
        extract_graph(mRes, loadedOptions.extrinsic, loadedOptions.rosy, loadedOptions.posy, adj_extr, O_extr, N_extr,
                      creaseIn, creaseOut, loadedOptions.deterministic);

        MatrixXu F_extr;
        bool pureQuad = !loadedOptions.dominant;
        extract_faces(adj_extr, O_extr, N_extr, Nf_extr, F_extr, loadedOptions.posy,
            mRes.scale(), creaseOut, true, pureQuad, bvh, loadedOptions.smoothIter);

        write_mesh(options.outputMeshPath, F_extr, O_extr, MatrixXf(), Nf_extr);
        if (bvh)
            delete bvh;
        return 0;
    }

    throw std::runtime_error("Unknown headless mode: " + options.mode);
}
