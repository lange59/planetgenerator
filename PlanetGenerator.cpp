#include "stdafx.h"
#include "PlanetGenerator.h"

PlanetGenerator::PlanetGenerator(std::string planet_name, size_t size) : planet_name(planet_name), size(size) {
    
    // Base Icosahedron, base 1

    planet_polyhedron.faces = std::vector<std::vector<uint32_t>>({
        {0, 1, 2} ,
        {0, 2, 3},
        {0, 3, 4},
        {0, 4, 5},
        {0, 5, 1},
        {1, 5, 7},
        {1, 7, 6},
        {1, 6, 2},
        {2, 6, 8},
        {2, 8, 3},
        {3, 8, 9},
        {3, 9, 4},
        {4, 9, 10},
        {4, 10, 5},
        {5, 10, 7},
        {6, 7, 11},
        {6, 11, 8},
        {7, 10, 11},
        {8, 11, 9},
        {9, 11, 10}
    });

    planet_polyhedron.vertices = std::vector<std::array<float, 3>>({
        {0, 0, 1.176 },
        {1.051, 0, 0.526},
        {0.324, 1.0, 0.525},
        {-0.851, 0.618, 0.526},
        {-0.851, -0.618, 0.526},
        {0.325, -1.0, 0.526},
        {0.851, 0.618, -0.526},
        {0.851, -0.618, -0.526},
        {-0.325, 1.0, -0.526},
        {-1.051, 0, -0.526},
        {-0.325, -1.0, -0.526},
        {0, 0, -1.176}
        });
}

void PlanetGenerator::generate() {
    
    std::cout << "Generating " << planet_name << " of size " << size;
    /*  First step apply conway operator: 
        u => subdivide
        I => inflate/spherize
        d => Dual
        k => Kis
        d => Dual
    */
    apply_subdivide(planet_polyhedron, size);

    apply_inflate(planet_polyhedron);

    apply_dual(planet_polyhedron);

    apply_kis(planet_polyhedron);

    apply_dual(planet_polyhedron);

    generate_ids();

    planet_polyhedron.to_OBJ("objs/" + planet_name + ".obj");
    //planet_polyhedron.to_DU("planets/" + planet_name + ".du");
}

/* ID Creation */

bool PlanetGenerator::is_red_valid(PolygonalMesh& mesh, OpenMesh::SmartFaceHandle face_handle,
    std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p0_p1,
    std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p0_p8,
    std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p1_p8, uint32_t size) {

    bool is_out = false;

    for (auto i = 0; i < ico_p1_p8.size(); i++) {
        if (ico_p1_p8[i].opp().next().next().opp().face() == face_handle) {
            is_out = true;
        }
    }

    return is_out;
}

bool PlanetGenerator::is_green_valid(PolygonalMesh& mesh, OpenMesh::SmartFaceHandle face_handle,
    std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p0_p1,
    std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p0_p8,
    std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p1_p8, uint32_t size) {

    bool is_out = false;

    for (auto i = 0; i < ico_p1_p8.size(); i++) {
        if (ico_p1_p8[i].prev().opp().next().next().next().opp().face() == face_handle) {
            is_out = true;
        }
    }

    return is_out;
}

bool PlanetGenerator::is_cyan_valid(PolygonalMesh& mesh, OpenMesh::SmartFaceHandle face_handle,
    std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p0_p1,
    std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p0_p8,
    std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p1_p8, uint32_t size) {

    bool is_same_as_edge_face = false;
    for (auto i = 0; i < size; i++) {
        if (ico_p1_p8[i].prev().opp().face() == face_handle) {
            is_same_as_edge_face = true;
            break;
        }
    }

    return is_same_as_edge_face;
}


std::array<float, 3> PlanetGenerator::get_center(uint32_t ix_face) {
    std::array<float, 3> fcenter = { 0, 0, 0 };
    
    auto ix_vertices = planet_polyhedron.faces[ix_face];
    // average vertex coords
    for (auto vidx = ix_vertices.begin(); vidx != ix_vertices.end(); vidx++) {
        fcenter = add(fcenter, planet_polyhedron.vertices[*vidx]);
    }
    auto tile_center = mult(1.0 / ix_vertices.size(), fcenter);

    return tile_center;
}

int32_t PlanetGenerator::get_farest_pentagon(std::vector<uint32_t> &pentagons, std::array<float,3> face_center, 
    std::map<uint32_t, int32_t> &p_index_to_ids) {
    auto distance = -1.0f;
    auto ix_f = -1;

    for (auto i = 0; i < pentagons.size(); i++) {
        auto t_center = get_center(pentagons[i]);
        auto d = calc_distance(face_center, t_center);
        if (d > distance && p_index_to_ids[pentagons[i]] == -1) {
            distance = d;
            ix_f = pentagons[i];
        }
    }

    return ix_f;
}

int32_t PlanetGenerator::get_nearest_pentagon(std::vector<uint32_t>& pentagons, std::array<float, 3>
    face_center, std::map<uint32_t, int32_t>& p_index_to_ids) {
    auto distance = 100000000000.0f;
    auto ix_f = -1;

    for (auto i = 0; i < pentagons.size(); i++) {
        auto t_center = get_center(pentagons[i]);
        auto d = calc_distance(face_center, t_center);
        if (d < distance && p_index_to_ids[pentagons[i]] == -1) {
            distance = d;
            ix_f = pentagons[i];
        }
    }

    return ix_f;
}


std::map<uint32_t, int32_t> PlanetGenerator::constructs_pentagons() {
    std::map<uint32_t, int32_t> p_index_to_ids;
    for (auto i = 0; i < planet_polyhedron.faces.size(); i++) {
        if (planet_polyhedron.faces[i].size() == 5) {
            p_index_to_ids[i] = -1;
        }
    }

    assert(p_index_to_ids.size() == 12);

    std::vector<uint32_t> pentagons;
    for (auto pair : p_index_to_ids) {
        pentagons.push_back(pair.first);
    }

    // The first face is ALWAYS a pentagon
    auto ix_tile_0 = pentagons[0];
    p_index_to_ids[ix_tile_0] = 0;

    auto ix_tile_2 = get_farest_pentagon(pentagons, get_center(ix_tile_0), p_index_to_ids);
    p_index_to_ids[ix_tile_2] = 2;

    auto ix_tile_1 = get_nearest_pentagon(pentagons, get_center(ix_tile_0), p_index_to_ids);
    p_index_to_ids[ix_tile_1] = 1;

    auto ix_tile_3 = get_farest_pentagon(pentagons, get_center(ix_tile_1), p_index_to_ids);
    p_index_to_ids[ix_tile_3] = 3;
    
    auto ix_tile_8 = -1;
    {
        auto t0_center = get_center(ix_tile_0);
        auto t1_center = get_center(ix_tile_1);
        auto distance = 100000000000.0f;
        auto ix_f = -1;
        for (auto i = 0; i < pentagons.size(); i++) {
            auto t_center = get_center(pentagons[i]);
            auto d = calc_distance(t0_center, t_center) + calc_distance(t1_center, t_center);

            if (d < distance && p_index_to_ids[pentagons[i]] == -1) {
                distance = d;
                ix_f = pentagons[i];
            }
        }
        ix_tile_8 = ix_f;
    }
    p_index_to_ids[ix_tile_8] = 8;

    auto ix_tile_10 = get_farest_pentagon(pentagons, get_center(ix_tile_8), p_index_to_ids);
    p_index_to_ids[ix_tile_10] = 10;

    auto ix_tile_11 = -1;
    {
        auto t0_center = get_center(ix_tile_0);
        auto t1_center = get_center(ix_tile_1);
        auto distance = 100000000000.0f;
        auto ix_f = -1;
        for (auto i = 0; i < pentagons.size(); i++) {
            auto t_center = get_center(pentagons[i]);
            auto d = calc_distance(t0_center, t_center) + calc_distance(t1_center, t_center);

            if (d < distance && p_index_to_ids[pentagons[i]] == -1) {
                distance = d;
                ix_f = pentagons[i];
            }
        }
        ix_tile_11 = ix_f;
    }
    p_index_to_ids[ix_tile_11] = 11;

    auto ix_tile_9 = get_farest_pentagon(pentagons, get_center(ix_tile_11), p_index_to_ids);
    p_index_to_ids[ix_tile_9] = 9;

    auto ix_tile_5 = -1;
    {
        auto t0_center = get_center(ix_tile_0);
        auto t10_center = get_center(ix_tile_10);
        auto t11_center = get_center(ix_tile_11);
        auto distance = 100000000000.0f;
        auto ix_f = -1;
        for (auto i = 0; i < pentagons.size(); i++) {
            auto t_center = get_center(pentagons[i]);
            auto d = calc_distance(t0_center, t_center) + calc_distance(t10_center, t_center) +
                calc_distance(t11_center, t_center);

            if (d < distance && p_index_to_ids[pentagons[i]] == -1) {
                distance = d;
                ix_f = pentagons[i];
            }
        }
        ix_tile_5 = ix_f;
    }
    p_index_to_ids[ix_tile_5] = 5;

    auto ix_tile_7 = get_farest_pentagon(pentagons, get_center(ix_tile_5), p_index_to_ids);
    p_index_to_ids[ix_tile_7] = 7;

    auto ix_tile_4 = -1;
    {
        auto t0_center = get_center(ix_tile_0);
        auto t5_center = get_center(ix_tile_5);
        auto distance = 100000000000.0f;
        auto ix_f = -1;
        for (auto i = 0; i < pentagons.size(); i++) {
            auto t_center = get_center(pentagons[i]);
            auto d = calc_distance(t0_center, t_center) + calc_distance(t5_center, t_center);

            if (d < distance && p_index_to_ids[pentagons[i]] == -1) {
                distance = d;
                ix_f = pentagons[i];
            }
        }
        ix_tile_4 = ix_f;
    }
    p_index_to_ids[ix_tile_4] = 4;

    auto ix_tile_6 = get_nearest_pentagon(pentagons, get_center(ix_tile_4), p_index_to_ids);
    p_index_to_ids[ix_tile_6] = 6;

    return p_index_to_ids;
}

void PlanetGenerator::generate_ids() {
    /* Pentagon Ids */
    std::map<uint32_t, int32_t> p_index_to_ids = constructs_pentagons();
    
    /* BEGIN CREATION MESH */
    PolygonalMesh mesh;
    PolygonalMesh::VertexHandle* vhandle = new PolygonalMesh::VertexHandle[planet_polyhedron.vertices.size()];

    std::map<int32_t, PolygonalMesh::FaceHandle> id_to_face_map;

    for (auto i = 0; i < planet_polyhedron.vertices.size(); i++) {
        vhandle[i] = mesh.add_vertex(PolygonalMesh::Point(planet_polyhedron.vertices[i][0], 
            planet_polyhedron.vertices[i][1], planet_polyhedron.vertices[i][2]));
    }

    std::vector<PolygonalMesh::VertexHandle>  face_vhandles;
    for (auto i = 0; i < planet_polyhedron.faces.size(); i++) {
        face_vhandles.clear();

        for (auto j = 0; j < planet_polyhedron.faces[i].size(); j++) {
            face_vhandles.push_back(vhandle[planet_polyhedron.faces[i][j]]);
        }
        auto face_handle = mesh.add_face(face_vhandles);
        auto& f = mesh.data(face_handle);

        // ADD PENTA TO id_to_face_map 
        if (planet_polyhedron.faces[i].size() == 5) {
            f.set_tile_id(p_index_to_ids[i]);
            id_to_face_map[p_index_to_ids[i]] = face_handle;
        }
        
    }
    delete[] vhandle;

    /* END CREATION MESH */


    // CREATE ALL PENTA-PENTA EDGES
    std::map<std::string, std::vector<OpenMesh::SmartHalfedgeHandle>> ico_edges;

    for (auto i = 0; i < 12; i++) {
        auto start_face_handle = id_to_face_map[i];
        auto edges = mesh.fe_iter(start_face_handle);

        uint32_t start_pentagon_tile_id = mesh.data(start_face_handle).tile_id();

        PolygonalMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(start_face_handle);

        // CORNERS

        for (; fh_it.is_valid(); ++fh_it) {
            // GET PENTA CORNER EDGE
            auto _corner_he = fh_it->opp().prev();
            auto _corner_face = _corner_he.prev().opp().face();

            std::vector<OpenMesh::SmartHalfedgeHandle> corners_he;

            auto i = 0;
            while (i < size) {
                corners_he.push_back(_corner_he);

                // ITERATE ON THE ICO EDGE TO SIZE-1
                _corner_he = _corner_he.prev().opp().next().next().opp().prev();
                _corner_face = _corner_he.prev().opp().face();

                i++;
            }

            // GET PENTA AT THE END
            _corner_face = _corner_he.next().opp().face();

            uint32_t end_pentagon_tile_id = mesh.data(_corner_face).tile_id();
            assert(end_pentagon_tile_id >= 0);

            ico_edges[format("p%d-p%d", start_pentagon_tile_id, end_pentagon_tile_id)] = corners_he;
        }
    }

    // 30 EDGES ON ICO, 0-1 // 1-0 SO 30 * 2 EDGES IN MAP
    assert(ico_edges.size() == 60);

    std::vector < std::array<uint32_t, 3>> ico_triangles;
    ico_triangles.push_back({ 0,1,8 });
    ico_triangles.push_back({ 0,1,11 });
    ico_triangles.push_back({ 0,5,11 });
    ico_triangles.push_back({ 0,4,5 });
    ico_triangles.push_back({ 0,4,8 });

    ico_triangles.push_back({ 4,8,9 });
    ico_triangles.push_back({ 7,8,9 });
    ico_triangles.push_back({ 1,7,8 });
    ico_triangles.push_back({ 1,6,7 });
    ico_triangles.push_back({ 1,6,11 });

    ico_triangles.push_back({ 6,10,11 });
    ico_triangles.push_back({ 5,10,11 });
    ico_triangles.push_back({ 3,5,10 });
    ico_triangles.push_back({ 3,4,5 });
    ico_triangles.push_back({ 3,4,9 });

    ico_triangles.push_back({ 2,3,9 });
    ico_triangles.push_back({ 2,7,9 });
    ico_triangles.push_back({ 2,6,7 });
    ico_triangles.push_back({ 2,6,10 });
    ico_triangles.push_back({ 2,3,10 });

    assert(ico_triangles.size() == 20);

    // CREATE PLANET ICO-TRIANGLE
    std::list<std::string> edges_list;
    for (auto i = 0; i < ico_triangles.size(); i++) {
        auto tri = ico_triangles[i];
        edges_list.push_back(format("p%d-p%d", tri[0], tri[1]));
        edges_list.push_back(format("p%d-p%d", tri[0], tri[2]));
        edges_list.push_back(format("p%d-p%d", tri[1], tri[2]));
    }
    assert(edges_list.size() == 60);

    uint32_t tile_id = 12;

    // ID for ico edges
    for (auto el_it = edges_list.begin(); el_it != edges_list.end(); el_it++) {
        auto corners_he = ico_edges[*el_it];
        auto face_handle = corners_he[0].prev().opp().face();
        // edge have id yet
        if (mesh.data(face_handle).tile_id() == 0) {
            for (auto i = 0; i < size - 1; i++) {
                auto face_handle = corners_he[i].prev().opp().face();
                mesh.data(face_handle).set_tile_id(tile_id);
                id_to_face_map[tile_id] = face_handle;
                tile_id++;
            }
        }
    }

    // FACE ON ICO EDGES = (SIZE - 1) * 30 + 12
    assert((tile_id) == ((size - 1) * 30 + 12));
    assert(ico_edges.size() == 60); // REASSERT, IF ERROR WILL BE 61

    std::vector<std::string> tri_edges_list;


    // TRIANGLE 1
    // tile_id = 192
    tri_edges_list.push_back("p0-p1");
    tri_edges_list.push_back("p0-p8");
    tri_edges_list.push_back("p1-p8");

    // TRIANGLE 2
    // tile_id = 256
    tri_edges_list.push_back("p1-p0");
    tri_edges_list.push_back("p1-p11");
    tri_edges_list.push_back("p0-p11");

    // TRIANGLE 3
    // tile_id = 320
    tri_edges_list.push_back("p0-p5");
    tri_edges_list.push_back("p0-p11");
    tri_edges_list.push_back("p5-p11");

    // TRIANGLE 4
    // tile_id = 384
    tri_edges_list.push_back("p0-p4");
    tri_edges_list.push_back("p0-p5");
    tri_edges_list.push_back("p4-p5");

    // TRIANGLE 5
    // tile_id = 448
    tri_edges_list.push_back("p0-p8");
    tri_edges_list.push_back("p0-p4");
    tri_edges_list.push_back("p8-p4");

    // TRIANGLE 6
    // tile_id = 512
    tri_edges_list.push_back("p4-p8");
    tri_edges_list.push_back("p4-p9");
    tri_edges_list.push_back("p8-p9");

    // TRIANGLE 7
    // tile_id = 576
    tri_edges_list.push_back("p7-p9");
    tri_edges_list.push_back("p7-p8");
    tri_edges_list.push_back("p9-p8");

    // TRIANGLE 8
    // tile_id = 640
    tri_edges_list.push_back("p1-p7");
    tri_edges_list.push_back("p1-p8");
    tri_edges_list.push_back("p7-p8");

    // TRIANGLE 9
    // tile_id = 704
    tri_edges_list.push_back("p1-p6");
    tri_edges_list.push_back("p1-p7");
    tri_edges_list.push_back("p6-p7");

    // TRIANGLE 10
    // tile_id = 768
    tri_edges_list.push_back("p1-p11");
    tri_edges_list.push_back("p1-p6");
    tri_edges_list.push_back("p11-p6");

    // TRIANGLE 11
    // tile_id = 832
    tri_edges_list.push_back("p6-p11");
    tri_edges_list.push_back("p6-p10");
    tri_edges_list.push_back("p11-p10");

    // TRIANGLE 12
    // tile_id = 896
    tri_edges_list.push_back("p10-p11");
    tri_edges_list.push_back("p5-p10");
    tri_edges_list.push_back("p11-p5");

    // TRIANGLE 13
    // tile_id = 960
    tri_edges_list.push_back("p3-p10");
    tri_edges_list.push_back("p3-p5");
    tri_edges_list.push_back("p10-p5");

    // TRIANGLE 14
    // tile_id = 1024
    tri_edges_list.push_back("p3-p5");
    tri_edges_list.push_back("p3-p4");
    tri_edges_list.push_back("p5-p4");

    // TRIANGLE 15
    // tile_id = 1088
    tri_edges_list.push_back("p3-p4");
    tri_edges_list.push_back("p3-p9");
    tri_edges_list.push_back("p4-p9");

    // TRIANGLE 16
    // tile_id = 1152
    tri_edges_list.push_back("p2-p3");
    tri_edges_list.push_back("p2-p9");
    tri_edges_list.push_back("p3-p9");

    // TRIANGLE 17
    // tile_id = 1216
    tri_edges_list.push_back("p2-p9");
    tri_edges_list.push_back("p2-p7");
    tri_edges_list.push_back("p9-p7");

    // TRIANGLE 18
    // tile_id = 1280
    tri_edges_list.push_back("p2-p7");
    tri_edges_list.push_back("p2-p6");
    tri_edges_list.push_back("p7-p6");

    // TRIANGLE 19
    // tile_id = 1344
    tri_edges_list.push_back("p2-p6");
    tri_edges_list.push_back("p2-p10");
    tri_edges_list.push_back("p6-p10");

    // TRIANGLE 20
    // tile_id = 1408
    tri_edges_list.push_back("p2-p10");
    tri_edges_list.push_back("p2-p3");
    tri_edges_list.push_back("p10-p3");

    assert(tri_edges_list.size() == 60);

    std::cout << tri_edges_list.size() / 3 << std::endl;

    for (auto i = 0; i < tri_edges_list.size() / 3; i++) {

        auto begin_tile_id = tile_id;

        auto ico_p0_p1 = ico_edges[tri_edges_list[i * 3]];
        auto ico_p0_p8 = ico_edges[tri_edges_list[i * 3 + 1]];
        auto ico_p1_p8 = ico_edges[tri_edges_list[i * 3 + 2]];

        for (auto i = 0; i < size; i++) {
            auto face_handle = ico_p0_p1[i].face();
            mesh.data(face_handle).set_tile_id(tile_id);
            id_to_face_map[tile_id] = face_handle;
            tile_id++;

            if (i != (size - 1)) {
                // DRAW THE LINE OF RED
                auto r_he = ico_p0_p1[i].prev().prev().opp().next().next().opp();
                auto r_f = r_he.face();

                auto r_max_iters = size;

                while (!is_red_valid(mesh, r_f, ico_p0_p1, ico_p0_p8, ico_p1_p8, size)) {

                    mesh.data(r_f).set_tile_id(tile_id);
                    id_to_face_map[tile_id] = r_f;
                    tile_id++;

                    r_max_iters--;
                    if (r_max_iters == 0) {
                        break;
                    }

                    r_he = r_he.prev().prev().opp().next().next().opp();
                    r_f = r_he.face();
                }



                // DRAW THE LINE OF GREEN
                auto g_he = ico_p0_p1[i].prev().prev().opp();
                auto g_f = g_he.face();
                auto g_max_iters = size;

                while (!is_green_valid(mesh, g_f, ico_p0_p1, ico_p0_p8, ico_p1_p8, size)) {
                    mesh.data(g_f).set_tile_id(tile_id);
                    id_to_face_map[tile_id] = g_f;
                    tile_id++;

                    g_max_iters--;
                    if (g_max_iters == 0) {
                        break;
                    }

                    g_he = g_he.next().next().opp().prev().prev().opp();
                    g_f = g_he.face();
                }


                // DRAW THE LINE OF CYAN
                auto c_e = ico_p0_p1[i].prev().prev().opp().prev().prev().prev().opp();
                auto c_f = c_e.face();
                auto c_max_iters = size;

                while (!is_cyan_valid(mesh, c_f, ico_p0_p1, ico_p0_p8, ico_p1_p8, size)) {
                    mesh.data(c_f).set_tile_id(tile_id);
                    id_to_face_map[tile_id] = c_f;
                    tile_id++;

                    c_max_iters--;
                    if (c_max_iters == 0) {
                        break;
                    }

                    c_e = c_e.next().next().opp().prev().prev().opp();
                    c_f = c_e.face();
                }
            }
        }

        std::cout << "Triangle " << i + 1 << std::endl;
        std::cout << "IDs: " << begin_tile_id << " - " << tile_id - 1 << std::endl;
    }

    /*
    // PRINTING PENTA ADJS FACES IDs
    for (auto i = 0; i < 12; i++) {
        auto penta_face_handle = id_to_face_map[i];
        auto ff_iter = mesh.ff_iter(penta_face_handle);

        std::cout << "P" << i;
        for (; ff_iter.is_valid(); ++ff_iter) {
            std::cout << " " << mesh.data(*ff_iter).tile_id();
        }
        std::cout << std::endl;
    }
    */
    Polyhedron new_planet_polyhedron;

    auto ix_v = 0;
    std::map<std::string, uint32_t> uniqVs;

    for (auto id = 0; id < planet_polyhedron.faces.size(); id++) {
        PolygonalMesh::FaceHandle face_handle = id_to_face_map[id];
        int32_t tile_id = mesh.data(face_handle).tile_id();
        assert(tile_id == id);

        auto fv_it = mesh.fv_iter(face_handle);
        std::vector<uint32_t> face_index;

        for (; fv_it.is_valid(); ++fv_it) {
            auto p = mesh.point(*fv_it);
            auto v = std::array<float, 3>({ p[0],p[1],p[2] });
            auto vK = format("v%f-v%f-v%f", v[0], v[1], v[2]);
            if (uniqVs.count(vK) == 0) {
                uniqVs[vK] = ix_v;
                new_planet_polyhedron.vertices.push_back(v);
                face_index.push_back(ix_v);
                ix_v++;
            }
            else {
                face_index.push_back(uniqVs[vK]);
            }
        }
        new_planet_polyhedron.faces.push_back(face_index);
    }

    planet_polyhedron.faces = new_planet_polyhedron.faces;
    planet_polyhedron.vertices = new_planet_polyhedron.vertices;
}


