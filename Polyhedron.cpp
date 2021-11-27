#include "stdafx.h"
#include "Polyhedron.h"


Polyhedron::Polyhedron()
{
}


std::vector< std::array<float, 3>> Polyhedron::centers() {
    std::vector< std::array<float, 3>> centers;

    for (auto face = faces.begin(); face != faces.end(); face++) {
        std::array<float, 3> fcenter = { 0, 0, 0 };
        // average vertex coords
        for (auto vidx = face->begin(); vidx != face->end(); vidx++) {
            fcenter = add(fcenter, vertices[*vidx]);
        }
        centers.push_back(mult(1.0 / face->size(), fcenter));
    }
    // return face-ordered array of centroids
    return centers;
}

std::vector<std::array<float, 3>> Polyhedron::normals() {
    if (_normals.size() == 0) {
        for (auto it = faces.begin(); it != faces.end(); it++) {
            std::vector<std::array<float, 3>> faces_vertices;
            for (auto i = it->begin(); i != it->end(); i++) {
                faces_vertices.push_back(vertices[*i]);
            }

            _normals.push_back(normal(faces_vertices));
        }
    }
    
    return _normals;
}


void Polyhedron::to_OBJ(std::string output_path)
{
    std::ofstream output_file;
    output_file.open(output_path);

    for (auto it = vertices.begin(); it != vertices.end(); it++) {
        output_file << format("v %f %f %f", (*it)[0], (*it)[1], (*it)[2]) << std::endl;
    }

    for (auto it = faces.begin(); it != faces.end(); it++) {
        std::vector<std::array<float, 3>> faces_vertices;
        for (auto i = it->begin(); i != it->end(); i++) {
            faces_vertices.push_back(vertices[*i]);
        }
        
        auto norm = normal(faces_vertices);

        output_file << format("vn %f %f %f", norm[0], norm[1], norm[2]) << std::endl;
    }

    for (auto i = 0; i < faces.size(); i++) {
        auto f = faces[i];
        output_file << "f ";
        for (auto v = f.begin(); v != f.end(); v++) {
            output_file << format("%d//%d ", (*v)+1, i+1);
        }
        output_file << std::endl;
    }


    output_file.close();
}


template <typename T, class... StreamArgs>
inline std::basic_ostream<StreamArgs...>&
operator <= (std::basic_ostream<StreamArgs...>& out, T const& data) {
    out.write(reinterpret_cast<char const*>(&data), sizeof(T));
    return out;
}

void Polyhedron::to_DU(std::string output_path) {
    // DU FORMAT v 0.1
    // HEADERS
    //
    // UINT32_T MAGIC_NUMBER 0x587478
    // UINT32_T FORMAT_VERSION
    // UINT64_T CREATE_TIMESTAMP
    // UINT32_T VERTICES_COUNT
    // UINT32_T FACES_COUNT
    // N VERTICES OF FLOAT/FLOAT/FLOAT
    // N FACES UINT8_T / UINT32_T / UINT32_T / UINT32_T / ... (UINT8_T count)

    std::ofstream file_stream;
    file_stream.open(output_path, std::ios::binary);

    uint32_t magic_number = 2021;
    uint32_t format_version = 1;
    uint64_t create_timestamp = 0;

    uint32_t vertices_count = vertices.size();
    uint32_t faces_count = faces.size();

    file_stream <= magic_number <= format_version <= create_timestamp;
    file_stream <= vertices_count <= faces_count;


    for (auto v_it = vertices.begin(); v_it != vertices.end(); v_it++) {
        file_stream <= v_it->at(0) <= v_it->at(1) <= v_it->at(2);
    }

    for (auto f_it = faces.begin(); f_it != faces.end(); f_it++) {
        uint8_t faces_vertices_count = f_it->size();
        file_stream <= faces_vertices_count;
        for (auto i = 0; i < faces_vertices_count; i++) {
            file_stream <= f_it->at(i);
        }
    }


    file_stream.close();
}

Polyhedron Polyhedron::from_OBJ(std::string input_path) {
    Polyhedron polyhedron;

    std::ifstream input_file;
    input_file.open(input_path);

    std::string line;
    while (std::getline(input_file, line))
    {
        std::istringstream iss(line);

        std::string lineType;
        iss >> lineType;

        if (lineType == "v") {
            float x, y, z;
            iss >> x >> y >> z;
            polyhedron.vertices.push_back(std::array<float, 3>({ x,y,z }));
        }
        else if (lineType == "vn") {
            float x, y, z;
            iss >> x >> y >> z;
            polyhedron._normals.push_back(std::array<float, 3>({ x,y,z }));
        }
        else if (lineType == "f") {
            //f 2877//1457 2874//1457 2860//1457 2856//1457 2522//1457 2521//1457
            std::string fStr;
            
            std::vector<uint32_t> ix_fs;
            while (iss >> fStr) {
                std::istringstream iss_face(fStr);
                uint32_t ix_f, ix_n;
                std::string t;

                iss_face >> ix_f >> t >> ix_n;
                ix_fs.push_back(ix_f-1);
            }  
            polyhedron.faces.push_back(ix_fs);
        }
    }

    input_file.close();

    return polyhedron;
}

Polyhedron Polyhedron::from_DU(std::string input_path) {
    Polyhedron polyhedron;
    
    std::ifstream file_stream;
    file_stream.open(input_path, std::ios::binary);

    uint32_t magic_number = 0;
    uint32_t format_version = 0;
    uint64_t create_timestamp = 0;

    uint32_t vertices_count = 0;
    uint32_t faces_count = 0;

    file_stream.read(reinterpret_cast<char *>(&magic_number), sizeof(magic_number));
    file_stream.read(reinterpret_cast<char*>(&format_version), sizeof(format_version));
    file_stream.read(reinterpret_cast<char*>(&create_timestamp), sizeof(create_timestamp));

    file_stream.read(reinterpret_cast<char*>(&vertices_count), sizeof(vertices_count));
    file_stream.read(reinterpret_cast<char*>(&faces_count), sizeof(faces_count));

    std::cout << "magic_number: " << magic_number << " format_version: "
        << format_version << " create_timestamp: " << create_timestamp << std::endl;
    std::cout << "vertices_count: " << vertices_count << " faces_count: "
        << faces_count << std::endl;

    polyhedron.vertices.resize(vertices_count);
    float v1, v2, v3;

    for (auto i = 0; i < polyhedron.vertices.size(); i++) {
        file_stream.read(reinterpret_cast<char*>(&v1), sizeof(float));
        file_stream.read(reinterpret_cast<char*>(&v2), sizeof(float));
        file_stream.read(reinterpret_cast<char*>(&v3), sizeof(float));

        polyhedron.vertices[i] = std::array<float, 3>({ v1,v2,v3 });
    }
    
    polyhedron.faces.resize(faces_count);
    for (auto i = 0; i < polyhedron.faces.size(); i++) {
        uint8_t faces_count = 0;
        file_stream.read(reinterpret_cast<char*>(&faces_count), sizeof(uint8_t));
        std::vector<uint32_t> ix_vertices(faces_count);
        
        file_stream.read(reinterpret_cast<char*>(ix_vertices.data()), sizeof(uint32_t)* faces_count);
        polyhedron.faces[i] = ix_vertices;
    }

    /*
    for (auto v_it = vertices.begin(); v_it != vertices.end(); v_it++) {
        file_stream <= v_it->at(0) <= v_it->at(1) <= v_it->at(2);
    }

    for (auto f_it = faces.begin(); f_it != faces.end(); f_it++) {
        uint8_t faces_count = f_it->size();
        file_stream <= faces_count;
        for (auto i = 0; i < faces_count; i++) {
            file_stream <= f_it->at(i);
        }
    }
    */

    file_stream.close();

    return polyhedron;
}






Polyflag::Polyflag()
{
}

void Polyflag::newV(std::string vertName, std::array<float, 3> coordinates) {
    if (vertidxs.count(vertName) == 0) {
        vertidxs[vertName] = 0;
        vertices[vertName] = coordinates;
    }
}

void Polyflag::newFlag(std::string faceName, std::string vertName1, std::string vertName2) {
    if (flags.count(faceName) == 0) {
        flags[faceName] = std::map<std::string, std::string>();
    }
    flags[faceName][vertName1] = vertName2;
}

Polyhedron Polyflag::topoly() {
    Polyhedron poly;

    uint32_t ctr = 0; // first number the vertices
    for (auto i = vertidxs.begin(); i != vertidxs.end(); i++) {
        auto idx = format("%s", i->first.c_str());
        auto vtx = vertices[idx];

        poly.vertices.push_back(vtx);
        vertidxs[idx] = ctr;
        ctr++;
    }


    ctr = 0;
    for (auto i = flags.begin(); i != flags.end(); i++) {
        auto idx = format("%s", i->first.c_str());
        auto face = flags[idx];

        std::string v0;

        poly.faces.push_back(std::vector<uint32_t>()); // new face
        // grab _any_ vertex as starting point
        for (auto j = face.begin(); j != face.end(); j++) {
            v0 = j->second;
            break;
        }
        // build face out of all the edge relations in the flag assoc array
        std::string v = v0; // v moves around face

        poly.faces[ctr].push_back(vertidxs[v]); //record index
        v = flags[idx][v]; // goto next vertex
        auto faceCTR = 0;
        while (v != v0) { // loop until back to start
            poly.faces[ctr].push_back(vertidxs[v]);
            v = flags[idx][v];
            faceCTR++;
        }
        ctr++;
    }

    return poly;
}




void apply_subdivide(Polyhedron &poly, uint32_t n)
{
    // Can only work with triangles
    for (auto fn = 0; fn < poly.faces.size(); fn++) {
        if (poly.faces[fn].size() != 3) {
            return;
        }
    }

    // Calculate redundant set of new vertices for subdivided mesh.
    std::vector<std::array<float, 3>> newVs;
    std::map<std::string, uint32_t> vmap;
    uint32_t pos = 0;

    for (auto fn = 0; fn < poly.faces.size(); fn++) {
        const std::vector<uint32_t> f = poly.faces[fn];

        std::array<float, 3> v1 = poly.vertices[f[0]];
        std::array<float, 3> v2 = poly.vertices[f[1]];
        std::array<float, 3> v3 = poly.vertices[f[2]];
        
        std::array<float, 3> v21 = sub(v2, v1);
        std::array<float, 3> v31 = sub(v3, v1);
        
        for (auto i = 0; i <= n; i++) {
            for (auto j = 0; j + i <= n; j++) {
                auto v = add(add(v1, mult(i * 1.0 / n, v21)), mult(j * 1.0 / n, v31));
                vmap[format("v${%d} - ${ %d } - ${ %d }", fn, i, j)] = pos++;
                newVs.push_back(v);
            }
        }
    }

    const float EPSILON_CLOSE = 1.0e-8;
    std::vector<std::array<float, 3>> uniqVs;
    uint32_t newpos = 0;
    std::map<uint32_t, uint32_t> uniqmap;

    uint32_t i = 0;

    for (auto it = newVs.begin(); it != newVs.end(); it++, i++) {  
        if (uniqmap.count(i)>0) { continue; } // already mapped

        uniqmap[i] = newpos;
        uniqVs.push_back(*it);
        
        for (auto j = i + 1; j < newVs.size(); j++) {
            std::array<float, 3> w = newVs[j];
            auto c = mag(sub(*it, w));
             if ( c < 0.000001) {
                uniqmap[j] = newpos;
            }
        }
        newpos++;
    }

    std::vector<std::vector<uint32_t>> faces;

    for (auto fn = 0; fn < poly.faces.size(); fn++) {
        for (auto i = 0; i < n; i++) {
            for (auto j = 0; j + i < n; j++) {
                faces.push_back({
                    uniqmap[vmap[format("v${%d} - ${ %d } - ${ %d }", fn, i, j)]],
                    uniqmap[vmap[format("v${%d} - ${ %d } - ${ %d }", fn, i + 1, j)]],
                    uniqmap[vmap[format("v${%d} - ${ %d } - ${ %d }", fn, i, j + 1)]]
                    });
            }
        }

        for (auto i = 1; i < n; i++) {
            for (auto j = 0; j + i < n; j++) {
                faces.push_back({
                        uniqmap[vmap[format("v${%d} - ${ %d } - ${ %d }", fn, i, j)]],
                        uniqmap[vmap[format("v${%d} - ${ %d } - ${ %d }", fn, i , j + 1)]],
                        uniqmap[vmap[format("v${%d} - ${ %d } - ${ %d }", fn, i - 1, j + 1)]]
                    });
            }
        }
    }

    poly.faces = faces;
    poly.vertices = uniqVs;
}

void apply_dual(Polyhedron& poly)
{
    Polyflag flag;

    std::vector<std::map<std::string, std::string>> face; // make table of face as fn of edge
    face.resize(poly.vertices.size());
    
    /*
    for (auto i = 0; i <= poly.vertices.size() - 1; i++) {
        face[i] = std::map<std::string, std::string>();
    } // create empty associative table
    */

    for (auto i = 0; i < poly.faces.size(); i++) {
        std::vector<uint32_t> f = poly.faces[i];
        uint32_t v1 = f[f.size() - 1]; //previous vertex
        for (auto v2 = f.begin(); v2 != f.end(); v2++) {
            // THIS ASSUMES that no 2 faces that share an edge share it in the same orientation!
            // which of course never happens for proper manifold meshes, so get your meshes right.
            face[v1][format("v%d", *v2)] = format("%d", i);
            v1 = *v2;
        }
    } // current becomes previous

    std::vector< std::array<float, 3>> centers = poly.centers();

    for (auto i = 0; i <= poly.faces.size() - 1; i++) {
        flag.newV(format("%d", i),centers[i]);
    }


    for (auto i = 0; i < poly.faces.size(); i++) {
        std::vector<uint32_t> f = poly.faces[i];
        uint32_t v1 = f[f.size() - 1]; //previous vertex
        for (auto v2 = f.begin(); v2 != f.end(); v2++) {
            flag.newFlag(format("%d",v1), face[*v2][format("v%d", v1)], format("%d", i));
            v1 = *v2;
        } // current becomes previous
    }

    Polyhedron dpoly = flag.topoly(); // build topological dual from flags

    /*
    // match F index ordering to V index ordering on dual
    std::vector<std::vector<uint32_t>> sortF;
    sortF.resize(dpoly.faces.size());

    for (auto f = dpoly.faces.begin(); f != dpoly.faces.end(); f++) {
        auto v1 = (*f)[0];
        auto v2 = (*f)[1];
        auto v3 = (*f)[2];

        uint32_t k = intersect(poly.faces[v1], poly.faces[v2], poly.faces[v3]);
        sortF[k] = *f;
    }
    */

    poly.vertices = dpoly.vertices;
    poly.faces = dpoly.faces;

} 

void apply_kis(Polyhedron& poly, uint32_t n, float apexdist)
{
    Polyflag flag;

    for (auto i = 0; i < poly.vertices.size(); i++) {
        // each old vertex is a new vertex
        auto p = poly.vertices[i];
        flag.newV(format("v%d", i), p);
    }


    std::vector< std::array<float, 3>> normals = poly.normals();
    std::vector< std::array<float, 3>> centers = poly.centers();

    bool foundAny = false;

    for (auto idx_face = 0; idx_face < poly.faces.size(); idx_face++) {
        auto f = poly.faces[idx_face];
        std::string v1 = format("v%d", f[f.size() - 1]);

        for (auto idx_v = 0; idx_v < f.size(); idx_v++) {
            std::string v2 = format("v%d", f[idx_v]);

            if (f.size() == n || n == 0) {
                foundAny = true;
                std::string apex = format("apex%d", idx_face);
                std::string fname = format("%d%s", idx_face, v1.c_str());
                flag.newV(apex, add(centers[idx_face], mult(apexdist, normals[idx_face])));
                flag.newFlag(fname, v1, v2); // the old edge of original face
                flag.newFlag(fname, v2, apex); // up to apex of pyramid
                flag.newFlag(fname, apex, v1); // and back down again
            }
            else {
                flag.newFlag(format("%d", idx_face), v1, v2);  // same old flag, if non-n
            }
            // current becomes previous
            v1 = v2;
        }
    }

    auto newpoly = flag.topoly();

    poly.faces = newpoly.faces;
    poly.vertices = newpoly.vertices;
}

void apply_inflate(Polyhedron& poly){
    
    std::vector<std::array<float, 3>> s;

    for (auto i = 0; i < poly.vertices.size(); i++) {
        auto vertex = poly.vertices[i];
        auto v = unit(vertex);
        s.push_back(v);
    }

    poly.vertices = s;
}


//missing string printf
//this is safe and convenient but not exactly efficient
std::string format(const char* fmt, ...) {
    int size = 512;
    char* buffer = 0;
    buffer = new char[size];
    va_list vl;
    va_start(vl, fmt);
    int nsize = vsnprintf(buffer, size, fmt, vl);
    if (size <= nsize) { //fail delete buffer and try again
        delete[] buffer;
        buffer = 0;
        buffer = new char[nsize + 1]; //+1 for /0
        nsize = vsnprintf(buffer, size, fmt, vl);
    }
    std::string ret(buffer);
    va_end(vl);
    delete[] buffer;
    return ret;
}

// 3d scalar multiplication
std::array<float, 3> mult(float c, std::array<float, 3> vec) {
    return std::array<float, 3>({ c * vec[0], c * vec[1], c * vec[2] });
}

// 3d vector addition
std::array<float, 3> add(std::array<float, 3> vec1, std::array<float, 3> vec2) {
    return std::array<float, 3>({ vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2] });
}

// 3d vector subtraction
std::array<float, 3> sub(std::array<float, 3> vec1, std::array<float, 3> vec2) {
    return std::array<float, 3>({ vec1[0] - vec2[0], vec1[1] - vec2[1], vec1[2] - vec2[2] });
}

// 3d dot product
float dot(std::array<float, 3> vec1, std::array<float, 3> vec2) {
    return (vec1[0] * vec2[0]) + (vec1[1] * vec2[1]) + (vec1[2] * vec2[2]);
}

float calc_distance(std::array<float, 3> vec1, std::array<float, 3> vec2) {
    //(Xa - Xb)² + (Ya - Yb)² + (Za - Zb)²
    
    return sqrt( pow((vec1[0] - vec2[0]),2) + pow((vec1[1] - vec2[1]), 2) + pow((vec1[2] - vec2[2]), 2));
}
// vector norm
float mag(std::array<float, 3> vec) {
    return sqrt(dot(vec, vec));
}

// vector magnitude squared
float mag2(std::array<float, 3> vec) {
    return dot(vec, vec);
}

std::array<float, 3> cross(std::array<float, 3> d1, std::array<float, 3> d2) {
    return std::array<float, 3>({
        (d1[1] * d2[2]) - (d1[2] * d2[1]),
        (d1[2] * d2[0]) - (d1[0] * d2[2]),
        (d1[0] * d2[1]) - (d1[1] * d2[0]) });
}

std::array<float, 3> orthogonal(std::array<float, 3> v1, std::array<float, 3> v2, std::array<float, 3> v3) {
    // adjacent edge vectors
    const std::array<float, 3> d1 = sub(v2, v1);
    const std::array<float, 3> d2 = sub(v3, v2);
    // cross product
    return cross(d1, d2);
};

//unit
std::array<float, 3> unit(std::array<float, 3> vec) {
    return mult(1 / sqrt(mag2(vec)), vec);
}

std::array<float, 3> normal(std::vector<std::array<float, 3>> vertices) {
    // running sum of normal vectors
    std::array<float, 3> normalV({ 0, 0, 0 });

    auto v1 = vertices[0];
    auto v2 = vertices[1];

    for (auto v3 = vertices.begin(); v3 != vertices.end(); v3++) {
        normalV = add(normalV, orthogonal(v1, v2, *v3));
        v1 = v2;
        v2 = *v3;
    } // shift over one
    return unit(normalV);
};


// find first element common to 3 sets by brute force search
uint32_t intersect(std::vector <uint32_t> set1, std::vector <uint32_t> set2, std::vector <uint32_t> set3) {
    for (auto s1 = set1.begin(); s1 != set1.end(); s1++) {
        for (auto s2 = set2.begin(); s2 != set2.end(); s2++) {
            if (*s1 == *s2) {
                for (auto s3 = set3.begin(); s3 != set3.end(); s3++) {
                    if (*s1 == *s3) {
                        return *s1;
                    }
                }
            }
        }
    }

    return 0;
};
