#pragma once

struct CustomTraits : public OpenMesh::DefaultTraits
{
    FaceTraits
    {
        private:
        int32_t tile_id_;

        public:
            const uint32_t tile_id() const { return tile_id_; }
            void set_tile_id(const uint32_t& tile_id) { tile_id_ = tile_id; }
    };
};


typedef OpenMesh::PolyMesh_ArrayKernelT<CustomTraits>  PolygonalMesh;

class PlanetGenerator
{
public:
	PlanetGenerator(std::string planet_name, size_t size);

	void generate();

private:
	std::string planet_name;
	size_t size;
	Polyhedron planet_polyhedron;

private:
    /* OPEN MESH PART */
    void generate_ids();

    std::array<float, 3> get_center(uint32_t ix_face);
    int32_t get_farest_pentagon(std::vector<uint32_t>& pentagons, std::array<float, 3> 
        face_center, std::map<uint32_t, int32_t>& p_index_to_ids);
    int32_t get_nearest_pentagon(std::vector<uint32_t>& pentagons, std::array<float, 3>
        face_center, std::map<uint32_t, int32_t>& p_index_to_ids);
    std::map<uint32_t, int32_t> constructs_pentagons();

    bool is_red_valid(PolygonalMesh& mesh, OpenMesh::SmartFaceHandle face_handle,
        std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p0_p1,
        std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p0_p8,
        std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p1_p8, uint32_t size);

    bool is_green_valid(PolygonalMesh& mesh, OpenMesh::SmartFaceHandle face_handle,
        std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p0_p1,
        std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p0_p8,
        std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p1_p8, uint32_t size);

    bool is_cyan_valid(PolygonalMesh& mesh, OpenMesh::SmartFaceHandle face_handle,
        std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p0_p1,
        std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p0_p8,
        std::vector<OpenMesh::SmartHalfedgeHandle>& ico_p1_p8, uint32_t size);

};

