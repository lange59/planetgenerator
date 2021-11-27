#pragma once

#include <vector>
#include <array>
#include <sstream>
#include <string>
#include <fstream>
#include <map>

#include <cstdarg>
#include <iostream>
#include <string.h>

class Polyhedron
{
public:
	Polyhedron();

	std::vector<std::vector<uint32_t>> faces;
	std::vector<std::array<float, 3>> vertices;
	std::vector<std::array<float, 3>> _normals;

	std::vector< std::array<float, 3>> centers();
	std::vector< std::array<float, 3>> normals();

	void to_OBJ(std::string output_path);
	void to_DU(std::string output_path);
	
	static Polyhedron from_OBJ(std::string input_path);
	static Polyhedron from_DU(std::string input_path);
};

class Polyflag
{
public:
	Polyflag();

	void newV(std::string vertName, std::array<float, 3> coordinates);
	void newFlag(std::string faceName, std::string vertName1, std::string vertName2);

	Polyhedron topoly();

	std::map<std::string, std::map<std::string, std::string>> flags;
	std::map<std::string, uint32_t> vertidxs;
	std::map<std::string, std::array<float, 3>> vertices;
};

// Conway Operators
void apply_subdivide(Polyhedron& poly, uint32_t n);
void apply_dual(Polyhedron& poly);
void apply_kis(Polyhedron& poly, uint32_t n = 0, float apexdist = 0.1);
void apply_inflate(Polyhedron& poly);


std::string format(const char* fmt, ...);

// 3d scalar multiplication
std::array<float, 3> mult(float c, std::array<float, 3> vec);

// 3d vector addition
std::array<float, 3> add(std::array<float, 3> vec1, std::array<float, 3> vec2);

// 3d vector subtraction
std::array<float, 3> sub(std::array<float, 3> vec1, std::array<float, 3> vec2);

// 3d dot product
float dot(std::array<float, 3> vec1, std::array<float, 3> vec2);

// vector norm
float mag(std::array<float, 3> vec);

// vector magnitude squared
float mag2(std::array<float, 3> vec);

float calc_distance(std::array<float, 3> vec1, std::array<float, 3> vec2);

std::array<float, 3> cross(std::array<float, 3> d1, std::array<float, 3> d2);

std::array<float, 3> orthogonal(std::array<float, 3> v1, std::array<float, 3> v2, std::array<float, 3> v3);

//unit
std::array<float, 3> unit(std::array<float, 3> vec);

std::array<float, 3> normal(std::vector<std::array<float, 3>> vertices);

// find first element common to 3 sets by brute force search
uint32_t intersect(std::vector <uint32_t> set1, std::vector <uint32_t> set2, std::vector <uint32_t> set3);