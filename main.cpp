#include "stdafx.h"
#include "Polyhedron.h"

#include "PlanetGenerator.h"

int main()
{
    auto planets = std::vector<std::pair<std::string, uint32_t>>();
    
    // Madis
    planets.push_back(std::make_pair("madis", 32));
    planets.push_back(std::make_pair("madis-moon-1", 7));
    planets.push_back(std::make_pair("madis-moon-2", 8));
    planets.push_back(std::make_pair("madis-moon-3", 11));

    // Alioth
    planets.push_back(std::make_pair("alioth", 93));
    planets.push_back(std::make_pair("alioth-moon-1", 22));
    planets.push_back(std::make_pair("alioth-moon-2", 22));
    planets.push_back(std::make_pair("sanctuary", 61));
    
    // Thades
    planets.push_back(std::make_pair("thades", 45));
    planets.push_back(std::make_pair("thades-moon-1", 10));
    planets.push_back(std::make_pair("thades-moon-2", 11));

    // Talemai
    planets.push_back(std::make_pair("talemai", 42));
    planets.push_back(std::make_pair("talemai-moon-1", 11));
    planets.push_back(std::make_pair("talemai-moon-2", 8));
    planets.push_back(std::make_pair("talemai-moon-3", 8));

    // Feli
    planets.push_back(std::make_pair("feli", 44));
    planets.push_back(std::make_pair("feli-moon-1", 10));

    // Sicari
    planets.push_back(std::make_pair("sicari", 37));

    // Sinnen
    planets.push_back(std::make_pair("sinnen", 40));
    planets.push_back(std::make_pair("sinnen-moon-1", 12));

    // Teoma
    planets.push_back(std::make_pair("teoma", 45));

    // Jago
    planets.push_back(std::make_pair("jago", 46));

    // Lacobus
    planets.push_back(std::make_pair("lacobus", 42));
    planets.push_back(std::make_pair("lacobus-moon-1", 13));
    planets.push_back(std::make_pair("lacobus-moon-2", 10));
    planets.push_back(std::make_pair("lacobus-moon-3", 11));

    // Symeon
    planets.push_back(std::make_pair("symeon", 36));
    
    // Ion
    planets.push_back(std::make_pair("ion", 42));
    planets.push_back(std::make_pair("ion-moon-1", 8));
    planets.push_back(std::make_pair("ion-moon-2", 11));

    std::cout << "Planets/moons count: " << planets.size() << std::endl;

/*
#ifndef _DEBUG
    for (auto planet_it = planets.begin(); planet_it != planets.end(); planet_it++) {
        std::cout << "Starting generating " << planet_it->first << std::endl;

        PlanetGenerator generator(planet_it->first, planet_it->second);
        generator.generate();
    }
#else
    PlanetGenerator generator(planets[1].first, planets[1].second);
    generator.generate();
#endif
*/
    for (auto i = 1; i < 94; i++) {
        PlanetGenerator generator(format("mn-%d", i), i);
        generator.generate();
    }

    return 0;
}