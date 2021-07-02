#include "master.h"

class Agent{
public:
    float x, y, z;

    Agent(float x, float y, float z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    ~Agent(){}
};

//struct GridCell{
//public:
//    int x, y, z;
//
//    GridCell(int x, int y, int z)
//    {
//        this->x = x;
//        this->y = y;
//        this->z = z;
//    }
//
//    ~GridCell(){}
//
//    bool operator==(const GridCell& other) const
//    {
//        return x == other.x && y == other.y && z == other.z;
//    }
//};
//
//
//
//struct hash_fn{
//public:
//    std::size_t operator()(const GridCell& other) const
//    {
//        std::size_t seed = 0;
//        boost::hash_combine(seed, other.x);
//        boost::hash_combine(seed, other.y);
//        boost::hash_combine(seed, other.z);
//        return seed;
//    }
//};

float randlfloat(float LO, float HI){
    return LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
}

int randint(float LO, float HI){
    return static_cast <int> (rand()) % static_cast <int>(HI - LO + 1) + LO;
}

void hashtable_test(Plate*& plate){
    std::unordered_multimap<GridCell, Agent, hash_fn> AgentMap;

    AgentMap.insert({{1,2,3},{1.1,2.3,3.5}});

    for(auto i=0;i<10;i++){
        GridCell* g = new GridCell(randint(0, plate->x), randint(0, plate->x), randint(plate->agar_height, plate->z));
        Agent* a = new Agent(randlfloat(0, plate->x), randlfloat(0, plate->x), randlfloat(plate->agar_height, plate->z));
        AgentMap.insert({*g, *a});
    }

//    for(auto it : AgentMap){
//        std::cout<<it.first.x<<" "<<it.first.y<<" "<<it.first.z<<" "<<it.second.x<<" "<<it.second.y<<" "<<it.second.z<<std::endl;
//    }

    std::cout<<AgentMap.size()<<std::endl;

    auto range = AgentMap.equal_range({1,2,3});
    for_each (
        range.first,
        range.second,
        [](auto it){std::cout<<it.first.x<<" "<<it.first.y<<" "<<it.first.z<<std::endl<<it.second.x<<" "<<it.second.y<<" "<<it.second.z<<std::endl;}
    );

    AgentMap.clear();
}
