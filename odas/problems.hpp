#ifndef __ODAS_PROBLEMS
#define __ODAS_PROBLEMS

#include <vector>

class Problems
{
private:
    std::vector<std::string> namelist;
    
public:
    Problems()
    {
        // heritage buildings
        namelist.push_back("data/HK_2424275_HKU_Main_Building.pcd");
        namelist.push_back("data/Dublin_25631284_Dublin_City_Hall.pcd");
        namelist.push_back("data/HK_181773526_HKU_HHY_Building.pcd");
        
        // modern buildings
        //namelist.push_back("data/HK_540115808_118_Connaught_Rd_W.pcd");
        namelist.push_back("data/Dublin_177562798_One_Georges_Quay_Plaza.pcd");
        namelist.push_back("data/Dublin_4560539_47_51_OConnell_Street_Upper.pcd");
        namelist.push_back("data/HK_26305856_Fruits_Wholesale_Market.pcd");
        
        // bridges & pier
        namelist.push_back("data/Dublin_46160694_Samuel_Beckett_Bridge.pcd");
        namelist.push_back("data/Dublin_4934685_Sean_OCasey_Bridge.pcd");
        namelist.push_back("data/HK_Two_Piers.pcd");
        
        // park - open space
        //namelist.push_back("data/Dublin_3885916_Trinity_College_Library_Square.pcd");
        
        // area of buildings
        //namelist.push_back("data/Goolwa_South_Australia_Gardiner_St.pcd");
        // namelist.push_back("data/Fremont_CA_Ondina_Drive.pcd");
        
        // natural landscape
        //namelist.push_back("data/Hawaii_Lua_Hohonu.pcd");
        //namelist.push_back("data/Willington_Kapiti_Island.pcd");
        
        // namelist.push_back("data/HK_320984029_Kwan_Yick_Building_Phase_I.pcd");
        // namelist.push_back("data/HK_155738272_Manhattan_Heights.pcd");
        // namelist.push_back("data/HK_320984071_72_The_Belchers_Blocks_1_2.pcd");
        // namelist.push_back("data/Dublin_113818887_Saint_Andrews_Church.pcd");
        // namelist.push_back("data/Dublin_294790295.pcd");
        // namelist.push_back("data/HK_26201751_Belcher_Bay_Park.pcd");
    }
    
    std::size_t get_problem_size()
    {
        return namelist.size();
    }
    
    std::string get_filename(const int id)
    {
        return namelist[id % namelist.size()];
    }
};

#endif

