#include "config.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <libpq-fe.h>
#include <map>
#include <pqxx/pqxx>
#include <set>
#include <unistd.h>
#include <vector>

std::string buildstr = "hostaddr=" + host + " dbname=" + dbname + " user=" + user + " password=" + password;

/**
 * 如何实现NN查询
 * 先根据mbb找到前100个，然后从最低的lod开始找，期间要涉及到排序，如果
 */
int N = 1;
int target = 929;
int number = 3;
std::string table1 = "nuclei";
std::string table2 = "vessel";
int distance = 100;
int sampleRate = 10;

class Range {
  public:
    int id;
    float mindis;
    float maxdis;
    Range(){};
    Range(int id_) : id(id_){};
    Range(int id_, float mindis_, float maxdis_) : id(id_), mindis(mindis_), maxdis(maxdis_) {}
    Range(float mindis_, float maxdis_) : mindis(mindis_), maxdis(maxdis_) {}
};

std::string buildQueryMbbSql(int id) {
    char sql[4096];
    if (table2 == "vessel") {
        sprintf(sql,
                "SELECT b.id as id , ST_3DDistance(a.geom, b.geom) as mindis, ST_3DMAXDistance(a.geom, b.geom) as "
                "maxdis FROM "
                "%s_box a, %s_box b WHERE a.id <> b.id AND a.id = %d;",
                table1.c_str(), table2.c_str(), id, distance);
    } else {
        sprintf(sql,
                "SELECT b.id as id, ST_3DDistance(a.geom, b.geom) as mindis, ST_3DMAXDistance(a.geom, b.geom) as "
                "maxdis FROM "
                "%s_box a, %s_box b WHERE a.id <> b.id AND a.id = %d AND ST_3DDWithin(a.geom, "
                "b.geom, %d);",
                table1.c_str(), table2.c_str(), id, distance);
    }
    return std::string(sql);
}
std::string buildIdList(const std::vector<int>& ids) {
    std::string idList;
    for (size_t i = 0; i < ids.size(); ++i) {
        idList += std::to_string(ids[i]);
        if (i < ids.size() - 1) {
            idList += ",";
        }
    }
    return idList;
}

/**
 * 返回的是当前lod下的distance和hausdorff距离
 */
std::string buildQueryOriginSql(int id, std::vector<int> ids) {
    char sql[4096];
    sprintf(sql,
            "SELECT b.id, ST_3DDistance(a.geom, b.geom) as dis FROM "
            "%s_100 a, %s_100 b WHERE a.id <> b.id AND a.id = '%d' AND b.id IN (%s);",
            table1.c_str(), table2.c_str(), id, buildIdList(ids).c_str());
    return std::string(sql);
}

// 根据query mbb返回的结果进行更新
void parseDistanceResult(pqxx::result& rows, std::map<int, Range>& candidates) {
    for (int i = 0; i < rows.size(); i++) {
        for (int j = 0; j < rows[i].size(); j++) {
            candidates[rows[i]["id"].as<int>()] = Range(rows[i]["id"].as<int>());
            candidates[rows[i]["id"].as<int>()].mindis = rows[i]["mindis"].as<float>();
            candidates[rows[i]["id"].as<int>()].maxdis = rows[i]["maxdis"].as<float>();
        }
    }
}

// 根据query lod返回的结果进行更新
void parseLodOriginResult(pqxx::result& rows, std::map<int, Range>& candidates) {
    for (int i = 0; i < rows.size(); i++) {
        for (int j = 0; j < rows[i].size(); j++) {
            candidates[rows[i]["id"].as<int>()].mindis = rows[i]["dis"].as<float>();
            candidates[rows[i]["id"].as<int>()].maxdis = rows[i]["dis"].as<float>();
        }
    }
}

// 根据query mbb返回的结果进行更新
void parseOriginResult(pqxx::result& rows, std::map<int, Range>& candidates) {
    for (int i = 0; i < rows.size(); i++) {
        for (int j = 0; j < rows[i].size(); j++) {
            candidates[rows[i]["id"].as<int>()].mindis = rows[i]["dis"].as<float>();
            candidates[rows[i]["id"].as<int>()].maxdis = rows[i]["dis"].as<float>();
        }
    }
}

std::vector<int> mapKeysToVector(const std::map<int, Range>& myMap) {
    std::vector<int> keys;
    for (const auto& pair : myMap) {
        keys.push_back(pair.first);
    }
    return keys;
}

std::vector<Range> mapValuesToVector(const std::map<int, Range>& myMap) {
    std::vector<Range> values;
    for (const auto& pair : myMap) {
        values.push_back(pair.second);
    }
    return values;
}

/**
 * 过滤掉不符合条件的distance，符合条件的直接假如，待确定的留在其中
 */
void filterByDistance(std::map<int, Range>& ranges, int number) {
    std::set<int> removable;
    std::set<float> window;
    std::vector<Range> rarr = mapValuesToVector(ranges);
    std::sort(rarr.begin(), rarr.end(), [](Range a, Range b) { return a.maxdis < b.maxdis; });
    for (int i = 0; i < number; i++) {
        window.insert(rarr[i].maxdis);
    }
    for (int i = 0; i < rarr.size(); i++) {
        if (rarr[i].mindis > *(window.rbegin())) {
            removable.insert(rarr[i].id);
        }
    }
    for (int id : removable) {
        ranges.erase(id);
    }
}

int main(int argc, char** argv) {
    int opt;
    while ((opt = getopt(argc, argv, "n:1:2:")) != -1) {
        switch (opt) {
            case 'n': N = std::stoi(optarg); break;
            case '1': table1 = optarg; break;
            case '2': table2 = optarg; break;
            default: std::cerr << "Usage: " << argv[0] << " -t <target> -1 <table1> -2 <table2>" << std::endl; return 1;
        }
    }

    pqxx::connection c(buildstr);
    pqxx::work w(c);
    std::vector<std::string> logs;
    logs.reserve(N);
    for (int i = 0; i < N; i += sampleRate) {
        std::string log = "";
        target = i;
        auto beforeTime = std::chrono::steady_clock::now();
        pqxx::result rows = w.exec(buildQueryMbbSql(target));
        std::map<int, Range> candidates;
        parseDistanceResult(rows, candidates);
        filterByDistance(candidates, number);
        auto afterTime = std::chrono::steady_clock::now();
        double duration_millsecond = std::chrono::duration<double, std::milli>(afterTime - beforeTime).count();
        log = log + std::to_string(target) + " " + std::to_string(duration_millsecond);
        rows = w.exec(buildQueryOriginSql(target, mapKeysToVector(candidates)));
        parseOriginResult(rows, candidates);
        filterByDistance(candidates, number);
        auto afterTime2 = std::chrono::steady_clock::now();
        duration_millsecond = std::chrono::duration<double, std::milli>(afterTime2 - afterTime).count();
        log = log + " " + std::to_string(duration_millsecond);
        logs.push_back(log);
    }
    for (int i = 0; i < logs.size(); i++) {
        std::cout << logs[i] << std::endl;
    }
    return 0;
}
