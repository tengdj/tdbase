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
int target = 6853;
int N = 10;
int number = 3;
std::string table1 = "nuclei";
std::string table2 = "nuclei";
int distance = 100;
std::vector<int> result;
std::vector<int> candidateNumber;
std::vector<double> iterTimes;
int sampleRate = 1;

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
std::string buildQueryLodSql(int lod, int id, std::vector<int> ids) {
    char sql[4096];
    sprintf(sql,
            "SELECT b.id, ST_3DDistance(a.geom, b.geom) as dis, "
            "b.hausdorff, b.phausdorff FROM "
            "%s_%d a, %s_%d b WHERE a.id <> b.id AND a.id = '%d' AND b.id IN (%s);",
            table1.c_str(), lod, table2.c_str(), lod, id, buildIdList(ids).c_str());
    return std::string(sql);
}

std::string buildQueryHausdorffSql(int lod, int id) {
    char sql[4096];
    sprintf(sql,
            "SELECT a.hausdorff, a.phausdorff "
            "FROM %s_%d a "
            "WHERE a.id = '%d' ",
            table1.c_str(), lod, id);
    return std::string(sql);
}

// 根据query mbb返回的结果进行更新
void parseDistanceResult(pqxx::result& rows, std::map<int, Range>& candidates) {
    for (int i = 0; i < rows.size(); i++) {
        for (int j = 0; j < rows[i].size(); j++) {
            candidates[rows[i]["id"].as<int>()] = Range(rows[i]["id"].as<int>(), rows[i]["mindis"].as<float>(), rows[i]["maxdis"].as<float>());
        }
    }
}

// 根据query lod返回的结果进行更新
void parseLodDistanceResult(pqxx::result& rows, std::map<int, Range>& candidates, std::pair<float, float> item) {
    for (int i = 0; i < rows.size(); i++) {
        for (int j = 0; j < rows[i].size(); j++) {
            candidates[rows[i]["id"].as<int>()].mindis = rows[i]["dis"].as<float>() - rows[i]["phausdorff"].as<float>() - item.second;
            candidates[rows[i]["id"].as<int>()].maxdis = rows[i]["dis"].as<float>() + rows[i]["hausdorff"].as<float>() + item.first;
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
void filterByDistance(std::map<int, Range>& ranges) {
    std::set<int> removable;
    std::set<float> window;
    std::vector<Range> rarr = mapValuesToVector(ranges);
    std::sort(rarr.begin(), rarr.end(), [](Range a, Range b) { return a.maxdis < b.maxdis; });
    // std::sort(rarr.begin(), rarr.end(), [](Range a, Range b) { return a.mindis < b.mindis; });
    // 如果有人的maxdis小于其他所有人的mindis，则说明可以直接确定下来
    for (int i = 0; i < number - result.size(); i++) {
        int flag = true;
        for (int j = i + 1; j < rarr.size(); j++) {
            if (rarr[i].maxdis > rarr[j].mindis) {
                flag = false;
                break;
            }
        }
        if (flag) {
            result.push_back(rarr[i].id);
            removable.insert(rarr[i].id);
        }
    }
    for (int id : removable) {
        for (auto it = rarr.begin(); it != rarr.end();) {
            if (id == (*it).id) {
                it = rarr.erase(it);
            } else {
                it++;
            }
        }
        ranges.erase(id);
    }
    if (result.size() == number) {
        return;
    }
    removable.clear();
    for (int i = 0; i < number - result.size(); i++) {
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
    std::vector<int> avgNumber(5, 0);
    std::vector<double> avgTime(5, 0);
    for (int i = 0; i < N; i++) {
        target = i;
        result.clear();
        result.shrink_to_fit();
        candidateNumber.clear();
        candidateNumber.shrink_to_fit();
        iterTimes.clear();
        iterTimes.shrink_to_fit();
        auto beforeTime = std::chrono::steady_clock::now();
        pqxx::result rows = w.exec(buildQueryMbbSql(target));
        std::map<int, Range> candidates;
        parseDistanceResult(rows, candidates);
        filterByDistance(candidates);
        std::string log = std::to_string(target) + " ";

        auto mbbTime = std::chrono::steady_clock::now();
        double duration_millsecond = std::chrono::duration<double, std::milli>(mbbTime - beforeTime).count();
        log = log + std::to_string(duration_millsecond) + " ";
        for (int lod = 20; lod <= 100; lod += 20) {
            candidateNumber.push_back(candidates.size());
            auto iterSt = std::chrono::steady_clock::now();
            rows = w.exec(buildQueryHausdorffSql(lod, target));
            std::pair<float, float> targetHausdorff = std::make_pair(rows[0]["hausdorff"].as<float>(), rows[0]["phausdorff"].as<float>());
            rows = w.exec(buildQueryLodSql(lod, target, mapKeysToVector(candidates)));
            parseLodDistanceResult(rows, candidates, targetHausdorff);
            filterByDistance(candidates);
            if (candidates.size() <= number - result.size()) {
                for (auto it : candidates) {
                    result.push_back(it.first);
                }
            }
            auto iterEd = std::chrono::steady_clock::now();
            double ts = std::chrono::duration<double, std::milli>(iterEd - iterSt).count();
            iterTimes.push_back(ts);
            if (result.size() == number) {
                // for (int i = 0; i < result.size(); i++) {
                //     std::cout << result[i] << " ";
                // }
                // std::cout << "\n";
                break;
            }
        }
        for (int i = 0; i < 5; i++) {
            if (i < candidateNumber.size()) {
                log = log + std::to_string(candidateNumber[i]) + " ";
                avgNumber[i] += candidateNumber[i];
            } else {
                log = log + std::to_string(0) + " ";
            }
        }
        for (int i = 0; i < 5; i++) {
            if (i < iterTimes.size()) {
                log = log + std::to_string(iterTimes[i]) + " ";
                avgTime[i] += iterTimes[i];
            } else {
                log = log + std::to_string(0) + " ";
            }
        }
        auto endTime = std::chrono::steady_clock::now();
        double allTime = std::chrono::duration<double, std::milli>(endTime - beforeTime).count();
        log += std::to_string(allTime);
        logs.push_back(log);
    }
    for (int i = 0; i < logs.size(); i++) {
        std::cout << logs[i] << std::endl;
    }
    for (int i = 0; i < avgNumber.size(); i++) {
        std::cout << (double)avgNumber[i] / (double)N << " ";
    }
    std::cout << "\n";
    for (int i = 0; i < avgTime.size(); i++) {
        std::cout << avgTime[i] / (double)N << " ";
    }
    std::cout << "\n";
    return 0;
}
