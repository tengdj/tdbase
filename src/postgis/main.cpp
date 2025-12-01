#include "config.h"
#include <boost/program_options.hpp>
#include <chrono>
#include <iostream>
#include <libpq-fe.h>
#include <pqxx/pqxx>
namespace po = boost::program_options;
std::string buildstr = "hostaddr=" + host + " dbname=" + dbname + " user=" + user + " password=" + password;

/**
 * 需要takein的参数
 * 1. 数据的数量
 * 2. table2
 * 3. sample rate // 针对的是nv all
 * 4. 查询的类型 within or nn
 * 5. progressive or origin
 */
// parameter
int N = 1;
int number = 3;
int distance = 50;
int initDistance = 100;
int sampleRate = 1;
std::string table1 = "nuclei";
std::string table2 = "nuclei";
std::string queryType = "within";
std::string queryMethod = "progressive";
//
pqxx::connection c(buildstr);
pqxx::work w(c);
std::vector<std::string> logs;

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

std::string buildQueryWithinMbbSql(int id) {
    char sql[4096];
    sprintf(sql,
            "SELECT b.id as id FROM %s_box a, %s_box b WHERE a.id <> b.id AND a.id = %d AND ST_3DDWithin(a.geom, "
            "b.geom, %d);",
            table1.c_str(), table2.c_str(), id, distance);
    return std::string(sql);
}

std::string buildQueryNNMbbSql(int id) {
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
                table1.c_str(), table2.c_str(), id, initDistance);
    }
    return std::string(sql);
}

std::string buildQueryOriginSql(int id, std::vector<int> ids) {
    char sql[4096];
    sprintf(sql,
            "SELECT b.id, ST_3DDistance(a.geom, b.geom) as dis FROM "
            "%s_100 a, %s_100 b WHERE a.id <> b.id AND a.id = '%d' AND b.id IN (%s);",
            table1.c_str(), table2.c_str(), id, buildIdList(ids).c_str());
    return std::string(sql);
}

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

void parseOriginResult(pqxx::result& rows, std::map<int, Range>& candidates) {
    for (int i = 0; i < rows.size(); i++) {
        for (int j = 0; j < rows[i].size(); j++) {
            candidates[rows[i]["id"].as<int>()].mindis = rows[i]["dis"].as<float>();
            candidates[rows[i]["id"].as<int>()].maxdis = rows[i]["dis"].as<float>();
        }
    }
}

void parseLodDistanceResult(pqxx::result& rows, std::map<int, Range>& candidates, std::pair<float, float> item) {
    for (int i = 0; i < rows.size(); i++) {
        for (int j = 0; j < rows[i].size(); j++) {
            candidates[rows[i]["id"].as<int>()].mindis = rows[i]["dis"].as<float>() - rows[i]["phausdorff"].as<float>() - item.second;
            candidates[rows[i]["id"].as<int>()].maxdis = rows[i]["dis"].as<float>() + rows[i]["hausdorff"].as<float>() + item.first;
        }
    }
}

void filterWithinByDistance(std::map<int, Range>& ranges, std::vector<int>& result, int distance) {
    for (auto it = ranges.begin(); it != ranges.end();) {
        if ((*it).second.maxdis <= distance) {
            result.push_back((*it).first);
            it = ranges.erase(it);
        } else if ((*it).second.mindis > distance) {
            it = ranges.erase(it);
        } else {
            it++;
        }
    }
}

void filterNNByDistance(std::map<int, Range>& ranges, std::vector<int>& result) {
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

void parseWithinMbbResult(pqxx::result& rows, std::map<int, Range>& candidates) {
    for (int i = 0; i < rows.size(); i++) {
        for (int j = 0; j < rows[i].size(); j++) {
            candidates[rows[i]["id"].as<int>()] = Range();
        }
    }
}

void parseNNMbbResult(pqxx::result& rows, std::map<int, Range>& candidates) {
    for (int i = 0; i < rows.size(); i++) {
        for (int j = 0; j < rows[i].size(); j++) {
            candidates[rows[i]["id"].as<int>()] = Range(rows[i]["id"].as<int>());
            candidates[rows[i]["id"].as<int>()].mindis = rows[i]["mindis"].as<float>();
            candidates[rows[i]["id"].as<int>()].maxdis = rows[i]["maxdis"].as<float>();
        }
    }
}

void withinProgressive(int target) {
    std::vector<int> candidateNumber;
    std::vector<double> iterTimes;
    auto beforeTime = std::chrono::steady_clock::now();
    pqxx::result rows = w.exec(buildQueryWithinMbbSql(target));
    std::map<int, Range> candidates;
    parseWithinMbbResult(rows, candidates);
    std::vector<int> result;
    std::string log = std::to_string(target) + " ";
    auto mbbTime = std::chrono::steady_clock::now();
    double duration_millsecond = std::chrono::duration<double, std::milli>(mbbTime - beforeTime).count();
    log = log + std::to_string(duration_millsecond) + " ";
    for (int lod = 20; lod <= 100; lod += 20) {
        if (candidates.empty()) {
            break;
        }
        candidateNumber.push_back(candidates.size());
        auto iterSt = std::chrono::steady_clock::now();
        rows = w.exec(buildQueryHausdorffSql(lod, target));
        std::pair<float, float> targetHausdorff = std::make_pair(rows[0]["hausdorff"].as<float>(), rows[0]["phausdorff"].as<float>());
        rows = w.exec(buildQueryLodSql(lod, target, mapKeysToVector(candidates)));
        parseLodDistanceResult(rows, candidates, targetHausdorff);
        filterWithinByDistance(candidates, result, distance);
        auto iterEd = std::chrono::steady_clock::now();
        double ts = std::chrono::duration<double, std::milli>(iterEd - iterSt).count();
        iterTimes.push_back(ts);
    }
    for (int i = 0; i < 5; i++) {
        if (i < candidateNumber.size()) {
            log = log + std::to_string(candidateNumber[i]) + " ";
        } else {
            log = log + std::to_string(0) + " ";
        }
    }
    for (int i = 0; i < 5; i++) {
        if (i < iterTimes.size()) {
            log = log + std::to_string(iterTimes[i]) + " ";
        } else {
            log = log + std::to_string(0) + " ";
        }
    }
    auto endTime = std::chrono::steady_clock::now();
    double allTime = std::chrono::duration<double, std::milli>(endTime - beforeTime).count();
    log += std::to_string(allTime);
    logs.push_back(log);
}

void withinOrigin(int target) {
    std::string log = "";
    auto beforeTime = std::chrono::steady_clock::now();
    pqxx::result rows = w.exec(buildQueryWithinMbbSql(target));
    std::map<int, Range> candidates;
    parseWithinMbbResult(rows, candidates);
    std::vector<int> result;
    auto afterTime = std::chrono::steady_clock::now();
    double duration_millsecond = std::chrono::duration<double, std::milli>(afterTime - beforeTime).count();
    log = log + std::to_string(target) + " " + std::to_string(duration_millsecond);

    if (candidates.empty()) {
        log = log + " " + std::to_string(0);
        logs.push_back(log);
        return;
    }
    rows = w.exec(buildQueryOriginSql(target, mapKeysToVector(candidates)));
    parseOriginResult(rows, candidates);
    filterWithinByDistance(candidates, result, distance);
    auto afterTime2 = std::chrono::steady_clock::now();
    duration_millsecond = std::chrono::duration<double, std::milli>(afterTime2 - afterTime).count();
    log = log + " " + std::to_string(duration_millsecond);
    logs.push_back(log);
}

void nnProgressive(int target) {
    std::vector<int> candidateNumber;
    std::vector<double> iterTimes;
    auto beforeTime = std::chrono::steady_clock::now();
    pqxx::result rows = w.exec(buildQueryNNMbbSql(target));
    std::map<int, Range> candidates;
    std::vector<int> result;
    parseNNMbbResult(rows, candidates);
    filterNNByDistance(candidates, result);
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
        filterNNByDistance(candidates, result);
        if (candidates.size() <= number - result.size()) {
            for (auto it : candidates) {
                result.push_back(it.first);
            }
        }
        auto iterEd = std::chrono::steady_clock::now();
        double ts = std::chrono::duration<double, std::milli>(iterEd - iterSt).count();
        iterTimes.push_back(ts);
        if (result.size() == number) {
            break;
        }
    }
    for (int i = 0; i < 5; i++) {
        if (i < candidateNumber.size()) {
            log = log + std::to_string(candidateNumber[i]) + " ";
        } else {
            log = log + std::to_string(0) + " ";
        }
    }
    for (int i = 0; i < 5; i++) {
        if (i < iterTimes.size()) {
            log = log + std::to_string(iterTimes[i]) + " ";
        } else {
            log = log + std::to_string(0) + " ";
        }
    }
    auto endTime = std::chrono::steady_clock::now();
    double allTime = std::chrono::duration<double, std::milli>(endTime - beforeTime).count();
    log += std::to_string(allTime);
    logs.push_back(log);
}

void nnOrigin(int target) {
    std::string log = "";
    std::vector<int> result;
    auto beforeTime = std::chrono::steady_clock::now();
    pqxx::result rows = w.exec(buildQueryNNMbbSql(target));
    std::map<int, Range> candidates;
    parseNNMbbResult(rows, candidates);
    filterNNByDistance(candidates, result);
    auto afterTime = std::chrono::steady_clock::now();
    double duration_millsecond = std::chrono::duration<double, std::milli>(afterTime - beforeTime).count();
    log = log + std::to_string(target) + " " + std::to_string(duration_millsecond);
    rows = w.exec(buildQueryOriginSql(target, mapKeysToVector(candidates)));
    parseOriginResult(rows, candidates);
    filterNNByDistance(candidates, result);
    auto afterTime2 = std::chrono::steady_clock::now();
    duration_millsecond = std::chrono::duration<double, std::milli>(afterTime2 - afterTime).count();
    log = log + " " + std::to_string(duration_millsecond);
    logs.push_back(log);
}

int main(int argc, char** argv) {
    po::options_description desc("exp setting");
    desc.add_options()("help,h", "produce help message")("num", po::value<int>(&N), "query number")(
        "samplerate", po::value<int>(&sampleRate), "used for nv origin")("table2", po::value<std::string>(&table2), "the second table")(
        "queryType", po::value<std::string>(&queryType), "within or nn")("queryMethod", po::value<std::string>(&queryMethod),
                                                                         "progressive or origin");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    logs.reserve(N);
    for (int i = 0; i < N; i += sampleRate) {
        if (queryType == "within" && queryMethod == "progressive") {
            withinProgressive(i);
        } else if (queryType == "within" && queryMethod == "origin") {
            withinOrigin(i);
        } else if (queryType == "nn" && queryMethod == "progressive") {
            nnProgressive(i);
        } else if (queryType == "nn" && queryMethod == "origin") {
            nnOrigin(i);
        }
        std::cout << *logs.rbegin() << std::endl;
    }
    return 0;
}
